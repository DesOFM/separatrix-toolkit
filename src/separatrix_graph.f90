module separatrix_graph
    use iso_fortran_env, only: int64
    use kinds, only: dp
    use separatrix_types, only: saddle_point, branch, separatrix_network
    use contour_marching, only: segment_set
    implicit none
    private

    public :: build_networks
    public :: append_network_array

    type :: graph_data
        integer :: nv=0, ne=0
        real(dp), allocatable :: x(:),y(:)
        integer, allocatable :: e1(:),e2(:),degree(:),saddle_id(:)
    end type graph_data

    type :: vertex_hash
        integer :: nbucket=0
        real(dp) :: cell=1.0_dp
        integer, allocatable :: head(:),next(:)
        integer(int64), allocatable :: qx(:),qy(:)
    end type vertex_hash

    type :: edge_hash
        integer :: nbucket=0
        integer, allocatable :: head(:),next(:)
    end type edge_hash

contains

    subroutine build_networks(seg,saddles,tol,networks)
        type(segment_set), intent(in) :: seg
        type(saddle_point), intent(in) :: saddles(:)
        real(dp), intent(in) :: tol
        type(separatrix_network), allocatable, intent(out) :: networks(:)
        type(graph_data) :: g
        integer, allocatable :: vcomp(:),ecomp(:)
        integer :: ncomp,c,nnet
        logical :: has_saddle
        type(separatrix_network), allocatable :: temp(:)

        call segments_to_graph(seg,saddles,tol,g)
        call graph_components(g,vcomp,ecomp,ncomp)
        allocate(temp(max(1,ncomp)));nnet=0
        do c=1,ncomp
            has_saddle=any(g%saddle_id(1:g%nv)>0.and.vcomp==c)
            if(.not.has_saddle)cycle
            nnet=nnet+1;temp(nnet)%id=nnet
            call extract_component_branches(g,c,vcomp,ecomp,temp(nnet))
        end do
        allocate(networks(nnet));if(nnet>0)networks=temp(1:nnet)
        deallocate(temp,vcomp,ecomp)
    end subroutine build_networks


    subroutine segments_to_graph(seg,saddles,tol,g)
        type(segment_set), intent(in) :: seg
        type(saddle_point), intent(in) :: saddles(:)
        real(dp), intent(in) :: tol
        type(graph_data), intent(out) :: g
        type(vertex_hash) :: vh
        type(edge_hash) :: eh
        integer :: maxv,s,v1,v2,k,bestv
        real(dp) :: d,best

        maxv=max(2,2*seg%n)
        allocate(g%x(maxv),g%y(maxv),g%e1(max(1,seg%n)),g%e2(max(1,seg%n)), &
                 g%degree(maxv),g%saddle_id(maxv))
        g%nv=0;g%ne=0;g%degree=0;g%saddle_id=0
        call init_vertex_hash(vh,maxv,tol)
        call init_edge_hash(eh,max(1,seg%n))

        do s=1,seg%n
            v1=find_or_add_vertex(g,vh,seg%x1(s),seg%y1(s),tol)
            v2=find_or_add_vertex(g,vh,seg%x2(s),seg%y2(s),tol)
            if(v1==v2)cycle
            if(edge_exists(eh,g,v1,v2))cycle
            g%ne=g%ne+1;g%e1(g%ne)=v1;g%e2(g%ne)=v2
            call insert_edge_hash(eh,v1,v2,g%ne)
        end do

        do s=1,g%ne
            g%degree(g%e1(s))=g%degree(g%e1(s))+1;g%degree(g%e2(s))=g%degree(g%e2(s))+1
        end do

        do k=1,size(saddles)
            best=huge(1.0_dp);bestv=0
            do s=1,g%nv
                d=hypot(g%x(s)-saddles(k)%x,g%y(s)-saddles(k)%y)
                if(d<best)then;best=d;bestv=s;end if
            end do
            if(bestv>0.and.best<=10.0_dp*tol)g%saddle_id(bestv)=saddles(k)%id
        end do
    end subroutine segments_to_graph


    subroutine init_vertex_hash(h,maxv,tol)
        type(vertex_hash), intent(out) :: h
        integer, intent(in) :: maxv
        real(dp), intent(in) :: tol
        h%nbucket=next_odd(4*maxv+1);h%cell=max(tol,tiny(1.0_dp))
        allocate(h%head(h%nbucket),h%next(maxv),h%qx(maxv),h%qy(maxv))
        h%head=0;h%next=0;h%qx=0_int64;h%qy=0_int64
    end subroutine init_vertex_hash

    subroutine init_edge_hash(h,maxe)
        type(edge_hash), intent(out) :: h
        integer, intent(in) :: maxe
        h%nbucket=next_odd(4*maxe+1)
        allocate(h%head(h%nbucket),h%next(maxe));h%head=0;h%next=0
    end subroutine init_edge_hash

    integer function find_or_add_vertex(g,h,x,y,tol) result(v)
        type(graph_data), intent(inout) :: g
        type(vertex_hash), intent(inout) :: h
        real(dp), intent(in) :: x,y,tol
        integer(int64) :: qx,qy
        integer :: dx,dy,b,i

        qx=int(floor(x/h%cell),int64);qy=int(floor(y/h%cell),int64)
        do dx=-1,1
            do dy=-1,1
                b=hash_cell(qx+dx,qy+dy,h%nbucket);i=h%head(b)
                do while(i/=0)
                    if((g%x(i)-x)**2+(g%y(i)-y)**2<=tol*tol)then;v=i;return;end if
                    i=h%next(i)
                end do
            end do
        end do

        g%nv=g%nv+1;if(g%nv>size(g%x))error stop 'Vertex buffer exceeded'
        v=g%nv;g%x(v)=x;g%y(v)=y;h%qx(v)=qx;h%qy(v)=qy
        b=hash_cell(qx,qy,h%nbucket);h%next(v)=h%head(b);h%head(b)=v
    end function find_or_add_vertex

    logical function edge_exists(h,g,v1,v2) result(found)
        type(edge_hash), intent(in) :: h
        type(graph_data), intent(in) :: g
        integer, intent(in) :: v1,v2
        integer :: b,e,a,c
        a=min(v1,v2);c=max(v1,v2);b=hash_edge(a,c,h%nbucket);e=h%head(b);found=.false.
        do while(e/=0)
            if(min(g%e1(e),g%e2(e))==a.and.max(g%e1(e),g%e2(e))==c)then;found=.true.;return;end if
            e=h%next(e)
        end do
    end function edge_exists

    subroutine insert_edge_hash(h,v1,v2,e)
        type(edge_hash), intent(inout) :: h
        integer, intent(in) :: v1,v2,e
        integer :: b
        b=hash_edge(min(v1,v2),max(v1,v2),h%nbucket);h%next(e)=h%head(b);h%head(b)=e
    end subroutine insert_edge_hash

    integer function hash_cell(qx,qy,n) result(h)
        integer(int64), intent(in) :: qx,qy
        integer, intent(in) :: n
        integer(int64) :: z
        z=qx*73856093_int64+qy*19349663_int64
        h=int(modulo(z,int(n,int64)))+1
    end function hash_cell

    integer function hash_edge(a,b,n) result(h)
        integer, intent(in) :: a,b,n
        integer(int64) :: z
        z=int(a,int64)*73856093_int64+int(b,int64)*19349663_int64
        h=int(modulo(z,int(n,int64)))+1
    end function hash_edge

    integer function next_odd(n) result(v)
        integer, intent(in) :: n
        v=max(17,n);if(mod(v,2)==0)v=v+1
    end function next_odd

    subroutine graph_components(g,vcomp,ecomp,ncomp)
        type(graph_data), intent(in) :: g
        integer, allocatable, intent(out) :: vcomp(:),ecomp(:)
        integer, intent(out) :: ncomp
        integer, allocatable :: queue(:)
        integer :: v,head,tail,e,w,nbr
        allocate(vcomp(g%nv),ecomp(g%ne),queue(g%nv));vcomp=0;ecomp=0;ncomp=0
        do v=1,g%nv
            if(vcomp(v)/=0)cycle
            ncomp=ncomp+1;head=1;tail=1;queue(1)=v;vcomp(v)=ncomp
            do while(head<=tail)
                w=queue(head);head=head+1
                do e=1,g%ne
                    if(g%e1(e)==w.or.g%e2(e)==w)then
                        ecomp(e)=ncomp
                        if(g%e1(e)==w)then;nbr=g%e2(e);else;nbr=g%e1(e);end if
                        if(vcomp(nbr)==0)then;tail=tail+1;queue(tail)=nbr;vcomp(nbr)=ncomp;end if
                    end if
                end do
            end do
        end do
    end subroutine graph_components

    subroutine extract_component_branches(g,comp,vcomp,ecomp,network)
        type(graph_data), intent(in) :: g
        integer, intent(in) :: comp,vcomp(:),ecomp(:)
        type(separatrix_network), intent(inout) :: network
        logical, allocatable :: used(:)
        type(branch), allocatable :: temp(:)
        integer :: v,e,nb
        allocate(used(g%ne));used=.false.;allocate(temp(max(1,count(ecomp==comp))));nb=0
        do v=1,g%nv
            if(vcomp(v)/=comp)cycle
            if(.not.(g%degree(v)/=2.or.g%saddle_id(v)>0))cycle
            do e=1,g%ne
                if(ecomp(e)/=comp.or.used(e))cycle
                if(g%e1(e)==v.or.g%e2(e)==v)then
                    nb=nb+1;call trace_branch(g,comp,ecomp,used,v,e,temp(nb));temp(nb)%id=nb
                end if
            end do
        end do
        do e=1,g%ne
            if(ecomp(e)==comp.and..not.used(e))then
                nb=nb+1;call trace_branch(g,comp,ecomp,used,g%e1(e),e,temp(nb));temp(nb)%id=nb
            end if
        end do
        network%nbranches=nb;allocate(network%branches(nb));if(nb>0)network%branches=temp(1:nb)
        deallocate(temp,used)
    end subroutine extract_component_branches

    subroutine trace_branch(g,comp,ecomp,used,start_v,start_e,b)
        type(graph_data), intent(in) :: g
        integer, intent(in) :: comp,ecomp(:),start_v,start_e
        logical, intent(inout) :: used(:)
        type(branch), intent(out) :: b
        integer, allocatable :: verts(:)
        integer :: current_v,current_e,next_v,next_e,nv_path,e
        allocate(verts(g%ne+2));nv_path=1;verts(1)=start_v;current_v=start_v;current_e=start_e
        do
            used(current_e)=.true.
            if(g%e1(current_e)==current_v)then;next_v=g%e2(current_e);else;next_v=g%e1(current_e);end if
            nv_path=nv_path+1;verts(nv_path)=next_v
            if(next_v==start_v.and.nv_path>2)exit
            if(g%saddle_id(next_v)>0.or.g%degree(next_v)/=2)exit
            next_e=0
            do e=1,g%ne
                if(ecomp(e)/=comp.or.used(e))cycle
                if(g%e1(e)==next_v.or.g%e2(e)==next_v)then;next_e=e;exit;end if
            end do
            if(next_e==0)exit
            current_v=next_v;current_e=next_e
        end do
        b%n=nv_path;allocate(b%x(nv_path),b%y(nv_path))
        do e=1,nv_path;b%x(e)=g%x(verts(e));b%y(e)=g%y(verts(e));end do
        b%saddle_start=g%saddle_id(verts(1));b%saddle_end=g%saddle_id(verts(nv_path));deallocate(verts)
    end subroutine trace_branch

    subroutine append_network_array(master,newnets)
        type(separatrix_network), allocatable, intent(inout) :: master(:)
        type(separatrix_network), intent(in) :: newnets(:)
        type(separatrix_network), allocatable :: temp(:)
        integer :: nold,nnew,i
        nold=0;if(allocated(master))nold=size(master);nnew=size(newnets);if(nnew==0)return
        allocate(temp(nold+nnew));if(nold>0)temp(1:nold)=master;temp(nold+1:nold+nnew)=newnets
        do i=1,nold+nnew;temp(i)%id=i;end do
        call move_alloc(temp,master)
    end subroutine append_network_array

end module separatrix_graph

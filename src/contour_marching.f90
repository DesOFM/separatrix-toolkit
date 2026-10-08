module contour_marching
    use kinds, only: dp
    use separatrix_types, only: grid_field
    implicit none
    private

    public :: segment_set
    public :: build_level_segments
    public :: clip_segments_outside_circle
    public :: add_segment

    real(dp), parameter :: tiny_value = 1.0e-14_dp

    type :: segment_set
        integer :: n = 0
        integer :: capacity = 0
        real(dp), allocatable :: x1(:), y1(:), x2(:), y2(:)
    end type segment_set

contains

    subroutine init_segments(seg, capacity)
        type(segment_set), intent(inout) :: seg
        integer, intent(in) :: capacity

        seg%n = 0
        seg%capacity = max(16, capacity)
        if (allocated(seg%x1)) deallocate(seg%x1, seg%y1, seg%x2, seg%y2)
        allocate(seg%x1(seg%capacity), seg%y1(seg%capacity), &
                 seg%x2(seg%capacity), seg%y2(seg%capacity))
    end subroutine init_segments


    subroutine ensure_capacity(seg, needed)
        type(segment_set), intent(inout) :: seg
        integer, intent(in) :: needed
        real(dp), allocatable :: a(:), b(:), c(:), d(:)
        integer :: newcap

        if (needed <= seg%capacity) return

        newcap = max(needed, 2*seg%capacity)
        allocate(a(newcap), b(newcap), c(newcap), d(newcap))

        if (seg%n > 0) then
            a(1:seg%n) = seg%x1(1:seg%n)
            b(1:seg%n) = seg%y1(1:seg%n)
            c(1:seg%n) = seg%x2(1:seg%n)
            d(1:seg%n) = seg%y2(1:seg%n)
        end if

        call move_alloc(a, seg%x1)
        call move_alloc(b, seg%y1)
        call move_alloc(c, seg%x2)
        call move_alloc(d, seg%y2)
        seg%capacity = newcap
    end subroutine ensure_capacity


    subroutine add_segment(seg, xa, ya, xb, yb)
        type(segment_set), intent(inout) :: seg
        real(dp), intent(in) :: xa, ya, xb, yb

        if (.not. allocated(seg%x1)) call init_segments(seg, 128)
        call ensure_capacity(seg, seg%n + 1)

        seg%n = seg%n + 1
        seg%x1(seg%n) = xa
        seg%y1(seg%n) = ya
        seg%x2(seg%n) = xb
        seg%y2(seg%n) = yb
    end subroutine add_segment


    subroutine build_level_segments(field, level, seg)
        type(grid_field), intent(in) :: field
        real(dp), intent(in) :: level
        type(segment_set), intent(out) :: seg

        integer :: i, j

        call init_segments(seg, 2*(field%nx-1)*(field%ny-1))

        do i = 1, field%nx-1
            do j = 1, field%ny-1
                call process_cell(field%x(i), field%x(i+1), &
                                  field%y(j), field%y(j+1), &
                                  field%psi(i,j), field%psi(i+1,j), &
                                  field%psi(i+1,j+1), field%psi(i,j+1), &
                                  level, seg)
            end do
        end do
    end subroutine build_level_segments


    subroutine process_cell(xa, xb, ya, yb, p00, p10, p11, p01, level, seg)
        real(dp), intent(in) :: xa, xb, ya, yb
        real(dp), intent(in) :: p00, p10, p11, p01, level
        type(segment_set), intent(inout) :: seg

        real(dp) :: f(4), xp(4), yp(4), fc
        integer :: np

        f = [p00-level, p10-level, p11-level, p01-level]
        np = 0

        call collect_edge(xa,ya,f(1), xb,ya,f(2), xp,yp,np,seg)
        call collect_edge(xb,ya,f(2), xb,yb,f(3), xp,yp,np,seg)
        call collect_edge(xb,yb,f(3), xa,yb,f(4), xp,yp,np,seg)
        call collect_edge(xa,yb,f(4), xa,ya,f(1), xp,yp,np,seg)

        select case (np)
        case (2)
            call add_segment(seg, xp(1),yp(1), xp(2),yp(2))
        case (4)
            fc = bilinear_decider(f)
            if (fc*f(1) >= 0.0_dp) then
                call add_segment(seg, xp(1),yp(1), xp(2),yp(2))
                call add_segment(seg, xp(3),yp(3), xp(4),yp(4))
            else
                call add_segment(seg, xp(1),yp(1), xp(4),yp(4))
                call add_segment(seg, xp(2),yp(2), xp(3),yp(3))
            end if
        case default
            ! np=0 is the normal empty-cell case.  np=1 or 3 can only
            ! occur in a degenerate cell with a contour lying exactly on
            ! an edge; that edge has already been emitted by collect_edge.
        end select
    end subroutine process_cell


    subroutine collect_edge(xa,ya,fa,xb,yb,fb,xp,yp,np,seg)
        real(dp), intent(in) :: xa,ya,fa,xb,yb,fb
        real(dp), intent(inout) :: xp(4),yp(4)
        integer, intent(inout) :: np
        type(segment_set), intent(inout) :: seg
        real(dp) :: xi,yi

        if (abs(fa)<=tiny_value .and. abs(fb)<=tiny_value) then
            call add_segment(seg,xa,ya,xb,yb)
            return
        end if

        if (abs(fa)<=tiny_value) then
            call add_unique_point(xp,yp,np,xa,ya)
        else if (abs(fb)<=tiny_value) then
            call add_unique_point(xp,yp,np,xb,yb)
        else if (fa*fb<0.0_dp) then
            call interpolate_edge(xa,ya,fa,xb,yb,fb,xi,yi)
            call add_unique_point(xp,yp,np,xi,yi)
        end if
    end subroutine collect_edge


    subroutine add_unique_point(xp,yp,np,x,y)
        real(dp), intent(inout) :: xp(4),yp(4)
        integer, intent(inout) :: np
        real(dp), intent(in) :: x,y
        integer :: k
        real(dp), parameter :: point_tol=1.0e-12_dp

        do k=1,np
            if ((xp(k)-x)**2+(yp(k)-y)**2 <= point_tol*point_tol) return
        end do

        if (np<4) then
            np=np+1
            xp(np)=x;yp(np)=y
        end if
    end subroutine add_unique_point


    logical function crosses(a,b)
        real(dp), intent(in) :: a,b

        if (abs(a) <= tiny_value .and. abs(b) <= tiny_value) then
            crosses = .false.
        else if (abs(a) <= tiny_value .or. abs(b) <= tiny_value) then
            crosses = .true.
        else
            crosses = (a*b < 0.0_dp)
        end if
    end function crosses


    subroutine interpolate_edge(xa,ya,fa, xb,yb,fb, xp,yp)
        real(dp), intent(in) :: xa,ya,fa,xb,yb,fb
        real(dp), intent(out) :: xp,yp
        real(dp) :: t

        if (abs(fa-fb) <= tiny_value) then
            t = 0.5_dp
        else
            t = fa/(fa-fb)
        end if

        t = max(0.0_dp, min(1.0_dp, t))
        xp = xa + t*(xb-xa)
        yp = ya + t*(yb-ya)
    end subroutine interpolate_edge


    real(dp) function bilinear_decider(f) result(value)
        real(dp), intent(in) :: f(4)
        real(dp) :: a,b,c,d,u,v

        ! f(u,v) = a + b*u + c*v + d*u*v
        a = f(1)
        b = f(2)-f(1)
        c = f(4)-f(1)
        d = f(3)-f(2)-f(4)+f(1)

        if (abs(d) > tiny_value) then
            u = -c/d
            v = -b/d
            if (u > 0.0_dp .and. u < 1.0_dp .and. &
                v > 0.0_dp .and. v < 1.0_dp) then
                value = a + b*u + c*v + d*u*v
                return
            end if
        end if

        value = 0.25_dp*sum(f)
    end function bilinear_decider


    subroutine clip_segments_outside_circle(seg, xc, yc, radius)
        type(segment_set), intent(inout) :: seg
        real(dp), intent(in) :: xc,yc,radius

        type(segment_set) :: out
        integer :: s, npiece
        real(dp) :: px(4), py(4)

        call init_segments(out, max(16,seg%n+8))

        do s = 1, seg%n
            call clip_one(seg%x1(s),seg%y1(s),seg%x2(s),seg%y2(s), &
                          xc,yc,radius,npiece,px,py)
            if (npiece >= 1) call add_segment(out,px(1),py(1),px(2),py(2))
            if (npiece >= 2) call add_segment(out,px(3),py(3),px(4),py(4))
        end do

        call move_alloc(out%x1, seg%x1)
        call move_alloc(out%y1, seg%y1)
        call move_alloc(out%x2, seg%x2)
        call move_alloc(out%y2, seg%y2)
        seg%n = out%n
        seg%capacity = out%capacity
    end subroutine clip_segments_outside_circle


    subroutine clip_one(x1,y1,x2,y2,xc,yc,radius,npiece,px,py)
        real(dp), intent(in) :: x1,y1,x2,y2,xc,yc,radius
        integer, intent(out) :: npiece
        real(dp), intent(out) :: px(4),py(4)

        real(dp) :: dx,dy,fx,fy,aa,bb,cc,disc,r1,r2
        real(dp) :: t(4),ta,tb,tm,q,tmp
        integer :: nt,i,j

        npiece = 0
        px = 0.0_dp
        py = 0.0_dp

        dx = x2-x1
        dy = y2-y1
        aa = dx*dx + dy*dy
        if (aa <= tiny_value) return

        fx = x1-xc
        fy = y1-yc
        bb = 2.0_dp*(fx*dx + fy*dy)
        cc = fx*fx + fy*fy - radius*radius
        disc = bb*bb - 4.0_dp*aa*cc

        nt = 2
        t(1)=0.0_dp
        t(2)=1.0_dp

        if (disc > tiny_value) then
            r1=(-bb-sqrt(disc))/(2.0_dp*aa)
            r2=(-bb+sqrt(disc))/(2.0_dp*aa)
            if (r1>0.0_dp .and. r1<1.0_dp) then
                nt=nt+1; t(nt)=r1
            end if
            if (r2>0.0_dp .and. r2<1.0_dp) then
                nt=nt+1; t(nt)=r2
            end if
        end if

        do i=1,nt-1
            do j=i+1,nt
                if (t(j)<t(i)) then
                    tmp=t(i); t(i)=t(j); t(j)=tmp
                end if
            end do
        end do

        do i=1,nt-1
            ta=t(i); tb=t(i+1)
            if (tb-ta <= 1.0e-13_dp) cycle
            tm=0.5_dp*(ta+tb)
            q=(x1+tm*dx-xc)**2 + (y1+tm*dy-yc)**2 - radius*radius
            if (q >= 0.0_dp) then
                npiece=npiece+1
                if (npiece==1) then
                    px(1)=x1+ta*dx; py(1)=y1+ta*dy
                    px(2)=x1+tb*dx; py(2)=y1+tb*dy
                else if (npiece==2) then
                    px(3)=x1+ta*dx; py(3)=y1+ta*dy
                    px(4)=x1+tb*dx; py(4)=y1+tb*dy
                end if
            end if
        end do
    end subroutine clip_one

end module contour_marching

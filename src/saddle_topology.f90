module saddle_topology
    use kinds, only: dp
    use separatrix_types, only: saddle_point
    use contour_marching, only: segment_set, clip_segments_outside_circle, add_segment
    implicit none
    private

    public :: reconstruct_level_saddles

contains

    subroutine reconstruct_level_saddles(seg, saddles, radius, connect_tol, min_cos)
        type(segment_set), intent(inout) :: seg
        type(saddle_point), intent(in) :: saddles(:)
        real(dp), intent(in) :: radius, connect_tol, min_cos

        integer :: k
        real(dp) :: local_radius

        do k=1,size(saddles)
            local_radius = safe_radius(saddles, k, radius)
            call clip_segments_outside_circle(seg, saddles(k)%x, saddles(k)%y, local_radius)
        end do

        do k=1,size(saddles)
            local_radius = safe_radius(saddles, k, radius)
            call add_saddle_spokes(seg, saddles(k), local_radius, connect_tol, min_cos)
        end do
    end subroutine reconstruct_level_saddles


    real(dp) function safe_radius(saddles, k, requested) result(r)
        type(saddle_point), intent(in) :: saddles(:)
        integer, intent(in) :: k
        real(dp), intent(in) :: requested

        integer :: j
        real(dp) :: d, nearest

        nearest = huge(1.0_dp)
        do j=1,size(saddles)
            if (j==k) cycle
            d = hypot2(saddles(k)%x-saddles(j)%x, saddles(k)%y-saddles(j)%y)
            nearest = min(nearest,d)
        end do

        r = requested
        if (nearest < huge(1.0_dp)/2.0_dp) r = min(r,0.25_dp*nearest)
    end function safe_radius


    subroutine add_saddle_spokes(seg, saddle, radius, tol, min_cos)
        type(segment_set), intent(inout) :: seg
        type(saddle_point), intent(in) :: saddle
        real(dp), intent(in) :: radius,tol,min_cos

        integer, parameter :: maxcand=64
        real(dp) :: cx(maxcand),cy(maxcand),dx(4),dy(4)
        logical :: used(maxcand)
        integer :: nc,s,e,k,best
        real(dp) :: r,vx,vy,dot,best_dot

        call saddle_directions(saddle%h11,saddle%h12,saddle%h22,dx,dy)

        nc=0
        do s=1,seg%n
            do e=1,2
                if (e==1) then
                    vx=seg%x1(s); vy=seg%y1(s)
                else
                    vx=seg%x2(s); vy=seg%y2(s)
                end if

                if (endpoint_degree(seg,vx,vy,tol) /= 1) cycle
                r=hypot2(vx-saddle%x,vy-saddle%y)
                if (abs(r-radius) <= max(10.0_dp*tol,1.0e-5_dp*max(1.0_dp,radius))) then
                    if (.not. candidate_exists(cx,cy,nc,vx,vy,10.0_dp*tol)) then
                        if (nc<maxcand) then
                            nc=nc+1; cx(nc)=vx; cy(nc)=vy
                        end if
                    end if
                end if
            end do
        end do

        if (nc < 4) then
            write(*,'(A,I0,A,I0)') 'WARNING: saddle ',saddle%id, &
                ' has fewer than four circle intersections: ',nc
        end if

        used=.false.
        do k=1,4
            best=0; best_dot=-huge(1.0_dp)
            do s=1,nc
                if (used(s)) cycle
                vx=cx(s)-saddle%x; vy=cy(s)-saddle%y
                r=hypot2(vx,vy)
                if (r<=tiny(1.0_dp)) cycle
                vx=vx/r; vy=vy/r
                dot=vx*dx(k)+vy*dy(k)
                if (dot>best_dot) then
                    best_dot=dot; best=s
                end if
            end do

            if (best>0 .and. best_dot>=min_cos) then
                used(best)=.true.
                call add_segment(seg,saddle%x,saddle%y,cx(best),cy(best))
            else
                write(*,'(A,I0,A,I0,A,F8.4)') 'WARNING: saddle ',saddle%id, &
                    ', direction ',k,' unmatched; best cosine=',best_dot
            end if
        end do
    end subroutine add_saddle_spokes


    subroutine saddle_directions(a,b,c,dx,dy)
        real(dp), intent(in) :: a,b,c
        real(dp), intent(out) :: dx(4),dy(4)
        real(dp) :: disc,r1,r2,norm
        integer :: k

        disc=b*b-a*c
        if (disc<=0.0_dp) error stop 'Hessian is not hyperbolic'

        if (abs(c)>=abs(a)) then
            r1=(-b+sqrt(disc))/c
            r2=(-b-sqrt(disc))/c
            dx=[1.0_dp,-1.0_dp,1.0_dp,-1.0_dp]
            dy=[r1,-r1,r2,-r2]
        else
            r1=(-b+sqrt(disc))/a
            r2=(-b-sqrt(disc))/a
            dx=[r1,-r1,r2,-r2]
            dy=[1.0_dp,-1.0_dp,1.0_dp,-1.0_dp]
        end if

        do k=1,4
            norm=hypot2(dx(k),dy(k))
            dx(k)=dx(k)/norm
            dy(k)=dy(k)/norm
        end do
    end subroutine saddle_directions


    integer function endpoint_degree(seg,x,y,tol) result(deg)
        type(segment_set), intent(in) :: seg
        real(dp), intent(in) :: x,y,tol
        integer :: s

        deg=0
        do s=1,seg%n
            if (hypot2(x-seg%x1(s),y-seg%y1(s))<=tol) deg=deg+1
            if (hypot2(x-seg%x2(s),y-seg%y2(s))<=tol) deg=deg+1
        end do
    end function endpoint_degree


    logical function candidate_exists(cx,cy,n,x,y,tol) result(found)
        real(dp), intent(in) :: cx(:),cy(:),x,y,tol
        integer, intent(in) :: n
        integer :: i
        found=.false.
        do i=1,n
            if (hypot2(cx(i)-x,cy(i)-y)<=tol) then
                found=.true.; return
            end if
        end do
    end function candidate_exists


    pure real(dp) function hypot2(x,y) result(r)
        real(dp), intent(in) :: x,y
        r=sqrt(x*x+y*y)
    end function hypot2

end module saddle_topology

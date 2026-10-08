module saddle_refinement
    use kinds, only: dp
    use separatrix_types, only: grid_field, saddle_point
    implicit none
    private

    public :: refine_saddles_quadratic

contains

    subroutine refine_saddles_quadratic(field, saddles)
        type(grid_field), intent(in) :: field
        type(saddle_point), allocatable, intent(inout) :: saddles(:)

        integer :: k, i0, j0

        do k = 1, size(saddles)
            call nearest_indices(field, saddles(k)%x, saddles(k)%y, i0, j0)
            call refine_one(field, i0, j0, saddles(k))
        end do

        call keep_hyperbolic(saddles)
        call deduplicate_close(field, saddles)
    end subroutine refine_saddles_quadratic


    subroutine nearest_indices(field, xs, ys, i0, j0)
        type(grid_field), intent(in) :: field
        real(dp), intent(in) :: xs, ys
        integer, intent(out) :: i0, j0

        integer :: i, j
        real(dp) :: best, d

        best = huge(1.0_dp)
        i0 = 1
        do i = 1, field%nx
            d = abs(field%x(i)-xs)
            if (d < best) then
                best = d
                i0 = i
            end if
        end do

        best = huge(1.0_dp)
        j0 = 1
        do j = 1, field%ny
            d = abs(field%y(j)-ys)
            if (d < best) then
                best = d
                j0 = j
            end if
        end do

        i0 = max(3, min(field%nx-2, i0))
        j0 = max(3, min(field%ny-2, j0))
    end subroutine nearest_indices


    subroutine refine_one(field, i0, j0, s)
        type(grid_field), intent(in) :: field
        integer, intent(in) :: i0, j0
        type(saddle_point), intent(inout) :: s

        integer, parameter :: ncoef = 6, npts = 25
        real(dp) :: A(npts,ncoef), b(npts), ATA(ncoef,ncoef), ATb(ncoef), c(ncoef)
        real(dp) :: X, Y, dxs, dys, det_h
        integer :: i, j, k, p, q, info

        k = 0
        do i = i0-2, i0+2
            do j = j0-2, j0+2
                k = k + 1
                X = field%x(i) - field%x(i0)
                Y = field%y(j) - field%y(j0)

                A(k,:) = [1.0_dp, X, Y, X*X, X*Y, Y*Y]
                b(k) = field%psi(i,j)
            end do
        end do

        ATA = 0.0_dp
        ATb = 0.0_dp

        do p = 1, ncoef
            ATb(p) = sum(A(:,p)*b(:))
            do q = 1, ncoef
                ATA(p,q) = sum(A(:,p)*A(:,q))
            end do
        end do

        call solve6(ATA, ATb, c, info)
        if (info /= 0) return

        s%h11 = 2.0_dp*c(4)
        s%h12 = c(5)
        s%h22 = 2.0_dp*c(6)
        det_h = s%h11*s%h22 - s%h12*s%h12
        s%det_h = det_h

        if (det_h >= 0.0_dp) return

        dxs = (-c(2)*s%h22 + s%h12*c(3)) / det_h
        dys = (-s%h11*c(3) + s%h12*c(2)) / det_h

        s%x = field%x(i0) + dxs
        s%y = field%y(j0) + dys

        s%psi = c(1) + c(2)*dxs + c(3)*dys + c(4)*dxs*dxs + &
                c(5)*dxs*dys + c(6)*dys*dys
    end subroutine refine_one


    subroutine solve6(A, b, x, info)
        real(dp), intent(inout) :: A(6,6)
        real(dp), intent(inout) :: b(6)
        real(dp), intent(out) :: x(6)
        integer, intent(out) :: info

        integer :: i, j, k, p
        real(dp) :: vmax, factor, tmp

        info = 0

        do k = 1,5
            p = k
            vmax = abs(A(k,k))

            do i = k+1,6
                if (abs(A(i,k)) > vmax) then
                    vmax = abs(A(i,k))
                    p = i
                end if
            end do

            if (vmax < 1.0e-20_dp) then
                info = 1
                return
            end if

            if (p /= k) then
                do j = 1,6
                    tmp = A(k,j)
                    A(k,j) = A(p,j)
                    A(p,j) = tmp
                end do
                tmp = b(k)
                b(k) = b(p)
                b(p) = tmp
            end if

            do i = k+1,6
                factor = A(i,k)/A(k,k)
                do j = k+1,6
                    A(i,j) = A(i,j) - factor*A(k,j)
                end do
                b(i) = b(i) - factor*b(k)
            end do
        end do

        if (abs(A(6,6)) < 1.0e-20_dp) then
            info = 1
            return
        end if

        x(6) = b(6)/A(6,6)

        do i = 5,1,-1
            tmp = b(i)
            do j = i+1,6
                tmp = tmp - A(i,j)*x(j)
            end do
            x(i) = tmp/A(i,i)
        end do
    end subroutine solve6


    subroutine keep_hyperbolic(saddles)
        type(saddle_point), allocatable, intent(inout) :: saddles(:)
        type(saddle_point), allocatable :: temp(:)
        logical, allocatable :: keep(:)

        allocate(keep(size(saddles)))
        keep = saddles%det_h < 0.0_dp
        allocate(temp(count(keep)))
        temp = pack(saddles, keep)
        call move_alloc(temp, saddles)
        deallocate(keep)
    end subroutine keep_hyperbolic


    subroutine deduplicate_close(field, saddles)
        type(grid_field), intent(in) :: field
        type(saddle_point), allocatable, intent(inout) :: saddles(:)
        type(saddle_point), allocatable :: temp(:)
        logical, allocatable :: keep(:)
        integer :: i,j,n
        real(dp) :: dx,dy,tol,d

        n=size(saddles)
        if (n<=1) return

        dx=abs(field%x(2)-field%x(1))
        dy=abs(field%y(2)-field%y(1))
        tol=sqrt(dx*dx+dy*dy)

        allocate(keep(n)); keep=.true.
        do i=1,n
            if (.not.keep(i)) cycle
            do j=i+1,n
                if (.not.keep(j)) cycle
                d=sqrt((saddles(i)%x-saddles(j)%x)**2 + &
                       (saddles(i)%y-saddles(j)%y)**2)
                if (d<=tol) keep(j)=.false.
            end do
        end do

        allocate(temp(count(keep)))
        temp=pack(saddles,keep)
        call move_alloc(temp,saddles)
        deallocate(keep)
    end subroutine deduplicate_close

end module saddle_refinement

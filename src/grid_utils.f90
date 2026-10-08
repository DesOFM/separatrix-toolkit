module grid_utils
    use kinds, only: dp
    use separatrix_types, only: grid_field
    implicit none
    private

    public :: validate_rectilinear_grid
    public :: characteristic_spacing
    public :: bilinear_interpolate

contains

    subroutine validate_rectilinear_grid(field)
        type(grid_field), intent(in) :: field
        integer :: i

        if (field%nx < 2 .or. field%ny < 2) error stop 'Grid needs at least 2x2 points'
        if (.not. allocated(field%x) .or. .not. allocated(field%y) .or. .not. allocated(field%psi)) &
            error stop 'Grid field is not allocated'
        if (size(field%psi,1) /= field%nx .or. size(field%psi,2) /= field%ny) &
            error stop 'psi dimensions do not match nx, ny'

        do i=2,field%nx
            if (field%x(i) <= field%x(i-1)) error stop 'x coordinates must be strictly increasing'
        end do
        do i=2,field%ny
            if (field%y(i) <= field%y(i-1)) error stop 'y coordinates must be strictly increasing'
        end do
    end subroutine validate_rectilinear_grid


    real(dp) function characteristic_spacing(field) result(h)
        type(grid_field), intent(in) :: field
        real(dp), allocatable :: d(:)
        integer :: n, i

        n=(field%nx-1)+(field%ny-1)
        allocate(d(n))
        do i=1,field%nx-1
            d(i)=field%x(i+1)-field%x(i)
        end do
        do i=1,field%ny-1
            d(field%nx-1+i)=field%y(i+1)-field%y(i)
        end do
        call sort_real(d)
        if (mod(n,2)==1) then
            h=d((n+1)/2)
        else
            h=0.5_dp*(d(n/2)+d(n/2+1))
        end if
        deallocate(d)
    end function characteristic_spacing


    real(dp) function bilinear_interpolate(field,xp,yp,inside) result(v)
        type(grid_field), intent(in) :: field
        real(dp), intent(in) :: xp,yp
        logical, intent(out), optional :: inside
        integer :: i,j
        real(dp) :: u,w
        logical :: ok

        ok = xp >= field%x(1) .and. xp <= field%x(field%nx) .and. &
             yp >= field%y(1) .and. yp <= field%y(field%ny)
        if (present(inside)) inside=ok
        if (.not.ok) then
            v=huge(1.0_dp)
            return
        end if

        i=locate_interval(field%x,xp)
        j=locate_interval(field%y,yp)
        u=(xp-field%x(i))/(field%x(i+1)-field%x(i))
        w=(yp-field%y(j))/(field%y(j+1)-field%y(j))

        v=(1.0_dp-u)*(1.0_dp-w)*field%psi(i,j) + &
          u*(1.0_dp-w)*field%psi(i+1,j) + &
          u*w*field%psi(i+1,j+1) + &
          (1.0_dp-u)*w*field%psi(i,j+1)
    end function bilinear_interpolate


    integer function locate_interval(a,x) result(i)
        real(dp), intent(in) :: a(:),x
        integer :: lo,hi,mid,n
        n=size(a)
        if (x<=a(1)) then
            i=1;return
        else if (x>=a(n)) then
            i=n-1;return
        end if
        lo=1;hi=n
        do while(hi-lo>1)
            mid=(lo+hi)/2
            if (a(mid)<=x) then
                lo=mid
            else
                hi=mid
            end if
        end do
        i=lo
    end function locate_interval


    subroutine sort_real(a)
        real(dp), intent(inout) :: a(:)
        integer :: i,j
        real(dp) :: key
        do i=2,size(a)
            key=a(i);j=i-1
            do while(j>=1)
                if (a(j)<=key) exit
                a(j+1)=a(j);j=j-1
            end do
            a(j+1)=key
        end do
    end subroutine sort_real

end module grid_utils

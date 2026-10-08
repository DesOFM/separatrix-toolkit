module field_provider
    use kinds, only: dp
    use separatrix_types, only: grid_field
    implicit none
    private

    public :: scalar_field
    public :: sample_function_on_grid

    abstract interface
        function scalar_field(x, y) result(value)
            import dp
            real(dp), intent(in) :: x
            real(dp), intent(in) :: y
            real(dp) :: value
        end function scalar_field
    end interface

contains

    subroutine sample_function_on_grid(field, f)
        type(grid_field), intent(inout) :: field
        procedure(scalar_field) :: f

        integer :: i, j

        if (.not. allocated(field%x)) error stop "field%x is not allocated"
        if (.not. allocated(field%y)) error stop "field%y is not allocated"

        field%nx = size(field%x)
        field%ny = size(field%y)

        if (allocated(field%psi)) deallocate(field%psi)
        allocate(field%psi(field%nx, field%ny))

        do i = 1, field%nx
            do j = 1, field%ny
                field%psi(i,j) = f(field%x(i), field%y(j))
            end do
        end do
    end subroutine sample_function_on_grid

end module field_provider

module grid_io
    use kinds, only: dp
    use separatrix_types, only: grid_field
    implicit none
    private

    public :: read_xyz_grid
    public :: write_xyz_grid

contains

    subroutine read_xyz_grid(filename, field)
        character(len=*), intent(in) :: filename
        type(grid_field), intent(out) :: field

        integer :: unit, ios, nrows, k, nx, ny, i, j
        real(dp), allocatable :: tx(:), ty(:), tp(:)
        real(dp) :: x0
        real(dp), parameter :: tol = 1.0e-12_dp

        nrows = 0
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop "Cannot open grid file"

        do
            read(unit, *, iostat=ios)
            if (ios /= 0) exit
            nrows = nrows + 1
        end do
        close(unit)

        if (nrows < 4) error stop "Grid file is too small"

        allocate(tx(nrows), ty(nrows), tp(nrows))

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) error stop "Cannot reopen grid file"

        do k = 1, nrows
            read(unit, *, iostat=ios) tx(k), ty(k), tp(k)
            if (ios /= 0) error stop "Error while reading x y psi"
        end do
        close(unit)

        x0 = tx(1)
        ny = 1
        do k = 2, nrows
            if (abs(tx(k) - x0) <= tol) then
                ny = ny + 1
            else
                exit
            end if
        end do

        if (mod(nrows, ny) /= 0) error stop "Input is not a rectangular x-major grid"
        nx = nrows / ny

        field%nx = nx
        field%ny = ny
        allocate(field%x(nx), field%y(ny), field%psi(nx,ny))

        k = 0
        do i = 1, nx
            do j = 1, ny
                k = k + 1
                field%x(i) = tx(k)
                if (i == 1) field%y(j) = ty(k)
                field%psi(i,j) = tp(k)
            end do
        end do

        deallocate(tx, ty, tp)
    end subroutine read_xyz_grid


    subroutine write_xyz_grid(filename, field)
        character(len=*), intent(in) :: filename
        type(grid_field), intent(in) :: field

        integer :: unit, i, j

        open(newunit=unit, file=filename, status='replace', action='write')

        do i = 1, field%nx
            do j = 1, field%ny
                write(unit,'(3ES24.15)') field%x(i), field%y(j), field%psi(i,j)
            end do
        end do

        close(unit)
    end subroutine write_xyz_grid

end module grid_io

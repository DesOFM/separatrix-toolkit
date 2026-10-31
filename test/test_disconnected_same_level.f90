program test_disconnected_same_level
    use kinds, only: dp
    use separatrix_types, only: grid_field, saddle_point, separatrix_network
    use field_provider, only: sample_function_on_grid
    use separatrix_builder, only: build_all_separatrices
    implicit none

    type(grid_field) :: f
    type(saddle_point), allocatable :: s(:)
    type(separatrix_network), allocatable :: n(:)

    integer :: i, k
    real(dp) :: xmax

    allocate(f%x(241), f%y(161))

    do i = 1, 241
        f%x(i) = -2.0_dp + 0.025_dp*real(i-1,dp)
    end do

    do i = 1, 161
        f%y(i) = -2.0_dp + 0.025_dp*real(i-1,dp)
    end do

    call sample_function_on_grid(f, psi)

    allocate(s(1))
    s(1)%id = 1
    s(1)%x = 0.0_dp
    s(1)%y = 0.0_dp
    s(1)%psi = 0.0_dp
    s(1)%h11 = 17.5_dp
    s(1)%h12 = 0.0_dp
    s(1)%h22 = -17.5_dp
    s(1)%det_h = s(1)%h11*s(1)%h22

    call build_all_separatrices(f, s, n, level_tol=1.0e-8_dp)

    if (size(n) /= 1) error stop 'disconnected-level test: expected one saddle-connected network'
    if (.not.allocated(n(1)%branches)) error stop 'disconnected-level test: no branches reconstructed'
    if (n(1)%nbranches <= 0) error stop 'disconnected-level test: expected at least one branch'

    xmax = -huge(1.0_dp)
    do k = 1, size(n(1)%branches)
        if (allocated(n(1)%branches(k)%x)) xmax = max(xmax, maxval(n(1)%branches(k)%x))
    end do

    ! The zero level also contains a disconnected circle centered at x=3 with radius 0.5.
    ! Saddle-connected branches reach the y-boundary near |x|=2 and must not absorb that circle.
    if (xmax > 2.2_dp) error stop 'disconnected-level test: unrelated contour component entered saddle network'

    print *, 'disconnected-level test: PASS'

contains

    function psi(x,y) result(v)
        real(dp), intent(in) :: x,y
        real(dp) :: v
        real(dp) :: saddle_factor, circle_factor

        saddle_factor = x*x - y*y
        circle_factor = (x-3.0_dp)**2 + y*y - 0.25_dp
        v = saddle_factor*circle_factor
    end function psi

end program test_disconnected_same_level

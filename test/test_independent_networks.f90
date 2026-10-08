program test_independent_networks
    use kinds, only: dp
    use separatrix_types, only: grid_field, saddle_point, separatrix_network
    use field_provider, only: sample_function_on_grid
    use critical_points, only: find_all_saddles
    use saddle_refinement, only: refine_saddles_quadratic
    use separatrix_builder, only: build_all_separatrices
    implicit none

    type(grid_field) :: f
    type(saddle_point), allocatable :: s(:)
    type(separatrix_network), allocatable :: n(:)

    integer :: i, ns

    allocate(f%x(241), f%y(161))

    do i = 1, 241
        f%x(i) = -3.0_dp + 0.025_dp*real(i-1,dp)
    end do

    do i = 1, 161
        f%y(i) = -2.0_dp + 0.025_dp*real(i-1,dp)
    end do

    call sample_function_on_grid(f, psi)
    call find_all_saddles(f, s, ns)
    call refine_saddles_quadratic(f, s)

    if (size(s) /= 2) error stop 'independent-networks test: expected two saddles'

    s(1)%id = 1
    s(2)%id = 2

    call build_all_separatrices(f, s, n, level_tol=1.0e-6_dp)

    if (size(n) /= 2) error stop 'independent-networks test: expected two networks'

    print *, 'independent-networks test: PASS'

contains

    function psi(x,y) result(v)
        real(dp), intent(in) :: x,y
        real(dp) :: v

        ! The linear tilt separates the two saddle energy levels.
        v = (x*x - 1.0_dp)**2 + 0.2_dp*x - y*y
    end function psi

end program test_independent_networks

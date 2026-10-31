program test_nonaxisymmetric_full_pipeline
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
    real(dp), parameter :: xs = -13.0_dp/54.0_dp
    real(dp), parameter :: ys = -4.0_dp/27.0_dp

    allocate(f%x(161), f%y(161))

    do i = 1, 161
        f%x(i) = -2.0_dp + 0.025_dp*real(i-1,dp)
        f%y(i) = -2.0_dp + 0.025_dp*real(i-1,dp)
    end do

    call sample_function_on_grid(f, psi)
    call find_all_saddles(f, s, ns)
    call refine_saddles_quadratic(f, s)

    if (size(s) /= 1) error stop 'nonaxisymmetric pipeline test: expected one saddle'
    if (abs(s(1)%x-xs) > 1.0e-6_dp) error stop 'nonaxisymmetric pipeline test: saddle x mismatch'
    if (abs(s(1)%y-ys) > 1.0e-6_dp) error stop 'nonaxisymmetric pipeline test: saddle y mismatch'
    if (s(1)%det_h >= 0.0_dp) error stop 'nonaxisymmetric pipeline test: point is not hyperbolic'
    if (abs(s(1)%h12) < 0.5_dp) error stop 'nonaxisymmetric pipeline test: mixed Hessian term was lost'

    s(1)%id = 1
    call build_all_separatrices(f, s, n)

    if (size(n) /= 1) error stop 'nonaxisymmetric pipeline test: expected one network'
    if (n(1)%nbranches /= 4) error stop 'nonaxisymmetric pipeline test: expected four boundary branches'

    print *, 'nonaxisymmetric full-pipeline test: PASS'

contains

    function psi(x,y) result(v)
        real(dp), intent(in) :: x,y
        real(dp) :: v

        ! Shifted, rotated indefinite quadratic. The xy term makes the saddle axes
        ! non-aligned with the computational grid, and the linear terms move it
        ! away from any symmetry center or grid node.
        v = x*x - 2.0_dp*y*y + 0.8_dp*x*y + 0.6_dp*x - 0.4_dp*y
    end function psi

end program test_nonaxisymmetric_full_pipeline

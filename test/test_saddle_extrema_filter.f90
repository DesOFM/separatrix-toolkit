program test_saddle_extrema_filter
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
    real(dp), parameter :: tol = 5.0e-3_dp

    allocate(f%x(241), f%y(161))

    do i = 1, 241
        f%x(i) = -3.0_dp + 0.025_dp*real(i-1,dp)
    end do

    do i = 1, 161
        f%y(i) = -2.0_dp + 0.025_dp*real(i-1,dp)
    end do

    call sample_function_on_grid(f, psi)
    call find_all_saddles(f, s, ns)

    ! The field has three critical points: minima at (+/-1,0)
    ! and one hyperbolic saddle at (0,0).  Candidate detection is
    ! intentionally topology-agnostic; refinement must remove extrema.
    if (ns < 3) error stop 'saddle-extrema filter test: expected at least three critical candidates'

    call refine_saddles_quadratic(f, s)

    if (size(s) /= 1) error stop 'saddle-extrema filter test: extrema were not removed'
    if (abs(s(1)%x) > tol .or. abs(s(1)%y) > tol) &
        error stop 'saddle-extrema filter test: wrong saddle position'
    if (s(1)%det_h >= 0.0_dp) error stop 'saddle-extrema filter test: retained point is not hyperbolic'

    s(1)%id = 1
    call build_all_separatrices(f, s, n, level_tol=1.0e-6_dp)

    if (size(n) /= 1) error stop 'saddle-extrema filter test: expected one separatrix network'
    if (n(1)%nbranches < 1) error stop 'saddle-extrema filter test: no branches reconstructed'

    print *, 'saddle-extrema filter test: PASS'

contains

    function psi(x,y) result(v)
        real(dp), intent(in) :: x,y
        real(dp) :: v

        ! Non-axisymmetric double-well field.
        ! Critical points:
        !   (0,0)    : Hessian diag(-4,2), saddle
        !   (+/-1,0) : Hessian diag(8,2), minima
        v = (x*x - 1.0_dp)**2 + y*y
    end function psi

end program test_saddle_extrema_filter

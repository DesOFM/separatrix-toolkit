program test_poincare_map
    use kinds, only: dp
    use poincare_map, only: map_jacobian_fd, refine_fixed_point, hyperbolic_eigendata
    implicit none

    real(dp) :: x(2), jac(2,2), lambda_s, lambda_u, vs(2), vu(2), residual
    integer :: status, iterations

    x = [0.0_dp, 0.0_dp]

    call refine_fixed_point(test_map, x, 1.0e-12_dp, 1.0e-6_dp, 12, status, residual, iterations)
    if (status /= 0) error stop 'fixed-point refinement failed'
    if (maxval(abs(x-[1.25_dp,-0.75_dp])) > 1.0e-10_dp) error stop 'wrong fixed point'
    if (residual > 1.0e-11_dp) error stop 'fixed-point residual too large'

    call map_jacobian_fd(test_map, x, 1.0e-6_dp, jac, status)
    if (status /= 0) error stop 'Jacobian evaluation failed'
    if (maxval(abs(jac-reshape([2.0_dp,0.0_dp,0.0_dp,0.5_dp],[2,2]))) > 1.0e-9_dp) &
        error stop 'wrong map Jacobian'

    call hyperbolic_eigendata(jac,lambda_s,lambda_u,vs,vu,status)
    if (status /= 0) error stop 'hyperbolic eigendata failed'
    if (abs(lambda_s-0.5_dp) > 1.0e-9_dp) error stop 'wrong stable multiplier'
    if (abs(lambda_u-2.0_dp) > 1.0e-9_dp) error stop 'wrong unstable multiplier'

    print *, 'test_poincare_map: PASS'

contains

    subroutine test_map(x0,x1,status)
        real(dp), intent(in) :: x0(2)
        real(dp), intent(out) :: x1(2)
        integer, intent(out) :: status
        real(dp), parameter :: xs(2)=[1.25_dp,-0.75_dp]

        x1(1)=xs(1)+2.0_dp*(x0(1)-xs(1))
        x1(2)=xs(2)+0.5_dp*(x0(2)-xs(2))
        status=0
    end subroutine test_map

end program test_poincare_map

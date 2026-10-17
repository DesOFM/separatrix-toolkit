module poincare_map
    use kinds, only: dp
    implicit none
    private

    public :: map2d
    public :: map_jacobian_fd
    public :: refine_fixed_point
    public :: hyperbolic_eigendata

    abstract interface
        subroutine map2d(x0, x1, status)
            import dp
            real(dp), intent(in) :: x0(2)
            real(dp), intent(out) :: x1(2)
            integer, intent(out) :: status
        end subroutine map2d
    end interface

contains

    subroutine map_jacobian_fd(map, x, delta, jac, status)
        procedure(map2d) :: map
        real(dp), intent(in) :: x(2)
        real(dp), intent(in) :: delta
        real(dp), intent(out) :: jac(2,2)
        integer, intent(out) :: status

        real(dp) :: h, xp(2), xm(2), fp(2), fm(2)
        integer :: j, istat

        if (delta <= 0.0_dp) then
            status = -1
            jac = 0.0_dp
            return
        end if

        do j = 1, 2
            h = delta * max(1.0_dp, abs(x(j)))
            xp = x
            xm = x
            xp(j) = xp(j) + h
            xm(j) = xm(j) - h

            call map(xp, fp, istat)
            if (istat /= 0) then
                status = istat
                jac = 0.0_dp
                return
            end if

            call map(xm, fm, istat)
            if (istat /= 0) then
                status = istat
                jac = 0.0_dp
                return
            end if

            jac(:,j) = (fp-fm)/(2.0_dp*h)
        end do

        status = 0
    end subroutine map_jacobian_fd


    subroutine refine_fixed_point(map, x, tol, fd_delta, max_iter, status, residual, iterations)
        procedure(map2d) :: map
        real(dp), intent(inout) :: x(2)
        real(dp), intent(in) :: tol, fd_delta
        integer, intent(in) :: max_iter
        integer, intent(out) :: status
        real(dp), intent(out), optional :: residual
        integer, intent(out), optional :: iterations

        real(dp) :: fx(2), r(2), jac(2,2), a(2,2), dx(2)
        real(dp) :: det_a, res
        integer :: k, istat

        if (tol <= 0.0_dp .or. fd_delta <= 0.0_dp .or. max_iter < 1) then
            status = -1
            if (present(residual)) residual = huge(1.0_dp)
            if (present(iterations)) iterations = 0
            return
        end if

        do k = 1, max_iter
            call map(x, fx, istat)
            if (istat /= 0) then
                status = istat
                if (present(iterations)) iterations = k-1
                return
            end if

            r = fx-x
            res = sqrt(sum(r*r))
            if (res <= tol) then
                status = 0
                if (present(residual)) residual = res
                if (present(iterations)) iterations = k-1
                return
            end if

            call map_jacobian_fd(map, x, fd_delta, jac, istat)
            if (istat /= 0) then
                status = istat
                if (present(iterations)) iterations = k-1
                return
            end if

            a = jac
            a(1,1) = a(1,1)-1.0_dp
            a(2,2) = a(2,2)-1.0_dp
            det_a = a(1,1)*a(2,2)-a(1,2)*a(2,1)

            if (abs(det_a) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a))**2)) then
                status = -2
                if (present(residual)) residual = res
                if (present(iterations)) iterations = k-1
                return
            end if

            ! Solve (DP-I) dx = -(P(x)-x).
            dx(1) = (-a(2,2)*r(1)+a(1,2)*r(2))/det_a
            dx(2) = ( a(2,1)*r(1)-a(1,1)*r(2))/det_a
            x = x+dx
        end do

        call map(x, fx, istat)
        if (istat == 0) then
            res = sqrt(sum((fx-x)**2))
        else
            res = huge(1.0_dp)
        end if

        status = 1
        if (present(residual)) residual = res
        if (present(iterations)) iterations = max_iter
    end subroutine refine_fixed_point


    subroutine hyperbolic_eigendata(jac, lambda_s, lambda_u, vec_s, vec_u, status)
        real(dp), intent(in) :: jac(2,2)
        real(dp), intent(out) :: lambda_s, lambda_u
        real(dp), intent(out) :: vec_s(2), vec_u(2)
        integer, intent(out) :: status

        real(dp) :: tr, detj, disc, l1, l2

        tr = jac(1,1)+jac(2,2)
        detj = jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
        disc = tr*tr-4.0_dp*detj

        if (disc < 0.0_dp) then
            status = -1
            lambda_s = 0.0_dp
            lambda_u = 0.0_dp
            vec_s = 0.0_dp
            vec_u = 0.0_dp
            return
        end if

        disc = sqrt(max(0.0_dp,disc))
        l1 = 0.5_dp*(tr+disc)
        l2 = 0.5_dp*(tr-disc)

        if (abs(l1) < abs(l2)) then
            lambda_s = l1
            lambda_u = l2
        else
            lambda_s = l2
            lambda_u = l1
        end if

        if (.not.(abs(lambda_s) < 1.0_dp .and. abs(lambda_u) > 1.0_dp)) then
            status = -2
            call eigenvector_2x2(jac,lambda_s,vec_s)
            call eigenvector_2x2(jac,lambda_u,vec_u)
            return
        end if

        call eigenvector_2x2(jac,lambda_s,vec_s)
        call eigenvector_2x2(jac,lambda_u,vec_u)
        status = 0
    end subroutine hyperbolic_eigendata


    subroutine eigenvector_2x2(a, lambda, v)
        real(dp), intent(in) :: a(2,2), lambda
        real(dp), intent(out) :: v(2)
        real(dp) :: nrm

        if (abs(a(1,2)) >= abs(a(2,1))) then
            v = [a(1,2), lambda-a(1,1)]
        else
            v = [lambda-a(2,2), a(2,1)]
        end if

        nrm = sqrt(sum(v*v))
        if (nrm <= tiny(1.0_dp)) then
            if (abs(a(1,1)-lambda) <= abs(a(2,2)-lambda)) then
                v = [1.0_dp,0.0_dp]
            else
                v = [0.0_dp,1.0_dp]
            end if
        else
            v = v/nrm
        end if
    end subroutine eigenvector_2x2

end module poincare_map

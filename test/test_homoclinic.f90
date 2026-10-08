program test_homoclinic
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

    integer :: i, ns, j, nhomoclinic

    allocate(f%x(201), f%y(161))

    do i = 1, 201
        f%x(i) = -2.5_dp + 0.025_dp*real(i-1,dp)
    end do

    do i = 1, 161
        f%y(i) = -2.0_dp + 0.025_dp*real(i-1,dp)
    end do

    call sample_function_on_grid(f, psi)
    call find_all_saddles(f, s, ns)
    call refine_saddles_quadratic(f, s)

    if (size(s) /= 1) error stop 'homoclinic test: expected one saddle'

    s(1)%id = 1
    call build_all_separatrices(f, s, n, level_tol=1.0e-8_dp)

    if (size(n) /= 1) error stop 'homoclinic test: expected one network'

    nhomoclinic = 0
    do j = 1, n(1)%nbranches
        if (n(1)%branches(j)%saddle_start > 0 .and. &
            n(1)%branches(j)%saddle_start == n(1)%branches(j)%saddle_end) then
            nhomoclinic = nhomoclinic + 1
        end if
    end do

    if (nhomoclinic /= 2) error stop 'homoclinic test: expected two homoclinic loops'

    print *, 'homoclinic test: PASS'

contains

    function psi(x,y) result(v)
        real(dp), intent(in) :: x,y
        real(dp) :: v

        ! Duffing Hamiltonian. H=0 contains two homoclinic loops.
        v = 0.5_dp*y*y - 0.5_dp*x*x + 0.25_dp*x**4
    end function psi

end program test_homoclinic

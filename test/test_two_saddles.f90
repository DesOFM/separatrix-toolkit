program test_two_saddles
    use kinds, only: dp
    use separatrix_types, only: grid_field,saddle_point,separatrix_network
    use field_provider, only: sample_function_on_grid
    use critical_points, only: find_all_saddles
    use saddle_refinement, only: refine_saddles_quadratic
    use separatrix_builder, only: build_all_separatrices
    implicit none

    type(grid_field) :: f
    type(saddle_point), allocatable :: s(:)
    type(separatrix_network), allocatable :: n(:)
    integer :: i,ns

    allocate(f%x(161),f%y(161))
    do i=1,161
        f%x(i)=-2.0_dp+0.025_dp*real(i-1,dp)
        f%y(i)=-2.0_dp+0.025_dp*real(i-1,dp)
    end do
    call sample_function_on_grid(f,psi)
    call find_all_saddles(f,s,ns)
    call refine_saddles_quadratic(f,s)
    ns=size(s)
    if(ns/=2) error stop 'two-saddle test: expected two saddles'
    call build_all_separatrices(f,s,n,level_tol=1.0e-8_dp)
    if(size(n)/=1) error stop 'two-saddle test: expected one connected network'
    if(n(1)%nbranches/=6) error stop 'two-saddle test: expected six branches'
    print *,'two-saddle heteroclinic test: PASS'
contains
    function psi(x,y) result(v)
        real(dp),intent(in)::x,y
        real(dp)::v
        v=(x*x-1.0_dp)**2-y*y
    end function psi
end program test_two_saddles

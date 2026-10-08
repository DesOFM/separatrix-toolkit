program test_single_saddle
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

    allocate(f%x(81),f%y(81))
    do i=1,81
        f%x(i)=-2.0_dp+0.05_dp*real(i-1,dp)
        f%y(i)=-2.0_dp+0.05_dp*real(i-1,dp)
    end do
    call sample_function_on_grid(f,psi)
    call find_all_saddles(f,s,ns)
    call refine_saddles_quadratic(f,s)
    ns=size(s)
    if(ns/=1) error stop 'single-saddle test: expected one saddle'
    call build_all_separatrices(f,s,n)
    if(size(n)/=1) error stop 'single-saddle test: expected one network'
    if(n(1)%nbranches/=4) error stop 'single-saddle test: expected four branches'
    print *,'single-saddle test: PASS'
contains
    function psi(x,y) result(v)
        real(dp),intent(in)::x,y
        real(dp)::v
        v=(x-0.017_dp)**2-(y+0.013_dp)**2
    end function psi
end program test_single_saddle

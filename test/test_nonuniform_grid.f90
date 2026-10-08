program test_nonuniform_grid
    use kinds, only: dp
    use separatrix_types, only: grid_field,saddle_point,separatrix_network
    use field_provider, only: sample_function_on_grid
    use critical_points, only: find_all_saddles
    use saddle_refinement, only: refine_saddles_quadratic
    use separatrix_builder, only: build_all_separatrices
    implicit none
    type(grid_field)::f;type(saddle_point),allocatable::s(:);type(separatrix_network),allocatable::n(:)
    integer::i,ns
    allocate(f%x(161),f%y(151))
    do i=1,161;f%x(i)=-2.0_dp+4.0_dp*(real(i-1,dp)/160.0_dp)**1.35_dp;end do
    do i=1,151;f%y(i)=-1.8_dp+3.6_dp*(real(i-1,dp)/150.0_dp)**1.20_dp;end do
    call sample_function_on_grid(f,psi);call find_all_saddles(f,s,ns);call refine_saddles_quadratic(f,s)
    if(size(s)/=1)error stop 'nonuniform-grid test: expected one saddle'
    if(hypot(s(1)%x-0.37_dp,s(1)%y+0.21_dp)>2.0e-3_dp)error stop 'nonuniform-grid test: saddle inaccurate'
    s(1)%id=1;call build_all_separatrices(f,s,n)
    if(size(n)/=1)error stop 'nonuniform-grid test: expected one network'
    print *,'nonuniform-grid test: PASS'
contains
    function psi(x,y) result(v)
        real(dp),intent(in)::x,y;real(dp)::v
        v=(x-0.37_dp)**2-(y+0.21_dp)**2
    end function psi
end program test_nonuniform_grid

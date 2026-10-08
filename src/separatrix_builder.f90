module separatrix_builder
    use kinds, only: dp
    use separatrix_types, only: grid_field,saddle_point,separatrix_network
    use contour_marching, only: segment_set,build_level_segments
    use saddle_topology, only: reconstruct_level_saddles
    use separatrix_graph, only: build_networks,append_network_array
    use grid_utils, only: characteristic_spacing
    use quality_metrics, only: evaluate_network_quality
    implicit none
    private

    public :: build_all_separatrices

contains

    subroutine build_all_separatrices(field,saddles,networks,level_tol,radius_factor,connect_tol,min_cos)
        type(grid_field), intent(in) :: field
        type(saddle_point), intent(in) :: saddles(:)
        type(separatrix_network), allocatable, intent(out) :: networks(:)
        real(dp), intent(in), optional :: level_tol,radius_factor,connect_tol,min_cos

        logical, allocatable :: assigned(:),ingroup(:)
        type(saddle_point), allocatable :: group(:)
        type(separatrix_network), allocatable :: localnets(:)
        type(segment_set) :: seg
        integer :: i,n,ng
        real(dp) :: lt,rf,ct,mc,radius,scale,level,h

        allocate(networks(0));n=size(saddles);if(n==0)return
        h=characteristic_spacing(field)
        rf=2.0_dp;if(present(radius_factor))then;if(radius_factor>0.0_dp)rf=radius_factor;end if
        ct=1.0e-7_dp*max(1.0_dp,maxval(abs(field%x)),maxval(abs(field%y)))
        if(present(connect_tol))then;if(connect_tol>0.0_dp)ct=connect_tol;end if
        mc=0.90_dp;if(present(min_cos))then;if(min_cos>0.0_dp.and.min_cos<=1.0_dp)mc=min_cos;end if
        scale=max(1.0_dp,maxval(abs(field%psi)))
        lt=1.0e-6_dp*scale;if(present(level_tol))then;if(level_tol>0.0_dp)lt=level_tol;end if
        radius=rf*h

        allocate(assigned(n),ingroup(n));assigned=.false.
        do i=1,n
            if(assigned(i))cycle
            level=saddles(i)%psi
            ingroup=(.not.assigned).and.(abs([(saddles(ng)%psi,ng=1,n)]-level)<=lt)
            ng=count(ingroup);allocate(group(ng));group=pack(saddles,ingroup);assigned=assigned.or.ingroup

            call build_level_segments(field,level,seg)
            call reconstruct_level_saddles(seg,group,radius,ct,mc)
            call build_networks(seg,group,ct,localnets)
            call evaluate_network_quality(field,level,localnets)
            call append_network_array(networks,localnets)

            if(allocated(localnets))deallocate(localnets)
            if(allocated(group))deallocate(group)
        end do
        deallocate(assigned,ingroup)
    end subroutine build_all_separatrices

end module separatrix_builder

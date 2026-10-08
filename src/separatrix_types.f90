module separatrix_types
    use kinds, only: dp
    implicit none
    private

    public :: grid_field
    public :: saddle_point
    public :: branch
    public :: separatrix_network

    type :: grid_field
        integer :: nx = 0
        integer :: ny = 0
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: psi(:,:)
    end type grid_field

    type :: saddle_point
        integer :: id = 0
        real(dp) :: x = 0.0_dp
        real(dp) :: y = 0.0_dp
        real(dp) :: psi = 0.0_dp
        real(dp) :: h11 = 0.0_dp
        real(dp) :: h12 = 0.0_dp
        real(dp) :: h22 = 0.0_dp
        real(dp) :: det_h = 0.0_dp
    end type saddle_point

    type :: branch
        integer :: id = 0
        integer :: saddle_start = 0
        integer :: saddle_end = 0
        integer :: n = 0
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: y(:)

        ! Quality metrics. They are filled by quality_metrics.f90.
        real(dp) :: arclength = 0.0_dp
        real(dp) :: rms_level_error = 0.0_dp
        real(dp) :: max_level_error = 0.0_dp
    end type branch

    type :: separatrix_network
        integer :: id = 0
        integer :: nbranches = 0
        real(dp) :: level = 0.0_dp
        real(dp) :: rms_level_error = 0.0_dp
        real(dp) :: max_level_error = 0.0_dp
        type(branch), allocatable :: branches(:)
    end type separatrix_network

end module separatrix_types

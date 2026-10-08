module quality_metrics
    use kinds, only: dp
    use separatrix_types, only: grid_field, separatrix_network
    use grid_utils, only: bilinear_interpolate
    implicit none
    private

    public :: evaluate_network_quality

contains

    subroutine evaluate_network_quality(field,level,networks)
        type(grid_field), intent(in) :: field
        real(dp), intent(in) :: level
        type(separatrix_network), intent(inout) :: networks(:)
        integer :: i,j,k,nall
        real(dp) :: value,err,sumsq,maxerr,total_len
        logical :: inside

        do i=1,size(networks)
            networks(i)%level=level
            sumsq=0.0_dp;maxerr=0.0_dp;nall=0

            do j=1,networks(i)%nbranches
                networks(i)%branches(j)%arclength=0.0_dp
                networks(i)%branches(j)%rms_level_error=0.0_dp
                networks(i)%branches(j)%max_level_error=0.0_dp
                total_len=0.0_dp

                do k=1,networks(i)%branches(j)%n
                    value=bilinear_interpolate(field,networks(i)%branches(j)%x(k), &
                                               networks(i)%branches(j)%y(k),inside)
                    if (inside) then
                        err=abs(value-level)
                        networks(i)%branches(j)%rms_level_error = &
                            networks(i)%branches(j)%rms_level_error + err*err
                        networks(i)%branches(j)%max_level_error = &
                            max(networks(i)%branches(j)%max_level_error,err)
                        sumsq=sumsq+err*err;maxerr=max(maxerr,err);nall=nall+1
                    end if
                    if (k>1) then
                        total_len=total_len+hypot(networks(i)%branches(j)%x(k)-networks(i)%branches(j)%x(k-1), &
                                                 networks(i)%branches(j)%y(k)-networks(i)%branches(j)%y(k-1))
                    end if
                end do

                networks(i)%branches(j)%arclength=total_len
                if (networks(i)%branches(j)%n>0) then
                    networks(i)%branches(j)%rms_level_error = sqrt( &
                        networks(i)%branches(j)%rms_level_error/real(networks(i)%branches(j)%n,dp))
                end if
            end do

            networks(i)%max_level_error=maxerr
            if (nall>0) networks(i)%rms_level_error=sqrt(sumsq/real(nall,dp))
        end do
    end subroutine evaluate_network_quality

end module quality_metrics

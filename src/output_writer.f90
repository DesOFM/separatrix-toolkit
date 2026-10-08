module output_writer
    use kinds, only: dp
    use separatrix_types, only: saddle_point, separatrix_network
    implicit none
    private

    public :: write_saddles
    public :: write_separatrices

contains

    subroutine write_saddles(filename, layer_id, saddles)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: layer_id
        type(saddle_point), intent(in) :: saddles(:)

        integer :: unit, i

        open(newunit=unit, file=filename, status='replace', action='write')
        write(unit,'(A)') '# layer saddle_id x y psi h11 h12 h22 det_h'

        do i = 1, size(saddles)
            write(unit,'(2I8,7ES24.15)') layer_id, saddles(i)%id, &
                saddles(i)%x, saddles(i)%y, saddles(i)%psi, &
                saddles(i)%h11, saddles(i)%h12, saddles(i)%h22, saddles(i)%det_h
        end do

        close(unit)
    end subroutine write_saddles


    subroutine write_separatrices(filename, layer_id, networks)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: layer_id
        type(separatrix_network), intent(in) :: networks(:)

        integer :: unit, isep, ib, ip

        open(newunit=unit, file=filename, status='replace', action='write')
        write(unit,'(A)') '# layer sep_id branch_id point_id x y saddle_start saddle_end'

        do isep = 1, size(networks)
            do ib = 1, networks(isep)%nbranches
                do ip = 1, networks(isep)%branches(ib)%n
                    write(unit,'(4I8,2ES24.15,2I8)') &
                        layer_id, networks(isep)%id, &
                        networks(isep)%branches(ib)%id, ip, &
                        networks(isep)%branches(ib)%x(ip), &
                        networks(isep)%branches(ib)%y(ip), &
                        networks(isep)%branches(ib)%saddle_start, &
                        networks(isep)%branches(ib)%saddle_end
                end do
            end do
        end do

        close(unit)
    end subroutine write_separatrices

end module output_writer

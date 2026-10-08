module diagnostics
    use kinds, only: dp
    use separatrix_types, only: grid_field, saddle_point, separatrix_network
    implicit none
    private
    public :: write_diagnostics_report, print_topology_summary
contains
    subroutine write_diagnostics_report(filename,field,saddles,networks)
        character(len=*), intent(in) :: filename
        type(grid_field), intent(in) :: field
        type(saddle_point), intent(in) :: saddles(:)
        type(separatrix_network), intent(in) :: networks(:)
        integer :: unit,i,j,nbranches,npoints,nhomo,nhetero,nboundary,nother
        call topology_counts(networks,nbranches,npoints,nhomo,nhetero,nboundary,nother)
        open(newunit=unit,file=trim(filename),status='replace',action='write')
        write(unit,'(A)')'SEPARATRIX TOOLKIT DIAGNOSTICS';write(unit,'(A)')'============================='
        write(unit,'(A,I0)')'nx = ',field%nx;write(unit,'(A,I0)')'ny = ',field%ny
        write(unit,'(A,2ES24.15)')'x range = ',minval(field%x),maxval(field%x)
        write(unit,'(A,2ES24.15)')'y range = ',minval(field%y),maxval(field%y)
        write(unit,'(A,2ES24.15)')'psi range = ',minval(field%psi),maxval(field%psi)
        write(unit,'(A,I0)')'saddles = ',size(saddles);write(unit,'(A,I0)')'networks = ',size(networks)
        write(unit,'(A,I0)')'branches = ',nbranches;write(unit,'(A,I0)')'points = ',npoints
        write(unit,'(A,I0)')'homoclinic branches = ',nhomo;write(unit,'(A,I0)')'heteroclinic branches = ',nhetero
        write(unit,'(A,I0)')'boundary branches = ',nboundary;write(unit,'(A,I0)')'other branches = ',nother
        write(unit,'(A)')'';write(unit,'(A)')'SADDLES';write(unit,'(A)')'-------'
        do i=1,size(saddles)
            write(unit,'(A,I0)')'saddle_id = ',saddles(i)%id
            write(unit,'(A,ES24.15)')'  x = ',saddles(i)%x;write(unit,'(A,ES24.15)')'  y = ',saddles(i)%y
            write(unit,'(A,ES24.15)')'  psi = ',saddles(i)%psi;write(unit,'(A,ES24.15)')'  det(H) = ',saddles(i)%det_h
        end do
        write(unit,'(A)')'';write(unit,'(A)')'NETWORKS';write(unit,'(A)')'--------'
        do i=1,size(networks)
            write(unit,'(A,I0,A,I0)')'network ',networks(i)%id,': branches = ',networks(i)%nbranches
            write(unit,'(A,ES14.6)')'  level = ',networks(i)%level
            write(unit,'(A,ES14.6)')'  rms |psi-level| = ',networks(i)%rms_level_error
            write(unit,'(A,ES14.6)')'  max |psi-level| = ',networks(i)%max_level_error
            do j=1,networks(i)%nbranches
                write(unit,'(A,I0,A,I0,A,I0,A,I0,A,A)')'  branch ',networks(i)%branches(j)%id, &
                    ': points=',networks(i)%branches(j)%n,', saddle_start=',networks(i)%branches(j)%saddle_start, &
                    ', saddle_end=',networks(i)%branches(j)%saddle_end,', type=', &
                    trim(branch_type(networks(i)%branches(j)%saddle_start,networks(i)%branches(j)%saddle_end))
                write(unit,'(A,ES14.6)')'    arclength = ',networks(i)%branches(j)%arclength
                write(unit,'(A,ES14.6)')'    rms |psi-level| = ',networks(i)%branches(j)%rms_level_error
                write(unit,'(A,ES14.6)')'    max |psi-level| = ',networks(i)%branches(j)%max_level_error
            end do
        end do
        close(unit)
    end subroutine write_diagnostics_report

    subroutine print_topology_summary(saddles,networks)
        type(saddle_point), intent(in) :: saddles(:);type(separatrix_network), intent(in) :: networks(:)
        integer :: nbranches,npoints,nhomo,nhetero,nboundary,nother,i
        call topology_counts(networks,nbranches,npoints,nhomo,nhetero,nboundary,nother)
        write(*,'(A,I0)')'Detected saddle points: ',size(saddles);write(*,'(A,I0)')'Detected separatrix networks: ',size(networks)
        write(*,'(A,I0)')'Total branches: ',nbranches;write(*,'(A,I0)')'  homoclinic: ',nhomo
        write(*,'(A,I0)')'  heteroclinic: ',nhetero;write(*,'(A,I0)')'  boundary-going: ',nboundary
        if(nother>0)write(*,'(A,I0)')'  other/unclassified: ',nother
        do i=1,size(networks)
            write(*,'(A,I0,A,ES10.3,A,ES10.3)')'  network ',networks(i)%id, &
                ': rms residual=',networks(i)%rms_level_error,', max=',networks(i)%max_level_error
        end do
    end subroutine print_topology_summary

    subroutine topology_counts(networks,nbranches,npoints,nhomo,nhetero,nboundary,nother)
        type(separatrix_network), intent(in) :: networks(:)
        integer,intent(out)::nbranches,npoints,nhomo,nhetero,nboundary,nother
        integer::i,j,s1,s2
        nbranches=0;npoints=0;nhomo=0;nhetero=0;nboundary=0;nother=0
        do i=1,size(networks);do j=1,networks(i)%nbranches
            nbranches=nbranches+1;npoints=npoints+networks(i)%branches(j)%n
            s1=networks(i)%branches(j)%saddle_start;s2=networks(i)%branches(j)%saddle_end
            if(s1>0.and.s2==s1)then;nhomo=nhomo+1
            else if(s1>0.and.s2>0.and.s1/=s2)then;nhetero=nhetero+1
            else if((s1>0.and.s2==0).or.(s2>0.and.s1==0))then;nboundary=nboundary+1
            else;nother=nother+1;end if
        end do;end do
    end subroutine topology_counts

    function branch_type(s1,s2) result(name)
        integer,intent(in)::s1,s2;character(len=32)::name
        if(s1>0.and.s2==s1)then;name='homoclinic'
        else if(s1>0.and.s2>0.and.s1/=s2)then;name='heteroclinic'
        else if((s1>0.and.s2==0).or.(s2>0.and.s1==0))then;name='boundary'
        else;name='other';end if
    end function branch_type
end module diagnostics

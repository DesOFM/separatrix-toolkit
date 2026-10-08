module critical_points
    use kinds, only: dp
    use separatrix_types, only: grid_field, saddle_point
    use grid_utils, only: validate_rectilinear_grid, characteristic_spacing
    implicit none
    private

    public :: find_all_saddles

contains

    subroutine find_all_saddles(field, saddles, nsaddles, root_tol, merge_tol)
        type(grid_field), intent(in) :: field
        type(saddle_point), allocatable, intent(out) :: saddles(:)
        integer, intent(out) :: nsaddles
        real(dp), intent(in), optional :: root_tol, merge_tol

        real(dp), allocatable :: gx(:,:), gy(:,:)
        type(saddle_point), allocatable :: candidates(:)
        integer :: i,j,ncand,maxcand
        real(dp) :: rt,mt,u,v,xr,yr,h
        logical :: ok

        call validate_rectilinear_grid(field)
        if (field%nx<4 .or. field%ny<4) error stop 'Grid is too small'

        h=characteristic_spacing(field)
        rt=1.0e-10_dp
        if (present(root_tol)) then
            if (root_tol > 0.0_dp) rt=root_tol
        end if
        mt=0.75_dp*h
        if (present(merge_tol)) then
            if (merge_tol > 0.0_dp) mt=merge_tol
        end if

        allocate(gx(field%nx,field%ny),gy(field%nx,field%ny))
        gx=0.0_dp;gy=0.0_dp

        do i=2,field%nx-1
            do j=1,field%ny
                gx(i,j)=first_derivative_3pt(field%x(i-1),field%x(i),field%x(i+1), &
                                              field%psi(i-1,j),field%psi(i,j),field%psi(i+1,j))
            end do
        end do
        do i=1,field%nx
            do j=2,field%ny-1
                gy(i,j)=first_derivative_3pt(field%y(j-1),field%y(j),field%y(j+1), &
                                              field%psi(i,j-1),field%psi(i,j),field%psi(i,j+1))
            end do
        end do

        maxcand=(field%nx-3)*(field%ny-3)
        allocate(candidates(maxcand));ncand=0

        do i=2,field%nx-2
            do j=2,field%ny-2
                if (.not.brackets_zero(gx(i,j),gx(i+1,j),gx(i+1,j+1),gx(i,j+1))) cycle
                if (.not.brackets_zero(gy(i,j),gy(i+1,j),gy(i+1,j+1),gy(i,j+1))) cycle

                call solve_bilinear_gradient( &
                    gx(i,j),gx(i+1,j),gx(i+1,j+1),gx(i,j+1), &
                    gy(i,j),gy(i+1,j),gy(i+1,j+1),gy(i,j+1),u,v,rt,ok)
                if (.not.ok) cycle
                if (u < -0.05_dp .or. u > 1.05_dp .or. v < -0.05_dp .or. v > 1.05_dp) cycle

                xr=field%x(i)+u*(field%x(i+1)-field%x(i))
                yr=field%y(j)+v*(field%y(j+1)-field%y(j))
                ncand=ncand+1
                candidates(ncand)%x=xr;candidates(ncand)%y=yr
                candidates(ncand)%psi=bilinear_value(field%psi(i,j),field%psi(i+1,j), &
                    field%psi(i+1,j+1),field%psi(i,j+1),u,v)
            end do
        end do

        call merge_candidates(candidates,ncand,mt,saddles,nsaddles)
        do i=1,nsaddles
            saddles(i)%id=i
        end do
        deallocate(gx,gy,candidates)
    end subroutine find_all_saddles


    pure real(dp) function first_derivative_3pt(xm,x0,xp,fm,f0,fp) result(df)
        real(dp), intent(in) :: xm,x0,xp,fm,f0,fp
        real(dp) :: hm,hp
        hm=x0-xm;hp=xp-x0
        df = -hp/(hm*(hm+hp))*fm + (hp-hm)/(hm*hp)*f0 + hm/(hp*(hm+hp))*fp
    end function first_derivative_3pt


    logical function brackets_zero(a,b,c,d) result(ok)
        real(dp), intent(in) :: a,b,c,d
        ok=(min(a,b,c,d)<=0.0_dp .and. max(a,b,c,d)>=0.0_dp)
    end function brackets_zero


    subroutine solve_bilinear_gradient(gx00,gx10,gx11,gx01,gy00,gy10,gy11,gy01,u,v,tol,ok)
        real(dp), intent(in) :: gx00,gx10,gx11,gx01,gy00,gy10,gy11,gy01,tol
        real(dp), intent(out) :: u,v
        logical, intent(out) :: ok
        real(dp) :: ax,bx,cx,dx,ay,by,cy,dy,fx,fy,j11,j12,j21,j22,det,du,dv,scale
        integer :: iter

        call bilinear_coeff(gx00,gx10,gx11,gx01,ax,bx,cx,dx)
        call bilinear_coeff(gy00,gy10,gy11,gy01,ay,by,cy,dy)
        u=0.5_dp;v=0.5_dp;ok=.false.
        scale=max(1.0_dp,abs(gx00),abs(gx10),abs(gx11),abs(gx01),abs(gy00),abs(gy10),abs(gy11),abs(gy01))

        do iter=1,30
            fx=ax+bx*u+cx*v+dx*u*v;fy=ay+by*u+cy*v+dy*u*v
            j11=bx+dx*v;j12=cx+dx*u;j21=by+dy*v;j22=cy+dy*u
            det=j11*j22-j12*j21
            if (abs(det)<=1.0e-18_dp*scale) return
            du=(-fx*j22+j12*fy)/det;dv=(-j11*fy+j21*fx)/det
            u=u+du;v=v+dv
            if (max(abs(du),abs(dv))<=tol) then
                fx=ax+bx*u+cx*v+dx*u*v;fy=ay+by*u+cy*v+dy*u*v
                ok=(max(abs(fx),abs(fy))<=1.0e-7_dp*scale);return
            end if
            if (abs(u)>4.0_dp .or. abs(v)>4.0_dp) return
        end do
    end subroutine solve_bilinear_gradient


    subroutine bilinear_coeff(f00,f10,f11,f01,a,b,c,d)
        real(dp), intent(in) :: f00,f10,f11,f01
        real(dp), intent(out) :: a,b,c,d
        a=f00;b=f10-f00;c=f01-f00;d=f11-f10-f01+f00
    end subroutine bilinear_coeff

    real(dp) function bilinear_value(f00,f10,f11,f01,u,v) result(f)
        real(dp), intent(in) :: f00,f10,f11,f01,u,v
        real(dp) :: a,b,c,d
        call bilinear_coeff(f00,f10,f11,f01,a,b,c,d);f=a+b*u+c*v+d*u*v
    end function bilinear_value

    subroutine merge_candidates(candidates,ncand,merge_tol,saddles,nsaddles)
        type(saddle_point), intent(in) :: candidates(:)
        integer, intent(in) :: ncand
        real(dp), intent(in) :: merge_tol
        type(saddle_point), allocatable, intent(out) :: saddles(:)
        integer, intent(out) :: nsaddles
        type(saddle_point), allocatable :: temp(:)
        integer, allocatable :: counts(:)
        integer :: i,j,best
        real(dp) :: d

        allocate(temp(max(1,ncand)),counts(max(1,ncand)));nsaddles=0;counts=0
        do i=1,ncand
            best=0
            do j=1,nsaddles
                d=hypot(candidates(i)%x-temp(j)%x,candidates(i)%y-temp(j)%y)
                if (d<=merge_tol) then;best=j;exit;end if
            end do
            if (best==0) then
                nsaddles=nsaddles+1;temp(nsaddles)=candidates(i);counts(nsaddles)=1
            else
                counts(best)=counts(best)+1
                temp(best)%x=temp(best)%x+(candidates(i)%x-temp(best)%x)/real(counts(best),dp)
                temp(best)%y=temp(best)%y+(candidates(i)%y-temp(best)%y)/real(counts(best),dp)
                temp(best)%psi=temp(best)%psi+(candidates(i)%psi-temp(best)%psi)/real(counts(best),dp)
            end if
        end do
        allocate(saddles(nsaddles));if(nsaddles>0)saddles=temp(1:nsaddles)
        deallocate(temp,counts)
    end subroutine merge_candidates

end module critical_points

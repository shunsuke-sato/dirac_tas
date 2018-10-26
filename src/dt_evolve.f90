subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer,parameter :: n_Runge_Kutta_4th = 0
  integer,parameter :: imethod = n_Runge_Kutta_4th

  select case(imethod)
  case(n_Runge_Kutta_4th)
    call dt_evolve_RK4(it)
  case default
    stop 'Error in dt_evolve'
  end select
  
end subroutine dt_evolve
!===================================================================
subroutine dt_evolve_RK4(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ik
  real(8) :: kx_t, ky_t, Ez_t
  complex(8) :: zrho_t(4,4),zrho_rk(4,4,0:4)

  do ik = nk_start, nk_end

    zrho_rk(:,:,0) = zrho_dm(:,:,ik)

! set fields t=t
    kx_t = kx(ik) + Act(1,it)
    ky_t = ky(ik) + Act(2,it)
    Ezt = Ezt(it)

! RK-1
    zrho_t(:,:) = zrho_rk(:,:,0)
    call apply_Lrho(zrho_t, zrho_rk(:,:,1), kx_t,ky_t,Ezt)

! set fields t=t + dt/2
    kx_t = kx(ik) + Act_dt2(1,it)
    ky_t = ky(ik) + Act_dt2(2,it)
    Ezt = Ezt_dt2(it)

! RK-2
    zrho_t(:,:) = zrho_rk(:,:,0) + 0.5d0*dt*zrho_rk(:,:,1)
    call apply_Lrho(zrho_t, zrho_rk(:,:,2), kx_t,ky_t,Ezt)

! RK-3
    zrho_t(:,:) = zrho_rk(:,:,0) + 0.5d0*dt*zrho_rk(:,:,2)
    call apply_Lrho(zrho_t, zrho_rk(:,:,3), kx_t,ky_t,Ezt)

! set fields t=t + dt
    kx_t = kx(ik) + Act(1,it+1)
    ky_t = ky(ik) + Act(2,it+1)
    Ezt = Ezt(it+1)

! RK-4
    zrho_t(:,:) = zrho_rk(:,:,0) + dt*zrho_rk(:,:,3)
    call apply_Lrho(zrho_t, zrho_rk(:,:,4), kx_t,ky_t,Ezt)




  end do



end subroutine dt_evolve_RK4
!===================================================================
!===================================================================

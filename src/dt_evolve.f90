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

    zrho_dm(:,:,ik) = zrho_dm(:,:,ik) + dt/6d0*( &
      zrho_rk(:,:,1) + 2d0*zrho_rk(:,:,2) + 2d0*zrho_rk(:,:,3) + zrho_rk(:,:,4))
      


  end do



end subroutine dt_evolve_RK4
!===================================================================
subroutine apply_Lrho(zrho_in, zLrho_out, kx_t, ky_t,Ezt)
  use global_variables
  implicit none
  complex(8),intent(in) :: zrho_in(4,4)
  complex(8),intent(out) :: zLrho_out(4,4)
  complex(8) :: zHmat_t(4,4)

  zhmat_t(1,1) = 0d0
  zhmat_t(2,1) = v_Fermi*(tau_chiral*kx_t+zI*kt_t)
  zhmat_t(3,1) = dip_core*Ezt
  zhmat_t(4,1) = 0d0

  zhmat_t(1,2) = v_Fermi*(tau_chiral*kx_t-zI*kt_t)
  zhmat_t(2,2) = 0d0
  zhmat_t(3,2) = 0d0
  zhmat_t(4,2) = dip_core*Ezt

  zhmat_t(1,3) = dip_core*Ezt
  zhmat_t(2,3) = 0d0
  zhmat_t(3,3) = eps_core
  zhmat_t(4,3) = 0d0

  zhmat_t(1,4) = 0d0
  zhmat_t(2,4) = dip_core*Ezt
  zhmat_t(3,4) = 0d0
  zhmat_t(4,4) = eps_core

  zLrho_out = -zI*( matmul(zhmat_t,zrho_in) - matmul(zrho_in,zhmat_t) )

end subroutine apply_Lrho


!===================================================================

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
  real(8) :: kx_t, ky_t, Ezt_t
  complex(8) :: zrho_t(4,4),zrho_rk(4,4,0:4)

  do ik = nk_start, nk_end

    zrho_rk(:,:,0) = zrho_dm(:,:,ik)

! set fields t=t
    kx_t = kx(ik) + Act(1,it)
    ky_t = ky(ik) + Act(2,it)
    Ezt_t = Ezt(it)

! RK-1
    zrho_t(:,:) = zrho_rk(:,:,0)
    call apply_Lrho(zrho_t, zrho_rk(:,:,1), kx_t,ky_t,Ezt_t)

! set fields t=t + dt/2
    kx_t = kx(ik) + Act_dt2(1,it)
    ky_t = ky(ik) + Act_dt2(2,it)
    Ezt_t = Ezt_dt2(it)

! RK-2
    zrho_t(:,:) = zrho_rk(:,:,0) + 0.5d0*dt*zrho_rk(:,:,1)
    call apply_Lrho(zrho_t, zrho_rk(:,:,2), kx_t,ky_t,Ezt_t)

! RK-3
    zrho_t(:,:) = zrho_rk(:,:,0) + 0.5d0*dt*zrho_rk(:,:,2)
    call apply_Lrho(zrho_t, zrho_rk(:,:,3), kx_t,ky_t,Ezt_t)

! set fields t=t + dt
    kx_t = kx(ik) + Act(1,it+1)
    ky_t = ky(ik) + Act(2,it+1)
    Ezt_t = Ezt(it+1)

! RK-4
    zrho_t(:,:) = zrho_rk(:,:,0) + dt*zrho_rk(:,:,3)
    call apply_Lrho(zrho_t, zrho_rk(:,:,4), kx_t,ky_t,Ezt_t)

    zrho_dm(:,:,ik) = zrho_dm(:,:,ik) + dt/6d0*( &
      zrho_rk(:,:,1) + 2d0*zrho_rk(:,:,2) + 2d0*zrho_rk(:,:,3) + zrho_rk(:,:,4))
      


  end do



end subroutine dt_evolve_RK4
!===================================================================
subroutine apply_Lrho(zrho_in, zLrho_out, kx_t, ky_t,Ezt_t)
  use global_variables
  implicit none
  complex(8),intent(in) :: zrho_in(4,4)
  complex(8),intent(out) :: zLrho_out(4,4)
  real(8),intent(in) :: kx_t, ky_t, Ezt_t
  complex(8) :: zHmat_t(4,4)
  complex(8) :: zeigv(4,4), zrho_col(4,4)
  real(8) :: occ_v, occ_c, occ_core, eps_v, eps_c, phi

  zhmat_t(1,1) = 0d0
  zhmat_t(2,1) = v_Fermi*(tau_chiral*kx_t+zI*ky_t)
  zhmat_t(3,1) = dip_core*Ezt_t
  zhmat_t(4,1) = 0d0

  zhmat_t(1,2) = v_Fermi*(tau_chiral*kx_t-zI*ky_t)
  zhmat_t(2,2) = 0d0
  zhmat_t(3,2) = 0d0
  zhmat_t(4,2) = dip_core*Ezt_t

  zhmat_t(1,3) = dip_core*Ezt_t
  zhmat_t(2,3) = 0d0
  zhmat_t(3,3) = eps_core
  zhmat_t(4,3) = 0d0

  zhmat_t(1,4) = 0d0
  zhmat_t(2,4) = dip_core*Ezt_t
  zhmat_t(3,4) = 0d0
  zhmat_t(4,4) = eps_core

  zLrho_out = -zI*( matmul(zhmat_t,zrho_in) - matmul(zrho_in,zhmat_t) )
!  call commutator4x4(zHmat_t,zrho_in,zLrho_out)
!  zLrho_out = -zI*zLrho_out

! relaxation
  eps_c = v_Fermi*sqrt(kx_t**2 + ky_t**2)
  eps_v = - eps_c
  occ_v = Fermi_Dirac_distribution(eps_v, mu_F, kbT)
  occ_c = Fermi_Dirac_distribution(eps_c, mu_F, kbT)
  occ_core = Fermi_Dirac_distribution(eps_core, mu_F, kbT)

  if(kx_t*tau_chiral > 0d0)then
    phi = atan(ky_t/kx_t*tau_chiral)
  else if(kx_t*tau_chiral < 0d0)then
    phi = atan(ky_t/kx_t*tau_chiral) + pi
  else
    if(ky_t > 0d0)then
      phi = 0.5d0*pi
    else
      phi = -0.5d0*pi
    end if
  end if

  zeigv = 0d0
  zeigv(1,1:2) = 1d0/sqrt(2d0)
  zeigv(2,1) = -exp(zI*phi)/sqrt(2d0)
  zeigv(2,2) =  exp(zI*phi)/sqrt(2d0)
  zeigv(3,3) = 1d0
  zeigv(4,4) = 1d0

  zrho_col = matmul(transpose(conjg(zeigv)),matmul(zrho_in,zeigv))
  zrho_col(1,1) = -(zrho_col(1,1)-occ_v)/T1_relax
  zrho_col(2:4,1) = -zrho_col(2:4,1)/T2_relax
  zrho_col(1,2) = -zrho_col(1,2)/T2_relax
  zrho_col(2,2) = -(zrho_col(2,2)-occ_c)/T1_relax
  zrho_col(3:4,2) = -zrho_col(3:4,2)/T2_relax
  zrho_col(1:2,3) = -zrho_col(1:2,3)/T2_relax
  zrho_col(3,3) = -(zrho_col(3,3)-occ_core)/T1_relax
  zrho_col(4,3) = -zrho_col(4,3)/T2_relax
  zrho_col(1:3,4) = -zrho_col(1:3,4)/T2_relax
  zrho_col(4,4) = -(zrho_col(4,4)-occ_core)/T1_relax

  zrho_col = matmul( matmul(zeigv, zrho_col), transpose(conjg(zeigv)))

  zLrho_out = zLrho_out + zrho_col

end subroutine apply_Lrho
subroutine commutator4x4(A,B,C)
  implicit none
  complex(8),intent(in)  :: A(4,4),B(4,4)
  complex(8),intent(out) :: C(4,4)

  C(1,1) = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1) + A(1,4)*B(4,1) &
          -B(1,1)*A(1,1) - B(1,2)*A(2,1) - B(1,3)*A(3,1) - B(1,4)*A(4,1)

  C(2,1) = A(2,1)*B(1,1) + A(2,2)*B(2,1) + A(2,3)*B(3,1) + A(2,4)*B(4,1) &
          -B(2,1)*A(1,1) - B(2,2)*A(2,1) - B(2,3)*A(3,1) - B(2,4)*A(4,1)

  C(3,1) = A(3,1)*B(1,1) + A(3,2)*B(2,1) + A(3,3)*B(3,1) + A(3,4)*B(4,1) &
          -B(3,1)*A(1,1) - B(3,2)*A(2,1) - B(3,3)*A(3,1) - B(3,4)*A(4,1)

  C(4,1) = A(4,1)*B(1,1) + A(4,2)*B(2,1) + A(4,3)*B(3,1) + A(4,4)*B(4,1) &
          -B(4,1)*A(1,1) - B(4,2)*A(2,1) - B(4,3)*A(3,1) - B(4,4)*A(4,1)

  C(1,2) = A(1,1)*B(1,2) + A(1,2)*B(2,2) + A(1,3)*B(3,2) + A(1,4)*B(4,2) &
          -B(1,1)*A(1,2) - B(1,2)*A(2,2) - B(1,3)*A(3,2) - B(1,4)*A(4,2)

  C(2,2) = A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2) + A(2,4)*B(4,2) &
          -B(2,1)*A(1,2) - B(2,2)*A(2,2) - B(2,3)*A(3,2) - B(2,4)*A(4,2)

  C(3,2) = A(3,1)*B(1,2) + A(3,2)*B(2,2) + A(3,3)*B(3,2) + A(3,4)*B(4,2) &
          -B(3,1)*A(1,2) - B(3,2)*A(2,2) - B(3,3)*A(3,2) - B(3,4)*A(4,2)

  C(4,2) = A(4,1)*B(1,2) + A(4,2)*B(2,2) + A(4,3)*B(3,2) + A(4,4)*B(4,2) &
          -B(4,1)*A(1,2) - B(4,2)*A(2,2) - B(4,3)*A(3,2) - B(4,4)*A(4,2)

  C(1,3) = A(1,1)*B(1,3) + A(1,2)*B(2,3) + A(1,3)*B(3,3) + A(1,4)*B(4,3) &
          -B(1,1)*A(1,3) - B(1,2)*A(2,3) - B(1,3)*A(3,3) - B(1,4)*A(4,3)

  C(2,3) = A(2,1)*B(1,3) + A(2,2)*B(2,3) + A(2,3)*B(3,3) + A(2,4)*B(4,3) &
          -B(2,1)*A(1,3) - B(2,2)*A(2,3) - B(2,3)*A(3,3) - B(2,4)*A(4,3)

  C(3,3) = A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3) + A(3,4)*B(4,3) &
          -B(3,1)*A(1,3) - B(3,2)*A(2,3) - B(3,3)*A(3,3) - B(3,4)*A(4,3)

  C(4,3) = A(4,1)*B(1,3) + A(4,2)*B(2,3) + A(4,3)*B(3,3) + A(4,4)*B(4,3) &
          -B(4,1)*A(1,3) - B(4,2)*A(2,3) - B(4,3)*A(3,3) - B(4,4)*A(4,3)

  C(1,4) = A(1,1)*B(1,4) + A(1,2)*B(2,4) + A(1,3)*B(3,4) + A(1,4)*B(4,4) &
          -B(1,1)*A(1,4) - B(1,2)*A(2,4) - B(1,3)*A(3,4) - B(1,4)*A(4,4)

  C(2,4) = A(2,1)*B(1,4) + A(2,2)*B(2,4) + A(2,3)*B(3,4) + A(2,4)*B(4,4) &
          -B(2,1)*A(1,4) - B(2,2)*A(2,4) - B(2,3)*A(3,4) - B(2,4)*A(4,4)

  C(3,4) = A(3,1)*B(1,4) + A(3,2)*B(2,4) + A(3,3)*B(3,4) + A(3,4)*B(4,4) &
          -B(3,1)*A(1,4) - B(3,2)*A(2,4) - B(3,3)*A(3,4) - B(3,4)*A(4,4)

  C(4,4) = A(4,1)*B(1,4) + A(4,2)*B(2,4) + A(4,3)*B(3,4) + A(4,4)*B(4,4) &
          -B(4,1)*A(1,4) - B(4,2)*A(2,4) - B(4,3)*A(3,4) - B(4,4)*A(4,4)


end subroutine commutator4x4

!===================================================================

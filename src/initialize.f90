subroutine initialize
  use global_variables
  implicit none
  integer :: ik, ikx, iky
    
  nk = nkx*nky
  dkx = 2d0*kx_max/nkx
  dky = 2d0*ky_max/nky
      
  if(nk < comm_nproc_global) call error_finalize('Error: nk < # of MPI processes')
  nk_average = nk/comm_nproc_global
  nk_remainder = mod(nk,comm_nproc_global)
  if(comm_id_global+1 <= nk_remainder)then
    nk_start = 1 + comm_id_global*(nk_average+1)
    nk_end   = nk_start + (nk_average + 1) -1
  else
    nk_start = 1 + nk_remainder*(nk_average+1) + nk_average*(comm_id_global - nk_remainder)
    nk_end    = nk_start + nk_average  -1
  end if
      
!    do ik = 1, comm_nproc_global
!      if(ik-1 == comm_id_global)then
!        write(*,*)ik-1,nk_start,nk_end,nk
!      end if
!    end do
  call read_basic_input('kx_shift',kx_shift,val_default = 0d0)
  call read_basic_input('ky_shift',ky_shift,val_default = 0d0)
      
  allocate(kx(nk),ky(nk))
  allocate(kx0(nk),ky0(nk))    
  allocate(zrho_dm(4,4,nk_start:nk_end))
  allocate(ik_table(nkx,nky))
      
      ! initialize k-grids
  ik = 0
  do ikx = 1, nkx
    do iky = 1, nky
      ik = ik + 1
      
      kx0(ik) = -kx_max + dkx*ikx -0.5d0*dkx
      ky0(ik) = -ky_max + dky*iky -0.5d0*dky
      
      ik_table(ikx,iky) = ik
      
    end do
  end do
  kx0 = kx0 + kx_shift
  ky0 = ky0 + ky_shift


  kx = kx0
  ky = ky0


  call initialize_density_matrix
  
end subroutine initialize
!-------------------------------------------------------------------------------
subroutine initialize_density_matrix
  use global_variables
  implicit none
  integer :: ik
  real(8) :: eps_v, eps_c, occ_v, occ_c, phi
  complex(8) :: zeig_vec(2,2), zrho_dm2x2(2,2)
  real(8) :: occ_core

  do ik = nk_start, nk_end
    eps_c = v_Fermi*sqrt(kx(ik)**2 + ky(ik)**2)
    eps_v = - eps_c
    
    if(kx(ik)*tau_chiral > 0d0)then
      phi = atan(ky(ik)/kx(ik)*tau_chiral)
    else if(kx(ik)*tau_chiral < 0d0)then
      phi = atan(ky(ik)/kx(ik)*tau_chiral) + pi
    else
      if(ky(ik) > 0d0)then
        phi = 0.5d0*pi
      else
        phi = -0.5d0*pi
      end if
    end if
      
    zeig_vec(1,1:2) = 1d0/sqrt(2d0)
    zeig_vec(2,1) = -exp(zI*phi)/sqrt(2d0)
    zeig_vec(2,2) =  exp(zI*phi)/sqrt(2d0)
    
    occ_v = Fermi_Dirac_distribution(eps_v, mu_F, kbT)
    occ_c = Fermi_Dirac_distribution(eps_c, mu_F, kbT)
    
    zrho_dm2x2(1,1) = occ_v*abs(zeig_vec(1,1))**2 &
                     +occ_c*abs(zeig_vec(1,2))**2
    
    zrho_dm2x2(1,2) = occ_v*zeig_vec(1,1)*conjg(zeig_vec(2,1)) &
                     +occ_c*zeig_vec(1,2)*conjg(zeig_vec(2,2))
    zrho_dm2x2(2,1) = conjg(zrho_dm2x2(1,2))
                     

    zrho_dm2x2(2,2) = occ_v*abs(zeig_vec(2,1))**2 &
                     +occ_c*abs(zeig_vec(2,2))**2

    zrho_dm(:,:,ik) = 0d0
    zrho_dm(1:2,1:2,ik) =  zrho_dm2x2(1:2,1:2)

    occ_core = Fermi_Dirac_distribution(eps_core, mu_F, kbT)
    zrho_dm(3,3,ik) =  occ_core
    zrho_dm(4,4,ik) =  occ_core

  end do
    
end subroutine initialize_density_matrix

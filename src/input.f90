subroutine input
  use global_variables
  implicit none
  real(8) :: T1_relax_fs,T2_relax_fs
  real(8) :: mu_F_ev, kbT_K
  real(8) :: T_propagation_fs
  real(8) :: E0_1_Vm, omega_1_ev, tpulse_1_fs
  real(8) :: E0_2_Vm, omega_2_ev, tpulse_2_fs

  if(if_root_global)then
    write(*,"(A)")"Start: input"
  end if

  call read_basic_input('kx_max',kx_max,val_default = -1d0)
  call read_basic_input('ky_max',ky_max,val_default = -1d0)

  call read_basic_input('nkx',nkx,val_default = -1)
  call read_basic_input('nky',nky,val_default = -1)

!  call read_basic_input('T_relax_fs',T_relax_fs,val_default = -1d0)
!  T_relax = T_relax_fs*fs

  call read_basic_input('T1_relax_fs',T1_relax_fs,val_default = -1d0)
  T1_relax = T1_relax_fs*fs

  call read_basic_input('T2_relax_fs',T2_relax_fs,val_default = -1d0)
  T2_relax = T2_relax_fs*fs

  v_Fermi = clight*1.12d6/299792458d0
  eps_gap = 0d0*ev
  eps_core = -300*ev
  dip_core = 1d-2
  

  call read_basic_input('mu_F_ev',mu_F_ev,val_default = 0d0)
  mu_F = mu_F_ev*ev

  call read_basic_input('kbT_K',kbT_K,val_default = 0d0)
  kbT  = kbT_K/11604.505d0*ev


  call read_basic_input('E0_1_Vm',E0_1_Vm,val_default = 0d0)  
  call read_basic_input('E0_2_Vm',E0_2_Vm,val_default = 0d0)  

  call read_basic_input('omega_1_ev',omega_1_ev,val_default = 200d-3)  
  call read_basic_input('omega_2_ev',omega_2_ev,val_default = 300d0)  

  call read_basic_input('tpulse_1_fs',tpulse_1_fs,val_default = 1d0)  
  call read_basic_input('tpulse_2_fs',tpulse_2_fs,val_default = 1d0)  


  f0_1 = E0_1_Vm*ev/angstrom*1d-10 !V/m
  f0_2 = E0_2_Vm*ev/angstrom*1d-10 !V/m
  omega_1 = omega_1_ev*ev
  omega_2 = omega_2_ev*ev
  tpulse_1 = tpulse_1_fs*fs !2d3*fs
  tpulse_2 = tpulse_2_fs*fs




  call read_basic_input('T_propagation_fs',T_propagation_fs,val_default = 0d0)
  T_propagation = T_propagation_fs*fs
  call read_basic_input('dt',dt,val_default = 0d0)
  if(if_root_global)write(*,"(A,2x,e26.16e3)")'input   dt =',dt
!  nt_probe_period = aint( (2d0*pi/omega_2)/dt) +1
!  dt = (2d0*pi/omega_2)/nt_probe_period
!  if(if_root_global)write(*,"(A,2x,e26.16e3)")'refined dt =',dt


  nt = aint(T_propagation/dt) + 1

  if(if_root_global)then
    write(*,"(A)")"Finish: input"
  end if

end subroutine input

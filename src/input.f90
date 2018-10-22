subroutine input
  use global_variables
  implicit none
  real(8) :: T_relax_fs,mu_F_ev, kbT_K
  real(8) :: T_propagation_fs



  call read_basic_input('kx_max',kx_max,val_default = -1d0)
  call read_basic_input('ky_max',ky_max,val_default = -1d0)

  call read_basic_input('nkx',nkx,val_default = -1)
  call read_basic_input('nky',nky,val_default = -1)

  call read_basic_input('T_relax_fs',T_relax_fs,val_default = -1d0)

  T_relax = T_relax_fs*fs

  v_Fermi = clight*1.12d6/299792458d0
  eps_gap = 0d0*ev
  eps_core = -300*ev
  

  call read_basic_input('mu_F_ev',mu_F_ev,val_default = 0d0)
  mu_F = mu_F_ev*ev

  call read_basic_input('kbT_K',kbT_K,val_default = 0d0)
  kbT  = kbT_K/11604.505d0*ev


  call read_basic_input('T_propagation_fs',T_propagation_fs,val_default = 0d0)
  T_propagation = T_propagation_fs*fs
  call read_basic_input('dt',dt,val_default = 0d0)
  nt = aint(T_propagation/dt) + 1


end subroutine input

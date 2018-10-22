module global_variables
  use parallel
  use communication
  use math
  use constants


! quantum
  real(8) :: kx_max, ky_max
  real(8) :: dkx, dky

! Dirac band parameter
  real(8) :: v_Fermi
  real(8) :: eps_gap = 0d0

! core-level
  real(8) :: eps_core
  real(8) :: dip_core

! relaxation parameters
  real(8) :: T1_relax
  real(8) :: T2_relax

! Fermi-Dirac distribution
  real(8) :: mu_F
  real(8) :: kbT


! band parameter
  integer,parameter :: nband_type_dirac_cone = 0
  integer,parameter :: nband_type_massive_dirac = 1
  real(8),parameter :: tau_chiral = 1d0 ! chirality of Dirac band



end module global_variables

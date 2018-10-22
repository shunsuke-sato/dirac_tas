module global_variables
  use parallel
  use communication
  use math
  use constants
  use inputoutput

! quantum
  real(8) :: kx_max, ky_max
  real(8) :: dkx, dky
  integer :: nk, nkx, nky
  integer :: nk_start, nk_end
  integer :: nk_average, nk_remainder
  real(8),allocatable :: kx0(:),ky0(:)
  real(8),allocatable :: kx(:),ky(:)
  real(8) :: kx_shift, ky_shift
  integer,allocatable :: ik_table(:,:)

! Dirac band parameter
  real(8) :: v_Fermi
  real(8) :: eps_gap = 0d0

! core-level
  real(8) :: eps_core
  real(8) :: dip_core

! relaxation parameters
  real(8) :: T_relax

! Fermi-Dirac distribution
  real(8) :: mu_F
  real(8) :: kbT

! density matrix
  complex(8),allocatable :: zrho_dm(:,:,:)


! time propagation
  real(8) :: dt, T_propagation
  integer :: nt 

! laser 
  real(8),allocatable :: Act(:,:),Act_dt2(:,:)
  real(8),allocatable :: Ezt(:),Ezt_dt2(:)


! band parameter
  integer,parameter :: nband_type_dirac_cone = 0
  integer,parameter :: nband_type_massive_dirac = 1
  real(8),parameter :: tau_chiral = 1d0 ! chirality of Dirac band



end module global_variables

subroutine impulsive_field(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ik, iexp
  integer,parameter :: nexp = 16
  real(8) :: dip_mat(4,4),mat_tmp(4,4)
  complex(8) :: impulse_mat(4,4), zfact

  if(.not.(if_impulse .and. it==it_impulse))return
  dip_mat = 0d0
  dip_mat(3,1) = dip_core
  dip_mat(1,3) = dip_core
  dip_mat(4,2) = dip_core
  dip_mat(2,4) = dip_core

  impulse_mat = 0d0
  impulse_mat(1,1) = 1d0
  impulse_mat(2,2) = 1d0
  impulse_mat(3,3) = 1d0
  impulse_mat(4,4) = 1d0

  mat_tmp = real(impulse_mat)
  zfact = 1d0
  do iexp = 1, nexp
    mat_tmp = matmul(mat_tmp,dip_mat)
    zfact = zfact*(-zI*kmom_impulse)/iexp
    impulse_mat = impulse_mat + zfact*mat_tmp
  end do

  do ik = nk_start, nk_end
    
    zrho_dm(:,:,ik) = matmul(&
      matmul(impulse_mat,zrho_dm(:,:,ik)),&
      conjg(transpose(impulse_mat)))

  end do
  

end subroutine impulsive_field

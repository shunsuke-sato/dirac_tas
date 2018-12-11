subroutine calc_dipole(dip_k)
  use global_variables
  implicit none
  real(8), intent(out) :: dip_k
  real(8) :: dip_mat(4,4)
  complex(8) :: ztmp_mat(4,4)
  integer :: ik

  dip_mat = 0d0
  dip_mat(1,3) = -dip_core
  dip_mat(2,4) = -dip_core
  dip_mat(3,1) = -dip_core
  dip_mat(4,2) = -dip_core

  dip_k = 0d0
  do ik = nk_start, nk_end

    ztmp_mat = matmul(dip_mat,zrho_dm(:,:,ik))
    dip_k = dip_k + ztmp_mat(1,1) &
                  + ztmp_mat(2,2) &
                  + ztmp_mat(3,3) &
                  + ztmp_mat(4,4) 

  end do

  dip_k = dip_k*dkx*dky*2d0/(2d0*pi)**2


end subroutine calc_dipole

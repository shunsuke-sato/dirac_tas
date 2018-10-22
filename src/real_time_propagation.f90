subroutine real_time_propagation
  use global_variables
  implicit none
  integer :: it

  
  do it = 0, nt

    call dt_evolve(it)
    call calc_dipole

  end do



  


end subroutine real_time_propagation

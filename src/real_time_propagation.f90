subroutine real_time_propagation
  use global_variables
  implicit none
  integer :: it, it_t
  real(8) :: dipole_t(0:nt+1)

  dipole_t = 0d0
  call calc_dipole(dipole_t(0))

  
  do it = 0, nt
    call impulsive_field(it)
    call dt_evolve(it)
    call calc_dipole(dipole_t(it+1))
  end do


  call comm_allreduce(dipole_t)
  if(if_root_global)then
    open(20,file='dipole_t.out')
    write(20,"(A,2x,I9)")"nt=",nt
    do it_t = 0, nt
      write(20,"(999e26.16)")dt*it_t,Ezt(it_t),dipole_t(it_t)
    end do
    close(20)
  end if
  


end subroutine real_time_propagation

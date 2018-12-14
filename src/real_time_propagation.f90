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
    write(20,"(A,2x,I9)")"#nt=",nt
    write(20,"(A,2x,I9)")"#it_impulse=",it_impulse
    do it_t = 0, nt
      write(20,"(999e26.16)")dt*it_t,Ezt(it_t),dipole_t(it_t)
    end do
    close(20)
  end if
  
  call calc_polarizability_with_cw_probe(dipole_t)


end subroutine real_time_propagation
!-----------------------------------------------------------------------------------------
subroutine calc_polarizability_with_cw_probe(dipole_t)
  use global_variables
  implicit none
  real(8),intent(in) :: dipole_t(0:nt+1)
  integer :: it
  complex(8) :: zPw, zEw, zfact
  
  if(.not. if_root_global)return
  
  zPw = 0d0
  zEw = 0d0
  do it = nt-nt_probe_period+1,nt
    zfact = exp(zI*omega_2*dt*it)
    zPw = zPw + zfact*dipole_t(it)
    zEw = zEw + zfact*Ezt(it)
  end do
  zPw = zPw*dt
  zEw = zEw*dt
  
  write(*,"(A)")"Polarizability"
  write(*,"(A,2x,e26.16e3)")"omega (a.u.) =",omega_2
  write(*,"(A,2x,2e26.16e3)")"chi (a.u.) =",zPw/zEw
    


end subroutine calc_polarizability_with_cw_probe


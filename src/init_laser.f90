subroutine init_laser
  use global_variables
  implicit none
  real(8) :: f0_1, f0_2
  real(8) :: omega_1, omega_2
  real(8) :: tpulse_1, tpulse_2
  real(8) :: tdelay
  integer :: it
  real(8) :: tt, xx
  real(8) :: t_impulse_fs, t_impulse

  if(if_root_global)then
    write(*,"(A)")"Start: init_laser"
  end if

  f0_1 = 2d7*ev/angstrom*1d-10 !V/m
  f0_2 = 0d4*ev/angstrom*1d-10 !V/m
  omega_1 = 200d-3*ev
  omega_2 = 300d0*ev
  tpulse_1 = 1d0*fs !2d3*fs
  tpulse_2 = 1d0*fs


  call read_basic_input('T_impulse_fs',T_impulse_fs,val_default = 0d0)  
  t_impulse = t_impulse_fs*fs
  it_impulse = nint(t_impulse/dt)
  kmom_impulse = 1d-2
  

  allocate(Act(2,-1:nt+1),Act_dt2(2,-1:nt+1))
  act = 0d0; act_dt2 = 0d0
  allocate(Ezt(-1:nt+1),Ezt_dt2(-1:nt+1))
  Ezt = 0d0; Ezt_dt2 = 0d0

! pump
  if(.true.)then
    do it = 0, nt
      tt = dt*it
      xx = tt
      Act(1,it) = -f0_1/omega_1*sin(omega_1*xx)
      Act(2,it) = -f0_1/omega_1*cos(omega_1*xx)

      xx = tt + 0.5d0*dt
      Act_dt2(1,it) = -f0_1/omega_1*sin(omega_1*xx)
      Act_dt2(2,it) = -f0_1/omega_1*cos(omega_1*xx)
    end do

  else
  do it = 0, nt
    tt = dt*it
    xx = tt - 0.5d0*tpulse_1
    if(abs(xx) < 0.5d0*tpulse_1)then
      Act(1,it) = -f0_1/omega_1*cos(pi*xx/tpulse_1)**4*sin(omega_1*xx)
      Act(2,it) = -f0_1/omega_1*cos(pi*xx/tpulse_1)**4*cos(omega_1*xx)
    end if

    xx = tt - 0.5d0*tpulse_1 + 0.5d0*dt
    if(abs(xx) < 0.5d0*tpulse_1)then
      Act_dt2(1,it) = -f0_1/omega_1*cos(pi*xx/tpulse_1)**4*sin(omega_1*xx)
      Act_dt2(2,it) = -f0_1/omega_1*cos(pi*xx/tpulse_1)**4*cos(omega_1*xx)
    end if

  end do
  end if

! probe
  do it = 0, nt
    tt = dt*it
    xx = tt - 0.5d0*tpulse_1-tdelay
    if(abs(xx) < 0.5d0*tpulse_2)then
      Ezt(it) = f0_2*cos(pi*xx/tpulse_2)**4*sin(omega_2*xx)
    end if

    xx = tt - 0.5d0*tpulse_1-tdelay + 0.5d0*dt
    if(abs(xx) < 0.5d0*tpulse_2)then
      Ezt_dt2(it) = f0_2*cos(pi*xx/tpulse_2)**4*sin(omega_2*xx)
    end if

  end do

  if(if_root_global)then
    write(*,"(A)")"Finish: init_laser"
  end if

end subroutine init_laser

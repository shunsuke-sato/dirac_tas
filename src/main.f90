program main
  use global_variables
  implicit none

  call init_parallel
  call init_input

!  write(*,*)'Hello world!!'
  call input
  call initialize
  call init_laser

  call real_time_propagation

  call fin_input
  call fin_parallel

end program main

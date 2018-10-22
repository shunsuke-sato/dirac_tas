program main
  use global_variables
  implicit none

  call init_parallel
  call init_input

  write(*,*)'Hello world!!'
  call input
  call initialize


  call fin_input
  call fin_parallel

end program main

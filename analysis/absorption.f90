program main
  implicit none
  complex(8),parameter :: zi = (0d0, 1d0)
  integer,parameter :: nt = 103353
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: wi = (300d0-2d0)*ev,wf=(300d0+2d0)*ev
  real(8),parameter :: dw = (1d-3)*ev
  integer,parameter :: nw = aint((wf-wi)/dw)+1
  real(8) :: tt(0:nt), pz(0:nt),dt
  real(8) :: ww
  integer :: iw, it
  real(8) :: f1
  complex(8) :: zd,zs,zw

  open(20,file='dipole_t.out')
  read(20,*)
  do it = 0, nt
     read(20,*)tt(it),f1,pz(it)
  end do
  close(20)
  dt = tt(1)-tt(0)

  open(20,file="zchi.out")
  do iw = 0,nw
     ww = dw*iw + wi

     zs = 0d0
     zw = 1d0
     zd = exp(-zI*ww*dt)
     do it = 0,nt
        zs = zs + pz(it)*zw
        zw = zw*zd
     end do
     zs=zs*dt
     write(20,"(999e26.16e3)")ww,zs
     
  end do
  close(20)


end program main

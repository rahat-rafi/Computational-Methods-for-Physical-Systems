program henon_heiles
  implicit none
  integer, parameter :: n=4, ncases=3
  integer :: istep, nsteps, nout, icase
  double precision :: y(n), t, dt, tmax
  double precision :: T0(ncases), px0, py0
  character(len=30) :: fname

  ! Integration parameters
  dt    = 0.002d0
  tmax  = 200.d0
  nsteps = int(tmax/dt)
  nout   = int(0.1d0/dt)   ! write every 0.1 time units

  ! Array of kinetic energies
  T0 = (/ 0.09d0, 0.125d0, 0.1482d0 /)

  ! Loop over cases of kinetic energies
  do icase=1,ncases
     px0 = sqrt(T0(icase))
     py0 = px0

     ! initial conditions: x=0, y=0, px=px0, py=py0
     y = (/ 0.d0, 0.d0, px0, py0 /)
     t = 0.d0

     write(fname,'(A,F0.4,A)') '4_henon_0', T0(icase), '.dat'
     open(10,file=trim(adjustl(fname)),status='unknown')
     write(10,'(a)') "# t x y px py"

     do istep=1,nsteps
        call rk4n(t,dt,y,n)
        if (mod(istep,nout)==0) then
           write(10,'(5e16.8)') t,y(1),y(2),y(3),y(4)
        end if
     end do

     close(10)
  end do

end program henon_heiles


! ================= SUBROUTINES BELOW =================

subroutine derivs(t,y,dydt)
  implicit none
  double precision, intent(in) :: t, y(4)
  double precision, intent(out):: dydt(4)

  ! updated for this problem
  dydt(1) = y(3)
  dydt(2) = y(4)
  dydt(3) = -y(1) - 2.d0*y(1)*y(2)
  dydt(4) = -y(2) - y(1)**2 + y(2)**2
end subroutine derivs


subroutine rk4n(t,dt,y,n)
  implicit none
  integer, intent(in) :: n
  double precision, intent(inout) :: y(n), t, dt
  double precision :: k1(n), k2(n), k3(n), k4(n), yt(n)
  integer :: i

  call derivs(t,y,k1)
  do i=1,n
     yt(i) = y(i) + 0.5d0*dt*k1(i)
  end do
  call derivs(t+0.5d0*dt,yt,k2)

  do i=1,n
     yt(i) = y(i) + 0.5d0*dt*k2(i)
  end do
  call derivs(t+0.5d0*dt,yt,k3)

  do i=1,n
     yt(i) = y(i) + dt*k3(i)
  end do
  call derivs(t+dt,yt,k4)

  do i=1,n
     y(i) = y(i) + dt/6.d0*(k1(i)+2.d0*k2(i)+2.d0*k3(i)+k4(i))
  end do
  t = t + dt
end subroutine rk4n


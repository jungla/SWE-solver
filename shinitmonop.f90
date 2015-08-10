module shinit
  use kinds
  implicit none

  contains

!***********************************************************************
! Initialize Fields, u,v, and zt for the monopole vortex problem
!***********************************************************************
subroutine uvh_initial(u,v,zt, xe, ye, nx,ny, time)
  implicit none
  integer, intent(in) :: nx              ! number of grid cells in x
  integer, intent(in) :: ny              ! number of grid cells in y
  real(r8), intent(in) :: xe(nx+1)       ! zonal edges
  real(r8), intent(in) :: ye(ny+1)       ! meriodional edges
  real(r8), intent(in), optional :: time ! timelevel
  real(r8), intent(out) :: u(nx+1,ny)    ! zonal velocity
  real(r8), intent(out) :: v(nx,ny+1)    ! meridional velocity
  real(r8), intent(out) :: zt(nx,ny)     ! pressure field

  integer :: i,j,nx1,ny1
  real(r8) :: xu,xv,xp,x
  real(r8) :: yu,yv,yp,y
  real(r8) :: r, r2, Vz, timei

  real(r8), parameter :: f=9.0e-5_r8       ! central Coriolis parameter(1/s)
  real(r8), parameter :: gravity=0.0810_r8 ! reduced gravity (m/s)
  real(r8), parameter :: xcent=900.0e3_r8  ! x-center of monopole (m)
  real(r8), parameter :: ycent=1400.0e3_r8 ! y-center of monopole (m)
  real(r8), parameter :: A=39.5061_r8      ! monopole pressure amplitude (m)
  real(r8), parameter :: L=200.0e3_r8      ! monopole decay scale (in km)
  real(r8), parameter :: pin=3.1415926535897932384_r8
  real(r8), parameter :: L2=L**2
  real(r8), parameter :: vamp=8._r8* gravity*A / ((f*L)**2)
  real(r8), parameter :: h_bar = 0.492508228253105_r8

! Q = (derf(27.d0/2.d0)+derf(9.d0/2.d0))*derf(7.d0)
! Q = Q*A*0.5_r8*L2*pin
! h_bar = Q/(3600.e3_r8*2800.e3_r8)
  write(6,*) 'h bar =',h_bar

  nx1 = nx+1
  ny1 = ny+1

  if (.not.present(time)) then
    timei = 0._r8
  else
    timei = time
  endif

  do j = 1,ny
    yu = 0.5_r8* (ye(j) + ye(j+1));   ! meridional center of cell j
    y  = yu - ycent;                  ! meriodonal distance from monopole
    do i = 1,nx1
      xu = xe(i);                     ! west edge of cell i
      x  = xu - xcent                 ! zonal distance from monopole
      r2 = x**2 + y**2
      r  = SQRT(r2)                   ! radial distance
      Vz = 0.5_r8*f*r* ( 1._r8 - SQRT( 1._r8 - vamp*EXP(-r2/L2)) )
      if (r > 0._r8) then
       u(i,j) = y*Vz/r
      else
       u(i,j) = 0._r8
      endif
    enddo
  enddo

  do j = 1,ny1
    yv = ye(j);                       ! south edge of cell j
    y  = yv - ycent;                  ! meriodonal distance from monopole
    do i = 1,nx
      xv = 0.5_r8*(xe(i)+xe(i+1))     ! zonal center of cell i
      x  = xv - xcent                 ! zonal distance from monopole
      r2 = x**2 + y**2
      r  = SQRT(r2)                   ! radial distance
      Vz = 0.5_r8*f*r* ( 1._r8 - SQRT( 1._r8 - vamp*EXP(-r2/L2)) )
      if ( r.ne.0._r8) then
       v(i,j) =-x*Vz/r
      else
       v(i,j) = 0.0_r8
      endif
    enddo
  enddo

  do j = 1,ny
    yp= 0.5_r8*(ye(j)+ye(j+1))
    y = yp - ycent
    do i = 1,nx
      xp= 0.5_r8*(xe(i)+xe(i+1))
      x = xp - xcent
      r2 = x**2 + y**2
      zt(i,j) = A*EXP(-r2/L2) - h_bar
    enddo
  enddo

  return
  end subroutine uvh_initial
end module shinit

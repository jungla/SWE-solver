module shinit
  use kinds
  use shparams, only : depth,gravity
  implicit none
  contains
!***********************************************************************
!***********************************************************************
subroutine uvh_initial(u,v,zt, xe,ye, Nx,Ny, time)
  implicit none
  integer, intent(in) :: Nx,Ny
  real(r8), intent(in) :: xe(Nx+1)
  real(r8), intent(in) :: ye(Ny+1)
  real(r8), intent(in), optional :: time ! timelevel
  real(r8), intent(out) :: u(Nx+1,Ny)
  real(r8), intent(out) :: v(Nx,Ny+1)
  real(r8), intent(out) :: zt(Nx,Ny)

  real(r8), parameter :: pin = 3.1415926535897932384_r8
  real(r8), parameter :: amp=1.0_r8 ! pressure amplitude at t=0
  real(r8) :: uamp     ! amplitude coeff of x-velocity
  real(r8) :: vamp     ! amplitude coeff of y-velocity
  real(r8) :: zamp     ! amplitude coeff of pressure
  real(r8) :: sigma    ! frequency of oscillations
  real(r8), parameter :: kx=2._r8 ! wave number in x
  real(r8), parameter :: ky=2._r8 ! wave number in y
  real(r8) :: k, C
  real(r8) :: xp,yp,timei
  integer :: i,j

  if (.not.present(time)) then
    timei = 0._r8
  else
    timei = time
  endif

  C = sqrt(depth*gravity)          ! gravity wave speed
  k = sqrt(kx*kx+ky*ky)            ! total wave number
  sigma = C * pin*k                ! frequency in a unit box
  uamp = amp * (gravity*kx/sigma) * sin(sigma*timei) ! x-velocity factor
  vamp = amp * (gravity*ky/sigma) * sin(sigma*timei) ! y-velocity factor
  zamp = amp * cos(sigma*timei)                      ! pressure factor

  do j = 1,Ny
    yp = 0.5_r8* (ye(j)+ye(j+1))
    do i = 1,Nx+1
      u(i,j)  = uamp * sin(pin*xe(i)) * cos(pin*yp)
    enddo
    do i = 1,Nx
      xp = 0.5_r8* (xe(i)+xe(i+1))
      zt(i,j) = amp  * cos(pin*xp) * cos(pin*yp)
    enddo
  enddo

  do j = 1,Ny+1
    do i = 1,Nx
      xp = 0.5_r8* (xe(i)+xe(i+1))
      v(i,j) = vamp * cos(pin*xp) * sin(pin*ye(j))
    enddo
  enddo

  return
end subroutine uvh_initial

end module shinit

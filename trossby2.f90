module shinit
  use kinds
  implicit none

  contains

!***********************************************************************
! Initialize Fields, u,v, and zt for the equatorial Rossby
! soliton. The IC is obtained from the asymptotic solution of
! J.P. Boyd.
!***********************************************************************
subroutine uvh_initial(u,v,zt, xe, ye, nx,ny, time)
  implicit none
  integer, intent(in) :: nx	         ! number of grid cells in x
  integer, intent(in) :: ny	         ! number of grid cells in y
  real(r8), intent(in) :: xe(nx+1)       ! zonal edges
  real(r8), intent(in) :: ye(ny+1)       ! meriodional edges
  real(r8), intent(in), optional :: time ! timelevel
  real(r8), intent(out) :: u(nx+1,ny)    ! zonal velocity
  real(r8), intent(out) :: v(nx,ny+1)    ! meridional velocity
  real(r8), intent(out) :: zt(nx,ny)     ! pressure field

  integer :: i,j, nx1,ny1
  real(r8), parameter :: B=0.5_r8
  real(r8), parameter :: A = 0.771_r8 * (B**2)
  real(r8), parameter :: c0 = -1.0_r8/3.0_r8
  real(r8), parameter :: c1 = -0.3950_r8*(B**2)
  real(r8), parameter :: cc = c0 + c1

  real(r8) :: xi,ey2,exi,eta,eta_xi
  real(r8) :: u0,v0,z0,xg,ulat1 
  real(r8) :: u1,v1,z1,yg,vlat1 
  real(r8) :: zlat1, timei

  nx1 = nx+1
  ny1 = ny+1

  if (.not.present(time)) then
    timei = 0._r8
  else
    timei = time
  endif

  do j = 1,ny
    yg = 0.5_r8* (ye(j) + ye(j+1));                 ! (j-0.5_r8)*dy + ymin
    ey2 = exp( -0.5_r8 * (yg**2) )	            ! e^{-y^2/2}
    ulat1 = first_order(yg, 1)
    do i = 1,nx1
      xg = xe(i);                                   ! (i-1)*dx + xmin
      xi = xg - cc*timei		            ! x - c t
      exi = exp( -B*xi )			    ! e^{-B\xi}
      eta = A* ( (2.0_r8*exi/(1.0_r8+exi**2)) **2 ) ! A \sech^2(B\xi)
      eta_xi =-2.0_r8* B* tanh(B*xi) * eta	    ! \eta_{\xi}
      u0 = 0.25_r8*eta* (6.0_r8*(yg**2)-9.0_r8) *ey2
      u1 = 0.5625_r8 * c1 * eta * (2.0_r8*(yg**2)+3.0_r8) * ey2 &
         + (eta**2) * ulat1
      u(i,j) = u0 + u1
    enddo
  enddo
  do j = 1,ny1
    yg = ye(j)                                      ! (j-1)*dy + ymin
    ey2 = exp( -0.5_r8 * (yg**2) )                  ! e^{-y^2/2}
    vlat1 = first_order(yg, 2)
    do i = 1,nx
      xg = 0.5_r8* (xe(i)+xe(i+1));                 ! (i-0.5_r8)*dx + xmin
      xi = xg - cc*timei		            ! x - c t
      exi = exp( -B*xi )			    ! e^{-B\xi}
      eta = A* ( (2.0_r8*exi/(1.0_r8+exi**2)) **2 ) ! A \sech^2(B\xi)
      eta_xi =-2.0_r8* B* tanh(B*xi) * eta	    ! \eta_{\xi}
      v0 = eta_xi * 2.0_r8*yg * ey2
      v1 = eta * eta_xi * vlat1
      v(i,j) = v0 + v1
    enddo
  enddo
  do j = 1,ny
    yg = 0.5_r8* (ye(j) + ye(j+1));                 ! (j-0.5_r8)*dy + ymin
      ey2 = exp( -0.50_r8* (yg**2) )	            ! e^{-y^2/2}
    zlat1 = first_order(yg, 3)
    do i = 1,nx
      xg = 0.5_r8* (xe(i)+xe(i+1));                 ! (i-0.5_r8)*dx + xmin
      xi = xg - cc*timei		            ! x - c t
      exi = exp( -B*xi )			    ! e^{-B\xi}
      eta = A* ( (2.0_r8*exi/(1.0_r8+exi**2)) **2 ) ! A \sech^2(B\xi)
      eta_xi =-2.0_r8* B* tanh(B*xi) * eta	    ! \eta_{\xi}
      z0 = 0.25_r8*eta* (6.0_r8*(yg**2)+3.0_r8) *ey2
      z1 = 0.5625_r8 * c1 * eta * (2.0_r8*(yg**2)-5.0_r8) * ey2 &
         + (eta**2) * zlat1
      zt(i,j) = z0 + z1
    enddo
  enddo

  return
end subroutine uvh_initial

!***********************************************************************
!		First Order Correction to Soliton Solution
! Compute latitudinal structure of first order correction.
!***********************************************************************
real(r8) function first_order(yv, flag)
  implicit none
  integer, intent(in) :: flag
  real(r8), intent(in) :: yv

  integer, parameter :: nmax=26
  real( r8 ) :: H(0:nmax)
  real( r8 ) :: sum
  integer :: n

  real(r8) :: un(0:nmax) = (/                                &
     1.789276e+00_r8, 0.0000000_r8, 0.1164146e+00_r8, 0.0000000_r8,&
   -0.3266961e-03_r8, 0.0000000_r8,-0.1274022e-02_r8, 0.0000000_r8,&
    0.4762876e-04_r8, 0.0000000_r8,-0.1120652e-05_r8, 0.0000000_r8,&
    0.1996333e-07_r8, 0.0000000_r8,-0.2891698e-09_r8, 0.0000000_r8,&
    0.3543594e-11_r8, 0.0000000_r8,-0.3770130e-13_r8, 0.0000000_r8,&
    0.3547600e-15_r8, 0.0000000_r8,-0.2994113e-17_r8, 0.0000000_r8,&
    0.2291658e-19_r8, 0.0000000_r8,-0.1178252e-21_r8/)
  real(r8) :: vn(0:nmax) = (/                                &
    0.0000000_r8, 0.0000000e+00_r8, 0.0000000_r8,-0.6697824e-01_r8,&
    0.0000000_r8,-0.2266569e-02_r8, 0.0000000_r8, 0.9228703e-04_r8,&
    0.0000000_r8,-0.1954691e-05_r8, 0.0000000_r8, 0.2925271e-07_r8,&
    0.0000000_r8,-0.3332983e-09_r8, 0.0000000_r8, 0.2916586e-11_r8,&
    0.0000000_r8,-0.1824357e-13_r8, 0.0000000_r8, 0.4920950e-16_r8,&
    0.0000000_r8, 0.6302640e-18_r8, 0.0000000_r8,-0.1289167e-19_r8,&
    0.0000000_r8, 0.1471189e-21_r8, 0.0000000_r8/)
  real(r8) :: zn(0:nmax) = (/                                &
    -3.071430e+00_r8, 0.0000000_r8,-0.3508384e-01_r8, 0.0000000_r8,&
   -0.1861060e-01_r8, 0.0000000_r8,-0.2496364e-03_r8, 0.0000000_r8,&
    0.1639537e-04_r8, 0.0000000_r8,-0.4410177e-06_r8, 0.0000000_r8,&
    0.8354759e-09_r8, 0.0000000_r8,-0.1254222e-09_r8, 0.0000000_r8,&
    0.1573519e-11_r8, 0.0000000_r8,-0.1702300e-13_r8, 0.0000000_r8,&
    0.1621976e-15_r8, 0.0000000_r8,-0.1382304e-17_r8, 0.0000000_r8,&
    0.1066277e-19_r8, 0.0000000_r8,-0.1178252e-21_r8/)

  call hermite(H,yv,nmax)
  sum = 0.0_r8
  select case (flag)
    case (1)
      do n = 0,nmax
        sum = sum + un(n) * H(n)
      enddo
    case (2)
      do n = 0,nmax
        sum = sum + vn(n) * H(n)
      enddo
    case (3)
      do n = 0,nmax
        sum = sum + zn(n) * H(n)
      enddo
  end select
  first_order = sum

  return
end function first_order
!***********************************************************************
!		Hermite Polynomials
!  Compute Hermite polynomials from order 0 to n for the argument x.
!***********************************************************************
subroutine hermite(H,x,n)
  implicit none
  integer, intent( in ) :: n
  real(r8), intent( in ) :: x

  real(r8), intent( out ) ::  H(0:n)

  integer :: i
  real*8 ::  two_x, ex2
  real(r8), parameter :: two=2.0_r8, one=1.0_r8

  ex2 = exp(-x*x/two)
  two_x = two * x
  H(0) = ex2
  if (n > 0) then
    H(1) = two_x * ex2
    do i = 2,n
      H(i) = two_x * H(i-1) - (two*i-two) * H(i-2)
    enddo
  endif

  return
end subroutine hermite

end module shinit

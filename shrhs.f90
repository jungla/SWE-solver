module shrhs

 use kinds
 use cgridops
 use shparams, only : gravity,depth,bdrag,viscosity,dx,dy,Nx,Ny,Nx1,Ny1
 implicit none

 real(r8), private :: fcorio(Nx1,Ny1)    ! coriolis parameter
 real(r8), private :: F(Nx,Ny)    ! coriolis parameter
 real(r8), private :: wsx(Ny)           ! x-wind stress
 real(r8) :: hxy(Nx1,Ny1),u2(Nx1,Ny), v2(Nx,Ny1), q(Nx1,Ny1)

 contains

!**********************************************************************
! RHS of shallow water equations
!**********************************************************************
subroutine swerhs(u,v,zt,ru,rv,rzt,h,dx,dy)
  implicit none
  real(r8), intent(in) :: u(Nx1,Ny)
  real(r8), intent(in) :: v(Nx,Ny1)
  real(r8), intent(in) :: zt(Nx,Ny)
  real(r8), intent(in) :: dx, dy, h(Nx,Ny) ! function values needed to compute rhs
  real(r8), intent(out) :: ru(Nx1,Ny)
  real(r8), intent(out) :: rv(Nx,Ny1)
  real(r8), intent(out) :: rzt(Nx,Ny)

  real(r8) :: a(2), hy(Nx,Ny1), hx(Nx1,Ny), hz(Nx,Ny) ! function values needed to compute rhs
  real(r8) :: vx(Nx1,Ny1), uy(Nx1,Ny1), u2x(Nx,Ny), v2y(Nx,Ny)
  real(r8) :: vhx(Nx1,Ny1), uhy(Nx1,Ny1)
  real(r8) :: vv(Nx,Ny1), uv(Nx1,Ny) ! viscosity
 
  integer :: M(2)


!!! first VELOCITIES as a function of the pressure gradient

 hz = h + zt

 u2 = u**2
 v2 = v**2

! define mass fluxes hx = u*hx, uy = v*hy

 a(1) = 0.5_r8
 a(2) = 0.5_r8
 
 M(1) = Nx
 M(2) = Ny
  
 call xOP_2D(hx, hz, 0._r8, a, Nx1, Nx, M, 0) ! from p to u-points: 0 
 hx(1,:) = 1.0_r8
 hx(Nx1,:) = 1.0_r8

 call yOP_2D(hy, hz, 0._r8, a, Nx, Ny1, Nx, Ny, M, 0) 
 hy(:,1) = 1.0_r8
 hy(:,Ny1) = 1.0_r8

! define hxy
 M(1) = Nx1
 M(2) = Ny

 call yOP_2D(hxy, hx, 0._r8, a, Nx1, Ny1, Nx1, Ny, M, 0) ! from u points to p-point...  
 hxy(:,1) = 1._r8
 hxy(:,Ny1)= 1._r8
 hxy(1,:) = 1._r8
 hxy(Nx1,:)= 1._r8

 hx = u*hx !hx mass flux
 hy = v*hy !hy mass flux

! compute PHI

 M(1) = Nx
 M(2) = Ny
 
 call yOP_2D(v2y, v2, 0._r8, a, Nx, Ny, Nx, Ny1, M, 1)  
 call xOP_2D(u2x, u2, 0._r8, a, Nx, Nx1, M, 1) 

 F = gravity*zt + v2y * 0.5_r8 + u2x * 0.5_r8

! compute q, relative vorticity

 M(1) = Nx
 M(2) = Ny1

 a(1) =  1._r8/dx
 a(2) = -1._r8/dx
 call xOP_2D(vx, v, 0._r8, a, Nx1, Nx, M, 0) ! from p to u-points: 0 
 vx(1,:) = 0._r8
 vx(Nx1,:) = 0._r8
 vx(:,1) = 0._r8
 vx(:,Ny1) = 0._r8

 M(1) = Nx1
 M(2) = Ny
 
 a(1) =  1._r8/dy
 a(2) = -1._r8/dy
 call yOP_2D(uy, u, 0._r8, a, Nx1, Ny1, Nx1, Ny, M, 0) ! from p tp u-points: 0
 uy(:,1) = 0._r8
 uy(:,Ny1) = 0._r8
 uy(1,:) = 0._r8
 uy(Nx1,:) = 0._r8

 q = (vx - uy + fcorio)/hxy ! add coriolis force too...

! RHS Momentum

! Fx and Fy (assign to ru and rv) ... I will add the vorticity component later
 M(1) = Nx
 M(2) = Ny

 a(1) = -1_r8/dx
 a(2) =  1_r8/dx
 call xOP_2D(ru, F, 0._r8, a, Nx1, Nx, M, 0) ! from p to u-points: 0 
 ru(Nx1,:) = 0._r8
 ru(1,:) = 0._r8

 a(1) = -1_r8/dy
 a(2) =  1_r8/dy
 call yOP_2D(rv, F, 0._r8, a, Nx, Ny1, Nx, Ny, M, 0) ! from p tp u-points: 0
 rv(:,1) = 0._r8
 rv(:,Ny1) = 0._r8


! rhs for ru and rv

! vhx uhy, averages of V and U in x and y 

 M(1) = Nx
 M(2) = Ny1
 a(1) = 0.5_r8
 a(2) = 0.5_r8
 call xOP_2D(vhx, hy, 0._r8, a, Nx1, Nx, M, 0) ! from p to u-points: 0 
 vhx(Nx1,:) = 0._r8
 vhx(1,:) = 0._r8
 vhx(:,Ny1) = 0._r8
 vhx(:,1) = 0._r8
 
 M(1) = Nx1
 M(2) = Ny
 a(1) = 0.5_r8
 a(2) = 0.5_r8
 call yOP_2D(uhy, hx, 0._r8, a, Nx1, Ny1, Nx1, Ny, M, 0) ! from p tp u-points: 0
 vhx(:,Ny1) = 0._r8
 vhx(:,1) = 0._r8
 vhx(Nx1,:) = 0._r8
 vhx(1,:) = 0._r8

 vhx = vhx*q
 uhy = uhy*q

! rv ... adding to the already present value of rv and ru
 M(1) = Nx  ! values used up to the end of the script
 M(2) = Ny

 a(1) = -0.5_r8
 a(2) = -0.5_r8
 call xOP_2D(rv, uhy, 1._r8, a, Nx, Nx1, M, 1) ! from p to u-points: 0 

!! viscosity for v
! M(1) = Nx
! M(2) = Ny1
! a(1) =  -1_r8/dy**2
! a(2) =  1_r8/dy**2
! call y2OP_2D(vv, v, 0._r8, a, Nx, Ny1, M)
! vv(:,1) = 0._r8 
! vv(:,Ny1) = 0._r8 

! viscosity for u
! M(1) = Nx1
! M(2) = Ny
! a(1) =  -1_r8/dx**2
! a(2) =  1_r8/dx**2
! call x2OP_2D(uv, u, 0._r8, a, Nx1, Nx1, M)
! uv(1,:) = 0._r8
! uv(Nx1,:) = 0._r8

 !ru = ru + viscosity*uv
 !rv = rv + viscosity*vv

! ru
 M(1) = Nx  ! values used up to the end of the script
 M(2) = Ny
 a(1) = +0.5_r8
 a(2) = +0.5_r8
 call yOP_2D(ru, vhx, 1._r8, a, Nx1, Ny, Nx1, Ny1, M, 1) ! from p tp u-points: 0

! RHS Continuity

 a(1) = -1_r8/dx
 a(2) =  1_r8/dx
 call xOP_2D(rzt, hx, 0._r8, a, Nx, Nx1, M, 1) ! div, from p to u points.
 
 a(1) = -1_r8/dy
 a(2) =  1_r8/dy
 call yOP_2D(rzt, hy, 1._r8, a, Nx, Ny, Nx, Ny1, M, 1) ! div, from p to u points.
!rzt =-(hx(2:Nx1,1:Ny)-hx(1:Nx,1:Ny))/dx-(hy(1:Nx,2:Ny1)-hy(1:Nx,1:Ny))/dy

 return
end subroutine swerhs

!***********************************************************************
! Perform some diagnostics calculations
!***********************************************************************
subroutine Diagnostics(u,v,zt,h,volume,pvbudget,TE,pv2budget)
  implicit none
  real(r8), intent(out) :: volume        ! volume budget
  real(r8), intent(out) :: pvbudget      ! vorticity budget
  real(r8), intent(out) :: TE            ! total energy budget
  real(r8), intent(out) :: pv2budget     ! enstropy budget
  real(r8), intent(in)  ::  u(Nx1,Ny )
  real(r8), intent(in)  ::  v(Nx ,Ny1)
  real(r8), intent(in)  :: zt(Nx ,Ny )
  real(r8), intent(in)  ::  h(Nx ,Ny )
  integer :: i,j,M(2)
  real(r8):: a(2),PE, KE(Nx,Ny)

  volume = 0._r8   ! Initialize Volume
  PE = 0._r8       ! potential energy
  KE = 0._r8       ! kinetic energy
  TE = 0._r8

! KE
  a(1) = 0.5_r8
  a(2) = 0.5_r8
  M(1) = Nx  ! values used up to the end of the script
  M(2) = Ny
  call xOP_2D(KE, u2, 0._r8, a, Nx, Nx1, M, 1) ! from u to p-points: 1

  M(1) = Nx  ! values used up to the end of the script
  M(2) = Ny
  call yOP_2D(KE, v2, 1._r8, a, Nx, Ny, Nx, Ny1, M, 1) ! from p tp u-points: 0

  KE = 0.5_r8*( (h + zt)*KE + gravity*zt**2)

!  TE        = PE+KE      ! total energy
  pvbudget = 0._r8      ! potential enstrophy 
  pv2budget = 0._r8      ! potential enstrophy 

  do j = 1,Ny
    do i = 1,Nx
      volume = volume + zt(i,j)
      TE = TE + KE(i,j)
      pv2budget = pv2budget + (q(i,j)**2)*0.5_r8*hxy(i,j)
      pvbudget = pvbudget + q(i,j)
    enddo
  enddo

  do j = 1,Ny1
   pvbudget = pvbudget + q(Nx1,j)
   pv2budget = pv2budget + q(Nx1,j)
  enddo 

  do i = 1,Nx1
   pvbudget = pvbudget + q(i,Ny1)**2*0.5_r8*hxy(i,Ny1) 
   pv2budget = pv2budget + q(i,Ny1)**2*0.5_r8*hxy(i,Ny1)
  enddo 

  TE = TE*dx*dy
  volume = volume*dx*dy


 
  
  return
end subroutine Diagnostics

!***********************************************************************
! Initialize the Coriolis force on a beta plane
!***********************************************************************
subroutine Set_Coriolis(ye,fc,beta)
  implicit none
  real(r8), intent(in) :: fc
  real(r8), intent(in) :: beta
  real(r8), intent(in) :: ye(Ny1)

  integer :: i,j
  real(r8) :: yc, yv, fco

  yc = 0.5_r8*(ye(Ny1)+ye(1))           ! center of domain
  do j = 1,Ny1                           ! why Ny1?? don't I want it defined on the the p-points??
    fco = fc + beta * (ye(j)-yc)
    do i = 1,Nx1                         ! same here... why Nx1?
      fcorio(i,j) = fco
    enddo
  enddo

  return
end subroutine Set_Coriolis
!***********************************************************************
subroutine Set_Wind(xe,ye)
  implicit none
  real(r8), intent(in) :: xe(:), ye(:)

!...Double Gyre Wind
! real(r8), parameter ::wmax=3.0e-5_r8     ! max-wind stress
! integer :: j
! real(r8) :: yv,pi,wk
! pi = 2._r8*asin(1._r8)
! wk = 2._r8*pi/(ye(Ny1)-ye(1))
! do j = 1,Ny
!   yv = 0.5_r8*(ye(j)+ye(j+1))
!   wsx(j) =-wmax*cos(wk*yv)
! enddo

  wsx = 0._r8

  return
end subroutine Set_Wind

end module shrhs

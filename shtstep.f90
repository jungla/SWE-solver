module shtstep

  use kinds
  use shparams, only : Nx,Ny  ! grid-dimension
  use shrhs   , only : swerhs ! to calculate time-tendencies
  implicit none

!
!    This module performs time-stepping on the shallow water
!    equations. The user can choose between RK4, RK3 and AB3,
!    or Leap-Frog Trapezoidal by calling RK4Step, RK3Step, 
!    AB3Step, or LFTStep. They all take the same input 
!    and return the same output.
!    Usage is RRRStep(u,v,zt,dt) where RRR is either RK4, RK3,
!    AB3, or LFT.
!    The past right hand sides for AB3 are stored within static
!    variables local to that subroutine. This is not the cleanest
!    of programming constructs but it is the easiest to
!    understand, and does not involve lugging variables around
!    the entire program. An alternative is to define the data
!    type below, AB3data, and pass it to AB3.
!   
!    Each subroutine does a single time-step, and update the
!    3 dynamic variables: 2 velocity components and pressure.
!    The schemes are explicit in time and the time-tendencies
!    are calculated in subroutine swerhs.
!    

  integer, parameter, private :: nab3=3
!
!...The data type below is necessary for AB3 time-stepping only
!---------------  Code that can work with GFORTRAN ------------------
! type AB3data
!   private
!   real(r8) ::  ru(Nx+1,Ny  ,nab3)
!   real(r8) ::  rv(Nx  ,Ny+1,nab3)
!   real(r8) :: rzt(Nx  ,Ny  ,nab3)
!   integer :: imod=-2
! end type AB3data
!---------------  Code that can work with GFORTRAN ------------------

!---------------  Code that can't work with GFORTRAN ------------------
! GFORTRAN cannot allocate components of a derived type-
! type AB3data
!   private
!   real(r8), allocatable ::  ru(:,:,:)
!   real(r8), allocatable ::  rv(:,:,:)
!   real(r8), allocatable :: rzt(:,:,:)
! end type AB3data
!---------------  Code that can't work with GFORTRAN ------------------
!

  contains

!**********************************************************************
!	Runge-Kutta Integration of Shallow Water Equations
!**********************************************************************
subroutine RK4Step(u,v,zt,dt,h,dx,dy)
  implicit none
! integer, intent(in) :: Nx,Ny                ! number of cells
  real(r8), intent(in) :: dt,dx,dy                  ! time step
  real(r8), intent(inout) ::  u(Nx+1,Ny  )    ! u-velocity
  real(r8), intent(inout) ::  v(Nx  ,Ny+1)    ! v-velocity
  real(r8), intent(inout) :: zt(Nx  ,Ny  )    ! pressure
  real(r8), intent(inout) :: h(Nx  ,Ny  )    ! pressure

  real(r8), parameter :: factor=1.0_r8/6.0_r8
  real(r8) :: dts
  real(r8) :: qu(Nx+1,Ny), qv(Nx,Ny+1), qzt(Nx,Ny) ! temp solution
  real(r8) :: su(Nx+1,Ny), sv(Nx,Ny+1), szt(Nx,Ny) ! acc slopes
  real(r8) :: ru(Nx+1,Ny), rv(Nx,Ny+1), rzt(Nx,Ny) ! stage slope
  integer  :: is

!             Stage 1 slope at "time" for forward Euler half step 
  dts = 0.5_r8*dt
  call swerhs(u,v,zt,su,sv,szt,h,dx,dy)
  call Update(qu,qv,qzt, u,v,zt, su,sv,szt, dts) ! stage 1 u,v,zt

  do is = 2,3
!.....slope estimates at "time+0.5*timestep"
    call swerhs(qu,qv,qzt, ru,rv,rzt, h,dx,dy)
    select case (is)
      case (2)
        dts = 0.5*dt ! Backward Euler half step Stage 2
      case (3)
        dts = dt     ! Midpoint rule full step Stage 3
    end select
    su =  su + 2.0_r8 *ru
    sv =  sv + 2.0_r8 *rv
    szt= szt + 2.0_r8 *rzt
    call Update(qu,qv,qzt, u,v,zt, ru,rv,rzt, dts) ! stage 2,3 u,v,zt
  enddo

!         Stage 4 slope at "time+timestep" for Simpson rule
  call swerhs(qu,qv,qzt,ru,rv,rzt,h,dx,dy)
  su  =  su + ru
  sv  =  sv + rv
  szt = szt + rzt
  dts= factor*dt
  call InplaceUpdate(u,v,zt, su,sv,szt, dts) ! final update  u,v,zt

  return
end subroutine RK4Step

!**********************************************************************
! Time integration is performed with a 3rd order Runge-Kutta scheme.
!**********************************************************************
subroutine RK3Step(u,v,zt,dt,h,dx,dy)
  implicit  none
! integer, intent(in)     :: Nx,Ny    ! number of cells
  real(r8), intent(in)    :: dt,dx,dy       ! time step size
  real(r8), intent(inout) ::  u(Nx+1,Ny  ) ! x-velocity
  real(r8), intent(inout) ::  v(Nx  ,Ny+1) ! y-velocity
  real(r8), intent(inout) ::  h(Nx  ,Ny) ! y-velocity
  real(r8), intent(inout) :: zt(Nx  ,Ny  ) ! pressure

  real(r8) :: a1,a2
  integer :: nst
  integer, parameter :: nstage=3
  real(r8), parameter :: rk3c(nstage)=(/0._r8,0.75_r8,1.0_r8/3.0_r8/)
  real(r8) :: qu(Nx+1,Ny),qv(Nx,Ny+1),qzt(Nx,Ny) ! store temporary values
  real(r8) :: ru(Nx+1,Ny),rv(Nx,Ny+1),rzt(Nx,Ny) ! store rhs

!                 stage 1 (time_rhs = time)
  call swerhs( u, v, zt,ru,rv,rzt,h,dx,dy)        ! Right hand side
  call Update(qu,qv,qzt, u,v,zt, ru,rv,rzt, dt)

!                 stage 2 (time_rhs = time + dt)
  nst = 2
  call swerhs(qu,qv,qzt,ru,rv,rzt,h,dx,dy)     ! Right hand side
  call InPlaceUpdate(qu,qv,qzt, ru,rv,rzt, dt)
  a1 = rk3c(nst);
  a2 = 1._r8-rk3c(nst)
  qu  =  a1*u + a2* qu     ! x-velocity estimate 2
  qv  =  a1*v + a2* qv     ! y-velocity estimate 2
  qzt = a1*zt + a2*qzt     ! pressure estimate 2

!                 stage 3 (time_rhs = time + 0.5_r8*dt)
  nst = 3
  call swerhs(qu,qv,qzt,ru,rv,rzt,h,dx,dy)     ! Right Hand Side
  call InPlaceUpdate(qu,qv,qzt, ru,rv,rzt, dt)
  a1 = rk3c(nst)
  a2 = 1._r8-rk3c(nst)
  u  =  a1*u + a2* qu       ! final x-velocity
  v  =  a1*v + a2* qv       ! final y-velocity
  zt = a1*zt + a2*qzt       ! final pressure

  return
end subroutine RK3Step
!**********************************************************************
! AB3 time-stepping
!**********************************************************************
subroutine AB3Step(u,v,zt,dt,h,dx,dy)
  implicit none
! type(AB3data), intent(inout) :: vd
! integer, intent(in) :: Nx,Ny
  real(r8),intent(in) :: dt,dx,dy                 ! time step
  real(r8),intent(inout) ::  u(Nx+1,Ny  )   ! dependent variables
  real(r8),intent(inout) ::  v(Nx  ,Ny+1)   ! dependent variables
  real(r8),intent(inout) ::  h(Nx  ,Ny)   ! dependent variables
  real(r8),intent(inout) :: zt(Nx  ,Ny  )   ! dependent variables

  real(r8), save ::  ru(Nx+1,Ny  ,nab3)     ! save ru between calls
  real(r8), save ::  rv(Nx  ,Ny+1,nab3)     ! save rv between calls
  real(r8), save :: rzt(Nx  ,Ny  ,nab3)     ! save rzt between calls
  integer, save :: imod=-2                  ! save imod between calls

!.The following is required for AB3 steps
  select case (imod)
   case (-2)
    call swerhs(u,v,zt, ru(1,1,1),rv(1,1,1),rzt(1,1,1),h,dx,dy)
    call RK3Step(u,v,zt, dt,h,dx,dy)         ! start-up first RK3 step
    imod = -1                     ! indicator for 2nd start-up
   case (-1)
    call swerhs(u,v,zt, ru(1,1,2),rv(1,1,2),rzt(1,1,2),h,dx,dy)
    call RK3Step(u,v,zt, dt,h,dx,dy)         ! start-up 2nd RK3 step
    imod = 0                      ! indicator to switch to regular AB3
   case (0)
    call AB3Work(u,v,zt,ru,rv,rzt,dt,(/1,2,3/),h,dx,dy)
    imod = 1
   case (1)
    call AB3Work(u,v,zt,ru,rv,rzt,dt,(/2,3,1/),h,dx,dy)
    imod = 2
   case default
    call AB3Work(u,v,zt,ru,rv,rzt,dt,(/3,1,2/),h,dx,dy)
    imod = 0
  end select

  return
end subroutine AB3Step
!**********************************************************************
! AB3 steps to evaluate right hand side once per step.
!**********************************************************************
subroutine AB3Work(u,v,zt,ru,rv,rzt,dt,iflip,h,dx,dy)
  implicit none
  integer, intent(in)   :: iflip(nab3) ! order of AB3 right hand side
! integer, intent(in)   :: Nx,Ny
  real(r8), intent(in)  :: dt,dx,dy
  real(r8), intent(inout) ::  ru(Nx+1,Ny  ,nab3)
  real(r8), intent(inout) ::  rv(Nx  ,Ny+1,nab3)
  real(r8), intent(inout) :: rzt(Nx  ,Ny  ,nab3)
  real(r8), intent(inout):: u(Nx+1,Ny),v(Nx,Ny+1),zt(Nx,Ny),h(Nx,Ny)

!...3rd order Adams-Bashforth factors
  real(r8), parameter :: alpha(nab3)=(/  5.0_r8/12.0_r8, &
                                       -16.0_r8/12.0_r8, &
                                        23.0_r8/12.0_r8/)
  integer :: k,if
  real(r8):: dts

  k = iflip(nab3)            ! last index can be overwritten
  call swerhs(u,v,zt, ru(1,1,k),rv(1,1,k),rzt(1,1,k),h,dx,dy)
  do if = 1,nab3
   k = iflip(if)
   dts= dt*alpha(if)
   call InPlaceUpdate(u,v,zt, ru(1,1,k),rv(1,1,k),rzt(1,1,k), dts)
  enddo

  return
end subroutine AB3Work
!**********************************************************************
! Leap frog trapezoidal time-stepping
!**********************************************************************
subroutine LFTStep(u,v,zt,dt,h,dx,dy)
  real(r8), intent(in)    :: dt,dx,dy       ! time step size and grid sizes
  real(r8), intent(inout) ::  u(Nx+1,Ny  ) ! x-velocity
  real(r8), intent(inout) ::  v(Nx  ,Ny+1) ! y-velocity
  real(r8), intent(inout) :: zt(Nx  ,Ny  ) ! pressure
  real(r8), intent(inout) :: h(Nx  ,Ny  ) ! pressure

  real(r8), save ::  qu(Nx+1,Ny  )     ! save  u @ t_{n-1} between calls
  real(r8), save ::  qv(Nx  ,Ny+1)     ! save  v @ t_{n-1} between calls
  real(r8), save :: qzt(Nx  ,Ny  )     ! save zt @ t_{n-1} between calls
  integer, save :: imod=-1             ! save imod for initialization
  real(r8) ::  ru(Nx+1,Ny  ), rv(Nx  ,Ny+1), rzt(Nx  ,Ny  )
  real(r8) ::  su(Nx+1,Ny  ), sv(Nx  ,Ny+1), szt(Nx  ,Ny  )
  real(r8) :: dts

  select case (imod)
   case (-1)
    qu = u; qv=v; qzt=zt;             ! store t_n for next call
    call RK3Step(u,v,zt, dt, h,dx, dy)          ! start-up first RK3 step
    imod = 0                          ! switch to LeapFrog Trapezoidal
   case default
    call swerhs(u,v,zt, ru,rv,rzt,h,dx,dy)    ! slope at t_n
    dts = 2._r8*dt                    ! leap-frog step
    call InPlaceUpdate(qu,qv,qzt, ru,rv,rzt, dts) ! leap-frog estimate
!.............This version averages u,v,zt @ t_{n+1/2}
!   qu =  0.5_r8* ( u+ qu)            ! u  estimate at t_{n+1/2}
!   qv =  0.5_r8* ( v+ qv)            ! v  estimate at t_{n+1/2}
!   qzt = 0.5_r8* (zt+qzt)            ! zt estimate at t_{n+1/2}
!   call swerhs(qu,qv,qzt, ru,rv,rzt) ! slope at t_{n+1/2}
!   call InPlaceUpdate(u,v,zt, ru,rv,rzt, dt) ! final update for t_{n+1}
!.............This version averages time-tendencies of u,v,zt @ t_{n+1/2}
    call swerhs(qu,qv,qzt, su,sv,szt,h,dx,dy) ! slope at t_{n+1}
    ru  = 0.5_r8*( ru+ su)
    rv  = 0.5_r8*( rv+ sv)
    rzt = 0.5_r8*(rzt+szt)
!.............Final stages
    qu = u; qv=v; qzt=zt              ! store t_n for next call
    call InPlaceUpdate( u, v, zt, ru,rv,rzt, dt) ! final update for t_{n+1}
  end select

  return
end subroutine LFTStep
!**********************************************************************
! Update values
!**********************************************************************
subroutine Update(qu,qv,qzt, u,v,zt, ru,rv,rzt, dts)
  implicit none
  real(r8), intent(in) :: dts
  real(r8), intent(out):: qu(Nx+1,Ny),qv(Nx,Ny+1),qzt(Nx,Ny) ! updated
  real(r8), intent(in) ::  u(Nx+1,Ny), v(Nx,Ny+1), zt(Nx,Ny) ! oldvalue
  real(r8), intent(in) :: ru(Nx+1,Ny),rv(Nx,Ny+1),rzt(Nx,Ny) ! rhs

  qu  =  u + dts*ru          ! temporary value for u
  qv  =  v + dts*rv          ! temporary value for v
  qzt = zt + dts*rzt         ! temporary value for zt

  return
end subroutine Update
!**********************************************************************
! Update values in place
!**********************************************************************
subroutine InPlaceUpdate(u,v,zt, ru,rv,rzt, dts)
  implicit none
  real(r8), intent(in) :: dts
  real(r8), intent(in) :: ru(Nx+1,Ny),rv(Nx,Ny+1),rzt(Nx,Ny) ! rhs
  real(r8), intent(inout) ::  u(Nx+1,Ny), v(Nx,Ny+1), zt(Nx,Ny) ! oldvalue

   u =  u + dts*ru          ! new value for u
   v =  v + dts*rv          ! new value for v
  zt = zt + dts*rzt         ! new value for zt

  return
end subroutine InPlaceUpdate
!!**********************************************************************
!subroutine AB3Init(vd,Nx,Ny)
!  implicit none
!  integer, intent(in) :: Nx,Ny
!  type(AB3data), intent(inout) :: vd
!
!  vd%imod =-2
!
!  return
!end subroutine AB3Init
!---------------  Code that can't work with GFORTRAN ------------------
!!**********************************************************************
!! If momentum and pressure tendencies are allocatable arrays
!! we need to allocate the memory for them at run time and
!! initialize the AB3 counter
!!**********************************************************************
!subroutine AB3Init(vd,Nx,Ny)
!  implicit none
!  integer, intent(in) :: Nx,Ny
!  type(AB3data), intent(inout) :: vd
!  allocate( vd%ru( Nx+1,Ny  ,nab3) )
!  allocate( vd%rv( Nx  ,Ny+1,nab3) )
!  allocate( vd%rzt(Nx  ,Ny  ,nab3) )
!  vd%imod =-2
!  if (allocated(vd%ru)) then
!    print *,'allocated in AB3Init'
!  endif
!!---------------  Code that can't work with GFORTRAN ------------------
!  return
!end subroutine AB3Init

end module shtstep

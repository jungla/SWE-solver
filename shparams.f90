module shparams

  use kinds
  implicit none

!                   Global Variables
!The variables declared here are global variable except for spval
!which is simply a ridiculous value.
!All variables are initially set to spval to make sure they are initialized
!properly, otherwise trouble would be spotted easily. The number of cells are set
!as parameters to allow the compiler to see the size of the arrays so it
!perform its cost benefit analysis on the basis of the actual array sizes.
!

!....Rossby soliton problem
! integer, parameter  :: Nx= 90   ! nb of x-grid cells
! integer, parameter  :: Ny= 70   ! nb of y-grid cells
 integer, parameter  :: Nx= 90*2   ! nb of x-grid cells
 integer, parameter  :: Ny= 70*2   ! nb of y-grid cells
!  integer, parameter  :: Nx=384   ! nb of x-grid cells
!  integer, parameter  :: Ny=128   ! nb of y-grid cells
!  integer, parameter  :: Nx=384*2   ! nb of x-grid cells
!  integer, parameter  :: Ny=128*2   ! nb of y-grid cells

!....Double Gyre Problem
! integer, parameter  :: Nx=204   ! nb of x-grid cells
! integer, parameter  :: Ny=204   ! nb of y-grid cells

!....Stading wave problem
! integer, parameter  :: Nx=101   ! nb of x-grid cells
! integer, parameter  :: Ny=101   ! nb of y-grid cells

  integer, parameter :: Nx1=Nx+1  ! nb of x-cell edges
  integer, parameter :: Ny1=Ny+1  ! nb of y-cell edges

  real(r8), parameter, private :: spval=-999999999! ludicrous value
  real(r8) :: dx=spval                 ! x-grid spacing
  real(r8) :: dy=spval                 ! y-grid spacing
  real(r8) :: gravity=spval            ! gravity
  real(r8) :: depth=spval              ! depth
  real(r8) :: bdrag=spval              ! Bottom stress
  real(r8) :: viscosity=spval          ! viscosity

  contains
!*********************************************************************
! Read the problem parameters
!*********************************************************************
subroutine Read_Params(xmin,xmax,ymin,ymax,fc,beta,dt,ntimestep,isnap)
  implicit none
  real(r8), intent(out) :: xmin,xmax
  real(r8), intent(out) :: ymin,ymax
  real(r8), intent(out) :: dt          ! time-step size
  real(r8), intent(out) :: fc          ! Coriolis for beta-plane
  real(r8), intent(out) :: beta        ! beta factor for beta-plane
  integer, intent(out)  :: ntimestep   ! number of time steps
  integer, intent(out)  :: isnap       ! snapshot interval

  read(5,*) xmin
  read(5,*) xmax
  read(5,*) ymin
  read(5,*) ymax
  read(5,*) gravity
  read(5,*) depth
  read(5,*) bdrag
  read(5,*) viscosity
  read(5,*) fc
  read(5,*) beta
  read(5,*) dt
  read(5,*) ntimestep
  read(5,*) isnap

  return
end subroutine Read_Params

end module shparams

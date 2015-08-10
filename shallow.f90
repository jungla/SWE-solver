program shallow
  use kinds
  use shparams
  use grid   , only: DefineCellEdges
!  use shtstep, only: LFTStep                   ! AB3 time stepping
! use shtstep, only: AB3Step                   ! AB3 time stepping
! use shtstep, only: RK3Step                   ! RK3 time stepping
  use shtstep, only: RK4Step                   ! RK4 time stepping
  use shinit , only: uvh_initial
  use shrhs  , only : Set_Coriolis, Set_Wind, Diagnostics

  implicit none

  real(r8) :: u(Nx1,Ny),ue(Nx1,Ny)                     ! zonal velocity on zonal edges
  real(r8) :: v(Nx,Ny1),ve(Nx,Ny1)                     ! meridional velocity on meridional edges
  real(r8) :: zt(Nx,Ny),zte(Nx,Ny)                     ! zonal cell edges
  real(r8) :: h(Nx,Ny)                                 ! zonal cell edges

!      Diagnostics variables
  real(r8) :: voli,volume                   ! volume
  real(r8) :: TEi,TE                        ! total energy budgets
  real(r8) :: pvbudgeti ,pvbudget           ! potential vorticity budgets
  real(r8) :: pv2budgeti,pv2budget          ! potential enstrophy budgets
  real(r8) :: fc,beta                       ! stuff for Coriolis force

!      Time steps variables
  integer  :: itime, ntimestep, isnap
  real(r8) :: time, dt

!      Grid information dx and dy are in module shparams
  real(r8) :: xmin,xmax,xe(Nx1)             ! zonal domain edges, and cell edges
  real(r8) :: ymin,ymax,ye(Nx1)             ! meridional domain edges, and cell edges

!................Preprocessing section
  call Read_Params(xmin,xmax,ymin,ymax,fc,beta,dt,ntimestep,isnap) ! physical & numerical parameters
  call DefineCellEdges(xe,dx,xmin,xmax,Nx) ! cell edges in x
  call DefineCellEdges(ye,dy,ymin,ymax,Ny) ! cell edges in y
  call Set_Coriolis(ye,fc,beta)            ! Set the coriolis parameter
  call Set_Wind(xe,ye)                     ! Set the Wind

  h = 1000.0_r8

!.Open output files
  open(unit=11,file='fort.11',form='unformatted',&
       status='unknown', action='write') ! file for u
  open(unit=12,file='fort.12',form='unformatted',&
       status='unknown', action='write') ! file for v
  open(unit=13,file='fort.13',form='unformatted',&
       status='unknown', action='write') ! file for zt

  time = 0._r8
!.initialization of variables
  call uvh_initial(u,v,zt,xe,ye,Nx,Ny,time)

!  call History(u,v,zt)

  itime = 1
  call RK4Step(u,v,zt,dt,h,dx,dy)    ! RK4 time-step
  call Diagnostics(u,v,zt,h,voli,pvbudgeti,TEi,pv2budgeti)
  
  do itime = 2,ntimestep
    call RK4Step(u,v,zt,dt,h,dx,dy)    ! RK4 time-step
!   call RK3Step(u,v,zt,dt,h,dx,dy)    ! RK3 time-step
!   call AB3Step(u,v,zt,dt)    ! AB3 time-step
!   call LFTStep(u,v,zt,dt,h,dx,dy)    ! LF-Trap time-step
    if (mod(itime,isnap) == 0) then
      time = itime*dt
!.....Compare to an analytic (exact) solution if one is available
!.....Comment it out otherwise
      call History(u,v,zt)     ! Output the data
      call Diagnostics(u,v,zt,h,volume,pvbudget,TE,pv2budget)          ! diagnostics
      write(6,*) 'Volume Error= '   , volume-voli, volume
!      write(6,*) 'Vorticity Error= ',  pvbudget -  pvbudgeti, pvbudget
!      write(6,*) 'Enstrophy Error= ',  pv2budget - pv2budgeti, pv2budget
      write(6,*) 'TE -error= ',       TE-TEi, TE

    endif
  enddo

  close(11)
  close(12)
  close(13)

  contains
!***********************************************************************
!***********************************************************************
subroutine History(u,v,zt)
  implicit none
  real(r8), intent(in) ::  u(Nx1,Ny )
  real(r8), intent(in) ::  v(Nx ,Ny1)
  real(r8), intent(in) :: zt(Nx ,Ny )
  write(11) u
  write(12) v
  write(13) zt
!  write(*,*) 'u:', u
!  write(*,*) 'v:', v
!  write(*,*) 'z:', zt
  return
end subroutine History

end program shallow

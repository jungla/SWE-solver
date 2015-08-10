module cgridops

  use kinds
  implicit none

! Not yet implemented
! integer, parameter :: nbcperiod=0   ! periodic BC flag
! integer, parameter :: nbcdirich=1   ! Dirichlet BC flag
! integer, parameter :: nbcneuman=2   ! Neumann BC flag

  contains

!
! This module contains utilities designed to operate on variables
! layed out on a C-Grid configuration. The operations are primarily
! for computing averages and differences on arrays defined on u-v-p-z-points.
! These operations are needed to reconstruct functions,
! compute gradient, divergence, and curl. For these routines to work properly
! it is important for the user to conform to the convention used in defining
! the variables, particularly the geometrical patterns of the grid.
!
! The following grid layout is assumed: 
! There are M cells in the i direction and K=M+1 edges
! There are N cells in the j direction and L=N+1 edges
! The origin for the numbering system is the lower left corner and coincides
! with z11.
!
!     p0L  u1L  p1L  u2L  p2L  u3L  p3L  uiL  piL  |......  uML  pML  uKL  pKL
!           |         |         |         |        |         |         |
!     v0L--z1L--v1L--z2L--v2L--z3L--v3L--ziL--viL-----------zML--vML--zKL--vKL
!           |         |         |         |        |         |         |
!     p0N  u1N  p1N  u2N  p2N  u3N  p3N  uiN  piN  |......  uMN  pMN  uKN  pKN
!           |         |         |         |        |         |         |
!     v0N--z1N--v1N--z2N--v2N--z3N--v3N--ziN--viN--|--------zMN--vMN--zKN--vKN
!           |         |         |         |        |         |         |
!     p0j  u1j  p1j  u2j  p2j  u3j  p3j  uij  pij  |......  uMj  pMj  uKj  pKj
!           |         |         |         |        |         |         |
!  ^  v0j--z1j--v1j--z2j--v2j--z3j--v3j--zij--vij---------- zMj--vMj--zKj--vKj
! j|        |                                      |                   |   
!  |  p03  u13  p13  u23  p23  u33  p33  ui3  pi3  |......  uM3  pM3  uK3  pK3
!           |         |         |         |        |         |         |
!     v03--z13--v13--z23--v23--z33--v33--zi3--vi3---------- zM3--vM3--zK3--vK3
!           |         |         |         |        |         |         |
!     p02--u12  p12  u22  p22  u32  p32  ui2  pi2  |......  uM2  pM2  uK2  pK2
!           |         |         |         |        |         |         |
!     v02--z12--v12--z22  v22--z32--v32--zi2--vi2---------- zM2--vM2--zK2--vK2
!           |                                      |                   |   
!     p01--u11  p11  u21  p21  u31  p31  ui1  pi1  |......  uM1  pM1  uK1  pK1
!           |         |         |         |        |         |         |
!     v01--z11--v11--z21--v21--z31--v31--zi1--vi1---------- zM1--vM1--zK1--vK1
!           |         |         |         |        |         |         |
!     p00  u10  p10  u20  p20  u30  p30  ui0  pi0  |......  uM0  pM0  uK0  pK0
!                           ----> 
!                             i
! The variables on the left and the bottom with a 0-index are halo points,
! and so are the variables on the right and top with a K or L index.
! The halo points can be used to enforce Boundary Conditions on the primary variables
! The halo points can also be used for parallel communication
! when the code is run on distributed memory processors. These are currently not enabled
! in the present version of this software.
!
! The main subroutines in this module are xOP_nD, and there
! are corresponding versions for operations defined in the y-direction when the grid
! is two-dimensional: yOP_nD.
! xOP_nD presume that there are NO halo points, and so operations
! on the edges (u-points) are not performed for lack of information.
! It is the users responsibility to take care of these edge values in the output
! array.
!
! This module is large because it contains a lot of code for testing. The actual
! number of useful subroutines for numerical computation is just 3!
!

!*******************************************************************************
! Performs the following operation
! g(i-igo,j) = a0*g(i-igo,j) + a(1)*f(i,j) + a(2)*f(i-1,j) for i = 2,...,M+igo
!                                                              j = 1,...,N
! For an x-derivative set a1=1/dx and a2=-1/dx
! For an x-average    set a1=1/2  and a2= 1/2
! Set a0=0._r8 to assign result to array g 
! Set a0=1._r8 to add result to array g 
! To perform an operation from p to u-points set igo=0
!    In this case the end points of the array g are untouched and must be set by BC
! To perform an operation from u to p-points set igo=1
!*******************************************************************************
subroutine xOP_1D(g, f, a0, a, Mg,Mf, M,igo)
  implicit none
  integer,  intent(in)    :: Mg         ! leading dimension of array g
  integer,  intent(in)    :: Mf         ! leading dimension of array f
  integer,  intent(in)    :: M          ! number of x-cells
  integer,  intent(in)    :: igo        ! indexing offset for array g
  real(r8), intent(inout) :: g(Mg)      ! output array
  real(r8), intent(in)    :: f(Mf)      ! input array
  real(r8), intent(in)    :: a0, a(2)   ! coefficients
  integer :: i
  do i = 2,M+igo
    g(i-igo) = a0*g(i-igo) + a(1)*f(i)+a(2)*f(i-1);
  enddo
  return
end subroutine xOP_1D

!*******************************************************************************
! Performs the following operation
! g(i-igo,j) = a0*g(i-igo,j) + a1*f(i,j) + a2*f(i-1,j) for i = 2,...,M+igo
!                                                          j = 1,...,N
! For an x-derivative set a1=1/dx and a2=-1/dx
! For an x-average    set a1=1/2  and a2= 1/2
! Set a0=0._r8 to assign result to array g 
! Set a0=1._r8 to add result to array g 
! To perform an operation from p to u-points set igo=0
!    In this case the end points of the array g are untouched and must be set by BC
! To perform an operation from u to p-points set igo=1
!*******************************************************************************
subroutine xOP_2D(g, f, a0, a, Mg,Mf, M,igo)
  implicit none
  integer,  intent(in)    :: Mg         ! leading dimension of array g
  integer,  intent(in)    :: Mf         ! leading dimension of array f
  integer,  intent(in)    :: M(2)       ! number of cells in each direction
  integer,  intent(in)    :: igo        ! indexing offset for array g
  real(r8), intent(inout) :: g(Mg,M(2))
  real(r8), intent(in)    :: f(Mf,M(2))
  real(r8), intent(in)    :: a0, a(2)
  integer :: i,j
  do j = 1,M(2)
    do i = 2,M(1)+igo
      g(i-igo,j) = a0*g(i-igo,j) + a(1)*f(i,j)+a(2)*f(i-1,j);
    enddo
  enddo
  return
end subroutine xOP_2D


!************************!
! Double derivative in x !
!************************!
subroutine x2OP_2D(g, f, a0, a, Mg, Mf, M)
   implicit none
   integer,  intent(in)    :: Mg         ! leading dimension of array g
   integer,  intent(in)    :: Mf         ! leading dimension of array f
   integer,  intent(in)    :: M(2)       ! number of cells in each direction
   real(r8), intent(inout) :: g(Mg,M(2))
   real(r8), intent(in)    :: f(Mf,M(2))
   real(r8), intent(in)    :: a0, a(2)
   integer :: i,j
   do j = 1,M(2)
     do i = 2,M(1)-1
       g(i,j) = a0*g(i,j) + a(1)*f(i+1,j) - 2._r8*f(i,j) + a(2)*f(i-1,j);
     enddo
   enddo
   return
end subroutine x2OP_2D

!*******************************************************************************
! Performs the following operation
! g(i,j-jgo) = a0*g(i,j-jgo) + a1*f(i,j) + a2*f(i,j-1)
! For an y-derivative set a1=1/dy and a2=-1/dy
! For an y-average    set a1=1/2  and a2= 1/2
! Set a0=0.0_r8 to assign result to array g 
! Set a0=1.0_r8 to add result to array g 
! To perform an operation from p to v-points set jgo=0
!    In this case the end points of the array g are untouched and must be set by BC
! To perform an operation from v to p-points set jgo=1
!*******************************************************************************
subroutine yOP_2D(g, f, a0, a, Mg,Ng, Mf,Nf, M,jgo)
  implicit none
  integer,  intent(in)    :: Mg,Mf      ! leading dimension of array g,f
  integer,  intent(in)    :: Ng,Nf      ! column dimension of array  g,f
  integer,  intent(in)    :: M(2)       ! number of actual points in x-y
  integer,  intent(in)    :: jgo       ! offset of array for the storage part
  real(r8), intent(inout) :: g(Mg,Ng)
  real(r8), intent(in)    :: f(Mf,Nf)
  real(r8), intent(in)    :: a0, a(2)
  integer :: i,j
  do j = 2,M(2)+jgo
    do i = 1,M(1)
      g(i,j-jgo) = a0*g(i,j-jgo) + a(1)*f(i,j)+a(2)*f(i,j-1);
    enddo
  enddo
  return
end subroutine yOP_2D

subroutine y2OP_2D(g, f, a0, a, Mg,Ng, M)
  implicit none
  integer,  intent(in)    :: Mg,Ng      ! leading dimension of array g,f
  integer,  intent(in)    :: M(2)       ! number of actual points in x-y
  real(r8), intent(inout) :: g(Mg,Ng)
  real(r8), intent(in)    :: f(Mg,Ng)
  real(r8), intent(in)    :: a0, a(2)
  integer :: i,j
  do j = 2,M(2)-1
    do i = 1,M(1)
      g(i,j) = a0*g(i,j) + a(1)*f(i,j)-2._r8*f(i,j)+a(2)*f(i,j-1);
    enddo
  enddo
  return
end subroutine y2OP_2D

!*******************************************************************************
! The following is a  test case to illustrate how to use the above subroutines
! The input/output arrays are fx/gx where x is u,p,v,or z
! The gx arrays are initialized to a special number -9.99 to spot entries
! in the output arrays that have been skipped.
! The input arrays are linear in x and y to facilitate the validation
! The gradients in x should be a constant=2, and the one in y should also be constant.
! The averaging operations should return the exact array as for the input
!*******************************************************************************
subroutine Test2DCgridOps(fp,fu,fv,fz, &
                        gp,gu,gv,gz, M,N)
  implicit none
  integer, intent(in) :: M,N
  real(r8), intent(in) :: fp(M,N), fu(M+1,N), fv(M  ,N+1), fz(M+1,N+1)
  real(r8), intent(out):: gp(M,N), gu(M+1,N), gv(M  ,N+1), gz(M+1,N+1)

  real(r8), parameter :: am(2) =(/0.5_r8, 0.5_r8/)      ! averaging coefficients
  real(r8), parameter :: ad(2) =(/1.0_r8,-1.0_r8/)      ! gradient coefficients

  integer :: Mu(2),Mv(2),Mp(2)

  Mp = (/M  ,N  /) ! number of cells to process in p-point op
  Mu = (/M+1,N  /) ! number of cells to process in u-point op
  Mv = (/M  ,N+1/) ! number of cells to process in v-point op

  print *,'============x-averaging with xOP============'
  gu =-9.99_r8; gz =-9.99_r8; gv =-9.99_r8; gp =-9.99_r8
  call xOP_2D(gu, fp, 0._r8, am, M+1, M  , Mp, 0) ! p to u x-average
  call PrettyPrint4(fz,fv,gu,fp,M,M,N,'x-mean p to u-points')

  call xOP_2D(gp, fu, 0._r8, am, M  , M+1, Mp, 1) ! u to p x-average
  call PrettyPrint4(fz,fv,fu,gp,M,M,N,'x-mean u to p-points')

  call xOP_2D(gz, fv, 0._r8, am, M+1, M  , Mv, 0) ! v to z x-average
  call PrettyPrint4(gz,fv,fu,fp,M,M,N,'x-mean v to z-points')

  call xOP_2D(gv, fz, 0._r8, am, M  ,M+1 , Mv, 1) ! z to v x-average
  call PrettyPrint4(fz,gv,fu,fp,M,M,N,'x-mean z to v-points')


  print *,'============x-gradient with xOP============'
  gu =-9.99_r8; gz =-9.99_r8; gv =-9.99_r8; gp =-9.99_r8
  call xOP_2D(gu, fp, 0._r8, ad, M+1, M  , Mp, 0) ! p to u x-gradient
  call PrettyPrint4(fz,fv,gu,fp,M,M,N,'x-gradient p to u-points')

  call xOP_2D(gp, fu, 0._r8, ad, M  , M+1, Mp, 1) ! u to p x-gradient
  call PrettyPrint4(fz,fv,fu,gp,M,M,N,'x-gradient u to p-points')

  call xOP_2D(gz, fv, 0._r8, ad, M+1, M  , Mv, 0) ! v to z x-gradient
  call PrettyPrint4(gz,fv,fu,fp,M,M,N,'x-gradient v to z-points')

  call xOP_2D(gv, fz, 0._r8, ad, M  ,M+1 , Mv, 1) ! z to v x-gradient
  call PrettyPrint4(fz,gv,fu,fp,M,M,N,'x-gradient z to v-points')

  print *,'============y-averaging with yOP============'
  gu =-9.99_r8; gz =-9.99_r8; gv =-9.99_r8; gp =-9.99_r8
  call yOP_2D(gv, fp, 0._r8, am, M  , N+1, M  ,N  , Mp , 0) ! p to v y-average
  call PrettyPrint4(fz,gv,fu,fp,M,M,N,'y-mean from p to v-points')

  call yOP_2D(gp, fv, 0._r8, am, M  , N  , M  ,N+1, Mp , 1) ! v to p y-average
  call PrettyPrint4(fz,fv,fu,gp,M,M,N,'y-mean from v to p-points')

  call yOP_2D(gz, fu, 0._r8, am, M+1, N+1, M+1,N  , Mu , 0) ! u to z y-average
  call PrettyPrint4(gz,fv,fu,fp,M,M,N,'y-mean from u to z-points')

  call yOP_2D(gu, fz, 0._r8, am, M+1, N, M+1,N+1  , Mu , 1) ! z to u y-average
  call PrettyPrint4(fz,fv,gu,fp,M,M,N,'y-mean from z to u-points')

  print *,'============y-gradients with yOP============'
  gu =-9.99_r8; gz =-9.99_r8; gv =-9.99_r8; gp =-9.99_r8
  call yOP_2D(gv, fp, 0._r8, ad, M  , N+1, M  ,N  , Mp , 0) ! p to v y-gradient
  call PrettyPrint4(fz,gv,fu,fp,M,M,N,'y-gradient from p to v-points')

  call yOP_2D(gp, fv, 0._r8, ad, M  , N  , M  ,N+1, Mp , 1) ! v to p y-gradient
  call PrettyPrint4(fz,fv,fu,gp,M,M,N,'y-gradient from v to p-points')

  call yOP_2D(gz, fu, 0._r8, ad, M+1, N+1, M+1,N  , Mu , 0) ! u to z y-gradient
  call PrettyPrint4(gz,fv,fu,fp,M,M,N,'y-gradient from u to z-points')

  call yOP_2D(gu, fz, 0._r8, ad, M+1, N, M+1,N+1  , Mu , 1) ! z to u y-gradient
  call PrettyPrint4(fz,fv,gu,fp,M,M,N,'y-gradient from z to u-points')

  return
end subroutine Test2DCgridOps

!*******************************************************************************
!    Test the xOP routines, the output arrays are gu and gp
! Untouched edge values are flagged with a -9.99
!*******************************************************************************
subroutine Test1DCgridOps(fp,fu, gp,gu, M)
  implicit none
  integer, intent(in) :: M
  real(r8), intent(in) :: fp(M), fu(M+1)
  real(r8), intent(out):: gp(M), gu(M+1)

  real(r8), parameter :: am(2) =(/0.5_r8, 0.5_r8/)      ! averaging coefficients
  real(r8), parameter :: ad(2) =(/1.0_r8,-1.0_r8/)      ! gradient coefficients

  gu =-9.99_r8
  gp =-9.99_r8

  print *,'============x-averaging with xOP============'
  call xOP_1D(gu, fp, 0._r8, am, M+1, M  , M, 0) ! p to u x-average
  call PrettyPrint3(gu,fp,M,'x-mean p to u-points')

  call xOP_1D(gp, fu, 0._r8, am, M  , M+1, M, 1) ! u to p x-average
  call PrettyPrint3(fu,gp,M,'x-mean u to p-points')

  gu =-9.99_r8
  gp =-9.99_r8

  print *,'============x-gradient with xOP============'
  call xOP_1D(gu, fp, 0._r8, ad, M+1, M  , M, 0) ! p to u x-gradient
  call PrettyPrint3(gu,fp,M,'x-gradient p to u-points')

  call xOP_1D(gp, fu, 0._r8, ad, M  , M+1, M, 1) ! u to p x-gradient
  call PrettyPrint3(fu,gp,M,'x-gradient u to p-points')

  return
end subroutine Test1DCgridOps

subroutine Test1DCgridSetMatrices(fp,fu,M)
  implicit none
  integer, intent(in) :: M
  real(r8), intent(out) :: fp(M), fu(M+1)
  integer :: i

  do i = 1,M
   fp(i) = (2*i-1)
  enddo
  do i = 1,M+1
   fu(i) = 2*(i-1)
  enddo

  print *,'================Original Arrays---============'
  call PrettyPrint3(fu,fp,M,'Original Layout of u-p arrays')

  return
end subroutine Test1DCgridSetMatrices

subroutine Test2DCgridSetMatrices(fp,fu,fv,fz,M,N)
  implicit none
  integer, intent(in) :: M,N
  real(r8), intent(out) :: fp(M,N), fu(M+1,N), fv(M  ,N+1), fz(M+1,N+1)
  integer :: i,j

  do j = 1,N
   do i = 1,M
    fp(i,j) = (2*j-1)*(2*M+1) + (2*i-1)
   enddo
  enddo
  do j = 1,N+1
   do i = 1,M
    fv(i,j) = 2*(j-1)*(2*M+1) + (2*i-1)
   enddo
  enddo
  do j = 1,N
   do i = 1,M+1
    fu(i,j) = (2*j-1)*(2*M+1) + 2*(i-1)
   enddo
  enddo
  do j = 1,N+1
   do i = 1,M+1
    fz(i,j) = 2*(j-1)*(2*M+1) + 2*(i-1)
   enddo
  enddo

  print *,'================Original Arrays---============'
  call PrettyPrint4(fz,fv,fu,fp,M,M,N,'Original Layout')

  return
end subroutine Test2DCgridSetMatrices
!*****************************************************************
! The following is a collection of subroutines to print arrays
! on a C-grid.
!*****************************************************************
subroutine PrettyPrint(f,Mf,M,N,string)
  integer, intent(in) :: Mf
  integer, intent(in) :: M,N
  real(r8), intent(in):: f(Mf,N)
  character(len=*), intent(in):: string
  integer :: i,j
  write(6,*) string
  do j = N,1,-1
   do i = 1,M
    write(6,fmt='(f8.2,$)') f(i,j)
   enddo
   write(6,*)
 enddo
 return
end subroutine PrettyPrint
!*****************************************************************
! Print a single line alternating u and v-points
!*****************************************************************
subroutine PrettyPrint2(fu,fp,M,c)
  integer, intent(in) :: M
  real(r8), intent(in):: fu(M+1)
  real(r8), intent(in):: fp(M)
  character(len=1), optional :: c
  integer :: i
  if (present(c)) then
    do i = 1,M
      write(6,fmt='(1x,a1,f7.2,1x,a1,f7.2,$)') c,fu(i),c,fp(i)
    enddo
    write(6,fmt='(1x,a1,f7.2,1x,a1)') c, fu(M+1), c
  else
    do i = 1,M
     write(6,fmt='(2f8.2,$)') fu(i),fp(i)
    enddo
    write(6,fmt='(f8.2)') fu(M+1)
  endif
  return
end subroutine PrettyPrint2
!*****************************************************************
! Print bars on the u-points
!*****************************************************************
subroutine PrintBars(M)
  implicit none
  integer, intent(in) :: M
  integer :: i
  write(6,fmt='(6x,a1,$)') '|'
  do i = 2,M+1
   write(6,fmt='(17x,a1,$)') '|'
  enddo
  write(6,*)
end subroutine PrintBars
!*****************************************************************
! Print a set of u and p points indicating the u-points by bars
!*****************************************************************
subroutine PrettyPrint3(fu,fp,M,string)
  implicit none
  integer, intent(in) :: M
  real(r8), intent(in):: fu(M+1)
  real(r8), intent(in):: fp(M  )
  character(len=*), intent(in):: string
  write(6,*) string
  call PrintBars(M)
  call PrettyPrint2(fu(1),fp(1),M,' ')
  call PrintBars(M)
  return
end subroutine PrettyPrint3
!*****************************************************************
!*****************************************************************
subroutine PrettyPrint4(fz,fv,fu,fp,Mf,M,N,string)
  integer, intent(in) :: Mf
  integer, intent(in) :: M,N
  real(r8), intent(in):: fz(Mf+1,N+1)
  real(r8), intent(in):: fv(Mf  ,N+1)
  real(r8), intent(in):: fu(Mf+1,N)
  real(r8), intent(in):: fp(Mf  ,N)
  character(len=*), intent(in):: string
  integer :: j
  write(6,*) string
  call PrintBars(M)
  call PrettyPrint2(fz(1,N+1),fv(1,N+1),M,' ')
  do j = N,1,-1
    call PrettyPrint2(fu(1,j),fp(1,j),M,' ')
    call PrettyPrint2(fz(1,j),fv(1,j),M,' ')
  enddo
  call PrintBars(M)
  return
end subroutine PrettyPrint4

subroutine maintests
  implicit none
  integer, parameter :: M=4,N=5
  real(r8) :: fp(M,N), fu(M+1,N), fv(M  ,N+1), fz(M+1,N+1)
  real(r8) :: gp(M,N), gu(M+1,N), gv(M  ,N+1), gz(M+1,N+1)
  real(r8) :: fp1d(M), fu1d(M+1)
  real(r8) :: gp1d(M), gu1d(M+1)
  integer :: i

  print *,'===============TESTS of 1D ROUTINES================'
  call Test1DCgridSetMatrices(fp1d,fu1d,M)
  call Test1DCgridOps(fp1d,fu1d, gp1d,gu1d, M)
  print *,'===============END TESTS of 1D ROUTINES================'
  do i = 1,3
    print *,'======================================================='
  enddo

  print *,'===============TESTS of 2D ROUTINES================'
  call Test2DCgridSetMatrices(fp,fu,fv,fz,M,N)
  call Test2DCgridOps(fp,fu,fv,fz, gp,gu,gv,gz, M,N)
  print *,'===============END TESTS of 2D ROUTINES================'

  stop
end subroutine maintests

end module cgridops


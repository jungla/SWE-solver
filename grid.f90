module grid
 use kinds
 implicit none

contains

subroutine DefineCellEdges(x,dx,xmin,xmax,M)
  implicit none
  integer, intent(in) :: M             ! number of cells
  real(r8), intent(in) :: xmin, xmax   ! extrema of interval
  real(r8), intent(out):: x(M+1)       ! abcissa of edges
  real(r8), intent(out):: dx           ! grid size
  integer :: i
  dx = (xmax-xmin)/real(M,r8)
  x(1) = xmin
  do i = 2,M
   x(i) = (i-1)*dx + xmin
  enddo
  x(M+1) = xmax
  return
end subroutine DefineCellEdges

end module grid

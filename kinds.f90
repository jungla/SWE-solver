! kinds.f
! this module retrieves the kind of a double precision variable
module kinds
  implicit none
  integer, parameter :: r8 = selected_real_kind(15, 307)
end module kinds

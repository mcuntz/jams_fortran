module mo_mainvar
  !
  use mo_kind, only: i4, sp, dp
  !
  implicit none
  !
  real(dp), dimension(:),     allocatable, target :: t
  real(sp), dimension(:,:,:), allocatable, target :: data
  real(dp), dimension(:,:),   allocatable, target :: Lat
  real(dp), dimension(:,:),   allocatable, target :: Lon
  !
end module mo_mainvar

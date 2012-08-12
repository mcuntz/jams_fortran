! ------------------------------------------------------------------------------
!
! module mainvar
!
! author: Stephan Thober
!
! created: 10 Aug 2012
!
! ------------------------------------------------------------------------------
module mo_mainvar
  !
  use mo_kind, only: i4, sp, dp
  !
  implicit none
  !
  character(256) :: OutFile
  !
  integer(i4), dimension(:), allocatable         :: src_add
  integer(i4), dimension(:), allocatable         :: dst_add
  integer(i4), dimension(:), allocatable, target :: x
  integer(i4), dimension(:), allocatable, target :: y
  !
  real(sp)                                        :: nodata
  real(sp), dimension(:,:,:), allocatable         :: src_data
  real(sp), dimension(:,:,:), allocatable, target :: dst_data
  !
  real(dp), dimension(:,:),   allocatable         :: wts
  real(dp), dimension(:),     allocatable, target :: time
  real(dp), dimension(:,:),   allocatable, target :: lon
  real(dp), dimension(:,:),   allocatable, target :: lat
  !
  !
end module mo_mainvar

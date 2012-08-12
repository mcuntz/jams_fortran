module mo_readdata
  ! ------------------------------------------------------------------------------
  !
  ! This module contains all Variables / Routines that interact with
  ! reading data.
  !
  ! author: Stephan Thober
  !
  ! created: 10 Aug 2012
  !
  ! ------------------------------------------------------------------------------

  implicit none

  private

  public :: readdata

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !
  ! author: Stephan Thober
  !
  ! created: 10 Aug 2011
  !
  ! ------------------------------------------------------------------------------

  subroutine readdata

    use mo_kind,      only: i4, dp
    use mo_mainvar,   only: dst_add, src_add, wts, OutFile, nodata, src_data, dst_data, time, lon, lat
    use mo_ncread,    only: Get_NcDim, Get_NcVar
    use mo_constants, only: nin, PI_D

    implicit none

    character(256) :: InFile
    character(256) :: WtsFile

    integer(i4), dimension(5)              :: dl ! dimension length
    integer(i4), dimension(:), allocatable :: src_grid_dims
    integer(i4), dimension(:), allocatable :: dst_grid_dims

    real(dp), dimension(:), allocatable    :: tlon
    real(dp), dimension(:), allocatable    :: tlat

    namelist /RemapNamelist/ InFile, WtsFile, OutFile, nodata

    open(unit=nin, file="../FORTRAN_chs_lib/test/test_mo_remap/namelist", status="old", action="read")
    read(nin,RemapNamelist)
    close(nin)

    ! read srcdata
    dl = Get_NcDim(InFile, 'Z')
    allocate(src_data(dl(1), dl(2), dl(3)))
    call Get_NcVar(InFile, 'Z', src_data)

    dl = Get_NcDim(InFile, 'time')
    allocate(time(dl(1)))
    call Get_NcVar(InFile, 'time', time)

    ! read Wts File
    dl = Get_NcDim(WtsFile, 'src_grid_dims')
    allocate(src_grid_dims(dl(1)))
    call Get_NcVar(WtsFile, 'src_grid_dims', src_grid_dims)

    dl = Get_NcDim(WtsFile, 'dst_grid_dims')
    allocate(dst_grid_dims(dl(1)))
    call Get_NcVar(WtsFile, 'dst_grid_dims', dst_grid_dims)

    dl = Get_NcDim(WtsFile, 'dst_grid_center_lon')
    allocate(tlon(dl(1)))
    call Get_NcVar(WtsFile, 'dst_grid_center_lon', tlon)

    dl = Get_NcDim(WtsFile, 'dst_grid_center_lat')
    allocate(tlat(dl(1)))
    call Get_NcVar(WtsFile, 'dst_grid_center_lat', tlat)

    dl = Get_NcDim(WtsFile, 'src_address')
    allocate(src_add(dl(1)))
    call Get_NcVar(WtsFile, 'src_address', src_add)

    dl = Get_NcDim(WtsFile, 'dst_address')
    allocate(dst_add(dl(1)))
    call Get_NcVar(WtsFile, 'dst_address', dst_add)

    dl = Get_NcDim(WtsFile, 'remap_matrix')
    allocate(wts(dl(1), dl(2)))
    call Get_NcVar(WtsFile, 'remap_matrix', wts)

    ! perform consistency checks
    if ((size(src_data,1) /= src_grid_dims(1)) .or. &
        (size(src_data,2) /= src_grid_dims(2))) stop 'ERROR *** size mismatch in src_data. sub readdata'
   
    if ((size(src_add) /= size(wts,2)) .or. &
        (size(dst_add) /= size(wts,2))) stop 'ERROR *** size mismatch in wts. sub readdata'

    allocate(dst_data(dst_grid_dims(1), dst_grid_dims(2), size(src_data,3)))

    ! convert lat / lon from radians to degrees
    tlon = tlon * 180._dp / PI_D
    tlat = tlat * 180._dp / PI_D

    allocate(lon(dst_grid_dims(1), dst_grid_dims(2)))
    allocate(lat(dst_grid_dims(1), dst_grid_dims(2)))

    lon = reshape(tlon, (/ dst_grid_dims(1), dst_grid_dims(2) /))
    lat = reshape(tlat, (/ dst_grid_dims(1), dst_grid_dims(2) /))

  end subroutine readdata

end module mo_readdata

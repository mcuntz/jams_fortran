! ------------------------------------------------------------------------------
!
! Test Program for writing nc files using the netcdf4 library.
!
! author: Stephan Thober & Matthias Cuntz
!
! created: 04.11.2011
! last update: 25.11.2016
!
! ------------------------------------------------------------------------------
program main

  use mo_kind,    only: i4, sp
  use mo_NcRead,  only: Get_NcVar, Get_NcDim
  use mo_setnc,   only: setnc
  use mo_NcWrite, only: create_netcdf, close_netcdf, write_static_netcdf, write_dynamic_netcdf, V
  use mo_mainvar, only: lat, lon, data, t
  use mo_utils,   only: notequal

  implicit none

  integer(i4)                             :: ncid
  integer(i4), dimension(5)               :: dimlen
  character(256)                          :: Filename
  character(256)                          :: Varname
  real(sp), dimension(:,:,:), allocatable :: data1
  logical :: isgood

  isgood   = .true.
  Filename = 'pr_1961-2000.nc'

  ! --------------------------------------------------------------------
  ! read all variables

  Varname  = 'pr'
  dimlen = Get_NcDim(Filename,Varname)

  allocate(data(dimlen(1),dimlen(2),dimlen(3)))
  allocate(lat(dimlen(1),dimlen(2)))
  allocate(lon(dimlen(1),dimlen(2)))
  allocate(t(dimlen(3)))

  call Get_NcVar(Filename,Varname, data)

  ! read lat and lon
  Varname = 'lat'
  call Get_NcVar(Filename, Varname, lat)

  Varname = 'lon'
  call Get_NcVar(Filename, Varname, lon)

  Varname = 'time'
  call Get_NcVar(Filename, Varname, t)

  ! --------------------------------------------------------------------
  ! write nc file
  !

  Filename = 'ncwrite_make_check_test_file'

  ! 1st set netcdf structure V
  call setnc

  ! 2nd create actual netcdf file
  call create_netcdf(Filename, ncid)

  ! 3rd write static Variable
  call write_static_netcdf(ncid)

  ! 4th write dynamic Variables, 1st record
  V(3)%G0_d => t(1)
  V(4)%G2_f => data(:,:,1)
  call write_dynamic_netcdf(ncid,1)
  V(3)%G0_d => t(2)
  V(4)%G2_f => data(:,:,2)
  call write_dynamic_netcdf(ncid,2)

  ! last close netcdf files
  call close_netcdf(ncid)

  ! Check
  Varname  = 'pr'
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data1(size(data,1),size(data,2),size(data,3)))
  call Get_NcVar(Filename,Varname,data1)
  ! if (any(abs(data-data1) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data,data1))) isgood = .false.

  ! --------------------------------------------------------------------
  ! write nc4 file
  !

  Filename = 'ncwrite_make_check_test_file'
  call create_netcdf(Filename, ncid, netcdf4=.true.)
  call write_static_netcdf(ncid)
  V(3)%G0_d => t(1)
  V(4)%G2_f => data(:,:,1)
  call write_dynamic_netcdf(ncid,1)
  V(3)%G0_d => t(2)
  V(4)%G2_f => data(:,:,2)
  call write_dynamic_netcdf(ncid,2)
  call close_netcdf(ncid)
  ! Check
  Varname  = 'pr'
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data1)
  ! if (any(abs(data-data1) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data,data1))) isgood = .false.

  
  if (isgood) then
     write(*,*) 'mo_ncwrite o.k.'
  else
     write(*,*) 'mo_ncwrite failed!'
  endif
  
end program main

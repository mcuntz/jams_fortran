
! ------------------------------------------------------------------------------
!
! Test Program for writing nc files using the netcdf4 library.
!
! author: Stephan Thober
!
! created: 04.11.2011
! last update: 21.11.2011
!
! ------------------------------------------------------------------------------
program ReadNc
!
use mo_kind,    only: i4, sp
use mo_NcRead,  only: Get_NcVar, Get_NcDim
use mo_setnc,   only: setnc
use mo_NcWrite, only: create_netcdf, close_netcdf, write_static_netcdf, write_dynamic_netcdf, V
use mo_mainvar, only: lat, lon, data, t
!
integer(i4)                             :: ncid
integer(i4), dimension(5)               :: dimlen
character(256)                          :: Filename
character(256)                          :: Varname
LOGICAL :: isgood
real(sp), dimension(:,:,:), allocatable, target :: data1
!
Filename = '../FORTRAN_chs_lib/test/test_mo_ncwrite/pr_1961-2000.nc'
!
! read all variables -------------------------------------------------
!
Varname  = 'pr'
dimlen = Get_NcDim(Filename,Varname)
!
allocate(data(dimlen(1),dimlen(2),dimlen(3)))
allocate(lat(dimlen(1),dimlen(2)))
allocate(lon(dimlen(1),dimlen(2)))
allocate(t(dimlen(3)))
!
call Get_NcVar(Filename,Varname, data)
!
! read lat and lon
Varname = 'lat'
call Get_NcVar(Filename, Varname, lat)
!
Varname = 'lon'
call Get_NcVar(Filename, Varname, lon)
!
Varname = 'time'
call Get_NcVar(Filename, Varname, t)
!
! WRITE nc file -------------------------------------------------------
Filename = 'Test.nc'
!
! 1st set netcdf structure V
call setnc
!
! 2nd create actual netcdf file
call create_netcdf(Filename, ncid)
!
! 3rd write static Variable
call write_static_netcdf(ncid)
!
! 4th write dynamic Variables, 1st record
V(3)%G0_d => t(1)
V(4)%G2_f => data(:,:,1)
call write_dynamic_netcdf(ncid,1)
V(3)%G0_d => t(2)
V(4)%G2_f => data(:,:,2)
call write_dynamic_netcdf(ncid,2)
!
! last close netcdf files
call close_netcdf(ncid)
!
! Check
isgood = .true.
Filename = 'Test.nc'
Varname  = 'pr'
dimlen = Get_NcDim(Filename,Varname)
allocate(data1(size(data,1),size(data,2),size(data,3)))
call Get_NcVar(Filename,Varname,data1)
if (any(abs(data-data1) > epsilon(1.0_sp))) isgood = .false.
if (isgood) then
   write(*,*) 'mo_ncwrite o.k.'
else
   write(*,*) 'mo_ncwrite failed!'
endif
!
end program ReadNc

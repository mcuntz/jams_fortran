
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
use mo_kind,    only: i4, sp, dp
use mo_NcRead,  only: Get_NcVar
use mo_setnc,   only: setnc
use mo_NcWrite, only: create_netcdf, close_netcdf, write_static_netcdf, write_dynamic_netcdf
use mo_mainvar, only: lat, lon, data, t
!
integer(i4)                             :: ncid
character(256)                          :: Filename
character(256)                          :: Varname
LOGICAL :: isgood
real(sp), dimension(:,:,:), allocatable, target :: data1
!
Filename = '../FORTRAN_chs_lib/test/test_mo_ncwrite/pr_1961-2000.nc'
!
! read all variables -------------------------------------------------
Varname  = 'pr'
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
print*, lon
print*, 'setnc'
call setnc
print*, 'finished setnc'
!
! 2nd create actual netcdf file
call create_netcdf(Filename, ncid)
print*, ncid
!
! 3rd write static Variable
call write_static_netcdf(ncid)
!
! 4th write dynamic Variables, 1st record
call write_dynamic_netcdf(ncid,1)
!
! last close netcdf files
call close_netcdf(ncid)
!
! Check
isgood = .true.
Filename = 'Test.nc'
Varname  = 'pr'
allocate(data1(size(data,1),size(data,2),size(data,3)))
call Get_NcVar(Filename,Varname,data1)
if (any(abs(data-data1) > epsilon(1.0_sp))) isgood = .false.
print*, data(10,10,10), data1(10,10,10)
if (isgood) then
   write(*,*) 'mo_ncwrite o.k.'
else
   write(*,*) 'mo_ncwrite failed!'
endif
!
end program ReadNc

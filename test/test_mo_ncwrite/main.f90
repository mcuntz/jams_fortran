
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
!
Filename = 'pr_1961-2000.nc'
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
call setnc
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
end program ReadNc

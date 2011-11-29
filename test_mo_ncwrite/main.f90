! ------------------------------------------------------------------------------
!
! Test Program for reading nc files using the netcdf4 library.
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
use mo_NcInter, only: set_NcVar, create_netcdf, close_netcdf, write_static
!
real(sp), dimension(:,:,:), allocatable :: data
integer(i4)                               :: ncid
character(256)                            :: Filename
character(256), dimension(1)              :: Varname
!
Filename = 'pr_1961-2000.nc'
!
! Variable name can be retrieved by a "ncdump -h <filename>"
Varname  = 'pr'
!
call Get_NcVar(Filename,Varname(1), data)
!
write(*,*) 'sum of data: ', sum(data)
! The sum of the data should be 0.1174308 in single precision
!
Filename= 'Attributes/'
!
call set_NcVar(VarName, Filename)
Filename = 'Test.nc'
call create_netcdf(Filename, ncid)
print*, ncid
call write_static(ncid,4,data)
call close_netcdf(ncid)
!
end program ReadNc

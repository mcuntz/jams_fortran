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
use mo_kind,   only: sp, dp
use mo_NcRead, only: Get_NcVar
!
real(sp), dimension(:,:,:,:), allocatable :: data
character(256)                            :: Filename
character(256)                            :: Varname
!
Filename = 'test_mo_NcRead/pr_1961-2000.nc'
!
! Variable name can be retrieved by a "ncdump -h <filename>"
Varname  = 'pr'
!
call Get_NcVar(Filename,Varname, data)
!
write(*,*) 'sum of data: ', sum(data)
! The sum of the data should be 0.1174308 in single precision
!
end program ReadNc

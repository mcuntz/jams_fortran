! ------------------------------------------------------------------------------
!
! Test Program for reading nc files using the netcdf4 library.
!
! author: Stephan Thober
!
! created: 04.11.2011
! last update: 09.12.2011
!
! ------------------------------------------------------------------------------
program ReadNc
!
use mo_kind,   only: i4, sp
use mo_NcRead, only: Get_NcVar, get_ncdim
!
real(sp), dimension(:,:,:), allocatable :: data
character(256)                            :: Filename
character(256)                            :: Varname
integer(i4), dimension(5)                 :: dl
LOGICAL :: isgood
!
Filename = '../FORTRAN_chs_lib/test/test_mo_ncread/pr_1961-2000.nc'
!
! Variable name can be retrieved by a "ncdump -h <filename>"
Varname  = 'pr'
!
dl = get_ncdim(Filename, Varname)
!
allocate(data(dl(1),dl(2),dl(3)))
!
call Get_NcVar(Filename,Varname, data)
!
! The sum of the data should be 0.1174308 in single precision
!write(*,*) 'sum of data: ', sum(data)
isgood = .true.
isgood = isgood .and. (anint(1e7_sp*sum(data)) == 1174308._sp)
if (isgood) then
   write(*,*) 'mo_ncread o.k.'
else
   write(*,*) 'mo_ncread failed!'
endif
!
deallocate(data)
!
end program ReadNc

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
use mo_NcRead, only: Get_NcVar, get_ncdim, NcOpen, NcClose
!
real(sp), dimension(:,:,:), allocatable   :: data
character(256)                            :: Filename
character(256)                            :: Varname
integer(i4)                               :: id
integer(i4)                               :: i
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
data = -9999._sp
!
! check dynamic read
id = NcOpen(trim(Filename)) ! open file and get file handle
do i = 1, size(data,3)
   !
   call Get_NcVar(Filename, Varname, data(:,:,i), (/1,1,i/),(/dl(1),dl(2),1/), id)
   !
end do
call NcClose(id)            ! close file
isgood = isgood .and. (anint(1e7_sp*sum(data)) == 1174308._sp)
!
if (isgood) then
   write(*,*) 'mo_ncread o.k.'
else
   write(*,*) 'mo_ncread failed!'
endif
!
deallocate(data)
!
end program ReadNc

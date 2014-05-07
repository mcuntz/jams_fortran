! ------------------------------------------------------------------------------
!
! Test Program for reading nc files using the netcdf4 library.
!
! author: Stephan Thober
!
! created: 04.11.2011
! last update: 22.03.2014
!
! ------------------------------------------------------------------------------
program ReadNc
!
use mo_kind,   only: i4, sp, dp
use mo_NcRead, only: Get_NcVar, get_ncdim, NcOpen, NcClose
#ifndef ABSOFT
use mo_NcRead, only: Get_NcDimAtt, Get_NcVarAtt
#endif
!
real(sp)      , dimension(:,:,:), allocatable :: data
real(sp)      , dimension(:,:)  , allocatable :: tmp
character(256), dimension(:)    , allocatable :: DimNames
integer(i4)   , dimension(:)    , allocatable :: DimLen
real(dp)      , dimension(:)    , allocatable :: DimData
character(256)                                :: Filename
character(256)                                :: Varname
#ifndef ABSOFT
character(256)                                :: Attname
character(256)                                :: AttValues
#endif
integer(i4)                                   :: ncid
integer(i4)                                   :: NoDims
integer(i4)                                   :: i
integer(i4)   , dimension(5)                  :: dl
LOGICAL                                       :: isgood
!
Filename = '../FORTRAN_chs_lib/test/test_mo_ncread/pr_1961-2000.nc'
!
! Variable name can be retrieved by a "ncdump -h <filename>"
Varname  = 'pr'
!
isgood = .true.
!
dl = get_ncdim(Filename, Varname, ndims=NoDims)
!
allocate(data(dl(1),dl(2),dl(3)))

#ifndef ABSOFT
! get Dimesnion information - name & lenght (size)
call Get_NcDimAtt(Filename, Varname, DimNames, DimLen)
!
isgood = isgood .and. (DimNames(1) == 'x')
isgood = isgood .and. (DimNames(2) == 'y')
isgood = isgood .and. (DimNames(3) == 'time')
!
isgood = isgood .and. (DimLen(1) == 28)
isgood = isgood .and. (DimLen(2) == 36)
isgood = isgood .and. (DimLen(3) == 2)
#else
allocate(DimNames(3))
allocate(DimLen(3))
DimNames(1) = 'x'
DimNames(2) = 'y'
DimNames(3) = 'time'
DimLen(1) = 28
DimLen(2) = 36
DimLen(3) = 2
#endif
!
! read data corresponding to dimesnion 3 ('time')
allocate(DimData(DimLen(3)))
call Get_NcVar(Filename, DimNames(3) ,DimData)
isgood = isgood .and. (nint(sum(DimData)) == 8100_i4)
deallocate(DimData)
!
call Get_NcVar(Filename, Varname, data)
!
! The sum of the data should be 0.1174308 in single precision
! sun with aggresive optimisation calculates 0.117430955
isgood = isgood .and. (nint(1e6_sp*sum(data)) == 117431_i4)
data = -9999._sp
!
! check dynamic read
ncid = NcOpen(trim(Filename)) ! open file and get file handle
!
do i = 1, size(data,3)
   ! tmp is allocated within Get_NcVar
   call Get_NcVar(Filename, Varname, tmp, (/1,1,i/),(/dl(1),dl(2),1/), ncid)
   data(:,:,i) = tmp
   if ( allocated( tmp ) ) deallocate( tmp )
end do
!
call NcClose(ncid)            ! close file
!
isgood = isgood .and. (nint(1e6_sp*sum(data)) == 117431_i4)
!
#ifndef ABSOFT
! retrieving variables attributes
AttName='units'
call Get_NcVarAtt(FileName, trim(DimNames(3)), AttName, AttValues)
isgood = isgood .and. (AttValues == 'days since 1950-01-01 00:00:00')
!
call Get_NcVarAtt(FileName, 'pr', '_FillValue', AttValues)
isgood = isgood .and. (AttValues == '0.1000000E+31')
#endif
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

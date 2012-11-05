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
use mo_NcRead,  only: Get_NcVar, Get_NcDim
use mo_setnc,   only: setnc
use mo_NcWrite, only: create_netcdf, close_netcdf, write_static_netcdf, write_dynamic_netcdf, V, dump_netcdf
use mo_mainvar, only: lat, lon, data, t
!
implicit none
!
integer(i4)                             :: ncid, i, j
integer(i4), dimension(5)               :: dimlen
character(256)                          :: Filename
character(256)                          :: Varname
LOGICAL :: isgood
real(sp), dimension(:,:,:),     allocatable :: data1, data2
real(sp), dimension(:,:,:,:),   allocatable :: data3, data4
real(sp), dimension(:,:),       allocatable :: data5, data6
real(sp), dimension(:),         allocatable :: data7, data8
real(sp), dimension(:,:,:,:,:), allocatable :: data9, data10
real(dp), dimension(:,:,:),     allocatable :: ddata2
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
!
! Dump nc file -------------------------------------------------------
Filename = 'Test.nc'
Varname  = 'var'
! 1D
allocate(data7(size(data,1)))
data7(:) = data(:,1,1)
call dump_netcdf(Filename, data7)
dimlen = Get_NcDim(Filename,Varname)
allocate(data8(dimlen(1)))
call Get_NcVar(Filename,Varname,data8)
if (any(abs(data7-data8) > epsilon(1.0_sp))) isgood = .false.

! 2D
allocate(data5(size(data,1),size(data,2)))
data5(:,:) = data(:,:,1)
call dump_netcdf(Filename, data5)
dimlen = Get_NcDim(Filename,Varname)
allocate(data6(dimlen(1),dimlen(2)))
call Get_NcVar(Filename,Varname,data6)
if (any(abs(data5-data6) > epsilon(1.0_sp))) isgood = .false.

! 3D
call dump_netcdf(Filename, data)
dimlen = Get_NcDim(Filename,Varname)
allocate(data2(dimlen(1),dimlen(2),dimlen(3)))
call Get_NcVar(Filename,Varname,data2)
if (any(abs(data-data2) > epsilon(1.0_sp))) isgood = .false.

! 4D
allocate(data3(size(data,1),size(data,2),10,size(data,3)))
!data3 = spread(data,3,10) ! does not work with gfortran 4.6.1 on Mac OS X 10.7; works with NAG
do i=1, 10
   data3(:,:,i,:) = data
end do
call dump_netcdf(Filename, data3)
dimlen = Get_NcDim(Filename,Varname)
allocate(data4(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
call Get_NcVar(Filename,Varname,data4)
if (any(abs(data3-data4) > epsilon(1.0_sp))) isgood = .false.

! 5D
allocate(data9(size(data,1),size(data,2),10,8,size(data,3)))
do i=1, 10
   do j=1, 8
      data9(:,:,i,j,:) = data
   end do
end do
call dump_netcdf(Filename, data9)
dimlen = Get_NcDim(Filename,Varname)
allocate(data10(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
call Get_NcVar(Filename,Varname,data10)
if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.

! 3D - dp
call dump_netcdf(Filename, real(data,dp))
dimlen = Get_NcDim(Filename,Varname)
allocate(ddata2(dimlen(1),dimlen(2),dimlen(3)))
call Get_NcVar(Filename,Varname,ddata2)
if (any(abs(real(data,dp)-ddata2) > epsilon(1.0_dp))) isgood = .false.
!
if (isgood) then
   write(*,*) 'mo_ncwrite o.k.'
else
   write(*,*) 'mo_ncwrite failed!'
endif
!
end program ReadNc

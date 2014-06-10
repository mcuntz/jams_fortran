! ------------------------------------------------------------------------------
!
! Test Program for writing nc files using the netcdf4 library.
!
! author: Stephan Thober & Matthias Cuntz
!
! created: 04.11.2011
! last update: 30.11.2012
!
! ------------------------------------------------------------------------------
program write_netcdf

use mo_kind,    only: i4, sp, dp
use mo_NcRead,  only: Get_NcVar, Get_NcDim
use mo_setnc,   only: setnc
use mo_NcWrite, only: create_netcdf, close_netcdf, write_static_netcdf, write_dynamic_netcdf, V, dump_netcdf, var2nc
use mo_mainvar, only: lat, lon, data, t
use mo_utils,   only: notequal

implicit none

integer(i4)                             :: ncid, i, j
integer(i4), dimension(5)               :: dimlen
character(256)                          :: Filename
character(256)                          :: Varname
character(256), dimension(5)            :: dimname
character(256), dimension(1)            :: tname
LOGICAL :: isgood
real(sp), dimension(:,:,:),     allocatable :: data1, data2, data11, data12
real(sp), dimension(:,:,:,:),   allocatable :: data3, data4
real(sp), dimension(:,:),       allocatable :: data5, data6, data13, data14
real(sp), dimension(:),         allocatable :: data7, data8
real(sp), dimension(:,:,:,:,:), allocatable :: data9, data10
real(dp), dimension(:,:,:),     allocatable :: ddata2
real(dp), dimension(:,:),       allocatable :: ddata5
integer(i4), dimension(:,:,:),  allocatable :: idata2

isgood   = .true.
Filename = '../FORTRAN_chs_lib/test/test_mo_ncwrite/pr_1961-2000.nc'

! read all variables -------------------------------------------------

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

! WRITE nc file -------------------------------------------------------
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
if (any(abs(data-data1) > epsilon(1.0_sp))) isgood = .false.

! WRITE nc4 file -------------------------------------------------------
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
if (any(abs(data-data1) > epsilon(1.0_sp))) isgood = .false.

! var2nc file --------------------------------------------------------
Filename = 'ncwrite_make_check_test_file'
dimname(1) = 'lat'
dimname(2) = 'lon'
dimname(3) = 'time'
dimname(4) = 'xx'
dimname(5) = 'yy'
tname(1)   = 'time'
! write static data
call var2nc( Filename, data(:,:,1), dimname(1:2), 'pre_static', &
     longname = 'precipitation', units = '[mm/d]', fill_value = -9999. )
Varname = 'pre_static'
call Get_NcVar(Filename,Varname,data1)
if (any(abs(data-data1) > epsilon(1.0_sp))) isgood = .false.

! write time - 1d unlimit
call var2nc( Filename, int(t,i4), tname, 'time', dim_unlimit = 1_i4, &
     units = 'days since 1984-08-28', fill_value=-9999, is_dim = .True. )
! write variable
call var2nc( Filename, data(14,14,:), tname, 'pre_1d', dim_unlimit = 1_i4 , &
     longname = 'precipitation', units = '[mm/d]' )
! read again
Varname = 'pre_1d'
dimlen  = Get_NcDim(Filename,Varname)
allocate( data7( dimlen(1) ) )
call Get_NcVar(Filename,Varname,data7)
if (any(abs(data(14,14,:)-data7) > epsilon(1.0_sp))) isgood = .false.

! write 2d - dp
call var2nc( Filename, real(data(14,:,:),dp), dimname(2:3), 'pre_2d', dim_unlimit = 2,  &
     longname = 'precipitation', units = '[mm/d]' )
Varname = 'pre_2d'
allocate( ddata5( size( data, 2), size( data, 3) ) )
call Get_NcVar( Filename, Varname, ddata5 )
if (any(notequal(ddata5,real(data(14,:,:), dp)))) isgood = .false.
deallocate( ddata5 )

! write 3d - sp - specify append, if save variable should not be used
call var2nc( Filename, data, dimname(1:3), 'pre_3d', dim_unlimit = 3_i4 , &
     longname = 'precipitation', units = '[mm/d]', f_exists = .true. )
Varname = 'pre_3d'
call Get_NcVar(Filename,Varname,data1)
if (any(abs(data-data1) > epsilon(1.0_sp))) isgood = .false.

! write 4d - sp
allocate(data3(size(data,1),size(data,2),size(data,3),10))
do i=1, 10
   data3(:,:,:,i) = data
end do
call var2nc( Filename, data3, dimname(1:4), 'pre_4d', dim_unlimit = 3_i4 , &
     longname = 'precipitation', units = '[mm/d]' )
allocate( data4( size( data3,1), size(data3,2), size(data3,3), size(data3,4) ) )
call Get_NcVar(Filename,'pre_4d',data4)
if (any(notequal(data3,data4))) isgood = .false.

! write 5d - sp
allocate(data9(size(data,1),size(data,2),size(data,3),10,8))
do i=1, 10
   do j=1, 8
      data9(:,:,:,i,j) = data
   end do
end do
call var2nc( Filename, data9, dimname, 'pre_5d', dim_unlimit = 3_i4 , &
     longname = 'precipitation', units = '[mm/d]' )
Varname = 'pre_5d'
dimlen = Get_NcDim(Filename,Varname)
allocate(data10(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
call Get_NcVar(Filename,Varname,data10)
if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.

! clean up
deallocate( data3, data4, data7, data9, data10 )

! append time - 1d unlimit
call var2nc( Filename, int(t,i4)+2_i4, tname, 'time', dim_unlimit = 1_i4, &
     units = 'days since 1984-08-28', is_dim = .True. )

 ! append 1d - sp
call var2nc( Filename, data(14,14,:), tname, 'pre_1d', dim_unlimit = 1_i4 , &
     longname = 'precipitation', units = '[mm/d]', f_exists = .true. )
! check
Varname = 'pre_1d'
dimlen  = Get_NcDim(Filename,Varname)
allocate( data7( dimlen(1) ) )
call Get_NcVar(Filename,Varname,data7)
if (any(abs(data(14,14,:)-data7(3:)) > epsilon(1.0_sp))) isgood = .false.

! append 2d - dp
call var2nc( Filename, real(data(14,:,:),dp), dimname(2:3), 'pre_2d', dim_unlimit = 2, &
     longname = 'precipitation', units = '[mm/d]' )
Varname = 'pre_2d'
allocate( ddata5( size( data, 2), 2 * size( data, 3) ) )
call Get_NcVar( Filename, Varname, ddata5 )
if (any(notequal(ddata5(:,3:),real(data(14,:,:), dp)))) isgood = .false.

! append 3d
call var2nc( Filename, data, dimname(1:3), 'pre_3d', dim_unlimit = 3_i4, &
     longname = 'precipitation', units = '[mm/d]'  )
deallocate( data1 )
allocate( data1( size(data,1), size(data,2), 2*size(data,3) ) )
Varname = 'pre_3d'
call Get_NcVar(Filename,Varname,data1)
if (any(abs(data-data1(:,:,3:)) > epsilon(1.0_sp))) isgood = .false.
deallocate( data1 )
allocate( data1( size( data, 1), size( data, 2 ), size( data, 3 ) ) )

! append 4d
allocate(data3(size(data,1),size(data,2),size(data,3),10))
do i=1, 10
   data3(:,:,:,i) = data
end do
call var2nc( Filename, data3, dimname(1:4), 'pre_4d', dim_unlimit = 3_i4 , &
     longname = 'precipitation', units = '[mm/d]' )
allocate( data4( size( data3,1), size(data3,2), 2*size(data3,3), size(data3,4) ) )
call Get_NcVar(Filename,'pre_4d',data4)
if (any(notequal(data3,data4(:,:,3:,:)))) isgood = .false.

! append 5d
allocate(data9(size(data,1),size(data,2),size(data,3),10,8))
do i=1, 10
   do j=1, 8
      data9(:,:,:,i,j) = data
   end do
end do
call var2nc( Filename, data9, dimname, 'pre_5d', dim_unlimit = 3_i4 , &
     longname = 'precipitation', units = '[mm/d]', f_exists = .true. )
Varname = 'pre_5d'
dimlen = Get_NcDim(Filename,Varname)
allocate(data10(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
call Get_NcVar(Filename,Varname,data10)
if (any(abs(data9-data10(:,:,3:,:,:)) > epsilon(1.0_sp))) isgood = .false.

! cleanup
deallocate( data7, data9, data10, ddata5, data3, data4 )

! Dump nc file -------------------------------------------------------
Filename = 'ncwrite_make_check_test_file'
Varname  = 'var'
! 1D
allocate(data7(size(data,1)))
data7(:) = data(:,1,1)
call dump_netcdf(Filename, data7)
dimlen = Get_NcDim(Filename,Varname)
allocate(data8(dimlen(1)))
call Get_NcVar(Filename,Varname,data8)
if (any(abs(data7-data8) > epsilon(1.0_sp))) isgood = .false.

allocate(data14(dimlen(1),dimlen(2)))
call Get_NcVar(Filename,Varname,data14)
if (any(abs(data7-data14(:,1)) > epsilon(1.0_sp))) isgood = .false.

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

! 3D - i4
call dump_netcdf(Filename, int(data,i4))
dimlen = Get_NcDim(Filename,Varname)
allocate(idata2(dimlen(1),dimlen(2),dimlen(3)))
call Get_NcVar(Filename,Varname,idata2)
if (any(abs(int(data,i4)-idata2) > 0_i4)) isgood = .false.

! 1D - append
call dump_netcdf(Filename, data7)
call dump_netcdf(Filename, data7, append=.true.)
call dump_netcdf(Filename, data7, append=.true.)
dimlen = Get_NcDim(Filename,Varname)
allocate(data13(dimlen(1),dimlen(2)))
call Get_NcVar(Filename,Varname,data13)
if (any(abs(data7-data13(:,1)) > epsilon(1.0_sp))) isgood = .false.
if (any(abs(data7-data13(:,2)) > epsilon(1.0_sp))) isgood = .false.
if (any(abs(data7-data13(:,3)) > epsilon(1.0_sp))) isgood = .false.

! 2D - append
call dump_netcdf(Filename, data5)
call dump_netcdf(Filename, data5, append=.true.)
dimlen = Get_NcDim(Filename,Varname)
allocate(data12(dimlen(1),dimlen(2),dimlen(3)))
call Get_NcVar(Filename,Varname,data12)
if (any(abs(data5-data12(:,:,1)) > epsilon(1.0_sp))) isgood = .false.
if (any(abs(data5-data12(:,:,2)) > epsilon(1.0_sp))) isgood = .false.

! 3D - append
call dump_netcdf(Filename, data(:,:,1:size(data,3)/2))
call dump_netcdf(Filename, data(:,:,size(data,3)/2+1:), append=.true.)
dimlen = Get_NcDim(Filename,Varname)
allocate(data11(dimlen(1),dimlen(2),dimlen(3)))
call Get_NcVar(Filename,Varname,data11)
if (any(abs(data-data11) > epsilon(1.0_sp))) isgood = .false.

! 5D - large file support
call dump_netcdf(Filename, data9, lfs=.true.)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,data10)
if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.

! 5D - netcdf4, deflate=default
call dump_netcdf(Filename, data9, netcdf4=.true.)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,data10)
if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.

! 5D - netcdf4, deflate=0
call dump_netcdf(Filename, data9, netcdf4=.true., deflate_level=0)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,data10)
if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.

! 5D - netcdf4, deflate=9
call dump_netcdf(Filename, data9, netcdf4=.true., deflate_level=9)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,data10)
if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.


! 1D - netcdf4
call dump_netcdf(Filename, data7, netcdf4=.true.)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,data8)
if (any(abs(data7-data8) > epsilon(1.0_sp))) isgood = .false.
! 2D - netcdf4
call dump_netcdf(Filename, data5, netcdf4=.true.)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,data6)
if (any(abs(data5-data6) > epsilon(1.0_sp))) isgood = .false.
! 3D - netcdf4
call dump_netcdf(Filename, data, netcdf4=.true.)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,data2)
if (any(abs(data-data2) > epsilon(1.0_sp))) isgood = .false.
! 4D - netcdf4
call dump_netcdf(Filename, data3, netcdf4=.true.)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,data4)
if (any(abs(data3-data4) > epsilon(1.0_sp))) isgood = .false.
! 5D - netcdf4
call dump_netcdf(Filename, data9, netcdf4=.true.)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,data10)
if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.
! 3D - dp - netcdf4
call dump_netcdf(Filename, real(data,dp), netcdf4=.true.)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,ddata2)
if (any(abs(real(data,dp)-ddata2) > epsilon(1.0_dp))) isgood = .false.
! 3D - i4 - netcdf4
call dump_netcdf(Filename, int(data,i4), netcdf4=.true.)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,idata2)
if (any(abs(int(data,i4)-idata2) > 0_i4)) isgood = .false.
! 1D - append - netcdf4
call dump_netcdf(Filename, data7, netcdf4=.true.)
call dump_netcdf(Filename, data7, append=.true.)
call dump_netcdf(Filename, data7, append=.true.)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,data13)
if (any(abs(data7-data13(:,1)) > epsilon(1.0_sp))) isgood = .false.
if (any(abs(data7-data13(:,2)) > epsilon(1.0_sp))) isgood = .false.
if (any(abs(data7-data13(:,3)) > epsilon(1.0_sp))) isgood = .false.
! 2D - append - netcdf4
call dump_netcdf(Filename, data5, netcdf4=.true.)
call dump_netcdf(Filename, data5, append=.true.)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,data12)
if (any(abs(data5-data12(:,:,1)) > epsilon(1.0_sp))) isgood = .false.
if (any(abs(data5-data12(:,:,2)) > epsilon(1.0_sp))) isgood = .false.
! 3D - append - netcdf4
call dump_netcdf(Filename, data(:,:,1:size(data,3)/2), netcdf4=.true.)
call dump_netcdf(Filename, data(:,:,size(data,3)/2+1:), append=.true.)
dimlen = Get_NcDim(Filename,Varname)
call Get_NcVar(Filename,Varname,data11)
if (any(abs(data-data11) > epsilon(1.0_sp))) isgood = .false.

! ! 3D - netcdf4, deflate=0
! Filename = 'ncwrite_make_check_test_file'
! !deallocate(ddata2)
! allocate(ddata2(1000,1000,100))
! ddata2 = 1.0_dp
! call random_number(ddata2(:,:,1:50))
! call dump_netcdf(Filename, ddata2, netcdf4=.false., lfs=.true., deflate_level=0)


if (isgood) then
   write(*,*) 'mo_ncwrite o.k.'
else
   write(*,*) 'mo_ncwrite failed!'
endif

end program Write_netcdf

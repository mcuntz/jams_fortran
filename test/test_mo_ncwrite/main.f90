! ------------------------------------------------------------------------------
!
! Test Program for writing nc files using the netcdf4 library.
!
! author: Stephan Thober & Matthias Cuntz
!
! created: 04.11.2011
! last update: 29.05.2015
!
! ------------------------------------------------------------------------------
program write_netcdf

  use mo_kind,    only: i4, sp, dp
  use mo_NcRead,  only: Get_NcVar, Get_NcDim, Get_ncVarAtt
  use mo_setnc,   only: setnc
  use mo_NcWrite, only: create_netcdf, close_netcdf, write_static_netcdf, write_dynamic_netcdf, V, dump_netcdf, var2nc
  use mo_mainvar, only: lat, lon, data, t
  use mo_utils,   only: notequal
  use mo_timer,   only: max_timers, timers_init, timer_clear, timer_start, timer_stop, timer_get

  implicit none

  integer(i4)                             :: ncid, i, j, n
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
  real(dp), dimension(:,:,:,:,:), allocatable :: ddata9
  integer(i4), dimension(:,:,:),  allocatable :: idata2
  real(dp),    dimension(:,:),       allocatable :: lon1, lat1
  integer(i4), dimension(:),         allocatable :: x1, y1
  real(dp),    dimension(:),         allocatable :: time1
  character(256)                                 :: att1, att2, oriFilename
  character(256), dimension(:,:), allocatable :: attributes

  isgood   = .true.
  Filename = '../FORTRAN_chs_lib/test/test_mo_ncwrite/pr_1961-2000.nc'
  oriFilename = trim(Filename)


  ! --------------------------------------------------------------------
  ! read all variables

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



  ! --------------------------------------------------------------------
  ! write nc file
  !

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
  ! if (any(abs(data-data1) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data,data1))) isgood = .false.



  ! --------------------------------------------------------------------
  ! write nc4 file
  !

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
  ! if (any(abs(data-data1) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data,data1))) isgood = .false.



  ! --------------------------------------------------------------------
  ! var2nc file
  !

  Filename = 'ncwrite_make_check_test_file'
  dimname(1) = 'lat'
  dimname(2) = 'lon'
  dimname(3) = 'time'
  dimname(4) = 'tile'
  dimname(5) = 'depth'
  tname(1)   = 'time' ! tname must be array
  ! create attributes
  allocate( attributes(4,2) )
  attributes(1,1) = 'long_name'
  attributes(1,2) = 'precipitation'
  attributes(2,1) = 'units'
  attributes(2,2) = '[mm/d]'
  attributes(3,1) = 'missing_value'
  attributes(3,2) = '-9999.'
  attributes(4,1) = 'scale_factor'
  attributes(4,2) = '1.'

  ! write static data
  call var2nc( Filename, data(:,:,1), dimname(1:2), 'pre_static', &
       long_name = 'precipitation', units = '[mm/d]', missing_value = -9999., create=.true. )
  Varname = 'pre_static'
  call Get_NcVar(Filename,Varname,data1)
  if (any(notequal(data,data1))) isgood = .false.

  ! write time - 1d unlimit
  call var2nc( Filename, int(t,i4), tname, 'time', dim_unlimited = 1_i4, &
       units = 'days since 1984-08-28', missing_value=-9999 )
  ! write variable
  call var2nc( Filename, data(14,14,:), tname, 'pre_1d', dim_unlimited = 1_i4 , &
       attributes = attributes(:2,:) )

  ! read again
  Varname = 'pre_1d'
  dimlen  = Get_NcDim(Filename,Varname)
  allocate( data7( dimlen(1) ) )
  call Get_NcVar(Filename,Varname,data7)
  ! if (any(abs(data(14,14,:)-data7) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data(14,14,:),data7))) isgood = .false.

  ! write 2d - dp
  call var2nc( Filename, real(data(14,:,:),dp), dimname(2:3), 'pre_2d', dim_unlimited = 2,  &
       long_name = 'precipitation', units = '[mm/d]' )
  Varname = 'pre_2d'
  allocate( ddata5( size( data, 2), size( data, 3) ) )
  call Get_NcVar( Filename, Varname, ddata5 )
  if (any(notequal(ddata5,real(data(14,:,:), dp)))) isgood = .false.
  deallocate( ddata5 )

  ! write 3d - sp - specify append, if save variable should not be used
  call var2nc( Filename, data, dimname(1:3), 'pre_3d', dim_unlimited = 3_i4 , &
       long_name = 'precipitation', units = '[mm/d]' )
  Varname = 'pre_3d'
  call Get_NcVar(Filename,Varname,data1)
  ! if (any(abs(data-data1) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data,data1))) isgood = .false.

  ! write 4d - sp
  allocate(data3(size(data,1),size(data,2),size(data,3),10))
  do i=1, 10
     data3(:,:,:,i) = data
  end do
  call var2nc( Filename, data3, dimname(1:4), 'pre_4d', dim_unlimited = 3_i4 , &
       long_name = 'precipitation', units = '[mm/d]' )
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
  call var2nc( Filename, data9, dimname, 'pre_5d', dim_unlimited = 3_i4 , &
       long_name = 'precipitation', units = '[mm/d]' )
  Varname = 'pre_5d'
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data10(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
  call Get_NcVar(Filename,Varname,data10)
  ! if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data9,data10))) isgood = .false.

  ! clean up
  deallocate( data3, data4, data7, data9, data10 )

  ! append time - 1d unlimit
  call var2nc( Filename, int(t,i4)+2_i4, tname, 'time', dim_unlimited = 1_i4, &
       units = 'days since 1984-08-28' )

  ! append 1d - sp
  call var2nc( Filename, data(14,14,:), tname, 'pre_1d', dim_unlimited = 1_i4 , &
       long_name = 'precipitation', units = '[mm/d]' )
  ! check
  Varname = 'pre_1d'
  dimlen  = Get_NcDim(Filename,Varname)
  allocate( data7( dimlen(1) ) )
  call Get_NcVar(Filename,Varname,data7)
  ! if (any(abs(data(14,14,:)-data7(3:)) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data(14,14,:),data7(3:)))) isgood = .false.

  ! append 2d - dp
  call var2nc( Filename, real(data(14,:,:),dp), dimname(2:3), 'pre_2d', dim_unlimited = 2, &
       long_name = 'precipitation', units = '[mm/d]' )
  Varname = 'pre_2d'
  allocate( ddata5( size( data, 2), 2 * size( data, 3) ) )
  call Get_NcVar( Filename, Varname, ddata5 )
  if (any(notequal(ddata5(:,3:),real(data(14,:,:), dp)))) isgood = .false.

  ! append 3d
  call var2nc( Filename, data, dimname(1:3), 'pre_3d', dim_unlimited = 3_i4, &
       long_name = 'precipitation', units = '[mm/d]'  )
  deallocate( data1 )
  allocate( data1( size(data,1), size(data,2), 2*size(data,3) ) )
  Varname = 'pre_3d'
  call Get_NcVar(Filename,Varname,data1)
  ! if (any(abs(data-data1(:,:,3:)) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data,data1(:,:,3:)))) isgood = .false.
  deallocate( data1 )
  allocate( data1( size( data, 1), size( data, 2 ), size( data, 3 ) ) )

  ! append 4d
  allocate(data3(size(data,1),size(data,2),size(data,3),10))
  do i=1, 10
     data3(:,:,:,i) = data
  end do
  call var2nc( Filename, data3, dimname(1:4), 'pre_4d', dim_unlimited = 3_i4 , &
       long_name = 'precipitation', units = '[mm/d]' )
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
  call var2nc( Filename, data9, dimname, 'pre_5d', dim_unlimited = 3_i4 , &
       long_name = 'precipitation', units = '[mm/d]' )
  Varname = 'pre_5d'
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data10(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
  call Get_NcVar(Filename,Varname,data10)
  ! if (any(abs(data9-data10(:,:,3:,:,:)) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data9,data10(:,:,3:,:,:)))) isgood = .false.
  ! cleanup
  deallocate( data7, data9, data10, ddata5, data3, data4 )

  ! Try to reproduce exactly input file
  Filename = 'ncwrite_make_check_test_file'
  dimname(1) = 'x'
  dimname(2) = 'y'
  dimname(3) = 'time'
  allocate(x1(size(lon,1)), y1(size(lon,2)))
  forall(i=1:size(lon,1)) x1(i) = i
  forall(i=1:size(lon,2)) y1(i) = i
  ! write static dimensions
  call var2nc(Filename, x1, dimname(1:1), dimname(1), create = .true.)
  call var2nc(Filename, y1, dimname(2:2), dimname(2) )
  ! write static data
  call Get_NcVarAtt(oriFileName, 'lon', 'long_name', att1)
  call Get_NcVarAtt(oriFileName, 'lon', 'units', att2)
  call var2nc(Filename, lon, dimname(1:2), 'lon', long_name=trim(att1), units=trim(att2))
  call Get_NcVarAtt(oriFileName, 'lat', 'long_name', att1)
  call Get_NcVarAtt(oriFileName, 'lat', 'units', att2)
  call var2nc(Filename, lat, dimname(1:2), 'lat', long_name=trim(att1), units=trim(att2))
  ! write first time step - time
  call Get_NcVarAtt(oriFileName, 'time', 'units', att2)
  call var2nc(Filename, t(1:1), dimname(3:3), dimname(3), &
       dim_unlimited=1_i4, units=trim(att2) )
  ! write first time step - data
  call Get_NcVarAtt(oriFileName, 'pr', 'long_name', att1)
  call Get_NcVarAtt(oriFileName, 'pr', 'units', att2)
  call var2nc(Filename, data(:,:,1:1), dimname(1:3), 'pr', &
       dim_unlimited=3_i4, long_name=trim(att1), units=trim(att2))
  ! write second to last time step - time and data
  do i=2, size(t,1)
     call var2nc(Filename, t(i:i), dimname(3:3), dimname(3), dim_unlimited=1_i4)
     call var2nc(Filename, data(:,:,i:i), dimname(1:3), 'pr', dim_unlimited=3_i4)
  end do
  !
  ! Check variables and attributes
  dimlen = Get_NcDim(Filename,'lon')
  allocate(lon1(dimlen(1),dimlen(2)))
  call Get_NcVar(Filename,'lon',lon1)
  if (any(notequal(lon,lon1))) isgood = .false.
  call Get_NcVarAtt(oriFileName, 'lon', 'long_name', att1)
  call Get_NcVarAtt(FileName, 'lon', 'long_name', att2)
  if (trim(att1) /= trim(att2)) isgood = .false.
  call Get_NcVarAtt(oriFileName, 'lon', 'units', att1)
  call Get_NcVarAtt(FileName, 'lon', 'units', att2)
  if (trim(att1) /= trim(att2)) isgood = .false.
  !
  dimlen = Get_NcDim(Filename,'lat')
  allocate(lat1(dimlen(1),dimlen(2)))
  call Get_NcVar(Filename,'lat',lat1)
  if (any(notequal(lon,lon1))) isgood = .false.
  call Get_NcVarAtt(oriFileName, 'lat', 'long_name', att1)
  call Get_NcVarAtt(FileName, 'lat', 'long_name', att2)
  if (trim(att1) /= trim(att2)) isgood = .false.
  call Get_NcVarAtt(oriFileName, 'lat', 'units', att1)
  call Get_NcVarAtt(FileName, 'lat', 'units', att2)
  if (trim(att1) /= trim(att2)) isgood = .false.
  !
  dimlen = Get_NcDim(Filename,'time')
  allocate(time1(dimlen(1)))
  call Get_NcVar(Filename,'time',time1)
  if (any(notequal(t,time1))) isgood = .false.
  call Get_NcVarAtt(oriFileName, 'time', 'units', att1)
  call Get_NcVarAtt(FileName, 'time', 'units', att2)
  if (trim(att1) /= trim(att2)) isgood = .false.
  !
  dimlen = Get_NcDim(Filename,'pr')
  allocate(data2(dimlen(1),dimlen(2),dimlen(3)))
  call Get_NcVar(Filename,'pr',data2)
  if (any(notequal(data,data2))) isgood = .false.
  call Get_NcVarAtt(oriFileName, 'pr', 'long_name', att1)
  call Get_NcVarAtt(FileName, 'pr', 'long_name', att2)
  if (trim(att1) /= trim(att2)) isgood = .false.
  call Get_NcVarAtt(oriFileName, 'pr', 'units', att1)
  call Get_NcVarAtt(FileName, 'pr', 'units', att2)
  if (trim(att1) /= trim(att2)) isgood = .false.
  !
  ! cleanup
  deallocate(x1, y1, lon1, lat1, time1, data2)

  ! --------------------------------------------------------------------
  ! Dump nc file
  !

  Filename = 'ncwrite_make_check_test_file'
  Varname  = 'var'
  ! 1D
  allocate(data7(size(data,1)))
  data7(:) = data(:,1,1)
  call dump_netcdf(Filename, data7)
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data8(dimlen(1)))
  call Get_NcVar(Filename,Varname,data8)
  ! if (any(abs(data7-data8) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data7,data8))) isgood = .false.

  allocate(data14(dimlen(1),dimlen(2)))
  call Get_NcVar(Filename,Varname,data14)
  ! if (any(abs(data7-data14(:,1)) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data7,data14(:,1)))) isgood = .false.

  ! 2D
  allocate(data5(size(data,1),size(data,2)))
  data5(:,:) = data(:,:,1)
  call dump_netcdf(Filename, data5)
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data6(dimlen(1),dimlen(2)))
  call Get_NcVar(Filename,Varname,data6)
  ! if (any(abs(data5-data6) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data5,data6))) isgood = .false.

  ! 3D
  call dump_netcdf(Filename, data)
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data2(dimlen(1),dimlen(2),dimlen(3)))
  call Get_NcVar(Filename,Varname,data2)
  ! if (any(abs(data-data2) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data,data2))) isgood = .false.

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
  ! if (any(abs(data3-data4) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data3,data4))) isgood = .false.

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
  ! if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data9,data10))) isgood = .false.

  ! 3D - dp
  call dump_netcdf(Filename, real(data,dp))
  dimlen = Get_NcDim(Filename,Varname)
  allocate(ddata2(dimlen(1),dimlen(2),dimlen(3)))
  call Get_NcVar(Filename,Varname,ddata2)
  ! if (any(abs(real(data,dp)-ddata2) > epsilon(1.0_dp))) isgood = .false.
  if (any(notequal(real(data,dp),ddata2))) isgood = .false.

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
  ! if (any(abs(data7-data13(:,1)) > epsilon(1.0_sp))) isgood = .false.
  ! if (any(abs(data7-data13(:,2)) > epsilon(1.0_sp))) isgood = .false.
  ! if (any(abs(data7-data13(:,3)) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data7,data13(:,1)))) isgood = .false.
  if (any(notequal(data7,data13(:,2)))) isgood = .false.
  if (any(notequal(data7,data13(:,3)))) isgood = .false.

  ! 2D - append
  call dump_netcdf(Filename, data5)
  call dump_netcdf(Filename, data5, append=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data12(dimlen(1),dimlen(2),dimlen(3)))
  call Get_NcVar(Filename,Varname,data12)
  ! if (any(abs(data5-data12(:,:,1)) > epsilon(1.0_sp))) isgood = .false.
  ! if (any(abs(data5-data12(:,:,2)) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data5,data12(:,:,1)))) isgood = .false.
  if (any(notequal(data5,data12(:,:,2)))) isgood = .false.

  ! 3D - append
  call dump_netcdf(Filename, data(:,:,1:size(data,3)/2))
  call dump_netcdf(Filename, data(:,:,size(data,3)/2+1:), append=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data11(dimlen(1),dimlen(2),dimlen(3)))
  call Get_NcVar(Filename,Varname,data11)
  ! if (any(abs(data-data11) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data,data11))) isgood = .false.

  ! 5D - large file support
  call dump_netcdf(Filename, data9, lfs=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data10)
  ! if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data9,data10))) isgood = .false.

  ! 5D - netcdf4, deflate=default
  call dump_netcdf(Filename, data9, netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data10)
  ! if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data9,data10))) isgood = .false.

  ! 5D - netcdf4, deflate=0
  call dump_netcdf(Filename, data9, netcdf4=.true., deflate_level=0)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data10)
  ! if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data9,data10))) isgood = .false.

  ! 5D - netcdf4, deflate=9
  call dump_netcdf(Filename, data9, netcdf4=.true., deflate_level=9)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data10)
  ! if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data9,data10))) isgood = .false.


  ! 1D - netcdf4
  call dump_netcdf(Filename, data7, netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data8)
  ! if (any(abs(data7-data8) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data7,data8))) isgood = .false.
  ! 2D - netcdf4
  call dump_netcdf(Filename, data5, netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data6)
  ! if (any(abs(data5-data6) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data5,data6))) isgood = .false.
  ! 3D - netcdf4
  call dump_netcdf(Filename, data, netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data2)
  ! if (any(abs(data-data2) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data,data2))) isgood = .false.
  ! 4D - netcdf4
  call dump_netcdf(Filename, data3, netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data4)
  ! if (any(abs(data3-data4) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data3,data4))) isgood = .false.
  ! 5D - netcdf4
  call dump_netcdf(Filename, data9, netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data10)
  ! if (any(abs(data9-data10) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data9,data10))) isgood = .false.
  ! 3D - dp - netcdf4
  call dump_netcdf(Filename, real(data,dp), netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,ddata2)
  ! if (any(abs(real(data,dp)-ddata2) > epsilon(1.0_dp))) isgood = .false.
  if (any(notequal(real(data,dp),ddata2))) isgood = .false.
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
  ! if (any(abs(data7-data13(:,1)) > epsilon(1.0_sp))) isgood = .false.
  ! if (any(abs(data7-data13(:,2)) > epsilon(1.0_sp))) isgood = .false.
  ! if (any(abs(data7-data13(:,3)) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data7,data13(:,1)))) isgood = .false.
  if (any(notequal(data7,data13(:,2)))) isgood = .false.
  if (any(notequal(data7,data13(:,3)))) isgood = .false.
  ! 2D - append - netcdf4
  call dump_netcdf(Filename, data5, netcdf4=.true.)
  call dump_netcdf(Filename, data5, append=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data12)
  ! if (any(abs(data5-data12(:,:,1)) > epsilon(1.0_sp))) isgood = .false.
  ! if (any(abs(data5-data12(:,:,2)) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data5,data12(:,:,1)))) isgood = .false.
  if (any(notequal(data5,data12(:,:,2)))) isgood = .false.
  ! 3D - append - netcdf4
  call dump_netcdf(Filename, data(:,:,1:size(data,3)/2), netcdf4=.true.)
  call dump_netcdf(Filename, data(:,:,size(data,3)/2+1:), append=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data11)
  ! if (any(abs(data-data11) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data,data11))) isgood = .false.

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

  ! ! --------------------------------------------------------------------
  ! ! Time var2nc
  ! !

  ! Filename = 'ncwrite_make_check_test_file'
  ! dimname(1) = 'lat'
  ! dimname(2) = 'lon'
  ! dimname(3) = 'depth'
  ! dimname(4) = 'tile'
  ! dimname(5) = 'time'
  ! tname(1)   = 'time'  ! tname must be array
  ! if (allocated(attributes)) deallocate(attributes)
  ! allocate(attributes(3,2))
  ! attributes(1,1) = 'long_name'
  ! attributes(1,2) = 'precipitation'
  ! attributes(2,1) = 'units'
  ! attributes(2,2) = '[mm/d]'
  ! attributes(3,1) = 'missing_value'
  ! attributes(3,2) = '-9999.'

  ! n = 42
  ! if (allocated(time1)) deallocate(time1)
  ! allocate(time1(n))
  ! if (allocated(ddata9)) deallocate(ddata9)
  ! allocate(ddata9(90,180,3,5,n))

  ! call timers_init
  ! do i=1, max_timers
  !    call timer_clear(i)
  ! end do

  ! forall(i=1:n) time1(i) = i
  ! ! normal var2nc
  ! call timer_start(1)
  ! call var2nc(Filename, time1(1:1), dimname(5:5), 'time', dim_unlimited=1_i4, create=.true.)
  ! do i=1, n
  !    ddata9 = real(i, dp)
  !    if (i /= 1) call var2nc(Filename, time1(i:i), dimname(5:5), 'time', dim_unlimited=1_i4)
  !    call var2nc(Filename, ddata9(:,:,:,:,i), dimname(1:5), 'prec5', dim_unlimited=5_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp)
  !    call var2nc(Filename, ddata9(1,:,:,:,i), dimname(2:5), 'prec4', dim_unlimited=4_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp)
  !    call var2nc(Filename, ddata9(1,1,:,:,i), dimname(3:5), 'prec3', dim_unlimited=3_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp)
  !    call var2nc(Filename, ddata9(1,1,1,:,i), dimname(4:5), 'prec2', dim_unlimited=2_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp)
  !    call var2nc(Filename, ddata9(1,1,1,1,i:i), dimname(5:5), 'prec1', dim_unlimited=1_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp)
  ! end do
  ! call timer_stop(1)
  ! write(*,*) 'var2nc without ncid [s]: ', timer_get(1)

  ! ncid = -1_i4
  ! ! var2nc leaving file open
  ! call timer_start(2)
  ! call var2nc(Filename, time1(1:1), dimname(5:5), 'time', dim_unlimited=1_i4, ncid=ncid, create=.true.)
  ! do i=1, n
  !    ddata9 = real(i, dp)
  !    if (i /= 1) call var2nc(Filename, time1(i:i), dimname(5:5), 'time', dim_unlimited=1_i4, ncid=ncid)
  !    call var2nc(Filename, ddata9(:,:,:,:,i), dimname(1:5), 'prec5', dim_unlimited=5_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp, ncid=ncid)
  !    call var2nc(Filename, ddata9(1,:,:,:,i), dimname(2:5), 'prec4', dim_unlimited=4_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp, ncid=ncid)
  !    call var2nc(Filename, ddata9(1,1,:,:,i), dimname(3:5), 'prec3', dim_unlimited=3_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp, ncid=ncid)
  !    call var2nc(Filename, ddata9(1,1,1,:,i), dimname(4:5), 'prec2', dim_unlimited=2_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp, ncid=ncid)
  !    call var2nc(Filename, ddata9(1,1,1,1,i:i), dimname(5:5), 'prec1', dim_unlimited=1_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp, ncid=ncid)
  ! end do
  ! call close_netcdf(ncid)
  ! call timer_stop(2)
  ! write(*,*) 'var2nc with ncid [s]: ', timer_get(2)

  ! ncid = -1_i4
  ! ! var2nc leaving file open and setting record directly
  ! call timer_start(3)
  ! call var2nc(Filename, time1(1:1), dimname(5:5), 'time', dim_unlimited=1_i4, ncid=ncid, create=.true.)
  ! do i=1, n
  !    ddata9 = real(i, dp)
  !    if (i /= 1) call var2nc(Filename, time1(i:i), dimname(5:5), 'time', dim_unlimited=1_i4, ncid=ncid, nrec=i)
  !    call var2nc(Filename, ddata9(:,:,:,:,i), dimname(1:5), 'prec5', dim_unlimited=5_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp, ncid=ncid, nrec=i)
  !    call var2nc(Filename, ddata9(1,:,:,:,i), dimname(2:5), 'prec4', dim_unlimited=4_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp, ncid=ncid, nrec=i)
  !    call var2nc(Filename, ddata9(1,1,:,:,i), dimname(3:5), 'prec3', dim_unlimited=3_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp, ncid=ncid, nrec=i)
  !    call var2nc(Filename, ddata9(1,1,1,:,i), dimname(4:5), 'prec2', dim_unlimited=2_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp, ncid=ncid, nrec=i)
  !    call var2nc(Filename, ddata9(1,1,1,1,i:i), dimname(5:5), 'prec1', dim_unlimited=1_i4, &
  !         attributes=attributes(:3,:), missing_value=-9999._dp, ncid=ncid, nrec=i)
  ! end do
  ! call close_netcdf(ncid)
  ! call timer_stop(3)
  ! write(*,*) 'var2nc with ncid and nrec [s]: ', timer_get(3)

  ! deallocate(attributes)
  ! deallocate(ddata9)
  
end program write_netcdf

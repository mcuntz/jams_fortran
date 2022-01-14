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
program main

  use mo_kind,   only: i4, sp, dp
  use mo_ansi_colors, only: color, c_red, c_green
  use mo_ncread, only: Get_NcVar, Get_NcDim, Get_ncVarAtt
  use mo_var2nc, only: var2nc
  use mo_utils,  only: notequal

  implicit none

  integer(i4)                             :: i, j
  ! integer(i4)                             :: n
  integer(i4), dimension(5)               :: dimlen
  character(256)                          :: Filename
  character(256)                          :: Varname
  character(256), dimension(5)            :: dimname
  character(256), dimension(1)            :: tname
  LOGICAL :: isgood
  real(sp), dimension(:,:,:),     allocatable :: data1, data2
  real(sp), dimension(:,:,:,:),   allocatable :: data3, data4
  real(sp), dimension(:,:),       allocatable :: data5
  real(sp), dimension(:),         allocatable :: data7
  real(sp), dimension(:,:,:,:,:), allocatable :: data9, data10
  real(dp), dimension(:,:),       allocatable :: ddata5
  ! real(dp), dimension(:,:,:,:,:), allocatable :: ddata9
  real(dp),    dimension(:,:),       allocatable :: lon1, lat1
  integer(i4), dimension(:),         allocatable :: x1, y1
  real(dp),    dimension(:),         allocatable :: time1
  character(256)                                 :: att1, att2, oriFilename
  character(256), dimension(:,:), allocatable :: attributes
  !
  real(dp), dimension(:),     allocatable, target :: t
  real(sp), dimension(:,:,:), allocatable, target :: data
  real(dp), dimension(:,:),   allocatable, target :: Lat
  real(dp), dimension(:,:),   allocatable, target :: Lon

  isgood   = .true.
  Filename = 'pr_1961-2000.nc'
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
  ! var2nc file
  !

  Filename = 'var2nc_make_check_test_file'
  dimname(1) = 'lat'
  dimname(2) = 'lon'
  dimname(3) = 'tile'
  dimname(4) = 'depth'
  dimname(5) = 'time'
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
  allocate(data5(size(data,1),size(data,2)))
  call Get_NcVar(Filename,Varname,data5)
  if (any(notequal(data(:,:,1),data5))) isgood = .false.

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
  call var2nc( Filename, real(data(14,:,:),dp), (/dimname(2), dimname(5)/), 'pre_2d', dim_unlimited = 2,  &
       long_name = 'precipitation', units = '[mm/d]' )
  Varname = 'pre_2d'
  allocate( ddata5( size( data, 2), size( data, 3) ) )
  call Get_NcVar( Filename, Varname, ddata5 )
  if (any(notequal(ddata5,real(data(14,:,:), dp)))) isgood = .false.
  deallocate( ddata5 )

  ! write 3d - sp - specify append, if save variable should not be used
  call var2nc( Filename, data, (/dimname(1), dimname(2), dimname(5)/), 'pre_3d', dim_unlimited = 3_i4 , &
       long_name = 'precipitation', units = '[mm/d]' )
  Varname = 'pre_3d'
  call Get_NcVar(Filename,Varname,data1)
  ! if (any(abs(data-data1) > epsilon(1.0_sp))) isgood = .false.
  if (any(notequal(data,data1))) isgood = .false.

  ! write 4d - sp
  allocate(data3(size(data,1),size(data,2),10,size(data,3)))
  do i=1, 10
     data3(:,:,i,:) = data
  end do
  call var2nc( Filename, data3, (/dimname(1), dimname(2), dimname(3), dimname(5)/), 'pre_4d', dim_unlimited = 4_i4 , &
       long_name = 'precipitation', units = '[mm/d]' )
  allocate( data4( size( data3,1), size(data3,2), size(data3,3), size(data3,4) ) )
  call Get_NcVar(Filename,'pre_4d',data4)
  if (any(notequal(data3,data4))) isgood = .false.

  ! write 5d - sp
  allocate(data9(size(data,1),size(data,2),10,8,size(data,3)))
  do i=1, 10
     do j=1, 8
        data9(:,:,i,j,:) = data
     end do
  end do
  call var2nc( Filename, data9, dimname, 'pre_5d', dim_unlimited = 5_i4 , &
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
  call var2nc( Filename, real(data(14,:,:),dp), (/dimname(2), dimname(5)/), 'pre_2d', dim_unlimited = 2, &
       long_name = 'precipitation', units = '[mm/d]' )
  Varname = 'pre_2d'
  allocate( ddata5( size( data, 2), 2 * size( data, 3) ) )
  call Get_NcVar( Filename, Varname, ddata5 )
  if (any(notequal(ddata5(:,3:),real(data(14,:,:), dp)))) isgood = .false.

  ! append 3d
  call var2nc( Filename, data, (/dimname(1), dimname(2), dimname(5)/), 'pre_3d', dim_unlimited = 3_i4, &
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
  allocate(data3(size(data,1),size(data,2),10,size(data,3)))
  do i=1, 10
     data3(:,:,i,:) = data
  end do
  call var2nc( Filename, data3, (/dimname(1), dimname(2), dimname(3), dimname(5)/), 'pre_4d', dim_unlimited = 4_i4 , &
       long_name = 'precipitation', units = '[mm/d]' )
  allocate( data4( size( data3,1), size(data3,2), size(data3,3), 2*size(data3,4) ) )
  call Get_NcVar(Filename,'pre_4d',data4)
  if (any(notequal(data3,data4(:,:,:,3:)))) isgood = .false.

  ! append 5d
  allocate(data9(size(data,1),size(data,2),10,8,size(data,3)))
  do i=1, 10
     do j=1, 8
        data9(:,:,i,j,:) = data
     end do
  end do
  call var2nc( Filename, data9, dimname, 'pre_5d', dim_unlimited = 5_i4 , &
       long_name = 'precipitation', units = '[mm/d]' )
  Varname = 'pre_5d'
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data10(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
  call Get_NcVar(Filename,Varname,data10)
  if (any(notequal(data9,data10(:,:,:,:,3:)))) isgood = .false.
  ! cleanup
  deallocate( data7, data9, data10, ddata5, data3, data4 )

  ! Try to reproduce exactly input file
  Filename = 'var2nc_make_check_test_file'
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

  if (isgood) then
     write(*,*) 'mo_var2nc ', color('o.k.', c_green)
  else
     write(*,*) 'mo_var2nc ', color('failed!', c_red)
  endif

end program main

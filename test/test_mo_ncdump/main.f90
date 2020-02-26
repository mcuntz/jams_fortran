! ------------------------------------------------------------------------------
!
! Test Program for dumping nc files using the netcdf3/4 library.
!
! author: Matthias Cuntz & Stephan Thober
!
! created: 04.11.2011
! last update: 25.11.2016
!
! ------------------------------------------------------------------------------
program main

  use mo_kind,    only: i4, sp, dp
  use mo_ansi_colors, only: color, c_red, c_green
  use mo_NcRead,  only: Get_NcVar, Get_NcDim
  use mo_ncdump,  only: dump_netcdf
  use mo_utils,   only: notequal

  implicit none

  integer(i4)                             :: i, j
  integer(i4), dimension(5)               :: dimlen
  character(256)                          :: Filename
  character(256)                          :: Varname
  LOGICAL :: isgood
  real(sp), dimension(:,:,:),     allocatable :: data2, data11, data12
  real(sp), dimension(:,:,:,:),   allocatable :: data3, data4
  real(sp), dimension(:,:),       allocatable :: data5, data6, data13, data14
  real(sp), dimension(:),         allocatable :: data7, data8
  real(sp), dimension(:,:,:,:,:), allocatable :: data9, data10
  real(dp), dimension(:,:,:),     allocatable :: ddata2
  integer(i4), dimension(:,:,:),  allocatable :: idata2

  real(dp), dimension(:),     allocatable, target :: t
  real(sp), dimension(:,:,:), allocatable, target :: data
  real(dp), dimension(:,:),   allocatable, target :: Lat
  real(dp), dimension(:,:),   allocatable, target :: Lon

  isgood   = .true.
  Filename = 'pr_1961-2000.nc'

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
  ! Dump nc file
  !

  Filename = 'ncdump_make_check_test_file'
  Varname  = 'var'
  ! 1D
  allocate(data7(size(data,1)))
  data7(:) = data(:,1,1)
  call dump_netcdf(Filename, data7)
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data8(dimlen(1)))
  call Get_NcVar(Filename,Varname,data8)
  if (any(notequal(data7,data8))) isgood = .false.

  allocate(data14(dimlen(1),dimlen(2)))
  call Get_NcVar(Filename,Varname,data14)
  if (any(notequal(data7,data14(:,1)))) isgood = .false.

  ! 2D
  allocate(data5(size(data,1),size(data,2)))
  data5(:,:) = data(:,:,1)
  call dump_netcdf(Filename, data5)
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data6(dimlen(1),dimlen(2)))
  call Get_NcVar(Filename,Varname,data6)
  if (any(notequal(data5,data6))) isgood = .false.

  ! 3D
  call dump_netcdf(Filename, data)
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data2(dimlen(1),dimlen(2),dimlen(3)))
  call Get_NcVar(Filename,Varname,data2)
  if (any(notequal(data,data2))) isgood = .false.

  ! 4D
  allocate(data3(size(data,1),size(data,2),10,size(data,3)))
  do i=1, 10
     data3(:,:,i,:) = data
  end do
  call dump_netcdf(Filename, data3)
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data4(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
  call Get_NcVar(Filename,Varname,data4)
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
  if (any(notequal(data9,data10))) isgood = .false.

  ! 3D - dp
  call dump_netcdf(Filename, real(data,dp))
  dimlen = Get_NcDim(Filename,Varname)
  allocate(ddata2(dimlen(1),dimlen(2),dimlen(3)))
  call Get_NcVar(Filename,Varname,ddata2)
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
  if (any(notequal(data7,data13(:,1)))) isgood = .false.
  if (any(notequal(data7,data13(:,2)))) isgood = .false.
  if (any(notequal(data7,data13(:,3)))) isgood = .false.

  ! 2D - append
  call dump_netcdf(Filename, data5)
  call dump_netcdf(Filename, data5, append=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data12(dimlen(1),dimlen(2),dimlen(3)))
  call Get_NcVar(Filename,Varname,data12)
  if (any(notequal(data5,data12(:,:,1)))) isgood = .false.
  if (any(notequal(data5,data12(:,:,2)))) isgood = .false.

  ! 3D - append
  call dump_netcdf(Filename, data(:,:,1:size(data,3)/2))
  call dump_netcdf(Filename, data(:,:,size(data,3)/2+1:), append=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  allocate(data11(dimlen(1),dimlen(2),dimlen(3)))
  call Get_NcVar(Filename,Varname,data11)
  if (any(notequal(data,data11))) isgood = .false.

  ! 5D - large file support
  call dump_netcdf(Filename, data9, lfs=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data10)
  if (any(notequal(data9,data10))) isgood = .false.

  ! 5D - netcdf4, deflate=default
  call dump_netcdf(Filename, data9, netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data10)
  if (any(notequal(data9,data10))) isgood = .false.

  ! 5D - netcdf4, deflate=0
  call dump_netcdf(Filename, data9, netcdf4=.true., deflate_level=0)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data10)
  if (any(notequal(data9,data10))) isgood = .false.

  ! 5D - netcdf4, deflate=9
  call dump_netcdf(Filename, data9, netcdf4=.true., deflate_level=9)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data10)
  if (any(notequal(data9,data10))) isgood = .false.


  ! 1D - netcdf4
  call dump_netcdf(Filename, data7, netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data8)
  if (any(notequal(data7,data8))) isgood = .false.
  ! 2D - netcdf4
  call dump_netcdf(Filename, data5, netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data6)
  if (any(notequal(data5,data6))) isgood = .false.
  ! 3D - netcdf4
  call dump_netcdf(Filename, data, netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data2)
  if (any(notequal(data,data2))) isgood = .false.
  ! 4D - netcdf4
  call dump_netcdf(Filename, data3, netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data4)
  if (any(notequal(data3,data4))) isgood = .false.
  ! 5D - netcdf4
  call dump_netcdf(Filename, data9, netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data10)
  if (any(notequal(data9,data10))) isgood = .false.
  ! 3D - dp - netcdf4
  call dump_netcdf(Filename, real(data,dp), netcdf4=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,ddata2)
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
  if (any(notequal(data7,data13(:,1)))) isgood = .false.
  if (any(notequal(data7,data13(:,2)))) isgood = .false.
  if (any(notequal(data7,data13(:,3)))) isgood = .false.
  ! 2D - append - netcdf4
  call dump_netcdf(Filename, data5, netcdf4=.true.)
  call dump_netcdf(Filename, data5, append=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data12)
  if (any(notequal(data5,data12(:,:,1)))) isgood = .false.
  if (any(notequal(data5,data12(:,:,2)))) isgood = .false.
  ! 3D - append - netcdf4
  call dump_netcdf(Filename, data(:,:,1:size(data,3)/2), netcdf4=.true.)
  call dump_netcdf(Filename, data(:,:,size(data,3)/2+1:), append=.true.)
  dimlen = Get_NcDim(Filename,Varname)
  call Get_NcVar(Filename,Varname,data11)
  if (any(notequal(data,data11))) isgood = .false.

  
  if (isgood) then
     write(*,*) 'mo_ncdump ', color('o.k.', c_green)
  else
     write(*,*) 'mo_ncdump ', color('failed!', c_red)
  endif
  
end program main

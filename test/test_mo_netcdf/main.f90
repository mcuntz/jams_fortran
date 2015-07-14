! ---------------------------------------------------------------------------------------
!
! mo_netcdf test program
!
! author: David Schaefer
!
! created: June 2015
!
! ----------------------------------------------------------------------------------------

program test_mo_netcdf

  use mo_kind  , only : i4, sp, dp
  use mo_netcdf, only : NcDataset, NcDimension, NcVariable
  use mo_utils , only : equal
  
  implicit none

  logical                  :: correct
  integer(i4),  parameter  :: ntime=30, ny=15, nx=25, nadd=12

  character(*), parameter  :: fname="netcdf_make_check_test_file"
  character(*), parameter  :: vname_time="time", vname_lat="lon", vname_lon="lat", vname_data="data"
  character(64)            :: wavalue, ravalue
  
  type(NcDataset)          :: nc
  type(NcDimension)        :: dim_x, dim_y, dim_time
  type(NcVariable)         :: var_lon, var_lat, var_time, var_data

  integer(i4)              :: wtime(ntime),  i
  integer(i4), allocatable :: rtime(:)
  real(sp)                 :: wlat(nx,ny), wlon(nx,ny)
  real(sp),    allocatable :: rlat(:,:),rlon(:,:)
  real(dp)                 :: wdata(nx,ny,ntime), wfvalue, rfvalue
  real(dp),    allocatable :: rdata(:,:,:)

  ! ------------------------------------------
  ! 0. Initialization ...
  ! ------------------------------------------
  
  ! assume the program is correct up to this point ...
  correct = .true.
  
  ! generate some dummy data
  wfvalue = -9999._dp
  wavalue = "David Schaefer"
  wtime = (/ (i, i = 1, ntime ) /)
  do i = 1, nx
     wlon(i,:) = real(i-.5_sp,sp)
  end do
  do i = 1, ny
     wlat(:,i) = real(i-.5_sp,sp)
  end do  
  do i = 1, ntime
     wdata(:,:,i) = (wlon + wlat) * i
  end do

  ! --------------------------------------------
  ! 1. Create a file and dynamically append data
  ! --------------------------------------------
  
  ! 1.1 create a file
  nc = NcDataset(fname,"w")

  ! create dimensions
  dim_x    = nc%setDimension("x",    nx) ! lenght < 0 -> unlimited dimension
  dim_y    = nc%setDimension("y",    ny)
  dim_time = nc%setDimension("time", -1)

  ! create variables
  var_time = nc%setVariable("time", "i32", (/dim_time/))
  var_lat  = nc%setVariable("lat" , "f32", (/dim_x,dim_y/))
  var_lon  = nc%setVariable("lon" , "f32", (/dim_x,dim_y/))
  var_data = nc%setVariable("data", "f64", (/dim_x, dim_y, dim_time/))

  ! set fill value
  call var_data%setFillValue(wfvalue)
  
  ! write data of static variables
  call var_lat%setData(wlat)
  call var_lon%setData(wlon)
    
  ! append data within a loop
  do i = 1,ntime
     call var_time%setData(wtime(i),     start=(/i/))
     call var_data%setData(wdata(:,:,i), start=(/1,1,i/))
  end do

  ! add some variable attributes
  call var_data%setAttribute("units",   "mm/d")
  call var_data%setAttribute("scaling", .1_dp)

  ! add global attributes
  call nc%setAttribute("Author", wavalue)
  call nc%setAttribute("Year",   2099_i4)

  ! close the file
  call nc%close()

  
  ! 1.2. Read the written data
  ! --------------------------

  ! open dataset
  nc = NcDataset(fname,"r")

  ! acces the variable
  var_time = nc%getVariable(vname_time)
  var_lat  = nc%getVariable(vname_lat)
  var_lat  = nc%getVariable(vname_lon)
  var_data = nc%getVariable(vname_data)
  
  ! read the data
  call var_time%getData(rtime)
  call var_lat%getData(rlat)
  call var_lon%getData(rlon)
  call var_data%getData(rdata)

  ! read the fill value
  call var_data%getFillValue(rfvalue)

  ! read a global attribute
  call nc%getAttribute("Author", ravalue)

  ! close dataset
  call nc%close()

  ! 1.3 Check
  ! ---------
  if (.not. all(rtime .eq. wtime))    correct = .false.
  if (.not. all(equal(rlat, wlat)))   correct = .false.
  if (.not. all(equal(rdata, wdata))) correct = .false.
  if (.not. equal(rfvalue, rfvalue))  correct = .false.
  if (.not. (ravalue .eq. wavalue))   correct = .false.

  ! --------------------------------
  ! 2. Append to an existing dataset
  ! --------------------------------
  !
  ! 2.1 Write data
  ! --------------
  !
  ! open dataset
  nc = NcDataset(fname,"a")

  ! acces variable
  var_time = nc%getVariable(vname_time)  
  var_data = nc%getVariable(vname_data)

  ! append data within a loop
  do i = ntime, ntime+nadd
     call var_time%setData(i,    start=(/1,1,i/))
     call var_data%setData(wlon, start=(/1,1,i/))
  end do

  ! close dataset
  call nc%close()
  
  ! 2.2 Read the appended data
  ! ---------------------------
  !
  ! open dataset
  nc = NcDataset(fname,"r")

  ! acces the variable
  var_data = nc%getVariable(vname_data)

  ! read the appended chunk of data
  call var_data%getData(rdata, start=(/1,1,ntime/), count=(/nx,ny,nadd/))

  ! close dataset
  call nc%close()

  ! 2.3 Check
  ! ---------
  do i=1,nadd     
     if (.not. all(equal(rdata(:,:,i), real(wlon,dp)))) then
        correct = .false.
     end if
  end do
  
  
  ! ----------------------------------
  ! 3. Dump some data - the short form
  ! ----------------------------------

  ! 3.1 Write data
  ! --------------
  ! open a file
  nc = NcDataset(fname,"w")

  ! create variable and dimensions
  var_data = nc%setVariable(vname_data, "f64", (/ &
       nc%setDimension("x", nx),              &
       nc%setDimension("y", ny),              &
       nc%setDimension("time", ntime) /)      &
       )
  
  ! write data
  call var_data%setData(wdata)

  ! write data
  call var_data%setData(wdata)

  ! close the file
  call nc%close()

  
  ! 3.2 Read the dumped data
  ! --------------------------
  ! open dataset
  nc = NcDataset(fname,"r")

  ! ! acces the variable
  var_data = nc%getVariable(vname_data)
  
  ! read the data
  call var_data%getData(rdata)

  ! close dataset
  call nc%close()

  ! 3.3 check
  ! ---------
  if (.not. all(equal(rdata, wdata))) correct = .false.
    
  ! --------------------------------------------------------------------------------------
  ! The moment of truth ...
  if (correct) then
     print*, "mo_netcdf is o.k."
  else
     print*, "mo_netcdf failed."
  endif

end program test_mo_netcdf

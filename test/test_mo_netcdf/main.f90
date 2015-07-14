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

  logical                 :: correct
  integer(i4), parameter  :: time=30, y=15, x=25, sup=10

  character(*), parameter :: fname="netcdf_make_check_test_file", dname_time="time", dname_y="y", dname_x="x"
  character(*), parameter :: vname_1d="time", vname_2d="xy", vname_3d="data"

  type(NcDataset)         :: nc
  type(NcDimension)       :: dim1, dim2, dim3
  type(NcVariable)        :: var_1d, var_2d, var_3d

  character(64)           :: wavalue, ravalue
  integer(i4)             :: wdata_1d(time), i
  real(sp)                :: wdata_2d(x,y)
  real(dp)                :: wdata_3d(x,y,time), wfvalue, rfvalue

  integer(i4), allocatable  :: rdata_1d(:)
  real(sp)   , allocatable  :: rdata_2d(:,:)
  real(dp)   , allocatable  :: rdata_3d(:,:,:)

  ! assume the program is correct up to this point ...
  correct = .true.
  
  ! generate some dummy data
  wfvalue = -9999_dp
  wavalue = "David Schaefer"
  wdata_1d = (/ (i, i = 1, time ) /)
  do i = 1,y
     wdata_2d(:,i) = real(i-.5_sp,sp)
  end do
  do i = 1, time
     wdata_3d(:,:,i) = real(i,dp)
  end do

  ! --------------------
  ! 1. Create a new file
  ! --------------------
  !
  ! 1.1 Write data
  ! --------------
  !
  ! open dataset
  nc = NcDataset(fname,"w")

  ! create dimensions
  dim1 = nc%setDimension(dname_x, x) ! lenght < 0 -> unlimited dimension
  dim2 = nc%setDimension(dname_y, y)
  dim3 = nc%setDimension(dname_time, -1)

  ! create variables
  var_1d = nc%setVariable(vname_1d, "i32", (/dim3/))
  var_2d = nc%setVariable(vname_2d, "f32", (/dim1,dim2/))
  var_3d = nc%setVariable(vname_3d, "f64", (/dim1,dim2,dim3/))

  ! set the variable fill value -> needs to be done BEFORE data is written
  call var_3d%setFillValue(wfvalue)

  ! write data
  call var_1d%setData(wdata_1d)
  call var_2d%setData(wdata_2d)
  call var_3d%setData(wdata_3d)

  ! ! create some variable attributes
  call var_1d%setAttribute("units", "days since the beginning of time")
  call var_3d%setAttribute("scaling_factor", .1)
  
  ! create some global attributes
  call nc%setAttribute("Author", wavalue)
  call nc%setAttribute("Iterations", 42)
  
  ! Don't forget to close the file, support for destructors 
  ! is pretty bad in current compilers ...
  call nc%close()
  
  ! 1.2. Read the written data
  ! --------------------------

  ! open dataset
  nc = NcDataset(fname,"r")

  ! acces the variable
  var_1d = nc%getVariable(vname_1d)
  var_2d = nc%getVariable(vname_2d)
  var_3d = nc%getVariable(vname_3d)
  
  ! read the data
  call var_1d%getData(rdata_1d)
  call var_2d%getData(rdata_2d)
  call var_3d%getData(rdata_3d)

  ! read the fill value
  call var_3d%getFillValue(rfvalue)

  ! read a global attribute
  call nc%getAttribute("Author", ravalue)

  ! close dataset
  call nc%close()

  ! 1.3 Check
  ! ---------
  if (.not. all(rdata_1d .eq. wdata_1d))    correct = .false.
  if (.not. all(equal(rdata_2d, wdata_2d))) correct = .false.
  if (.not. all(equal(rdata_3d, wdata_3d))) correct = .false.
  if (.not. equal(rfvalue, rfvalue))        correct = .false.
  if (.not. (ravalue .eq. wavalue))         correct = .false.

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
  var_3d = nc%getVariable(vname_3d)

  ! append one timestep
  call var_3d%setData(wdata_2d+1,start=(/1,1,time+1/))

  ! append data within a loop
  do i = 2,sup
     call var_3d%setData(wdata_2d+i,start=(/1,1,time+i/))
  end do

  ! close dataset
  call nc%close()
  
  ! 2.2 Read the appended data
  ! ---------------------------
  !
  ! open dataset
  nc = NcDataset(fname,"r")

  ! acces the variable
  var_3d = nc%getVariable(vname_3d)

  ! read the data
  call var_3d%getData(rdata_3d,start=(/1,1,time+1/),count=(/x,y,sup/))

  ! close dataset
  call nc%close()

  ! 2.3 Check
  ! ---------
  do i=1,sup
     if (.not. all(equal(rdata_3d(:,:,i),real(wdata_2d+i,dp)))) then
        correct = .false.
     end if
  end do
  
  ! --------------------------------
  ! 3. Write into existing variables
  ! --------------------------------
  !
  ! 3.1 Write data
  ! --------------
  !
  ! open dataset
  nc = NcDataset(fname,"a")
  
  ! acces the variables
  var_1d = nc%getVariable(vname_1d)
  var_3d = nc%getVariable(vname_3d)

  ! write
  call var_1d%setData((/42,42,42,42/),stride=(/3/))
  
  ! close dataset
  call nc%close()

  ! 3.2 Read data
  ! --------------
  !
  ! open dataset
  nc = NcDataset(fname,"r")
  
  ! acces the variables
  var_1d = nc%getVariable(vname_1d)

  ! read data
  call var_1d%getData(rdata_1d,stride=(/3/),count=(/4/))
  
  ! close dataset
  call nc%close()

  ! 3.3 Check
  ! ---------
  if (.not. all(rdata_1d .eq. (/42,42,42,42/))) correct = .false.  

  
  ! --------------------------------------------------------------------------------------
  ! The moment of truth ...
  if (correct) then
     print*, "mo_netcdf is o.k."
  else
     print*, "mo_netcdf failed."
  endif

end program test_mo_netcdf

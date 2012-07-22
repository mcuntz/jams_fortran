Module mo_setnc
  !
  ! This module sets the netcdf structure V, which means it sets the attributes,
  ! the dimensions and the data arrays
  !
  ! number precision
  use mo_kind, only: i4, sp
  !
  ! use user defined variables
  use mo_MainVar, only: lat, lon
  !
  ! use netcdfstruct V
  use mo_NCWrite, only: V, Dnc, Gatt, ndims, nVars
  !
  ! use constants of netcdf library
  use netcdf, only: NF90_CHAR, NF90_FLOAT, NF90_UNLIMITED, NF90_DOUBLE
  !
  implicit none
  !
  private

  !
  public :: setnc
  !
  ! HISTORY
  !     WRITTEN,  LUIS SAMANIEGO, Feb 2011
  !     MODIFIED, STEPHAN THOBER, NOV 2011, COMMENTED AND GENERALIZED
  !
contains
  !
  subroutine setnc
    !
    implicit none
    !
    ! local variables  
    integer(i4) :: i
    !
    ! define parameters
    nDims  = 3                          ! nr. dim types
    nVars  = 4                          ! total nr. var to print (including dimensions)
    !
    ! allocate arrays
    allocate ( Dnc(nDims)  )
    allocate ( V(nVars)    )
    !
    ! define dimensions --------------------------------------------------------
    Dnc(1)%name      = "easting"
    Dnc(1)%len       = 28
    !
    Dnc(2)%name      = "northing"
    Dnc(2)%len       = 36
    !
    Dnc(3)%name      = "time"
    Dnc(3)%len       = NF90_UNLIMITED
    
    !
    ! DIMENSION VARIABLES ------------------------------------------------------
    
    ! EASTING
    i                =  1
    V(i)%name        =  Dnc(i)%name
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  0
    V(i)%nSubs       =  0
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    V(i)%wFlag       = .true.
    !
    ! pointer to actual data      
    V(i)%G2_d        => Lon
    !
    ! attributes (other possibilities: add_offset, valid_min, valid_max)  
    V(i)%nAtt          = 3
    !
    V(i)%att(1)%name   = "units"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values  = "m"
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "x-coordinate in cartesian coordinates GK4"   
    !
    V(i)%att(3)%name   = "valid_range"
    V(i)%att(3)%xType  = NF90_FLOAT
    V(i)%att(3)%nValues= 2
    !
    write( V(i)%att(3)%values, '(2f15.2)')  Lon(1,1), Lon(size(Lon,1),size(Lon,2))
    !
    ! NORTHING
    i                =  2
    V(i)%name        =  Dnc(i)%name
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  0
    V(i)%nSubs       =  0
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    V(i)%wFlag       = .true.
    !
    ! pointer      
    V(i)%G2_d        => Lat
    !
    ! attributes
    V(i)%nAtt          = 3
    !
    V(i)%att(1)%name   = "units"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values  = "m"
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "y-coordinate in cartesian coordinates GK4"
    !
    V(i)%att(3)%name   = "valid_range"
    V(i)%att(3)%xType  = NF90_FLOAT
    V(i)%att(3)%nValues= 2
    write( V(i)%att(3)%values, '(2f15.2)')  Lat(1,1), Lat(size(Lat,1),size(Lat,2))
    !
    ! TIME
    i                =  3
    V(i)%name        =  Dnc(i)%name
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  0
    V(i)%nSubs       =  0
    V(i)%nDims       =  1
    V(i)%dimTypes    =  (/i,0,0,0,0/)
    V(i)%wFlag       = .true.
    !
    ! pointer will be attributed during run time
    !
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "units"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "time_unit_since"
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "time"
    !  
    ! FIELD VARIABLES ----------------------------------------------------------
    !
    ! pr
    i                =  4
    V(i)%name        =  "pr"
    V(i)%xType       =  NF90_FLOAT
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  3
    V(i)%dimTypes    =  (/1,2,3,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    !                   during running time 
    ! attributes 
    V(i)%nAtt          = 7
    !
    V(i)%att(1)%name   = "units"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "mm"
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "precipitation"
    !
    V(i)%att(3)%name   = "valid_range"
    V(i)%att(3)%xType  = NF90_FLOAT
    V(i)%att(3)%nValues= 2
    write( V(i)%att(3)%values, '(2f15.2)')  0._sp, 100._sp
    !
    V(i)%att(4)%name   = "scale_factor"
    V(i)%att(4)%xType  = NF90_FLOAT
    V(i)%att(4)%nValues= 1
    V(i)%att(4)%values = "1."
    !
    V(i)%att(5)%name   = "_FillValue"
    V(i)%att(5)%xType  = NF90_FLOAT
    V(i)%att(5)%nValues= 1
    V(i)%att(5)%values = "-9999."
    !
    V(i)%att(6)%name   = "missing_value"
    V(i)%att(6)%xType  = NF90_FLOAT
    V(i)%att(6)%nValues= 1
    V(i)%att(6)%values = "-9999."
    !
    V(i)%att(7)%name   = "coordinates"
    V(i)%att(7)%xType  = NF90_CHAR
    V(i)%att(7)%nValues= 1
    V(i)%att(7)%values = "lon lat"
    !
    ! global attributes ---------------------------------------------------
    Gatt(1)%name       = "title"
    Gatt(1)%values     = "This is a test run!"
    !
    Gatt(2)%name       = "history"
    Gatt(2)%values     = "Stephan Thober, Nov 2011"
    !
  end subroutine  setnc
end Module mo_setnc

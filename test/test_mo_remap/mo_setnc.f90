module mo_setnc

  ! -----------------------------------------------------------------------
  !
  ! This module sets the netcdf structure V, which means it sets the attributes,
  ! the dimensions and the data arrays

  ! HISTORY
  !     Written,  Luis Samaniego, Feb 2011
  !     Modified, Stephan Thober, Nov 2011 - commented and generalized

  ! -----------------------------------------------------------------------

  implicit none

  private

  public :: setnc

  ! -----------------------------------------------------------------------

contains

  ! -----------------------------------------------------------------------

  subroutine setnc

    use mo_kind,    only: i4
    use mo_mainvar, only: dst_data, x, y, lat, lon
    use mo_ncwrite, only: V, Dnc, Gatt, ndims, nVars
    use netcdf,     only: NF90_CHAR, NF90_FLOAT, NF90_INT, NF90_UNLIMITED, NF90_DOUBLE

    implicit none

    integer(i4) :: i

    ! define parameters
    nDims  = 3                          ! nr. dim types
    nVars  = 3 + nDims                  ! total nr. var to print (including dimensions)

    ! allocate arrays
    allocate(Dnc(nDims))
    allocate(V(nVars))

    ! define dimensions --------------------------------------------------------
    Dnc(1)%name      = "x"
    Dnc(1)%len       = size(dst_data,1)
    !
    Dnc(2)%name      = "y"
    Dnc(2)%len       = size(dst_data,2)
    !
    Dnc(3)%name      = "time"
    Dnc(3)%len       = NF90_UNLIMITED
    !
    allocate(x(size(dst_data,1)))
    forall(i=1:size(dst_data,1)) x(i) = i
    allocate(y(size(dst_data,2)))
    forall(i=1:size(dst_data,2)) y(i) = i

    ! DIMENSION VARIABLES ------------------------------------------------------
    ! Number of Cells
    i                =  1
    V(i)%name        =  Dnc(i)%name
    V(i)%xType       =  NF90_INT
    V(i)%nLvls       =  0
    V(i)%nSubs       =  0
    V(i)%nDims       =  1
    V(i)%dimTypes    =  (/1,0,0,0,0/)
    V(i)%wFlag       = .true.
    ! pointer to actual data
    V(i)%G1_i        => x
    ! attributes (other possibilities: add_offset, valid_min, valid_maxval)  
    V(i)%nAtt          = 1
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Number of cells in x direction"   

    ! NORTHING
    i                =  2
    V(i)%name        =  Dnc(i)%name
    V(i)%xType       =  NF90_INT
    V(i)%nLvls       =  0
    V(i)%nSubs       =  0
    V(i)%nDims       =  1
    V(i)%dimTypes    =  (/2,0,0,0,0/)
    V(i)%wFlag       = .true.
    ! pointer    
    V(i)%G1_i        => y
    ! attributes
    V(i)%nAtt          = 1
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Number of cells in x direction"

    ! Time
    i                =  3
    V(i)%name        =  Dnc(i)%name
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  0
    V(i)%nSubs       =  0
    V(i)%nDims       =  1
    V(i)%dimTypes    =  (/3,0,0,0,0/)
    V(i)%wFlag       = .true.
    ! pointer will be attributed during runtime
    !V(i)%G1_d        => time
    ! attributes
    V(i)%nAtt          = 3
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Number of timesteps"
    !
    V(i)%att(2)%name   = "units"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "days since 0001-01-01 00:00:00"
    !
    V(i)%att(3)%name   = "calendar"
    V(i)%att(3)%xType  = NF90_CHAR
    V(i)%att(3)%nValues= 1
    V(i)%att(3)%values = "proleptic_gregorian"

    ! FIELD VARIABLES ----------------------------------------------------------
    i                =  4
    V(i)%name        =  "Z"
    V(i)%xType       =  NF90_FLOAT
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  3
    V(i)%dimTypes    =  (/1,2,3,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer will be attributed during runtime
    !V(i)%G3_f        => dst_data
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "reflectivity"
    !
    V(i)%att(2)%name   = "_FillValue"
    V(i)%att(2)%xType  = NF90_FLOAT
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "-9999."

    ! latitude
    i                =  5
    V(i)%name        =  "lat"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => lat
    ! attributes 
    V(i)%nAtt          = 1
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Latitude"

    ! longitude
    i                =  6
    V(i)%name        =  "lon"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! pointer  
    V(i)%G2_d        => lon
    ! attributes 
    V(i)%nAtt          = 1
    !
    V(i)%att(1)%name   = "long_name"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "Longitude"

    ! global attributes ---------------------------------------------------
    Gatt(1)%name       = "Title"
    Gatt(1)%values     = "P-values of Generated Gaussian Field"
    !
    Gatt(2)%name       = "history"
    Gatt(2)%values     = "Stephan Thober, Aug 2012"

  end subroutine setnc

end module mo_setnc

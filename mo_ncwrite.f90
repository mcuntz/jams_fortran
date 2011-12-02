module mo_NCinter

  ! This module provides routines for writing data in a nc file using the netcdf4 library.

  ! Please make sure you include the netcdf4 library in your makefile in order to use this module
  ! see this webpage for netcdf documentation:
  ! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html

  ! numerical precision
  use mo_kind, only: i4, sp, dp 
  
  ! functions and constants of netcdf4 library
  use netcdf, only: NF90_FLOAT, NF90_CHAR, NF90_UNLIMITED, NF90_CLOBBER, NF90_INT, nf90_def_var, &
                    nf90_put_att, nf90_enddef, nf90_create, nf90_def_dim, nf90_put_var, nf90_close, &
                    nf90_noerr, nf90_strerror
  !
  private
  !
  public :: set_NcVar     ! set the variables and attributes which shall be written in the file
  public :: create_netcdf ! write the variables and attributes in the file
  public :: Write_Static  ! write the actual data in the file statically ( at once )
  public :: Write_Dynamic ! write the actual data in the file dynamically ( one record after the other )
  public :: close_netcdf  ! close the netcdf file
  !
  ! Definition of PRIVATE Variables ------------------------------------------------------------------------
  !
  ! definition of parameters
  integer(i4), parameter                    :: nMaxDim = 5         ! nr. max dimensions
  integer(i4), parameter                    :: nMaxAtt = 20        ! nr. max attributes
  integer(i4), parameter                    :: maxLen  = 256       ! nr. string lenght
  integer(i4), parameter                    :: nGAtt   = 2         ! nr. global attributes
  integer(i4), parameter                    :: nAttDim = 2         ! dim array of attribute values
  !
  ! definition of variables
  integer(i4)                               :: nVars               ! nr. variables
  integer(i4)                               :: nDims               ! nr. dimensions 
  integer(i4)                               :: nLats               ! latitude  => nRows
  integer(i4)                               :: nLons               ! longitude => nCols
  integer(i4)                               :: nLvls               ! nr. levels
  integer(i4)                               :: nSubs               ! nr. of sub units / tiles/ fractions
  integer(i4)                               :: nRecs               ! time intervals
  !
  type attribute
    character (len=maxLen)                  :: name                ! attribute name
    integer(i4)                             :: xType               ! attribute of the values
    integer(i4)                             :: nValues             ! number of attributes       
    character (len=maxLen)                  :: values               ! numbers or "characters" separed by spaces
  end type attribute
  !
  type variable
     character (len=maxLen)                 :: name                ! short name
     integer(i4)                            :: xType               ! NF90 var. type
     integer(i4)                            :: nLvls               ! number of levels  
     integer(i4)                            :: nSubs               ! number of subparts     
     logical                                :: unlimited           ! time limited
     integer(i4)                            :: variD               ! Id
     integer(i4)                            :: nDims               ! field dimension
     integer(i4), dimension(nMaxDim)        :: dimIds              ! passing var. dimensions
     integer(i4), dimension(nMaxDim)        :: dimTypes            ! type of dimensions
     integer(i4)                            :: nAtt                ! nr. attributes     
     type(attribute), dimension(:), allocatable    :: att                 ! var. attributes   
     integer(i4), dimension(nMaxDim)        :: start               ! starting indices for netCDF
     integer(i4), dimension(nMaxDim)        :: count               ! counter          for netCDF
  end type variable
  type(variable), dimension(:), allocatable :: V                   ! variable list
  !
  !type(attribute), dimension(nGAtt)         :: globalAtt           ! global attributes for netCDF
  !
  type dims
     character (len=maxLen)                 :: name                ! dim. name
     integer(i4)                            :: len                 ! dim. lenght, undefined time => NF90_UNLIMITED
     integer(i4)                            :: dimId               ! dim. Id
  end type dims
  type (dims), dimension(:), allocatable    :: Dnc                 ! dimesions netCDF e.g. lon, lat, level, time
  !
  interface Write_Static
     module procedure Write_Static_i_2d, write_Static_i_3d, write_static_i_4d, &
                      Write_Static_s_2d, write_Static_s_3d, write_static_s_4d, &
                      Write_Static_d_2d, Write_Static_d_3d, Write_Static_d_4d
  end interface Write_Static
  !
  interface Write_Dynamic
     module procedure Write_Dynamic_i_0d, Write_Dynamic_i_1d, Write_Dynamic_i_2d, &
                      write_Dynamic_i_3d, write_Dynamic_i_4d
     module procedure Write_Dynamic_s_0d, Write_Dynamic_s_1d, Write_Dynamic_s_2d, &
                      write_Dynamic_s_3d, write_Dynamic_s_4d
     module procedure Write_Dynamic_d_0d, Write_Dynamic_d_1d, Write_Dynamic_d_2d, &
                      write_Dynamic_d_3d, write_Dynamic_d_4d
  end interface Write_Dynamic
  !
  contains
!******************************************************************************

!    NAME
!        set_NcVar

!    PURPOSE
!        This subroutine reads the attributes of the dimensions and variables
!        and stores these in the private data structure of this module. The
!        variables have to be specified in the variable VarList. Furthermore
!        the Path of the folder, where the files with the attributes are
!        written, have to be specified.

!    CALLING SEQUENCE
!        call set_NcVar(VariableList, AttributePath)

!    INTENT(IN)
!        character(256), dimension(:) :: VariableList      List of the variable names

!    INTENT(IN)
!        character(256)               :: AttributePath     Path of the attribute files

!    RESTRICTIONS
!        Can read attributes which are written in files in the <AttributePath>.
!        The following files have to be given:
!          - dim.txt
!                This file contains the dimensions to be written. It has to be of
!                the following form.
!             >>>Number_of_Dimensions: 3
!                <dim1_name> <dim1_len>
!                <dim2_name> <dim2_len>
!                ...
!                <dimn_name> <dimn_len>
!                <<<EOF
!                The dimension names as well as the dimension lengths do not contain
!                any spaces.
!          For each dimension and each variable a file containing the attributes has to 
!          be specified. These have to look like this:
!          - <name>_att.txt
!                <name> can be either Variable name or dimension name. These files have to
!                be written like this
!             >>>xtype= (NF90_FLOAT or NF90_INT or NF90_CHAR)
!                nLvls= 1
!                nSubs= 1
!                nDims= 3
!                dimTypes= 1 2 3 0 0 
!                Number_of_Attributes: <Att_Number>
!                name= <Att_name>
!                xType= (NF90_FLOAT or NF90_INT or NF90_CHAR
!                nValues= 1
!                values= m
!                name= long_name
!xtype= NF90_CHAR
!nValues= 1
!values= "x-coordinate in cartesian coordinates GK4"   

!    LITERATURE
!        http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html

!    HISTORY
!        Written,  Luis Samaniego, Feb 2011
!        Modified, Stephan Thober, Nov 2011 - restructured
!
!******************************************************************************
subroutine set_NcVar(VarList, AttPath)
  !
  implicit none
  !
  ! input variables
  character(256), dimension(:), intent(in)     :: VarList   ! List of  Variable names
  character(256),               intent(in)     :: AttPath
  !
  ! local variables  
  character(256)                               :: dummy, dummy2
  character(256)                               :: filename
  integer(i4)                                  :: i, j
  integer(i4)                                  :: nDims
  !
  ! read dimensions
  open (unit = 110, file=trim(AttPath)//"dim.txt", status="old", &
                                    action="read")
  !
  read(110,*) dummy, nDims
  !
  allocate(Dnc(ndims))
  !
  do i = 1, ndims
     !
     read(110,*) Dnc(i)%name, Dnc(i)%len
     if (Dnc(i)%len == -1) Dnc(i)%len = NF90_UNLIMITED
     !
  end do
  !
  close(110)
  !
  ! read and set all Attributes
  allocate ( V(size(Dnc,1)+size(VarList,1)))
  !
  ! 1st read dimension attributes
  Dimloop: do i = 1, size(Dnc,1)
     !
     V(i)%name = trim(Dnc(i)%name)
     write(filename,*) trim(AttPath)//trim(V(i)%name)//'_att.txt'
     open(unit = 120, file = filename, status = "old", &
                                       action = "read" )
     !
     read(120,*) dummy, dummy2
     select case (trim(dummy2))
        case ('NF90_FLOAT')
           V(i)%xType = NF90_FLOAT
        case ('NF90_CHAR')
           V(i)%xType = NF90_CHAR
     end select
     read(120,*) dummy, V(i)%nLvls
     read(120,*) dummy, V(i)%nSubs
     read(120,*) dummy, V(i)%nDims
     read(120,*) dummy, V(i)%dimTypes
     !
     read(120,*) dummy, V(i)%natt
     allocate(V(i)%att(V(i)%natt))
     !
     ! read attributes of ith Variable
     do j = 1, V(i)%natt
        !
        read(120,*) dummy, V(i)%att(j)%name
        read(120,*) dummy, dummy2
        select case (trim(dummy2))
        case ('NF90_FLOAT')
           V(i)%att(j)%xType = NF90_FLOAT
        case ('NF90_CHAR')
           V(i)%att(j)%xType = NF90_CHAR
        end select
        read(120,*) dummy, V(i)%att(j)%nValues
        read(120,*) dummy, V(i)%att(j)%values
        !
     end do
     !
     close(120)
     !
     V(i)%nLvls       =  0
     V(i)%nSubs       =  0
     V(i)%nDims       =  1
     V(i)%dimTypes    =  (/i,0,0,0,0/)
     !
  end do Dimloop
  !
  ! 2nd read variable attributes
  VarLoop: do i = size(Dnc,1) + 1, size(V,1)
     !
     ! get Attributes and properties of variable from file
     V(i)%name = trim(Varlist(i-size(Dnc,1)))
     write(filename,*) trim(AttPath)//trim(V(i)%name)//'_att.txt'
     open(unit = 120, file = filename, status = "old", &
                                       action = "read" )
     !
     read(120,*) dummy, dummy2
     select case (trim(dummy2))
        case ('NF90_FLOAT')
           V(i)%xType = NF90_FLOAT
        case ('NF90_CHAR')
           V(i)%xType = NF90_CHAR
     end select
     read(120,*) dummy, V(i)%nLvls
     read(120,*) dummy, V(i)%nSubs
     read(120,*) dummy, V(i)%nDims
     read(120,*) dummy, V(i)%dimTypes
     !
     read(120,*) dummy, V(i)%natt
     allocate(V(i)%att(V(i)%natt))
     !
     ! read attributes of ith Variable
     do j = 1, V(i)%natt
        !
        read(120,*) dummy, V(i)%att(j)%name
        read(120,*) dummy, dummy2
        select case (trim(dummy2))
        case ('NF90_FLOAT')
           V(i)%att(j)%xType = NF90_FLOAT
        case ('NF90_CHAR')
           V(i)%att(j)%xType = NF90_CHAR
        end select
        read(120,*) dummy, V(i)%att(j)%nValues
        read(120,*) dummy, V(i)%att(j)%values
        !
     end do
     !
     close(120)
     !
  end do VarLoop
  !
end subroutine  set_NcVar
!******************************************************************************
!  CREATE netCDF                                                              *
!******************************************************************************
subroutine create_netCDF(fname, ncid)
  !
  implicit none
  !
  ! netcdf related variables
  character(len=maxLen), intent(in)          :: fname
  integer(i4), intent(out)                   :: ncid
  integer(i4)                                :: i, j, k
  integer(i4), dimension (nAttDim)           :: att_INT  
  real(sp), dimension (nAttDim)              :: att_FLOAT
  character(len=maxLen), dimension (nAttDim) :: att_CHAR
  !
  ! 1  Create netCDF dataset: enter define mode     ->  get ncId  
  call check( nf90_create ( trim(fname), NF90_CLOBBER, ncId ))
  !
  ! 2  Define dimensions                                 -> get dimId
  do i=1, size(Dnc,1)
    call check( nf90_def_dim ( ncId, Dnc(i)%name , Dnc(i)%len, Dnc(i)%dimId ))
  end do
  !
  ! 3 Define dimids array, which is used to pass the dimids of the dimensions of
  ! the netCDF variables
  nVars = size(V,1)
  do i = 1, nVars
     print*, 'ndims',V(i)%nDims,'name',trim(V(i)%name)
    V(i)%unlimited = .false.
    V(i)%dimids  = 0
    V(i)%start   = 1
    V(i)%count   = 1
    do k = 1, V(i)%nDims
      if ( Dnc( V(i)%dimTypes(k) )%len  == NF90_UNLIMITED ) V(i)%unlimited = .true.
      V(i)%dimids(k) = Dnc( V(i)%dimTypes(k) )%dimId
    end do
    if ( V(i)%unlimited ) then
      ! set counts for unlimited files (time is always the last dimension
      if (V(i)%nDims == 1) cycle
      do k = 1, V(i)%nDims - 1
         V(i)%count(k)  = Dnc( V(i)%dimTypes(k) )%len
      end do
    end if
    print*, V(i)%dimids
  end do
  !
  ! 4 Define the netCDF variables and attributes                            -> get varId
  print*, 'set attributes...'
  do i=1, nVars
    call check( nf90_def_var ( ncId, V(i)%name, V(i)%xtype, V(i)%dimids(1:V(i)%nDims), V(i)%varId ))
    do k = 1, V(i)%nAtt
      select case ( V(i)%att(k)%xType )
        case (NF90_CHAR)
          ! read( V(i)%att(k)%values, *) ( att_CHAR(j), j =1, V(i)%att(k)%nValues )
          read( V(i)%att(k)%values, '(a)')  att_CHAR(1) 
          call check( nf90_put_att ( ncId, V(i)%varId, V(i)%att(k)%name, att_CHAR(1) ))     
        case (NF90_INT)
          read( V(i)%att(k)%values, *) ( att_INT(j), j =1, V(i)%att(k)%nValues )
          call check( nf90_put_att ( ncId, V(i)%varId, V(i)%att(k)%name, att_INT(1:V(i)%att(k)%nValues) ))
        case (NF90_FLOAT)
          read( V(i)%att(k)%values, *) ( att_FLOAT(j), j =1, V(i)%att(k)%nValues )
          call check( nf90_put_att ( ncId, V(i)%varId, V(i)%att(k)%name, att_FLOAT(1:V(i)%att(k)%nValues) ))
      end select
    end do
  end do
  !
  ! 5 Global attributes
  do k = 1, nGAtt
!    call check( nf90_put_att ( ncId, NF90_GLOBAL, globalAtt(k)%name, globalAtt(k)%values ))
  end do
  !
  ! 6 end definitions: leave define mode
  call check ( nf90_enddef ( ncId ))
  print *, "NetCDF file was created", ncId
  ! 
end subroutine create_netCDF

!******************************************************************************
!  WRITE STATIC                                                               *
!                                                                             *
!  NOTES: 1) netCDF file must be on *** data mode ***                         *
!         2) any numeric data type is allowed but NOT character               *
!         3) writes the whole array                                           *
!                                                                             *
!******************************************************************************
subroutine write_static_i_2d(ncId, i, data)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)                :: ncId
  integer(i4), intent(in)                :: i    ! ~ Var(i)
  integer(i4), dimension(:,:),intent(in) :: data
 
  ! NOTES: 1) netCDF file must be on *** data mode ***
  !
  ! write all static variables
  if ( V(i)%unlimited     ) stop 'ERROR*** This Variable is unlimited!!!'
  print *, "Var. ", trim(V(i)%name) ," is  static"
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data ))
  !
end subroutine write_static_i_2d
!
subroutine write_static_i_3d(ncId, i, data)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)                  :: ncId
  integer(i4), intent(in)                  :: i    ! ~ Var(i)
  integer(i4), dimension(:,:,:),intent(in) :: data
 
  ! NOTES: 1) netCDF file must be on *** data mode ***
  !
  ! write all static variables
  if ( V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  print *, "Var. ", trim(V(i)%name) ," is  static"
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data ))
  !
end subroutine write_static_i_3d
!
subroutine write_static_i_4d(ncId, i, data)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)                     :: ncId
  integer(i4), intent(in)                     :: i    ! ~ Var(i)
  integer(i4), dimension(:,:,:,:), intent(in) :: data
 
  ! NOTES: 1) netCDF file must be on *** data mode ***
  !
  ! write all static variables
  if ( V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  print *, "Var. ", trim(V(i)%name) ," is  static"
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data ))
  !
end subroutine write_static_i_4d
!
subroutine write_static_s_2d(ncId, i, data)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)             :: ncId
  integer(i4), intent(in)             :: i    ! ~ Var(i)
  real(sp), dimension(:,:),intent(in) :: data
 
  ! NOTES: 1) netCDF file must be on *** data mode ***
  !
  ! write all static variables
  if ( V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  print *, "Var. ", trim(V(i)%name) ," is  static"
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data ))
  !
end subroutine write_static_s_2d
!
subroutine write_static_s_3d(ncId, i, data)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)               :: ncId
  integer(i4), intent(in)               :: i    ! ~ Var(i)
  real(sp), dimension(:,:,:),intent(in) :: data
 
  ! NOTES: 1) netCDF file must be on *** data mode ***
  !
  ! write all static variables
  if ( V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  print *, "Var. ", trim(V(i)%name) ," is  static"
  print*,size(data,1), size(data,2), size(data,3)
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data ))
  !
end subroutine write_static_s_3d
!
subroutine write_static_s_4d(ncId, i, data)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)                  :: ncId
  integer(i4), intent(in)                  :: i    ! ~ Var(i)
  real(sp), dimension(:,:,:,:), intent(in) :: data
 
  ! NOTES: 1) netCDF file must be on *** data mode ***
  !
  ! write all static variables
  if ( V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  print *, "Var. ", trim(V(i)%name) ," is  static"
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data ))
  !
end subroutine write_static_s_4d
!
subroutine write_static_d_2d(ncId, i, data)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)             :: ncId
  integer(i4), intent(in)             :: i    ! ~ Var(i)
  real(dp), dimension(:,:),intent(in) :: data
 
  ! NOTES: 1) netCDF file must be on *** data mode ***
  !
  ! write all static variables
  if ( V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  print *, "Var. ", trim(V(i)%name) ," is  static"
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data ))
  !
end subroutine write_static_d_2d
!
subroutine write_static_d_3d(ncId, i, data)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)               :: ncId
  integer(i4), intent(in)               :: i    ! ~ Var(i)
  real(dp), dimension(:,:,:),intent(in) :: data
 
  ! NOTES: 1) netCDF file must be on *** data mode ***
  !
  ! write all static variables
  if ( V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  print *, "Var. ", trim(V(i)%name) ," is  static"
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data ))
  !
end subroutine write_static_d_3d
!
subroutine write_static_d_4d(ncId, i, data)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)                  :: ncId
  integer(i4), intent(in)                  :: i    ! ~ Var(i)
  real(dp), dimension(:,:,:,:), intent(in) :: data
 
  ! NOTES: 1) netCDF file must be on *** data mode ***
  !
  ! write all static variables
  if ( V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  print *, "Var. ", trim(V(i)%name) ," is  static"
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data ))
  !
end subroutine write_static_d_4d
!******************************************************************************
!  WRITE netCDF                                                               *
!                                                                             *
!  NOTES: 1) netCDF file must be on *** data mode ***                         *
!         2) any numeric data type is allowed but NOT character               *
!         3) writes only one specified record                                 *
!                                                                             *
!******************************************************************************
subroutine write_dynamic_i_0d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)           :: ncId
  integer(i4), intent(in)           :: i
  integer(i4), intent(in), optional :: iRec
  integer(i4), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start )) 
  !
end subroutine write_dynamic_i_0d
!
subroutine write_dynamic_i_1d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4),               intent(in)           :: ncId
  integer(i4),               intent(in)           :: i
  integer(i4),               intent(in), optional :: iRec
  integer(i4), dimension(:), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start, V(i)%count ))
  !
end subroutine write_dynamic_i_1d
!
subroutine write_dynamic_i_2d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4),                 intent(in)           :: ncId
  integer(i4),                 intent(in)           :: i
  integer(i4),                 intent(in), optional :: iRec
  integer(i4), dimension(:,:), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start, V(i)%count ))
  !
end subroutine write_dynamic_i_2d
!
subroutine write_dynamic_i_3d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4),                   intent(in)           :: ncId
  integer(i4),                   intent(in)           :: i
  integer(i4),                   intent(in), optional :: iRec
  integer(i4), dimension(:,:,:), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start, V(i)%count ))
  !
end subroutine write_dynamic_i_3d
!
subroutine write_dynamic_i_4d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4),                     intent(in)           :: ncId
  integer(i4),                     intent(in)           :: i
  integer(i4),                     intent(in), optional :: iRec
  integer(i4), dimension(:,:,:,:), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start, V(i)%count ))
  !
end subroutine write_dynamic_i_4d
!
subroutine write_dynamic_s_0d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)           :: ncId
  integer(i4), intent(in)           :: i
  integer(i4), intent(in), optional :: iRec
  real(sp),    intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start )) 
  !
end subroutine write_dynamic_s_0d
!
subroutine write_dynamic_s_1d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4),            intent(in)           :: ncId
  integer(i4),            intent(in)           :: i
  integer(i4),            intent(in), optional :: iRec
  real(sp), dimension(:), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start, V(i)%count ))
  !
end subroutine write_dynamic_s_1d
!
subroutine write_dynamic_s_2d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4),              intent(in)           :: ncId
  integer(i4),              intent(in)           :: i
  integer(i4),              intent(in), optional :: iRec
  real(sp), dimension(:,:), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start, V(i)%count ))
  !
end subroutine write_dynamic_s_2d
!
subroutine write_dynamic_s_3d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4),                intent(in)           :: ncId
  integer(i4),                intent(in)           :: i
  integer(i4),                intent(in), optional :: iRec
  real(sp), dimension(:,:,:), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start, V(i)%count ))
  !
end subroutine write_dynamic_s_3d
!
subroutine write_dynamic_s_4d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4),                  intent(in)           :: ncId
  integer(i4),                  intent(in)           :: i
  integer(i4),                  intent(in), optional :: iRec
  real(sp), dimension(:,:,:,:), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start, V(i)%count ))
  !
end subroutine write_dynamic_s_4d
!
subroutine write_dynamic_d_0d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4), intent(in)           :: ncId
  integer(i4), intent(in)           :: i
  integer(i4), intent(in), optional :: iRec
  real(dp),    intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start )) 
  !
end subroutine write_dynamic_d_0d
!
subroutine write_dynamic_d_1d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4),            intent(in)           :: ncId
  integer(i4),            intent(in)           :: i
  integer(i4),            intent(in), optional :: iRec
  real(dp), dimension(:), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start, V(i)%count ))
  !
end subroutine write_dynamic_d_1d
!
subroutine write_dynamic_d_2d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4),              intent(in)           :: ncId
  integer(i4),              intent(in)           :: i
  integer(i4),              intent(in), optional :: iRec
  real(dp), dimension(:,:), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start, V(i)%count ))
  !
end subroutine write_dynamic_d_2d
!
subroutine write_dynamic_d_3d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4),                intent(in)           :: ncId
  integer(i4),                intent(in)           :: i
  integer(i4),                intent(in), optional :: iRec
  real(dp), dimension(:,:,:), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start, V(i)%count ))
  !
end subroutine write_dynamic_d_3d
!
subroutine write_dynamic_d_4d(ncId, i, data, irec)
  !
  implicit none
  !
  ! netcdf related variables
  integer(i4),                  intent(in)           :: ncId
  integer(i4),                  intent(in)           :: i
  integer(i4),                  intent(in), optional :: iRec
  real(dp), dimension(:,:,:,:), intent(in)           :: data
  !
  ! set values for variables (one scalar or grid at a time)
  if ( .not. V(i)%unlimited     ) stop 'ERROR*** This Variable is not unlimited!!!'
  !
  if (iRec ==1) print *, "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
  V(i)%start ( V(i)%nDims ) = iRec
  !
  call check( nf90_put_var ( ncId,  V(i)%varId, data, V(i)%start, V(i)%count ))
  !
end subroutine write_dynamic_d_4d
!
!******************************************************************************
!  CLOSE netCDF                                                               *
!******************************************************************************
subroutine close_netCDF(ncId)
  !
  implicit none
  ! 
  integer(i4), intent(in)        :: ncId
  !
  ! close: save new netCDF dataset
  call check( nf90_close( ncId ))
  print *, "NetCDF file was saved"
  !
end subroutine close_netCDF

!******************************************************************************
!  ERROR STATUS                                                              *
!******************************************************************************
subroutine check(status)
  !
  implicit none
  !
  integer(i4), intent(in) :: status
  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop
  end if
end subroutine check  
!
end module mo_NCinter



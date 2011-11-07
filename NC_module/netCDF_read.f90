!******************************************************************************
!  NetCDF_Interface v1                                                        *
!  PURPOSE:                                                                   *
!           This subroutines reads data from a netcdf file. Data array is     *
!           NOT allowed to have more than 4 dimensions.                       *
!                                                                             *
!           Call Get_netcdf_Var(filename,varname,array)                       *
!                                                                             *
!           filename: intent(in)                                              *
!                     character(256)                                          *
!                     Filename of nc file to be read                          *
!                                                                             *
!           varname:  intent(in)                                              *
!                     character(256)                                          *
!                     variable name as EXACTLY specified in the file          *
!                                                                             *
!           array:    intent(inout)                                           *
!                     single or double precision                              *
!                     Array is an allocatable array (has not been allocated   *
!                     when calling). It can have 2, 3 or 4 dimensions. When   *
!                     the actual data has less dimensions than the array, the *
!                     dimensions of the array are filled with ones. The data  *
!                     will be stored in this array.                           *
!                                                                             *
!  MISSING TESTS: TESTS FOR DIFFERENT COMPILERS THAN INTEL11.1.075            *
!                                                                             *
!  FURTHER ENHANCEMENTS: PACK DATA                                            *
!                                                                             *
!  AUTHOR:  Stephan Thober UFZ 2011                                           *
!                                                                             *
!  UPDATES:                                                                   *
!           created  Thober     04.11.2011                                    *
!-----------------------------------------------------------------------------*
!  REQUIREMENTS AND SUGESTIONS                                                *
!           Matthias Cuntz UFZ                                                *
!******************************************************************************
module mo_netcdf_read
  !
  use mo_kind, only: i4, sp, dp
  use netcdf,  only: nf90_open, nf90_get_var, nf90_close, NF90_MAX_NAME , &
                     nf90_inquire, nf90_inq_varids, nf90_inquire_variable, &
                     nf90_inquire_dimension, NF90_NOWRITE
  ! netCDF configuration for mHM
  ! see: http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html
  !
  private
  !
  public :: Get_netcdf_Var
  !
  interface Get_netcdf_Var
     module procedure Get_netcdf_Var_3d_sp, Get_netcdf_Var_3d_dp, Get_netcdf_Var_2d_sp, &
                      Get_netcdf_Var_2d_dp, Get_netcdf_Var_4d_sp, Get_netcdf_Var_4d_dp
  end interface
  !
contains

! ------------------------------------------------------------------------------
! SUBROUTINE GET_NETCDF_VAR
! ------------------------------------------------------------------------------
subroutine Get_netcdf_Var_4d_sp(Filename, VarName, Data)
  !
  implicit none
  !
  character(256), intent(in)                :: Filename
  character(256), intent(in)                :: VarName   ! Variable name == Longname in netcdf dataset
  integer(i4)                               :: varid     ! id of variable to be read
  integer(i4)                               :: vartype   ! type of variable
  integer(i4), dimension(4)                 :: dl        ! length of dimensions
  real(sp), dimension(:,:,:,:), allocatable, &
                              intent(inout) :: Data      ! array where values should be stored
  !
  ! further Variables
  integer(i4) :: ncid   ! id of input stream
  integer(i4) :: status ! status of read stream
  !
  if (allocated(data)) stop 'ERROR*** data is already allocated. subroutine Get_netcdf_var'
  !
  ! Open NetCDF filename
  status = nf90_open(trim(Filename),NF90_NOWRITE, ncid)
  if ( status /= 0) stop 'ERROR*** nc file could not be opened!'
  !
  ! Inquire file and check if VarName exists in the dataset
  ! Get also the id and the length of the dimensions
  call get_Info(Varname,ncid,varid,dl,vartype)
  !
  ! check variable type ( 5 equals float type )
  if ( vartype /= 5 ) stop 'ERROR*** type of variable does not match argument type. subroutine get_Info'
  !
  ! get values
  allocate(data(dl(1),dl(2),dl(3),dl(4)))
  status = nf90_get_var(ncid, varid, data)
  if ( status /= 0) stop 'ERROR*** data could not be read!'
  !
  ! close File
  status = nf90_close(ncid)
  if ( status /= 0) stop 'ERROR*** nc file could not be closed!'
  !
end subroutine Get_netcdf_Var_4d_sp

! ------------------------------------------------------------------------------
! SUBROUTINE GET_NETCDF_VAR
! ------------------------------------------------------------------------------
subroutine Get_netcdf_Var_4d_dp(Filename, VarName, Data)
  !
  implicit none
  !
  character(256), intent(in)                :: Filename
  character(256), intent(in)                :: VarName   ! Variable name == Longname in netcdf dataset
  integer(i4)                               :: varid     ! id of variable to be read
  integer(i4)                               :: vartype   ! type of variable
  integer(i4), dimension(4)                 :: dl        ! length of dimensions
  real(dp), dimension(:,:,:,:), allocatable, &
                              intent(inout) :: Data      ! array where values should be stored
  !
  ! further Variables
  integer(i4) :: ncid   ! id of input stream
  integer(i4) :: status ! status of read stream
  !
  if (allocated(data)) stop 'ERROR*** data is already allocated. subroutine Get_netcdf_var'
  !
  ! Open NetCDF filename
  status = nf90_open(trim(Filename),NF90_NOWRITE, ncid)
  if ( status /= 0) stop 'ERROR*** nc file could not be opened!'
  !
  ! Inquire file and check if VarName exists in the dataset
  ! Get also the id and the length of the dimensions
  call get_Info(Varname,ncid,varid,dl,vartype)
  !
  ! check variable type ( 5 equals float type )
  if ( vartype /= 6 ) stop 'ERROR*** type of variable does not match argument type. subroutine get_Info'
  !
  ! get values
  allocate(data(dl(1),dl(2),dl(3),dl(4)))
  status = nf90_get_var(ncid, varid, data)
  if ( status /= 0) stop 'ERROR*** data could not be read!'
  !
  ! close File
  status = nf90_close(ncid)
  if ( status /= 0) stop 'ERROR*** nc file could not be closed!'
  !
end subroutine Get_netcdf_Var_4d_dp

! ------------------------------------------------------------------------------
! SUBROUTINE GET_NETCDF_VAR
! ------------------------------------------------------------------------------
subroutine Get_netcdf_Var_2d_dp(Filename, VarName, Data)
  !
  implicit none
  !
  character(256), intent(in)                :: Filename
  character(256), intent(in)                :: VarName   ! Variable name == Longname in netcdf dataset
  integer(i4)                               :: varid     ! id of variable to be read
  integer(i4)                               :: vartype   ! type of variable
  integer(i4), dimension(2)                 :: dl        ! length of dimensions
  real(dp), dimension(:,:), allocatable, &
                              intent(inout) :: Data      ! array where values should be stored
  !
  ! further Variables
  integer(i4) :: ncid   ! id of input stream
  integer(i4) :: status ! status of read stream
  !
  if (allocated(data)) stop 'ERROR*** data is already allocated. subroutine Get_netcdf_var'
  !
  ! Open NetCDF filename
  status = nf90_open(trim(Filename),NF90_NOWRITE, ncid)
  if ( status /= 0) stop 'ERROR*** nc file could not be opened!'
  !
  ! Inquire file and check if VarName exists in the dataset
  ! Get also the id and the length of the dimensions
  call get_Info(Varname,ncid,varid,dl,vartype)
  !
  ! check variable type ( 5 equals float type )
  if ( vartype /= 6 ) stop 'ERROR*** type of variable does not match argument type. subroutine get_Info'
  !
  ! get values
  allocate(data(dl(1),dl(2)))
  status = nf90_get_var(ncid, varid, data)
  if ( status /= 0) stop 'ERROR*** data could not be read!'
  !
  ! close File
  status = nf90_close(ncid)
  if ( status /= 0) stop 'ERROR*** nc file could not be closed!'
  !
end subroutine Get_netcdf_Var_2d_dp

! ------------------------------------------------------------------------------
! SUBROUTINE GET_NETCDF_VAR
! ------------------------------------------------------------------------------
subroutine Get_netcdf_Var_2d_sp(Filename, VarName, Data)
  !
  implicit none
  !
  character(256), intent(in)                :: Filename
  character(256), intent(in)                :: VarName   ! Variable name == Longname in netcdf dataset
  integer(i4)                               :: varid     ! id of variable to be read
  integer(i4)                               :: vartype   ! type of variable
  integer(i4), dimension(2)                 :: dl        ! length of dimensions
  real(sp), dimension(:,:), allocatable, &
                              intent(inout) :: Data      ! array where values should be stored
  !
  ! further Variables
  integer(i4) :: ncid   ! id of input stream
  integer(i4) :: status ! status of read stream
  !
  if (allocated(data)) stop 'ERROR*** data is already allocated. subroutine Get_netcdf_var'
  !
  ! Open NetCDF filename
  status = nf90_open(trim(Filename),NF90_NOWRITE, ncid)
  if ( status /= 0) stop 'ERROR*** nc file could not be opened!'
  !
  ! Inquire file and check if VarName exists in the dataset
  ! Get also the id and the length of the dimensions
  call get_Info(Varname,ncid,varid,dl,vartype)
  !
  ! check variable type ( 5 equals float type )
  if ( vartype /= 5 ) stop 'ERROR*** type of variable does not match argument type. subroutine get_Info'
  !
  ! get values
  allocate(data(dl(1),dl(2)))
  status = nf90_get_var(ncid, varid, data)
  if ( status /= 0) stop 'ERROR*** data could not be read!'
  !
  ! close File
  status = nf90_close(ncid)
  if ( status /= 0) stop 'ERROR*** nc file could not be closed!'
  !
end subroutine Get_netcdf_Var_2d_sp

! ------------------------------------------------------------------------------
! SUBROUTINE GET_NETCDF_VAR
! ------------------------------------------------------------------------------
subroutine Get_netcdf_Var_3d_dp(Filename, VarName, Data)
  !
  implicit none
  !
  character(256), intent(in)                :: Filename
  character(256), intent(in)                :: VarName   ! Variable name == Longname in netcdf dataset
  integer(i4)                               :: varid     ! id of variable to be read
  integer(i4)                               :: vartype   ! type of variable
  integer(i4), dimension(3)                 :: dl        ! length of dimensions
  real(dp), dimension(:,:,:), allocatable, &
                              intent(inout) :: Data      ! array where values should be stored
  !
  ! further Variables
  integer(i4) :: ncid   ! id of input stream
  integer(i4) :: status ! status of read stream
  !
  if (allocated(data)) stop 'ERROR*** data is already allocated. subroutine Get_netcdf_var'
  !
  ! Open NetCDF filename
  status = nf90_open(trim(Filename),NF90_NOWRITE, ncid)
  if ( status /= 0) stop 'ERROR*** nc file could not be opened!'
  !
  ! Inquire file and check if VarName exists in the dataset
  ! Get also the id and the length of the dimensions
  call get_Info(Varname,ncid,varid,dl,vartype)
  !
  ! check variable type ( 5 equals float type )
  if ( vartype /= 6 ) stop 'ERROR*** type of variable does not match argument type. subroutine get_Info'
  !
  ! get values
  allocate(data(dl(1),dl(2),dl(3)))
  status = nf90_get_var(ncid, varid, data)
  if ( status /= 0) stop 'ERROR*** data could not be read!'
  !
  ! close File
  status = nf90_close(ncid)
  if ( status /= 0) stop 'ERROR*** nc file could not be closed!'
  !
end subroutine Get_netcdf_Var_3d_dp

! ------------------------------------------------------------------------------
! SUBROUTINE GET_NETCDF_VAR
! ------------------------------------------------------------------------------
subroutine Get_netcdf_Var_3d_sp(Filename, VarName, Data)
  !
  implicit none
  !
  character(256), intent(in)                :: Filename
  character(256), intent(in)                :: VarName   ! Variable name == Longname in netcdf dataset
  integer(i4)                               :: varid     ! id of variable to be read
  integer(i4)                               :: vartype   ! type of variable
  integer(i4), dimension(3)                 :: dl        ! length of dimensions
  real(sp), dimension(:,:,:), allocatable, &
                              intent(inout) :: Data      ! array where values should be stored
  !
  ! further Variables
  integer(i4) :: ncid   ! id of input stream
  integer(i4) :: status ! status of read stream
  !
  if (allocated(data)) stop 'ERROR*** data is already allocated. subroutine Get_netcdf_var'
  !
  ! Open NetCDF filename
  status = nf90_open(trim(Filename),NF90_NOWRITE, ncid)
  if ( status /= 0) stop 'ERROR*** nc file could not be opened!'
  !
  ! Inquire file and check if VarName exists in the dataset
  ! Get also the id and the length of the dimensions
  call get_Info(Varname,ncid,varid,dl,vartype)
  !
  ! check variable type ( 5 equals float type )
  if ( vartype /= 5 ) stop 'ERROR*** type of variable does not match argument type. subroutine get_Info'
  !
  ! get values
  allocate(data(dl(1),dl(2),dl(3)))
  status = nf90_get_var(ncid, varid, data)
  if ( status /= 0) stop 'ERROR*** data could not be read!'
  !
  ! close File
  status = nf90_close(ncid)
  if ( status /= 0) stop 'ERROR*** nc file could not be closed!'
  !
end subroutine Get_netcdf_Var_3d_sp

! ------------------------------------------------------------------------------
!
! SUBROUTINE GET_INFO
!
! This function inquires the nc file.
! See: http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90/ -> NF90-INQUIRE
! for detailed information on output.
! Moreover, check if the Variable name does really exist in the file.
!
! ------------------------------------------------------------------------------
subroutine get_Info(Varname, ncid, varid, dl, kind)
  !
  implicit none
  !
  character(256),            intent(in)  :: Varname
  integer(i4),               intent(in)  :: ncid 
  integer(i4), dimension(:), intent(inout)  :: dl
  integer(i4),               intent(out) :: varid   ! variable id of data to be read
  integer(i4),               intent(out) :: kind    ! type of the variable
  !
  integer(i4)                            :: nVars   ! Number of Variables
  integer(i4), dimension(:), allocatable :: varids  ! Ids of all variables
  integer(i4), dimension(:), allocatable :: DimID   ! Id of dimension
  integer(i4)                            :: NumDims ! Number of Dimensions for specific variable
  character(NF90_MAX_NAME)               :: name    ! name of Variables in the file
  integer(i4)                            :: n       ! loop index
  integer(i4)                            :: status  ! control variable
  logical                                :: flag    ! flag if Variable is in file
  !
  status = nf90_inquire(ncid,nVariables = nVars)
  if ( status /= 0 ) stop 'ERROR*** Could NOT inquire Number of Variables'
  !
  allocate(varids(nVars))
  !
  ! get the varids
  status = nf90_inq_varids(ncid, nVars, varids)
  if ( status /= 0 ) stop 'ERROR*** Could NOT inquire IDs of Variables'
  !
  flag = .false.
  !
  varidloop: do n = 1, size(varids)
     !
     ! check whether Varname does exist in the file
     status = nf90_inquire_variable(ncid, varids(n), name, ndims = NumDims)
     if ( status /= 0 ) stop 'ERROR*** Could NOT inquire Number of Dimensions'
     !
     if (trim (Varname) == trim(name)) then
        flag  = .true.
        varid = varids(n)
        exit
     end if
     !
  end do varidloop
  !
  if ( .not. flag ) stop 'ERROR*** Varname does not exist in file. subroutine get_Info'
  !
  ! check consistency
  if ( NumDims > size(dl) ) &
       stop 'ERROR*** Dimension size of Variable is greater than dims of array. subroutine get_Info'
  ! if the data has less dimensions than the array, fill the dimension with ones
  if ( size(dl) > NumDims ) dl(NumDims+1:size(dl)) = 1
  !
  ! get the dimensions of this variable
  allocate(DimId(NumDims))
  status = nf90_inquire_variable(ncid, varids(n), name, xtype=kind, dimids=DimId)
  if ( status /= 0 ) stop 'ERROR*** Could NOT inquire Ids of Dimensions'
  !
  dimloop: do n = 1, NumDims
     !
     status = nf90_inquire_dimension(ncid, DimId(n), name, dl(n))
     if ( status /= 0 ) stop 'ERROR*** Could NOT inquire Length of Dimension'
     !
     write(*,'(a10,i1,a4,a<len(trim(name))>)') 'Dimension ',n,' is ', trim(name)
     write(*,'(a14,i5.5)') 'The Length is ',dl(n)
     !
  end do dimloop
  !
end subroutine get_Info

! ------------------------------------------------------------------------------
!
! SUBROUTINE GET_DIM
!
! This subroutine gets the dimension ids.
!
! ------------------------------------------------------------------------------
subroutine get_dim
end subroutine get_dim
!
end module mo_netcdf_read


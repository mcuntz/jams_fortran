module mo_NcRead
  
  ! This module provides subroutines for reading arrays from nc file using the netcdf4 library
  
  ! Please make sure you include the netcdf4 library in your makefile in order to use this module


!  AUTHOR:  Stephan Thober UFZ 2011                                           *
!                                                                             *
!  UPDATES:                                                                   *
!           created  Thober     04.11.2011                                    *

  ! numerical precision
  use mo_kind, only: i4, sp, dp
  !
  ! functions and constants of netcdf4 library
  use netcdf,  only: nf90_open, nf90_get_var, nf90_close, NF90_MAX_NAME , &
                     nf90_inquire, nf90_inq_varids, nf90_inquire_variable, &
                     nf90_inquire_dimension, NF90_NOWRITE
  !
  implicit none
  !
  private
  !
  public :: Get_NcVar
  !
  interface Get_NcVar
     module procedure Get_NcVar_3d_sp, Get_NcVar_3d_dp, Get_NcVar_2d_sp, &
                      Get_NcVar_2d_dp, Get_NcVar_4d_sp, Get_NcVar_4d_dp
  end interface
  !
contains
  !
  ! ------------------------------------------------------------------------------

  !    NAME
  !        Get_NcVar

  !    PURPOSE
  !        Reads a 2 - 4 dimensional array from a nc file given
  !        the variable name EXACTLY as specified in the file.
  !        When calling, the array has to have the allocatable attribute, but is
  !        not yet allocated. This will be done by Get_Nc_Var. If the dimension
  !        of the actual data is less than the ones of the array, then the dimension
  !        lengths of the array will be filled with ones.

  !    CALLING SEQUENCE
  !        call Get_NcVar(File, Variable_Name, array)

  !    INTENT(IN)
  !        character(256) :: File                                         Name of the nc file

  !    INTENT(IN)
  !        character(256) :: Variable_Name                                Name of the Variable in the nc file
  
  !    INTENT(INOUT)
  !        real(sp/dp), dimension(:,:[,:[,:]]), allocatable :: array      array where data will be read
  
  !    RESTRICTIONS
  !        Output array is a floating point of 2-4 dimensions.
  !        NOT yet tested for different compilers than intel11.1.075
  !        CANNOT read packed data

  !    EXAMPLE
  !        see test program in directory test_mo_NcRead

  !    LITERATURE
  !        http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html

  !    HISTORY
  !        Written,  Stephan Thober, Nov 2011
  !        Modified, Stephan Thober, Nov 2011 - added comments

  ! ------------------------------------------------------------------------------
  subroutine Get_NcVar_4d_sp(Filename, VarName, Data)
    !
    implicit none
    !
    character(256), intent(in)                :: Filename
    character(256), intent(in)                :: VarName ! Variable name
    integer(i4)                               :: ncid    ! id of input stream
    integer(i4)                               :: status  ! status of read stream
    integer(i4)                               :: varid   ! id of variable to be read
    integer(i4)                               :: vartype ! type of variable
    integer(i4), dimension(4)                 :: dl      ! length of dimensions
    real(sp), dimension(:,:,:,:), allocatable, &
         intent(inout) :: Data    ! array where values should be stored
    !
    if (allocated(data)) stop 'ERROR*** data is already allocated. subroutine Get_Ncvar'
    !
    ! Open NetCDF filename
    status = nf90_open(trim(Filename),NF90_NOWRITE, ncid)
    if ( status /= 0) stop 'ERROR*** nc file could not be opened!'
    !
    ! Inquire file and check if VarName exists in the dataset
    ! Get also the id and the length of the dimensions
    call get_Info(Varname,ncid,varid,dl,vartype)
    !
    ! check variable type ( 5 equals single precision )
    if ( vartype /= 5 ) stop 'ERROR*** type of variable does not match argument type. subroutine get_Info'
    !
    ! get values by varid
    allocate(data(dl(1),dl(2),dl(3),dl(4)))
    status = nf90_get_var(ncid, varid, data)
    if ( status /= 0) stop 'ERROR*** data could not be read!'
    !
    ! close File
    status = nf90_close(ncid)
    if ( status /= 0) stop 'ERROR*** nc file could not be closed!'
    !
  end subroutine Get_NcVar_4d_sp

  ! ------------------------------------------------------------------------------
  ! SUBROUTINE GET_NCVAR
  ! ------------------------------------------------------------------------------
  subroutine Get_NcVar_4d_dp(Filename, VarName, Data)
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
    if (allocated(data)) stop 'ERROR*** data is already allocated. subroutine Get_Ncvar'
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
  end subroutine Get_NcVar_4d_dp

  ! ------------------------------------------------------------------------------
  ! SUBROUTINE GET_NCVAR
  ! ------------------------------------------------------------------------------
  subroutine Get_NcVar_2d_dp(Filename, VarName, Data)
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
    if (allocated(data)) stop 'ERROR*** data is already allocated. subroutine Get_Ncvar'
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
  end subroutine Get_NcVar_2d_dp

  ! ------------------------------------------------------------------------------
  ! SUBROUTINE GET_NCVAR
  ! ------------------------------------------------------------------------------
  subroutine Get_NcVar_2d_sp(Filename, VarName, Data)
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
    if (allocated(data)) stop 'ERROR*** data is already allocated. subroutine Get_Ncvar'
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
  end subroutine Get_NcVar_2d_sp

  ! ------------------------------------------------------------------------------
  ! SUBROUTINE GET_NCVAR
  ! ------------------------------------------------------------------------------
  subroutine Get_NcVar_3d_dp(Filename, VarName, Data)
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
    if (allocated(data)) stop 'ERROR*** data is already allocated. subroutine Get_Ncvar'
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
  end subroutine Get_NcVar_3d_dp

  ! ------------------------------------------------------------------------------
  ! SUBROUTINE GET_NCVAR
  ! ------------------------------------------------------------------------------
  subroutine Get_NcVar_3d_sp(Filename, VarName, Data)
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
    if (allocated(data)) stop 'ERROR*** data is already allocated. subroutine Get_Ncvar'
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
  end subroutine Get_NcVar_3d_sp

  ! ------------------------------------------------------------------------------
  !
  ! SUBROUTINE GET_INFO
  !
  ! This subroutine is PRIVATE and therefore does not exist outside of this module.
  !
  ! This subroutine inquires the nc file. Given the Variable name and the stream 
  ! of the nc file, this subroutine determines the variable id, the kind of the
  ! variable and the length of the dimensions in the file.
  !
  ! See: http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90/ -> NF90-INQUIRE
  ! for detailed information.
  !
  ! ------------------------------------------------------------------------------
  subroutine get_Info(Varname, ncid, varid, dl, kind)
    !
    implicit none
    !
    character(256),            intent(in)     :: Varname
    integer(i4),               intent(in)     :: ncid 
    integer(i4), dimension(:), intent(inout)  :: dl
    integer(i4),               intent(out)    :: varid   ! variable id of data to be read
    integer(i4),               intent(out)    :: kind    ! type of the variable
    !
    integer(i4)                               :: nVars   ! Number of Variables
    integer(i4), dimension(:), allocatable    :: varids  ! Ids of all variables
    integer(i4), dimension(:), allocatable    :: DimID   ! Id of dimension
    integer(i4)                               :: NumDims ! Number of Dimensions for specific variable
    character(NF90_MAX_NAME)                  :: name    ! name of Variables in the file
    integer(i4)                               :: n       ! loop index
    integer(i4)                               :: status  ! control variable
    logical                                   :: flag    ! flag if Variable is in file
    !
    ! inquire number of variables
    status = nf90_inquire(ncid,nVariables = nVars)
    if ( status /= 0 ) stop 'ERROR*** Could NOT inquire Number of Variables'
    !
    allocate(varids(nVars))
    !
    ! get the ids of the variable
    status = nf90_inq_varids(ncid, nVars, varids)
    if ( status /= 0 ) stop 'ERROR*** Could NOT inquire IDs of Variables'
    !
    flag = .false.
    !
    varidloop: do n = 1, size(varids)
       !
       ! check whether Varname equals the name in the file
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
    ! check whether variable name was found in the file
    if ( .not. flag ) stop 'ERROR*** Varname does not exist in file. subroutine get_Info'
    !
    ! check consistency of dimensions
    if ( NumDims > size(dl) ) &
         stop 'ERROR*** Dimension size of Variable is greater than dims of array. subroutine get_Info'
    !
    ! if the data has less dimensions than the array, fill the dimension with ones
    if ( size(dl) > NumDims ) dl(NumDims+1:size(dl)) = 1
    !
    ! get the dimension ids ant the type of this variable
    allocate(DimId(NumDims))
    status = nf90_inquire_variable(ncid, varids(n), name, xtype=kind, dimids=DimId)
    if ( status /= 0 ) stop 'ERROR*** Could NOT inquire Ids of Dimensions'
    !
    ! go through dimension ids and get the dimesnion length
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
  !
end module mo_NcRead


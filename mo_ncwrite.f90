module mo_ncwrite
  !
  ! This module provides a structure and subroutines for writing netcdf files.

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! The UFZ Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The UFZ Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the UFZ Fortran library. If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2011-2012 Luis Samaniego, Stephan Thober, Matthias Cuntz

  use mo_kind,         only: i4, sp, dp
  use mo_string_utils, only: nonull

  ! functions and constants of netcdf4 library
  use netcdf,  only: nf90_create, nf90_def_dim, NF90_UNLIMITED, nf90_def_var, &
       NF90_CHAR, nf90_put_att, NF90_INT, NF90_INT, NF90_GLOBAL, &
       nf90_enddef, nf90_put_var, NF90_FLOAT, NF90_DOUBLE, &
       NF90_close, nf90_noerr, nf90_strerror, NF90_CLOBBER, &
       NF90_MAX_NAME, NF90_WRITE, nf90_inq_varid, nf90_inquire_variable, &
       nf90_inquire_dimension, nf90_open

  ! public routines -------------------------------------------------------------------
  public :: close_netcdf         ! save and close the netcdf file
  public :: create_netcdf        ! create the nc file with variables and their attributes, after they were set
  public :: dump_netcdf          ! simple dump of variable into a netcdf file
  public :: write_dynamic_netcdf ! write dynamically (one record after the other) in the file
  public :: write_static_netcdf  ! write static data in the file
  ! public types -------------------------------------------------------------------
  public :: dims
  public :: variable
  public :: attribute

  ! public parameters
  integer(i4), public, parameter :: nMaxDim = 5         ! nr. max dimensions
  integer(i4), public, parameter :: nMaxAtt = 20        ! nr. max attributes
  integer(i4), public, parameter :: maxLen  = 256       ! nr. string length
  integer(i4), public, parameter :: nGAtt   = 20        ! nr. global attributes
  integer(i4), public, parameter :: nAttDim = 2         ! dim array of attribute values

  ! public types -----------------------------------------------------------------
  type dims
     character (len=maxLen)                 :: name                ! dim. name
     integer(i4)                            :: len                 ! dim. lenght, undefined time => NF90_UNLIMITED
     integer(i4)                            :: dimId               ! dim. Id
  end type dims

  type attribute
     character (len=maxLen)                  :: name                ! attribute name
     integer(i4)                             :: xType               ! attribute of the values
     integer(i4)                             :: nValues             ! number of attributes
     character (len=maxLen)                  :: values              ! numbers or "characters" separed by spaces
  end type attribute

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
     type(attribute), dimension(nMaxAtt)    :: att                 ! var. attributes
     integer(i4), dimension(nMaxDim)        :: start               ! starting indices for netcdf
     integer(i4), dimension(nMaxDim)        :: count               ! counter          for netcdf
     logical                                :: wFlag               ! write flag
     integer(i4),                     pointer :: G0_i              ! array pointing model variables
     integer(i4), dimension(:      ), pointer :: G1_i              ! array pointing model variables
     integer(i4), dimension(:,:    ), pointer :: G2_i              ! array pointing model variables
     integer(i4), dimension(:,:,:  ), pointer :: G3_i              ! array pointing model variables
     integer(i4), dimension(:,:,:,:), pointer :: G4_i              ! array pointing model variables
     real(sp),                        pointer :: G0_f              ! array pointing model variables
     real(sp),    dimension(:    ),   pointer :: G1_f              ! array pointing model variables
     real(sp),    dimension(:,:    ), pointer :: G2_f              ! array pointing model variables
     real(sp),    dimension(:,:,:  ), pointer :: G3_f              ! array pointing model variables
     real(sp),    dimension(:,:,:,:), pointer :: G4_f              ! array pointing model variables
     real(dp),                        pointer :: G0_d              ! array pointing model variables
     real(dp),    dimension(:      ), pointer :: G1_d              ! array pointing model variables
     real(dp),    dimension(:,:    ), pointer :: G2_d              ! array pointing model variables
     real(dp),    dimension(:,:,:  ), pointer :: G3_d              ! array pointing model variables
     real(dp),    dimension(:,:,:,:), pointer :: G4_d              ! array pointing model variables
  end type variable

  ! public variables -----------------------------------------------------------------
  integer(i4),     public                            :: nVars   ! nr. variables
  integer(i4),     public                            :: nDims   ! nr. dimensions
  type (dims),     public, dimension(:), allocatable :: Dnc     ! dimensions list
  type(variable),  public, dimension(:), allocatable :: V       ! variable list, THIS STRUCTURE WILL BE WRITTEN IN THE FILE
  type(attribute), public, dimension(nGAtt)          :: gatt    ! global attributes for netcdf

  ! ------------------------------------------------------------------

  !     NAME
  !         dump_netcdf

  !     PURPOSE
  !         Simple write of a variable in a netcdf file.

  !         The variabel can be 1 to 5 dimensional and single or double precision.

  !         1D and 2D are dumped as static variables. From 3 to 5 dimension, the last
  !         dimension will be defined as time.
  !         The Variable will be called var.

  !     CALLING SEQUENCE
  !         call dump_netcdf(filename, arr)

  !     INTENT(IN)
  !         character(len=*) :: filename                  name of netcdf output file
  !         real(sp/dp) ::      arr(:[,:[,:[,:[,:]]]])    1D to 5D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         If dimension

  !     EXAMPLE
  !         call dump_netcdf('test.nc, myarray)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Nov 2012
  !         Modified, Stephan Thober, Nov 2012 - added functions for i4 variables
  !         Modified, Matthias Cuntz and Juliane Mai
  !                                   Nov 2012 - append
  !                                            - fake time dimension for 1D and 2D
  !                                            - make i4 behave exactly as sp and dp
  interface dump_netcdf
     module procedure dump_netcdf_1d_sp, dump_netcdf_2d_sp, dump_netcdf_3d_sp, &
          dump_netcdf_4d_sp, dump_netcdf_5d_sp, &
          dump_netcdf_1d_dp, dump_netcdf_2d_dp, dump_netcdf_3d_dp, &
          dump_netcdf_4d_dp, dump_netcdf_5d_dp, &
          dump_netcdf_1d_i4, dump_netcdf_2d_i4, dump_netcdf_3d_i4, &
          dump_netcdf_4d_i4, dump_netcdf_5d_i4
  end interface dump_netcdf

  ! ----------------------------------------------------------------------------

  private

  ! ----------------------------------------------------------------------------

contains

  ! ----------------------------------------------------------------------------
  !
  ! NAME
  !     close_netcdf
  !
  ! PURPOSE
  !     closes a stream of an open netcdf file and saves the file.
  !
  ! CALLING SEQUENCE
  !     call close_netcdf(nc, Info=Info)
  !
  ! INTENT(IN)
  !     integer(i4) :: nc - stream id of an open netcdf file which shall be closed
  !
  ! INTENT(IN), OPTIONAL
  !     logical :: Info - if given and true, progress information will be written
  !
  ! RESTRICTIONS
  !     closes only an already open stream
  !
  ! EXAMPLE
  !     see test_mo_ncwrite
  !
  ! LITERATURE
  !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html
  !
  ! HISTORY
  !     Written  Luis Samaniego  Feb 2011
  !     Modified Stephan Thober  Dec 2011 - added comments and generalized
  !              Matthias Cuntz  Jan 2012 - Info

  subroutine close_netcdf(ncId,Info)
    !
    implicit none
    !
    integer(i4), intent(in)           :: ncId
    logical,     intent(in), optional :: Info
    !
    ! close: save new netcdf dataset
    call check(nf90_close( ncId ))
    if (present(Info)) then
       if (Info) write(*,*) "NetCDF file was saved."
    end if
    !
  end subroutine close_netcdf

  ! ------------------------------------------------------------------------------
  !
  ! NAME
  !     create_netcdf
  !
  ! PURPOSE
  !     This subroutine will open a new netcdf file and write the variable
  !     attributes stored in the structure V in the file. Therefore V
  !     has to be already set. See the file set_netcdf in test_mo_ncwrite
  !     for an example.
  !
  ! CALLING SEQUENCE
  !     call create_netcdf(File, nc, Info=Info)
  !
  ! INTENT(IN)
  !     character(len=maxLen) :: File ! Filename of File to be written
  !
  ! INTENT(IN), OPTIONAL
  !     logical               :: Info ! if given and true, progress information will be written
  !
  ! INTENT(OUT)
  !     integer(i4)           :: nc   ! integer value of the stream for the opened file
  !
  ! RESTRICTIONS
  !     This routine only writes attributes and variables which have been stored in V
  !     nothing else.
  !
  ! EXAMPLE
  !     see test_mo_ncwrite
  !
  ! LITERATURE
  !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html
  !
  ! HISTORY
  !     Written  Luis Samaniego  Feb 2011
  !     Modified Stephan Thober  Dec 2011 - added comments and generalized
  !              Matthias Cuntz  Jan 2012 - Info

  subroutine create_netcdf(Filename,ncid,Info)
    !
    implicit none
    !
    ! netcdf related variables
    character(len=*), intent(in)           :: Filename
    integer(i4),      intent(out)          :: ncid
    logical,          intent(in), optional :: Info
    !
    integer(i4)                               :: i, j, k
    integer(i4), dimension(nAttDim)           :: att_INT
    real(sp),    dimension(nAttDim)           :: att_FLOAT
    character(len=maxLen), dimension(nAttDim) :: att_CHAR
    !
    ! 1  Create netcdf dataset: enter define mode     ->  get ncId
    call check(nf90_create(trim(Filename), NF90_CLOBBER, ncId ))
    !
    ! 2  Define dimensions                                 -> get dimId
    do i=1, nDims
       call check(nf90_def_dim(ncId, Dnc(i)%name , Dnc(i)%len, Dnc(i)%dimId ))
    end do
    !
    ! 3 Define dimids array, which is used to pass the dimids of the dimensions of
    ! the netcdf variables
    do i = 1, nVars
       V(i)%unlimited = .false.
       V(i)%dimids  = 0
       V(i)%start   = 1
       V(i)%count   = 1
       do k = 1, V(i)%nDims
          if ( Dnc( V(i)%dimTypes(k) )%len  == NF90_UNLIMITED ) V(i)%unlimited = .true.
          V(i)%dimids(k) = Dnc( V(i)%dimTypes(k) )%dimId
       end do
       if ( V(i)%unlimited ) then
          ! set counts for unlimited files (time is always the last dimension)
          if (V(i)%nDims == 1) cycle
          do k = 1, V(i)%nDims - 1
             V(i)%count(k)  = Dnc( V(i)%dimTypes(k) )%len
          end do
       end if
    end do
    !
    ! 4 Define the netcdf variables and atributes                            -> get varId
    do i=1, nVars
       if ( .not. V(i)%wFlag ) cycle
       call check(nf90_def_var(ncId, V(i)%name, V(i)%xtype, V(i)%dimids(1:V(i)%nDims), V(i)%varId ))
       do k = 1, V(i)%nAtt
          select case (V(i)%att(k)%xType)
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
       if (nonull(Gatt(k)%name)) then
          call check(nf90_put_att(ncId, NF90_GLOBAL, Gatt(k)%name, Gatt(k)%values))
       endif
    end do
    !
    ! 6 end definitions: leave define mode
    call check(nf90_enddef( ncId ))
    if (present(Info)) then
       if (Info) write(*,*) "NetCDF file was created", ncId
    end if
    !
  end subroutine create_netcdf

  ! ------------------------------------------------------------------

  subroutine dump_netcdf_1d_sp(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:),     intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 1 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_1d_sp: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim+1
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim+1) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_1d_sp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_1d_sp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_1d_sp: time name problem."
          endif
       enddo
       !
       ! append
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = dims(ndim+1) + i
          call check(nf90_put_var(ncid, varid(ndim+1), (/dims(ndim+1)+i/), (/dims(ndim+1)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))
       !
       ! define dim variables
       do i=1, ndim
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+2)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = i
          call check(nf90_put_var(ncid, varid(ndim+1), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_1d_sp


  subroutine dump_netcdf_2d_sp(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:,:),   intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 2 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_2d_sp: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim+1
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim+1) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_2d_sp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_2d_sp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_2d_sp: time name problem."
          endif
       enddo
       !
       ! append
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = dims(ndim+1) + i
          call check(nf90_put_var(ncid, varid(ndim+1), (/dims(ndim+1)+i/), (/dims(ndim+1)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))
       !
       ! define dim variables
       do i=1, ndim
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+2)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = i
          call check(nf90_put_var(ncid, varid(ndim+1), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_2d_sp


  subroutine dump_netcdf_3d_sp(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 3 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_3d_sp: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_3d_sp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_3d_sp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_3d_sp: time name problem."
          endif
       enddo
       !
       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,i), start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))
       !
       ! define dim variables
       do i=1, ndim-1
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,i), start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_3d_sp


  subroutine dump_netcdf_4d_sp(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 4 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_4d_sp: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_4d_sp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_4d_sp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_4d_sp: time name problem."
          endif
       enddo
       !
       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,i), start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))
       !
       ! define dim variables
       do i=1, ndim-1
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,i), start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_4d_sp


  subroutine dump_netcdf_5d_sp(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:,:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 5 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_5d_sp: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_5d_sp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_5d_sp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_5d_sp: time name problem."
          endif
       enddo
       !
       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,:,i), start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))
       !
       ! define dim variables
       do i=1, ndim-1
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,:,i), start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_5d_sp


  subroutine dump_netcdf_1d_dp(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:),     intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 1 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_1d_dp: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim+1
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim+1) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_1d_dp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_1d_dp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_1d_dp: time name problem."
          endif
       enddo
       !
       ! append
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = dims(ndim+1) + i
          call check(nf90_put_var(ncid, varid(ndim+1), (/dims(ndim+1)+i/), (/dims(ndim+1)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))
       !
       ! define dim variables
       do i=1, ndim
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+2)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = i
          call check(nf90_put_var(ncid, varid(ndim+1), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_1d_dp


  subroutine dump_netcdf_2d_dp(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:,:),   intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 2 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_2d_dp: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim+1
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim+1) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_2d_dp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_2d_dp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_2d_dp: time name problem."
          endif
       enddo
       !
       ! append
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = dims(ndim+1) + i
          call check(nf90_put_var(ncid, varid(ndim+1), (/dims(ndim+1)+i/), (/dims(ndim+1)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))
       !
       ! define dim variables
       do i=1, ndim
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+2)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = i
          call check(nf90_put_var(ncid, varid(ndim+1), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_2d_dp


  subroutine dump_netcdf_3d_dp(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 3 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_3d_dp: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_3d_dp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_3d_dp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_3d_dp: time name problem."
          endif
       enddo
       !
       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,i), start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))
       !
       ! define dim variables
       do i=1, ndim-1
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,i), start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_3d_dp


  subroutine dump_netcdf_4d_dp(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 4 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_4d_dp: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_4d_dp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_4d_dp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_4d_dp: time name problem."
          endif
       enddo
       !
       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,i), start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))
       !
       ! define dim variables
       do i=1, ndim-1
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,i), start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_4d_dp


  subroutine dump_netcdf_5d_dp(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:,:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 5 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_5d_dp: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_5d_dp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_5d_dp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_5d_dp: time name problem."
          endif
       enddo
       !
       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,:,i), start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))
       !
       ! define dim variables
       do i=1, ndim-1
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,:,i), start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_5d_dp


  subroutine dump_netcdf_1d_i4(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),         dimension(:),     intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 1 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_1d_i4: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim+1
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim+1) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_1d_i4: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_1d_i4: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_1d_i4: time name problem."
          endif
       enddo
       !
       ! append
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = dims(ndim+1) + i
          call check(nf90_put_var(ncid, varid(ndim+1), (/dims(ndim+1)+i/), (/dims(ndim+1)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))
       !
       ! define dim variables
       do i=1, ndim
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+2)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = i
          call check(nf90_put_var(ncid, varid(ndim+1), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_1d_i4


  subroutine dump_netcdf_2d_i4(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),         dimension(:,:),   intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 2 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_2d_i4: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim+1
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim+1) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_2d_i4: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_2d_i4: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_2d_i4: time name problem."
          endif
       enddo
       !
       ! append
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = dims(ndim+1) + i
          call check(nf90_put_var(ncid, varid(ndim+1), (/dims(ndim+1)+i/), (/dims(ndim+1)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))
       !
       ! define dim variables
       do i=1, ndim
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+2)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = i
          call check(nf90_put_var(ncid, varid(ndim+1), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_2d_i4


  subroutine dump_netcdf_3d_i4(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),      dimension(:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 3 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_3d_i4: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_3d_i4: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_3d_i4: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_3d_i4: time name problem."
          endif
       enddo
       !
       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,i), start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))
       !
       ! define dim variables
       do i=1, ndim-1
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,i), start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_3d_i4


  subroutine dump_netcdf_4d_i4(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),         dimension(:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 4 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_4d_i4: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_4d_i4: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_4d_i4: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_4d_i4: time name problem."
          endif
       enddo
       !
       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,i), start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))
       !
       ! define dim variables
       do i=1, ndim-1
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,i), start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_4d_i4


  subroutine dump_netcdf_5d_i4(filename, arr, append)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),      dimension(:,:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file

    integer(i4),      parameter         :: ndim = 5 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)
    !
    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))
       !
       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_5d_i4: number of variable dimensions /= number of file variable dimensions."
       !
       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_5d_i4: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_5d_i4: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_5d_i4: time name problem."
          endif
       enddo
       !
       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,:,i), start, counter))
       end do
    else
       ! open file
       call check(nf90_create(trim(filename), NF90_CLOBBER, ncid))
       !
       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))
       !
       ! define dim variables
       do i=1, ndim-1
          call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
       end do
       ! define time variable
       call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       !
       ! define variable
       call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1)))
       !
       ! end define mode
       call check(nf90_enddef(ncid))
       !
       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do
       !
       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,:,i), start, counter))
       end do
    endif
    !
    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_5d_i4

  ! ----------------------------------------------------------------------------
  !
  ! NAME
  !     write_dynamic_netcdf
  !
  ! PURPOSE
  !     This routine writes data, where one dimension has the unlimited attribute.
  !     Therefore, the number of the record which should be written has to be
  !     specified.
  !
  ! CALLING SEQUENCE
  !     call write_dynamic_netcdf(nc, rec, Info=Info)
  !
  ! INTENT(IN)
  !     integer(i4) :: nc - stream id of an open netcdf file where data should be written
  !                         can be obtained by an create_netcdf call
  !
  ! INTENT(IN)
  !     integer(i4) :: rec - record id of record which will be written in the file
  !
  ! INTENT(IN), OPTIONAL
  !     logical :: Info - if given and true, progress information will be written
  !
  ! RESTRICTIONS
  !     Writes only data, where the data pointers of the structure V are assigned
  !     and where one dimension has the unlimited attribute. Moreover only one
  !     record will be written.
  !     Writes only 1 to 4 dim arrays, integer, single or double precision.
  !
  ! EXAMPLE
  !     see test_mo_ncwrite
  !
  ! LITERATURE
  !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html
  !
  ! HISTORY
  !     Written  Luis Samaniego  Feb 2011
  !     Modified Stephan Thober  Dec 2011 - added comments and generalized
  !              Matthias Cuntz  Jan 2012 - Info
  !              Stephan Thober  Jan 2012 - iRec is not optional

  subroutine write_dynamic_netcdf(ncId, irec, Info)
    !
    implicit none
    !
    ! netcdf related variables
    integer(i4), intent(in)           :: ncId
    integer(i4), intent(in)           :: iRec
    logical,     intent(in), optional :: Info
    !
    integer(i4)                       :: i
    ! NOTES: 1) netcdf file must be on *** data mode ***
    !        2) start and end of the data chuck is controled by
    !           V(:)%start and  V(:)%count
    !
    ! set values for variables (one scalar or grid at a time)
    !
    do i = 1, nVars
       if ( .not. V(i)%unlimited     ) cycle
       if ( .not. V(i)%wFlag   ) cycle
       if (present(Info)) then
          if ((iRec ==1) .and. Info) write(*,*) "Var. ",i , trim(V(i)%name) ," is  dynamic" !, iRec
       end if
       V(i)%start ( V(i)%nDims ) = iRec
       select case (V(i)%xtype)
       case (NF90_INT)
          select case (V(i)%nDims-1)
          case (0)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G0_i, V(i)%start ))
          case (1)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G1_i, V(i)%start, V(i)%count ))
          case (2)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G2_i, V(i)%start, V(i)%count ))
          case (3)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G3_i, V(i)%start, V(i)%count ))
          case (4)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G4_i, V(i)%start, V(i)%count ))
          end select
       case (NF90_FLOAT)
          select case (V(i)%nDims-1)
          case (0)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G0_f, V(i)%start ))
          case (1)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G1_f, V(i)%start, V(i)%count ))
          case (2)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G2_f, V(i)%start, V(i)%count ))
          case (3)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G3_f, V(i)%start, V(i)%count ))
          case (4)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G4_f, V(i)%start, V(i)%count ))
          end select
       case (NF90_DOUBLE)
          select case (V(i)%nDims-1)
          case (0)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G0_d, V(i)%start ))
          case (1)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G1_d, V(i)%start, V(i)%count ))
          case (2)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G2_d, V(i)%start, V(i)%count ))
          case (3)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G3_d, V(i)%start, V(i)%count ))
          case (4)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G4_d, V(i)%start, V(i)%count ))
          end select
       end select
    end do

  end subroutine write_dynamic_netcdf

  ! ----------------------------------------------------------------------------
  !
  ! NAME
  !     write_static_netcdf
  !
  ! PURPOSE
  !     This routines writes static data in the netcdf file that is data
  !     where no dimension has the unlimited attribute.
  !
  ! CALLING SEQUENCE
  !     call write_static_netcdf(nc, Info=Info)
  !
  ! INTENT(IN)
  !     integer(i4) :: nc - stream id of an open netcdf file where data should be written
  !                         can be obtained by an create_netcdf call
  !
  ! INTENT(IN), OPTIONAL
  !     logical :: Info - if given and true, progress information will be written
  !
  ! RESTRICTIONS
  !     Writes only data, where the data pointers of the structure V are assigned.
  !     Writes only 1 to 4 dim arrays, integer, single or double precision.
  !
  ! EXAMPLE
  !     see test_mo_ncwrite
  !
  ! LITERATURE
  !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html
  !
  ! HISTORY
  !     Written  Luis Samaniego  Feb 2011
  !     Modified Stephan Thober  Dec 2011 - added comments and generalized
  !              Matthias Cuntz  Jan 2012 - Info

  subroutine write_static_netcdf(ncId, Info)
    !
    implicit none
    !
    ! netcdf related variables
    integer(i4), intent(in)           :: ncId
    logical,     intent(in), optional :: Info
    !
    integer(i4)                       :: i

    ! write all static variables
    do i = 1, nVars
       if ( V(i)%unlimited     ) cycle
       if ( .not. V(i)%wFlag   ) cycle
       if (present(Info)) then
          if (Info) write(*,*) "Var. ", trim(V(i)%name) ," is  static"
       end if
       select case (V(i)%xtype)
       case (NF90_INT)
          select case (V(i)%nDims)
          case (0)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G0_i ))
          case (1)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G1_i ))
          case (2)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G2_i ))
          case (3)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G3_i ))
          case (4)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G4_i ))
          end select
       case (NF90_FLOAT)
          select case (V(i)%nDims)
          case (0)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G0_f ))
          case (1)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G1_f ))
          case (2)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G2_f ))
          case (3)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G3_f ))
          case (4)
             call check( nf90_put_var( ncId,  V(i)%varId, V(i)%G4_f ))
          end select
       case (NF90_DOUBLE)
          select case (V(i)%nDims)
          case (0)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G0_d ))
          case (1)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G1_d ))
          case (2)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G2_d ))
          case (3)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G3_d ))
          case (4)
             call check(nf90_put_var( ncId,  V(i)%varId, V(i)%G4_d ))
          end select
       end select
    end do

  end subroutine write_static_netcdf

  ! -----------------------------------------------------------------------------
  !  private error checking routine
  subroutine check(status)
    !
    implicit none
    !
    integer(i4), intent(in) :: status
    !
    if (status /= nf90_noerr) then
       write(*,*) trim(nf90_strerror(status))
       stop
    end if
    !
  end subroutine check

end module mo_ncwrite

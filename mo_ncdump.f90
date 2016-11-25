module mo_ncdump

  ! This module provides subroutines for dumping variables into netcdf files.

  ! License
  ! -------
  ! This file is part of the JAMS Fortran library.

  ! The JAMS Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The JAMS Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the JAMS Fortran library (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2012-2016 Matthias Cuntz, Stephan Thober

  use mo_kind, only: i4, sp, dp

  ! functions and constants of netcdf4 library
  use netcdf,  only: nf90_create, nf90_def_dim, NF90_UNLIMITED, nf90_def_var, &
       NF90_INT, nf90_enddef, nf90_put_var, NF90_FLOAT, NF90_DOUBLE, &
       NF90_close, nf90_noerr, nf90_strerror, NF90_CLOBBER, &
       NF90_MAX_NAME, NF90_WRITE, nf90_inq_varid, nf90_inquire_variable, &
       nf90_inquire_dimension, nf90_open, NF90_64BIT_OFFSET, &
       nf90_inq_varid
#ifndef NETCDF3
  use netcdf, only: NF90_NETCDF4
#endif

  ! public routines -------------------------------------------------------------------
  public :: dump_netcdf          ! simple dump of variable into a netcdf file

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
  !         call dump_netcdf(filename, arr, append, lfs, netcdf4, deflate_level)

  !     INTENT(IN)
  !         character(len=*) :: filename                  name of netcdf output file
  !         real(sp/dp)      :: arr(:[,:[,:[,:[,:]]]])    1D to 5D-array with input numbers

  !     INTENT(IN), OPTIONAL
  !         logical          :: lfs                       True: enable netcdf3 large file support, i.e. 64-bit offset
  !         logical          :: netcdf4,                  True: use netcdf4 format
  !         integer(i4)      :: deflate_level             Compression level in netcdf4 (default: 1)

  !     RESTRICTIONS
  !         If dimension

  !     EXAMPLE
  !         call dump_netcdf('test.nc', myarray)
  !         call dump_netcdf('test.nc', myarray, netcdf4=.true.)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Nov 2012
  !         Modified, Stephan Thober, Nov 2012 - added functions for i4 variables
  !                   Matthias Cuntz and Juliane Mai
  !                                   Nov 2012 - append
  !                                            - fake time dimension for 1D and 2D
  !                                            - make i4 behave exactly as sp and dp
  !                                   Mar 2013 - lfs, netcdf4, deflate_level
  !                                   Nov 2016 - netcdf4 and deflate_level if compiled with -DNETCDF3
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

#ifdef NETCDF3
  INTEGER, PARAMETER :: NF90_NETCDF4 = NF90_64BIT_OFFSET ! to be available for compilation
#endif

  ! ----------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------

  subroutine dump_netcdf_1d_sp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:),     intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 1 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim+1) :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_1d_sp: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))

       ! define dim variables
       do i=1, ndim
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim) = dims(1:ndim)
          chunksizes(ndim+1) = 1
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+2) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+2)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_1d_sp


  subroutine dump_netcdf_2d_sp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:,:),   intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 2 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim+1) :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_2d_sp: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))

       ! define dim variables
       do i=1, ndim
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim) = dims(1:ndim)
          chunksizes(ndim+1) = 1
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+2) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+2)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_2d_sp


  subroutine dump_netcdf_3d_sp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 3 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_3d_sp: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_3d_sp


  subroutine dump_netcdf_4d_sp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 4 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_4d_sp: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_4d_sp


  subroutine dump_netcdf_5d_sp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:,:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 5 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_5d_sp: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_5d_sp


  subroutine dump_netcdf_1d_dp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:),     intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 1 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim+1) :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_1d_dp: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))

       ! define dim variables
       do i=1, ndim
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim) = dims(1:ndim)
          chunksizes(ndim+1) = 1
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+2) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+2)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_1d_dp


  subroutine dump_netcdf_2d_dp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:,:),   intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 2 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim+1) :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_2d_dp: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))

       ! define dim variables
       do i=1, ndim
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim) = dims(1:ndim)
          chunksizes(ndim+1) = 1
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+2) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+2)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_2d_dp


  subroutine dump_netcdf_3d_dp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 3 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_3d_dp: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
          ! call check(nf90_set_fill(ncid, NF90_NOFILL, old_fill_mode))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(Filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_3d_dp


  subroutine dump_netcdf_4d_dp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 4 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_4d_dp: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_4d_dp


  subroutine dump_netcdf_5d_dp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:,:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 5 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_5d_dp: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_5d_dp


  subroutine dump_netcdf_1d_i4(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),         dimension(:),     intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 1 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim+1) :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_1d_i4: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))

       ! define dim variables
       do i=1, ndim
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim) = dims(1:ndim)
          chunksizes(ndim+1) = 1
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+2) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+2)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_1d_i4


  subroutine dump_netcdf_2d_i4(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),         dimension(:,:),   intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 2 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim+1) :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_2d_i4: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))

       ! define dim variables
       do i=1, ndim
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim) = dims(1:ndim)
          chunksizes(ndim+1) = 1
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+2) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+2)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_2d_i4


  subroutine dump_netcdf_3d_i4(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),      dimension(:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 3 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_3d_i4: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_3d_i4


  subroutine dump_netcdf_4d_i4(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),         dimension(:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 4 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_4d_i4: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_4d_i4


  subroutine dump_netcdf_5d_i4(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),      dimension(:,:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 5 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

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
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
#ifndef NETCDF3
    if (present(netcdf4)) inetcdf4 = netcdf4
#endif
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_5d_i4: number of variable dimensions /= number of file variable dimensions."

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
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1) &
#ifndef NETCDF3
               , chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       else
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

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

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_5d_i4

  
  ! -----------------------------------------------------------------------------
  ! PRIVATE PART
  !

  ! -----------------------------------------------------------------------------

  !  private error checking routine
  subroutine check(status)

    implicit none

    integer(i4), intent(in) :: status

    if (status /= nf90_noerr) then
       write(*,*) 'mo_ncdump.check error: ', trim(nf90_strerror(status))
       stop
    end if

  end subroutine check

  ! -----------------------------------------------------------------------------

end module mo_ncdump

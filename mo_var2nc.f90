module mo_var2nc

  ! This module provides a routines for writing netcdf files.

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2014-2016 Stephan Thober, Matthias Cuntz - mc (at) macu (dot) de
  !
  ! Permission is hereby granted, free of charge, to any person obtaining a copy
  ! of this software and associated documentation files (the "Software"), to deal
  ! in the Software without restriction, including without limitation the rights
  ! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  ! copies of the Software, and to permit persons to whom the Software is
  ! furnished to do so, subject to the following conditions:
  !
  ! The above copyright notice and this permission notice shall be included in all
  ! copies or substantial portions of the Software.
  !
  ! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  ! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  ! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  ! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  ! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  ! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  ! SOFTWARE.

  use mo_kind,  only: i4, sp, dp
  use mo_utils, only: ne

  ! functions and constants of netcdf4 library
  use netcdf, only: &
       NF90_INT, NF90_FLOAT, NF90_DOUBLE, NF90_UNLIMITED, NF90_WRITE, NF90_NOERR, &
       nf90_create, nf90_def_dim, nf90_def_var, nf90_put_att, &
       nf90_enddef, nf90_put_var, nf90_close, nf90_strerror, &
       nf90_inq_varid, nf90_inquire_variable, nf90_inquire_dimension, nf90_open, &
       nf90_inq_varid, nf90_inq_dimid, nf90_inquire, nf90_get_var, nf90_fill_float, &
       nf90_fill_double, nf90_fill_int, nf90_redef
#ifndef __NETCDF3__
  use netcdf, only: NF90_NETCDF4
#else
  use netcdf, only: NF90_64BIT_OFFSET
#endif

  ! public routines -------------------------------------------------------------------
  public :: var2nc               ! simple dump of multiple variables with attributes into one netcdf file


  ! ------------------------------------------------------------------
  !
  !     NAME
  !         var2nc
  !
  !     PURPOSE
  !         Provide an intermediate way to write netcdf files
  !         between dump_netcdf and create_netcdf
  !
  !>        \brief Extended dump_netcdf for multiple variables
  !
  !>        \details Write different variables including attributes to netcdf
  !>        file. The attributes are restricted to long_name, units,
  !>        and missing_value. It is also possible to append variables
  !>        when an unlimited dimension is specified.
  !
  !     CALLING SEQUENCE
  !>        call var2nc(f_name, arr, dnames, v_name, dim_unlimited,
  !>                    long_name, units, missing_value, create)
  !
  !     INTENT(IN)
  !>        \param[in] "character(*) :: f_name"                             filename
  !>        \param[in] "integer(i4)/real(sp,dp) :: arr(:[,:[,:[,:[,:]]]])"  array to write
  !>        \param[in] "character(*) :: dnames(:)"                          dimension names
  !>        \param[in] "character(*) :: v_name"                             variable name
  !
  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4),             optional :: dim_unlimited"     index of unlimited dimension
  !>        \param[in] "character(*),            optional :: long_name"         attribute
  !>        \param[in] "character(*),            optional :: units"             attribute
  !>        \param[in] "integer(i4)/real(sp,dp), optional :: missing_value"     attribute
  !>        \param[in] "character(256), dimension(:,:), optional :: attributes" two dimensional array of attributes
  !>           size of first dimension equals number of attributes
  !>           first entry of second dimension equals attribute name (e.g. long_name)
  !>           second entry of second dimension equals attribute value (e.g. precipitation)
  !>           every attribute is written as string with the exception of missing_value
  !>        \param[in] "logical,                 optional :: create"            flag - specify whether a
  !>                                                                            output file should be
  !>                                                                            created, default
  !>        \param[inout] "integer(i4)/real(sp,dp), optional :: ncid"   if not given filename will be opened and closed
  !>                                                                    if given and <0 then file will be opened
  !>                                                                    and ncid will return the file unit.
  !>                                                                    if given and >0 then file is assumed open and
  !>                                                                    ncid is used as file unit.
  !>        \param[in] "integer(i4),                optional :: nrec"   if given: start point on unlimited dimension.
  !
  !     RESTRICTIONS
  !>        \note It is not allowed to write the folloing numbers for the indicated type\n
  !>                        number | kind\n
  !>                -2.1474836E+09 | integer(i4)\n
  !>                 9.9692100E+36 | real(sp)\n
  !>        9.9692099683868690E+36 | real(dp)\n
  !>        These numbers are netcdf fortran 90 constants! They are used to determine the
  !>        chunksize of the already written variable. Hence, this routine cannot append
  !>        correctly to variables when these numbers are used. Only five dimensional
  !>        variables can be written, only one unlimited dimension can be defined.
  !>
  !>        The unlimited dimension can be any dimension in netcdf4 files but must be the
  !>        last dimension when NETCDF3 is used (compiled with -DNETCDF3).
  !
  !     EXAMPLE
  !         Let <field> be some three dimensional array
  !           dnames(1) = 'x'
  !           dnames(2) = 'y'
  !           dnames(3) = 'time'
  !
  !         The simplest call to write <field> to a file is
  !           call var2nc('test.nc', field, dnames, 'h')
  !
  !         With attributes it looks like
  !           call var2nc('test.nc', field, dnames, 'h', &
  !                       long_name = 'height', units = '[m]', missing_value = -9999)
  !         or alternatively
  !         character(256), dimension(3,2) :: attributes
  !         attributes(1,1) = 'long_name'
  !         attributes(1,2) = 'precipitation'
  !         attributes(2,1) = 'units'
  !         attributes(2,2) = '[mm/d]'
  !         attributes(3,1) = 'missing_value'
  !         attributes(3,2) = '-9999.'
  !         call var2nc('test.nc', field, dnames, 'h', attributes = attributes, create = .true. )

  !         To be able to dynamically write <field>, an unlimited dimension
  !         needs to be specified (in this example field could also be only two
  !         dimensional)
  !         call var2nc('test.nc', field(:,:,1), dnames, 'h', dim_unlimited=3)
  !         Now one can append an arbitrary number of time steps, e.g., the next 9
  !         and the time has to be added again before
  !           call var2nc('test.nc', (/20,...,100/), dnames(3:3), 'time',
  !                       dim_unlimited = 1 )
  !           call var2nc('test.nc', field(:,:,2:10, dnames, 'h', dim_unlimited=3)
  !
  !         You can also write another variable sharing the same dimensions
  !           call var2nc('test.nc', field_2, dnames(1:2), 'h_2')
  !
  !         The netcdf file can stay open after the first call and subsequent calls can use the file unit
  !           ncid = -1_i4
  !           call var2nc('test.nc', field_1, dnames(1:1), 'h_1', ncid=ncid) ! opens file
  !           call var2nc('test.nc', field_2, dnames(1:2), 'h_2', ncid=ncid) ! uses ncid from last call
  !           call close_netcdf(ncid)
  !
  !         One can also give the start record number (on the unlimited dimension)
  !           ncid = -1_i4
  !           call var2nc('test.nc', time1, dnames(3:3), 'time', dim_unlimited=1_i4, ncid=ncid, create=.true.)
  !           do i=1, n
  !               call var2nc('test.nc', field_2(:,:,i), dnames(1:3), 'h_2', dim_unlimited=3_i4, ncid=ncid, nrec=i)
  !           end do
  !           call close_netcdf(ncid)
  !
  !         That's it, enjoy!
  !           -> see also example in test program
  !
  !    LITERATURE
  !        The manual of the used netcdf fortran library can be found in
  !        Robert Pincus & Ross Rew, The netcdf Fortran 90 Interface Guide
  !
  !    HISTORY
  !>       \author Stephan Thober & Matthias Cuntz
  !>       \date May 2014
  !        Modified Stephan Thober - Jun 2014 added deflate, shuffle, and chunksizes
  !                 Stephan Thober - Jun 2014 automatically append variable at the end,
  !                                           renamed _FillValue to missing_value
  !                 Stephan Thober - Jul 2014 add attributes array, introduced unlimited dimension
  !                                           that is added to the dimensions of the given array
  !                 Stephan Thober - Jan 2015 changed chunk_size convention to one chunk per unit in
  !                                           unlimited dimension (typically time)
  !                 Matthias Cuntz - Feb 2015 d_unlimit was not set in 5d cases
  !                                           use ne from mo_utils for fill value comparisons
  !                                           dummy(1) was sp instead of i4 in var2nc_1d_i4
  !                 Matthias Cuntz - May 2015 ncid for opening the file only once
  !                                           nrec for writing a specific record
  !                 Matthias Cuntz - Nov 2016 reopen definition mode when not create and file exists
  !                                           netCDF3: unlimited dimension must be last dimension when compiled with -DNETCDF3

  interface var2nc
     module procedure var2nc_1d_i4, var2nc_1d_sp, var2nc_1d_dp, &
          var2nc_2d_i4, var2nc_2d_sp, var2nc_2d_dp, &
          var2nc_3d_i4, var2nc_3d_sp, var2nc_3d_dp, &
          var2nc_4d_i4, var2nc_4d_sp, var2nc_4d_dp, &
          var2nc_5d_i4, var2nc_5d_sp, var2nc_5d_dp
  end interface var2nc

  ! ----------------------------------------------------------------------------

  private

#ifdef __NETCDF3__
  INTEGER, PARAMETER :: NF90_NETCDF4 = NF90_64BIT_OFFSET ! to be available for compilation
#endif

  ! ----------------------------------------------------------------------------

contains


  ! ------------------------------------------------------------------


  subroutine var2nc_1d_i4( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 1
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    integer(i4),      dimension(:),             intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! attributes
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    integer(i4),                      optional, intent(in) :: missing_value
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4)                            :: f_handle
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    integer(i4), dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size(dnames, 1)
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( size( dnames, 1) .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( size( dnames, 1) .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1) /)
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1
    counter(:)   = dims
    dummy        = nf90_fill_int
    dummy_count  = 1
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim .ne. ndim) stop "var2nc_1d_i4: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_1d_i4: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_1d_i4: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( dummy(1) /= nf90_fill_int ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimension
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_INT, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
       ))
#endif
       !
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(I6)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! inquire dimensions
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_1d_i4: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_1d_i4: variable dimension /= file variable dimension."
    enddo
    ! write time and variable
    call check(nf90_put_var(f_handle, varid(ndim+1), arr, start, counter ) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_1d_i4

  subroutine var2nc_1d_sp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 1
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(sp),         dimension(:),             intent(in) :: arr
    character(len=*), dimension(:),             intent(in) :: dnames
    character(len=*),                           intent(in) :: v_name
    ! attributes
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    real(sp),                         optional, intent(in) :: missing_value
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4)                            :: f_handle
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(sp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size(dnames, 1)
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( size( dnames, 1) .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( size( dnames, 1) .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1) /)
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1
    counter(:)   = dims
    dummy_count  = 1
    dummy        = nf90_fill_float
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim .ne. ndim) stop "var2nc_1d_sp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_1d_sp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_1d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_float) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimension
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_FLOAT, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! inquire dimensions
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_1d_sp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_1d_sp: variable dimension /= file variable dimension."
    enddo
    ! write time and variable
    call check(nf90_put_var(f_handle, varid(ndim+1), arr, start, counter ) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_1d_sp

  subroutine var2nc_1d_dp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 1
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(dp),         dimension(:),             intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! attributes
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    real(dp),                         optional, intent(in) :: missing_value
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4)                            :: f_handle
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(dp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size(dnames, 1)
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( size( dnames, 1) .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( size( dnames, 1) .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1) /)
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1
    counter(:)   = dims
    dummy_count  = 1
    dummy        = nf90_fill_double
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim .ne. ndim) stop "var2nc_1d_dp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_1d_dp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_1d_dp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_double) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimension
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_DOUBLE, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! inquire dimensions
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_1d_dp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_1d_dp: variable dimension /= file variable dimension."
    enddo
    ! write time and variable
    call check(nf90_put_var(f_handle, varid(ndim+1), arr, start, counter ) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_1d_dp

  subroutine var2nc_2d_i4( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 2
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    integer(i4),      dimension(:,:),           intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    integer(i4),                      optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    integer(i4), dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
         ( d_unlimit .lt. 0_i4 ) ) then
       write(*,*) '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       write(*,*) '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_int
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_2d_i4: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_2d_i4: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_2d_i4: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( dummy(1) /= nf90_fill_int ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_INT, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(I6)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_2d_i4: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_2d_i4: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_2d_i4

  subroutine var2nc_2d_sp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 2
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(sp),         dimension(:,:),           intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(sp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(sp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       write(*,*) '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       write(*,*) '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_float
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_2d_sp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_2d_sp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_2d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_float) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_FLOAT, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_2d_sp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_2d_sp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_2d_sp

  subroutine var2nc_2d_dp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 2
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(dp),         dimension(:,:),           intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(dp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(dp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       write(*,*) '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       write(*,*) '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1
    dummy        = nf90_fill_double
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_2d_dp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_2d_dp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_2d_dp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_double) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_DOUBLE, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value',missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_2d_dp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_2d_dp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_2d_dp

  subroutine var2nc_3d_i4( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 3
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    integer(i4),      dimension(:,:,:),         intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    integer(i4),                      optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    integer(i4), dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       write(*,*) '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       write(*,*) '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), size( arr, 3), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2), size( arr, 3) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_int
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_3d_i4: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_3d_i4: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_3d_i4: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( dummy(1) /= nf90_fill_int ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_INT, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(I6)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_3d_i4: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_3d_i4: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_3d_i4

  subroutine var2nc_3d_sp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 3
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(sp),         dimension(:,:,:),         intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(sp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(sp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       write(*,*) '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       write(*,*) '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    ! set chunk sizes and dimension names
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), size( arr, 3), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2), size( arr, 3) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_float
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_3d_sp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_3d_sp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_3d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_float) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_FLOAT, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
            ))
#endif
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_3d_sp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_3d_sp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_3d_sp

  subroutine var2nc_3d_dp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 3
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(dp),         dimension(:,:,:),         intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(dp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(dp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       write(*,*) '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       write(*,*) '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), size( arr, 3), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2), size( arr, 3) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1
    dummy        = nf90_fill_double
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_3d_dp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_3d_dp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_3d_dp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_double) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_DOUBLE, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value',missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_3d_dp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_3d_dp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_3d_dp

  subroutine var2nc_4d_i4( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 4
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    integer(i4),      dimension(:,:,:,:),       intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    integer(i4),                      optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    integer(i4), dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       write(*,*) '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       write(*,*) '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), &
            size( arr, 3), size( arr, 4 ),  1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2), &
            size( arr, 3), size( arr, 4) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_int
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_4d_i4: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_4d_i4: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_4d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( dummy(1) /= nf90_fill_int ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_INT, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(I6)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_4d_i4: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_4d_i4: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_4d_i4

  subroutine var2nc_4d_sp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 4
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(sp),         dimension(:,:,:,:),       intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(sp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(sp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       write(*,*) '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       write(*,*) '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), &
            size( arr, 3), size( arr, 4 ),  1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2), &
            size( arr, 3), size( arr, 4) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_float
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_4d_sp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_4d_sp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_4d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1), nf90_fill_float) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_FLOAT, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_4d_sp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_4d_sp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_4d_sp

  subroutine var2nc_4d_dp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 4
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(dp),         dimension(:,:,:,:),       intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(dp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(dp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       write(*,*) '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       write(*,*) '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! Initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), &
            size( arr, 3), size( arr, 4 ),  1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2), &
            size( arr, 3), size( arr, 4) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1
    dummy        = nf90_fill_double
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_4d_dp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_4d_dp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_4d_dp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_double) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_DOUBLE, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value',missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_4d_dp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_4d_dp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_4d_dp

  subroutine var2nc_5d_i4( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 5
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    integer(i4),      dimension(:,:,:,:,:),     intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    integer(i4),                      optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    integer(i4), dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    d_unlimit = 0_i4
    if (present(dim_unlimited)) d_unlimit = dim_unlimited
    ! consistency checks
    if ( ndim .gt. ndim_const ) then
       print *, '***ERROR more than five dimension names given'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       write(*,*) '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       write(*,*) '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    chunksizes   = (/ size( arr, 1), size( arr, 2), size( arr, 3), size( arr, 4), size( arr, 5) /)
    if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1
    dims(1:ndim) = shape( arr )
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_int
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_5d_i4: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_5d_i4: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_5d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( dummy(1) /= nf90_fill_int ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_INT, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(I6)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_5d_i4: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_5d_i4: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_5d_i4

  subroutine var2nc_5d_sp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 5
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(sp),         dimension(:,:,:,:,:),     intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(sp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(sp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    d_unlimit = 0_i4
    if (present(dim_unlimited)) d_unlimit = dim_unlimited
    ! consistency checks
    if ( ndim .gt. ndim_const ) then
       print *, '***ERROR more than five dimension names given'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       write(*,*) '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       write(*,*) '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    chunksizes   = (/ size( arr, 1), size( arr, 2), size( arr, 3), size( arr, 4), size( arr, 5) /)
    if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1
    dims(1:ndim) = shape( arr )
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_float
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_5d_sp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_5d_sp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_5d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_float) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_FLOAT, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_5d_sp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_5d_sp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_5d_sp

  subroutine var2nc_5d_dp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 5
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(dp),         dimension(:,:,:,:,:),     intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(dp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(dp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    d_unlimit = 0_i4
    if (present(dim_unlimited)) d_unlimit = dim_unlimited
    ! consistency checks
    if ( ndim .gt. ndim_const ) then
       print *, '***ERROR more than five dimension names given'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       write(*,*) '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       write(*,*) '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    chunksizes   = (/ size( arr, 1), size( arr, 2), size( arr, 3), size( arr, 4), size( arr, 5) /)
    if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1
    dims(1:ndim) = shape( arr )
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1
    dummy        = nf90_fill_double
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    else
       create_loc = .false.
    endif
    ! check whether variable exists
    if ( NF90_NOERR .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_5d_dp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_5d_dp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_5d_dp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_double) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! reopen definition
       if (.not. create_loc) call check(nf90_redef(f_handle))
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( NF90_NOERR .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_DOUBLE, dimid, varid(ndim+1) &
#ifndef __NETCDF3__
            , chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
#else
          ))
#endif
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value',missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_5d_dp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_5d_dp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_5d_dp


  ! -----------------------------------------------------------------------------
  ! PRIVATE PART
  !


  ! ----------------------------------------------------------------------------

  ! NAME
  !     close_netcdf

  ! PURPOSE
  !     closes a stream of an open netcdf file and saves the file.

  ! CALLING SEQUENCE
  !     call close_netcdf(nc)

  ! INTENT(IN)
  !     integer(i4) :: ncid - stream id of an open netcdf file which shall be closed

  ! RESTRICTIONS
  !     Closes only an already open stream

  ! EXAMPLE
  !     See test_mo_var2nc

  ! LITERATURE
  !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html

  ! HISTORY
  !     Written  Luis Samaniego  Feb 2011
  !     Modified Stephan Thober  Dec 2011 - added comments and generalized
  !              Matthias Cuntz  Jan 2012 - Info
  !              Matthias Cuntz  Mar 2013 - removed Info

  subroutine close_netcdf(ncId)

    implicit none

    integer(i4), intent(in)           :: ncId

    ! close: save new netcdf dataset
    call check(nf90_close(ncId))

  end subroutine close_netcdf

  ! -----------------------------------------------------------------------------

  ! private open netcdf function - returns file handle
  function open_netcdf( f_name, create )
    implicit none
    ! input variables
    character(len=*),           intent(in) :: f_name
    logical,                    intent(in) :: create ! flag indicates that file exists
    ! output variables
    integer(i4)   :: open_netcdf
    !
    if ( create ) then
       ! create file
#ifndef __NETCDF3__
       call check( nf90_create( trim(f_name), NF90_NETCDF4, open_netcdf ) )
#else
       call check( nf90_create( trim(f_name), NF90_64BIT_OFFSET, open_netcdf ) )
#endif
    else
       ! open file
       call check( nf90_open( trim(f_name), NF90_WRITE, open_netcdf ) )
    end if
  end function open_netcdf


  ! -----------------------------------------------------------------------------

  !  private error checking routine
  subroutine check(status)

    implicit none

    integer(i4), intent(in) :: status

    if (status /= NF90_NOERR) then
       write(*,*) 'mo_var2nc.check error: ', trim(nf90_strerror(status))
       stop
    end if

  end subroutine check

  ! -----------------------------------------------------------------------------

end module mo_var2nc

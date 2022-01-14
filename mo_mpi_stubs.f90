!> \file mo_mpi_stubs.f90

!> \brief Module provides dummy functions for common MPI routines.

!> \details This module provides dummy functions for the most commonly called MPI routines.\n
!>          Most of the stub routines do not do anything. In a few cases, where it makes sense,
!>          they do some simple action or return a value that is appropriate for the serial processing case.\n
!>          mo_mpi_stubs can be used as a convenience, when a real MPI implementation is not available,
!>          and the user simply wants to test-compile a code. It may also be useful in those occasions\n
!>          when a code has been so carefully written that it will still execute correctly on a single processor.\n
!>
!>          The original code was provided by J. Burkardt who based the code on a similar package
!>          supplied as part of the LAMMPS program.\n
!>          Burkardt's code can be downloaded at http://people.sc.fsu.edu/~jburkardt/f_src/mpi_stubs/mpi_stubs.html.
!>
!> \authors LAMMPS project, John Burkardt (F90), Matthias Cuntz (JAMS)
!> \date March 2015
module mo_mpi_stubs

  ! This module part of the JAMS Fortran library.

  ! Written Matthias Cuntz, March 2015

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2015 LAMMPS, John Burkardt, Matthias Cuntz - mc (at) macu (dot) de
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

  use mo_kind, only: i4, i8, sp, dp

  implicit none

  integer(i4), parameter :: mpi_comm_world = 0
  !
  !  Return values
  integer(i4), parameter :: mpi_failure = 1
  integer(i4), parameter :: mpi_success = 0
  !
  !  recv message status
  integer(i4), parameter :: mpi_status_size = 3
  integer(i4), parameter :: mpi_source = 1
  integer(i4), parameter :: mpi_tag = 2
  integer(i4), parameter :: mpi_count = 3
  !
  !  recv flags
  integer(i4), parameter :: mpi_any_source = -1
  integer(i4), parameter :: mpi_any_tag = -1
  !
  !  data types and sizes
  integer(i4), parameter :: mpi_integer = 1
  integer(i4), parameter :: mpi_real = 2
  integer(i4), parameter :: mpi_double_precision = 3
  integer(i4), parameter :: mpi_logical = 4
  integer(i4), parameter :: mpi_character = 5
  !
  !  allreduce operations
  integer(i4), parameter :: mpi_sum = 1
  integer(i4), parameter :: mpi_max = 2
  integer(i4), parameter :: mpi_min = 3
  integer(i4), parameter :: mpi_product = 4

  interface mpi_allgather
     module procedure mpi_allgather_i4, mpi_allgather_sp, mpi_allgather_dp
  end interface mpi_allgather

  interface mpi_allgatherv
     module procedure mpi_allgatherv_i4, mpi_allgatherv_sp, mpi_allgatherv_dp
  end interface mpi_allgatherv

  interface mpi_allreduce
     module procedure mpi_allreduce_i4, mpi_allreduce_sp, mpi_allreduce_dp
  end interface mpi_allreduce

  interface mpi_reduce
     module procedure mpi_reduce_0d_i4, mpi_reduce_1d_i4, mpi_reduce_sp, mpi_reduce_dp
  end interface mpi_reduce

  interface mpi_reduce_scatter
     module procedure mpi_reduce_scatter_i4, mpi_reduce_scatter_sp, mpi_reduce_scatter_dp
  end interface mpi_reduce_scatter

contains

  subroutine mpi_abort( comm, errorcode, ierror )

    !*****************************************************************************80
    !
    !! MPI_ABORT shuts down the processes in a given communicator.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    08 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Input, integer ERRORCODE, the error code to be returned.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: comm
    integer(i4) :: errorcode
    integer(i4) :: ierror

    ierror = MPI_SUCCESS

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_ABORT:'
    write ( *, '(a,i12)' ) '  Shut down with error code = ', errorcode

    stop

  end subroutine mpi_abort


  subroutine mpi_allgather_i4( data1, nsend, sendtype,data2, nrecv, recvtype, &
       comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_ALLGATHER gathers data from all the processes in a communicator.
    !
    !  Discussion:
    !
    !    Copy values from DATA1 to DATA2.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: nsend

    integer(i4) :: comm
    integer(i4) :: data1(nsend)
    integer(i4) :: data2(nsend)
    integer(i4) :: ierror
    integer(i4) :: nrecv
    integer(i4) :: recvtype
    integer(i4) :: sendtype

    ierror = MPI_SUCCESS

    if ( sendtype == mpi_integer ) then
       call mpi_copy_integer ( data1, data2, nsend, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_allgather_i4


  subroutine mpi_allgather_sp( data1, nsend, sendtype,data2, nrecv, recvtype, &
       comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_ALLGATHER gathers data from all the processes in a communicator.
    !
    !  Discussion:
    !
    !    Copy values from DATA1 to DATA2.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: nsend

    integer(i4) :: comm
    real(sp) :: data1(nsend)
    real(sp) :: data2(nsend)
    integer(i4) :: ierror
    integer(i4) :: nrecv
    integer(i4) :: recvtype
    integer(i4) :: sendtype

    ierror = MPI_SUCCESS

    if ( sendtype == mpi_real ) then
       call mpi_copy_real ( data1, data2, nsend, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_allgather_sp


  subroutine mpi_allgather_dp( data1, nsend, sendtype,data2, nrecv, recvtype, &
       comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_ALLGATHER gathers data from all the processes in a communicator.
    !
    !  Discussion:
    !
    !    Copy values from DATA1 to DATA2.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: nsend

    integer(i4) :: comm
    real(dp) :: data1(nsend)
    real(dp) :: data2(nsend)
    integer(i4) :: ierror
    integer(i4) :: nrecv
    integer(i4) :: recvtype
    integer(i4) :: sendtype

    ierror = MPI_SUCCESS

    if ( sendtype == mpi_double_precision ) then
       call mpi_copy_double_precision ( data1, data2, nsend, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_allgather_dp


  subroutine mpi_allgatherv_i4( data1, nsend, sendtype, data2, nrecv, ndispls, &
       recvtype, comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_ALLGATHERV gathers data from all the processes in a communicator.
    !
    !  Discussion:
    !
    !    Copy values from DATA1 to DATA2.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: nsend

    integer(i4) :: comm
    integer(i4) :: data1(nsend)
    integer(i4) :: data2(nsend)
    integer(i4) :: ierror
    integer(i4) :: ndispls
    integer(i4) :: nrecv
    integer(i4) :: recvtype
    integer(i4) :: sendtype

    ierror = MPI_SUCCESS

    if ( sendtype == mpi_integer ) then
       call mpi_copy_integer ( data1, data2, nsend, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_allgatherv_i4


  subroutine mpi_allgatherv_sp( data1, nsend, sendtype, data2, nrecv, ndispls, &
       recvtype, comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_ALLGATHERV gathers data from all the processes in a communicator.
    !
    !  Discussion:
    !
    !    Copy values from DATA1 to DATA2.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: nsend

    integer(i4) :: comm
    real(sp) :: data1(nsend)
    real(sp) :: data2(nsend)
    integer(i4) :: ierror
    integer(i4) :: ndispls
    integer(i4) :: nrecv
    integer(i4) :: recvtype
    integer(i4) :: sendtype

    ierror = MPI_SUCCESS

    if ( sendtype == mpi_real ) then
       call mpi_copy_real ( data1, data2, nsend, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_allgatherv_sp


  subroutine mpi_allgatherv_dp( data1, nsend, sendtype, data2, nrecv, ndispls, &
       recvtype, comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_ALLGATHERV gathers data from all the processes in a communicator.
    !
    !  Discussion:
    !
    !    Copy values from DATA1 to DATA2.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: nsend

    integer(i4) :: comm
    real(dp) :: data1(nsend)
    real(dp) :: data2(nsend)
    integer(i4) :: ierror
    integer(i4) :: ndispls
    integer(i4) :: nrecv
    integer(i4) :: recvtype
    integer(i4) :: sendtype

    ierror = MPI_SUCCESS

    if ( sendtype == mpi_double_precision ) then
       call mpi_copy_double_precision ( data1, data2, nsend, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_allgatherv_dp


  subroutine mpi_allreduce_i4( data1, data2, n, datatype, operation, comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_ALLREDUCE carries out a reduction operation.
    !
    !  Discussion:
    !
    !    The reduction operations are MAXIMUM, MINIMUM, PRODUCT and SUM.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !    Thanks to Simppa Akaslompolo for correcting this routine!
    !    12 January 2012.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    12 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, DATATYPE DATA1(N), the data to be processed.
    !
    !    Output, DATATYPE DATA2(N), the value of the reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    integer(i4) :: data1(n)
    integer(i4) :: data2(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: operation

    ierror = MPI_SUCCESS

    if ( datatype == mpi_integer ) then
       call mpi_reduce_integer ( data1, data2, n, operation, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_allreduce_i4


  subroutine mpi_allreduce_sp( data1, data2, n, datatype, operation, comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_ALLREDUCE carries out a reduction operation.
    !
    !  Discussion:
    !
    !    The reduction operations are MAXIMUM, MINIMUM, PRODUCT and SUM.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !    Thanks to Simppa Akaslompolo for correcting this routine!
    !    12 January 2012.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    12 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, DATATYPE DATA1(N), the data to be processed.
    !
    !    Output, DATATYPE DATA2(N), the value of the reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    real(sp) :: data1(n)
    real(sp) :: data2(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: operation

    ierror = MPI_SUCCESS

    if ( datatype == mpi_real ) then
       call mpi_reduce_real ( data1, data2, n, operation, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_allreduce_sp


  subroutine mpi_allreduce_dp( data1, data2, n, datatype, operation, comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_ALLREDUCE carries out a reduction operation.
    !
    !  Discussion:
    !
    !    The reduction operations are MAXIMUM, MINIMUM, PRODUCT and SUM.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !    Thanks to Simppa Akaslompolo for correcting this routine!
    !    12 January 2012.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    12 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, DATATYPE DATA1(N), the data to be processed.
    !
    !    Output, DATATYPE DATA2(N), the value of the reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    real(dp) :: data1(n)
    real(dp) :: data2(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: operation

    ierror = MPI_SUCCESS

    if ( datatype == mpi_double_precision ) then
       call mpi_reduce_double_precision ( data1, data2, n, operation, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_allreduce_dp


  subroutine mpi_barrier( comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_BARRIER forces processes within a communicator to wait together.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: comm
    integer(i4) :: ierror

    ierror = MPI_SUCCESS

    return

  end subroutine mpi_barrier


  subroutine mpi_bcast( data, n, datatype, node, comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_BCAST broadcasts data from one process to all others.
    !
    !  Discussion:
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, datatype DATA(N), the data to be broadcast.
    !
    !    Input, integer N, the number of items of data.
    !
    !    Input, integer DATATYPE, the MPI code for the datatype of the data.
    !
    !    Input, integer NODE, the rank of the sending process within the
    !    given communicator.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    integer(i4) :: data(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: node

    ierror = MPI_SUCCESS

    return

  end subroutine mpi_bcast


  subroutine mpi_bsend( data, n, datatype, iproc, itag, comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_BSEND sends data from one process to another, using buffering.
    !
    !  Discussion:
    !
    !    Warn against sending message to self, since no data copy is done.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, datatype DATA(N), the data to be sent.
    !
    !    Input, integer N, the number of data items to send.
    !
    !    Input, integer DATAYTPE, the MPI code for the datatype.
    !
    !    Input, integer IPROC, the rank of the process within the communicator
    !    that is to receive the message.
    !
    !    Input, integer ITAG, a tag for the message.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    integer(i4) :: data(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: iproc
    integer(i4) :: itag

    ierror = MPI_FAILURE

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_BSEND - Error!'
    write ( *, '(a)' )  '  Should not send message to self.'

    return

  end subroutine mpi_bsend


  subroutine mpi_cart_create( comm, ndims, dims, periods, reorder, comm_cart, &
       ierror )

    !*****************************************************************************80
    !
    !! MPI_CART_CREATE creates a communicator for a Cartesian topology.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: ndims

    integer(i4) :: comm
    integer(i4) :: comm_cart
    integer(i4) :: dims(*)
    integer(i4) :: ierror
    logical :: periods(*)
    logical :: reorder

    ierror = MPI_SUCCESS

    return

  end subroutine mpi_cart_create


  subroutine mpi_cart_get( comm, ndims, dims, periods, coords, ierror )

    !*****************************************************************************80
    !
    !! MPI_CART_GET returns the "Cartesian coordinates" of the calling process.
    !
    !  Discussion:
    !
    !    Set all coordinates to 0.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: ndims

    integer(i4) :: comm
    integer(i4) :: coords(*)
    integer(i4) :: dims(*)
    integer(i4) :: ierror
    logical :: periods(*)

    ierror = MPI_SUCCESS

    coords(1:ndims) = 0

    return

  end subroutine mpi_cart_get


  subroutine mpi_cart_shift( comm, idir, idisp, isource, idest, ierror )

    !*****************************************************************************80
    !
    !! MPI_CART_SHIFT finds the destination and source for Cartesian shifts.
    !
    !  Discussion:
    !
    !    Set ISOURCE = IDEST = SELF = 0.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: comm
    integer(i4) :: idest
    integer(i4) :: idir
    integer(i4) :: idisp
    integer(i4) :: ierror
    integer(i4) :: isource

    ierror = MPI_SUCCESS
    isource = 0
    idest = 0

    return

  end subroutine mpi_cart_shift


  subroutine mpi_comm_dup( comm, comm_out, ierror )

    !*****************************************************************************80
    !
    !! MPI_COMM_DUP duplicates a communicator.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: comm
    integer(i4) :: comm_out
    integer(i4) :: ierror

    ierror = MPI_SUCCESS

    return

  end subroutine mpi_comm_dup


  subroutine mpi_comm_free( comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_COMM_FREE "frees" a communicator.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: comm
    integer(i4) :: ierror

    ierror = MPI_SUCCESS

    return

  end subroutine mpi_comm_free


  subroutine mpi_comm_rank( comm, me, ierror )

    !*****************************************************************************80
    !
    !! MPI_COMM_RANK reports the rank of the calling process.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: comm
    integer(i4) :: ierror
    integer(i4) :: me

    ierror = MPI_SUCCESS
    me = 0

    return

  end subroutine mpi_comm_rank


  subroutine mpi_comm_size( comm, nprocs, ierror )

    !*****************************************************************************80
    !
    !! MPI_COMM_SIZE reports the number of processes in a communicator.
    !
    !  Discussion:
    !
    !    The routine simply returns NPROCS = 1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: comm
    integer(i4) :: ierror
    integer(i4) :: nprocs

    ierror = MPI_SUCCESS
    nprocs = 1

    return

  end subroutine mpi_comm_size


  subroutine mpi_comm_split( comm, icolor, ikey, comm_new, ierror )

    !*****************************************************************************80
    !
    !! MPI_COMM_SPLIT splits up a communicator based on a key.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: comm
    integer(i4) :: comm_new
    integer(i4) :: icolor
    integer(i4) :: ierror
    integer(i4) :: ikey

    ierror = MPI_SUCCESS

    return

  end subroutine mpi_comm_split


  subroutine mpi_copy_double_precision( data1, data2, n, ierror )

    !*****************************************************************************80
    !
    !! MPI_COPY_DOUBLE copies a double precision vector.
    !
    !  Discussion:
    !
    !    This routine is not part of the MPI standard.  However, it is
    !    needed by other routines which do emulate standard MPI routines.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, double precision DATA1(N), the data to be copied.
    !
    !    Output, double precision DATA2(N), the copied data.
    !
    !    Input, integer N, the number of items of data.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    real(dp) :: data1(n)
    real(dp) :: data2(n)
    integer(i4) :: ierror

    ierror = MPI_SUCCESS

    data2(1:n) = data1(1:n)

    return

  end subroutine mpi_copy_double_precision


  subroutine mpi_copy_integer( data1, data2, n, ierror )

    !*****************************************************************************80
    !
    !! MPI_COPY_INTEGER copies an integer vector.
    !
    !  Discussion:
    !
    !    This routine is not part of the MPI standard.  However, it is
    !    needed by other routines which do emulate standard MPI routines.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer DATA1(N), the data to be copied.
    !
    !    Output, integer DATA2(N), the copied data.
    !
    !    Input, integer N, the number of items of data.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: data1(n)
    integer(i4) :: data2(n)
    integer(i4) :: ierror

    ierror = MPI_SUCCESS

    data2(1:n) = data1(1:n)

    return

  end subroutine mpi_copy_integer


  subroutine mpi_copy_real( data1, data2, n, ierror )

    !*****************************************************************************80
    !
    !! MPI_COPY_REAL copies a real vector.
    !
    !  Discussion:
    !
    !    This routine is not part of the MPI standard.  However, it is
    !    needed by other routines which do emulate standard MPI routines.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, real DATA1(N), the data to be copied.
    !
    !    Output, real DATA2(N), the copied data.
    !
    !    Input, integer N, the number of items of data.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    real(sp) :: data1(n)
    real(sp) :: data2(n)
    integer(i4) :: ierror

    ierror = MPI_SUCCESS

    data2(1:n) = data1(1:n)

    return

  end subroutine mpi_copy_real


  subroutine mpi_finalize( ierror )

    !*****************************************************************************80
    !
    !! MPI_FINALIZE shuts down the MPI library.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: ierror

    ierror = MPI_SUCCESS

    return

  end subroutine mpi_finalize


  subroutine mpi_get_count( istatus, datatype, icount, ierror )

    !*****************************************************************************80
    !
    !! MPI_GET_COUNT reports the actual number of items transmitted.
    !
    !  Discussion:
    !
    !    Warn against querying message from self, since no data copy is done.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: datatype
    integer(i4) :: icount
    integer(i4) :: ierror
    integer(i4) :: istatus

    ierror = MPI_FAILURE

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_GET_COUNT - Error!'
    write ( *, '(a)' ) '  Should not query message from self.'

    return

  end subroutine mpi_get_count


  subroutine mpi_init( ierror )

    !*****************************************************************************80
    !
    !! MPI_INIT initializes the MPI library.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: ierror

    ierror = MPI_SUCCESS

    return

  end subroutine mpi_init


  subroutine mpi_irecv( data, n, datatype, iproc, itag, comm, irequest, ierror )

    !*****************************************************************************80
    !
    !! MPI_IRECV receives data from another process.
    !
    !  Discussion:
    !
    !    Warn against receiving message from self, since no data copy is done.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    integer(i4) :: data(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: iproc
    integer(i4) :: irequest
    integer(i4) :: itag

    ierror = MPI_FAILURE

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_IRECV - Error!'
    write ( *, '(a)' ) '  Should not "recv" message from self.'

    return

  end subroutine mpi_irecv


  subroutine mpi_isend( data, n, datatype, iproc, itag, comm, request, ierror )

    !*****************************************************************************80
    !
    !! MPI_ISEND sends data to another process using nonblocking transmission.
    !
    !  Discussion:
    !
    !    Warn against sending message to self, since no data copy is done.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    15 August 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, datatype DATA(N), the data to be sent.
    !
    !    Input, integer N, the number of data items to send.
    !
    !    Input, integer DATAYTPE, the MPI code for the datatype.
    !
    !    Input, integer IPROC, the rank of the process within the communicator
    !    that is to receive the message.
    !
    !    Input, integer ITAG, a tag for the message.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer REQUEST, a handle.  To determine if the data has been
    !    received yet, call MPI_Test or MPI_Wait, including the value of REQUEST.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    integer(i4) :: data(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: iproc
    integer(i4) :: itag
    integer(i4) :: request

    request = 0
    ierror = MPI_FAILURE

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_ISEND - Error!'
    write ( *, '(a)' )  '  Should not send message to self.'

    return

  end subroutine mpi_isend


  subroutine mpi_recv( data, n, datatype, iproc, itag, comm, istatus, ierror )

    !*****************************************************************************80
    !
    !! MPI_RECV receives data from another process within a communicator.
    !
    !  Discussion:
    !
    !    Warn against receiving message from self, since no data copy is done.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    integer(i4) :: data(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: iproc
    integer(i4) :: istatus
    integer(i4) :: itag

    ierror = MPI_FAILURE

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_RECV - Error!'
    write ( *, '(a)' ) '  Should not "recv" message from self.'

    return

  end subroutine mpi_recv


  subroutine mpi_reduce_0d_i4( data1, data2, n, datatype, operation, receiver, &
       comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_REDUCE carries out a reduction operation.
    !
    !  Discussion:
    !
    !    The reduction operations are sum, maximum, minimum, product.
    !
    !    The first two arguments must not overlap or share memory in any way.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !    Thanks to Simppa Akaslompolo for correcting this routine!
    !    11 January 2012.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, DATATYPE DATA1(N), the data to be processed.
    !
    !    Output (to RECEIVER only), DATATYPE DATA2(N), the value of the
    !    reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Input, integer RECEIVER, the process that is to receive the
    !    result.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    integer(i4) :: data1
    integer(i4) :: data2
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: operation
    integer(i4) :: receiver

    if ( datatype == mpi_integer ) then
       call mpi_reduce_0d_integer ( data1, data2, n, operation, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_reduce_0d_i4


  subroutine mpi_reduce_1d_i4( data1, data2, n, datatype, operation, receiver, &
       comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_REDUCE carries out a reduction operation.
    !
    !  Discussion:
    !
    !    The reduction operations are sum, maximum, minimum, product.
    !
    !    The first two arguments must not overlap or share memory in any way.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !    Thanks to Simppa Akaslompolo for correcting this routine!
    !    11 January 2012.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, DATATYPE DATA1(N), the data to be processed.
    !
    !    Output (to RECEIVER only), DATATYPE DATA2(N), the value of the
    !    reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Input, integer RECEIVER, the process that is to receive the
    !    result.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    integer(i4) :: data1(n)
    integer(i4) :: data2(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: operation
    integer(i4) :: receiver

    if ( datatype == mpi_integer ) then
       call mpi_reduce_integer ( data1, data2, n, operation, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_reduce_1d_i4


  subroutine mpi_reduce_sp( data1, data2, n, datatype, operation, receiver, &
       comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_REDUCE carries out a reduction operation.
    !
    !  Discussion:
    !
    !    The reduction operations are sum, maximum, minimum, product.
    !
    !    The first two arguments must not overlap or share memory in any way.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !    Thanks to Simppa Akaslompolo for correcting this routine!
    !    11 January 2012.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, DATATYPE DATA1(N), the data to be processed.
    !
    !    Output (to RECEIVER only), DATATYPE DATA2(N), the value of the
    !    reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Input, integer RECEIVER, the process that is to receive the
    !    result.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    real(sp) :: data1(n)
    real(sp) :: data2(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: operation
    integer(i4) :: receiver

    if ( datatype == mpi_real ) then
       call mpi_reduce_real ( data1, data2, n, operation, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_reduce_sp


  subroutine mpi_reduce_dp( data1, data2, n, datatype, operation, receiver, &
       comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_REDUCE carries out a reduction operation.
    !
    !  Discussion:
    !
    !    The reduction operations are sum, maximum, minimum, product.
    !
    !    The first two arguments must not overlap or share memory in any way.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !    Thanks to Simppa Akaslompolo for correcting this routine!
    !    11 January 2012.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, DATATYPE DATA1(N), the data to be processed.
    !
    !    Output (to RECEIVER only), DATATYPE DATA2(N), the value of the
    !    reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Input, integer RECEIVER, the process that is to receive the
    !    result.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    real(dp) :: data1(n)
    real(dp) :: data2(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: operation
    integer(i4) :: receiver

    if ( datatype == mpi_double_precision ) then
       call mpi_reduce_double_precision ( data1, data2, n, operation, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_reduce_dp


  subroutine mpi_reduce_double_precision( data1, data2, n, operation, ierror )

    !*****************************************************************************80
    !
    !! MPI_REDUCE_DOUBLE_PRECISION: reduction operation on double precision values.
    !
    !  Discussion:
    !
    !    The reduction operations are sum, maximum, minimum, product.
    !
    !    Thanks to Simppa Akaslompolo for correcting this routine!
    !    11 January 2012.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, double precision DATA1(N), the data to be processed.
    !
    !    Output, double precision DATA2(N), the value of the reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    real(dp) :: data1(n)
    real(dp) :: data2(n)
    integer(i4) :: ierror
    integer(i4) :: operation

    ierror = MPI_SUCCESS

    if ( operation == mpi_max ) then

       data2(1:n) = data1(1:n)

    else if ( operation == mpi_min ) then

       data2(1:n) = data1(1:n)

    else if ( operation == mpi_product ) then

       data2(1:n) = data1(1:n)

    else if ( operation == mpi_sum ) then

       data2(1:n) = data1(1:n)

    else

       ierror = MPI_FAILURE

    end if

    return

  end subroutine mpi_reduce_double_precision


  subroutine mpi_reduce_0d_integer( data1, data2, n, operation, ierror )

    !*****************************************************************************80
    !
    !! MPI_REDUCE_0d_INTEGER carries out a reduction operation on integers.
    !
    !  Discussion:
    !
    !    The reduction operations are sum, maximum, minimum, product.
    !
    !    Thanks to Simppa Akaslompolo for correcting this routine!
    !    11 January 2012.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer DATA1, the data to be processed.
    !
    !    Output, integer DATA2, the value of the reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: data1
    integer(i4) :: data2
    integer(i4) :: ierror
    integer(i4) :: operation

    ierror = MPI_SUCCESS

    if ( operation == mpi_max ) then

       data2 = data1

    else if ( operation == mpi_min ) then

       data2 = data1

    else if ( operation == mpi_product ) then

       data2 = data1

    else if ( operation == mpi_sum ) then

       data2 = data1

    else

       ierror = MPI_FAILURE

    end if

    return

  end subroutine mpi_reduce_0d_integer


  subroutine mpi_reduce_integer( data1, data2, n, operation, ierror )

    !*****************************************************************************80
    !
    !! MPI_REDUCE_INTEGER carries out a reduction operation on integers.
    !
    !  Discussion:
    !
    !    The reduction operations are sum, maximum, minimum, product.
    !
    !    Thanks to Simppa Akaslompolo for correcting this routine!
    !    11 January 2012.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer DATA1(N), the data to be processed.
    !
    !    Output, integer DATA2(N), the value of the reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: data1(n)
    integer(i4) :: data2(n)
    integer(i4) :: ierror
    integer(i4) :: operation

    ierror = MPI_SUCCESS

    if ( operation == mpi_max ) then

       data2(1:n) = data1(1:n)

    else if ( operation == mpi_min ) then

       data2(1:n) = data1(1:n)

    else if ( operation == mpi_product ) then

       data2(1:n) = data1(1:n)

    else if ( operation == mpi_sum ) then

       data2(1:n) = data1(1:n)

    else

       ierror = MPI_FAILURE

    end if

    return

  end subroutine mpi_reduce_integer


  subroutine mpi_reduce_real( data1, data2, n, operation, ierror )

    !*****************************************************************************80
    !
    !! MPI_REDUCE_REAL carries out a reduction operation on reals.
    !
    !  Discussion:
    !
    !    The reduction operations are sum, maximum, minimum, product.
    !
    !    Thanks to Simppa Akaslompolo for correcting this routine!
    !    11 January 2012.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, real DATA1(N), the data to be processed.
    !
    !    Output, real DATA2(N), the value of the reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    real(sp) :: data1(n)
    real(sp) :: data2(n)
    integer(i4) :: ierror
    integer(i4) :: operation

    ierror = MPI_SUCCESS

    if ( operation == mpi_max ) then

       data2(1:n) = data1(1:n)

    else if ( operation == mpi_min ) then

       data2(1:n) = data1(1:n)

    else if ( operation == mpi_product ) then

       data2(1:n) = data1(1:n)

    else if ( operation == mpi_sum ) then

       data2(1:n) = data1(1:n)

    else

       ierror = MPI_FAILURE

    end if

    return

  end subroutine mpi_reduce_real


  subroutine mpi_reduce_scatter_i4( data1, data2, n, datatype, operation, comm, &
       ierror )

    !*****************************************************************************80
    !
    !! MPI_REDUCE_SCATTER collects a message of the same length from each process.
    !
    !  Discussion:
    !
    !    Copy values from DATA1 to DATA2.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, DATATYPE DATA1(N), the data to be processed.
    !
    !    Output, DATATYPE DATA2, the value of the reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    integer(i4) :: data1(n)
    integer(i4) :: data2(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: operation

    ierror = MPI_SUCCESS

    if ( datatype == mpi_integer ) then
       call mpi_copy_integer ( data1, data2, n, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_reduce_scatter_i4


  subroutine mpi_reduce_scatter_sp( data1, data2, n, datatype, operation, comm, &
       ierror )

    !*****************************************************************************80
    !
    !! MPI_REDUCE_SCATTER collects a message of the same length from each process.
    !
    !  Discussion:
    !
    !    Copy values from DATA1 to DATA2.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, DATATYPE DATA1(N), the data to be processed.
    !
    !    Output, DATATYPE DATA2, the value of the reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    real(sp) :: data1(n)
    real(sp) :: data2(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: operation

    ierror = MPI_SUCCESS

    if ( datatype == mpi_real ) then
       call mpi_copy_real ( data1, data2, n, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_reduce_scatter_sp


  subroutine mpi_reduce_scatter_dp( data1, data2, n, datatype, operation, comm, &
       ierror )

    !*****************************************************************************80
    !
    !! MPI_REDUCE_SCATTER collects a message of the same length from each process.
    !
    !  Discussion:
    !
    !    Copy values from DATA1 to DATA2.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, DATATYPE DATA1(N), the data to be processed.
    !
    !    Output, DATATYPE DATA2, the value of the reduction operation.
    !
    !    Input, integer N, the number of items in DATA1.
    !
    !    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
    !
    !    Input, integer OPERATION, should have the value of one of the symbolic
    !    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    real(dp) :: data1(n)
    real(dp) :: data2(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: operation

    ierror = MPI_SUCCESS

    if ( datatype == mpi_double_precision ) then
       call mpi_copy_double_precision ( data1, data2, n, ierror )
    else
       ierror = MPI_FAILURE
    end if

    return

  end subroutine mpi_reduce_scatter_dp


  subroutine mpi_rsend( data, n, datatype, iproc, itag, comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_RSEND "ready sends" data from one process to another.
    !
    !  Discussion:
    !
    !    Warn against sending message to self, since no data copy is done.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    integer(i4) :: data(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: iproc
    integer(i4) :: itag

    ierror = MPI_FAILURE

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_RSEND - Error!'
    write ( *, '(a)' ) '  Should not send message to self.'

    return

  end subroutine mpi_rsend


  subroutine mpi_send( data, n, datatype, iproc, itag, comm, ierror )

    !*****************************************************************************80
    !
    !! MPI_SEND sends data from one process to another.
    !
    !  Discussion:
    !
    !    Warn against sending message to self, since no data copy is done.
    !
    !    The data to be transferred can be integer, real, or double precision.
    !    In this routine, it is declared and documented as INTEGER type,
    !    but using the other types should generally not cause a problem.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Input, datatype DATA(N), the data to be sent.
    !
    !    Input, integer N, the number of data items to send.
    !
    !    Input, integer DATAYTPE, the MPI code for the datatype.
    !
    !    Input, integer IPROC, the rank of the process within the communicator
    !    that is to receive the message.
    !
    !    Input, integer ITAG, a tag for the message.
    !
    !    Input, integer COMM, the MPI communicator.
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: n

    integer(i4) :: comm
    integer(i4) :: data(n)
    integer(i4) :: datatype
    integer(i4) :: ierror
    integer(i4) :: iproc
    integer(i4) :: itag

    ierror = MPI_FAILURE

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_SEND - Error!'
    write ( *, '(a)' )  '  Should not send message to self.'

    return

  end subroutine mpi_send


  subroutine mpi_wait( irequest, istatus, ierror )

    !*****************************************************************************80
    !
    !! MPI_WAIT waits for an I/O request to complete.
    !
    !  Discussion:
    !
    !    Warn against waiting on message from self, since no data copy is done.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 December 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: ierror
    integer(i4) :: irequest
    integer(i4) :: istatus

    ierror = MPI_FAILURE

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_WAIT - Error!'
    write ( *, '(a)' ) '  Should not wait on message from self.'

    return

  end subroutine mpi_wait


  subroutine mpi_waitall( icount, irequest, istatus, ierror )

    !*****************************************************************************80
    !
    !! MPI_WAITALL waits until all I/O requests have completed.
    !
    !  Discussion:
    !
    !    Warn against waiting on message from self, since no data copy is done.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: icount
    integer(i4) :: ierror
    integer(i4) :: irequest
    integer(i4) :: istatus

    ierror = MPI_FAILURE

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_WAITALL - Error!'
    write ( *, '(a)' ) '  Should not wait on message from self.'

    return

  end subroutine mpi_waitall


  subroutine mpi_waitany( icount, array_of_requests, index, istatus, ierror )

    !*****************************************************************************80
    !
    !! MPI_WAITANY waits until one I/O requests has completed.
    !
    !  Discussion:
    !
    !    Warn against waiting on message from self, since no data copy is done.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Output, integer IERROR, is nonzero if an error occurred.
    !
    implicit none

    integer(i4) :: array_of_requests(*)
    integer(i4) :: icount
    integer(i4) :: ierror
    integer(i4) :: index
    integer(i4) :: istatus

    ierror = MPI_FAILURE

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_WAITANY - Error!'
    write ( *, '(a)' ) '  Should not wait on message from self.'

    return

  end subroutine mpi_waitany


  function mpi_wtick( )

    !*****************************************************************************80
    !
    !! MPI_WTICK returns the number of seconds per clock tick.
    !
    !  Discussion:
    !
    !    The value returned here is simply a dummy value.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 October 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Output, real(dp) :: MPI_WTICK, the number of seconds per clock tick.
    !
    implicit none

    real(dp) :: mpi_wtick

    mpi_wtick = 1.0D+00

    return

  end function mpi_wtick


  function mpi_wtime( )

    !*****************************************************************************80
    !
    !! MPI_WTIME returns the elapsed wall clock time.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 October 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    William Gropp, Ewing Lusk, Anthony Skjellum,
    !    Using MPI: Portable Parallel Programming with the
    !    Message-Passing Interface,
    !    Second Edition,
    !    MIT Press, 1999,
    !    ISBN: 0262571323.
    !
    !  Parameters:
    !
    !    Output, real(dp) :: MPI_WTIME, the elapsed wall clock time.
    !
    implicit none

    integer(i8) :: count
    integer(i8) :: count_max
    integer(i8) :: count_rate
    real(dp) :: mpi_wtime

    call system_clock ( count, count_rate, count_max )

    mpi_wtime = real(count,dp)/real(count_rate,dp)

    return

  end function mpi_wtime


end module mo_mpi_stubs

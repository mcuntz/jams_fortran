!> \file mo_file_utils.f90

!> \brief General file utilities

!> \details This module provides general utilities for handling files.

!> \authors Matthias Cuntz
!> \date Dec 2011
module mo_file_utils

  ! This module provides general utilities for handling files

  ! Written  Matthias Cuntz, Dec 2011
  ! Modified Matthias Cuntz, Jan 2017 - istart,istop optional in find_next_unit
  !          Matthias Cuntz, Jan 2017 - lines_in_file

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

  ! Copyright 2011-2017 Matthias Cuntz

  use mo_kind, only: i4

  implicit none

  private

  public :: find_next_unit    ! find file handle that is not used yet
  public :: lines_in_file     ! count number of lines in file

  ! ------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------

  !     NAME
  !         find_next_unit

  !     PURPOSE
  !>        \brief Find unused file unit.

  !>        \details Starting from a given number this routine checks in a range
  !>        if the numbers are already assigned to a file handle.
  !>        It reports the first number that is not associated yet.

  !     CALLING SEQUENCE
  !         iout = find_next_unit(istart,istop)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4) :: istart"     Starting unit (default: 20)
  !>        \param[in] "integer(i4) :: istop"      Last checked unit (default: 1000)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: find_next_unit ! Free unit in interval [istart,istop], returns -1 if no free unit.

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         iout = find_next_unit(100)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Written, Matthias Cuntz - modified from Echam5 mo_filename.f90/find_next_free_unit
  !>                                          Echam5, (C) MPI-MET, Hamburg, Germany
  !>        \date Dec 2011
  !         Modified, Matthias Cuntz, Jan 2017 - istart, istop optinal

  function find_next_unit(istart,istop)

    implicit none

    integer(i4), intent(in), optional :: istart, istop
    integer(i4)                       :: find_next_unit

    logical     :: opened
    integer(i4) :: i, ifirst, ilast

    ifirst = 20
    if (present(istart)) ifirst = istart
    ilast = 1000
    if (present(istop)) ilast = istop
    
    opened         = .true.
    find_next_unit = -1
    do i=ifirst, ilast
       inquire(unit=i, opened=opened)
       if (.not. opened) then
          find_next_unit = i
          exit
       end if
    end do

  end function find_next_unit

  ! ------------------------------------------------------------------

  !     NAME
  !         lines_in_file

  !     PURPOSE
  !>        \brief Count number of lines in a file.

  !>        \details Count number of lines in file,
  !>        optionally excluding blank lines and lines starting with comment character.

  !     CALLING SEQUENCE
  !         num = lines_in_file(filename, comment, noblank)

  !     INTENT(IN)
  !>        \param[in] "character(*) :: filename"     Name of file to count lines

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "character(*) :: comment"     Exclude lines starting with comment.
  !>                                                 Leading blanks are ignored.
  !>        \param[in] "logical :: noblank"          .false.:  count blank lines,
  !>                                                           i.e. lines, where trim(adjustl(line)) is epmty.
  !>                                                 .true.:   do not count blank lines,
  !>                                                           i.e. where trim(adjustl(line)) is epmty.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: lines_in_file    ! number of lines, without commented or blank lines (optional)

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         iout = lines_in_file(filename, comment='!', noblank=.true.)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Jan 2017

  function lines_in_file(filename, comment, noblank)

    implicit none

    character(len=*),           intent(in) :: filename
    character(len=*), optional, intent(in) :: comment
    logical,          optional, intent(in) :: noblank
    integer                                :: lines_in_file

    ! ToDo

    integer(i4), parameter :: maxlen = 8192
    logical     :: iexist, inoblank
    integer(i4) :: nun, ierr
    character(len=maxlen) :: iread, icomment

    ! options
    inoblank = .false.
    if (present(noblank)) inoblank = noblank
    
    ! file exists
    inquire(file=filename, exist=iexist)
    if (.not. iexist) then
       lines_in_file = -1
       return
    endif
    
    ! open
    nun = find_next_unit()
    open(nun, file=filename, status='old', action='read', form="formatted", recl=maxlen)
    ! count
    lines_in_file = 0
    ierr = 0 
    do while (ierr==0)
       read(nun, '(a)', iostat=ierr) iread
       if (ierr==0) then
          ! exclude commented lines
          if (present(comment)) then
             icomment = adjustl(iread)
             if (icomment(1:len(comment)) == comment) cycle
          endif
          ! exclude blank lines
          if (inoblank) then
             if (trim(adjustl(iread)) == '') cycle
          endif
          lines_in_file = lines_in_file + 1
       endif
    end do

    close(nun)

  end function lines_in_file

  ! ------------------------------------------------------------------

END MODULE mo_file_utils

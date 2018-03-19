!> \file mo_tee.f90

!> \brief Writes message to standard output and into file; similar to Unix tee utility.

!> \details Mimicks the Unix tee utility. Writes message on standard out and
!> into a file. File can be given via a unit or a filename, in which case appends
!> to the file or optionally overwrites it.

!> \authors Matthias Cuntz
!> \date Mar 2018
module mo_tee

  ! This module provides the tee utility

  ! Written  Matthias Cuntz, Mar 2018

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

  ! Copyright 2018 Matthias Cuntz

  implicit none

  public :: tee ! write to standard out and into file

  ! ------------------------------------------------------------------

  !     NAME
  !         tee

  !     PURPOSE
  !>        \brief Writes message to standard out and into file; similar to Unix tee utility.

  !>        \details Mimicks the Unix tee utility. Writes message on standard out and
  !>        into a file. File can be given as a unit or as a filename, in which case tee
  !>        appends to the file, or optionally overwrites the file.

  !     CALLING SEQUENCE
  !         call tee(file, string, overwrite)

  !     INTENT(IN)
  !>        \param[in] "integer(i4)/character(len=*) :: file"   Unit or filename
  !>        \param[in] "character(len=*) :: string"   String to be written to standard out and file

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical :: overwrite"   .true.: overwrite existing file if file is filename;
  !>                                            .false.: append to existing file if file is filename,
  !>                                                     otherwise create a new file
  !>                                            (default: .false.)

  !     EXAMPLE
  !         call tee(logfile, 'I am a message.')
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Mar 2018
  interface tee
     module procedure tee_unit, tee_filename
  end interface tee

  private

  ! ------------------------------------------------------------------

contains

  subroutine tee_unit(unit, string)

    use mo_kind, only: i4

    implicit none

    integer(i4),      intent(in) :: unit
    character(len=*), intent(in) :: string

    write(unit,*) string
    write(*,*) string

  end subroutine tee_unit


  subroutine tee_filename(filename, string, overwrite)

    use mo_kind,       only: i4
    use mo_file_utils, only: find_next_unit

    implicit none

    character(len=*),  intent(in) :: filename
    character(len=*),  intent(in) :: string
    logical, optional, intent(in) :: overwrite

    integer(i4) :: nun
    logical     :: over

    ! optionals
    over = .false.
    if (present(overwrite)) over = overwrite

    ! open file accordinf to optionals
    nun = find_next_unit()
    if (over) then
       open(unit=nun, file=filename, action='write', form="formatted", &
            status='replace')
    else
       open(unit=nun, file=filename, action='write', form="formatted", &
            status='unknown', position='append')
    endif

    ! tee
    write(nun,*) string
    write(*,*) string

    close(nun)

  end subroutine tee_filename

  ! ------------------------------------------------------------------

END MODULE mo_tee

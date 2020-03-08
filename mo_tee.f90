!> \file mo_tee.f90

!> \brief Write out concatenated strings on standard out and to a given file or unit.

!> \details Write out several strings concatenated on standard out and to a given file or unit,
!> advancing or not, overwriting the file or appending.
!> Similar to the Unix tee utility.

!> \authors Matthias Cuntz
!> \date Mar 2018
module mo_tee

  ! This module provides the tee utility.

  ! Written  Matthias Cuntz, Mar 2018

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2018 Matthias Cuntz - mc (at) macu (dot) de
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

  implicit none

  public :: tee ! write to standard out and into file

  ! ------------------------------------------------------------------

  !     NAME
  !         tee

  !     PURPOSE
  !>        \brief Write out concatenated strings on standard out and to a given file or unit.

  !>        \details Write out several strings concatenated on standard out and to a given file or unit.
  !>        It is similar to the Unix utility tee.
  !>        File can be given as a unit or as a filename, in which case tee appends to the file,
  !>        or optionally overwrites the file.
  !>        The concatenated string can be written advancing or not.

  !     CALLING SEQUENCE
  !         call tee(file, t01=' ', t02='', t03='', t04='', t05='', t06='', t07='', &
  !                  t08='', t09='', t10='', advance='yes', overwrite=.false.)

  !     INTENT(IN)
  !>        \param[in] "integer(i4)/character(len=*) :: file"     Unit or filename

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "character(len=*), optional :: t01"        1st string (default: ' ')
  !>        \param[in] "character(len=*), optional :: t02"        2nd string (default: '')
  !>        \param[in] "character(len=*), optional :: t03"        3rd string (default: '')
  !>        \param[in] "character(len=*), optional :: t04"        4th string (default: '')
  !>        \param[in] "character(len=*), optional :: t05"        5th string (default: '')
  !>        \param[in] "character(len=*), optional :: t06"        6th string (default: '')
  !>        \param[in] "character(len=*), optional :: t07"        7th string (default: '')
  !>        \param[in] "character(len=*), optional :: t08"        8th string (default: '')
  !>        \param[in] "character(len=*), optional :: t09"        9th string (default: '')
  !>        \param[in] "character(len=*), optional :: t10"        10th string (default: '')
  !>        \param[in] "integer         , optional :: unit"       Unit to write to (default: *)
  !>        \param[in] "character(len=*), optional :: advance"    advance keyword of write (default: 'yes')\n
  !>                                                              'yes': newline will be written after message\n
  !>                                                              'no':  no newline at end of message
  !>        \param[in] "logical, optional :: overwrite"           .true.: overwrite existing file if file is filename;\n
  !>                                                              .false.: append to existing file if file is filename,\n
  !>                                                              otherwise create a new file (default: .false.)

  !     EXAMPLE
  !         call tee(logfile, 'I am ', 'a message.')
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

  subroutine tee_unit(unit, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, advance)

    use mo_kind,    only: i4
    use mo_message, only: message

    implicit none

    integer(i4),      intent(in)           :: unit
    character(len=*), intent(in), optional :: t01
    character(len=*), intent(in), optional :: t02
    character(len=*), intent(in), optional :: t03
    character(len=*), intent(in), optional :: t04
    character(len=*), intent(in), optional :: t05
    character(len=*), intent(in), optional :: t06
    character(len=*), intent(in), optional :: t07
    character(len=*), intent(in), optional :: t08
    character(len=*), intent(in), optional :: t09
    character(len=*), intent(in), optional :: t10
    character(len=*), intent(in), optional :: advance

    call message(t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, advance=advance)
    call message(t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, unit, advance)

  end subroutine tee_unit


  subroutine tee_filename(filename, t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, advance, overwrite)

    use mo_kind,       only: i4
    use mo_file_utils, only: find_next_unit
    use mo_message,    only: message

    implicit none

    character(len=*), intent(in)           :: filename
    character(len=*), intent(in), optional :: t01
    character(len=*), intent(in), optional :: t02
    character(len=*), intent(in), optional :: t03
    character(len=*), intent(in), optional :: t04
    character(len=*), intent(in), optional :: t05
    character(len=*), intent(in), optional :: t06
    character(len=*), intent(in), optional :: t07
    character(len=*), intent(in), optional :: t08
    character(len=*), intent(in), optional :: t09
    character(len=*), intent(in), optional :: t10
    character(len=*), intent(in), optional :: advance
    logical,          intent(in), optional :: overwrite

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
    call message(t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, advance=advance)
    call message(t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, nun, advance)

    close(nun)

  end subroutine tee_filename

  ! ------------------------------------------------------------------

END MODULE mo_tee

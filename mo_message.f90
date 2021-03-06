!> \file mo_message.f90

!> \brief Write out concatenated strings

!> \details Write out several strings concatenated on standard out or a given unit, either advancing or not.

!> \author Matthias Cuntz
!> \date Jul 2011

module mo_message

  ! This module supplies routines to write out text

  ! Written  Jul 2011, Matthias Cuntz - Inspired from Echam5 mo_exception.f90
  ! Modified Mar 2018, Matthias Cuntz - Use * instead of nout=6 for standard out

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2011-2018 Matthias Cuntz - mc (at) macu (dot) de
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

  private

  public :: message_text    ! dummy string to use in subroutines
  public :: message         ! versatile routine to write out strings in file or on screen

  character(len=1024) :: message_text = ''

  ! ------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------

  !     NAME
  !         message

  !     PURPOSE
  !>        \brief Write out several string concatenated either on screen or in a file.

  !     CALLING SEQUENCE
  !         call message(t01=' ', t02='', t03='', t04='', t05='', t06='', t07='', &
  !                      t08='', t09='', t10='', unit=unit, advance='yes')

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

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

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         "Only" 10 input strings

  !     EXAMPLE
  !         call message('A=',advance='no')
  !         call message(num2str(a))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz - modified from Echam5, (C) MPI-MET, Hamburg, Germany
  !>        \date Dec 2011
  !         Modified Matthias Cuntz, Mar 2018 - use * for standard out

  subroutine message(t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, uni, advance)

    use mo_kind, only: i4
    
    implicit none

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
    integer,          intent(in), optional :: uni
    character(len=*), intent(in), optional :: advance

    integer(i4)          :: iout
    character(len=32000) :: out
    character(len=3)     :: iadv
#ifdef __GFORTRAN__
    character(len=32000) :: nold
#endif

    iout = -1
    if (present(uni)) iout = uni
    if (present(advance)) then
       iadv = ''
       iadv(1:min(len(advance),3)) = advance(1:min(len(advance),3))
    else
       iadv = 'yes'
    end if

    out = ''
    ! start from back so that trim does not remove user desired blanks
#ifdef __GFORTRAN__
    ! GFORTRAN has problems with concatenation operator //
    ! It is also weird in write:
    !    write(out,'(A,A)') t10, trim(out)
    ! writes t10 twice into out.
    nold = out
    if (present(t10)) write(out,'(a,a)') t10, trim(nold)
    nold = out
    if (present(t09)) write(out,'(a,a)') t09, trim(nold)
    nold = out
    if (present(t08)) write(out,'(a,a)') t08, trim(nold)
    nold = out
    if (present(t07)) write(out,'(a,a)') t07, trim(nold)
    nold = out
    if (present(t06)) write(out,'(a,a)') t06, trim(nold)
    nold = out
    if (present(t05)) write(out,'(a,a)') t05, trim(nold)
    nold = out
    if (present(t04)) write(out,'(a,a)') t04, trim(nold)
    nold = out
    if (present(t03)) write(out,'(a,a)') t03, trim(nold)
    nold = out
    if (present(t02)) write(out,'(a,a)') t02, trim(nold)
    nold = out
    if (present(t01)) write(out,'(a,a)') t01, trim(nold)
    ! output at least one space otherwise some compilers get confused on Mac (empty assembler statement)
    if ((lle(trim(out),'') .and. lge(trim(out),''))) then
       nold = out
       write(out,'(a,a)') trim(nold), ' '
    endif
    if (iout > 0) then
       write(iout,'(a)',advance=iadv) trim(out)
    else
       write(*,'(a)',advance=iadv) trim(out)
    endif
#else
    if (present(t10)) out = t10//trim(out)
    if (present(t09)) out = t09//trim(out)
    if (present(t08)) out = t08//trim(out)
    if (present(t07)) out = t07//trim(out)
    if (present(t06)) out = t06//trim(out)
    if (present(t05)) out = t05//trim(out)
    if (present(t04)) out = t04//trim(out)
    if (present(t03)) out = t03//trim(out)
    if (present(t02)) out = t02//trim(out)
    if (present(t01)) out = t01//trim(out)
    ! output at least one space otherwise some compilers get confused on Mac (empty assembler statement)
    if ((lle(trim(out),'') .and. lge(trim(out),''))) then
       if (iout > 0) then
          write(iout,'(a)',advance=iadv) trim(out)//' '
       else
          write(*,'(a)',advance=iadv) trim(out)//' '
       endif
    else
       if (iout > 0) then
          write(iout,'(a)',advance=iadv) trim(out)
       else
          write(*,'(a)',advance=iadv) trim(out)
       endif
    endif
#endif

  end subroutine message

end module mo_message

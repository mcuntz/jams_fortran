!> \file mo_finish.f90

!> \brief Finish a program gracefully

!> \details This module supplies a routine that writes out final comments and then stops.

!> \authors Original of Echam5, (C) MPI-MET, Hamburg, Germany.\n
!> Modified Matthias Cuntz
!> \date Jan 2011

MODULE mo_finish

  ! This module supplies a routine to finish a program gracefully.

  ! Written Jan 2011, Matthias Cuntz - Adapted from Echam5, (C) MPI-MET, Hamburg, Germany
  ! Modified May 2016, Matthias Cuntz - copied separator from mo_string_utils

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2011-2016 Matthias Cuntz, MPI-MET Hamburg Germany - mc (at) macu (dot) de
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

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: finish     ! Write out error message and stop program

  ! ------------------------------------------------------------------

  CHARACTER(len=*), PARAMETER :: separator = repeat('-',70)

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         finish

  !     PURPOSE
  !>        \brief Finish a program gracefully

  !>        \details Stop a program but writing out a message first that is separated
  !>        from earlier output by -------------- (i.e. the separator)

  !     CALLING SEQUENCE
  !         call finish(name, text=text, unit=unit)

  !     INTENT(IN)
  !>        \param[in] "character(len=*) :: name"         First string separated from otional second by :

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "character(len=*), optional :: text"         Second string separated by :
  !>        \param[in] "integer, optional          :: unit"         File unit for write (default: *)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         call finish('main','End gracefully')
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Written, Matthias Cuntz - modified from Echam5, (C) MPI-MET, Hamburg, Germany
  !>        \date Dec 2011

  SUBROUTINE finish(name, text, unit)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN)           :: name
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: text
    INTEGER,          INTENT(IN), OPTIONAL :: unit

    IF (PRESENT(unit)) THEN
       WRITE (unit,'(a)') separator

       IF (PRESENT(text)) THEN
          WRITE (unit,'(a,a,a)') name, ': ', text
       ELSE
          WRITE (unit,'(a)') name
       END IF

       WRITE (unit,'(a)') separator
    ELSE
       WRITE (*,'(a)') separator

       IF (PRESENT(text)) THEN
          WRITE (*,'(a,a,a)') name, ': ', text
       ELSE
          WRITE (*,'(a)') name
       END IF

       WRITE (*,'(a)') separator
    ENDIF

    STOP

  END SUBROUTINE finish

END MODULE mo_finish

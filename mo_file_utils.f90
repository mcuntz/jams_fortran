MODULE mo_file_utils

  ! This module provides general utilities for handling files

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: find_next_unit    ! find file handle that is not used yet

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         find_next_unit

  !     PURPOSE
  !         Starting from a given number this routine checks in a range if the numbers
  !         are already assigned to a file handle. It reports the first number that is
  !         not associated yet.

  !     CALLING SEQUENCE
  !         iout = find_next_unit(istart,istop)
  
  !     INDENT(IN)
  !         integer :: istart         Starting unit
  !         integer :: istop          Last checked unit

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         integer :: iout           Free unit in interval [istart,istop], returns -1 if no free unit.

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         iout = find_next_unit(100,200)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011 - modified from Echam5 mo_filename.f90/find_next_free_unit
  !                                              Echam5, (C) MPI-MET, Hamburg, Germany

  FUNCTION find_next_unit(istart,istop)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: istart, istop
    INTEGER             :: find_next_unit

    LOGICAL :: opened
    INTEGER :: i

    opened         = .true.
    find_next_unit = -1
    DO i=istart, istop
       INQUIRE(unit=i, opened=opened)
       IF (.NOT. opened) THEN
          find_next_unit = i
          EXIT
       END IF
    END DO

  END FUNCTION find_next_unit

  ! ------------------------------------------------------------------

END MODULE mo_file_utils

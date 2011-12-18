MODULE mo_finish

  ! This module supplies a routine to finish a program gracefully.

  ! Written Jan 2011, Matthias Cuntz - Adapted from Echam5, (C) MPI-MET, Hamburg, Germany

  USE mo_constants,    ONLY: nerr
  USE mo_string_utils, ONLY: separator

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: finish     ! Write out error message and stop program

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         finish

  !     PURPOSE
  !         Stop a program but writing out a message first that is separated
  !         from earlier output by -------------- (i.e. the separator of mo_sting_utils)

  !     CALLING SEQUENCE
  !         call finish(name, text=text, unit=unit)
  
  !     INDENT(IN)
  !         character(len=*) :: name         First string separated from otional second by :

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         character(len=*) :: text         Second string separated by :
  !         integer          :: unit         File unit for write (default: nerr)

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         call finish('main','End gracefully')
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011 - modified from Echam5, (C) MPI-MET, Hamburg, Germany

  SUBROUTINE finish(name, text, unit)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN)           :: name
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: text
    INTEGER,          INTENT(IN), OPTIONAL :: unit

    INTEGER :: nunit
    
    IF (PRESENT(unit)) THEN
       nunit = unit
    ELSE
       nunit = nerr
    ENDIF
    
    WRITE (nunit,'(a)') separator

    IF (PRESENT(text)) THEN
      WRITE (nunit,'(a,a,a)') name, ': ', text
    ELSE
      WRITE (nunit,'(a)') name
    END IF

    WRITE (nunit,'(a)') separator

    STOP

  END SUBROUTINE finish

END MODULE mo_finish

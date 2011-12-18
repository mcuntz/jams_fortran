MODULE mo_message

  ! This module supplies routines to write out text

  ! Written Jul 2011, Matthias Cuntz - Inspired from Echam5 mo_exception.f90

  USE mo_constants, ONLY: nout

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: message_text    ! dummy string to use in subroutines
  PUBLIC :: message         ! versatile routine to write out strings in file or on screen

  CHARACTER(len=1024) :: message_text = ''

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         message

  !     PURPOSE
  !         Write out several string concatenated either on screen or in a file.

  !     CALLING SEQUENCE
  !         call message(t01=t01, t02=t02, t03=t03, t04=t04, t05=t05, t06=t06, t07=t07, &
  !                      t08=t08, t09=t09, t10=t10, unit=unit, advance=advance)
  
  !     INDENT(IN)
  !         None

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         character(len=*) :: t01        1st string
  !         character(len=*) :: t02        2nd string
  !         character(len=*) :: t03        3rd string
  !         character(len=*) :: t04        4th string
  !         character(len=*) :: t05        5th string
  !         character(len=*) :: t06        6th string
  !         character(len=*) :: t07        7th string
  !         character(len=*) :: t08        8th string
  !         character(len=*) :: t09        9th string
  !         character(len=*) :: t10        10th string
  !         integer          :: unit       Unit to write to (default: nout)
  !         character(len=*) :: advance    WRITE advance keyword (default: 'yes')
  !                                        yes: newline will be written after message
  !                                        no:  no newline at end of message

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
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
  !         Written,  Matthias Cuntz, Dec 2011 - modified from Echam5, (C) MPI-MET, Hamburg, Germany

  SUBROUTINE message(t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, uni, advance)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t01
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t02
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t03
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t04
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t05
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t06
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t07
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t08
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t09
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t10
    INTEGER,          INTENT(IN), OPTIONAL :: uni
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: advance

    INTEGER              :: iout
    CHARACTER(len=32000) :: out
    CHARACTER(len=3)     :: iadv


    if (present(uni)) then
       iout = uni
    else
       iout = nout
    end if
    if (present(advance)) then
       iadv = ''
       iadv(1:min(len(advance),3)) = advance(1:min(len(advance),3))
    else
       iadv = 'yes'
    end if

    out = ''
    ! start from back so that trim does not remove user desired blanks
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
    write(iout,'(a)',advance=iadv) trim(out)//' '

  END SUBROUTINE message

END MODULE mo_message

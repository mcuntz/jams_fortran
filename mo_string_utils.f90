MODULE mo_string_utils

  ! This module holds string conversion utilities

  USE mo_kind, ONLY: i4, i8, sp, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nonull    ! Check if string is still NULL
  PUBLIC :: num2str   ! Convert a number to a string
  PUBLIC :: separator ! Format string: '-----...-----'
  PUBLIC :: tolower   ! Conversion   : 'ABCXYZ' -> 'abcxyz'   
  PUBLIC :: toupper   ! Conversion   : 'abcxyz' -> 'ABCXYZ'

  INTERFACE num2str
     MODULE PROCEDURE i42str, i82str, sp2str, dp2str, log2str
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: separator = repeat('-',70)

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         nonull

  !     PURPOSE
  !         Checks if string was already used, i.e. does not contain NULL character anymore.

  !     CALLING SEQUENCE
  !         used = nonull(str)
  
  !     INDENT(IN)
  !         character(len=*) :: str    String

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         logical :: used    .true.: string was already set; .false.: string still in initialised state

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         if (nonull(str)) write(*,*) trim(str)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Jan 2012

  FUNCTION nonull(str)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(in) :: str
    LOGICAL                      :: nonull

    if (scan(str, achar(0)) == 0) then
       nonull = .true.
    else
       nonull = .false.
    endif

  END FUNCTION nonull

  ! ------------------------------------------------------------------

  !     NAME
  !         num2str

  !     PURPOSE
  !         Convert a number or logical to a string with an optional format.

  !     CALLING SEQUENCE
  !         str = num2str(num,form=form)
  
  !     INDENT(IN)
  !         integer(i4/i8)/real(sp/dp)/logical :: num    Number or logical

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         character(len=X) :: str              String of formatted input number or logical
  !                                              Ouput length X is:
  !                                              i4    - 10
  !                                              i8    - 20
  !                                              sp/dp - 32
  !                                              log   - 10

  !     INDENT(IN), OPTIONAL
  !         character(len=*) :: form             Format string
  !                                              Defaults are:
  !                                              i4    - '(I10)'
  !                                              i8    - '(I20)'
  !                                              sp/dp - '(G32.5)'
  !                                              log   - '(L10)'

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Uses WRITE to write into string. Recursive write is not permitted before Fortran 2003
  !         so that one cannot use
  !             write(*,*) 'A='//num2str(a)
  !         Use 'call message' from mo_messages.f90
  !             use mo_messages, only message
  !             call message('A=', trim(num2str(a)))
  !         or write into another string first:
  !             str = 'A='//num2str(a)
  !             write(*,*) trim(str)

  !     EXAMPLE
  !         str = num2str(3.1415217_i4,'(F3.1)')
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011 - modified from Echam5, (C) MPI-MET, Hamburg, Germany

  PURE FUNCTION i42str(nn,form)
    ! returns integer nn as a string (often needed in printing messages)
    IMPLICIT NONE
    INTEGER(i4),      INTENT(IN)           :: nn
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=10) :: i42str

    if (present(form)) then
       write(i42str,form) nn
    else
       write(i42str,'(I10)') nn
    end if
    !i42str = adjustl(i42str)

  END FUNCTION i42str


  PURE FUNCTION i82str(nn,form)
    ! returns integer nn as a string (often needed in printing messages)
    IMPLICIT NONE
    INTEGER(i8),      INTENT(IN)           :: nn
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=20) :: i82str

    if (present(form)) then
       write(i82str,form) nn
    else
       write(i82str,'(I20)') nn
    end if
    !i82str = adjustl(i82str)

  END FUNCTION i82str


  PURE FUNCTION sp2str(rr,form)
    ! returns real rr as a string (often needed in printing messages)
    IMPLICIT NONE
    REAL(sp),         INTENT(IN)           :: rr
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=32) :: sp2str

    if (present(form)) then
       write(sp2str,form) rr
    else
       write(sp2str,'(G32.5)') rr
    end if
    !sp2str = adjustl(sp2str)

  END FUNCTION sp2str


  PURE FUNCTION dp2str(rr,form)
    ! returns real rr as a string (often needed in printing messages)
    IMPLICIT NONE
    REAL(dp),         INTENT(IN)           :: rr
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=32) :: dp2str

    if (present(form)) then
       write(dp2str,form) rr
    else
       write(dp2str,'(G32.5)') rr
    end if
    !dp2str = adjustl(dp2str)

  END FUNCTION dp2str


  PURE FUNCTION log2str(ll,form)
    ! returns logical ll as a string (often needed in printing messages)
    IMPLICIT NONE
    LOGICAL,          INTENT(in)           :: ll
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=10) :: log2str

    if (present(form)) then
       write(log2str,form) ll
    else
       write(log2str,'(L10)') ll
    end if
    !log2str = adjustl(log2str)

  END FUNCTION log2str

  ! ------------------------------------------------------------------

  !     NAME
  !         tolower

  !     PURPOSE
  !         Convert all upper case letters in string to lower case letters.

  !     CALLING SEQUENCE
  !         low = tolower(upper)
  
  !     INDENT(IN)
  !         character(len=*) :: upper    String

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         character(len=len_trim(upper)) :: low    String where all uppercase in input is converted to lowercase

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! Returns 'hallo'
  !         low = tolower('Hallo')
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011 - modified from Echam5, (C) MPI-MET, Hamburg, Germany

  FUNCTION tolower(upper)

    IMPLICIT NONE

    CHARACTER(LEN=*)              ,INTENT(in) :: upper
    CHARACTER(LEN=LEN_TRIM(upper))            :: tolower

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('a')-ICHAR('A')

    DO i=1,LEN_TRIM(upper)
      IF (ICHAR(upper(i:i)) >= ICHAR('A') .AND. &
          ICHAR(upper(i:i)) <= ICHAR('Z')) THEN
        tolower(i:i) = CHAR( ICHAR(upper(i:i)) + idel )
      ELSE
        tolower(i:i) = upper(i:i)
      END IF
    END DO

  END FUNCTION tolower

  ! ------------------------------------------------------------------

  !     NAME
  !         toupper

  !     PURPOSE
  !         Convert all lower case letters in string to upper case letters.

  !     CALLING SEQUENCE
  !         up = toupper(lower)
  
  !     INDENT(IN)
  !         character(len=*) :: lower    String

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         character(len=len_trim(lower)) :: up    String where all lowercase in input is converted to uppercase

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! Returns 'HALLO'
  !         up = toupper('Hallo')
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011 - modified from Echam5, (C) MPI-MET, Hamburg, Germany

  FUNCTION toupper (lower)

    IMPLICIT NONE

    CHARACTER(LEN=*)              ,INTENT(in) :: lower
    CHARACTER(LEN=LEN_TRIM(lower))            :: toupper

    INTEGER            :: i
    INTEGER, PARAMETER :: idel = ICHAR('A')-ICHAR('a')

    DO i=1,LEN_TRIM(lower)
      IF (ICHAR(lower(i:i)) >= ICHAR('a') .AND. &
          ICHAR(lower(i:i)) <= ICHAR('z')) THEN
        toupper(i:i) = CHAR( ICHAR(lower(i:i)) + idel )
      ELSE
        toupper(i:i) = lower(i:i)
      END IF
    END DO

  END FUNCTION toupper

END MODULE mo_string_utils


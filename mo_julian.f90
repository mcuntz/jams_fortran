MODULE mo_julian

  ! This module provides Julian day conversion routines

  ! Written  Matthias Cuntz, Dec 2011

  USE mo_kind, ONLY: i4, sp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: caldat          ! Day , month and year from Julian day
  PUBLIC :: julday          ! Julian day from day, month and year
  PUBLIC :: ndays           ! IMSL Julian day from day, month and year
!   PUBLIC :: ndyin           ! Day, month and year from IMSL Julian day

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         caldat

  !     PURPOSE
  !         Inverse of the function julday. Here julian is input as a Julian Day Number,
  !         and the routine outputs id, mm, and iyyy as the day, month, and year on which the specified 
  !         Julian Day started at noon.

  !     CALLING SEQUENCE
  !         call caldat(Julday, dd, mm, yy)
  
  !     INDENT(IN)
  !         integer(i4) :: Julday     Julian day

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         integer(i4) :: dd         Day in month of Julian day
  !         integer(i4) :: mm         Month in year of Julian day
  !         integer(i4) :: yy         Year of Julian day

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Note: day and month are reversed compared to original Numerical recipes routine.

  !     EXAMPLE
  !         ! 2415021 is 01.01.1900
  !         call caldat(2415021, dd, mm, yy)
  !         -> see also example in test directory

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011 - modified caldat from Numerical Recipes

  SUBROUTINE caldat(julian,id,mm,iyyy)

    use mo_kind, only: i4

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN)  :: julian
    INTEGER(i4), INTENT(OUT) :: id, mm, iyyy

    INTEGER(i4) :: ja, jalpha, jb, jc, jd, je
    INTEGER(i4), PARAMETER :: IGREG = 2299161_i4

    if (julian >= IGREG) then
       jalpha=int(((julian-1867216_i4)-0.25_sp)/36524.25_sp)
       ja=julian+1+jalpha-int(0.25_sp*jalpha)
    else
       ja=julian
    end if
    jb=ja+1524_i4
    jc=int(6680.0_sp+((jb-2439870_i4)-122.1_sp)/365.25_sp)
    jd=365*jc+int(0.25_sp*jc)
    je=int((jb-jd)/30.6001_sp)
    id=jb-jd-int(30.6001_sp*je)
    mm=je-1
    if (mm > 12) mm=mm-12
    iyyy=jc-4715_i4
    if (mm > 2) iyyy=iyyy-1
    if (iyyy <= 0) iyyy=iyyy-1
    
  END SUBROUTINE caldat

  ! ------------------------------------------------------------------

  !     NAME
  !         julday

  !     PURPOSE
  !         In this routine julday returns the Julian Day Number that begins at noon of the calendar
  !         date specified by month mm, day id, and year iyyy, all integer variables. Positive year
  !         signifies A.D.; negative, B.C. Remember that the year after 1 B.C. was 1 A.D.

  !     CALLING SEQUENCE
  !         julian = julday(dd, mm, yy)
  
  !     INDENT(IN)
  !         integer(i4) :: dd         Day in month of Julian day
  !         integer(i4) :: mm         Month in year of Julian day
  !         integer(i4) :: yy         Year of Julian day

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         integer(i4) :: julian     Julian day

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Note: day and month are reversed compared to original Numerical recipes routine.

  !     EXAMPLE
  !         ! 2415021 is 01.01.1900
  !         julian = julday(01,01,1990)
  !         -> see also example in test directory

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011 - modified julday from Numerical Recipes

  FUNCTION julday(id,mm,iyyy)

    use mo_kind,   only: i4, sp
    
    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: id, mm, iyyy
    INTEGER(i4) :: julday
    INTEGER(i4), PARAMETER :: IGREG=15+31*(10+12*1582)
    INTEGER(i4) :: ja,jm,jy

    jy=iyyy
    if (jy == 0) stop 'julday: there is no year zero'
    if (jy < 0) jy=jy+1
    if (mm > 2) then
       jm=mm+1
    else
       jy=jy-1
       jm=mm+13
    end if
    julday=int(365.25_sp*jy)+int(30.6001_sp*jm)+id+1720995_i4
    if (id+31*(mm+12*iyyy) >= IGREG) then
       ja=int(0.01_sp*jy)
       julday=julday+2-ja+int(0.25_sp*ja)
    end if
    
  END FUNCTION julday

  ! ------------------------------------------------------------------

  !     NAME
  !         ndays

  !     PURPOSE
  !         In this routine ndays returns the IMSL Julian Day Number. Julian days begin at noon of the calendar
  !         date specified by month mm, day id, and year iyyy, all integer variables. IMSL treats 01.01.1900
  !         as a reference and assigns a Julian day 0 to it.
  !         ndays = julday(dd,mm,yy) - julday(01,01,1900)

  !     CALLING SEQUENCE
  !         julian = ndays(dd, mm, yy)
  
  !     INDENT(IN)
  !         integer(i4) :: dd         Day in month of IMSL Julian day
  !         integer(i4) :: mm         Month in year of IMSL Julian day
  !         integer(i4) :: yy         Year of IMSL Julian day

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         integer(i4) :: julian     IMSL Julian day, i.e. days before or after 01.01.1900

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! 2415021 is 01.01.1900
  !         julian = ndays(01,01,1990)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011

  FUNCTION ndays(id,mm,iyyy)

    use mo_kind,   only: i4, sp
    
    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: id, mm, iyyy
    INTEGER(i4) :: ndays

    INTEGER(i4), PARAMETER :: IMSL = 2415021_i4

    ndays = julday(id, mm, iyyy) - IMSL
    
  END FUNCTION ndays

  ! ------------------------------------------------------------------

END MODULE mo_julian

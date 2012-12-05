MODULE mo_julian

  ! This module provides Julian day conversion routines
  !
  ! Written  Matthias Cuntz, Dec 2011

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! The UFZ Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The UFZ Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the UFZ Fortran library. If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2011 Matthias Cuntz


  ! Note on Numerical Recipes License
  ! ---------------------------------
  ! Be aware that some code is under the Numerical Recipes License 3rd
  ! edition <http://www.nr.com/aboutNR3license.html>

  ! The Numerical Recipes Personal Single-User License lets you personally
  ! use Numerical Recipes code ("the code") on any number of computers,
  ! but only one computer at a time. You are not permitted to allow anyone
  ! else to access or use the code. You may, under this license, transfer
  ! precompiled, executable applications incorporating the code to other,
  ! unlicensed, persons, providing that (i) the application is
  ! noncommercial (i.e., does not involve the selling or licensing of the
  ! application for a fee), and (ii) the application was first developed,
  ! compiled, and successfully run by you, and (iii) the code is bound
  ! into the application in such a manner that it cannot be accessed as
  ! individual routines and cannot practicably be unbound and used in
  ! other programs. That is, under this license, your application user
  ! must not be able to use Numerical Recipes code as part of a program
  ! library or "mix and match" workbench.

  ! Businesses and organizations that purchase the disk or code download,
  ! and that thus acquire one or more Numerical Recipes Personal
  ! Single-User Licenses, may permanently assign those licenses, in the
  ! number acquired, to individual employees. Such an assignment must be
  ! made before the code is first used and, once made, it is irrevocable
  ! and can not be transferred. 

  ! If you do not hold a Numerical Recipes License, this code is only for
  ! informational and educational purposes but cannot be used.

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: caldat          ! Day , month and year from Julian day
  PUBLIC :: julday          ! Julian day from day, month and year
  PUBLIC :: ndays           ! IMSL Julian day from day, month and year
  PUBLIC :: ndyin           ! Day, month and year from IMSL Julian day

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         caldat

  !     PURPOSE
  !         Inverse of the function julday. Here julian is input as a Julian Day Number,
  !         and the routine outputs id, mm, and yy as the day, month, and year on which the specified 
  !         Julian Day started at noon.

  !         The zeroth Julian Day is 01.01.-4713, i.e. the 1st January 4713 BC.

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

  SUBROUTINE caldat(julian,dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN)  :: julian
    INTEGER(i4), INTENT(OUT) :: dd, mm, yy

    INTEGER(i4) :: ja, jalpha, jb, jc, jd, je
    INTEGER(i4), PARAMETER :: IGREG = 2299161_i4

    if (julian >= IGREG) then
       jalpha = int((real(julian-1867216_i4,dp)-0.25_dp)/36524.25_dp, i4)
       ja     = julian + 1 + jalpha - int(0.25_dp*real(jalpha,dp), i4)
    else
       ja     = julian
    end if
    jb = ja + 1524_i4
    jc = int(6680.0_dp+(real(jb-2439870_i4,dp)-122.1_dp)/365.25_dp, i4)
    jd = 365*jc + int(0.25_dp*real(jc,dp), i4)
    je = int(real(jb-jd,dp)/30.6001_dp, i4)
    dd = jb - jd - int(30.6001_dp*real(je,dp), i4)
    mm = je - 1
    if (mm > 12) mm = mm - 12
    yy = jc - 4715_i4
    if (mm > 2)  yy = yy - 1
    if (yy <= 0) yy = yy - 1
    
  END SUBROUTINE caldat

  ! ------------------------------------------------------------------

  !     NAME
  !         julday

  !     PURPOSE
  !         In this routine julday returns the Julian Day Number that begins at noon of the calendar
  !         date specified by month mm, day dd, and year yy, all integer variables. Positive year
  !         signifies A.D.; negative, B.C. Remember that the year after 1 B.C. was 1 A.D.

  !         The zeroth Julian Day is 01.01.-4713, i.e. the 1st January 4713 BC.

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

  FUNCTION julday(dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: dd, mm, yy
    INTEGER(i4) :: julday
    INTEGER(i4), PARAMETER :: IGREG=15+31*(10+12*1582)
    INTEGER(i4) :: ja,jm,jy

    jy=yy
    if (jy == 0) stop 'julday: there is no year zero'
    if (jy < 0) jy=jy+1
    if (mm > 2) then
       jm=mm+1
    else
       jy=jy-1
       jm=mm+13
    end if

    !Bug julday=int(365.25_dp*jy)+int(30.6001_dp*jm)+dd+1720995_i4
    julday = 365*jy + int(0.25_dp*real(jy,dp)+2000._dp, i4) + int(30.6001_dp*real(jm,dp), i4) + dd + 1718995_i4
    if (dd+31*(mm+12*yy) >= IGREG) then
       ja = int(0.01_dp*jy, i4)
       julday=julday + 2 - ja + int(0.25_dp*ja, i4)
    end if
    
  END FUNCTION julday

  ! ------------------------------------------------------------------

  !     NAME
  !         ndays

  !     PURPOSE
  !         In this routine ndays returns the IMSL Julian Day Number. Julian days begin at noon of the calendar
  !         date specified by month mm, day dd, and year yy, all integer variables. IMSL treats 01.01.1900
  !         as a reference and assigns a Julian day 0 to it.
  !             ndays = julday(dd,mm,yy) - julday(01,01,1900)

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
  !         ! 0 is 01.01.1900
  !         julian = ndays(01,01,1990)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011

  FUNCTION ndays(dd,mm,yy)
    
    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: dd, mm, yy
    INTEGER(i4) :: ndays

    INTEGER(i4), PARAMETER :: IMSLday = 2415021_i4

    ndays = julday(dd, mm, yy) - IMSLday
    
  END FUNCTION ndays

  ! ------------------------------------------------------------------

  !     NAME
  !         ndyin

  !     PURPOSE
  !         Inverse of the function ndys. Here ISML Julian is input as a Julian Day Number
  !         minus the Julian Day Number of 01.01.1900, and the routine outputs id, mm, and yy
  !         as the day, month, and year on which the specified Julian Day started at noon.
  !           ndyin is caldat(IMSLJulian + 2415021, dd, mm, yy)

  !     CALLING SEQUENCE
  !         call ndyin(julian, dd, mm, yy)
  
  !     INDENT(IN)
  !         integer(i4) :: julian     IMSL Julian day, i.e. days before or after 01.01.1900

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         integer(i4) :: dd         Day in month of IMSL Julian day
  !         integer(i4) :: mm         Month in year of IMSL Julian day
  !         integer(i4) :: yy         Year of IMSL Julian day

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! 0 is 01.01.1900
  !         call ndyin(0,dd,mm,yy)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Dec 2011

  SUBROUTINE ndyin(julian,dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN)  :: julian
    INTEGER(i4), INTENT(OUT) :: dd, mm, yy

    INTEGER(i4), PARAMETER :: IMSLday = 2415021_i4

    call caldat(julian+IMSLday, dd, mm, yy)
    
  END SUBROUTINE ndyin

  ! ------------------------------------------------------------------

END MODULE mo_julian

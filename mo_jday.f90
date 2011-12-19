!
!**************************************************************************
!
! PROGRAM    :: Julian_Day_Conversion
! PURPOSE    :: Convert a given julian day into day, month and year 
!               and Vice-versa 
! Created by :: R. Kumar (UFZ) 
!               Sources from where basic *.f90 files were taken are 
!               appropriately acknowledged. 
!
!
! NOTE ::=================================================================
!  module mo_kind
!  implicit none
!  ! Integer and floating point section
!  ! see:   http://fortranwiki.org/fortran/show/Real+precision
!  integer, parameter                        :: i4 = selected_int_kind(9)
!  integer, parameter                        :: sp = selected_real_kind(6,37)
!  integer, parameter                        :: dp = selected_real_kind(15,307)
!  !
!end module mo_kind
!=============================================================================
!
!**************************************************************************
MODULE mo_jday
  !
  use mo_kind, only: i4
  !
  implicit none
  !
  private
  !
  public :: NDYIN ! calculates day, month and year of a given julian day
  public :: NDAYS ! calculates julian day of a given day, month and year
  !
 CONTAINS
   !
   !**************************************************************************
   !
   !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      
   !
   !**************************************************************************
   !
   ! ** ROUTINE: ODS_CalDat() --- Convert Julian day to calendar date
   ! 
   ! ** DESCRIPTION:
   !
   !     The routine returns the calendar day based on the Julian
   !     day number. The algorithm was adopted from Press et al.
   !
   !     The Julian day number on May 23, 1968 is 2,440,000. 
   !
   !     The return values are day, month, & year of a julian day 
   !
   !     A negative number implies that the year is B.C.
   !
   !     Refrence: Press, William H., Saul A. Teukolsky, William T.
   !          Vetterling and Brian P. Flannery, 1992: Numerical
   !          Recipes in Fortran, 2nd Ed. Cambridge University
   !          Press, New York, NY, 963pp.
   !-------------------------------------------------------------------------
   ! *** CHANGE MADE BY R. KUMAR ---> TO MAKE IT COMPATABLE WITH IMSL LIBRARY
   !     IMSL treats 1.1.1900 as a base and assigins for this day, the julian
   !     day value as 0.
   !
   !     The present algorithm is based on the Gregorian Calendar
   !     adopted Oct 12, 1582
   !
   !     According to this algorithm the julian day of 1.1.1900 is 2415021
   !
   !     The number 2415021 will be ADDED to the Julian day in the present  
   !     algorithm so that it is compatable with IMSL produced day month and
   !      year for a given julian day
   !-------------------------------------------------------------------------
   subroutine NDYIN(JulDay, day1, month1, year1)
     use mo_kind
     implicit none
     !
     ! INPUT & OUTPUT PARAMETERS:
     integer(i4), intent(in)          :: JulDay                      ! Julian day number
     integer(i4), intent(out)         :: Day1, Month1, Year1            ! day, month & year
     !
     ! Other variables
     integer(i4)                      :: Day, Month, Year            ! day, month & year
     integer(i4), parameter           :: REF_to_IMSL_jday = 2415021  !>> SEE ABOVE FOR CLARIFICATION
     integer(i4)                      :: JulDay_IMSL_REF_FREE        ! Julian day FREE from IMSL refrence
     !
     integer(i4), parameter           :: iGreg = 2299161             ! The Julian day number of the Gregorian Calendar adopted Oct 12, 1582
     integer(i4)                      :: Alpha
     integer(i4)                      :: ja, jb, jc, jd, je
     !
     ! 
     ! Convert the imsl refrenced julian day to the GREGORIAN CALENDAR Julian day
     ! on which this algorithm is based
     JulDay_IMSL_REF_FREE = JulDay + REF_to_IMSL_jday
     !
     !Cross-over to Gregorian Calendar produces this correction
     if ( JulDay_IMSL_REF_FREE >= iGreg ) then
       Alpha = int((( real(JulDay_IMSL_REF_FREE, sp) - 1867216.0_sp ) - 0.25_sp ) / 36524.25_sp )
       ja    = JulDay_IMSL_REF_FREE + 1 + Alpha - int( 0.25_sp * real(Alpha, sp) )
     else 
       ! no correction
       ja    = JulDay
     end if
     !
     jb    = ja + 1524
     jc    = int( 6680.0_sp + ((real(jb, sp) - 2439870) - 122.1_sp ) / 365.25_sp )
     jd    = 365 * jc + int( 0.25_sp * real(jc, sp) )
     je    = int( ((jb - jd) / 30.6001_sp ), i4)
     Day   = jb - jd - int( 30.6001_sp * real(je, sp) )
     Month = je - 1
     !
     if( Month > 12 ) Month = Month - 12
     Year  = jc - 4715
     !
     if( Month  > 2 ) Year = Year - 1
     if( Year  <= 0 ) Year = Year - 1
     !
     ! OUTPUT
     day1   = day
     month1 = month
     year1  = year
     !
   end subroutine NDYIN
   !
   !
   !**************************************************************************
   !
   !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      
   !
   !**************************************************************************
   !
   ! ** ROUTINE: ODS_Julian() --- Returns the Julian day
   ! 
   ! ** DESCRIPTION:
   !
   !     The routine returns the Julian day number based on the
   !     calendar date. The algorithm was adopted from Press et al.
   !
   !     The Julian day number on May 23, 1968 is 2,440,000.  
   !
   !     The year zero is treated as the year 1 since the year 0 does
   !     not exist.
   !
   !     Refrence: Press, William H., Saul A. Teukolsky, William T.
   !          Vetterling and Brian P. Flannery, 1992: Numerical
   !          Recipes in Fortran, 2nd Ed. Cambridge University
   !          Press, New York, NY, 963pp.
   !-------------------------------------------------------------------------
   ! *** CHANGE MADE BY R. KUMAR ---> TO MAKE IT COMPATABLE WITH IMSL LIBRARY
   !     IMSL treats 1.1.1900 as a base and assigins for this day, the julian
   !     day value as 0.
   !
   !     The present algorithm is based on the Gregorian Calendar
   !     adopted Oct 12, 1582
   !
   !     According to this algorithm the julian day of 1.1.1900 is 2415021
   !
   !     The number 2415021 will be deducted from the Julian day estimated  
   !     with the present algorithm so that it is compatable with IMSL 
   !     calculated julian day. 
   !-------------------------------------------------------------------------
   !
   integer(i4) function NDAYS(day, month, year)
     use mo_kind
     implicit NONE
     !
     ! INPUT PARAMETERS:
     integer(i4)            :: Day, Month, Year
     !
     integer(i4), parameter :: REF_to_IMSL_jday = 2415021 !>> SEE ABOVE FOR CLARIFICATION
     integer(i4)            :: iGreg  ! Gregorian Calendar adopted Oct 12, 1582
     integer(i4)            :: JulDay
     integer(i4)            :: jy, jm, ja
     !
     iGreg = 15 + 31 * (10 + 12 * 1582) 
     !
     ! Change year 0 to year 1
     if(Year == 0) Year = 1
     !
     ! Account for the nonexisting year 0
     if(Year  < 0 ) Year = Year + 1
     if(Month > 2 ) then
       jy = Year
       jm = Month + 1
     else
       jy = Year  - 1
       jm = Month + 13
     endif
     !
     JulDay =  int ( 365.25_sp  * real(jy, sp) ) &
             + int ( 30.6001_sp * real(jm, sp) ) &
             + Day + 1720995
     ! Test whether to change to Gregorian Celendar
     if( Day + 31 * ( Month + 12 * Year ) >= iGreg) then
       ja     = int( 0.01_sp * real(jy, sp) )
       Julday = JulDay + 2 - ja + int( 0.25_sp * real(ja, sp) )
     endif
     !
     ! Refrence JDAY based on the julian day of 01.01.1900 == 0
     JulDay = JulDay - REF_to_IMSL_jday   
     !
     NDAYS = JulDay
     !
   end function NDAYS
   !
 END MODULE Mo_jday

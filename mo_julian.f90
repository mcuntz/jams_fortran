!> \file mo_julian.f90

!> \brief Julian date conversion routines

!> \details Julian date to and from day, month, year, and also from day, month, year, hour, minute, and second.\n
!> Different calendars provided:\n
!>     gregorian/standard, lilian, 360day, 365day, 366day, 360_day, 365_day/noleap, 366_day/all_leap.\n
!> There are differences between 365day and 365_day calendars, fo example; see below.
!> There is no julian calendar.\n
!> Convenience routines for Julian dates of IMSL are provided (start at 01.01.1990).
!> Also relative dates can be given in dec2date with a unit indicating the reference date.

!> \note Julian day definition starts at noon of the 1st January 4713 BC.\n
!> Here, the astronomical definition is used,
!> i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
!> This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n
!> \n
!> julday and caldat start at midnight of the 1st January 4713 BC.
!> So date2dec and julday as well as dec2date and caldat are shifted by half a day.\n
!> Use date2dec with dec2date together for fractional Julian dates
!> and use julday with caldat together for integer Julian days.

!> \note Lilian dates start at midnight of 15.10.1582, when Pope Gregory XIII introduced the Gregorian calendar.

!> \note 360day, 365day and 366day calendars start at 01.01.0000 00:00:00,\n
!> while 360_day, 365_day and 366_day calendars are Gregorian calendars starting at 01.01.-4712 noon,
!> as given by the CF-convention. 365_day and 366_day calendars are also called noleap and all_leap.

!> \note Units in dec2date can only be "days/minutes/hours/seconds since YYYY-MM-DD hh:mm:ss".
!> Any precision after reference time giving the time zone (Z, +hh:mm) will be ignored.
!> Numbers in reference date do not have to have leading zeros.

!> \author Matthias Cuntz
!> \date Dec 2011

MODULE mo_julian

  ! This module provides Julian day conversion routines

  ! Written  Matthias Cuntz, Dec 2011
  ! Modified Matthias Cuntz, Jan 2013 - added date2dec and dec2date
  ! Modified Matthias Cuntz, May 2014 - changed to new algorithm with astronomical units
  !                                     removed numerical recipes
  ! Modified David Schaefer, Oct 2015 - addded 360 day calendar procedures
  ! Modified David Schaefer, Jan 2016 - addded 365 day calendar procedures
  ! Modified David Schaefer, Feb 2016 - implemented wrapper function and the module calendar state
  ! Modified Matthias Cuntz, Dec 2016 - remove save variable calendar and associated routines setCalendar
  !                                     pass calendar to conversion routines
  ! Modified Matthias Cuntz, Dec 2016 - Lilian date, fracday
  ! Modified Matthias Cuntz, Mar 2018 - units in dec2date
  ! Modified Matthias Cuntz, Oct 2018 - calendars 360_day, 365_day/noleap, 366_day/all_leap
  !                                   - complete all missing routines of 36*day
  ! Modified Matthias Cuntz, Apr 2019 - allow non-leading zero notation in units

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2011-2019 Matthias Cuntz, David Schaefer - mc (at) macu (dot) de
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

  USE mo_kind, ONLY: i4, i8, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: caldat
  PUBLIC :: date2dec     ! Fractional Julian day from day, month, year, hour, minute, and second
  PUBLIC :: dec2date     ! Day, month, year, hour, minute, and second from fractional Julian day
  PUBLIC :: julday       ! Julian day from day, month and year
  PUBLIC :: ndays        ! IMSL Julian day from day, month and year
  PUBLIC :: ndyin        ! Day, month and year from IMSL Julian day

  public :: caldat360day, caldat365day, caldat366day, caldat360_day, &
       caldat365_day, caldat366_day, caldatJulian, caldatLilian
  public :: date2dec360day, date2dec365day, date2dec366day, date2dec360_day, &
       date2dec365_day, date2dec366_day, date2decJulian, date2decLilian
  public :: dec2date360day, dec2date365day, dec2date366day, dec2date360_day, &
       dec2date365_day, dec2date366_day, dec2dateJulian, dec2dateLilian
  public :: julday360day, julday365day, julday366day, julday360_day, &
       julday365_day, julday366_day, juldayJulian, juldayLilian

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         caldat

  !     PURPOSE
  !>        \brief Day, month and year from Julian day in the current or given calendar

  !>        \details Wrapper around the calendar specific caldat procedures.
  !>        Inverse of the function julday. Here julian is input as a Julian Day Number,
  !>        and the routine outputs d0d, mm, and yy as the day, month, and year on which the specified
  !>        Julian Day started at noon.

  !>        The zeroth Julian Day depends on the called procedure. See their documentation for details.

  !     CALLING SEQUENCE
  !         call caldat(julday, dd, mm, yy, calendar)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: julday"     Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: calendar"    The calendar to use. Available calendars:\n
  !>                    gregorian/standard, lilian, 360day, 365day, 366day, 360_day, 365_day/noleap, 366_day/all_leap\n
  !>                    All other strings are taken as 'gregorian'.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, David Schaefer
  !>        \date Jan 2015
  !>        Modified, Matthias Cuntz, Oct 2018 - added more calendars
  elemental subroutine caldat(julian, dd, mm, yy, calendar)

    implicit none

    integer(i4),  intent(in)           :: julian
    integer(i4),  intent(out)          :: dd, mm, yy
    character(*), intent(in), optional :: calendar

    if (present(calendar)) then
       select case(calendar)
       case("360day")
          call caldat360day(julian,dd,mm,yy)
       case("365day")
          call caldat365day(julian,dd,mm,yy)
       case("366day")
          call caldat366day(julian,dd,mm,yy)
       case("360_day")
          call caldat360_day(julian,dd,mm,yy)
       case("365_day", "noleap")
          call caldat365_day(julian,dd,mm,yy)
       case("366_day", "all_leap")
          call caldat366_day(julian,dd,mm,yy)
       case("lilian")
          call caldatLilian(julian,dd,mm,yy)
       case default
          call caldatJulian(julian,dd,mm,yy)
       end select
    else
       call caldatJulian(julian,dd,mm,yy)
    endif

  end subroutine caldat


  ! ------------------------------------------------------------------

  !     NAME
  !         date2dec

  !     PURPOSE
  !>        \brief Fractional Julian day from day, month, year, hour, minute, second in the current calendar

  !>        \details Wrapper around the calendar specific date2dec procedures.
  !>        In this routine date2dec returns the fractional Julian Day that begins at noon
  !>        of the calendar date specified by month mm, day dd, and year yy, all integer variables.

  !>        The zeroth Julian Day depends on the called procedure. See their documentation for details.

  !     CALLING SEQUENCE
  !         date2dec = date2dec(dd, mm, yy, hh, nn, ss, fracday, calendar)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: dd"         Day in month of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: mm"         Month in year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: yy"         Year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: hh"         Hours of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: nn"         Minutes of hour of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: ss"         Secondes of minute of hour of Julian day (default: 0)
  !>        \param[in] "real(dp),    optional :: fracday"    Fractional day, only used if hh,nn,ss not given (default: 0)
  !>        \param[in] "integer(i4), optional :: calendar"   The calendar to use. Available calendars:\n
  !>                    gregorian/standard, lilian, 360day, 365day, 366day, 360_day, 365_day/noleap, 366_day/all_leap\n
  !>                    All other strings are taken as 'julian'.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: date2dec  !     Fractional Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, David Schaefer
  !>        \date Jan 2015
  !>        Modified, Matthias Cuntz, Oct 2018 - added more calendars
  elemental function date2dec(dd, mm, yy, hh, nn, ss, fracday, calendar)

    implicit none

    integer(i4),  intent(in), optional :: dd, mm, yy
    integer(i4),  intent(in), optional :: hh, nn, ss
    real(dp),     intent(in), optional :: fracday
    character(*), intent(in), optional :: calendar
    real(dp)                           :: date2dec

    if (present(calendar)) then
       select case(calendar)
       case("360day")
          date2dec = date2dec360day(dd, mm, yy, hh, nn, ss, fracday)
       case("365day")
          date2dec = date2dec365day(dd, mm, yy, hh, nn, ss, fracday)
       case("366day")
          date2dec = date2dec366day(dd, mm, yy, hh, nn, ss, fracday)
       case("360_day")
          date2dec = date2dec360_day(dd, mm, yy, hh, nn, ss, fracday)
       case("365_day", "noleap")
          date2dec = date2dec365_day(dd, mm, yy, hh, nn, ss, fracday)
       case("366_day", "all_leap")
          date2dec = date2dec366_day(dd, mm, yy, hh, nn, ss, fracday)
       case("lilian")
          date2dec = date2decLilian(dd, mm, yy, hh, nn, ss, fracday)
       case default
          date2dec = date2decJulian(dd, mm, yy, hh, nn, ss, fracday)
       end select
    else
       date2dec = date2decJulian(dd, mm, yy, hh, nn, ss, fracday)
    endif

  end function date2dec


  ! ------------------------------------------------------------------

  !     NAME
  !         dec2date

  !     PURPOSE
  !>        \brief Day, month, year, hour, minute, and second from fractional Julian day in a given calendar.

  !>        \details Wrapper around the calendar specific dec2date procedures.
  !>        Inverse of the function date2dec. Here dec2date is input as a fractional Julian Day.
  !>        The routine outputs dd, mm, yy, hh, nn, ss as the day, month, year, hour, minute, and second
  !>        on which the specified Julian Day started at noon.\n
  !>        Fractional day can be output optionally.

  !>        The zeroth Julian Day depends on the called procedure. See their documentation for details.

  !     CALLING SEQUENCE
  !         call dec2date(fJulian, dd, mm, yy, hh, nn, ss, fracday, calendar, units)

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: fJulian"     fractional Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "character(len=*), optional :: calendar"    The calendar to use. Available calendars:\n
  !>                    gregorian/standard, lilian, 360day, 365day, 366day, 360_day, 365_day/noleap, 366_day/all_leap\n
  !>                    All other strings are taken as 'julian'.
  !>        \param[in] "character(len=*), optional :: units"       Units of decimal date in ISO 8601 format.\n
  !>                   Can only be "days/minutes/hours/seconds since YYYY-MM-DD hh:mm:ss".\n
  !>                   Any precision after reference time giving the time zone
  !>                   (http://www.cl.cam.ac.uk/%7Emgk25/iso-time.html -> Z, +hh:mm) will be ignored, e.g.
  !>                   YYYY-MM-DD hh:mm:ssZ for UTC or YYYY-MM-DD hh:mm:ss+hh:mm for a specific time zone.
  !>
  !>                   Allow exception to ISO 8601 that numbers must not have leading zeros, i.e.
  !>                   Y-M-D h:m, which is for example produced by the climate data operators (cdo).

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "integer(i4), optional :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4), optional :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4), optional :: yy"         Year of Julian day
  !>        \param[out] "integer(i4), optional :: hh"         Hour of Julian day
  !>        \param[out] "integer(i4), optional :: nn"         Minute in hour of Julian day
  !>        \param[out] "integer(i4), optional :: ss"         Second in minute of hour of Julian day
  !>        \param[out] "real(dp),    optional :: fracday"    Fractional day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, David Schaefer
  !>        \date Jan 2015
  !>        Modified, Matthias Cuntz, Mar 2018 - units
  !>        Modified, Matthias Cuntz, Oct 2018 - added more calendars
  !>        Modified, Matthias Cuntz, Apr 2019 - allow non-leading zero notation in units
  elemental subroutine dec2date(julian, dd, mm, yy, hh, nn, ss, fracday, calendar, units)

    implicit none

    real(dp),         intent(in)            :: julian
    integer(i4),      intent(out), optional :: dd, mm, yy, hh, nn, ss
    real(dp),         intent(out), optional :: fracday
    character(len=*), intent(in),  optional :: calendar
    character(len=*), intent(in),  optional :: units

    character(64) :: icalendar, iunit, idate0
    integer(i4)   :: year0, month0, day0, hour0, minute0, second0
    real(dp)      :: jdate0, jdate, eps

    icalendar = 'julian'
    if (present(calendar)) icalendar = calendar

    jdate = julian
    ! No write in pure routines so no error handling (commented)
    if (present(units)) then
       ! get reference date
       hour0   = 0
       minute0 = 0
       second0 = 0
       iunit  = units(1:index(trim(units),'since')-1) ! days/minutes/hours/seconds
       idate0 = units(index(trim(units),'since')+6:)  ! YYYY-MM-DD hh:mm:ss
       read(idate0(1:index(trim(idate0),'-')-1),*) year0
       idate0 = idate0(index(trim(idate0),'-')+1:) ! MM-DD hh:mm:ss
       read(idate0(1:index(trim(idate0),'-')-1),*) month0
       idate0 = idate0(index(trim(idate0),'-')+1:) ! DD hh:mm:ss
       if (index(trim(idate0),' ') > 0) then
          read(idate0(1:index(trim(idate0),' ')-1),*) day0
          idate0 = idate0(index(trim(idate0),' ')+1:) ! hh:mm:ss
          if (index(trim(idate0),':') > 0) then
             read(idate0(1:index(trim(idate0),':')-1),*) hour0
             idate0 = idate0(index(trim(idate0),':')+1:) ! mm:ss
             if (index(trim(idate0),':') > 0) then
                read(idate0(1:index(trim(idate0),':')-1),*) minute0
                idate0 = idate0(index(trim(idate0),':')+1:) ! ss
                if (len_trim(idate0) > 0) read(idate0(1:len_trim(idate0)),*) second0
             else
                read(idate0(1:len_trim(idate0)),*) minute0
             endif
          else
             read(idate0(1:len_trim(idate0)),*) hour0
          endif
       else
          read(idate0(1:len_trim(idate0)),*) day0
       endif
       jdate0 = date2dec(day0, month0, year0, hour0, minute0, second0, calendar=icalendar)
       ! Julian date should have a small offset proportional to julian date for correct re-conversion.
       ! Add it to final Julian date so substract it from jdate0
       eps = epsilon(1.0_dp)
       if (eps * abs(jdate0) > eps) then
          if (jdate0 > 0._dp) then
             jdate0 = jdate0 / (1._dp+eps)
          else
             jdate0 = jdate0 / (1._dp-eps)
          endif
       else
          jdate0 = jdate0 - eps
       endif
       ! Add time to reference date
       select case(trim(iunit))
       case("days")
          jdate = jdate0 + julian
       case("hours")
          jdate = jdate0 + julian/24._dp
       case("minutes")
          jdate = jdate0 + julian/1440._dp
       case("seconds")
          jdate = jdate0 + julian/86400._dp
       ! case default
       !    write(*,*) 'Error dec2date: time unit not days, hours, minutes, or seconds.'
       !    write(*,*) '    Units must be "days/minutes/hours/seconds since YYYY-MM-DD hh:mm:ss".'
       !    write(*,*) '    hh:mm:ss can be omitted. Units given: ', trim(units)
       !    stop
       end select
       ! Now add the offset for the numerics
       eps = max(eps * abs(jdate), eps)
       jdate = jdate + eps
    endif

    ! dec2date
    select case(icalendar)
    case("360day")
       call dec2date360day(jdate, dd, mm, yy, hh, nn, ss, fracday)
    case("365day")
       call dec2date365day(jdate, dd, mm, yy, hh, nn, ss, fracday)
    case("366day")
       call dec2date366day(jdate, dd, mm, yy, hh, nn, ss, fracday)
    case("360_day")
       call dec2date360_day(jdate, dd, mm, yy, hh, nn, ss, fracday)
    case("365_day", "noleap")
       call dec2date365_day(jdate, dd, mm, yy, hh, nn, ss, fracday)
    case("366_day", "all_leap")
       call dec2date366_day(jdate, dd, mm, yy, hh, nn, ss, fracday)
    case("lilian")
       call dec2dateLilian(jdate, dd, mm, yy, hh, nn, ss, fracday)
    case default
       call dec2dateJulian(jdate, dd, mm, yy, hh, nn, ss, fracday)
    end select

  end subroutine dec2date


  ! ------------------------------------------------------------------

  !     NAME
  !         julday

  !     PURPOSE
  !>        \brief Julian day from day, month and year in the current or given calendar

  !>        \details Wrapper around the calendar specific julday procedures.
  !>        In this routine julday returns the Julian Day Number that begins at noon of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables.

  !>        The zeroth Julian Day depends on the called procedure. See their documentation for details.

  !     CALLING SEQUENCE
  !         julian = julday(dd, mm, yy, calendar)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[in] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[in] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: calendar"    The calendar to use. Available calendars:\n
  !>                    gregorian/standard, lilian, 360day, 365day, 366day, 360_day, 365_day/noleap, 366_day/all_leap\n
  !>                    All other strings are taken as 'julian'.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: julian  !     Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, David Schaefer
  !>        \date Jan 2015
  !>        Modified, Matthias Cuntz, Oct 2018 - added more calendars
  elemental function julday(dd, mm, yy, calendar)

    implicit none

    integer(i4),  intent(in)           :: dd, mm, yy
    character(*), intent(in), optional :: calendar
    integer(i4)                        :: julday

    if (present(calendar)) then
       select case(calendar)
       case("360day")
          julday = julday360day(dd, mm, yy)
       case("365day")
          julday = julday365day(dd, mm, yy)
       case("366day")
          julday = julday366day(dd, mm, yy)
       case("360_day")
          julday = julday360_day(dd, mm, yy)
       case("365_day", "noleap")
          julday = julday365_day(dd, mm, yy)
       case("366_day", "all_leap")
          julday = julday366_day(dd, mm, yy)
       case("lilian")
          julday = juldayLilian(dd, mm, yy)
       case default
          julday = juldayJulian(dd, mm, yy)
       end select
    else
       julday = juldayJulian(dd, mm, yy)
    endif

  end function julday


  ! ------------------------------------------------------------------

  !     NAME
  !         ndays

  !     PURPOSE
  !>        \brief IMSL Julian day from day, month and year

  !>        \details In this routine ndays returns the IMSL Julian Day Number. Julian days begin at noon of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables. IMSL treats 01.01.1900
  !>        as a reference and assigns a Julian day 0 to it.
  !>            ndays = julday(dd,mm,yy) - julday(01,01,1900)

  !     CALLING SEQUENCE
  !         julian = ndays(dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: dd"         Day in month of IMSL Julian day
  !>        \param[in] "integer(i4) :: mm"         Month in year of IMSL Julian day
  !>        \param[in] "integer(i4) :: yy"         Year of IMSL Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: julian &mdash;     IMSL Julian day, i.e. days before or after 01.01.1900

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! 0 is 01.01.1900
  !         julian = ndays(01,01,1990)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Dec 2011

  ELEMENTAL FUNCTION ndays(dd,mm,yy)

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
  !>        \brief Day, month and year from IMSL Julian day

  !>        \details Inverse of the function ndys. Here ISML Julian is input as a Julian Day Number
  !>        minus the Julian Day Number of 01.01.1900, and the routine outputs id, mm, and yy
  !>        as the day, month, and year on which the specified Julian Day started at noon.
  !>          ndyin is caldat(IMSLJulian + 2415021, dd, mm, yy)

  !     CALLING SEQUENCE
  !         call ndyin(julian, dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: julian"     IMSL Julian day, i.e. days before or after 01.01.1900

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4) :: dd"         Day in month of IMSL Julian day
  !>        \param[out] "integer(i4) :: mm"         Month in year of IMSL Julian day
  !>        \param[out] "integer(i4) :: yy"         Year of IMSL Julian day

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
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
  !>        \author Written, Matthias Cuntz
  !>        \date Dec 2011

  ELEMENTAL SUBROUTINE ndyin(julian,dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN)  :: julian
    INTEGER(i4), INTENT(OUT) :: dd, mm, yy

    INTEGER(i4), PARAMETER :: IMSLday = 2415021_i4

    call caldat(julian+IMSLday, dd, mm, yy)

  END SUBROUTINE ndyin

  
  ! ------------------------------------------------------------------
  ! PRIVATE
  ! ------------------------------------------------------------------

  
  ! ------------------------------------------------------------------
  ! CALDAT
  ! ------------------------------------------------------------------


  ! ------------------------------------------------------------------

  !     NAME
  !         caldat360day

  !     PURPOSE
  !>        \brief Day, month and year from Julian day in a 360 day calendar.

  !>        \details Inverse of the function julday360day. Here julian is input as a Julian Day Number,
  !>        and the routine outputs dd, mm, and yy as the day, month, and year on which the specified
  !>        Julian Day started at noon.

  !>        The zeroth Julian Day here is 01.01.0000

  !     CALLING SEQUENCE
  !         call caldat360day(julday, dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: julday"     Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, David Schaefer
  !>        \date Oct 2015
  elemental subroutine caldat360day(julian,dd,mm,yy)

    implicit none

    integer(i4), intent(in)  :: julian
    integer(i4), intent(out) :: dd, mm, yy
    integer(i4), parameter   :: year=360, month=30
    integer(i4)              :: remainder

    yy = julian/year
    remainder = mod(abs(julian),year)
    mm = remainder/month + 1
    dd = mod(abs(julian), month) + 1

  end subroutine caldat360day


  ! ------------------------------------------------------------------

  !     NAME
  !         caldat365day

  !     PURPOSE
  !>        \brief Day, month and year from Julian day in a 365 day calendar.

  !>        \details Inverse of the function julday365day. Here julian is input as a Julian Day Number,
  !>        and the routine outputs dd, mm, and yy as the day, month, and year on which the specified
  !>        Julian Day started at noon.

  !>        The zeroth Julian Day here is 01.01.0000

  !     CALLING SEQUENCE
  !         call caldat365day(julday, dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: julday"     Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, David Schaefer
  !>        \date Dec 2015
  elemental subroutine caldat365day(julian,dd,mm,yy)

    implicit none

    integer(i4), intent(in)  :: julian
    integer(i4), intent(out) :: dd, mm, yy
    integer(i4), parameter   :: year=365
    integer(i4), dimension(12), parameter :: months=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    integer(i4)              :: remainder

    yy = julian/year
    remainder = mod(abs(julian),year) + 1

    do mm = 1,size(months)
       if (remainder .le. months(mm)) then
          exit
       end if
       remainder = remainder - months(mm)
    end do

    dd = remainder

  end subroutine caldat365day


  ! ------------------------------------------------------------------

  !     NAME
  !         caldat366day

  !     PURPOSE
  !>        \brief Day, month and year from Julian day in a 366 day calendar.

  !>        \details Inverse of the function julday366day. Here julian is input as a Julian Day Number,
  !>        and the routine outputs dd, mm, and yy as the day, month, and year on which the specified
  !>        Julian Day started at noon.

  !>        The zeroth Julian Day here is 01.01.0000

  !     CALLING SEQUENCE
  !         call caldat366day(julday, dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: julday"     Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental subroutine caldat366day(julian,dd,mm,yy)

    implicit none

    integer(i4), intent(in)  :: julian
    integer(i4), intent(out) :: dd, mm, yy
    integer(i4), parameter   :: year=366
    integer(i4), dimension(12), parameter :: months=(/31,29,31,30,31,30,31,31,30,31,30,31/)
    integer(i4)              :: remainder

    yy = julian/year
    remainder = mod(abs(julian),year) + 1

    do mm = 1,size(months)
       if (remainder .le. months(mm)) then
          exit
       end if
       remainder = remainder - months(mm)
    end do

    dd = remainder

  end subroutine caldat366day


  ! ------------------------------------------------------------------

  !     NAME
  !         caldat360_day

  !     PURPOSE
  !>        \brief Day, month and year from Julian day in a Gregorian calendar with only 30-day months.

  !>        \details Inverse of the function julday360_day. Here julian is input as a Julian Day Number,
  !>        and the routine outputs dd, mm, and yy as the day, month, and year on which the specified
  !>        Julian Day started at noon.

  !>        The zeroth Julian Day here is 01.01.0000

  !     CALLING SEQUENCE
  !         call caldat360_day(julday, dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: julday"     Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental subroutine caldat360_day(julian,dd,mm,yy)

    implicit none

    integer(i4), intent(in)  :: julian
    integer(i4), intent(out) :: dd, mm, yy
    integer(i4), parameter   :: year=360, month=30
    integer(i4)              :: remainder

    yy = julian/year - 4712
    remainder = mod(abs(julian),year)
    mm = remainder/month + 1
    dd = mod(abs(julian), month) + 1

  end subroutine caldat360_day


  ! ------------------------------------------------------------------

  !     NAME
  !         caldat365_day

  !     PURPOSE
  !>        \brief Day, month and year from Julian day in a Gregorian calendar without leap years.

  !>        \details Inverse of the function julday365_day. Here julian is input as a Julian Day Number,
  !>        and the routine outputs dd, mm, and yy as the day, month, and year on which the specified
  !>        Julian Day started at noon.

  !>        The zeroth Julian Day here is 01.01.-4712.

  !     CALLING SEQUENCE
  !         call caldat365_day(julday, dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: julday"     Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental subroutine caldat365_day(julian,dd,mm,yy)

    implicit none

    integer(i4), intent(in)  :: julian
    integer(i4), intent(out) :: dd, mm, yy
    integer(i4), parameter   :: year=365
    integer(i4), dimension(12), parameter   :: months=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    integer(i4)              :: remainder

    yy = julian/year - 4712
    remainder = mod(abs(julian),year) + 1

    do mm = 1,size(months)
       if (remainder .le. months(mm)) then
          exit
       end if
       remainder = remainder - months(mm)
    end do

    dd = remainder

  end subroutine caldat365_day


  ! ------------------------------------------------------------------

  !     NAME
  !         caldat366_day

  !     PURPOSE
  !>        \brief Day, month and year from Julian day in a Gregorian calendar without leap years.

  !>        \details Inverse of the function julday366_day. Here julian is input as a Julian Day Number,
  !>        and the routine outputs dd, mm, and yy as the day, month, and year on which the specified
  !>        Julian Day started at noon.

  !>        The zeroth Julian Day here is 01.01.-4712.

  !     CALLING SEQUENCE
  !         call caldat366_day(julday, dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: julday"     Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental subroutine caldat366_day(julian,dd,mm,yy)

    implicit none

    integer(i4), intent(in)  :: julian
    integer(i4), intent(out) :: dd, mm, yy
    integer(i4), parameter   :: year=366
    integer(i4), dimension(12), parameter   :: months=(/31,29,31,30,31,30,31,31,30,31,30,31/)
    integer(i4)              :: remainder

    yy = julian/year - 4712
    remainder = mod(abs(julian),year) + 1

    do mm = 1,size(months)
       if (remainder .le. months(mm)) then
          exit
       end if
       remainder = remainder - months(mm)
    end do

    dd = remainder

  end subroutine caldat366_day

  
  ! ------------------------------------------------------------------

  !     NAME
  !         caldatJulian

  !     PURPOSE
  !>        \brief Day, month and year from Julian day

  !>        \details Inverse of the function juldayJulian. Here julian is input as a Julian Day Number,
  !>        and the routine outputs id, mm, and yy as the day, month, and year on which the specified
  !>        Julian Day started at noon.

  !>        The zeroth Julian Day is 01.01.-4712, i.e. the 1st January 4713 BC.

  !>        Julian day definition starts at 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n

  !     CALLING SEQUENCE
  !         call caldatJulian(Julday, dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: Julday"     Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !>        \note Julian day definition starts at noon of the 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n
  !>        \n
  !>        julday and caldat start at midnight of the 1st January 4713 BC.
  !>        So date2decJulian and juldayJulian as well as dec2dateJulian and caldatJulian are shifted by half a day.\n
  !>        Use date2decJulian with dec2dateJulian together for fractional Julian dates
  !>        and use juldayJulian with caldatJulian together for integer Julian days.

  !     EXAMPLE
  !         ! 2415021 is 01.01.1900
  !         call caldatJulian(2415021, dd, mm, yy)
  !         -> see also example in test directory

  !     LITERATURE
  !         http://de.wikipedia.org/wiki/Julianisches_Datum
  !         which is different to the english Wiki
  !             http://en.wikipedia.org/wiki/Julian_day
  !         It is essentially the same as Numerical Recipes but uses astronomical instead of historical units.

  !     HISTORY
  !>        \author Written, Matthias Cuntz - modified julday from Numerical Recipes
  !>        \date Dec 2011
  !>        Modified Matthias Cuntz, May 2014 - changed to new algorithm with astronomical units
  !>                                            removed numerical recipes
  !>                 David Schaefer, Jan 2016 - renamed procedure
  elemental subroutine caldatJulian(julian,dd,mm,yy)

    implicit none

    integer(i4), intent(in)  :: julian
    integer(i4), intent(out) :: dd, mm, yy

    integer(i8) :: a, b, c, d, e, g
    integer(i4), parameter :: IGREG = 2299161_i4

    if (julian < IGREG) then
       A = int(julian, i8) ! julian
    else
       g = int((real(julian,dp)-1867216.25_dp)/36524.25_dp, i8) ! gregorian
       A = julian + 1_i8 + g - g/4_i8
    endif

    B = A + 1524_i8
    C = int((real(B,dp)-122.1_dp) / 365.25_dp, i8)
    D = int(365.25_dp * real(C,dp), i8)
    E = int(real(B-D,dp) / 30.6001_dp, i8)

    dd = int(B - D - int(30.6001_dp*real(E,dp), i8), i4)

    if (E<14_i8) then
       mm = int(E-1_i8, i4)
    else
       mm = int(E-13_i8, i4)
    endif

    if (mm > 2) then
       yy = int(C - 4716_i8, i4)
    else
       yy = int(C - 4715_i8, i4)
    endif

  end subroutine caldatJulian


  ! ------------------------------------------------------------------

  !     NAME
  !         caldatLilian

  !     PURPOSE
  !>        \brief Day, month and year from Lilian date

  !>        \details Inverse of the function juldayLilian. Here lilian is input as a Lilian Date,
  !>        and the routine outputs dd, mm, and yy as the day, month, and year on which the specified
  !>        Lilian Date started at midnight.

  !>        The first Lilian Day is 15.10.1582,
  !>        i.e. the day Pope Gregory XIII introduced the Gregorian calendar.

  !     CALLING SEQUENCE
  !         call caldatLilian(lilian, dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: lilian"     Lilian date

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4) :: dd"         Day in month of Lilian day
  !>        \param[out] "integer(i4) :: mm"         Month in year of Lilian day
  !>        \param[out] "integer(i4) :: yy"         Year of Lilian day

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! 115861 is 01.01.1900
  !         call caldatLilian(115861, dd, mm, yy)
  !         -> see also example in test directory

  !     LITERATURE
  !         https://en.wikipedia.org/wiki/Lilian_date

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Dec 2016
  elemental subroutine caldatLilian(lilian,dd,mm,yy)

    implicit none

    integer(i4), intent(in)  :: lilian
    integer(i4), intent(out) :: dd, mm, yy

    integer(i4), parameter :: iglil = 2299160_i4

    call caldatJulian(lilian+iglil,dd,mm,yy)

  end subroutine caldatLilian

  
  ! ------------------------------------------------------------------
  ! DATE2DEC
  ! ------------------------------------------------------------------


  ! ------------------------------------------------------------------

  !     NAME
  !         date2dec360day

  !     PURPOSE
  !>        \brief Fractional Julian day from day, month, year, hour, minute, second in 360 day calendar.

  !>        \details In this routine date2dec360day returns the fractional Julian Day that begins at noon
  !>        of the calendar date specified by month mm, day dd, and year yy, all integer variables.\n
  !>        Hours hh, minutes nn, seconds ss, or optinially fractional day can be given.

  !>        The zeroth Julian Day is 01.01.0000 at noon.

  !     CALLING SEQUENCE
  !         date2dec360day = date2dec360day(dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: dd"         Day in month of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: mm"         Month in year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: yy"         Year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: hh"         Hours of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: nn"         Minutes of hour of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: ss"         Secondes of minute of hour of Julian day (default: 0)
  !>        \param[in] "real(dp),    optional :: fracday"    Fractional day if hh,nn,ss not given (default: 0)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: date2dec360day  !     Fractional Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, David Schaefer
  !>        \date Oct 2015
  !>        Modified Matthias Cuntz, May 2014 - fractional day
  elemental function date2dec360day(dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    integer(i4), intent(in), optional :: dd, mm, yy
    integer(i4), intent(in), optional :: hh, nn, ss
    real(dp),    intent(in), optional :: fracday
    real(dp)                          :: date2dec360day, eps
    integer(i4)                       :: idd, imm, iyy
    real(dp)                          :: ihh, inn, iss
    real(dp)                          :: H, hour

    ! Presets
    idd = 1
    if (present(dd)) idd = dd
    imm = 1
    if (present(mm)) imm = mm
    iyy = 1
    if (present(yy)) iyy = yy
    if (present(hh) .or. present(nn) .or. present(ss)) then
       ihh = 0.0_dp
       if (present(hh)) ihh = real(hh,dp)
       inn = 0.0_dp
       if (present(nn)) inn = real(nn,dp)
       iss = 0.0_dp
       if (present(ss)) iss = real(ss,dp)
       H = ihh/24._dp + inn/1440._dp + iss/86400._dp
    else
       if (present(fracday)) then
          H = fracday
       else
          H = 0.0_dp
       endif
    endif

    hour = H - .5_dp

    ! Fractional Julian day starts at noon
    date2dec360day = real(julday360day(idd,imm,iyy),dp) + hour

    ! Add a small offset (proportional to julian date) for correct re-conversion.
    eps = epsilon(1.0_dp)
    eps = max(eps * abs(date2dec360day), eps)
    date2dec360day = date2dec360day + eps

  end function date2dec360day


  ! ------------------------------------------------------------------

  !     NAME
  !         date2dec365day

  !     PURPOSE
  !>        \brief Fractional Julian day from day, month, year, hour, minute, second in 365 day calendar.

  !>        \details In this routine date2dec365day returns the fractional Julian Day that begins at noon
  !>        of the calendar date specified by month mm, day dd, and year yy, all integer variables.\n
  !>        Hours hh, minutes nn, seconds ss, or optinially fractional day can be given.

  !>        The zeroth Julian Day is 01.01.0000 at noon.

  !     CALLING SEQUENCE
  !         date2dec365day = date2dec365day(dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: dd"         Day in month of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: mm"         Month in year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: yy"         Year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: hh"         Hours of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: nn"         Minutes of hour of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: ss"         Secondes of minute of hour of Julian day (default: 0)
  !>        \param[in] "real(dp),    optional :: fracday"    Fractional day if hh,nn,ss not given (default: 0)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: date2dec365day  !     Fractional Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, David Schaefer
  !>        \date Dec 2015
  !>        Modified Matthias Cuntz, May 2014 - fractional day
  elemental function date2dec365day(dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    integer(i4), intent(in), optional :: dd, mm, yy
    integer(i4), intent(in), optional :: hh, nn, ss
    real(dp),    intent(in), optional :: fracday
    real(dp)                          :: date2dec365day, eps
    integer(i4)                       :: idd, imm, iyy
    real(dp)                          :: ihh, inn, iss
    real(dp)                          :: H, hour

    ! Presets
    idd = 1
    if (present(dd)) idd = dd
    imm = 1
    if (present(mm)) imm = mm
    iyy = 1
    if (present(yy)) iyy = yy
    if (present(hh) .or. present(nn) .or. present(ss)) then
       ihh = 0.0_dp
       if (present(hh)) ihh = real(hh,dp)
       inn = 0.0_dp
       if (present(nn)) inn = real(nn,dp)
       iss = 0.0_dp
       if (present(ss)) iss = real(ss,dp)
       H = ihh/24._dp + inn/1440._dp + iss/86400._dp
    else
       if (present(fracday)) then
          H = fracday
       else
          H = 0.0_dp
       endif
    endif

    hour = H - .5_dp

    ! Fractional Julian day starts at noon
    date2dec365day = real(julday365day(idd,imm,iyy),dp) + hour

    ! Add a small offset (proportional to julian date) for correct re-conversion.
    eps = epsilon(1.0_dp)
    eps = max(eps * abs(date2dec365day), eps)
    date2dec365day = date2dec365day + eps

  end function date2dec365day


  ! ------------------------------------------------------------------

  !     NAME
  !         date2dec366day

  !     PURPOSE
  !>        \brief Fractional Julian day from day, month, year, hour, minute, second in 366 day calendar.

  !>        \details In this routine date2dec366day returns the fractional Julian Day that begins at noon
  !>        of the calendar date specified by month mm, day dd, and year yy, all integer variables.\n
  !>        Hours hh, minutes nn, seconds ss, or optinially fractional day can be given.

  !>        The zeroth Julian Day is 01.01.0000 at noon.

  !     CALLING SEQUENCE
  !         date2dec366day = date2dec366day(dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: dd"         Day in month of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: mm"         Month in year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: yy"         Year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: hh"         Hours of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: nn"         Minutes of hour of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: ss"         Secondes of minute of hour of Julian day (default: 0)
  !>        \param[in] "real(dp),    optional :: fracday"    Fractional day if hh,nn,ss not given (default: 0)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: date2dec366day  !     Fractional Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental function date2dec366day(dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    integer(i4), intent(in), optional :: dd, mm, yy
    integer(i4), intent(in), optional :: hh, nn, ss
    real(dp),    intent(in), optional :: fracday
    real(dp)                          :: date2dec366day, eps
    integer(i4)                       :: idd, imm, iyy
    real(dp)                          :: ihh, inn, iss
    real(dp)                          :: H, hour

    ! Presets
    idd = 1
    if (present(dd)) idd = dd
    imm = 1
    if (present(mm)) imm = mm
    iyy = 1
    if (present(yy)) iyy = yy
    if (present(hh) .or. present(nn) .or. present(ss)) then
       ihh = 0.0_dp
       if (present(hh)) ihh = real(hh,dp)
       inn = 0.0_dp
       if (present(nn)) inn = real(nn,dp)
       iss = 0.0_dp
       if (present(ss)) iss = real(ss,dp)
       H = ihh/24._dp + inn/1440._dp + iss/86400._dp
    else
       if (present(fracday)) then
          H = fracday
       else
          H = 0.0_dp
       endif
    endif

    hour = H - .5_dp

    ! Fractional Julian day starts at noon
    date2dec366day = real(julday366day(idd,imm,iyy),dp) + hour

    ! Add a small offset (proportional to julian date) for correct re-conversion.
    eps = epsilon(1.0_dp)
    eps = max(eps * abs(date2dec366day), eps)
    date2dec366day = date2dec366day + eps

  end function date2dec366day


  ! ------------------------------------------------------------------

  !     NAME
  !         date2dec360_day

  !     PURPOSE
  !>        \brief Fractional Julian day from day, month, year, hour, minute, second
  !>        in a Gregorian calendar with only 30-day months.

  !>        \details In this routine date2dec360_day returns the fractional Julian Day that begins at noon
  !>        of the calendar date specified by month mm, day dd, and year yy, all integer variables.\n
  !>        Hours hh, minutes nn, seconds ss, or optinially fractional day can be given.

  !>        The zeroth Julian Day is 01.01.-4712 at noon.

  !     CALLING SEQUENCE
  !         date2dec360_day = date2dec360_day(dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: dd"         Day in month of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: mm"         Month in year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: yy"         Year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: hh"         Hours of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: nn"         Minutes of hour of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: ss"         Secondes of minute of hour of Julian day (default: 0)
  !>        \param[in] "real(dp),    optional :: fracday"    Fractional day if hh,nn,ss not given (default: 0)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: date2dec360_day  !     Fractional Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental function date2dec360_day(dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    integer(i4), intent(in), optional :: dd, mm, yy
    integer(i4), intent(in), optional :: hh, nn, ss
    real(dp),    intent(in), optional :: fracday
    real(dp)                          :: date2dec360_day, eps
    integer(i4)                       :: idd, imm, iyy
    real(dp)                          :: ihh, inn, iss
    real(dp)                          :: H, hour

    ! Presets
    idd = 1
    if (present(dd)) idd = dd
    imm = 1
    if (present(mm)) imm = mm
    iyy = 1
    if (present(yy)) iyy = yy
    if (present(hh) .or. present(nn) .or. present(ss)) then
       ihh = 0.0_dp
       if (present(hh)) ihh = real(hh,dp)
       inn = 0.0_dp
       if (present(nn)) inn = real(nn,dp)
       iss = 0.0_dp
       if (present(ss)) iss = real(ss,dp)
       H = ihh/24._dp + inn/1440._dp + iss/86400._dp
    else
       if (present(fracday)) then
          H = fracday
       else
          H = 0.0_dp
       endif
    endif

    hour = H - .5_dp

    ! Fractional Julian day starts at noon
    date2dec360_day = real(julday360_day(idd,imm,iyy),dp) + hour

    ! Add a small offset (proportional to julian date) for correct re-conversion.
    eps = epsilon(1.0_dp)
    eps = max(eps * abs(date2dec360_day), eps)
    date2dec360_day = date2dec360_day + eps

  end function date2dec360_day


  ! ------------------------------------------------------------------

  !     NAME
  !         date2dec365_day

  !     PURPOSE
  !>        \brief Fractional Julian day from day, month, year, hour, minute, second
  !>        in a Gregorian calendar without leap years.

  !>        \details In this routine date2dec365_day returns the fractional Julian Day that begins at noon
  !>        of the calendar date specified by month mm, day dd, and year yy, all integer variables.\n
  !>        Hours hh, minutes nn, seconds ss, or optinially fractional day can be given.

  !>        The zeroth Julian Day is 01.01.0000 at noon.

  !     CALLING SEQUENCE
  !         date2dec365_day = date2dec365_day(dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: dd"         Day in month of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: mm"         Month in year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: yy"         Year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: hh"         Hours of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: nn"         Minutes of hour of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: ss"         Secondes of minute of hour of Julian day (default: 0)
  !>        \param[in] "real(dp),    optional :: fracday"    Fractional day if hh,nn,ss not given (default: 0)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: date2dec365_day  !     Fractional Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental function date2dec365_day(dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    integer(i4), intent(in), optional :: dd, mm, yy
    integer(i4), intent(in), optional :: hh, nn, ss
    real(dp),    intent(in), optional :: fracday
    real(dp)                          :: date2dec365_day, eps
    integer(i4)                       :: idd, imm, iyy
    real(dp)                          :: ihh, inn, iss
    real(dp)                          :: H, hour

    ! Presets
    idd = 1
    if (present(dd)) idd = dd
    imm = 1
    if (present(mm)) imm = mm
    iyy = 1
    if (present(yy)) iyy = yy
    if (present(hh) .or. present(nn) .or. present(ss)) then
       ihh = 0.0_dp
       if (present(hh)) ihh = real(hh,dp)
       inn = 0.0_dp
       if (present(nn)) inn = real(nn,dp)
       iss = 0.0_dp
       if (present(ss)) iss = real(ss,dp)
       H = ihh/24._dp + inn/1440._dp + iss/86400._dp
    else
       if (present(fracday)) then
          H = fracday
       else
          H = 0.0_dp
       endif
    endif

    hour = H - .5_dp

    ! Fractional Julian day starts at noon
    date2dec365_day = real(julday365_day(idd,imm,iyy),dp) + hour

    ! Add a small offset (proportional to julian date) for correct re-conversion.
    eps = epsilon(1.0_dp)
    eps = max(eps * abs(date2dec365_day), eps)
    date2dec365_day = date2dec365_day + eps

  end function date2dec365_day


  ! ------------------------------------------------------------------

  !     NAME
  !         date2dec366_day

  !     PURPOSE
  !>        \brief Fractional Julian day from day, month, year, hour, minute, second
  !>        in a Gregorian calendar with only leap years.

  !>        \details In this routine date2dec366_day returns the fractional Julian Day that begins at noon
  !>        of the calendar date specified by month mm, day dd, and year yy, all integer variables.\n
  !>        Hours hh, minutes nn, seconds ss, or optinially fractional day can be given.

  !>        The zeroth Julian Day is 01.01.0000 at noon.

  !     CALLING SEQUENCE
  !         date2dec366_day = date2dec366_day(dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: dd"         Day in month of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: mm"         Month in year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: yy"         Year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: hh"         Hours of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: nn"         Minutes of hour of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: ss"         Secondes of minute of hour of Julian day (default: 0)
  !>        \param[in] "real(dp),    optional :: fracday"    Fractional day if hh,nn,ss not given (default: 0)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: date2dec366_day  !     Fractional Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental function date2dec366_day(dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    integer(i4), intent(in), optional :: dd, mm, yy
    integer(i4), intent(in), optional :: hh, nn, ss
    real(dp),    intent(in), optional :: fracday
    real(dp)                          :: date2dec366_day, eps
    integer(i4)                       :: idd, imm, iyy
    real(dp)                          :: ihh, inn, iss
    real(dp)                          :: H, hour

    ! Presets
    idd = 1
    if (present(dd)) idd = dd
    imm = 1
    if (present(mm)) imm = mm
    iyy = 1
    if (present(yy)) iyy = yy
    if (present(hh) .or. present(nn) .or. present(ss)) then
       ihh = 0.0_dp
       if (present(hh)) ihh = real(hh,dp)
       inn = 0.0_dp
       if (present(nn)) inn = real(nn,dp)
       iss = 0.0_dp
       if (present(ss)) iss = real(ss,dp)
       H = ihh/24._dp + inn/1440._dp + iss/86400._dp
    else
       if (present(fracday)) then
          H = fracday
       else
          H = 0.0_dp
       endif
    endif

    hour = H - .5_dp

    ! Fractional Julian day starts at noon
    date2dec366_day = real(julday366_day(idd,imm,iyy),dp) + hour

    ! Add a small offset (proportional to julian date) for correct re-conversion.
    eps = epsilon(1.0_dp)
    eps = max(eps * abs(date2dec366_day), eps)
    date2dec366_day = date2dec366_day + eps

  end function date2dec366_day


  ! ------------------------------------------------------------------

  !     NAME
  !         date2decJulian

  !     PURPOSE
  !>        \brief Fractional Julian day from day, month, year, hour, minute, second

  !>        \details In this routine date2decJulian returns the fractional Julian Day that begins at noon
  !>        of the calendar date specified by month mm, day dd, and year yy, all integer variables.\n
  !>        Hours hh, minutes nn, seconds ss, or optinially fractional day can be given.

  !>        The zeroth Julian Day is 01.01.-4712 at noon, i.e. the 1st January 4713 BC 12:00:00 h.

  !>        Julian day definition starts at noon of the 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n

  !     CALLING SEQUENCE
  !         date2dec = date2decJulian(dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: dd"         Day in month of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: mm"         Month in year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: yy"         Year of Julian day (default: 1)
  !>        \param[in] "integer(i4), optional :: hh"         Hours of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: nn"         Minutes of hour of Julian day (default: 0)
  !>        \param[in] "integer(i4), optional :: ss"         Secondes of minute of hour of Julian day (default: 0)
  !>        \param[in] "real(dp),    optional :: fracday"    Fractional day if hh,nn,ss not given (default: 0)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: date2dec &mdash;     Fractional Julian day

  !     RESTRICTIONS
  !>        \note Julian day definition starts at noon of the 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n
  !>        \n
  !>        juldayJulian and caldatJulian start at midnight of the 1st January 4713 BC.
  !>        So date2decJulian and juldayJulian as well as dec2dateJulian and caldatJulian are shifted by half a day.\n
  !>        Use date2decJulian with dec2dateJulian together for fractional Julian dates
  !>        and use juldayJulian with caldatJulian together for integer Julian days.

  !     EXAMPLE
  !         ! 2415020.5 is 01.01.1900 00:00
  !         julian = date2decJulian(01,01,1990)
  !         ! 2415021.0 is 01.01.1900 12:00
  !         julian = date2decJulian(01,01,1990,12,00)
  !         -> see also example in test directory

  !     LITERATURE
  !         http://de.wikipedia.org/wiki/Julianisches_Datum
  !         which is different to the english Wiki
  !             http://en.wikipedia.org/wiki/Julian_day
  !         Numerical regulation of fractions is after
  !             IDL routine julday.pro. Copyright (c) 1988-2011, ITT Visual Information Solutions.

  !     HISTORY
  !>        \author Written,  Matthias Cuntz
  !>        \date Jan 2013
  !>        Modified Matthias Cuntz, May 2014 - changed to new algorithm with astronomical units
  !>                                            removed numerical recipes
  !>                 David Schaefer, Jan 2016 - renamed procedure
  !>                 Matthias Cuntz, Dec 2016 - fractional day
  elemental function date2decJulian(dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    integer(i4), intent(in), optional :: dd, mm, yy
    integer(i4), intent(in), optional :: hh, nn, ss
    real(dp),    intent(in), optional :: fracday
    real(dp)                          :: date2decjulian

    integer(i4), parameter :: IGREG2 = 15 + 31*(10+12*1582)
    integer(i4), parameter :: IGREG1 =  4 + 31*(10+12*1582)
    integer(i4) :: idd, imm, iyy
    real(dp)    :: ihh, inn, iss
    integer(i8) :: jm, jy
    real(dp)    :: jd, H, eps
    integer(i8) :: A, B

    ! Presets
    idd = 1
    if (present(dd)) idd = dd
    imm = 1
    if (present(mm)) imm = mm
    iyy = 1
    if (present(yy)) iyy = yy
    if (present(hh) .or. present(nn) .or. present(ss)) then
       ihh = 0.0_dp
       if (present(hh)) ihh = real(hh,dp)
       inn = 0.0_dp
       if (present(nn)) inn = real(nn,dp)
       iss = 0.0_dp
       if (present(ss)) iss = real(ss,dp)
       H = ihh/24._dp + inn/1440._dp + iss/86400._dp
    else
       if (present(fracday)) then
          H = fracday
       else
          H = 0.0_dp
       endif
    endif

    if (imm > 2) then
       jm = int(imm, i8)
       jy = int(iyy, i8)
    else
       jm = int(imm+12, i8)
       jy = int(iyy-1, i8)
    endif

    jd = real(idd, dp)

    if (dd+31*(mm+12*yy) >= IGREG2) then ! gregorian
       A = jy/100_i8
       B = 2_i8 - A + A/4_i8
    else if (dd+31*(mm+12*yy) <= IGREG1) then ! julian
       B = 0_i8
       ! else
       !    stop 'No Gregorian dates between 04.10.1582 and 15.10.1582'
    endif

    ! Fractional Julian day starts at noon
    date2decJulian = floor(365.25_dp*real(jy+4716_i8,dp)) + floor(30.6001_dp*real(jm+1_i8,dp)) + jd + H + real(B,dp) - 1524.5_dp

    ! Add a small offset (proportional to julian date) for correct re-conversion.
    eps = epsilon(1.0_dp)
    eps = max(eps * abs(date2decJulian), eps)
    date2decJulian = date2decJulian + eps

  end function date2decJulian


  ! ------------------------------------------------------------------

  !     NAME
  !         date2decLilian

  !     PURPOSE
  !>        \brief Fractional Lilian day from day, month, year, hour, minute, second

  !>        \details In this routine date2decLilian returns the fractional Lilian Day that begins at midnight
  !>        of the calendar date specified by month mm, day dd, and year yy, all integer variables.
  !>        Hours hh, minutes nn, seconds ss, or optinially fractional day can be given.

  !>        The first Lilian Day is 15.10.1582,
  !>        i.e. the day Pope Gregory XIII introduced the Gregorian calendar.

  !     CALLING SEQUENCE
  !         date2dec = date2decLilian(dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: dd"         Day in month of Lilian day (default: 1)
  !>        \param[in] "integer(i4), optional :: mm"         Month in year of Lilian day (default: 1)
  !>        \param[in] "integer(i4), optional :: yy"         Year of Lilian day (default: 1)
  !>        \param[in] "integer(i4), optional :: hh"         Hours of Lilian day (default: 0)
  !>        \param[in] "integer(i4), optional :: nn"         Minutes of hour of Lilian day (default: 0)
  !>        \param[in] "integer(i4), optional :: ss"         Secondes of minute of hour of Lilian day (default: 0)
  !>        \param[in] "real(dp),    optional :: fracday"    Fractional day if hh,nn,ss not given (default: 0)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: date2dec &mdash;     Fractional Lilian day

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! 115861.5 is 01.01.1900 12:00
  !         Lilian = date2decLilian(01,01,1990,12,00)
  !         -> see also example in test directory

  !     LITERATURE
  !         https://en.wikipedia.org/wiki/Lilian_date

  !     HISTORY
  !>        \author Written,  Matthias Cuntz
  !>        \date Dec 2016
  elemental function date2decLilian(dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    integer(i4), intent(in), optional :: dd, mm, yy
    integer(i4), intent(in), optional :: hh, nn, ss
    real(dp),    intent(in), optional :: fracday
    real(dp)                          :: date2decLilian

    real(dp), parameter :: IGLIL = 2299159.5_dp

    date2decLilian = date2decJulian(dd,mm,yy,hh,nn,ss,fracday) - IGLIL

  end function date2decLilian

  
  ! ------------------------------------------------------------------
  ! DEC2DATE
  ! ------------------------------------------------------------------


  ! ------------------------------------------------------------------

  !     NAME
  !         dec2date360day

  !     PURPOSE
  !>        \brief Day, month, year, hour, minute, and second from fractional Julian day in a 360 day calendar.

  !>        \details Inverse of the function date2dec360day. Here dec2date360day is input as a fractional Julian Day.
  !>        The routine outputs dd, mm, yy, hh, nn, ss as the day, month, year, hour, minute, and second
  !>        on which the specified Julian Day started at noon.\n
  !>        Fractional day can be output optionally.
  
  !>        The zeroth Day is 01.01.0000 at noon.

  !     CALLING SEQUENCE
  !         call dec2date360day(fJulian, dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: fJulian"     fractional Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "integer(i4), optional :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4), optional :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4), optional :: yy"         Year of Julian day
  !>        \param[out] "integer(i4), optional :: hh"         Hour of Julian day
  !>        \param[out] "integer(i4), optional :: nn"         Minute in hour of Julian day
  !>        \param[out] "integer(i4), optional :: ss"         Second in minute of hour of Julian day
  !>        \param[out] "real(dp),    optional :: fracday"    Fractional day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, David Schaefer
  !>        \date Oct 2015
  !>        Modified Matthias Cuntz, May 2014 - fractional day
  elemental subroutine dec2date360day(julian, dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    real(dp),    intent(in)            :: julian
    integer(i4), intent(out), optional :: dd, mm, yy
    integer(i4), intent(out), optional :: hh, nn, ss
    real(dp),    intent(out), optional :: fracday
    integer(i4)                        :: day, month, year
    real(dp)                           :: fraction, fJulian
    integer(i4)                        :: hour, minute, second

    fJulian = julian + .5_dp
    call caldat360day(int(floor(fJulian),i4),day,month,year)

    fraction = fJulian - floor(fJulian)
    hour     = min(max(floor(fraction * 24.0_dp), 0), 23)
    fraction = fraction - real(hour,dp)/24.0_dp
    minute   = min(max(floor(fraction*1440.0_dp), 0), 59)
    second   = max(nint((fraction - real(minute,dp)/1440.0_dp)*86400.0_dp), 0)

    if (second==60) then
       second = 0
       minute = minute + 1
       if (minute==60) then
          minute = 0
          hour   = hour + 1
          if (hour==24) then
             hour = 0
             call caldat360day(julday360day(day, month, year) + 1, day, month, year)
          endif
       endif
    endif

    if (present(dd)) dd = day
    if (present(mm)) mm = month
    if (present(yy)) yy = year
    if (present(hh)) hh = hour
    if (present(nn)) nn = minute
    if (present(ss)) ss = second
    if (present(fracday)) fracday = real(hour,dp)/24._dp + real(minute,dp)/1440._dp + real(second,dp)/86400._dp

  end subroutine dec2date360day


  ! ------------------------------------------------------------------

  !     NAME
  !         dec2date365day

  !     PURPOSE
  !>        \brief Day, month, year, hour, minute, and second from fractional Julian day in a 365_day calendar

  !>        \details Inverse of the function date2dec. Here dec2date365day is input as a fractional Julian Day.
  !>        The routine outputs dd, mm, yy, hh, nn, ss as the day, month, year, hour, minute, and second
  !>        on which the specified Julian Day started at noon.\n
  !>        Fractional day can be output optionally.

  !>        The zeroth Julian Day is 01.01.0000 at noon.

  !     CALLING SEQUENCE
  !         call dec2date365day(fJulian, dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: fJulian"     fractional Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "integer(i4), optional :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4), optional :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4), optional :: yy"         Year of Julian day
  !>        \param[out] "integer(i4), optional :: hh"         Hour of Julian day
  !>        \param[out] "integer(i4), optional :: nn"         Minute in hour of Julian day
  !>        \param[out] "integer(i4), optional :: ss"         Second in minute of hour of Julian day
  !>        \param[out] "real(dp),    optional :: fracday"    Fractional day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, David Schaefer
  !>        \date Dec 2015
  !>        Modified Matthias Cuntz, May 2014 - fractional day
  elemental subroutine dec2date365day(julian, dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    real(dp),    intent(in)            :: julian
    integer(i4), intent(out), optional :: dd, mm, yy
    integer(i4), intent(out), optional :: hh, nn, ss
    real(dp),    intent(out), optional :: fracday
    integer(i4)                        :: day, month, year
    real(dp)                           :: fraction, fJulian
    integer(i4)                        :: hour, minute, second

    fJulian = julian + .5_dp
    call caldat365day(int(floor(fJulian),i4),day,month,year)

    fraction = fJulian - floor(fJulian)
    hour     = min(max(floor(fraction * 24.0_dp), 0), 23)
    fraction = fraction - real(hour,dp)/24.0_dp
    minute   = min(max(floor(fraction*1440.0_dp), 0), 59)
    second   = max(nint((fraction - real(minute,dp)/1440.0_dp)*86400.0_dp), 0)

    ! If seconds==60
    if (second==60) then
       second = 0
       minute = minute + 1
       if (minute==60) then
          minute = 0
          hour   = hour + 1
          if (hour==24) then
             hour = 0
             call caldat365day(julday365day(day, month, year) + 1, day, month, year)
          endif
       endif
    endif

    if (present(dd)) dd = day
    if (present(mm)) mm = month
    if (present(yy)) yy = year
    if (present(hh)) hh = hour
    if (present(nn)) nn = minute
    if (present(ss)) ss = second
    if (present(fracday)) fracday = real(hour,dp)/24._dp + real(minute,dp)/1440._dp + real(second,dp)/86400._dp

  end subroutine dec2date365day


  ! ------------------------------------------------------------------

  !     NAME
  !         dec2date366day

  !     PURPOSE
  !>        \brief Day, month, year, hour, minute, and second from fractional Julian day in a 366_day calendar

  !>        \details Inverse of the function date2dec. Here dec2date366day is input as a fractional Julian Day.
  !>        The routine outputs dd, mm, yy, hh, nn, ss as the day, month, year, hour, minute, and second
  !>        on which the specified Julian Day started at noon.\n
  !>        Fractional day can be output optionally.

  !>        The zeroth Julian Day is 01.01.0000 at noon.

  !     CALLING SEQUENCE
  !         call dec2date366day(fJulian, dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: fJulian"     fractional Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "integer(i4), optional :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4), optional :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4), optional :: yy"         Year of Julian day
  !>        \param[out] "integer(i4), optional :: hh"         Hour of Julian day
  !>        \param[out] "integer(i4), optional :: nn"         Minute in hour of Julian day
  !>        \param[out] "integer(i4), optional :: ss"         Second in minute of hour of Julian day
  !>        \param[out] "real(dp),    optional :: fracday"    Fractional day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental subroutine dec2date366day(julian, dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    real(dp),    intent(in)            :: julian
    integer(i4), intent(out), optional :: dd, mm, yy
    integer(i4), intent(out), optional :: hh, nn, ss
    real(dp),    intent(out), optional :: fracday
    integer(i4)                        :: day, month, year
    real(dp)                           :: fraction, fJulian
    integer(i4)                        :: hour, minute, second

    fJulian = julian + .5_dp
    call caldat366day(int(floor(fJulian),i4),day,month,year)

    fraction = fJulian - floor(fJulian)
    hour     = min(max(floor(fraction * 24.0_dp), 0), 23)
    fraction = fraction - real(hour,dp)/24.0_dp
    minute   = min(max(floor(fraction*1440.0_dp), 0), 59)
    second   = max(nint((fraction - real(minute,dp)/1440.0_dp)*86400.0_dp), 0)

    ! If seconds==60
    if (second==60) then
       second = 0
       minute = minute + 1
       if (minute==60) then
          minute = 0
          hour   = hour + 1
          if (hour==24) then
             hour = 0
             call caldat366day(julday366day(day, month, year) + 1, day, month, year)
          endif
       endif
    endif

    if (present(dd)) dd = day
    if (present(mm)) mm = month
    if (present(yy)) yy = year
    if (present(hh)) hh = hour
    if (present(nn)) nn = minute
    if (present(ss)) ss = second
    if (present(fracday)) fracday = real(hour,dp)/24._dp + real(minute,dp)/1440._dp + real(second,dp)/86400._dp

  end subroutine dec2date366day


  ! ------------------------------------------------------------------

  !     NAME
  !         dec2date360_day

  !     PURPOSE
  !>        \brief Day, month, year, hour, minute, and second from fractional Julian day
  !>        in a Gregorian calendar with only 30-day months.
  !>        \details Inverse of the function date2dec360day. Here dec2date360_day is input as a fractional Julian Day.
  !>        The routine outputs dd, mm, yy, hh, nn, ss as the day, month, year, hour, minute, and second
  !>        on which the specified Julian Day started at noon.\n
  !>        Fractional day can be output optionally.
  
  !>        The zeroth Day is 01.01.0000 at noon.

  !     CALLING SEQUENCE
  !         call dec2date360_day(fJulian, dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: fJulian"     fractional Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "integer(i4), optional :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4), optional :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4), optional :: yy"         Year of Julian day
  !>        \param[out] "integer(i4), optional :: hh"         Hour of Julian day
  !>        \param[out] "integer(i4), optional :: nn"         Minute in hour of Julian day
  !>        \param[out] "integer(i4), optional :: ss"         Second in minute of hour of Julian day
  !>        \param[out] "real(dp),    optional :: fracday"    Fractional day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental subroutine dec2date360_day(julian, dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    real(dp),    intent(in)            :: julian
    integer(i4), intent(out), optional :: dd, mm, yy
    integer(i4), intent(out), optional :: hh, nn, ss
    real(dp),    intent(out), optional :: fracday
    integer(i4)                        :: day, month, year
    real(dp)                           :: fraction, fJulian
    integer(i4)                        :: hour, minute, second

    fJulian = julian + .5_dp
    call caldat360_day(int(floor(fJulian),i4),day,month,year)

    fraction = fJulian - floor(fJulian)
    hour     = min(max(floor(fraction * 24.0_dp), 0), 23)
    fraction = fraction - real(hour,dp)/24.0_dp
    minute   = min(max(floor(fraction*1440.0_dp), 0), 59)
    second   = max(nint((fraction - real(minute,dp)/1440.0_dp)*86400.0_dp), 0)

    if (second==60) then
       second = 0
       minute = minute + 1
       if (minute==60) then
          minute = 0
          hour   = hour + 1
          if (hour==24) then
             hour = 0
             call caldat360_day(julday360_day(day, month, year) + 1, day, month, year)
          endif
       endif
    endif

    if (present(dd)) dd = day
    if (present(mm)) mm = month
    if (present(yy)) yy = year
    if (present(hh)) hh = hour
    if (present(nn)) nn = minute
    if (present(ss)) ss = second
    if (present(fracday)) fracday = real(hour,dp)/24._dp + real(minute,dp)/1440._dp + real(second,dp)/86400._dp

  end subroutine dec2date360_day


  ! ------------------------------------------------------------------

  !     NAME
  !         dec2date365_day

  !     PURPOSE
  !>        \brief Day, month, year, hour, minute, and second from fractional Julian day
  !>        in a Gregorian calendar without leap years.

  !>        \details Inverse of the function date2dec. Here dec2date365_day is input as a fractional Julian Day.
  !>        The routine outputs dd, mm, yy, hh, nn, ss as the day, month, year, hour, minute, and second
  !>        on which the specified Julian Day started at noon.\n
  !>        Fractional day can be output optionally.

  !>        The zeroth Julian Day is 01.01.0000 at noon.

  !     CALLING SEQUENCE
  !         call dec2date365_day(fJulian, dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: fJulian"     fractional Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "integer(i4), optional :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4), optional :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4), optional :: yy"         Year of Julian day
  !>        \param[out] "integer(i4), optional :: hh"         Hour of Julian day
  !>        \param[out] "integer(i4), optional :: nn"         Minute in hour of Julian day
  !>        \param[out] "integer(i4), optional :: ss"         Second in minute of hour of Julian day
  !>        \param[out] "real(dp),    optional :: fracday"    Fractional day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental subroutine dec2date365_day(julian, dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    real(dp),    intent(in)            :: julian
    integer(i4), intent(out), optional :: dd, mm, yy
    integer(i4), intent(out), optional :: hh, nn, ss
    real(dp),    intent(out), optional :: fracday
    integer(i4)                        :: day, month, year
    real(dp)                           :: fraction, fJulian
    integer(i4)                        :: hour, minute, second

    fJulian = julian + .5_dp
    call caldat365_day(int(floor(fJulian),i4),day,month,year)

    fraction = fJulian - floor(fJulian)
    hour     = min(max(floor(fraction * 24.0_dp), 0), 23)
    fraction = fraction - real(hour,dp)/24.0_dp
    minute   = min(max(floor(fraction*1440.0_dp), 0), 59)
    second   = max(nint((fraction - real(minute,dp)/1440.0_dp)*86400.0_dp), 0)

    ! If seconds==60
    if (second==60) then
       second = 0
       minute = minute + 1
       if (minute==60) then
          minute = 0
          hour   = hour + 1
          if (hour==24) then
             hour = 0
             call caldat365_day(julday365_day(day, month, year) + 1, day, month, year)
          endif
       endif
    endif

    if (present(dd)) dd = day
    if (present(mm)) mm = month
    if (present(yy)) yy = year
    if (present(hh)) hh = hour
    if (present(nn)) nn = minute
    if (present(ss)) ss = second
    if (present(fracday)) fracday = real(hour,dp)/24._dp + real(minute,dp)/1440._dp + real(second,dp)/86400._dp

  end subroutine dec2date365_day


  ! ------------------------------------------------------------------

  !     NAME
  !         dec2date366_day

  !     PURPOSE
  !>        \brief Day, month, year, hour, minute, and second from fractional Julian day
  !>        in a Gregorian calendar with only leap years.

  !>        \details Inverse of the function date2dec. Here dec2date366_day is input as a fractional Julian Day.
  !>        The routine outputs dd, mm, yy, hh, nn, ss as the day, month, year, hour, minute, and second
  !>        on which the specified Julian Day started at noon.\n
  !>        Fractional day can be output optionally.

  !>        The zeroth Julian Day is 01.01.0000 at noon.

  !     CALLING SEQUENCE
  !         call dec2date366_day(fJulian, dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: fJulian"     fractional Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "integer(i4), optional :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4), optional :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4), optional :: yy"         Year of Julian day
  !>        \param[out] "integer(i4), optional :: hh"         Hour of Julian day
  !>        \param[out] "integer(i4), optional :: nn"         Minute in hour of Julian day
  !>        \param[out] "integer(i4), optional :: ss"         Second in minute of hour of Julian day
  !>        \param[out] "real(dp),    optional :: fracday"    Fractional day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental subroutine dec2date366_day(julian, dd, mm, yy, hh, nn, ss, fracday)

    implicit none

    real(dp),    intent(in)            :: julian
    integer(i4), intent(out), optional :: dd, mm, yy
    integer(i4), intent(out), optional :: hh, nn, ss
    real(dp),    intent(out), optional :: fracday
    integer(i4)                        :: day, month, year
    real(dp)                           :: fraction, fJulian
    integer(i4)                        :: hour, minute, second

    fJulian = julian + .5_dp
    call caldat366_day(int(floor(fJulian),i4),day,month,year)

    fraction = fJulian - floor(fJulian)
    hour     = min(max(floor(fraction * 24.0_dp), 0), 23)
    fraction = fraction - real(hour,dp)/24.0_dp
    minute   = min(max(floor(fraction*1440.0_dp), 0), 59)
    second   = max(nint((fraction - real(minute,dp)/1440.0_dp)*86400.0_dp), 0)

    ! If seconds==60
    if (second==60) then
       second = 0
       minute = minute + 1
       if (minute==60) then
          minute = 0
          hour   = hour + 1
          if (hour==24) then
             hour = 0
             call caldat366_day(julday366_day(day, month, year) + 1, day, month, year)
          endif
       endif
    endif

    if (present(dd)) dd = day
    if (present(mm)) mm = month
    if (present(yy)) yy = year
    if (present(hh)) hh = hour
    if (present(nn)) nn = minute
    if (present(ss)) ss = second
    if (present(fracday)) fracday = real(hour,dp)/24._dp + real(minute,dp)/1440._dp + real(second,dp)/86400._dp

  end subroutine dec2date366_day


  ! ------------------------------------------------------------------

  !     NAME
  !         dec2dateJulian

  !     PURPOSE
  !>        \brief Day, month, year, hour, minute, and second from fractional Julian day

  !>        \details Inverse of the function date2decJulian. Here dec2dateJulian is input as a fractional Julian Day,
  !>        which starts at noon of the 1st January 4713 BC, i.e. 01.01.-4712.
  !>        The routine outputs dd, mm, yy, hh, nn, ss as the day, month, year, hour, minute, and second
  !>        on which the specified Julian Day started at noon.\n
  !>        Fractional day can be output optionally.

  !>        The zeroth Julian Day is 01.01.-4712 at noon, i.e. the 1st January 4713 BC at noon.

  !>        Julian day definition starts at 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n

  !     CALLING SEQUENCE
  !         call dec2dateJulian(fJulian, dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: fJulian"     fractional Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "integer(i4), optional :: dd"         Day in month of Julian day
  !>        \param[out] "integer(i4), optional :: mm"         Month in year of Julian day
  !>        \param[out] "integer(i4), optional :: yy"         Year of Julian day
  !>        \param[out] "integer(i4), optional :: hh"         Hour of Julian day
  !>        \param[out] "integer(i4), optional :: nn"         Minute in hour of Julian day
  !>        \param[out] "integer(i4), optional :: ss"         Second in minute of hour of Julian day
  !>        \param[out] "real(dp),    optional :: fracday"    Fractional day

  !     RESTRICTIONS
  !>        \note Julian day definition starts at noon of the 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n
  !>        \n
  !>        juldayJulian and caldatJulian start at midnight of the 1st January 4713 BC.
  !>        So date2decJulian and juldayJulian as well as dec2dateJulian and caldatJulian are shifted by half a day.\n
  !>        Use date2decJulian with dec2dateJulian together for fractional Julian dates
  !>        and use juldayJulian with caldatJulian together for integer Julian days.

  !     EXAMPLE
  !         ! 2415020.5 is 01.01.1900 00:00
  !         call caldatJulian(2415020.5, dd, mm, yy, hh, nn)
  !         ! 2415021.0 is 01.01.1900 12:00
  !         call caldatJulian(2415021., dd, mm, yy, hh, nn, ss)
  !         -> see also example in test directory

  !     LITERATURE
  !         http://de.wikipedia.org/wiki/Julianisches_Datum
  !         which is different to the english Wiki
  !             http://en.wikipedia.org/wiki/Julian_day
  !         It is essentially the same as Numerical Recipes but uses astronomical instead of historical units.
  !         Here the sometimes 60 sec as output are corrected at the end.

  !     HISTORY
  !>        \author Written,  Matthias Cuntz
  !>        \date Jan 2013
  !>        Modified Matthias Cuntz, May 2014 - changed to new algorithm with astronomical units
  !>                                            removed numerical recipes
  !>                 David Schaefer, Jan 2016 - renamed procedure
  !>                 Matthias Cuntz, Dec 2016 - fractional day
  ELEMENTAL SUBROUTINE dec2dateJulian(julian, dd, mm, yy, hh, nn, ss, fracday)

    IMPLICIT NONE

    REAL(dp),    INTENT(IN)            :: julian
    INTEGER(i4), INTENT(OUT), OPTIONAL :: dd, mm, yy
    INTEGER(i4), INTENT(OUT), OPTIONAL :: hh, nn, ss
    REAL(dp),    INTENT(OUT), OPTIONAL :: fracday

    INTEGER(i4) :: day, month, year, hour, minute, second
    REAL(dp)    :: fraction

    INTEGER(i8) :: A, B, C, D, E, g, Z
    INTEGER(i4), PARAMETER :: IGREG = 2299161_i4

    Z = int(julian + 0.5, i8)

    if (Z < IGREG) then
       A = Z ! julian
    else
       g = int((real(Z,dp)-1867216.25_dp)/36524.25_dp, i8) ! gregorian
       A = Z + 1_i8 + g - g/4_i8
    endif

    B = A + 1524_i8
    C = int((real(B,dp)-122.1_dp) / 365.25_dp, i8)
    D = int(365.25_dp * real(C,dp), i8)
    E = int(real(B-D,dp) / 30.6001_dp, i8)

    day = int(B - D - int(30.6001_dp*real(E,dp), i8), i4)

    if (E<14_i8) then
       month = int(E-1_i8, i4)
    else
       month = int(E-13_i8, i4)
    endif

    if (month > 2) then
       year = int(C - 4716_i8, i4)
    else
       year = int(C - 4715_i8, i4)
    endif

    ! ! Fractional part
    ! eps = 1e-12_dp ! ~ 5000*epsilon(1.0_dp)
    ! eps = max(eps * abs(real(Z,dp)), eps)
    ! fraction = julian + 0.5_dp - real(Z,dp)
    ! hour     = min(max(floor(fraction * 24.0_dp + eps), 0), 23)
    ! fraction = fraction - real(hour,dp)/24.0_dp
    ! minute   = min(max(floor(fraction*1440.0_dp + eps), 0), 59)
    ! second   = max(nint((fraction - real(minute,dp)/1440.0_dp)*86400.0_dp), 0)

    ! Fractional part
    fraction = julian + 0.5_dp - real(Z,dp)
    hour     = min(max(floor(fraction * 24.0_dp), 0), 23)
    fraction = fraction - real(hour,dp)/24.0_dp
    minute   = min(max(floor(fraction*1440.0_dp), 0), 59)
    second   = max(nint((fraction - real(minute,dp)/1440.0_dp)*86400.0_dp), 0)

    ! If seconds==60
    if (second==60) then
       second = 0
       minute = minute + 1
       if (minute==60) then
          minute = 0
          hour   = hour + 1
          if (hour==24) then
             hour = 0
             call caldat(julday(day, month, year) + 1, day, month, year)
          endif
       endif
    endif

    if (present(dd)) dd = day
    if (present(mm)) mm = month
    if (present(yy)) yy = year
    if (present(hh)) hh = hour
    if (present(nn)) nn = minute
    if (present(ss)) ss = second
    if (present(fracday)) fracday = real(hour,dp)/24._dp + real(minute,dp)/1440._dp + real(second,dp)/86400._dp

  END SUBROUTINE dec2dateJulian


  ! ------------------------------------------------------------------

  !     NAME
  !         dec2dateLilian

  !     PURPOSE
  !>        \brief Day, month, year, hour, minute, and second from fractional Lilian day

  !>        \details Inverse of the function date2decLilian. Here dec2dateLilian is input as a fractional Lilian Day,
  !>        which starts at midnight of the 15th October 1582.
  !>        The routine outputs dd, mm, yy, hh, nn, ss as the day, month, year, hour, minute, and second
  !>        on which the specified Lilian Day started at midnight.\n
  !>        Fractional day can be output optionally.

  !>        The first Lilian Day is 15.10.1582,
  !>        i.e. the day Pope Gregory XIII introduced the Gregorian calendar.

  !     CALLING SEQUENCE
  !         call dec2dateLilian(fLilian, dd, mm, yy, hh, nn, ss, fracday)

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: fLilian"     fractional Lilian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "integer(i4), optional :: dd"         Day in month of Lilian day
  !>        \param[out] "integer(i4), optional :: mm"         Month in year of Lilian day
  !>        \param[out] "integer(i4), optional :: yy"         Year of Lilian day
  !>        \param[out] "integer(i4), optional :: hh"         Hour of Lilian day
  !>        \param[out] "integer(i4), optional :: nn"         Minute in hour of Lilian day
  !>        \param[out] "integer(i4), optional :: ss"         Second in minute of hour of Lilian day
  !>        \param[out] "real(dp),    optional :: fracday"    Fractional day

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! 115861.5 is 01.01.1900 12:00
  !         call caldatLilian(115861.5, dd, mm, yy, hh, nn, ss)
  !         -> see also example in test directory

  !     LITERATURE
  !         https://en.wikipedia.org/wiki/Lilian_date

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Dec 2016
  ELEMENTAL SUBROUTINE dec2dateLilian(lilian, dd, mm, yy, hh, nn, ss, fracday)

    IMPLICIT NONE

    REAL(dp),    INTENT(IN)            :: lilian
    INTEGER(i4), INTENT(OUT), OPTIONAL :: dd, mm, yy
    INTEGER(i4), INTENT(OUT), OPTIONAL :: hh, nn, ss
    REAL(dp),    INTENT(OUT), OPTIONAL :: fracday

    REAL(dp), PARAMETER :: IGLIL = 2299159.5_dp

    call dec2dateJulian(lilian+IGLIL,dd,mm,yy,hh,nn,ss,fracday)

  END SUBROUTINE dec2dateLilian

  
  ! ------------------------------------------------------------------
  ! JULDAY
  ! ------------------------------------------------------------------


  ! ------------------------------------------------------------------

  !     NAME
  !         julday360day

  !     PURPOSE
  !>        \brief Julian day from day, month and year in a 360 day calendar.

  !>        \details In this routine julday360day returns the Julian Day Number that begins at noon of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables.

  !>        The zeroth Julian Day is 01.01.0000

  !     CALLING SEQUENCE
  !         julian = julday360day(dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[in] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[in] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: julian  !     Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, David Schaefer
  !>        \date Oct 2015
  elemental function julday360day(dd,mm,yy)

    implicit none

    integer(i4), intent(in) :: dd, mm, yy
    integer(i4)             :: julday360day
    integer(i4), parameter  :: year=360, month=30

    julday360day = abs(yy)*year + (mm-1)*month + (dd-1)
    if (yy < 0) julday360day = julday360day * (-1)

  end function julday360day


  ! ------------------------------------------------------------------

  !     NAME
  !         julday365day

  !     PURPOSE
  !>        \brief Julian day from day, month and year in a 365 day calendar.

  !>        \details In this routine julday365day returns the Julian Day Number that begins at noon of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables.

  !>        The zeroth Julian Day is 01.01.0000

  !     CALLING SEQUENCE
  !         julian = julday365day(dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[in] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[in] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: julian  !     Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, David Schaefer
  !>        \date Dec 2015
  elemental function julday365day(dd,mm,yy)

    implicit none

    integer(i4), intent(in) :: dd, mm, yy
    integer(i4)             :: julday365day
    integer(i4), parameter  :: year=365
    integer(i4),dimension(12),parameter  :: months=(/31,28,31,30,31,30,31,31,30,31,30,31/)

    julday365day = abs(yy)*year +  sum(months(1:mm-1)) + (dd-1)

    if (yy < 0) julday365day = julday365day * (-1)

  end function julday365day


  ! ------------------------------------------------------------------

  !     NAME
  !         julday366day

  !     PURPOSE
  !>        \brief Julian day from day, month and year in a 366 day calendar.

  !>        \details In this routine julday366day returns the Julian Day Number that begins at noon of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables.

  !>        The zeroth Julian Day is 01.01.0000

  !     CALLING SEQUENCE
  !         julian = julday366day(dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[in] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[in] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: julian  !     Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental function julday366day(dd,mm,yy)

    implicit none

    integer(i4), intent(in) :: dd, mm, yy
    integer(i4)             :: julday366day
    integer(i4), parameter  :: year=366
    integer(i4),dimension(12),parameter  :: months=(/31,29,31,30,31,30,31,31,30,31,30,31/)

    julday366day = abs(yy)*year +  sum(months(1:mm-1)) + (dd-1)

    if (yy < 0) julday366day = julday366day * (-1)

  end function julday366day


  ! ------------------------------------------------------------------

  !     NAME
  !         julday360_day

  !     PURPOSE
  !>        \brief Julian day from day, month and year in a Gregorian calendar with only 30-day months.

  !>        \details In this routine julday360_day returns the Julian Day Number that begins at noon of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables.

  !>        The zeroth Julian Day is 01.01.0000

  !     CALLING SEQUENCE
  !         julian = julday360_day(dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[in] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[in] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: julian  !     Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental function julday360_day(dd,mm,yy)

    implicit none

    integer(i4), intent(in) :: dd, mm, yy
    integer(i4)             :: julday360_day
    integer(i4), parameter  :: year=360, month=30

    julday360_day = (yy+4712)*year + (mm-1)*month + (dd-1)

  end function julday360_day


  ! ------------------------------------------------------------------

  !     NAME
  !         julday365_day

  !     PURPOSE
  !>        \brief Julian day from day, month and year in a Gregorian calendar without leap years.
  
  !>        \details In this routine julday365_day returns the Julian Day Number that begins at noon of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables.

  !>        The zeroth Julian Day is 01.01.-4712.

  !     CALLING SEQUENCE
  !         julian = julday365_day(dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[in] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[in] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: julian  !     Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental function julday365_day(dd,mm,yy)

    implicit none

    integer(i4), intent(in) :: dd, mm, yy
    integer(i4)             :: julday365_day
    integer(i4), parameter  :: year=365
    integer(i4),dimension(12),parameter :: months=(/31,28,31,30,31,30,31,31,30,31,30,31/)

    julday365_day = (yy+4712)*year + sum(months(1:mm-1)) + (dd-1)

  end function julday365_day


  ! ------------------------------------------------------------------

  !     NAME
  !         julday366_day

  !     PURPOSE
  !>        \brief Julian day from day, month and year in a Gregorian calendar with only leap years.
  
  !>        \details In this routine julday366_day returns the Julian Day Number that begins at noon of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables.

  !>        The zeroth Julian Day is 01.01.-4712.

  !     CALLING SEQUENCE
  !         julian = julday366_day(dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[in] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[in] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: julian  !     Julian day

  !     EXAMPLE
  !         -> see example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Oct 2018
  elemental function julday366_day(dd,mm,yy)

    implicit none

    integer(i4), intent(in) :: dd, mm, yy
    integer(i4)             :: julday366_day
    integer(i4), parameter  :: year=366
    integer(i4),dimension(12),parameter :: months=(/31,29,31,30,31,30,31,31,30,31,30,31/)

    julday366_day = (yy+4712)*year + sum(months(1:mm-1)) + (dd-1)

  end function julday366_day


  ! ------------------------------------------------------------------

  !     NAME
  !         juldayJulian

  !     PURPOSE
  !>        \brief Julian day from day, month and year

  !>        \details In this routine juldayJulian returns the Julian Day Number that begins at noon of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables.

  !>        The zeroth Julian Day is 01.01.-4712 at noon, i.e. the 1st January 4713 BC 12:00:00 h.

  !>        Julian day definition starts at noon of the 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n

  !     CALLING SEQUENCE
  !         julian = juldayJulian(dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: dd"         Day in month of Julian day
  !>        \param[in] "integer(i4) :: mm"         Month in year of Julian day
  !>        \param[in] "integer(i4) :: yy"         Year of Julian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: julian &mdash;     Julian day

  !     RESTRICTIONS
  !>        \note Julian day definition starts at noon of the 1st January 4713 BC.\n
  !>        Here, the astronomical definition is used,
  !>        i.e. the year 1 BC (historic) is counted as 0 (astronomic), 2 BC is -1, etc.\n
  !>        This means that Julian day definition starts as 01.01.-4712 in astronomical units.\n
  !>        \n
  !>        juldayJulian and caldatJulian start at midnight of the 1st January 4713 BC.
  !>        So date2decJulian and juldayJulian as well as dec2dateJulian and caldatJulian are shifted by half a day.\n
  !>        Use date2decJulian with dec2dateJulian together for fractional Julian dates
  !>        and use juldayJulian with caldatJulian together for integer Julian days.

  !     EXAMPLE
  !         ! 2415021 is 01.01.1900
  !         ! 2440588 is 01.01.1970
  !         julian = juldayJulian(01,01,1990)
  !         -> see also example in test directory

  !     LITERATURE
  !         http://de.wikipedia.org/wiki/Julianisches_Datum
  !         which is different to the english Wiki
  !             http://en.wikipedia.org/wiki/Julian_day
  !         It is essentially the same as Numerical Recipes but uses astronomical instead of historical units.

  !     HISTORY
  !>        \author Written, Matthias Cuntz - modified julday from Numerical Recipes
  !>        \date Dec 2011
  !>        Modified Matthias Cuntz, May 2014 - changed to new algorithm with astronomical units
  !>                                            removed numerical recipes
  !>                 David Schaefer, Jan 2016 - renamed procedure
  ELEMENTAL FUNCTION juldayJulian(dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: dd, mm, yy
    INTEGER(i4) :: juldayJulian

    INTEGER(i4), PARAMETER :: IGREG2 = 15 + 31*(10+12*1582)
    INTEGER(i4), PARAMETER :: IGREG1 =  4 + 31*(10+12*1582)
    INTEGER(i8) :: jd, jm, jy
    INTEGER(i8) :: A, B

    if (mm > 2) then
       jm = int(mm, i8)
       jy = int(yy, i8)
    else
       jm = int(mm+12, i8)
       jy = int(yy-1, i8)
    endif

    jd = int(dd, i8)

    if (dd+31*(mm+12*yy) >= IGREG2) then ! gregorian
       A = jy/100_i8
       B = 2_i8 - A + A/4_i8
    else if (dd+31*(mm+12*yy) <= IGREG1) then ! julian
       B = 0_i8
       ! else
       !    stop 'No Gregorian dates between 04.10.1582 and 15.10.1582'
    endif

    ! add 0.5 to Wiki formula because formula was for fractional day
    ! juldayJulian = int(365.25_dp*real(jy+4716_i8,dp) + real(int(30.6001*real(jm+1_i8,dp),i8),dp) + real(jd+B,dp) - 1524.5_dp, i4)
    juldayJulian = int(365.25_dp*real(jy+4716_i8,dp) + real(int(30.6001*real(jm+1_i8,dp),i8),dp) &
         + real(jd+B,dp) - 1524.5_dp + 0.5_dp, i4)

  END FUNCTION juldayJulian


  ! ------------------------------------------------------------------

  !     NAME
  !         juldayLilian

  !     PURPOSE
  !>        \brief Lilian date from day, month and year

  !>        \details In this routine juldayLilian returns the Lilian Date that begins at midnight of the calendar
  !>        date specified by month mm, day dd, and year yy, all integer variables.

  !>        The first Lilian Day is 15.10.1582,
  !>        i.e. the day Pope Gregory XIII introduced the Gregorian calendar.

  !     CALLING SEQUENCE
  !         Lilian = juldayLilian(dd, mm, yy)

  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: dd"         Day in month of Lilian day
  !>        \param[in] "integer(i4) :: mm"         Month in year of Lilian day
  !>        \param[in] "integer(i4) :: yy"         Year of Lilian day

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return integer(i4) :: Lilian &mdash;     Lilian date

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! 115861 is 01.01.1900
  !         Lilian = juldayLilian(01,01,1990)
  !         -> see also example in test directory

  !     LITERATURE
  !         https://en.wikipedia.org/wiki/Lilian_date

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Dec 2016
  ELEMENTAL FUNCTION juldayLilian(dd,mm,yy)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: dd, mm, yy
    INTEGER(i4) :: juldayLilian

    INTEGER(i4), PARAMETER :: IGLIL = 2299160_i4

    juldayLilian = juldayJulian(dd,mm,yy) - IGLIL

  END FUNCTION juldayLilian

  ! ------------------------------------------------------------------

END MODULE mo_julian

!> \file mo_interpol.f90

!> \brief Module for linear and cubic B spline interpolation

!> \details The module provides one routine for linear interpolation of irregular grids
!>          and one routine for interpolation with cubic B splines.

!> \authors Matthias Cuntz
!> \date 2011-2014
MODULE mo_interpol

  ! This module provides a linear interpolation routine for irregular grids
  ! and a routine for cubic B spline interpolation

  ! Literature
  !   Inspired by interpolation routines of Numerical Recipes
  !       WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !           Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !           Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996
  !   and IDL routine interpol.pro Copyright (C) 1982-2004, Research Systems, Inc.
  !   The cubic B spline is a Fortran 90 version of the Fortran 77 version of John Burkardt.

  ! Written  Mar 2011, Matthias Cuntz - only linear interpolation
  ! Modified Jul 2013, Matthias Cuntz - added cubic B splines
  !          May 2014, Matthias Cuntz - removed numerical recipes

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! It is NOT released under the GNU Lesser General Public License, yet.

  ! If you use this routine, please contact Matthias Cuntz.

  ! Copyright 2011-2014 Matthias Cuntz

  USE mo_kind,  ONLY: i4, sp, dp
  USE mo_utils, only: le, locate

  Implicit NONE

  PUBLIC :: interpol        ! Linearly interpolate vectors with an irregular grid
  PUBLIC :: spline_b        ! Cubic B spline interpolation

  ! ------------------------------------------------------------------

  !     NAME
  !         interpol

  !     PURPOSE
  !>        \brief Linear interpolation on irregular grid.

  !>        \details Linearly interpolate vectors with an irregular grid.
  !>           yout = interpol(yin, xin, xout)
  !>           be z = xin(i) then
  !>           yout(i) = yin(i) + (z - floor(z)) * (yin(i+1) - yin(i))

  !     CALLING SEQUENCE
  !         yout = interpol(yin, xin, xout)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: yin(:)"      1D-array with input values to interpolate
  !>        \param[in] "real(sp/dp) :: xin(:)"      1D-array with position of yin values
  !>        \param[in] "real(sp/dp) :: xout(:)"     1D-array with wanted positions of interpolated values

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
  !>       \return     real(sp/dp) :: yout &mdash; interpolated values at xout

  !     RESTRICTIONS
  !>       \note Index arrays are not working with all compilers. Compile with -DABSOFT if problems occur.

  !     EXAMPLE
  !         ! Fill masked values of vector w with interpolated values
  !         w = interpol(pack(w,wmaske), pack(date,wmaske), date)
  !         -> see also example in test directory

  !     LITERATURE
  !       WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !           Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !           Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996
  !       IDL routine interpol.pro Copyright (C) 1982-2004, Research Systems, Inc.

  !     HISTORY
  !>        \author Written,  Matthias Cuntz
  !>        \date Nov 2011
  INTERFACE interpol
     MODULE PROCEDURE interpol_sp, interpol_dp
  END INTERFACE interpol

  ! ------------------------------------------------------------------

  !     NAME
  !         spline_b

  !     PURPOSE
  !>        \brief Evaluates a cubic B spline approximant.

  !>        \details The cubic B spline will approximate the data, but is not
  !>                 designed to interpolate it.
  !>
  !>                 In effect, two "phantom" data values are appended to the data,
  !>                 so that the spline will interpolate the first and last data values.

  !     CALLING SEQUENCE
  !         yout = spline_b(yin, xin, xout)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: yin(:)"      1D-array with input values to interpolate
  !>        \param[in] "real(sp/dp) :: xin(:)"      1D-array with position of yin values
  !>        \param[in] "real(sp/dp) :: xout(:)"     1D-array with wanted positions of interpolated values

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
  !>       \return     real(sp/dp) :: yout &mdash; spline-interpolated values at xout

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! Fill masked values of vector w with interpolated values
  !         w = spline_b(pack(w,wmaske), pack(date,wmaske), date)
  !         -> see also example in test directory

  !     LITERATURE
  !         Carl deBoor, A Practical Guide to Splines, Springer, 2001, ISBN: 0387953663.

  !     HISTORY
  !>        \author Written, John Burkardt (F77, Feb 2004). Modified Matthias Cuntz (F90, Jul 2013)
  !>        \date 1999-2013
  INTERFACE spline_b
     MODULE PROCEDURE spline_b_dp, spline_b_sp
  END INTERFACE spline_b

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

  ! For splines
  INTERFACE bracket
     MODULE PROCEDURE bracket_dp, bracket_sp
  END INTERFACE bracket

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  ! BRACKET searches a sorted vector for successive brackets of a value.
  !
  !  Discussion:
  !    A vector is an array of double precision real values.
  !    If the values in the vector are thought of as defining intervals
  !    on the real line, then this routine searches for the interval
  !    nearest to or containing the given value.
  !
  !  Licensing:
  !    This code is distributed under the GNU LGPL license.
  !
  !  Parameters:
  !    Input, real(dp) X(N), an array sorted into ascending order.
  !    Input, real(dp) XVAL(M), values to be bracketed.
  !    Output, integer(i4) LEFT(M), RIGHT(M), the results of the search.
  !    Either:
  !      XVAL < X(1), when LEFT = 1, RIGHT = 2;
  !      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
  !    or
  !      X(LEFT) <= XVAL <= X(RIGHT).

  !     HISTORY
  !        \author Written, John Burkardt (F77, Apr 1999). Modified Matthias Cuntz (F90, Jul 2013)
  !        \date 1999-2013
  subroutine bracket_dp(x, xval, left, right)

    implicit none

    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(in) :: xval
    integer(i4), dimension(size(xval)), intent(out) :: left
    integer(i4), dimension(size(xval)), intent(out) :: right

    integer(i4) :: n, m
    integer(i4) :: i, j
    logical :: istheend

    n = size(x)
    m = size(xval)
    do j=1, m
       istheend = .true.
       do i = 2, n - 1
          if ( xval(j) < x(i) ) then
             left(j)  = i - 1
             right(j) = i
             istheend = .false.
             exit
          end if
       end do
       if (istheend) then
          left(j)  = n - 1
          right(j) = n
       endif
    end do

    return

  end subroutine bracket_dp


  subroutine bracket_sp(x, xval, left, right)

    implicit none

    real(sp), dimension(:), intent(in) :: x
    real(sp), dimension(:), intent(in) :: xval
    integer(i4), dimension(size(xval)), intent(out) :: left
    integer(i4), dimension(size(xval)), intent(out) :: right

    integer(i4) :: n, m
    integer(i4) :: i, j
    logical :: istheend

    n = size(x)
    m = size(xval)
    do j=1, m
       istheend = .true.
       do i = 2, n - 1
          if ( xval(j) < x(i) ) then
             left(j)  = i - 1
             right(j) = i
             istheend = .false.
             exit
          end if
       end do
       if (istheend) then
          left(j)  = n - 1
          right(j) = n
       endif
    end do

    return

  end subroutine bracket_sp

  ! ------------------------------------------------------------------

  FUNCTION interpol_dp(v, x, u)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:), INTENT(IN) :: v
    REAL(dp),    DIMENSION(:), INTENT(IN) :: x
    REAL(dp),    DIMENSION(:), INTENT(IN) :: u
    REAL(dp),    DIMENSION(size(u))       :: interpol_dp

    INTEGER(i4), DIMENSION(size(u)) :: s
#ifdef ABSOFT
    REAL(dp)                        :: ums
#else
    REAL(dp),    DIMENSION(size(u)) :: ums
#endif
    INTEGER(i4) :: m
#ifdef ABSOFT
    INTEGER(i4) :: n
    INTEGER(i4) :: i
#endif

    m = size(v)
    if (size(x) /= m) stop 'Error interpol_dp: size(v) /= size(x)'

    s = min(max(locate(x, u), 1_i4), m-1_i4)  ! Subscript intervals

#ifdef ABSOFT
    n = size(u)
    do i=1, n
       ums = u(i) - x(s(i))
       if (abs(ums) < epsilon(1.0_dp)) ums = 0.0_dp
       interpol_dp(i) = v(s(i)) + ums*(v(s(i)+1)-v(s(i)))/(x(s(i)+1)-x(s(i)))
    end do
#else
    ums = (u-x(s)) ! distance from point before
    where (abs(ums) < epsilon(1.0_dp)) ums = 0.0_dp ! for numerical stability
    interpol_dp = v(s) + ums*(v(s+1)-v(s))/(x(s+1)-x(s))
#endif

  END FUNCTION interpol_dp


  FUNCTION interpol_sp(v, x, u)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:), INTENT(IN) :: v
    REAL(sp),    DIMENSION(:), INTENT(IN) :: x
    REAL(sp),    DIMENSION(:), INTENT(IN) :: u
    REAL(sp),    DIMENSION(size(u))       :: interpol_sp

    INTEGER(i4), DIMENSION(size(u)) :: s
#ifdef ABSOFT
    REAL(sp)                        :: ums
#else
    REAL(sp),    DIMENSION(size(u)) :: ums
#endif
    INTEGER(i4) :: m
#ifdef ABSOFT
    INTEGER(i4) :: n
    INTEGER(i4) :: i
#endif

    m = size(v)
    if (size(x) /= m) stop 'Error interpol_sp: size(v) /= size(x)'

    s = min(max(locate(x, u), 1_i4), m-1_i4)  ! Subscript intervals

#ifdef ABSOFT
    n = size(u)
    do i=1, n
       ums = u(i) - x(s(i))
       if (abs(ums) < epsilon(1.0_sp)) ums = 0.0_sp
       interpol_sp(i) = v(s(i)) + ums*(v(s(i)+1)-v(s(i)))/(x(s(i)+1)-x(s(i)))
    end do
#else
    ums = (u-x(s)) ! distance from point before
    where (abs(ums) < epsilon(1.0_sp)) ums = 0.0_sp ! for numerical stability
    interpol_sp = v(s) + ums*(v(s+1)-v(s))/(x(s+1)-x(s))
#endif

  END FUNCTION interpol_sp

  ! ------------------------------------------------------------------

  !    The cubic B spline will approximate the data, but is not
  !    designed to interpolate it.
  !
  !    In effect, two "phantom" data values are appended to the data,
  !    so that the spline will interpolate the first and last data values.
  function spline_b_dp(ydata, tdata, tval)

    implicit none

    real(dp), dimension(:), intent(in) :: ydata
    real(dp), dimension(:), intent(in) :: tdata
    real(dp), dimension(:), intent(in) :: tval
    real(dp), dimension(size(tval))    :: spline_b_dp
    
    integer(i4) :: ndata
    real(dp),    dimension(size(tval)) :: bval
    integer(i4), dimension(size(tval)) :: left
    integer(i4), dimension(size(tval)) :: right
    real(dp),    dimension(size(tval)) :: u

    ndata = size(tdata)
    !
    !  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
    !
    call bracket(tdata, tval, left, right)
    !
    !  Evaluate the 5 nonzero B spline basis functions in the interval,
    !  weighted by their corresponding data values.
    !
    u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
    spline_b_dp = 0.0_dp
    !
    !  B function associated with node LEFT - 1, (or "phantom node"),
    !  evaluated in its 4th interval.
    !
    bval = ( ( (     - 1.0_dp   &
         * u + 3.0_dp ) &
         * u - 3.0_dp ) &
         * u + 1.0_dp ) / 6.0_dp

    where ( 0 < left-1 )
       spline_b_dp = spline_b_dp + ydata(left-1) * bval
    elsewhere
       spline_b_dp = spline_b_dp + ( 2.0_dp * ydata(1) - ydata(2) ) * bval
    end where
    !
    !  B function associated with node LEFT,
    !  evaluated in its third interval.
    !
    bval = ( ( (       3.0_dp   &
         * u - 6.0_dp ) &
         * u + 0.0_dp ) &
         * u + 4.0_dp ) / 6.0_dp

    spline_b_dp = spline_b_dp + ydata(left) * bval
    !
    !  B function associated with node RIGHT,
    !  evaluated in its second interval.
    !
    bval = ( ( (     - 3.0_dp   &
         * u + 3.0_dp ) &
         * u + 3.0_dp ) &
         * u + 1.0_dp ) / 6.0_dp

    spline_b_dp = spline_b_dp + ydata(right) * bval
    !
    !  B function associated with node RIGHT+1, (or "phantom node"),
    !  evaluated in its first interval.
    !
    bval = u**3 / 6.0_dp

    where ( right+1 <= ndata )
       spline_b_dp = spline_b_dp + ydata(right+1) * bval
    elsewhere
       spline_b_dp = spline_b_dp + ( 2.0_dp * ydata(ndata) - ydata(ndata-1) ) * bval
    end where

    return

  end function spline_b_dp


  function spline_b_sp(ydata, tdata, tval)

    implicit none

    real(sp), dimension(:), intent(in) :: ydata
    real(sp), dimension(:), intent(in) :: tdata
    real(sp), dimension(:), intent(in) :: tval
    real(sp), dimension(size(tval))    :: spline_b_sp
    
    integer(i4) :: ndata
    real(sp),    dimension(size(tval)) :: bval
    integer(i4), dimension(size(tval)) :: left
    integer(i4), dimension(size(tval)) :: right
    real(sp),    dimension(size(tval)) :: u

    ndata = size(tdata)
    !
    !  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
    !
    call bracket(tdata, tval, left, right)
    !
    !  Evaluate the 5 nonzero B spline basis functions in the interval,
    !  weighted by their corresponding data values.
    !
    u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
    spline_b_sp = 0.0_sp
    !
    !  B function associated with node LEFT - 1, (or "phantom node"),
    !  evaluated in its 4th interval.
    !
    bval = ( ( (     - 1.0_sp   &
         * u + 3.0_sp ) &
         * u - 3.0_sp ) &
         * u + 1.0_sp ) / 6.0_sp

    where ( 0 < left-1 )
       spline_b_sp = spline_b_sp + ydata(left-1) * bval
    elsewhere
       spline_b_sp = spline_b_sp + ( 2.0_sp * ydata(1) - ydata(2) ) * bval
    end where
    !
    !  B function associated with node LEFT,
    !  evaluated in its third interval.
    !
    bval = ( ( (       3.0_sp   &
         * u - 6.0_sp ) &
         * u + 0.0_sp ) &
         * u + 4.0_sp ) / 6.0_sp

    spline_b_sp = spline_b_sp + ydata(left) * bval
    !
    !  B function associated with node RIGHT,
    !  evaluated in its second interval.
    !
    bval = ( ( (     - 3.0_sp   &
         * u + 3.0_sp ) &
         * u + 3.0_sp ) &
         * u + 1.0_sp ) / 6.0_sp

    spline_b_sp = spline_b_sp + ydata(right) * bval
    !
    !  B function associated with node RIGHT+1, (or "phantom node"),
    !  evaluated in its first interval.
    !
    bval = u**3 / 6.0_sp

    where ( right+1 <= ndata )
       spline_b_sp = spline_b_sp + ydata(right+1) * bval
    elsewhere
       spline_b_sp = spline_b_sp + ( 2.0_sp * ydata(ndata) - ydata(ndata-1) ) * bval
    end where

    return

  end function spline_b_sp

  ! ------------------------------------------------------------------

END MODULE mo_interpol

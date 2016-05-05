!> \file mo_utils.f90

!> \brief General utilities for the CHS library

!> \details This module provides general utilities such as comparisons of two reals.

!> \authors Matthias Cuntz, Juliane Mai
!> \date Feb 2014
MODULE mo_utils

  ! Written  Matthias Cuntz, Juliane Mai, Feb 2014
  ! Modified Matthias Cuntz, Juliane Mai, Feb 2014 - equal, notequal
  !          Matthias Cuntz,              May 2014 - swap
  !          Matthias Cuntz,              May 2014 - is_finite, is_nan, is_normal, special_value

  ! License
  ! -------
  ! This file is part of the JAMS Fortran library.

  ! The JAMS Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The JAMS Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the JAMS Fortran library (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2014 Matthias Cuntz, Juliane Mai

  USE mo_kind,         only: sp, dp, i4, spc, dpc

  IMPLICIT NONE

  PUBLIC :: equal         ! a == b, a .eq. b
  PUBLIC :: greaterequal  ! a >= b, a .ge. b
  PUBLIC :: lesserequal   ! a <= b, a .le. b
  PUBLIC :: notequal      ! a /= b, a .ne. b
  PUBLIC :: eq            ! a == b, a .eq. b
  PUBLIC :: ge            ! a >= b, a .ge. b
  PUBLIC :: le            ! a <= b, a .le. b
  PUBLIC :: ne            ! a /= b, a .ne. b
  PUBLIC :: is_finite     ! .true. if not IEEE Inf and not IEEE NaN
  PUBLIC :: is_nan        ! .true. if IEEE NaN
  PUBLIC :: is_normal     ! .true. if not IEEE Inf and not IEEE NaN
  PUBLIC :: locate        ! Find closest values in a monotonic series
  PUBLIC :: swap          ! swaps arrays or elements of an array
  PUBLIC :: special_value ! Special IEEE values

  ! ------------------------------------------------------------------

  !     NAME
  !         equal / notequal / greaterequal / lesserequal

  !     PURPOSE
  !         Elemental function returning .true. or .false. depending if the reals are equal or not.
  !
  !>        \brief Comparison of real values.
  !
  !>        \details Compares two reals if they are numerically equal or not, i.e.
  !>        equal: \f[ |\frac{a-b}{b}| < \epsilon \f]
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: a"        First number to compare
  !>        \param[in] "real(sp/dp) :: b"        Second number to compare
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>       \return     real(sp/dp) :: equal &mdash; \f$ a == b \f$ logically true or false
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         vec1 = (/ 1., 2., 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 1., 3., -999., 10., 6. /)
  !         isequal = equal(vec1, vec2)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \authors Matthias Cuntz, Juliane Mai
  !>        \date Feb 2014
  !         Modified, Matthias Cuntz, Juliane Mai, Feb 2014 - sp, dp
  INTERFACE equal
     MODULE PROCEDURE equal_sp, equal_dp
  END INTERFACE equal

  INTERFACE notequal
     MODULE PROCEDURE notequal_sp, notequal_dp
  END INTERFACE notequal

  INTERFACE greaterequal
     MODULE PROCEDURE greaterequal_sp, greaterequal_dp
  END INTERFACE greaterequal

  INTERFACE lesserequal
     MODULE PROCEDURE lesserequal_sp, lesserequal_dp
  END INTERFACE lesserequal

  INTERFACE eq
     MODULE PROCEDURE equal_sp, equal_dp
  END INTERFACE eq

  INTERFACE ne
     MODULE PROCEDURE notequal_sp, notequal_dp
  END INTERFACE ne

  INTERFACE ge
     MODULE PROCEDURE greaterequal_sp, greaterequal_dp
  END INTERFACE ge

  INTERFACE le
     MODULE PROCEDURE lesserequal_sp, lesserequal_dp
  END INTERFACE le


  ! ------------------------------------------------------------------

  !     NAME
  !         is_finite / is_nan / is_normal

  !     PURPOSE
  !         Elemental inquiry functions returning .true. if the argument has a value
  !         implied by the name of the function.
  !
  !>        \brief .true. if not IEEE Inf, IEEE NaN, nor IEEE Inf nor IEEE NaN, respectively.
  !
  !>        \details Checks for IEEE Inf and IEEE NaN, i.e. Infinity and Not-a-Number.\n
  !>                 Wraps to functions of the intrinsic module ieee_arithmetic
  !>                 but gives alternatives for gfortran, which does not provide ieee_arithmetic.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x"        Number to check
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>       \return logical :: is_finite/is_nan/is_normal &mdash; \f$ a /= Inf, a == NaN, a /= Inf and a == NaN \f$,
  !>                                                             logically true or false
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         vec1 = (/ NaN, 2., 3., Inf, 5., 6. /)
  !         is_finite = equal(vec1)
  !         is_nan    = equal(vec1)
  !         is_normal = equal(vec1)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \authors Matthias Cuntz
  !>        \date Mar 2015
  INTERFACE is_finite
     MODULE PROCEDURE is_finite_sp, is_finite_dp
  END INTERFACE is_finite  

  INTERFACE is_nan
     MODULE PROCEDURE is_nan_sp, is_nan_dp
  END INTERFACE is_nan
  
  INTERFACE is_normal
     MODULE PROCEDURE is_normal_sp, is_normal_dp
  END INTERFACE is_normal  


  ! ------------------------------------------------------------------

  !     NAME
  !         locate

  !     PURPOSE
  !         Find closest values in a monotonic series
  !
  !>         \brief Find closest values in a monotonic series, returns the indexes.
  !
  !>         \details Given an array x(1:n), and given a value y,
  !>         returns a value j such that y is between
  !>         x(j) and x(j+1).\n
  !
  !>         x must be monotonically increasing.\n
  !>         j=0 or j=N is returned to indicate that x is out of range.
  !
  !     INTENT(IN)
  !>        \param[in] "real(dp/sp)    :: x(:)"            Sorted array
  !>        \param[in] "real(dp/sp)    :: y[(:)]"          Value(s) of which the closest match in x(:) is wanted
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>       \return     integer(i4) :: index[(:)] &mdash; index(es) of x so that y is between x(index) and x(index+1)
  !
  !     RESTRICTIONS
  !>       \note x must be monotonically increasing.\n
  !
  !     EXAMPLE
  !         x = (/ 1., 2., 3., -999., 5., 6. /)
  !         y = (/ 1.1, 5.6 /)
  !         ii = locate(x, y)
  !         -> ii == (/ 1, 5 /)
  !         y = 1.1
  !         ii = locate(x, y)
  !         -> ii == 1
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date May 2014
  INTERFACE locate
     MODULE PROCEDURE locate_0d_dp, locate_0d_sp, locate_1d_dp, locate_1d_sp
  END INTERFACE locate

  
  ! ------------------------------------------------------------------

  !     NAME
  !         swap

  !     PURPOSE
  !         Swap two values/arrays or two elements in 1D-array.
  !
  !>        \brief Swap two values or exchange two elements in array.
  !
  !>        \details Swaps either two entities, i.e. scalars, vectors, matrices,
  !>                 or exchanges two elements in a vector.
  !>                 If an optinal mask is given, the only elements with mask==.true. will be exchanged.\n
  !>                 The call is either \n
  !>                   call swap(x, y, mask=mask) \n
  !>                 or \n
  !>                   call swap(vec, i, j, mask=mask)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4)    :: i"   Index of first element to be swapped with second [case swap(vec,i,j)]
  !>        \param[in] "integer(i4)    :: j"   Index of second element to be swapped with first [case swap(vec,i,j)]
  !
  !     INTENT(INOUT)
  !>        \param[inout] "real(sp/dp)/integer(i4)/complex(spc/dpc) :: x[(:,...)]"
  !>                       First scalar or array to swap with second [case swap(x,y)]
  !>        \param[inout] "real(sp/dp)/integer(i4)/complex(spc/dpc) :: y[(:[,:])]"
  !>                       Second scalar or array to swap with first [case swap(x,y)]
  !>
  !>        \param[inout] "real(sp/dp)/integer(i4)/complex(spc/dpc) :: x(:)"
  !>                       Vector of which to elements are swapped [case swap(vec,i,j)]
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: mask[(:,...)]" scalar or array logical mask\n
  !>                                                 If present, only those elements will be swapped
  !>                                                 where mask==.true.
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         vec1 = (/ 1., 2., 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 1., 3., -999., 10., 6. /)
  !         call swap(vec1, vec2, mask=(vec==-999.))
  !         call swap(vec1, 1, 3)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date May 2014
  INTERFACE swap
     MODULE PROCEDURE &
          swap_xy_dp,       swap_xy_sp,       swap_xy_i4,       swap_xy_dpc,       swap_xy_spc, &
          swap_xy_mask_dp,  swap_xy_mask_sp,  swap_xy_mask_i4,  swap_xy_mask_dpc,  swap_xy_mask_spc, &
          swap_vec_dp,      swap_vec_sp,      swap_vec_i4,      swap_vec_dpc,      swap_vec_spc,&
          swap_vec_mask_dp, swap_vec_mask_sp, swap_vec_mask_i4, swap_vec_mask_dpc, swap_vec_mask_spc
  END INTERFACE swap

  
  ! ------------------------------------------------------------------

  !     NAME
  !         special_value

  !     PURPOSE
  !         Mimics the function ieee_value of the intrinsic module ieee_arithmetic.
  !
  !>        \brief Special IEEE values.
  !
  !>        \details Returns special IEEE values such as Infinity or Not-a-Number.\n
  !>                 Wraps to function ieee_value of the intrinsic module ieee_arithmetic
  !>                 but gives alternatives for gfortran, which does not provide ieee_arithmetic.\n
  !>                 Quiet and signaling NaN are the same in case of gfortran;\n
  !>                 also denormal values are the same as inf.
  !>
  !>                 Current special values are:\n
  !>                 IEEE_SIGNALING_NAN\n
  !>                 IEEE_QUIET_NAN\n
  !>                 IEEE_NEGATIVE_INF\n
  !>                 IEEE_POSITIVE_INF\n
  !>                 IEEE_NEGATIVE_DENORMAL\n
  !>                 IEEE_POSITIVE_DENORMAL\n
  !>                 IEEE_NEGATIVE_NORMAL\n
  !>                 IEEE_POSITIVE_NORMAL\n
  !>                 IEEE_NEGATIVE_ZERO\n
  !>                 IEEE_POSITIVE_ZERO
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x"         dummy for kind of output
  !>        \param[in] "character(le=*) :: name   ieee signal nanme
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>       \return real(sp/dp) :: special_value &mdash; IEEE special value\n
  !>                 IEEE_SIGNALING_NAN\n
  !>                 IEEE_QUIET_NAN (==IEEE_SIGNALING_NAN for gfortran)\n
  !>                 IEEE_NEGATIVE_INF\n
  !>                 IEEE_POSITIVE_INF\n
  !>                 IEEE_NEGATIVE_DENORMAL (==-0.0 for gfortran)\n
  !>                 IEEE_POSITIVE_DENORMAL (==0.0 for gfortran)\n
  !>                 IEEE_NEGATIVE_NORMAL (==-1.0 for gfortran)\n
  !>                 IEEE_POSITIVE_NORMAL (==1.0 for gfortran)\n
  !>                 IEEE_NEGATIVE_ZERO\n
  !>                 IEEE_POSITIVE_ZERO\n

  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         NaN = special_value(1.0, 'IEEE_QUIET_NAN')
  !         nan = special_value(1.0_dp, 'ieee_quiet_nan')
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \authors Matthias Cuntz
  !>        \date Mar 2015
  INTERFACE special_value
     MODULE PROCEDURE special_value_sp, special_value_dp
  END INTERFACE special_value

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         toupper_

  !     PURPOSE
  !         \brief Convert to upper case

  !         \details Convert all lower case letters in string to upper case letters.

  !     CALLING SEQUENCE
  !         up = toupper_(lower)

  !     INTENT(IN)
  !         \param[in] "character(len=*) :: lower"    String

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
  !         \return character(len=len_trim(lower)) :: up  &mdash;  String where all lowercase in input is converted to uppercase

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! Returns 'HALLO'
  !         up = toupper_('Hallo')
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         \author Matthias Cuntz - modified from Echam5, (C) MPI-MET, Hamburg, Germany
  !         \date Dec 2011

  FUNCTION toupper_(lower)

    IMPLICIT NONE

    CHARACTER(LEN=*)              ,INTENT(in) :: lower
    CHARACTER(LEN=LEN_TRIM(lower))            :: toupper_

    INTEGER            :: i
    INTEGER, PARAMETER :: idel = ICHAR('A')-ICHAR('a')

    DO i=1,LEN_TRIM(lower)
       IF (ICHAR(lower(i:i)) >= ICHAR('a') .AND. &
            ICHAR(lower(i:i)) <= ICHAR('z')) THEN
          toupper_(i:i) = CHAR( ICHAR(lower(i:i)) + idel )
       ELSE
          toupper_(i:i) = lower(i:i)
       END IF
    END DO

  END FUNCTION toupper_

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION equal_dp(a, b)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: a
    REAL(dp), INTENT(IN) :: b
    LOGICAL              :: equal_dp

    if ((epsilon(1.0_dp)*abs(b) - abs(a-b)) < 0.0_dp) then
       equal_dp = .false.
    else
       equal_dp = .true.
    endif

  END FUNCTION equal_dp


  ELEMENTAL PURE FUNCTION equal_sp(a, b)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: a
    REAL(sp), INTENT(IN) :: b
    LOGICAL              :: equal_sp

    if ((epsilon(1.0_sp)*abs(b) - abs(a-b)) < 0.0_sp) then
       equal_sp = .false.
    else
       equal_sp = .true.
    endif

  END FUNCTION equal_sp

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION greaterequal_dp(a, b)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: a
    REAL(dp), INTENT(IN) :: b
    LOGICAL              :: greaterequal_dp

    greaterequal_dp = .true.
    ! 1st part is /=, 2nd part is the a<b
    if ( ((epsilon(1.0_dp)*abs(b) - abs(a-b)) < 0.0_dp) .and. (a < b) ) greaterequal_dp = .false.

  END FUNCTION greaterequal_dp


  ELEMENTAL PURE FUNCTION greaterequal_sp(a, b)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: a
    REAL(sp), INTENT(IN) :: b
    LOGICAL              :: greaterequal_sp

    greaterequal_sp = .true.
    ! 1st part is /=, 2nd part is the a<b
    if ( ((epsilon(1.0_sp)*abs(b) - abs(a-b)) < 0.0_sp) .and. (a < b) ) greaterequal_sp = .false.

  END FUNCTION greaterequal_sp

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION lesserequal_dp(a, b)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: a
    REAL(dp), INTENT(IN) :: b
    LOGICAL              :: lesserequal_dp

    lesserequal_dp = .true.
    ! 1st part is /=, 2nd part is the a>b
    if ( ((epsilon(1.0_dp)*abs(b) - abs(a-b)) < 0.0_dp) .and. (a > b) ) lesserequal_dp = .false.

  END FUNCTION lesserequal_dp


  ELEMENTAL PURE FUNCTION lesserequal_sp(a, b)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: a
    REAL(sp), INTENT(IN) :: b
    LOGICAL              :: lesserequal_sp

    lesserequal_sp = .true.
    ! 1st part is /=, 2nd part is the a>b
    if ( ((epsilon(1.0_sp)*abs(b) - abs(a-b)) < 0.0_sp) .and. (a > b) ) lesserequal_sp = .false.

  END FUNCTION lesserequal_sp

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION notequal_dp(a, b)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: a
    REAL(dp), INTENT(IN) :: b
    LOGICAL              :: notequal_dp

    if ((epsilon(1.0_dp)*abs(b) - abs(a-b)) < 0.0_dp) then
       notequal_dp = .true.
    else
       notequal_dp = .false.
    endif

  END FUNCTION notequal_dp


  ELEMENTAL PURE FUNCTION notequal_sp(a, b)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: a
    REAL(sp), INTENT(IN) :: b
    LOGICAL              :: notequal_sp

    if ((epsilon(1.0_sp)*abs(b) - abs(a-b)) < 0.0_sp) then
       notequal_sp = .true.
    else
       notequal_sp = .false.
    endif

  END FUNCTION notequal_sp
  
  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION is_finite_dp(a)

#ifndef GFORTRAN
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#endif

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: a
    LOGICAL              :: is_finite_dp

#ifndef GFORTRAN
    is_finite_dp = ieee_is_finite(a)
#else
    is_finite_dp = (.not. ((a > huge(a)) .or. (a < -huge(a)))) .and. (.not. is_nan(a))
#endif

  END FUNCTION is_finite_dp

  ELEMENTAL PURE FUNCTION is_finite_sp(a)

#ifndef GFORTRAN
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#endif

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: a
    LOGICAL              :: is_finite_sp

#ifndef GFORTRAN
    is_finite_sp = ieee_is_finite(a)
#else
    is_finite_sp = (.not. ((a > huge(a)) .or. (a < -huge(a)))) .and. (.not. is_nan(a))
#endif

  END FUNCTION is_finite_sp


  ELEMENTAL PURE FUNCTION is_nan_dp(a)

#ifndef GFORTRAN
  use, intrinsic :: ieee_arithmetic, only: isnan => ieee_is_nan
#endif

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: a
    LOGICAL              :: is_nan_dp

    ! isnan introduced in gfortran rev 4.2
#ifdef GFORTRAN41
    is_nan_dp = a /= a
#else
    is_nan_dp = isnan(a)
#endif

  END FUNCTION is_nan_dp

  ELEMENTAL PURE FUNCTION is_nan_sp(a)

#ifndef GFORTRAN
  use, intrinsic :: ieee_arithmetic, only: isnan => ieee_is_nan
#endif

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: a
    LOGICAL              :: is_nan_sp

    ! isnan introduced in gfortran rev 4.2
#ifdef GFORTRAN41
    is_nan_sp = a /= a
#else
    is_nan_sp = isnan(a)
#endif

  END FUNCTION is_nan_sp


  ELEMENTAL PURE FUNCTION is_normal_dp(a)

#ifndef GFORTRAN
  use, intrinsic :: ieee_arithmetic, only: ieee_is_normal
#endif

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: a
    LOGICAL              :: is_normal_dp

#ifndef GFORTRAN
    is_normal_dp = ieee_is_normal(a)
#else
    is_normal_dp = is_finite(a)
#endif

  END FUNCTION is_normal_dp

  ELEMENTAL PURE FUNCTION is_normal_sp(a)

#ifndef GFORTRAN
  use, intrinsic :: ieee_arithmetic, only: ieee_is_normal
#endif

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: a
    LOGICAL              :: is_normal_sp

#ifndef GFORTRAN
    is_normal_sp = ieee_is_normal(a)
#else
    is_normal_sp = is_finite(a)
#endif

  END FUNCTION is_normal_sp

  ! ------------------------------------------------------------------

  ! Given an array x(1:N), and given a value y, returns a value j such that y is between
  !  x(j) and x(j+1). x must be monotonically increasing.
  !  j=0 or j=N is returned to indicate that x is out of range.

  FUNCTION locate_0d_dp(x,y)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: x
    REAL(dp),               INTENT(IN) :: y
    INTEGER(i4)                        :: locate_0d_dp

    INTEGER(i4), dimension(1) :: c

    c = minloc(abs(x-y))
    if (le(x(c(1)),y)) then
       locate_0d_dp = c(1)
    else
       locate_0d_dp = c(1)-1
    endif

  END FUNCTION locate_0d_dp

  FUNCTION locate_0d_sp(x,y)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN) :: x
    REAL(sp),               INTENT(IN) :: y
    INTEGER(i4)                        :: locate_0d_sp

    INTEGER(i4), dimension(1) :: c

    c = minloc(abs(x-y))
    if (le(x(c(1)),y)) then
       locate_0d_sp = c(1)
    else
       locate_0d_sp = c(1)-1
    endif

  END FUNCTION locate_0d_sp

  FUNCTION locate_1d_dp(x,y)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: x
    REAL(dp), DIMENSION(:), INTENT(IN) :: y
    INTEGER(i4), DIMENSION(:), allocatable :: locate_1d_dp

    INTEGER(i4) :: ny, i
    INTEGER(i4), dimension(1) :: c


    ny = size(y)
    if (.not. allocated(locate_1d_dp)) allocate(locate_1d_dp(ny))

    do i=1, ny
       c = minloc(abs(x-y(i)))
       if (le(x(c(1)),y(i))) then
          locate_1d_dp(i) = c(1)
       else
          locate_1d_dp(i) = c(1)-1
       endif
    end do

  END FUNCTION locate_1d_dp

  FUNCTION locate_1d_sp(x,y)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN) :: x
    REAL(sp), DIMENSION(:), INTENT(IN) :: y
    INTEGER(i4), DIMENSION(:), allocatable :: locate_1d_sp

    INTEGER(i4) :: ny, i
    INTEGER(i4), dimension(1) :: c


    ny = size(y)
    if (.not. allocated(locate_1d_sp)) allocate(locate_1d_sp(ny))

    do i=1, ny
       c = minloc(abs(x-y(i)))
       if (le(x(c(1)),y(i))) then
          locate_1d_sp(i) = c(1)
       else
          locate_1d_sp(i) = c(1)-1
       endif
    end do

  END FUNCTION locate_1d_sp

  ! ------------------------------------------------------------------

  elemental pure subroutine swap_xy_dp(x, y)

    implicit none

    real(dp),          intent(inout) :: x
    real(dp),          intent(inout) :: y

    real(dp) :: z

    z = x
    x = y
    y = z

  end subroutine swap_xy_dp

  elemental pure subroutine swap_xy_sp(x, y)

    implicit none

    real(sp),          intent(inout) :: x
    real(sp),          intent(inout) :: y

    real(sp) :: z

    z = x
    x = y
    y = z

  end subroutine swap_xy_sp

  elemental pure subroutine swap_xy_i4(x, y)

    implicit none

    integer(i4),       intent(inout) :: x
    integer(i4),       intent(inout) :: y

    integer(i4) :: z

    z = x
    x = y
    y = z

  end subroutine swap_xy_i4

  elemental pure subroutine swap_xy_dpc(x, y)

    implicit none

    complex(dpc),          intent(inout) :: x
    complex(dpc),          intent(inout) :: y

    complex(dpc) :: z

    z = x
    x = y
    y = z

  end subroutine swap_xy_dpc

  elemental pure subroutine swap_xy_spc(x, y)

    implicit none

    complex(spc),          intent(inout) :: x
    complex(spc),          intent(inout) :: y

    complex(spc) :: z

    z = x
    x = y
    y = z

  end subroutine swap_xy_spc

  elemental pure subroutine swap_xy_mask_dp(x, y, mask)

    implicit none

    real(dp), intent(inout) :: x
    real(dp), intent(inout) :: y
    logical,  intent(in)    :: mask

    real(dp) :: z

    if (mask) then
       z = x
       x = y
       y = z
    endif

  end subroutine swap_xy_mask_dp

  elemental pure subroutine swap_xy_mask_sp(x, y, mask)

    implicit none

    real(sp), intent(inout) :: x
    real(sp), intent(inout) :: y
    logical,  intent(in)    :: mask

    real(sp) :: z

    if (mask) then
       z = x
       x = y
       y = z
    endif

  end subroutine swap_xy_mask_sp

  elemental pure subroutine swap_xy_mask_i4(x, y, mask)

    implicit none

    integer(i4), intent(inout) :: x
    integer(i4), intent(inout) :: y
    logical,     intent(in)    :: mask

    integer(i4) :: z

    if (mask) then
       z = x
       x = y
       y = z
    endif

  end subroutine swap_xy_mask_i4

  elemental pure subroutine swap_xy_mask_dpc(x, y, mask)

    implicit none

    complex(dpc), intent(inout) :: x
    complex(dpc), intent(inout) :: y
    logical,  intent(in)    :: mask

    complex(dpc) :: z

    if (mask) then
       z = x
       x = y
       y = z
    endif

  end subroutine swap_xy_mask_dpc

  elemental pure subroutine swap_xy_mask_spc(x, y, mask)

    implicit none

    complex(spc), intent(inout) :: x
    complex(spc), intent(inout) :: y
    logical,  intent(in)    :: mask

    complex(spc) :: z

    if (mask) then
       z = x
       x = y
       y = z
    endif

  end subroutine swap_xy_mask_spc


  subroutine swap_vec_dp(x,i1,i2)

    implicit none

    real(dp),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2

    real(dp) :: z

    z     = x(i1)
    x(i1) = x(i2)
    x(i2) = z

  end subroutine swap_vec_dp

  subroutine swap_vec_sp(x,i1,i2)

    implicit none

    real(sp),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2

    real(sp) :: z

    z     = x(i1)
    x(i1) = x(i2)
    x(i2) = z

  end subroutine swap_vec_sp

  subroutine swap_vec_i4(x,i1,i2)

    implicit none

    integer(i4),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2

    integer(i4) :: z

    z     = x(i1)
    x(i1) = x(i2)
    x(i2) = z

  end subroutine swap_vec_i4

  subroutine swap_vec_dpc(x,i1,i2)

    implicit none

    complex(dpc),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2

    complex(dpc) :: z

    z     = x(i1)
    x(i1) = x(i2)
    x(i2) = z

  end subroutine swap_vec_dpc

  subroutine swap_vec_spc(x,i1,i2)

    implicit none

    complex(spc),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2

    complex(spc) :: z

    z     = x(i1)
    x(i1) = x(i2)
    x(i2) = z

  end subroutine swap_vec_spc

  subroutine swap_vec_mask_dp(x, i1, i2, mask)

    implicit none

    real(dp),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2
    logical,                   intent(in)    :: mask

    real(dp) :: z

    if (mask) then
       z     = x(i1)
       x(i1) = x(i2)
       x(i2) = z
    endif

  end subroutine swap_vec_mask_dp

  subroutine swap_vec_mask_sp(x, i1, i2, mask)

    implicit none

    real(sp),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2
    logical,                   intent(in)    :: mask

    real(sp) :: z

    if (mask) then
       z     = x(i1)
       x(i1) = x(i2)
       x(i2) = z
    endif

  end subroutine swap_vec_mask_sp

  subroutine swap_vec_mask_i4(x, i1, i2, mask)

    implicit none

    integer(i4),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2
    logical,                   intent(in)    :: mask

    integer(i4) :: z

    if (mask) then
       z     = x(i1)
       x(i1) = x(i2)
       x(i2) = z
    endif

  end subroutine swap_vec_mask_i4

  subroutine swap_vec_mask_dpc(x, i1, i2, mask)

    implicit none

    complex(dpc),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2
    logical,                   intent(in)    :: mask

    complex(dpc) :: z

    if (mask) then
       z     = x(i1)
       x(i1) = x(i2)
       x(i2) = z
    endif

  end subroutine swap_vec_mask_dpc

  subroutine swap_vec_mask_spc(x, i1, i2, mask)

    implicit none

    complex(spc),    dimension(:), intent(inout) :: x
    integer(i4),               intent(in)    :: i1
    integer(i4),               intent(in)    :: i2
    logical,                   intent(in)    :: mask

    complex(spc) :: z

    if (mask) then
       z     = x(i1)
       x(i1) = x(i2)
       x(i2) = z
    endif

  end subroutine swap_vec_mask_spc

  ! ------------------------------------------------------------------

  function special_value_dp(x, ieee)

#ifndef GFORTRAN
    use, intrinsic :: ieee_arithmetic, only: ieee_value, &
         IEEE_SIGNALING_NAN, &
         IEEE_QUIET_NAN, &
         IEEE_NEGATIVE_INF, &
         IEEE_POSITIVE_INF, &
         IEEE_NEGATIVE_DENORMAL, &
         IEEE_POSITIVE_DENORMAL, &
         IEEE_NEGATIVE_NORMAL, &
         IEEE_POSITIVE_NORMAL, &
         IEEE_NEGATIVE_ZERO, &
         IEEE_POSITIVE_ZERO
#endif

    implicit none
    
    real(dp),         intent(in) :: x
    character(len=*), intent(in) :: ieee
    real(dp)                     :: special_value_dp

    ! local
    character(len=21) :: ieee_up
#ifdef GFORTRAN
    real(dp) :: tmp
#endif

    ieee_up = toupper_(ieee)
#ifndef GFORTRAN
    select case(trim(ieee_up))
    case('IEEE_SIGNALING_NAN')
       special_value_dp = ieee_value(x, IEEE_SIGNALING_NAN)
    case('IEEE_QUIET_NAN')
       special_value_dp = ieee_value(x, IEEE_QUIET_NAN)
    case('IEEE_NEGATIVE_INF')
       special_value_dp = ieee_value(x, IEEE_NEGATIVE_INF)
    case('IEEE_POSITIVE_INF')
       special_value_dp = ieee_value(x, IEEE_POSITIVE_INF)
    case('IEEE_NEGATIVE_DENORMAL')
       special_value_dp = ieee_value(x, IEEE_NEGATIVE_DENORMAL)
    case('IEEE_POSITIVE_DENORMAL')
       special_value_dp = ieee_value(x, IEEE_POSITIVE_DENORMAL)
    case('IEEE_NEGATIVE_NORMAL')
       special_value_dp = ieee_value(x, IEEE_NEGATIVE_NORMAL)
    case('IEEE_POSITIVE_NORMAL')
       special_value_dp = ieee_value(x, IEEE_POSITIVE_NORMAL)
    case('IEEE_NEGATIVE_ZERO')
       special_value_dp = ieee_value(x, IEEE_NEGATIVE_ZERO)
    case('IEEE_POSITIVE_ZERO')
       special_value_dp = ieee_value(x, IEEE_POSITIVE_ZERO)
    case default
       special_value_dp = 0.0_dp
    end select
#else
    select case(ieee_up)
    case('IEEE_SIGNALING_NAN')
       tmp = 0.0_dp
       special_value_dp = tmp/tmp
    case('IEEE_QUIET_NAN')
       tmp = 0.0_dp
       special_value_dp = tmp/tmp
    case('IEEE_NEGATIVE_INF')
       tmp = huge(x)
       special_value_dp = -tmp*tmp
    case('IEEE_POSITIVE_INF')
       tmp = huge(x)
       special_value_dp = tmp*tmp
    case('IEEE_NEGATIVE_DENORMAL')
       special_value_dp = -0.0_dp
    case('IEEE_POSITIVE_DENORMAL')
       special_value_dp = 0.0_dp
    case('IEEE_NEGATIVE_NORMAL')
       special_value_dp = -1.0_dp
    case('IEEE_POSITIVE_NORMAL')
       special_value_dp = 1.0_dp
    case('IEEE_NEGATIVE_ZERO')
       special_value_dp = -0.0_dp
    case('IEEE_POSITIVE_ZERO')
       special_value_dp = 0.0_dp
    case default
       special_value_dp = 0.0_dp
    end select
#endif

  end function special_value_dp

  function special_value_sp(x, ieee)

#ifndef GFORTRAN
    use, intrinsic :: ieee_arithmetic, only: ieee_value, &
         IEEE_SIGNALING_NAN, &
         IEEE_QUIET_NAN, &
         IEEE_NEGATIVE_INF, &
         IEEE_POSITIVE_INF, &
         IEEE_NEGATIVE_DENORMAL, &
         IEEE_POSITIVE_DENORMAL, &
         IEEE_NEGATIVE_NORMAL, &
         IEEE_POSITIVE_NORMAL, &
         IEEE_NEGATIVE_ZERO, &
         IEEE_POSITIVE_ZERO
#endif

    implicit none
    
    real(sp),         intent(in) :: x
    character(len=*), intent(in) :: ieee
    real(sp)                     :: special_value_sp

    ! local
    character(len=21) :: ieee_up
#ifdef GFORTRAN
    real(sp) :: tmp
#endif

    ieee_up = toupper_(ieee)
#ifndef GFORTRAN
    select case(trim(ieee_up))
    case('IEEE_SIGNALING_NAN')
       special_value_sp = ieee_value(x, IEEE_SIGNALING_NAN)
    case('IEEE_QUIET_NAN')
       special_value_sp = ieee_value(x, IEEE_QUIET_NAN)
    case('IEEE_NEGATIVE_INF')
       special_value_sp = ieee_value(x, IEEE_NEGATIVE_INF)
    case('IEEE_POSITIVE_INF')
       special_value_sp = ieee_value(x, IEEE_POSITIVE_INF)
    case('IEEE_NEGATIVE_DENORMAL')
       special_value_sp = ieee_value(x, IEEE_NEGATIVE_DENORMAL)
    case('IEEE_POSITIVE_DENORMAL')
       special_value_sp = ieee_value(x, IEEE_POSITIVE_DENORMAL)
    case('IEEE_NEGATIVE_NORMAL')
       special_value_sp = ieee_value(x, IEEE_NEGATIVE_NORMAL)
    case('IEEE_POSITIVE_NORMAL')
       special_value_sp = ieee_value(x, IEEE_POSITIVE_NORMAL)
    case('IEEE_NEGATIVE_ZERO')
       special_value_sp = ieee_value(x, IEEE_NEGATIVE_ZERO)
    case('IEEE_POSITIVE_ZERO')
       special_value_sp = ieee_value(x, IEEE_POSITIVE_ZERO)
    case default
       special_value_sp = 0.0_sp
    end select
#else
    select case(ieee_up)
    case('IEEE_SIGNALING_NAN')
       tmp = 0.0_sp
       special_value_sp = tmp/tmp
    case('IEEE_QUIET_NAN')
       tmp = 0.0_sp
       special_value_sp = tmp/tmp
    case('IEEE_NEGATIVE_INF')
       tmp = huge(x)
       special_value_sp = -tmp*tmp
    case('IEEE_POSITIVE_INF')
       tmp = huge(x)
       special_value_sp = tmp*tmp
    case('IEEE_NEGATIVE_DENORMAL')
       special_value_sp = -0.0_sp
    case('IEEE_POSITIVE_DENORMAL')
       special_value_sp = 0.0_sp
    case('IEEE_NEGATIVE_NORMAL')
       special_value_sp = -1.0_sp
    case('IEEE_POSITIVE_NORMAL')
       special_value_sp = 1.0_sp
    case('IEEE_NEGATIVE_ZERO')
       special_value_sp = -0.0_sp
    case('IEEE_POSITIVE_ZERO')
       special_value_sp = 0.0_sp
    case default
       special_value_sp = 0.0_sp
    end select
#endif

  end function special_value_sp

END MODULE mo_utils

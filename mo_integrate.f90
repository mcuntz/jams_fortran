!> \file mo_integrate.f90

!> \brief Provides integration routines

!> \details This module provides routine for numerical integration such a Newton-Cotes formulas, etc.

!> \authors Matthias Cuntz
!> \date Mar 2013
MODULE mo_integrate

  ! Provide numerical integration.

  ! Written  Matthias Cuntz, Mar 2013

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
  ! along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2013 Matthias Cuntz

  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE

  PUBLIC :: int_regular ! Integrate regularily spaced data

  ! ------------------------------------------------------------------

  !     NAME
  !         int_regular

  !     PURPOSE
  !>        \brief Integrate regularily spaced data.

  !>        \details Integrates regularily spaced data with a 5-point Newton-Cotes formula:
  !>        \f[ \int y = \frac{2}{45} dx \sum_{i=5,n-4,4} 7 y_{i-4} + 32 y_{i-3} + 12 y_{i-2} + 32 y_{i-1} + 7 y_{i} \f]

  !>        dx=1 if not given.

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: dat(:)"        \f$ y_i \f$ 1D-array with y-values.

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>       \param[in] "real(sp/dp) :: dx"             x-spacing (default=1.)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return     real(sp/dp) :: integral &mdash; \f$ \int y \f$ integral of y values

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., 4., 5., 6., 7., 8., 9. /)
  !         m   = int_regular(vec)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date Mar 2013
  INTERFACE int_regular
     MODULE PROCEDURE int_regular_sp, int_regular_dp
  END INTERFACE int_regular

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  ! integrates tabulated function with equal distant points usign 5-point Newton-Cotes formula
  FUNCTION int_regular_dp(dat, dx)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)  :: dat
    REAL(dp),               OPTIONAL, INTENT(IN)  :: dx
    REAL(dp)                                      :: int_regular_dp

    INTEGER(i4) :: n, n0
    REAL(dp)    :: ddx

    if (size(dat,1) < 5) stop 'Error int_regular_dp: size(dat) < 5'

    if (present(dx)) then
       ddx = dx*2.0_dp/45.0_dp
    else
       ddx = 2.0_dp/45.0_dp
    endif

    !print*, 'dat ',dat
    !print*, 'ddx ',ddx
    !pause

    n0 = 5
    n  = size(dat,1)

    if (ddx .gt. 0.0_dp) then
       int_regular_dp = sum( &
            (7.0_dp*(dat(n0-4:n-4:4) + dat(n0:n:4)) + &
            32.0_dp*(dat(n0-3:n-3:4) + dat(n0-1:n-1:4)) + &
            12.0_dp*dat(n0-2:n-2:4)) )
       ! to avoid underflow issues
       if ( ddx .lt. 1.0_dp ) then
          if ( int_regular_dp .gt. tiny(1.0_dp)/ddx ) then
             int_regular_dp = ddx * int_regular_dp
          else 
             int_regular_dp = tiny(1.0_dp)
          end if
       else
          int_regular_dp = ddx * int_regular_dp
       end if
    else
       int_regular_dp = 0.0_dp
    end if

  END FUNCTION int_regular_dp

  FUNCTION int_regular_sp(dat, dx)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)  :: dat
    REAL(sp),               OPTIONAL, INTENT(IN)  :: dx
    REAL(sp)                                      :: int_regular_sp

    INTEGER(i4) :: n, n0
    REAL(sp)    :: ddx

    if (size(dat,1) < 5) stop 'Error int_regular_sp: size(dat) < 5'

    if (present(dx)) then
       ddx = dx*2.0_sp/45.0_sp
    else
       ddx = 2.0_sp/45.0_sp
    endif

    n0 = 5
    n  = size(dat,1)
    if (ddx .gt. 0.0_sp) then
       int_regular_sp = sum( &
            (7.0_sp*(dat(n0-4:n-4:4) + dat(n0:n:4)) + &
            32.0_sp*(dat(n0-3:n-3:4) + dat(n0-1:n-1:4)) + &
            12.0_sp*dat(n0-2:n-2:4)) )
       ! to avoid underflow issues
       if ( ddx .lt. 1.0_sp ) then
          if ( int_regular_sp .gt. tiny(1.0_sp)/ddx ) then
             int_regular_sp = ddx * int_regular_sp
          else 
             int_regular_sp = tiny(1.0_sp)
          end if
       else
          int_regular_sp = ddx * int_regular_sp
       end if
    else
       int_regular_sp = 0.0_sp
    end if

  END FUNCTION int_regular_sp

END MODULE mo_integrate

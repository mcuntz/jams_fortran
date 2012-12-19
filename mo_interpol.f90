MODULE mo_interpol

  ! This module provides a linear interpolation routine for irregular grids

  ! Literature
  !   Inspired by interpolation routines of Numerical Recipes
  !       WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !           Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !           Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996
  !   and IDL routine interpol.pro Copyright (C) 1982-2004, Research Systems, Inc.

  ! Written March 2011, Matthias Cuntz

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

  USE mo_kind, ONLY: i4, sp, dp

  Implicit NONE

  PUBLIC :: interpol        ! Linearly interpolate vectors with an irregular grid

  ! ------------------------------------------------------------------

  !     NAME
  !         interpol

  !     PURPOSE
  !         Linearly interpolate vectors with an irregular grid.
  !            yout = interpol(yin, xin, xout)
  !            be z = xin(i) then
  !            yout(i) = yin(i) + (z - floor(z)) * (yin(i+1) - yin(i))

  !     CALLING SEQUENCE
  !         yout = interpol(yin, xin, xout)

  !     INTENT(IN)
  !         real(sp/dp) :: yin(:)      1D-array with input values to interpolate
  !         real(sp/dp) :: xin(:)      1D-array with position of yin values
  !         real(sp/dp) :: xout(:)     1D-array with wanted positions of interpolated values

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: yout        yin interpolated at xout

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Index arrays are not working with all compilers. Compile with -DABSOFT if problems occur.

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
  !         Written,  Matthias Cuntz, Nov 2011
  INTERFACE interpol
     MODULE PROCEDURE interpol_sp, interpol_dp
  END INTERFACE interpol

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

  ! Private function, from numerical recipes
  INTERFACE locate
     MODULE PROCEDURE locate_0d_sp, locate_1d_sp, locate_0d_dp, locate_1d_dp
  END INTERFACE locate

  ! ------------------------------------------------------------------

CONTAINS

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

  ! From numerical recipes documentation
  ! Given an array xx(1:N), and given a value x, returns a value j such that x is between
  !  xx(j) and xx(j+1). xx must be monotonic, either increasing or decreasing. j=0 or
  ! j=N is returned to indicate that x is out of range.
  FUNCTION locate_0d_dp(xx,x)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: xx
    REAL(dp),               INTENT(IN) :: x
    INTEGER(i4)                        :: locate_0d_dp

    INTEGER(i4) :: n, jl, jm, ju
    LOGICAL     :: ascnd

    n = size(xx)
    ascnd = (xx(n) >= xx(1))
    jl = 0
    ju = n+1
    do
       if (ju-jl <= 1) exit
       jm = (ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       end if
    end do
    if (x == xx(1)) then
       locate_0d_dp = 1
    else if (x == xx(n)) then
       locate_0d_dp = n-1
    else
       locate_0d_dp = jl
    end if

  END FUNCTION locate_0d_dp


  FUNCTION locate_1d_dp(xx,x)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: xx
    REAL(dp), DIMENSION(:), INTENT(IN) :: x
    INTEGER(i4), DIMENSION(size(x))    :: locate_1d_dp

    INTEGER(i4) :: n, jl, jm, ju
    LOGICAL     :: ascnd
    INTEGER(i4) :: nx, i

    n = size(xx)
    ascnd = (xx(n) >= xx(1))

    nx = size(x)
    do i=1, nx
       jl = 0
       ju = n+1
       do
          if (ju-jl <= 1) exit
          jm = (ju+jl)/2
          if (ascnd .eqv. (x(i) >= xx(jm))) then
             jl = jm
          else
             ju = jm
          end if
       end do
       if (x(i) == xx(1)) then
          locate_1d_dp(i) = 1
       else if (x(i) == xx(n)) then
          locate_1d_dp(i) = n-1
       else
          locate_1d_dp(i) = jl
       end if
    end do

  END FUNCTION locate_1d_dp


  FUNCTION locate_0d_sp(xx,x)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN) :: xx
    REAL(sp),               INTENT(IN) :: x
    INTEGER(i4)                        :: locate_0d_sp

    INTEGER(i4) :: n, jl, jm, ju
    LOGICAL     :: ascnd

    n = size(xx)
    ascnd = (xx(n) >= xx(1))
    jl = 0
    ju = n+1
    do
       if (ju-jl <= 1) exit
       jm = (ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       end if
    end do
    if (x == xx(1)) then
       locate_0d_sp = 1
    else if (x == xx(n)) then
       locate_0d_sp = n-1
    else
       locate_0d_sp = jl
    end if

  END FUNCTION locate_0d_sp


  FUNCTION locate_1d_sp(xx,x)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN) :: xx
    REAL(sp), DIMENSION(:), INTENT(IN) :: x
    INTEGER(i4), DIMENSION(size(x))    :: locate_1d_sp

    INTEGER(i4) :: n, jl, jm, ju
    LOGICAL     :: ascnd
    INTEGER(i4) :: nx, i

    n = size(xx)
    ascnd = (xx(n) >= xx(1))

    nx = size(x)
    do i=1, nx
       jl = 0
       ju = n+1
       do
          if (ju-jl <= 1) exit
          jm = (ju+jl)/2
          if (ascnd .eqv. (x(i) >= xx(jm))) then
             jl = jm
          else
             ju = jm
          end if
       end do
       if (x(i) == xx(1)) then
          locate_1d_sp(i) = 1
       else if (x(i) == xx(n)) then
          locate_1d_sp(i) = n-1
       else
          locate_1d_sp(i) = jl
       end if
    end do

  END FUNCTION locate_1d_sp

  ! ------------------------------------------------------------------

END MODULE mo_interpol

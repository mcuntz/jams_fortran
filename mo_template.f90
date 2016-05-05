!> \file mo_template.f90

!> \brief Template for future module development.

!> \details This module serves as a template for future model developments.
!> It shows the module structure, coding style, and documentation.

!> \authors Matthias Cuntz, Juliane Mai, Christoph Schneider
!> \date Nov 2011
MODULE mo_template

  ! This module is a template for the JAMS Fortran library.

  ! The module provides an example procedures and demonstrates the agreed coding standard:
  !     - Use mo_kind, only: i4, sp, dp, lgt, etc.
  !     - Always make procedures for single and double precision (sp, dp), i.e. include interface to module procedures.
  !     - Include optional mask argument, if possible.
  !     - Make the module private by default.
  !     - Make routines available explicitly, i.e. public.
  !     - Break lines at column 130 at most, including comments.
  !     - Do not use tabs in files.
  !     - Give 1-line descriptions after the public definition.
  !     - Documentation:
  !       * document before the individual routines;
  !       * but document before the module interface not the separate routines,
  !         i.e. only one documentation per interface, no separate docu for sp and dp;
  !       * follow the documentation structure before the interface of mean in mo_template.f90;
  !       * break comment lines at column 130 at most as well.
  !     - Sort routines alphabetically in the file and in the public definitions.
  !     - The modules should be tested with at least two different compilers (of different vendors)
  !     - Use a subdirectory test_mo_xxx for testing where you simply link your modules:
  !       This means you do on the command prompt in the test directory:
  !         ln -s ../mo_kind.f90
  !         ln -s ../mo_xxx.f90
  !     - Always include the license terms of the JAMS Fortran library (below)
  !     - Add the "Note on Numerical Recipes License" if you use Numerical Recipes code in you module
  !     - Add a Copyright with the developping years and all names of developpers (see License below)

  ! Written  Matthias Cuntz, Nov 2011
  ! Modified Matthias Cuntz, Nov 2011 - add private
  !                          Nov 2011 - add public
  !          Matthias Cuntz & Christoph Schneider, Nov 2011 - doxygen

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

  ! Copyright 2011-2016 Matthias Cuntz, Juliane Mai


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

  ! Always use the number precisions of mo_kind
  USE mo_kind, ONLY: i4, sp, dp

  ! Of course
  IMPLICIT NONE

  ! Explicitly make public only the routines, parameters, etc. that shall be provided
  ! Sort alphabetically and give 1-line descriptions.
  PUBLIC :: circum  ! Circumference of circle
  PUBLIC :: mean    ! Average
  PUBLIC :: PI_dp   ! Constant Pi in double precision
  PUBLIC :: PI_sp   ! Constant Pi in single precision

  ! Documentation of routines is in front of the routines.
  ! But documentation of module interface routines must be in front of the interface definition.

  ! Interfaces for single and double precision routines; sort alphabetically
  ! Interfaces have to be in a public section for doxygen
  ! ------------------------------------------------------------------

  !     NAME
  !         mean

  !     PURPOSE
  !         Calculates the average value of a vector, i.e. the first moment of a series of numbers:
  !             mean = sum(x)/n
  !
  !>        \brief The average.
  !
  !>        \details Calculates the average value of a vector, i.e. the first moment of a series of numbers:
  !>        \f[ \bar{x} = \frac{1}{N} \sum_{i=1}^N x_i \f]
  !
  !>        If an optinal mask is given, the mean is only over those locations that correspond
  !>        to true values in the mask.\n
  !>        x can be single or double precision. The result will have the same numerical precision.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: dat(:)"        \f$ x_i \f$ 1D-array with input numbers
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: mask(:)" 1D-array with input mask\n
  !>                                                 If present, only those locations in dat corresponding
  !>                                                 to the true values in mask are used.
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>       \return     real(sp/dp) :: mean &mdash; \f$ \bar{x} \f$ average of all elements in vec
  !
  !     RESTRICTIONS
  !>       \note Input values must be floating points.
  !
  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = mean(vec, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         Sokal RR & Rohlf FJ - Biometry: the principle and practice of statistics in biological research,
  !             Freeman & Co., ISBN 0-7167-2411-1
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date Nov 2011
  !         Modified, Matthias Cuntz, Nov 2011 - include mask
  !                   Matthias Cuntz, Nov 2011 - test size(mask) == size(dat)
  INTERFACE mean
     MODULE PROCEDURE mean_sp, mean_dp
  END INTERFACE mean

  ! ------------------------------------------------------------------

  ! Make everything private by default.
  ! Do it after public interface definitions (for doxygen)
  PRIVATE

  ! ------------------------------------------------------------------

  ! Public parameters
  !> Constant Pi in double precision
  REAL(dp), PARAMETER :: PI_dp = 3.141592653589793238462643383279502884197_dp
  !> Constant Pi in single precision
  REAL(sp), PARAMETER :: PI_sp = 3.141592653589793238462643383279502884197_sp

  ! Private global parameters (not used, only for demonstration)
  INTEGER(i4), PARAMETER :: iTest=1

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         circum
  !
  !     PURPOSE
  !>        \brief Circumference of a circle
  !
  !>        \details Calculates the circumference of a circle
  !>        \f[ c = 2 \pi r \f]
  !
  !     CALLING SEQUENCE
  !         out = circum(radius)
  !
  !     INTENT(IN)
  !>        \param[in] "real(dp) :: radius"        Radius
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
  !>       \return     real(dp) :: circum &mdash; circumference of circle.
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         r = (/ 1., 2, 3., 5., 6. /)
  !         c = circum(r)
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date Dec 2012

  elemental pure function circum(radius)

    implicit none

    real(dp), intent(in) :: radius
    real(dp)             :: circum

    circum = 2.0_dp * pi_dp * radius

  end function circum

  ! ------------------------------------------------------------------

  FUNCTION mean_dp(dat, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)  :: dat
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                      :: mean_dp

    REAL(dp) :: n

    LOGICAL, DIMENSION(size(dat)) :: maske

    if (present(mask)) then
       if (size(mask) /= size(dat)) stop 'Error mean_dp: size(mask) /= size(dat)'
       maske = mask
       n = real(count(maske),dp)
    else
       maske(:) = .true.
       n = real(size(dat),dp)
    endif
    if (n <= (1.0_dp+tiny(1.0_dp))) stop 'mean_dp: n must be at least 2'

    ! Mean
    mean_dp  = sum(dat(:), mask=maske)/n

  END FUNCTION mean_dp


  FUNCTION mean_sp(dat, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)  :: dat
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                      :: mean_sp

    REAL(sp) :: n

    LOGICAL, DIMENSION(size(dat)) :: maske

    if (present(mask)) then
       if (size(mask) /= size(dat)) stop 'Error mean_sp: size(mask) /= size(dat)'
       maske = mask
       n = real(count(maske),sp)
    else
       maske(:) = .true.
       n = real(size(dat),sp)
    endif
    if (n <= (1.0_sp+tiny(1.0_sp))) stop 'mean_sp: n must be at least 2'

    ! Mean
    mean_sp  = sum(dat(:), mask=maske)/n

  END FUNCTION mean_sp

  ! ------------------------------------------------------------------

END MODULE mo_template

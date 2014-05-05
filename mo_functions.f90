!> \file mo_functions.f90

!> \brief Special functions such as gamma function.

!> \details This module provides special functions such as the gamma function.

!> \authors Matthias Cuntz
!> \date May 2014
MODULE mo_functions

  ! Provide special functions.

  ! Written  Matthias Cuntz, May 2014

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

  ! Copyright 2014 Matthias Cuntz, Juliane Mai

  USE mo_kind, ONLY: i4, i8, sp, dp

  IMPLICIT NONE

  PUBLIC :: factln    ! Ln(n!)
  PUBLIC :: factorial ! n!
  PUBLIC :: gammln    ! Ln(gamma)
  PUBLIC :: gamma     ! gamma

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         factln
  !
  !     PURPOSE
  !         Logarithm of factorial.
  !
  !>        \brief \f$ \ln(n!) \f$.
  !
  !>        \details Elemental function of the logarithm of the factorial:
  !>        \f[ \ln(n!) \f]
  !>        with
  !>        \f[ n! = 1 \dot 2 \dot ... \dot n \f]
  !>
  !>         Uses n! = \f$ \gamma(n+1) \f$.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8) :: n"        input
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
  !     RETURNS
  !>       \return     real(dp) :: factln &mdash; \f$ \ln(n!) \f$
  !
  !     RESTRICTIONS
  !>       \note return is always double precision.
  !
  !     EXAMPLE
  !         fac = exp(factln(5))
  !         fac -> 120
  !         -> see also example in test directory
  !
  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date May 2014
  INTERFACE factln
     MODULE PROCEDURE factln_i4, factln_i8
  END INTERFACE factln

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         factorial
  !
  !     PURPOSE
  !         Factorial.
  !
  !>        \brief \f$ n! \f$.
  !
  !>        \details Elemental function of the factorial:
  !>        \f[ n! = 1 \dot 2 \dot ... \dot n \f]
  !>
  !>         Uses n! = \f$ \gamma(n+1) \f$.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8) :: n"        input
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
  !     RETURNS
  !>       \return     real(dp) :: factln &mdash; \f$ \ln(n!) \f$
  !
  !     RESTRICTIONS
  !>       \note return is always double precision.
  !
  !     EXAMPLE
  !         fac = factorial(5)
  !         fac -> 120
  !         -> see also example in test directory
  !
  !     HISTORY
  !>        \author Matthias Cuntz, Juliane May
  !>        \date Feb 2013
  ! ------------------------------------------------------------------  
  INTERFACE factorial
     MODULE PROCEDURE factorial_i4, factorial_i8
  END INTERFACE factorial

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         gamma
  !
  !     PURPOSE
  !         The gamma function.
  !
  !>        \brief \f$ \gamma(z) \f$.
  !
  !>        \details Elemental function of the gamma function:
  !>        \f[ \gamma = \Int_0^\Inf t^{z-1} e^{-t} dt \f]
  !>
  !>         Uses Lanczos-type approximation to ln(gamma) for z > 0.
  !>         The function calculates gamma from abs(z).
  !>         Accuracy: About 14 significant digits except for small regions in the vicinity of 1 and 2.
  !
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: z"        input
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
  !     RETURNS
  !>       \return     real(sp/dp) :: gamma &mdash; \f$ \gamma(z) \f$
  !
  !     RESTRICTIONS
  !>       \note z must be >0. gammaln is calculated from abs(z).
  !
  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -9., 5., 6. /)
  !         gl  = gamma(vec)
  !         -> see also example in test directory
  !
  !     LITERATURE
  !         Lanczos, C. ''A precision approximation of the gamma function'', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
  !
  !     HISTORY
  !>        \author Alan Miller, 1 Creswick Street, Brighton, Vic. 3187, Australia, e-mail: amiller @ bigpond.net.au
  !>        \date Oct 1996
  !         Modified, Matthias Cuntz, May 2014
  INTERFACE gamma
     MODULE PROCEDURE gamma_sp, gamma_dp
  END INTERFACE gamma

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         gammln
  !
  !     PURPOSE
  !         Logarithm of the gamma function.
  !
  !>        \brief \f$ \ln(\gamma(z)) \f$.
  !
  !>        \details Elemental function of the logarithm of the gamma function:
  !>        \f[ \gamma = \Int_0^\Inf t^{z-1} e^{-t} dt \f]
  !>
  !>         Uses Lanczos-type approximation to ln(gamma) for z > 0.
  !>         The function calculates ln(gamma) from abs(z).
  !>         Accuracy: About 14 significant digits except for small regions in the vicinity of 1 and 2.
  !
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: z"        input
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
  !     RETURNS
  !>       \return     real(sp/dp) :: gammln &mdash; \f$ \ln\gamma(z) \f$
  !
  !     RESTRICTIONS
  !>       \note z must be >0. gammaln is calculated from abs(z).
  !
  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -9., 5., 6. /)
  !         gl  = gammln(vec)
  !         -> see also example in test directory
  !
  !     LITERATURE
  !         Lanczos, C. ''A precision approximation of the gamma function'', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
  !
  !     HISTORY
  !>        \author Alan Miller, 1 Creswick Street, Brighton, Vic. 3187, Australia, e-mail: amiller @ bigpond.net.au
  !>        \date Oct 1996
  !         Modified, Matthias Cuntz, May 2014
  INTERFACE gammln
     MODULE PROCEDURE gammln_sp, gammln_dp
  END INTERFACE gammln

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION factln_i4(n)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: n
    REAL(dp)                :: factln_i4

    factln_i4 = gammln(real(n,dp)+1.0_dp)

  end function factln_i4

  ELEMENTAL PURE FUNCTION factln_i8(n)

    IMPLICIT NONE

    INTEGER(i8), INTENT(IN) :: n
    REAL(dp)                :: factln_i8

    factln_i8 = gammln(real(n,dp)+1.0_dp)

  end function factln_i8

  ! ------------------------------------------------------------------

  elemental pure function factorial_i4(n)
    ! Returns the factorial n!
    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: n
    INTEGER(i4)             :: factorial_i4

    factorial_i4 = nint(exp(factln(n)),i4)

  end function factorial_i4

  elemental pure function factorial_i8(n)
    ! Returns the factorial n!
    IMPLICIT NONE

    INTEGER(i8), INTENT(IN) :: n
    INTEGER(i8)             :: factorial_i8

    factorial_i8 = nint(exp(factln(n)),i8)

  end function factorial_i8

  ! ------------------------------------------------------------------

  elemental pure function gamma_dp(z)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: z
    REAL(dp)             :: gamma_dp

    gamma_dp = exp(gammln(z))
    
  end function gamma_dp

  elemental pure function gamma_sp(z)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: z
    REAL(sp)             :: gamma_sp

    gamma_sp = exp(gammln(z))
    
  end function gamma_sp

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION gammln_dp(z)
    !  Uses Lanczos-type approximation to ln(gamma) for z > 0.
    !  Reference:
    !       Lanczos, C. ''A precision approximation of the gamma
    !               function'', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
    !  Accuracy: About 14 significant digits except for small regions
    !            in the vicinity of 1 and 2.
    !  Programmer: Alan Miller
    !              1 Creswick Street, Brighton, Vic. 3187, Australia
    !  e-mail: amiller @ bigpond.net.au
    !  Latest revision - 14 October 1996

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: z
    REAL(dp)             :: gammln_dp

    ! Local variables
    REAL(dp), parameter :: a(9) = (/ 0.9999999999995183_dp, 676.5203681218835_dp, &
         -1259.139216722289_dp, 771.3234287757674_dp, &
         -176.6150291498386_dp, 12.50734324009056_dp, &
         -0.1385710331296526_dp, 0.9934937113930748E-05_dp, &
         0.1659470187408462E-06_dp /)
    REAL(dp), parameter :: zero = 0.0_dp,   &
         one = 1.0_dp, &
         lnsqrt2pi =  0.9189385332046727_dp, &
         half = 0.5_dp, &
         sixpt5 = 6.5_dp, &
         seven = 7.0_dp
    REAL(dp)  :: tmp
    INTEGER   :: j

    ! if (z <= 0.0_dp) stop 'Error gammln_dp: z <= 0'
    gammln_dp = zero
    tmp = abs(z) + seven
    DO j = 9, 2, -1
       gammln_dp = gammln_dp + a(j)/tmp
       tmp = tmp - one
    END DO
    gammln_dp = gammln_dp + a(1)
    gammln_dp = LOG(gammln_dp) + lnsqrt2pi - (abs(z) + sixpt5) + (abs(z) - half)*LOG(abs(z) + sixpt5)

  END FUNCTION gammln_dp

  ELEMENTAL PURE FUNCTION gammln_sp(z)
    !  Uses Lanczos-type approximation to ln(gamma) for z > 0.
    !  Reference:
    !       Lanczos, C. ''A precision approximation of the gamma
    !               function'', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
    !  Accuracy: About 14 significant digits except for small regions
    !            in the vicinity of 1 and 2.
    !  Programmer: Alan Miller
    !              1 Creswick Street, Brighton, Vic. 3187, Australia
    !  e-mail: amiller @ bigpond.net.au
    !  Latest revision - 14 October 1996

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: z
    REAL(sp)             :: gammln_sp

    ! Local variables
    REAL(sp), parameter :: a(9) = (/ 0.9999999999995183_sp, 676.5203681218835_sp, &
         -1259.139216722289_sp, 771.3234287757674_sp, &
         -176.6150291498386_sp, 12.50734324009056_sp, &
         -0.1385710331296526_sp, 0.9934937113930748E-05_sp, &
         0.1659470187408462E-06_sp /)
    REAL(sp), parameter :: zero = 0.0_sp,   &
         one = 1.0_sp, &
         lnsqrt2pi =  0.9189385332046727_sp, &
         half = 0.5_sp, &
         sixpt5 = 6.5_sp, &
         seven = 7.0_sp
    REAL(sp)  :: tmp
    INTEGER   :: j

    ! if (z <= 0.0_sp) stop 'Error gammln_sp: z <= 0'
    gammln_sp = zero
    tmp = abs(z) + seven
    DO j = 9, 2, -1
       gammln_sp = gammln_sp + a(j)/tmp
       tmp = tmp - one
    END DO
    gammln_sp = gammln_sp + a(1)
    gammln_sp = LOG(gammln_sp) + lnsqrt2pi - (abs(z) + sixpt5) + (abs(z) - half)*LOG(abs(z) + sixpt5)

  END FUNCTION gammln_sp

  ! ------------------------------------------------------------------

END MODULE mo_functions

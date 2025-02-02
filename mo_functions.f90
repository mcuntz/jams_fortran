!> \file mo_functions.f90

!> \brief Special functions such as gamma function.

!> \details This module provides special functions such as the gamma function.

!> \authors Matthias Cuntz
!> \date May 2014
MODULE mo_functions

  ! Provide special functions.

  ! Written  Matthias Cuntz, May 2014
  ! Modified Matthias Cuntz, Dec 2017 - morris
  !          Stephan Thober, Aug 2020 - added beta function

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2014-2017 Matthias Cuntz - mc (at) macu (dot) de
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

  USE mo_kind, ONLY: i4, i8, sp, dp

  IMPLICIT NONE

  PUBLIC :: factln    ! ln(n!)
  PUBLIC :: factorial ! n!
  PUBLIC :: gammln    ! ln(gamma)
  PUBLIC :: gamm      ! gamm
  PUBLIC :: morris    ! elementary effects test function after Morris (1991)
  PUBLIC :: beta      ! beta function

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
  !>         Uses n! = \f$ \Gamma(n+1) \f$.
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
  !>         Uses n! = \f$ \Gamma(n+1) \f$.
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
  !>        \author Matthias Cuntz
  !>        \date Feb 2013
  ! ------------------------------------------------------------------
  INTERFACE factorial
     MODULE PROCEDURE factorial_i4, factorial_i8
  END INTERFACE factorial

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         gamm
  !
  !     PURPOSE
  !         The gamma function.
  !
  !>        \brief \f$ \Gamma(z) \f$.
  !
  !>        \details Elemental function of the gamma function:
  !>        \f[ \Gamma = \Int_0^\Inf t^{z-1} e^{-t} dt \f]
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
  !>       \return     real(sp/dp) :: gamm &mdash; \f$ \Gamma(z) \f$
  !
  !     RESTRICTIONS
  !>       \note z must be >0. gammln is calculated from abs(z).
  !
  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -9., 5., 6. /)
  !         gl  = gamm(vec)
  !         -> see also example in test directory
  !
  !     LITERATURE
  !         Lanczos, C. ''A precision approximation of the gamma function'', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
  !
  !     HISTORY
  !>        \author Alan Miller, 1 Creswick Street, Brighton, Vic. 3187, Australia, e-mail: amiller @ bigpond.net.au
  !>        \date Oct 1996
  !         Modified, Matthias Cuntz, May 2014
  INTERFACE gamm
     MODULE PROCEDURE gamm_sp, gamm_dp
  END INTERFACE gamm

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         gammln
  !
  !     PURPOSE
  !         Logarithm of the gamma function.
  !
  !>        \brief \f$ \ln(\Gamma(z)) \f$.
  !
  !>        \details Elemental function of the logarithm of the gamma function:
  !>        \f[ \Gamma = \Int_0^\Inf t^{z-1} e^{-t} dt \f]
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
  !>       \return     real(sp/dp) :: gammln &mdash; \f$ \ln\Gamma(z) \f$
  !
  !     RESTRICTIONS
  !>       \note z must be >0. gammln is calculated from abs(z).
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
  !
  !     NAME
  !         morris
  !
  !     PURPOSE
  !>        \details 20-dimension test function for elementary effects after Morris (1991).
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp), dimension(20[,npoints]) :: x"        number of points
  !>        \param[in] "real(sp/dp)                          :: beta0"    parameters
  !>        \param[in] "real(sp/dp), dimension(20)           :: beta1"    parameters
  !>        \param[in] "real(sp/dp), dimension(20,20)        :: beta2"    parameters
  !>        \param[in] "real(sp/dp), dimension(20,20,20)     :: beta3"    parameters
  !>        \param[in] "real(sp/dp), dimension(20,20,20,20)  :: beta4"    parameters
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
  !>       \return     real(sp/dp)[, dimension(npoints)] :: morris &mdash; result of Morris function
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         -> see example in test directory
  !
  !     LITERATURE
  !         Morris (1991), Factorial sampling plans for preliminary computational experiments,
  !         Technometrics 33, 161-174.
  !
  !     HISTORY
  !>        \author Matthias Cuntz, Juliane Mai, Mar 2015 in Python
  !>        \date Mar 2015
  !         Modified, Matthias Cuntz, Dec 2017 - ported to Fortran
  INTERFACE morris
     MODULE PROCEDURE morris_0d_dp, morris_1d_dp, &
          morris_0d_sp, morris_1d_sp
  END INTERFACE morris

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         beta
  !
  !     PURPOSE
  !>        \details calculates beta function
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp)                          :: alpha"    parameters
  !>        \param[in] "real(sp/dp)                          :: beta"     parameters
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
  !>       \return     real(sp/dp)[, dimension(npoints)] :: beta &mdash; result of beta function
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         -> see example in test directory
  !
  !     LITERATURE
  !         https://en.wikipedia.org/wiki/Beta_function
  !
  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Aug 2020
  !         Modified
  INTERFACE beta
     MODULE PROCEDURE beta_sp, beta_dp
  END INTERFACE beta

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

  elemental pure function gamm_dp(z)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: z
    REAL(dp)             :: gamm_dp

    gamm_dp = exp(gammln(z))

  end function gamm_dp

  elemental pure function gamm_sp(z)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: z
    REAL(sp)             :: gamm_sp

    gamm_sp = exp(gammln(z))

  end function gamm_sp

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

  FUNCTION morris_0d_dp(x, beta0, beta1, beta2, beta3, beta4)

    use mo_kind, only: dp, i4

    implicit none

    real(dp), dimension(20),          intent(in) :: x
    real(dp),                         intent(in) :: beta0
    real(dp), dimension(20),          intent(in) :: beta1
    real(dp), dimension(20,20),       intent(in) :: beta2
    real(dp), dimension(20,20,20),    intent(in) :: beta3
    real(dp), dimension(20,20,20,20), intent(in) :: beta4
    real(dp)                                     :: morris_0d_dp

    ! local
    real(dp),    dimension(20) :: om
    integer(i4), dimension(3) :: ii
    integer(i4) :: i, j, l, s, nn

    om = 2._dp*(x - 0.5_dp)
    ii = (/2, 4, 6/) + 1 ! +1 Python -> Fortran
    om(ii) = 2._dp*(1.1_dp*x(ii)/(x(ii)+0.1_dp) - 0.5_dp)

    nn = size(x)

    morris_0d_dp = beta0
    do i=1, nn
       morris_0d_dp = morris_0d_dp + beta1(i)*om(i)
       do j=i+1, nn
          morris_0d_dp = morris_0d_dp + beta2(i,j)*om(i)*om(j)
          do l=j+1, nn
             morris_0d_dp = morris_0d_dp + beta3(i,j,l)*om(i)*om(j)*om(l)
             do s=l+1, nn
                morris_0d_dp = morris_0d_dp + beta4(i,j,l,s)*om(i)*om(j)*om(l)*om(s)
             end do
          end do
       end do
    end do

  END FUNCTION morris_0d_dp

  FUNCTION morris_1d_dp(x, beta0, beta1, beta2, beta3, beta4)

    use mo_kind, only: dp, i4

    implicit none

    real(dp), dimension(:,:),         intent(in) :: x
    real(dp),                         intent(in) :: beta0
    real(dp), dimension(20),          intent(in) :: beta1
    real(dp), dimension(20,20),       intent(in) :: beta2
    real(dp), dimension(20,20,20),    intent(in) :: beta3
    real(dp), dimension(20,20,20,20), intent(in) :: beta4
    real(dp), dimension(size(x,2))               :: morris_1d_dp

    ! local
    real(dp),    dimension(20,size(x,2)) :: om
    integer(i4), dimension(3)            :: ii
    integer(i4) :: i, j, l, s, nn

    om = 2._dp*(x - 0.5_dp)
    ii = (/2, 4, 6/)
    om(ii,:) = 2._dp*(1.1_dp*x(ii,:)/(x(ii,:)+0.1_dp) - 0.5_dp)

    nn = size(x)

    morris_1d_dp(:) = beta0
    do i=1, nn
       morris_1d_dp(:) = morris_1d_dp(:) + beta1(i)*om(i,:)
       do j=i+1, nn
          morris_1d_dp(:) = morris_1d_dp(:) + beta2(i,j)*om(i,:)*om(j,:)
          do l=j+1, nn
             morris_1d_dp(:) = morris_1d_dp(:) + beta3(i,j,l)*om(i,:)*om(j,:)*om(l,:)
             do s=l+1, nn
                morris_1d_dp(:) = morris_1d_dp(:) + beta4(i,j,l,s)*om(i,:)*om(j,:)*om(l,:)*om(s,:)
             end do
          end do
       end do
    end do

  END FUNCTION morris_1d_dp


  FUNCTION morris_0d_sp(x, beta0, beta1, beta2, beta3, beta4)

    use mo_kind, only: sp, i4

    implicit none

    real(sp), dimension(20),          intent(in) :: x
    real(sp),                         intent(in) :: beta0
    real(sp), dimension(20),          intent(in) :: beta1
    real(sp), dimension(20,20),       intent(in) :: beta2
    real(sp), dimension(20,20,20),    intent(in) :: beta3
    real(sp), dimension(20,20,20,20), intent(in) :: beta4
    real(sp)                                     :: morris_0d_sp

    ! local
    real(sp),    dimension(20) :: om
    integer(i4), dimension(3) :: ii
    integer(i4) :: i, j, l, s, nn

    om = 2._sp*(x - 0.5_sp)
    ii = (/2, 4, 6/)
    om(ii) = 2._sp*(1.1_sp*x(ii)/(x(ii)+0.1_sp) - 0.5_sp)

    nn = size(x)

    morris_0d_sp = beta0
    do i=1, nn
       morris_0d_sp = morris_0d_sp + beta1(i)*om(i)
       do j=i+1, nn
          morris_0d_sp = morris_0d_sp + beta2(i,j)*om(i)*om(j)
          do l=j+1, nn
             morris_0d_sp = morris_0d_sp + beta3(i,j,l)*om(i)*om(j)*om(l)
             do s=l+1, nn
                morris_0d_sp = morris_0d_sp + beta4(i,j,l,s)*om(i)*om(j)*om(l)*om(s)
             end do
          end do
       end do
    end do

  END FUNCTION morris_0d_sp

  FUNCTION morris_1d_sp(x, beta0, beta1, beta2, beta3, beta4)

    use mo_kind, only: sp, i4

    implicit none

    real(sp), dimension(:,:),         intent(in) :: x
    real(sp),                         intent(in) :: beta0
    real(sp), dimension(20),          intent(in) :: beta1
    real(sp), dimension(20,20),       intent(in) :: beta2
    real(sp), dimension(20,20,20),    intent(in) :: beta3
    real(sp), dimension(20,20,20,20), intent(in) :: beta4
    real(sp), dimension(size(x,2))               :: morris_1d_sp

    ! local
    real(sp),    dimension(20,size(x,2)) :: om
    integer(i4), dimension(3)            :: ii
    integer(i4) :: i, j, l, s, nn

    om = 2._sp*(x - 0.5_sp)
    ii = (/2, 4, 6/)
    om(ii,:) = 2._sp*(1.1_sp*x(ii,:)/(x(ii,:)+0.1_sp) - 0.5_sp)

    nn = size(x)

    morris_1d_sp(:) = beta0
    do i=1, nn
       morris_1d_sp(:) = morris_1d_sp(:) + beta1(i)*om(i,:)
       do j=i+1, nn
          morris_1d_sp(:) = morris_1d_sp(:) + beta2(i,j)*om(i,:)*om(j,:)
          do l=j+1, nn
             morris_1d_sp(:) = morris_1d_sp(:) + beta3(i,j,l)*om(i,:)*om(j,:)*om(l,:)
             do s=l+1, nn
                morris_1d_sp(:) = morris_1d_sp(:) + beta4(i,j,l,s)*om(i,:)*om(j,:)*om(l,:)*om(s,:)
             end do
          end do
       end do
    end do

  END FUNCTION morris_1d_sp

  ! ------------------------------------------------------------------

  elemental pure function beta_sp(alpha, beta)

    use mo_kind, only: sp

    real(sp), intent(in) :: alpha
    real(sp), intent(in) :: beta
    real(sp)             :: beta_sp

    beta_sp = (gamm(alpha) * gamm(beta)) / gamm(alpha + beta)

  end function beta_sp

  elemental pure function beta_dp(alpha, beta)

    use mo_kind, only: dp

    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: beta
    real(dp)             :: beta_dp

    beta_dp = (gamm(alpha) * gamm(beta)) / gamm(alpha + beta)

  end function beta_dp

END MODULE mo_functions

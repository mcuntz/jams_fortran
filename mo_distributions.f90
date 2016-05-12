!> \file mo_distributions.f90

!> \brief Statistical distributions

!> \details Continuous distribution density functions such as the normal distribution.

!> \author Matthias Cuntz
!> \date Apr 2016

module mo_distributions

  ! This module provides distribution density functions such as the normal distribution.

  ! Written  Matthias Cuntz, Apr 2016

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

  ! Copyright 2016 Matthias Cuntz

  implicit none

  private

  public :: ep          ! Exponential Power
  public :: ep01        ! Exponential Power with mean=0, stddev=1
  public :: laplace     ! Laplace = Double Exponential
  public :: laplace01   ! Laplace = Double Exponential with mean=0, stddev=1
  public :: normal      ! Normal = Gauss
  public :: normal01    ! Normal = Gauss with mean=0, stddev=1
  public :: sep         ! Skew Exponential Power
  public :: sep01       ! Skew Exponential Power with mean=0, stddev=1
  public :: studentt    ! Student t
  public :: studentt01  ! Student t with mean=0, stddev=1
  public :: sstudentt   ! Skewed Student t
  public :: sstudentt01 ! Skewed Student t with mean=0, stddev=1

  ! ------------------------------------------------------------------

  !     NAME
  !         ep

  !     PURPOSE
  !>        \brief Exponential power distribution

  !>        \details The exponential power distribution with given mean \f$ \mu \f$,
  !>        standard deviation \f$ \sigma \f$ and kurtosis \f$ \beta \f$ parameter:
  !>        \f[ EP(x;\mu,\sigma,\beta) = SEP((x-\mu)/\sigma;0,1,\beta) / \sigma \f]
  !>        \f[ EP(x;0,1,\beta) = \omega_\beta
  !>            \exp \left\{ -c_\beta |x|^\frac{2}{1+\beta} \right\} \f]
  !>        with
  !>        \f[ c_\beta = \left( \frac{\Gamma(3[1+\beta]/2)}{\Gamma([1+\beta]/2)} \right)^\frac{1}{1+\beta} \f]
  !>        \f[ \omega_\beta = \frac{\sqrt(\Gamma(3[1+\beta]/2))}{(1+\beta)(\Gamma([1+\beta]/2))^{3/2}} \f]

  !     CALLING SEQUENCE
  !         out = ep(x, mu=mu, sig=sig, kurt=kurt)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: mu        mean (default: 0)"
  !>        \param[in] "real(sp/dp) :: sig       standard deviation (default: 1)"
  !>        \param[in] "real(sp/dp) :: kurt      kurtosis parameter (default: 0)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return     real(sp/dp) :: ep &mdash; \f$ EP(x;\mu,\sigma,\beta) \f$

  !     RESTRICTIONS
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = ep(vec, mean, sig, kurt)
  !         -> see also example in test directory

  !     LITERATURE
  !         Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and
  !             predictive inference of hydrologic models with correlated, heteroscedastic,
  !             and non-Gaussian errors. Water Resources Research 46, W10531.

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface ep
     module procedure ep_sp, ep_dp
  end interface ep

  ! ------------------------------------------------------------------

  !     NAME
  !         ep01

  !     PURPOSE
  !>        \brief Exponential power distribution with mean=0 and stddev=1

  !>        \details The exponential power distribution with zero mean and unit
  !>        standard deviation given kurtosis parameter \f$ \beta \f$:
  !>        \f[ EP(x;0,1,\beta) = \omega_\beta
  !>            \exp \left\{ -c_\beta |x|^\frac{2}{1+\beta} \right\} \f]
  !>        with
  !>        \f[ c_\beta = \left( \frac{\Gamma(3[1+\beta]/2)}{\Gamma([1+\beta]/2)} \right)^\frac{1}{1+\beta} \f]
  !>        \f[ \omega_\beta = \frac{\sqrt(\Gamma(3[1+\beta]/2))}{(1+\beta)(\Gamma([1+\beta]/2))^{3/2}} \f]

  !     CALLING SEQUENCE
  !         out = ep01(x, kurt=kurt)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: kurt      kurtosis parameter (default: 0)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return     real(sp/dp) :: sep01 &mdash; \f$ SEP(x;0,1,\xi,\beta) \f$

  !     RESTRICTIONS
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = ep01((vec-mean)/sig, 0.1)/sig
  !         -> see also example in test directory

  !     LITERATURE
  !         Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and
  !             predictive inference of hydrologic models with correlated, heteroscedastic,
  !             and non-Gaussian errors. Water Resources Research 46, W10531.

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface ep01
     module procedure ep01_sp, ep01_dp
  end interface ep01

  ! ------------------------------------------------------------------

  !     NAME
  !         laplace

  !     PURPOSE
  !>        \brief The Laplace or double exponential distribution

  !>        \details Samples the probability density function of a the Laplace distribution
  !>        given mean and standard deviation:\n
  !>        \f[ L(x;\mu,\sigma) = L((x-\mu)/\sigma;0,1) / \sigma \f]
  !>        \f[ L(z;0,1) = \frac{1}{2} \exp (-|z|) \f]

  !     CALLING SEQUENCE
  !         out = laplace(x, mu=mu, sig=sig)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: mu      mean (default: 0)"
  !>        \param[in] "real(sp/dp) :: sig     standard deviation (default: 1)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return     real(sp/dp) :: laplace &mdash; \f$ L(x;\mu,\sigma) \f$

  !     RESTRICTIONS
  !         Standard deviation must be greater than zero: sig > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = laplace(vec, mu=2, sig=1.5)
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface laplace
     module procedure laplace_sp, laplace_dp
  end interface laplace

  ! ------------------------------------------------------------------

  !     NAME
  !         laplace01

  !     PURPOSE
  !>        \brief The Laplace or double exponential distribution with mean=0 and stddev=1

  !>        \details Samples the probability density function of the Laplace distribution
  !>        with zero mean and unit standard deviation:\n
  !>        \f[ L(z;0,1) = \frac{1}{2} \exp (-|z|) \f]

  !     CALLING SEQUENCE
  !         out = laplace01(x)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

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
  !>       \return     real(sp/dp) :: laplace &mdash; \f$ L(x;0,1) \f$

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = laplace01((vec-mu)/sig)/sig
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface laplace01
     module procedure laplace01_sp, laplace01_dp
  end interface laplace01

  ! ------------------------------------------------------------------

  !     NAME
  !         normal

  !     PURPOSE
  !>        \brief The normal (gaussian) distribution

  !>        \details Samples the probability density function of a the normal distribution
  !>        given mean and standard deviation:\n
  !>        \f[ N(x;\mu,\sigma) = N((x-\mu)/\sigma;0,1) / \sigma \f]
  !>        \f[ N(z;0,1) = \frac{1}{\sqrt{2\pi}} \exp (-z^2/2) \f]

  !     CALLING SEQUENCE
  !         out = normal(x, mu=mu, sig=sig)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: mu      mean (default: 0)"
  !>        \param[in] "real(sp/dp) :: sig     standard deviation (default: 1)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return     real(sp/dp) :: normal &mdash; \f$ N(x;\mu,\sigma) \f$

  !     RESTRICTIONS
  !         Standard deviation must be greater than zero: sig > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = normal(vec, mu=2, sig=1.5)
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface normal
     module procedure normal_sp, normal_dp
  end interface normal

  ! ------------------------------------------------------------------

  !     NAME
  !         normal01

  !     PURPOSE
  !>        \brief The normal (gaussian) distribution with mean=0 and stddev=1

  !>        \details Samples the probability density function of the normal distribution
  !>        with zero mean and unit standard deviation:\n
  !>        \f[ N(z;0,1) = \frac{1}{\sqrt{2\pi}} \exp (-z^2/2) \f]

  !     CALLING SEQUENCE
  !         out = normal01(x)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

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
  !>       \return     real(sp/dp) :: normal &mdash; \f$ N(x;0,1) \f$

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = normal01((vec-mu)/sig)/sig
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface normal01
     module procedure normal01_sp, normal01_dp
  end interface normal01

  ! ------------------------------------------------------------------

  !     NAME
  !         sep

  !     PURPOSE
  !>        \brief Skew exponential power distribution

  !>        \details The skew exponential power distribution with given mean \f$ \mu \f$,
  !>        standard deviation \f$ \sigma \f$, skewness \f$ \xi \f$ and kurtosis \f$ \beta \f$ parameters:
  !>        \f[ SEP(x;\mu,\sigma,\xi,\beta) = SEP((x-\mu)/\sigma;0,1,\xi,\beta) / \sigma \f]
  !>        \f[ SEP(x;0,1,\xi,\beta) = \frac{2\sigma_\xi}{\xi+\xi^{-1}} \omega_\beta
  !>            \exp \left\{ -c_\beta |x_\xi|^\frac{2}{1+\beta} \right\} \f]
  !>        with
  !>        \f[ x_\xi = \xi^{\mathit{sign}(\mu_\xi + \sigma_\xi x)} \f]
  !>        \f[ c_\beta = \left( \frac{\Gamma(3[1+\beta]/2)}{\Gamma([1+\beta]/2)} \right)^\frac{1}{1+\beta} \f]
  !>        \f[ \omega_\beta = \frac{\sqrt(\Gamma(3[1+\beta]/2))}{(1+\beta)(\Gamma([1+\beta]/2))^{3/2}} \f]
  !>        \f[ \sigma_\xi = \sqrt{(M_2-M_1^2)(\xi^2+\xi^{-2})+2M_1^2-M_2} \f]
  !>        \f[ \mu_\xi = M_1(\xi-\xi^{-1}) \f]
  !>        \f[ M_1 = \frac{\Gamma(1+\beta)}{\sqrt{\Gamma(3[1+\beta]/2)\Gamma([1+\beta]/2)}} \f]
  !>        \f[ M_2 = 1 \f]
  !>        \f[ M_1 = \frac{\Gamma(1+\beta)}{\sqrt{\Gamma(3[1+\beta]/2)\Gamma([1+\beta]/2)}} \f]

  !     CALLING SEQUENCE
  !         out = sep(x, mu=mu, sig=sig, skew=skew, kurt=kurt)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: mu        mean (default: 0)"
  !>        \param[in] "real(sp/dp) :: sig       standard deviation (default: 1)"
  !>        \param[in] "real(sp/dp) :: skew      skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: kurt      kurtosis parameter (default: 0)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return     real(sp/dp) :: sep &mdash; \f$ SEP(x;\mu,\sigma,\xi,\beta) \f$

  !     RESTRICTIONS
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = sep(vec, mean, sig, skew, kurt)
  !         -> see also example in test directory

  !     LITERATURE
  !         Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and
  !             predictive inference of hydrologic models with correlated, heteroscedastic,
  !             and non-Gaussian errors. Water Resources Research 46, W10531.

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface sep
     module procedure sep_sp, sep_dp
  end interface sep

  ! ------------------------------------------------------------------

  !     NAME
  !         sep01

  !     PURPOSE
  !>        \brief Skew exponential power distribution with mean=0 and stddev=1

  !>        \details The skew exponential power distribution with zero mean and unit
  !>        standard deviation given skewness and kurtosis parameters \f$ \xi \f$ and \f$ \beta \f$:
  !>        \f[ SEP(x;0,1,\xi,\beta) = \frac{2\sigma_\xi}{\xi+\xi^{-1}} \omega_\beta
  !>            \exp \left\{ -c_\beta |x_\xi|^\frac{2}{1+\beta} \right\} \f]
  !>        with
  !>        \f[ x_\xi = \xi^{\mathit{sign}(\mu_\xi + \sigma_\xi x)} \f]
  !>        \f[ c_\beta = \left( \frac{\Gamma(3[1+\beta]/2)}{\Gamma([1+\beta]/2)} \right)^\frac{1}{1+\beta} \f]
  !>        \f[ \omega_\beta = \frac{\sqrt(\Gamma(3[1+\beta]/2))}{(1+\beta)(\Gamma([1+\beta]/2))^{3/2}} \f]
  !>        \f[ \sigma_\xi = \sqrt{(M_2-M_1^2)(\xi^2+\xi^{-2})+2M_1^2-M_2} \f]
  !>        \f[ \mu_\xi = M_1(\xi-\xi^{-1}) \f]
  !>        \f[ M_1 = \frac{\Gamma(1+\beta)}{\sqrt{\Gamma(3[1+\beta]/2)\Gamma([1+\beta]/2)}} \f]
  !>        \f[ M_2 = 1 \f]
  !>        \f[ M_1 = \frac{\Gamma(1+\beta)}{\sqrt{\Gamma(3[1+\beta]/2)\Gamma([1+\beta]/2)}} \f]

  !     CALLING SEQUENCE
  !         out = sep01(x, skew=skew, kurt=kurt)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: skew      skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: kurt      kurtosis parameter (default: 0)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return     real(sp/dp) :: sep01 &mdash; \f$ SEP(x;0,1,\xi,\beta) \f$

  !     RESTRICTIONS
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = sep01((vec-mean)/sig, 1.1, 0.1)/sig
  !         -> see also example in test directory

  !     LITERATURE
  !         Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and
  !             predictive inference of hydrologic models with correlated, heteroscedastic,
  !             and non-Gaussian errors. Water Resources Research 46, W10531.

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface sep01
     module procedure sep01_sp, sep01_dp
  end interface sep01

  ! ------------------------------------------------------------------

  !     NAME
  !         sstudentt

  !     PURPOSE
  !>        \brief Skewed Student's t distribution

  !>        \details The Skewed Student t-distribution given mean, standard deviation and degrees of freedom:
  !>        \f[ sstudentt(x;\nu,\mu,\sigma,\xi) = sstudentt((x-\mu)/\sigma;\nu,0,1,\xi) / \sigma \f]
  !>        \f[ sstudentt(x;\nu,0,1,\xi) = \frac{2\sigma_\xi}{\xi+\xi^{-1}}
  !>        \frac{\Gamma \left( \frac{\nu+1}{2} \right) \sqrt{\frac{\nu}{\nu-2}}}{
  !>        \Gamma \left( \frac{\nu}{2}\right) \sqrt{\pi \nu} }
  !>        \left( 1 + \frac{x_\xi^2}{\nu-2} \right) ^{-\frac{\nu+1}{2}} \f]
  !>        with
  !>        \f[ x_\xi = \xi^{\mathit{sign}(\mu_\xi + \sigma_\xi x)} \f]
  !>        \f[ \mu_\xi = \frac{2(\xi^2-\xi^{-2})}{\xi+\xi^{-1}}
  !>        \frac{\Gamma \left( \frac{\nu+1}{2} \right) \sqrt{\frac{\nu}{\nu-2}}}{
  !>        \Gamma \left( \frac{\nu}{2}\right) \sqrt{\pi \nu} } \frac{\nu-2}{\nu-1} \f]
  !>        \f[ \sigma_\xi = \sqrt{-\mu_\xi^2 + \frac{\xi^3+\xi^{-3}}{\xi+\xi^{-1}}} \f]

  !     CALLING SEQUENCE
  !         out = sstudentt(x, nu, mu=mu, sig=sig, skew=skew)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: mu      mean (default: 0)"
  !>        \param[in] "real(sp/dp) :: sig     standard deviation (default: 1)"
  !>        \param[in] "real(sp/dp) :: skew      skewness parameter (default: 1)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return     real(sp/dp) :: sstudentt &mdash; \f$ sstudentt(x;\nu,\mu,\sigma,\xi) \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than two: \f$\nu\f$ > 2.
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = sstudentt(vec, 3., mean, sig)
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface sstudentt
     module procedure sstudentt_sp, sstudentt_dp
  end interface sstudentt

  ! ------------------------------------------------------------------

  !     NAME
  !         sstudentt01

  !     PURPOSE
  !>        \brief Student's t distribution with mean=0 and stddev=1

  !>        \details The Student t-distribution with zero mean and unit
  !>        standard deviation:
  !>        \f[ sstudentt(x;\nu,0,1,\xi) = \frac{\Gamma \left( \frac{\nu+1}{2} \right) \sqrt{\frac{\nu}{\nu-2}}}{
  !>        \Gamma \left( \frac{\nu}{2}\right) \sqrt{\pi \nu} }
  !>        \left( 1 + \frac{x^2}{\nu-2} \right) ^{-\frac{\nu+1}{2}} \f]

  !     CALLING SEQUENCE
  !         out = sstudentt01(x, nu, skew=skew)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: skew      skewness parameter (default: 1)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return     real(sp/dp) :: sstudentt01 &mdash; \f$ sstudentt(x;\nu,0,1,\xi) \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than two: \f$\nu\f$ > 2.
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = sstudentt01((vec-mean)/sig, 3.)/sig
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface sstudentt01
     module procedure sstudentt01_sp, sstudentt01_dp
  end interface sstudentt01

  ! ------------------------------------------------------------------

  !     NAME
  !         studentt

  !     PURPOSE
  !>        \brief Student's t distribution

  !>        \details The Student t-distribution given mean, standard deviation and degrees of freedom:
  !>        \f[ studentt(x;\nu,\mu,\sigma) = studentt((x-\mu)/\sigma;\nu,0,1) / \sigma \f]
  !>        \f[ studentt(x;\nu,0,1) = \frac{\Gamma \left( \frac{\nu+1}{2} \right)}{
  !>        \Gamma \left( \frac{\nu}{2}\right) \sqrt{\pi \nu} }
  !>        \left( 1 + \frac{x^2}{\nu} \right) ^{-\frac{\nu+1}{2}} \f]

  !     CALLING SEQUENCE
  !         out = studentt(x, nu, mu=mu, sig=sig, skew)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: mu      mean (default: 0)"
  !>        \param[in] "real(sp/dp) :: sig     standard deviation (default: 1)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return     real(sp/dp) :: studentt &mdash; \f$ studentt(x;\nu,\mu,\sigma) \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than two: \f$\nu\f$ > 2.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = studentt(vec, 3., mean, sig)
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface studentt
     module procedure studentt_sp, studentt_dp
  end interface studentt

  ! ------------------------------------------------------------------

  !     NAME
  !         studentt01

  !     PURPOSE
  !>        \brief Student's t distribution with mean=0 and stddev=1

  !>        \details The Student t-distribution with zero mean and unit
  !>        standard deviation:
  !>        \f[ studentt(x;\nu,0,1) = \frac{\Gamma \left( \frac{\nu+1}{2} \right)}{
  !>        \Gamma \left( \frac{\nu}{2}\right) \sqrt{\pi \nu} }
  !>        \left( 1 + \frac{x^2}{\nu} \right) ^{-\frac{\nu+1}{2}} \f]

  !     CALLING SEQUENCE
  !         out = studentt01(x, nu)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

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
  !>       \return     real(sp/dp) :: studentt01 &mdash; \f$ studentt(x;\nu,0,1) \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than two: \f$\nu\f$ > 2.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = studentt01((vec-mean)/sig, 3.)/sig
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface studentt01
     module procedure studentt01_sp, studentt01_dp
  end interface studentt01

CONTAINS

  ! ------------------------------------------------------------------

  elemental pure function ep_dp(x, mu, sig, kurt)

    use mo_kind,      only: dp

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: mu
    real(dp), optional, intent(in) :: sig
    real(dp), optional, intent(in) :: kurt
    real(dp)                       :: ep_dp

    real(dp) :: imu, isig

    imu = 0.0_dp
    if (present(mu)) imu = mu
    isig = 1.0_dp
    if (present(sig)) isig = sig

    ep_dp = ep01((x-imu)/isig, kurt)/isig

  end function ep_dp

  elemental pure function ep_sp(x, mu, sig, kurt)

    use mo_kind,      only: sp

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: mu
    real(sp), optional, intent(in) :: sig
    real(sp), optional, intent(in) :: kurt
    real(sp)                       :: ep_sp

    real(sp) :: imu, isig

    imu = 0.0_sp
    if (present(mu)) imu = mu
    isig = 1.0_sp
    if (present(sig)) isig = sig

    ep_sp = ep01((x-imu)/isig, kurt)/isig

  end function ep_sp

  ! ------------------------------------------------------------------

  elemental pure function ep01_dp(x, kurt)

    use mo_kind,      only: dp
    use mo_utils,     only: ne
    use mo_functions, only: gamma

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: kurt
    real(dp)                       :: ep01_dp

    real(dp) :: beta
    real(dp) :: b1, b3
    real(dp) :: g1, g3
    real(dp) :: c_beta, om_beta

    beta = 0.0_dp
    if (present(kurt)) beta = kurt

    ! helpers
    if (ne(beta, -1.0_dp)) then
       b1 = 0.5_dp*(1.0_dp + beta)
       b3 = 1.5_dp*(1.0_dp + beta)
       g1 = gamma(b1)
       g3 = gamma(b3)
       ! -> 0
       c_beta = (g3/g1)**(1.0_dp/(1.0_dp+beta))
       ! -> sqrt(1/12)
       om_beta = sqrt(g3)/((1.0_dp+beta)*sqrt(g1**3))
    else
       c_beta  = 0.0_dp
       om_beta = sqrt(1.0_dp/12.0_dp)
    endif
    ! pdf
    if (abs(beta+1.0_dp) < 0.003_dp) then ! 2/(1-0.997) ~ 666
       ep01_dp = om_beta
    else
       ep01_dp = om_beta * exp(-c_beta*abs(x)**(2.0_dp/(1.0_dp+beta)))
    endif

  end function ep01_dp

  elemental pure function ep01_sp(x, kurt)

    use mo_kind,      only: sp
    use mo_utils,     only: ne
    use mo_functions, only: gamma

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: kurt
    real(sp)                       :: ep01_sp

    real(sp) :: beta
    real(sp) :: b1, b3
    real(sp) :: g1, g3
    real(sp) :: c_beta, om_beta

    beta = 0.0_sp
    if (present(kurt)) beta = kurt

    ! helpers
    if (ne(beta, -1.0_sp)) then
       b1 = 0.5_sp*(1.0_sp + beta)
       b3 = 1.5_sp*(1.0_sp + beta)
       g1 = gamma(b1)
       g3 = gamma(b3)
       ! -> 0
       c_beta = (g3/g1)**(1.0_sp/(1.0_sp+beta))
       ! -> sqrt(1/12)
       om_beta = sqrt(g3)/((1.0_sp+beta)*sqrt(g1**3))
    else
       c_beta  = 0.0_sp
       om_beta = sqrt(1.0_sp/12.0_sp)
    endif
    ! pdf
    if (abs(beta+1.0_sp) < 0.003_sp) then ! 2/(1-0.997) ~ 666
       ep01_sp = om_beta
    else
       ep01_sp = om_beta * exp(-c_beta*abs(x)**(2.0_sp/(1.0_sp+beta)))
    endif

  end function ep01_sp

  ! ------------------------------------------------------------------

  elemental pure function laplace_dp(x, mu, sig)

    use mo_kind,      only: dp

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: mu
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: laplace_dp

    real(dp) :: imu, isig

    imu = 0.0_dp
    if (present(mu)) imu = mu
    isig = 1.0_dp
    if (present(sig)) isig = sig

    laplace_dp = laplace01((x-imu)/isig)/isig

  end function laplace_dp

  elemental pure function laplace_sp(x, mu, sig)

    use mo_kind,      only: sp

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: mu
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: laplace_sp

    real(sp) :: imu, isig

    imu = 0.0_sp
    if (present(mu)) imu = mu
    isig = 1.0_sp
    if (present(sig)) isig = sig

    laplace_sp = laplace01((x-imu)/isig)/isig

  end function laplace_sp

  ! ------------------------------------------------------------------

  elemental pure function laplace01_dp(x)

    use mo_kind,      only: dp

    real(dp), intent(in) :: x
    real(dp)             :: laplace01_dp

    laplace01_dp = 0.5_dp * exp(-abs(x))

  end function laplace01_dp

  elemental pure function laplace01_sp(x)

    use mo_kind,      only: sp

    real(sp), intent(in) :: x
    real(sp)             :: laplace01_sp

    laplace01_sp = 0.5_sp * exp(-abs(x))

  end function laplace01_sp

  ! ------------------------------------------------------------------

  elemental pure function normal_dp(x, mu, sig)

    use mo_kind,      only: dp

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: mu
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: normal_dp

    real(dp) :: imu, isig

    imu = 0.0_dp
    if (present(mu)) imu = mu
    isig = 1.0_dp
    if (present(sig)) isig = sig

    normal_dp = normal01((x-imu)/isig)/isig

  end function normal_dp

  elemental pure function normal_sp(x, mu, sig)

    use mo_kind,      only: sp

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: mu
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: normal_sp

    real(sp) :: imu, isig

    imu = 0.0_sp
    if (present(mu)) imu = mu
    isig = 1.0_sp
    if (present(sig)) isig = sig

    normal_sp = normal01((x-imu)/isig)/isig

  end function normal_sp

  ! ------------------------------------------------------------------

  elemental pure function normal01_dp(x)

    use mo_kind,      only: dp
    use mo_constants, only: pi_dp

    real(dp), intent(in) :: x
    real(dp)             :: normal01_dp

    normal01_dp = 1.0_dp/sqrt(2.0_dp*pi_dp) * exp(-0.5_dp*x*x)

  end function normal01_dp

  elemental pure function normal01_sp(x)

    use mo_kind,      only: sp
    use mo_constants, only: pi_sp

    real(sp), intent(in) :: x
    real(sp)             :: normal01_sp

    normal01_sp = 1.0_sp/sqrt(2.0_sp*pi_sp) * exp(-0.5_sp*x*x)

  end function normal01_sp

  ! ------------------------------------------------------------------

  elemental pure function sep_dp(x, mu, sig, skew, kurt)

    use mo_kind,      only: dp

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: mu
    real(dp), optional, intent(in) :: sig
    real(dp), optional, intent(in) :: skew
    real(dp), optional, intent(in) :: kurt
    real(dp)                       :: sep_dp

    real(dp) :: imu, isig

    imu = 0.0_dp
    if (present(mu)) imu = mu
    isig = 1.0_dp
    if (present(sig)) isig = sig

    sep_dp = sep01((x-imu)/isig, skew, kurt)/isig

  end function sep_dp

  elemental pure function sep_sp(x, mu, sig, skew, kurt)

    use mo_kind,      only: sp

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: mu
    real(sp), optional, intent(in) :: sig
    real(sp), optional, intent(in) :: skew
    real(sp), optional, intent(in) :: kurt
    real(sp)                       :: sep_sp

    real(sp) :: imu, isig

    imu = 0.0_sp
    if (present(mu)) imu = mu
    isig = 1.0_sp
    if (present(sig)) isig = sig

    sep_sp = sep01((x-imu)/isig, skew, kurt)/isig

  end function sep_sp

  ! ------------------------------------------------------------------

  elemental pure function sep01_dp(x, skew, kurt)

    use mo_kind,      only: dp
    use mo_utils,     only: ne
    use mo_functions, only: gamma

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: skew
    real(dp), optional, intent(in) :: kurt
    real(dp)                       :: sep01_dp

    real(dp) :: xi, beta
    real(dp) :: b1, b2, b3
    real(dp) :: g1, g2, g3
    real(dp) :: M1, M2, xi1
    real(dp) :: mu_xi, sig_xi, c_beta, om_beta, a_xi

    xi = 1.0_dp
    if (present(skew)) xi = skew
    beta = 0.0_dp
    if (present(kurt)) beta = kurt

    ! helpers
    if (ne(beta, -1.0_dp)) then
       b1 = 0.5_dp*(1.0_dp + beta)
       b2 =         1.0_dp + beta
       b3 = 1.5_dp*(1.0_dp + beta)
       g1 = gamma(b1)
       g2 = gamma(b2)
       g3 = gamma(b3)
       ! -> sqrt(3/4)
       M1 = g2 / sqrt(g3*g1)
       ! -> 0
       c_beta = (g3/g1)**(1.0_dp/(1.0_dp+beta))
       ! -> sqrt(1/12)
       om_beta = sqrt(g3)/((1.0_dp+beta)*sqrt(g1**3))
    else
       M1      = sqrt(0.75_dp)
       c_beta  = 0.0_dp
       om_beta = sqrt(1.0_dp/12.0_dp)
    endif
    M2  = 1.0_dp
    xi1 = 1.0_dp/xi
    ! coefficients
    mu_xi   = M1*(xi-xi1)
    sig_xi  = (M2-M1*M1)*(xi*xi+xi1*xi1) + 2.0_dp*M1*M1 - M2
    if (sig_xi > 0.0_dp) then
       sig_xi  = sqrt(sig_xi)
    else
       sig_xi  = 0.0_dp
    endif
    a_xi    = mu_xi+sig_xi*x
    if (a_xi < 0.0_dp) then
       a_xi = a_xi * xi
    else
       a_xi = a_xi * xi1
    endif
    ! pdf
    if (abs(beta+1.0_dp) < 0.003_dp) then ! 2/(1-0.997) ~ 666
       sep01_dp = 2.0_dp*sig_xi/(xi+xi1) * om_beta
    else
       sep01_dp = 2.0_dp*sig_xi/(xi+xi1) * om_beta * exp(-c_beta*abs(a_xi)**(2.0_dp/(1.0_dp+beta)))
    endif

  end function sep01_dp

  elemental pure function sep01_sp(x, skew, kurt)

    use mo_kind,      only: sp
    use mo_utils,     only: ne
    use mo_functions, only: gamma

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: skew
    real(sp), optional, intent(in) :: kurt
    real(sp)                       :: sep01_sp

    real(sp) :: xi, beta
    real(sp) :: b1, b2, b3
    real(sp) :: g1, g2, g3
    real(sp) :: M1, M2, xi1
    real(sp) :: mu_xi, sig_xi, c_beta, om_beta, a_xi

    xi = 1.0_sp
    if (present(skew)) xi = skew
    beta = 0.0_sp
    if (present(kurt)) beta = kurt

    ! helpers
    if (ne(beta, -1.0_sp)) then
       b1 = 0.5_sp*(1.0_sp + beta)
       b2 =         1.0_sp + beta
       b3 = 1.5_sp*(1.0_sp + beta)
       g1 = gamma(b1)
       g2 = gamma(b2)
       g3 = gamma(b3)
       ! -> sqrt(3/4)
       M1 = g2 / sqrt(g3*g1)
       ! -> 0
       c_beta = (g3/g1)**(1.0_sp/(1.0_sp+beta))
       ! -> sqrt(1/12)
       om_beta = sqrt(g3)/((1.0_sp+beta)*sqrt(g1**3))
    else
       M1      = sqrt(0.75_sp)
       c_beta  = 0.0_sp
       om_beta = sqrt(1.0_sp/12.0_sp)
    endif
    M2  = 1.0_sp
    xi1 = 1.0_sp/xi
    ! coefficients
    mu_xi   = M1*(xi-xi1)
    sig_xi  = (M2-M1*M1)*(xi*xi+xi1*xi1) + 2.0_sp*M1*M1 - M2
    if (sig_xi > 0.0_sp) then
       sig_xi  = sqrt(sig_xi)
    else
       sig_xi  = 0.0_sp
    endif
    a_xi    = mu_xi+sig_xi*x
    if (a_xi < 0.0_sp) then
       a_xi = a_xi * xi
    else
       a_xi = a_xi * xi1
    endif
    ! pdf
    if (abs(beta+1.0_sp) < 0.003_sp) then ! 2/(1-0.997) ~ 666
       sep01_sp = 2.0_sp*sig_xi/(xi+xi1) * om_beta
    else
       sep01_sp = 2.0_sp*sig_xi/(xi+xi1) * om_beta * exp(-c_beta*abs(a_xi)**(2.0_sp/(1.0_sp+beta)))
    endif

  end function sep01_sp

  ! ------------------------------------------------------------------

  elemental pure function sstudentt_dp(x, nu, mu, sig, skew)

    use mo_kind,      only: dp

    real(dp),           intent(in) :: x
    real(dp),           intent(in) :: nu
    real(dp), optional, intent(in) :: mu
    real(dp), optional, intent(in) :: sig
    real(dp), optional, intent(in) :: skew
    real(dp)                       :: sstudentt_dp

    real(dp) :: imu, isig

    imu = 0.0_dp
    if (present(mu)) imu = mu
    isig = 1.0_dp
    if (present(sig)) isig = sig

    sstudentt_dp = sstudentt01((x-imu)/isig, nu, skew)/isig

  end function sstudentt_dp

  elemental pure function sstudentt_sp(x, nu, mu, sig, skew)

    use mo_kind,      only: sp

    real(sp),           intent(in) :: x
    real(sp),           intent(in) :: nu
    real(sp), optional, intent(in) :: mu
    real(sp), optional, intent(in) :: sig
    real(sp), optional, intent(in) :: skew
    real(sp)                       :: sstudentt_sp

    real(sp) :: imu, isig

    imu = 0.0_sp
    if (present(mu)) imu = mu
    isig = 1.0_sp
    if (present(sig)) isig = sig

    sstudentt_sp = sstudentt01((x-imu)/isig, nu, skew)/isig

  end function sstudentt_sp

  ! ------------------------------------------------------------------

  elemental pure function sstudentt01_dp(x, nu, skew)

    use mo_kind,      only: dp
    use mo_constants, only: pi_dp
    use mo_functions, only: gamma

    real(dp),           intent(in) :: x
    real(dp),           intent(in) :: nu
    real(dp), optional, intent(in) :: skew
    real(dp)                       :: sstudentt01_dp

    real(dp) :: xi, xi1
    real(dp) :: c
    real(dp) :: mu_xi, sig_xi, a_xi

    xi = 0.0_dp
    if (present(skew)) xi = skew

    xi1    = 1.0_dp/xi
    ! c      = gamma(0.5_dp*(nu+1.0_dp)) * sqrt(nu/(nu-2.0_dp)) / (gamma(0.5_dp*nu) * sqrt(pi_dp*nu))
    c      = gamma(0.5_dp*(nu+1.0_dp)) / (gamma(0.5_dp*nu) * sqrt(pi_dp*nu))
    mu_xi  = 2.0_dp * (xi*xi-xi1*xi)/(xi+xi1) * c * (nu-2.0_dp)/(nu-1.0_dp)
    sig_xi = -mu_xi*mu_xi + (xi**3+xi**3)/(xi+xi1)
    if (sig_xi > 0.0_dp) then
       sig_xi = sqrt(sig_xi)
    else
       sig_xi = 0.0_dp
    endif
    a_xi = mu_xi+sig_xi*x
    if (a_xi < 0.0_dp) then
       a_xi = a_xi * xi
    else
       a_xi = a_xi * xi1
    endif
    ! sstudentt01_dp = 2.0_dp * sig_xi / (xi+xi1) * c * (1.0_dp + a_xi*a_xi/(nu-2.0_dp))**(-0.5_dp*(nu+1.0_dp))
    sstudentt01_dp = 2.0_dp * sig_xi / (xi+xi1) * c * (1.0_dp + a_xi*a_xi/nu)**(-0.5_dp*(nu+1.0_dp))

  end function sstudentt01_dp

  elemental pure function sstudentt01_sp(x, nu, skew)

    use mo_kind,      only: sp
    use mo_constants, only: pi_sp
    use mo_functions, only: gamma

    real(sp),           intent(in) :: x
    real(sp),           intent(in) :: nu
    real(sp), optional, intent(in) :: skew
    real(sp)                       :: sstudentt01_sp

    real(sp) :: xi, xi1
    real(sp) :: c
    real(sp) :: mu_xi, sig_xi, a_xi

    xi = 0.0_sp
    if (present(skew)) xi = skew

    xi1    = 1.0_sp/xi
    c      = gamma(0.5_sp*(nu+1.0_sp)) * sqrt(nu/(nu-2.0_sp)) / (gamma(0.5_sp*nu) * sqrt(pi_sp*nu))
    mu_xi  = 2.0_sp * (xi*xi-xi1*xi)/(xi+xi1) * c * (nu-2.0_sp)/(nu-1.0_sp)
    sig_xi = -mu_xi*mu_xi + (xi**3+xi**3)/(xi+xi1)
    if (sig_xi > 0.0_sp) then
       sig_xi = sqrt(sig_xi)
    else
       sig_xi = 0.0_sp
    endif
    a_xi = mu_xi+sig_xi*x
    if (a_xi < 0.0_sp) then
       a_xi = a_xi * xi
    else
       a_xi = a_xi * xi1
    endif
    sstudentt01_sp = 2.0_sp * sig_xi / (xi+xi1) * c * (1.0_sp + a_xi*a_xi/(nu-2.0_sp))**(-0.5_sp*(nu+1.0_sp))

  end function sstudentt01_sp

  ! ------------------------------------------------------------------

  elemental pure function studentt_dp(x, nu, mu, sig)

    use mo_kind,      only: dp

    real(dp),           intent(in) :: x
    real(dp),           intent(in) :: nu
    real(dp), optional, intent(in) :: mu
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: studentt_dp

    real(dp) :: imu, isig

    imu = 0.0_dp
    if (present(mu)) imu = mu
    isig = 1.0_dp
    if (present(sig)) isig = sig

    studentt_dp = studentt01((x-imu)/isig, nu)/isig

  end function studentt_dp

  elemental pure function studentt_sp(x, nu, mu, sig)

    use mo_kind,      only: sp

    real(sp),           intent(in) :: x
    real(sp),           intent(in) :: nu
    real(sp), optional, intent(in) :: mu
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: studentt_sp

    real(sp) :: imu, isig

    imu = 0.0_sp
    if (present(mu)) imu = mu
    isig = 1.0_sp
    if (present(sig)) isig = sig

    studentt_sp = studentt01((x-imu)/isig, nu)/isig

  end function studentt_sp

  ! ------------------------------------------------------------------

  elemental pure function studentt01_dp(x, nu)

    use mo_kind,      only: dp
    use mo_constants, only: pi_dp
    use mo_functions, only: gamma

    real(dp), intent(in) :: x
    real(dp), intent(in) :: nu
    real(dp)             :: studentt01_dp

    real(dp) :: c

    c = gamma(0.5_dp*(nu+1.0_dp)) / (gamma(0.5_dp*nu) * sqrt(pi_dp*nu))
    studentt01_dp = c * (1.0_dp + x*x/nu)**(-0.5_dp*(nu+1.0_dp))

  end function studentt01_dp

  elemental pure function studentt01_sp(x, nu)

    use mo_kind,      only: sp
    use mo_constants, only: pi_sp
    use mo_functions, only: gamma

    real(sp), intent(in) :: x
    real(sp), intent(in) :: nu
    real(sp)             :: studentt01_sp

    real(sp) :: c

    c = gamma(0.5_sp*(nu+1.0_sp)) / (gamma(0.5_sp*nu) * sqrt(pi_sp*nu))
    studentt01_sp = c * (1.0_sp + x*x/nu)**(-0.5_sp*(nu+1.0_sp))

  end function studentt01_sp

  ! ------------------------------------------------------------------

end module mo_distributions

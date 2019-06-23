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
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2016 Matthias Cuntz - mc (at) macu (dot) de
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

  implicit none

  private

  public :: ep            ! Exponential Power of Box and Tiao (1992)
  public :: ep01          ! Exponential Power of Box and Tiao (1992) with location=0, scale=1
  public :: laplace       ! Laplace = Double Exponential
  public :: laplace01     ! Laplace = Double Exponential with location=0, scale=1
  public :: normal        ! Normal = Gauss
  public :: normal01      ! Normal = Gauss with location=0, scale=1
  public :: sep           ! Skew Exponential Power
  public :: sep01         ! Skew Exponential Power with location=0, scale=1
  public :: sep_fs        ! Skew Exponential Power after Fernandez and Steel (1998)
  public :: sep01_fs      ! Skew Exponential Power after Fernandez and Steel (1998) with location=0, scale=1
  ! public :: sep_fs_mean   ! Mean of skew Exponential Power after Fernandez and Steel (1998)
  ! public :: sep01_fs_mean ! Mean of skew Exponential Power after Fernandez and Steel (1998) with location=0, scale=1
  ! public :: sep_fs_std    ! Standard dev of skew Exponential Power after Fernandez and Steel (1998)
  ! public :: sep01_fs_std  ! Standard dev of skew Exponential Power after Fernandez and Steel (1998) with location=0, scale=1
  public :: st            ! Skew Student t
  public :: st01          ! Skew Student t with location=0, scale=1
  public :: st_fs         ! Skew Student t after Fernandez and Steel (1998)
  public :: st01_fs       ! Skew Student t after Fernandez and Steel (1998) with location=0, scale=1
  public :: st_fs_mean    ! mean of skew Student t after Fernandez and Steel (1998)
  public :: st01_fs_mean  ! mean of skew Student t after Fernandez and Steel (1998) with location=0, scale=1
  public :: st_fs_std     ! Standard dev of skew Student t after Fernandez and Steel (1998)
  public :: st01_fs_std   ! Standard dev of skew Student t after Fernandez and Steel (1998) with location=0, scale=1
  public :: t             ! Student t
  public :: t01           ! Student t with location=0, scale=1

  ! ------------------------------------------------------------------

  !     NAME
  !         ep

  !     PURPOSE
  !>        \brief Exponential power distribution of Box and Tiao (1992)

  !>        \details The exponential power distribution of Box and Tiao (1992) with given location \f$ loc \f$,
  !>        scale  \f$ sca \f$ or standard deviation \f$ \sigma \f$ and kurtosis \f$ \beta \f$ parameter:
  !>        \f[ EP(x;loc,sca,\beta) = EP((x-loc)/sac;0,1,\beta) / sca \f]
  !>        \f[ EP(x;0,1,\beta) = \omega_\beta
  !>            \exp \left\{ -c_\beta |x|^\frac{2}{1+\beta} \right\} \f]
  !>        with
  !>        \f[ c_\beta = \left( \frac{\Gamma(3[1+\beta]/2)}{\Gamma([1+\beta]/2)} \right)^\frac{1}{1+\beta} \f]
  !>        \f[ \omega_\beta = \frac{\sqrt(\Gamma(3[1+\beta]/2))}{(1+\beta)(\Gamma([1+\beta]/2))^{3/2}} \f]

  !     CALLING SEQUENCE
  !         out = ep(x, loc=loc, sca=sca, beta=beta, sig=sig)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: loc       location (default: 0)"
  !>        \param[in] "real(sp/dp) :: sca       scale (default: 1)"
  !>        \param[in] "real(sp/dp) :: beta      kurtosis parameter (default: 0)"
  !>        \param[in] "real(sp/dp) :: sig       standard deviation (default: ???), overwrites scale"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: ep &mdash; \f$ EP(x;loc,sca,\beta) \f$

  !     RESTRICTIONS
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = ep(vec, loc, sca, beta)
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
  !>        \brief Exponential power distribution of Box and Tiao (1992) with location=0 and scale=1

  !>        \details The exponential power distribution of Box and Tiao (1992) at location zero and with scale one
  !>        given kurtosis parameter \f$ \beta \f$:
  !>        \f[ EP(x;0,1,\beta) = \omega_\beta
  !>            \exp \left\{ -c_\beta |x|^\frac{2}{1+\beta} \right\} \f]
  !>        with
  !>        \f[ c_\beta = \left( \frac{\Gamma(3[1+\beta]/2)}{\Gamma([1+\beta]/2)} \right)^\frac{1}{1+\beta} \f]
  !>        \f[ \omega_\beta = \frac{\sqrt(\Gamma(3[1+\beta]/2))}{(1+\beta)(\Gamma([1+\beta]/2))^{3/2}} \f]

  !     CALLING SEQUENCE
  !         out = ep01(x, beta=beta)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: beta      kurtosis parameter (default: 0)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: ep01 &mdash; \f$ EP(x;0,1,\xi,\beta) \f$

  !     RESTRICTIONS
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = ep01((vec-loc)/sca, 0.1)/sca
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
  !>        given location \f$ loc \f$, scale \f$ sca \f$ or standard deviation \f$ \sigma \f$:\n
  !>        \f[ L(x;loc,sca) = L((x-loc)/sca;0,1) / sca \f]
  !>        \f[ L(z;0,1) = \frac{1}{2} \exp (-|z|) \f]

  !     CALLING SEQUENCE
  !         out = laplace(x, loc=loc, sca=sca, sig=sig)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: loc     location (default: 0)"
  !>        \param[in] "real(sp/dp) :: sca     scale (default: 1)"
  !>        \param[in] "real(sp/dp) :: sig     standard deviation (default: ???), overwrites scale"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: laplace &mdash; \f$ L(x;loc,sca) \f$

  !     RESTRICTIONS
  !         Standard deviation must be greater than zero: sca or sig > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = laplace(vec, loc=2, sig=1.5)
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
  !>        \brief The Laplace or double exponential distribution with location=0 and scale=1

  !>        \details Samples the probability density function of the Laplace distribution
  !>       at location zero and with scale one:\n
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
  !>       \return real(sp/dp) :: laplace &mdash; \f$ L(x;0,1) \f$

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = laplace01((vec-loc)/sca)/sca
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
  !>        given location \f$ loc \f$, scale  \f$ sca \f$ or standard deviation \f$ \sigma \f$:\n
  !>        \f[ N(x;loc,sca) = N((x-loc)/sca;0,1) / sca \f]
  !>        \f[ N(z;0,1) = \frac{1}{\sqrt{2\pi}} \exp (-z^2/2) \f]

  !     CALLING SEQUENCE
  !         out = normal(x, loc=loc, sca=sca, sig=sig)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: loc     location (default: 0)"
  !>        \param[in] "real(sp/dp) :: sca     scale (default: 1)"
  !>        \param[in] "real(sp/dp) :: sig     standard deviation (default: ???), overwrites scale"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: normal &mdash; \f$ N(x;loc,scale) \f$

  !     RESTRICTIONS
  !         Standard deviation must be greater than zero: sca or sig > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = normal(vec, loc=2, sig=1.5)
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
  !>        \brief The normal (gaussian) distribution with location=0 and scale=1

  !>        \details Samples the probability density function of the normal distribution
  !>        at location zero and with scale one:\n
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
  !>       \return real(sp/dp) :: normal &mdash; \f$ N(x;0,1) \f$

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = normal01((vec-loc)/sca)/sca
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

  !>        \details The skew exponential power distribution with given location \f$ loc \f$,
  !>        scale \f$ sca \f$ or standard deviation \f$ \sigma \f$, skewness \f$ \xi \f$ and
  !>        kurtosis \f$ \beta \f$ parameters:
  !>        \f[ SEP(x;loc,sca,\xi,\beta) = SEP((x-loc)/sca;0,1,\xi,\beta) / sca \f]

  !     CALLING SEQUENCE
  !         out = sep(x, loc=loc, sca=sca, xi=xi, beta=beta, sig=sig)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: loc       location (default: 0)"
  !>        \param[in] "real(sp/dp) :: sca       scale (default: 1)"
  !>        \param[in] "real(sp/dp) :: xi      skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: beta      kurtosis parameter (default: 0)"
  !>        \param[in] "real(sp/dp) :: sig       standard deviation (default: ???), owerwrites scale"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: sep &mdash; \f$ SEP(x;loc,scale,\xi,\beta) \f$

  !     RESTRICTIONS
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = sep(vec, loc, sca, xi, beta)
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
  !>        \brief Skew exponential power distribution with location=0 and scale=1

  !>        \details The skew exponential power distribution at location zero and with scale one
  !>        given skewness and kurtosis parameters \f$ \xi \f$ and \f$ \beta \f$:
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
  !         out = sep01(x, xi=xi, beta=beta)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: xi      skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: beta      kurtosis parameter (default: 0)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: sep01 &mdash; \f$ SEP(x;0,1,\xi,\beta) \f$

  !     RESTRICTIONS
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = sep01((vec-loc)/sca, 1.1, 0.1)/sca
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
  !         sep_fs

  !     PURPOSE
  !>        \brief Skew exponential power distribution after Fernandez and Steel (1998)

  !>        \details The skew exponential power distribution after Fernandez and Steel (1998)
  !>        with given location \f$ loc \f$, scale \f$ sca \f$ or standard deviation \f$ \sigma \f$,
  !>        skewness \f$ \xi \f$, and kurtosis \f$ \beta \f$ parameters:
  !>        \f[ SEP_{FS}(x;loc,sca,\xi,\beta) = SEP_{FS}((x-loc)/sca;0,1,\xi,\beta) / sca \f]

  !     CALLING SEQUENCE
  !         out = sep_fs(x, loc=loc, sca=sca, xi=xi, beta=beta, sig=sig)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: loc       location (default: 0)"
  !>        \param[in] "real(sp/dp) :: sca       scale (default: 1)"
  !>        \param[in] "real(sp/dp) :: xi        skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: beta      kurtosis parameter (default: 0)"
  !>        \param[in] "real(sp/dp) :: sig       standard deviation (default: 1), owerwrites scale"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: sep_fs &mdash; \f$ SEP_{FS}(x;loc,scale,\xi,\beta) \f$

  !     RESTRICTIONS
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = sep_fs(vec, loc, sca, xi, beta)
  !         -> see also example in test directory

  !     LITERATURE
  !         Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and
  !             predictive inference of hydrologic models with correlated, heteroscedastic,
  !             and non-Gaussian errors. Water Resources Research 46, W10531.

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface sep_fs
     module procedure sep_fs_sp, sep_fs_dp
  end interface sep_fs

  ! ------------------------------------------------------------------

  !     NAME
  !         sep01_fs

  !     PURPOSE
  !>        \brief Skew exponential power distribution after Fernandez and Steel (1998)
  !>        with location=0 and scale=1

  !>        \details The skew exponential power distribution after Fernandez and Steel (1998)
  !>        at location zero and with scale one, given skewness and kurtosis parameters \f$ \xi \f$
  !>        and \f$ \beta \f$:
  !>        \f[ SEP_fs(x;0,1,\xi,\beta) = \frac{2}{\xi+\xi^{-1}}
  !>            EP(xi^{-sign(x)}x;0,1,\beta) \f]

  !     CALLING SEQUENCE
  !         out = sep01_fs(x, xi=xi, beta=beta)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: xi        skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: beta      kurtosis parameter (default: 0)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: sep01_fs &mdash; \f$ SEP_fs(x;0,1,\xi,\beta) \f$

  !     RESTRICTIONS
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = sep01_fs((vec-loc)/sca, 1.1, 0.1)/sca
  !         -> see also example in test directory

  !     LITERATURE
  !         Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and
  !             predictive inference of hydrologic models with correlated, heteroscedastic,
  !             and non-Gaussian errors. Water Resources Research 46, W10531.

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface sep01_fs
     module procedure sep01_fs_sp, sep01_fs_dp
  end interface sep01_fs

  ! ------------------------------------------------------------------

  !     NAME
  !         sep_fs_mean

  !     PURPOSE
  !>        \brief Mean of skew exponential power distribution after Fernandez and Steel (1998)

  !>        \details The mean of the skew exponential power distribution after Fernandez and Steel (1998)
  !>        with given location \f$ loc \f$, scale \f$ sca \f$ or standard deviation \f$ \sigma \f$,
  !>        skewness \f$ \xi \f$, and kurtosis \f$ \beta \f$ parameters:
  !>        \f[ sep_fs_mean = \int_{-\infty}^{\infty} x SEP_{FS}(x;loc,sca,\xi,\beta) dx \f]

  !     CALLING SEQUENCE
  !         out = sep_fs_mean(loc=loc, sca=sca, xi=xi, beta=beta, sig=sig)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: loc       location (default: 0)"
  !>        \param[in] "real(sp/dp) :: sca       scale (default: 1)"
  !>        \param[in] "real(sp/dp) :: xi        skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: beta      kurtosis parameter (default: 0)"
  !>        \param[in] "real(sp/dp) :: sig       standard deviation (default: 1), owerwrites scale"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: sep_fs_mean &mdash; \f$ \int_{-\infty}^{\infty} x SEP_{FS}(x;loc,scale,\xi,\beta) dx \f$

  !     RESTRICTIONS
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = sep_fs_mean(loc, sca, xi, beta)
  !         -> see also example in test directory

  !     LITERATURE
  !         Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and
  !             predictive inference of hydrologic models with correlated, heteroscedastic,
  !             and non-Gaussian errors. Water Resources Research 46, W10531.

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  ! interface sep_fs_mean
  !    module procedure sep_fs_mean_sp, sep_fs_mean_dp
  ! end interface sep_fs_mean

  ! ------------------------------------------------------------------

  !     NAME
  !         sep01_fs_mean

  !     PURPOSE
  !>        \brief Mean of skew exponential power distribution after Fernandez and Steel (1998)
  !>        with location=0 and scale=1

  !>        \details The mean of the skew exponential power distribution after Fernandez and Steel (1998)
  !>        at location zero and with scale one, given skewness and kurtosis parameters \f$ \xi \f$
  !>        and \f$ \beta \f$:
  !>        \f[ sep01_fs_mean = \int_{-\infty}^{\infty} x SEP_{FS}(x;0,1,\xi,\beta) dx \f]

  !     CALLING SEQUENCE
  !         out = sep01_fs_mean(xi=xi, beta=beta)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: xi        skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: beta      kurtosis parameter (default: 0)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: sep01_fs_mean &mdash; \f$ \int_{-\infty}^{\infty} x SEP_{FS}(x;0,1,\xi,\beta) dx \f$

  !     RESTRICTIONS
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = sep01_fs_mean(1.1, 0.1)/sca
  !         -> see also example in test directory

  !     LITERATURE
  !         Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and
  !             predictive inference of hydrologic models with correlated, heteroscedastic,
  !             and non-Gaussian errors. Water Resources Research 46, W10531.

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  ! interface sep01_fs_mean
  !    module procedure sep01_fs_mean_sp, sep01_fs_mean_dp
  ! end interface sep01_fs_mean

  ! ------------------------------------------------------------------

  !     NAME
  !         sep_fs_std

  !     PURPOSE
  !>        \brief Standard deviation of the skew exponential power distribution after Fernandez and Steel (1998)

  !>        \details The standard deviation of the skew exponential power distribution
  !>        after Fernandez and Steel (1998)
  !>        with given location \f$ loc \f$, scale \f$ sca \f$ or standard deviation \f$ \sigma \f$,
  !>        skewness \f$ \xi \f$, and kurtosis \f$ \beta \f$ parameters:
  !>        \f[ sep_fs_std^2 = \int_{-\infty}^{\infty} x^2 SEP_{FS}(x;loc,sca,\xi,\beta) dx \f]

  !     CALLING SEQUENCE
  !         out = sep_fs_std(sca=sca, xi=xi, beta=beta, sig=sig)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: sca       scale (default: 1)"
  !>        \param[in] "real(sp/dp) :: xi        skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: beta      kurtosis parameter (default: 0)"
  !>        \param[in] "real(sp/dp) :: sig       standard deviation (default: 1), owerwrites scale"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: sep_fs_std &mdash; \f$ \sqrt \int_{-\infty}^{\infty} x^2 SEP_{FS}(x;loc,sca,\xi,\beta) dx \f$

  !     RESTRICTIONS
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = sep_fs_std(sca, xi, beta)
  !         -> see also example in test directory

  !     LITERATURE
  !         Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and
  !             predictive inference of hydrologic models with correlated, heteroscedastic,
  !             and non-Gaussian errors. Water Resources Research 46, W10531.

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  ! interface sep_fs_std
  !    module procedure sep_fs_std_sp, sep_fs_std_dp
  ! end interface sep_fs_std

  ! ------------------------------------------------------------------

  !     NAME
  !         sep01_fs_std

  !     PURPOSE
  !>        \brief Standard deviation of the skew exponential power distribution after Fernandez and Steel (1998)
  !>        with location=0 and scale=1

  !>        \details The standard deviation of the skew exponential power distribution
  !>        after Fernandez and Steel (1998)
  !>        at location zero and with scale one, given skewness and kurtosis parameters \f$ \xi \f$
  !>        and \f$ \beta \f$:
  !>        \f[ sep01_fs_std^2 = \int_{-\infty}^{\infty} x^2 SEP_{FS}(x;0,1,\xi,\beta) dx \f]

  !     CALLING SEQUENCE
  !         out = sep01_fs_std(xi=xi, beta=beta)

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: xi        skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: beta      kurtosis parameter (default: 0)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: sep01_fs_std &mdash; \f$ \sqrt \int_{-\infty}^{\infty} x^2 SEP_{FS}(x;0,1,\xi,\beta) dx \f$

  !     RESTRICTIONS
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.
  !         Kurtosis factor must be greater than -1 and lower +1: -1 < \f$\beta\f$ < 1.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = sep01_fs_std(1.1, 0.1)
  !         -> see also example in test directory

  !     LITERATURE
  !         Schoups G & Vrugt JA (2010) A formal likelihood function for parameter and
  !             predictive inference of hydrologic models with correlated, heteroscedastic,
  !             and non-Gaussian errors. Water Resources Research 46, W10531.

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  ! interface sep01_fs_std
  !    module procedure sep01_fs_std_sp, sep01_fs_std_dp
  ! end interface sep01_fs_std

  ! ------------------------------------------------------------------

  !     NAME
  !         st

  !     PURPOSE
  !>        \brief Skewed Student's t distribution

  !>        \details The Skewed Student t-distribution given location \f$ loc \f$,
  !>        scale \f$ sca \f$ or standard deviation \f$ \sigma \f$ and degrees of freedom \f$ \nu \f$:
  !>        \f[ st(x;\nu,loc,sca,\xi) = st((x-loc)/sca;\nu,0,1,\xi) / sca \f]
  !>        \f[ st(x;\nu,0,1,\xi) = \frac{2\sigma_\xi}{\xi+\xi^{-1}}
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
  !         out = st(x, nu, loc=loc, sca=sca, xi=xi, sig=sig)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: loc      location (default: 0)"
  !>        \param[in] "real(sp/dp) :: sca      scale (default: 1)"
  !>        \param[in] "real(sp/dp) :: xi     skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: sig      standard deviation (default: ???), overwrites scale"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: st &mdash; \f$ st(x;\nu,loc,sca,\xi) \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than one: \f$\nu\f$ > 1.
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = st(vec, 3., loc, sca)
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface st
     module procedure st_sp, st_dp
  end interface st

  ! ------------------------------------------------------------------

  !     NAME
  !         st01

  !     PURPOSE
  !>        \brief Skewed Student's t distribution at location zero and scale 1.

  !>        \details The Skewed Student t-distribution at location zero and scale one,
  !>        given degrees of freedom \f$ \nu \f$:
  !>        \f[ st(x;\nu,0,1,\xi) = \frac{2\sigma_\xi}{\xi+\xi^{-1}}
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
  !         out = st01(x, nu, xi=xi)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: xi      skewness parameter (default: 1)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: st01 &mdash; \f$ st(x;\nu,0,1,\xi) \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than one, better greater than two: \f$\nu\f$ > 1 (2).
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = st01((vec-loc)/sca, 3.)/sca
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface st01
     module procedure st01_sp, st01_dp
  end interface st01

  ! ------------------------------------------------------------------

  !     NAME
  !         st_fs

  !     PURPOSE
  !>        \brief Skewed Student's t distribution after Fernandez and Steel (1998)

  !>        \details The Skewed Student t-distribution after Fernandez and Steel (1998)
  !>        given location \f$ loc \f$, scale \f$ sca \f$ or standard deviation \f$ \sigma \f$
  !>        and degrees of freedom \f$ \nu \f$:
  !>        \f[ t_{FS}(x;loc,sca,\xi,\beta) = \frac{2}{\xi+\xi^{-1}}
  !>            t(xi^{-sign(x)}x;loc,sca,\beta) \f]

  !     CALLING SEQUENCE
  !         out = st_fs(x, nu, loc=loc, sca=sca, xi=xi, sig=sig)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: loc      location (default: 0)"
  !>        \param[in] "real(sp/dp) :: sca      scale (default: 1)"
  !>        \param[in] "real(sp/dp) :: xi     skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: sig      standard deviation (default: ???), overwrites scale"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: st &mdash; \f$ st_fs(x;\nu,loc,sca,\xi) \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than one: \f$\nu\f$ > 1.
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = st_fs(vec, 3., loc, sca)
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface st_fs
     module procedure st_fs_sp, st_fs_dp
  end interface st_fs

  ! ------------------------------------------------------------------

  !     NAME
  !         st01_fs

  !     PURPOSE
  !>        \brief Skewed Student's t distribution after Fernandez and Steel (1998)
  !>        at location zero and scale one.

  !>        \details The Skewed Student t-distribution after Fernandez and Steel (1998)
  !>        at location zero and scale one, given degrees of freedom \f$ \nu \f$:
  !>        \f[ t_{FS}(x;0,1,\xi,\beta) = \frac{2}{\xi+\xi^{-1}}
  !>            t(xi^{-sign(x)}x;0,1,\beta) \f]

  !     CALLING SEQUENCE
  !         out = st01_fs(x, nu, xi=xi)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: xi      skewness parameter (default: 1)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: st01_fs &mdash; \f$ st_fs(x;\nu,0,1,\xi) \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than one, better greater than two: \f$\nu\f$ > 1 (2).
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = st01_fs((vec-loc)/sca, 3.)/sca
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface st01_fs
     module procedure st01_fs_sp, st01_fs_dp
  end interface st01_fs
  
  ! ------------------------------------------------------------------

  !     NAME
  !         st_fs_mean

  !     PURPOSE
  !>        \brief Mean of skewed Student's t distribution after Fernandez and Steel (1998)

  !>        \details The mean of the skew Student t-distribution given location \f$ loc \f$,
  !>        scale \f$ sca \f$ or standard deviation \f$ \sigma \f$ and degrees of freedom \f$ \nu \f$:
  !>        \f[ st_fs_mean = \int_{-\infty}^{\infty} x st_{FS}(x;loc,sca,\xi,\beta) dx \f]

  !     CALLING SEQUENCE
  !         out = st_fs_mean(nu, loc=loc, sca=sca, xi=xi, sig=sig)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: loc      location (default: 0)"
  !>        \param[in] "real(sp/dp) :: sca      scale (default: 1)"
  !>        \param[in] "real(sp/dp) :: xi       skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: sig      standard deviation (default: ???), overwrites scale"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: st_fs_mean &mdash; \f$ \int_{-\infty}^{\infty} x st_{FS}(x;loc,sca,\xi,\beta) dx \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than one: \f$\nu\f$ > 1.
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = st_fs_mean(3., loc, sca)
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface st_fs_mean
     module procedure st_fs_mean_sp, st_fs_mean_dp
  end interface st_fs_mean

  ! ------------------------------------------------------------------

  !     NAME
  !         st01_fs_mean

  !     PURPOSE
  !>        \brief Mean of skewed Student's t distribution after Fernandez and Steel (1998)
  !>        at location zero and scale one.

  !>        \details The mean of the skew Student t-distribution at location zero and scale one
  !>        given degrees of freedom \f$ \nu \f$:
  !>        \f[ st01_fs_mean = \int_{-\infty}^{\infty} x st_{FS}(x;0,1,\xi,\beta) dx \f]

  !     CALLING SEQUENCE
  !         out = st01_fs_mean(nu, xi=xi)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: xi      skewness parameter (default: 1)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: st01_fs_mean &mdash; \f$ \int_{-\infty}^{\infty} x st_{FS}(x;0,1,\xi,\beta) dx \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than one, better greater than two: \f$\nu\f$ > 1 (2).
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = st01_fs_mean((vec-loc)/sca, 3.)/sca
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface st01_fs_mean
     module procedure st01_fs_mean_sp, st01_fs_mean_dp
  end interface st01_fs_mean
  
  ! ------------------------------------------------------------------

  !     NAME
  !         st_fs_std

  !     PURPOSE
  !>        \brief Standard deviation of skewed Student's t distribution after Fernandez and Steel (1998)

  !>        \details The standard deviation of the skew Student t-distribution given location \f$ loc \f$,
  !>        scale \f$ sca \f$ or standard deviation \f$ \sigma \f$ and degrees of freedom \f$ \nu \f$:
  !>        \f[ st_fs_std^2 = \int_{-\infty}^{\infty} x^2 st_{FS}(x;loc,sca,\xi,\beta) dx \f]

  !     CALLING SEQUENCE
  !         out = st_fs_std(sca=sca, xi=xi, sig=sig)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: sca      scale (default: 1)"
  !>        \param[in] "real(sp/dp) :: xi       skewness parameter (default: 1)"
  !>        \param[in] "real(sp/dp) :: sig      standard deviation (default: ???), overwrites scale"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: st_fs_std &mdash; \f$ \sqrt \int_{-\infty}^{\infty} x^2 st_{FS}(x;loc,sca,\xi,\beta) dx \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than one: \f$\nu\f$ > 1.
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = st_fs_std(3., loc, sca)
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface st_fs_std
     module procedure st_fs_std_sp, st_fs_std_dp
  end interface st_fs_std

  ! ------------------------------------------------------------------

  !     NAME
  !         st01_fs_std

  !     PURPOSE
  !>        \brief Standard deviation of skewed Student's t distribution after Fernandez and Steel (1998)
  !>        at location zero and scale one.

  !>        \details The standard deviation of the skew Student t-distribution at location zero and scale one
  !>        given degrees of freedom \f$ \nu \f$:
  !>        \f[ st01_fs_std^2 = \int_{-\infty}^{\infty} x^2 st_{FS}(x;0,1,\xi,\beta) dx \f]

  !     CALLING SEQUENCE
  !         out = st01_fs_std(nu, xi=xi)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: xi      skewness parameter (default: 1)"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: st01_fs_std &mdash; \f$ \sqrt \int_{-\infty}^{\infty} x^2 st_{FS}(x;0,1,\xi,\beta) dx \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than one, better greater than two: \f$\nu\f$ > 1 (2).
  !         Skewness factor must be greater than zero: \f$\xi\f$ > 0.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = st01_fs_std(3., xi)
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface st01_fs_std
     module procedure st01_fs_std_sp, st01_fs_std_dp
  end interface st01_fs_std

  ! ------------------------------------------------------------------

  !     NAME
  !         t

  !     PURPOSE
  !>        \brief Student's t distribution

  !>        \details The Student t-distribution given location \f$ loc \f$,
  !>        scale \f$ sca \f$ or standard deviation \f$ \sigma \f$, and degrees of freedom \f$ \nu \f$:
  !>        \f[ t(x;\nu,loc,sca) = t((x-loc)/sca;\nu,0,1) / sca \f]
  !>        \f[ t(x;\nu,0,1) = \frac{\Gamma \left( \frac{\nu+1}{2} \right)}{
  !>        \Gamma \left( \frac{\nu}{2}\right) \sqrt{\pi \nu} }
  !>        \left( 1 + \frac{x^2}{\nu} \right) ^{-\frac{\nu+1}{2}} \f]

  !     CALLING SEQUENCE
  !         out = t(x, nu, loc=loc, sca=sca, xi=xi, sig=sig)

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x       input (array) number"
  !>        \param[in] "real(sp/dp) :: nu      degrees of freedom"

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(sp/dp) :: loc     location (default: 0)"
  !>        \param[in] "real(sp/dp) :: sca     scale (default: 1)"
  !>        \param[in] "real(sp/dp) :: sig     standard deviation (default: ???), overwrites scale"

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return real(sp/dp) :: t &mdash; \f$ t(x;\nu,loc,sca) \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than two: \f$\nu\f$ > 2.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = t(vec, 3., loc, sca)
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface t
     module procedure t_sp, t_dp
  end interface t

  ! ------------------------------------------------------------------

  !     NAME
  !         t01

  !     PURPOSE
  !>        \brief Student's t distribution with location=0 and scale=1

  !>        \details The Student t-distribution at location zero and with scale one:
  !>        \f[ t(x;\nu,0,1) = \frac{\Gamma \left( \frac{\nu+1}{2} \right)}{
  !>        \Gamma \left( \frac{\nu}{2}\right) \sqrt{\pi \nu} }
  !>        \left( 1 + \frac{x^2}{\nu} \right) ^{-\frac{\nu+1}{2}} \f]

  !     CALLING SEQUENCE
  !         out = t01(x, nu)

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
  !>       \return real(sp/dp) :: t01 &mdash; \f$ t(x;\nu,0,1) \f$

  !     RESTRICTIONS
  !         Degrees of freedom must be greater than two: \f$\nu\f$ > 2.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = t01((vec-loc)/sca, 3.)/sca
  !         -> see also example in test directory

  !     HISTORY
  !>        \author Written, Matthias Cuntz
  !>        \date Apr 2016
  interface t01
     module procedure t01_sp, t01_dp
  end interface t01

CONTAINS

  ! ------------------------------------------------------------------

  elemental pure function ep_dp(x, loc, sca, beta, sig)

    use mo_kind, only: dp

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: loc
    real(dp), optional, intent(in) :: sca
    real(dp), optional, intent(in) :: beta
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: ep_dp

    real(dp) :: iloc, isca

    iloc = 0.0_dp
    if (present(loc)) iloc = loc
    isca = 1.0_dp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig

    ep_dp = ep01((x-iloc)/isca, beta)/isca

  end function ep_dp

  elemental pure function ep_sp(x, loc, sca, beta, sig)

    use mo_kind, only: sp

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: loc
    real(sp), optional, intent(in) :: sca
    real(sp), optional, intent(in) :: beta
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: ep_sp

    real(sp) :: iloc, isca

    iloc = 0.0_sp
    if (present(loc)) iloc = loc
    isca = 1.0_sp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig

    ep_sp = ep01((x-iloc)/isca, beta)/isca

  end function ep_sp

  ! ------------------------------------------------------------------

  elemental pure function ep01_dp(x, beta)

    use mo_kind,      only: dp
    use mo_utils,     only: ne
    use mo_functions, only: gamm

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: beta
    real(dp)                       :: ep01_dp

    real(dp) :: ibeta
    real(dp) :: b1, b3
    real(dp) :: g1, g3
    real(dp) :: c_beta, om_beta
    real(dp) :: height, x1

    ibeta = 0.0_dp
    if (present(beta)) ibeta = beta

    ! helpers
    if (ne(ibeta, -1.0_dp)) then
       b1 = 0.5_dp*(1.0_dp + ibeta)
       b3 = 1.5_dp*(1.0_dp + ibeta)
       g1 = gamm(b1)
       g3 = gamm(b3)
       ! -> 0
       c_beta = (g3/g1)**(1.0_dp/(1.0_dp+ibeta))
       ! -> sqrt(1/12)
       om_beta = sqrt(g3)/((1.0_dp+ibeta)*sqrt(g1**3))
    else
       c_beta  = 0.0_dp
       om_beta = sqrt(1.0_dp/12.0_dp)
    endif
    ! pdf
    if (abs(ibeta+1.0_dp) < 0.003_dp) then ! 2/(1-0.997) ~ 666
       ! Uniform between [-x1,x1]
       height = om_beta
       x1 = 0.5_dp/height ! int(pdf) = 1 = 2*x1*height
       if ((x > (-x1)) .and. (x < x1)) then
          ep01_dp = height
       else
          ep01_dp = 0.0_dp
       endif
    else
       ep01_dp = om_beta * exp(-c_beta*abs(x)**(2.0_dp/(1.0_dp+ibeta)))
    endif

  end function ep01_dp

  elemental pure function ep01_sp(x, beta)

    use mo_kind,      only: sp
    use mo_utils,     only: ne
    use mo_functions, only: gamm

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: beta
    real(sp)                       :: ep01_sp

    real(sp) :: ibeta
    real(sp) :: b1, b3
    real(sp) :: g1, g3
    real(sp) :: c_beta, om_beta
    real(sp) :: height, x1

    ibeta = 0.0_sp
    if (present(beta)) ibeta = beta

    ! helpers
    if (ne(ibeta, -1.0_sp)) then
       b1 = 0.5_sp*(1.0_sp + ibeta)
       b3 = 1.5_sp*(1.0_sp + ibeta)
       g1 = gamm(b1)
       g3 = gamm(b3)
       ! -> 0
       c_beta = (g3/g1)**(1.0_sp/(1.0_sp+ibeta))
       ! -> sqrt(1/12)
       om_beta = sqrt(g3)/((1.0_sp+ibeta)*sqrt(g1**3))
    else
       c_beta  = 0.0_sp
       om_beta = sqrt(1.0_sp/12.0_sp)
    endif
    ! pdf
    if (abs(ibeta+1.0_sp) < 0.003_sp) then ! 2/(1-0.997) ~ 666
       ! Uniform between [-x1,x1]
       height = om_beta
       x1 = 0.5_sp/height ! int(pdf) = 1 = 2*x1*height
       if ((x > (-x1)) .and. (x < x1)) then
          ep01_sp = height
       else
          ep01_sp = 0.0_sp
       endif
    else
       ep01_sp = om_beta * exp(-c_beta*abs(x)**(2.0_sp/(1.0_sp+ibeta)))
    endif

  end function ep01_sp

  ! ------------------------------------------------------------------

  elemental pure function laplace_dp(x, loc, sca, sig)

    use mo_kind, only: dp

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: loc
    real(dp), optional, intent(in) :: sca
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: laplace_dp

    real(dp) :: iloc, isca

    iloc = 0.0_dp
    if (present(loc)) iloc = loc
    isca = 1.0_dp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig/sqrt(2.0_dp)

    laplace_dp = laplace01((x-iloc)/isca)/isca

  end function laplace_dp

  elemental pure function laplace_sp(x, loc, sca, sig)

    use mo_kind, only: sp

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: loc
    real(sp), optional, intent(in) :: sca
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: laplace_sp

    real(sp) :: iloc, isca

    iloc = 0.0_sp
    if (present(loc)) iloc = loc
    isca = 1.0_sp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig/sqrt(2.0_sp)

    laplace_sp = laplace01((x-iloc)/isca)/isca

  end function laplace_sp

  ! ------------------------------------------------------------------

  elemental pure function laplace01_dp(x)

    use mo_kind, only: dp

    real(dp), intent(in) :: x
    real(dp)             :: laplace01_dp

    laplace01_dp = 0.5_dp * exp(-abs(x))

  end function laplace01_dp

  elemental pure function laplace01_sp(x)

    use mo_kind, only: sp

    real(sp), intent(in) :: x
    real(sp)             :: laplace01_sp

    laplace01_sp = 0.5_sp * exp(-abs(x))

  end function laplace01_sp

  ! ------------------------------------------------------------------

  elemental pure function normal_dp(x, loc, sca, sig)

    use mo_kind, only: dp

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: loc
    real(dp), optional, intent(in) :: sca
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: normal_dp

    real(dp) :: iloc, isca

    iloc = 0.0_dp
    if (present(loc)) iloc = loc
    isca = 1.0_dp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig

    normal_dp = normal01((x-iloc)/isca)/isca

  end function normal_dp

  elemental pure function normal_sp(x, loc, sca, sig)

    use mo_kind, only: sp

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: loc
    real(sp), optional, intent(in) :: sca
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: normal_sp

    real(sp) :: iloc, isca

    iloc = 0.0_sp
    if (present(loc)) iloc = loc
    isca = 1.0_sp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig

    normal_sp = normal01((x-iloc)/isca)/isca

  end function normal_sp

  ! ------------------------------------------------------------------

  elemental pure function normal01_dp(x)

    use mo_kind,      only: dp
    use mo_constants, only: pi_dp

    real(dp), intent(in) :: x
    real(dp)             :: normal01_dp

    normal01_dp = 1.0_dp/sqrt(2.0_dp*pi_dp) * exp(-0.5_dp*x**2)

  end function normal01_dp

  elemental pure function normal01_sp(x)

    use mo_kind,      only: sp
    use mo_constants, only: pi_sp

    real(sp), intent(in) :: x
    real(sp)             :: normal01_sp

    normal01_sp = 1.0_sp/sqrt(2.0_sp*pi_sp) * exp(-0.5_sp*x**2)

  end function normal01_sp

  ! ------------------------------------------------------------------

  elemental pure function sep_dp(x, loc, sca, xi, beta, sig)

    use mo_kind, only: dp

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: loc
    real(dp), optional, intent(in) :: sca
    real(dp), optional, intent(in) :: xi
    real(dp), optional, intent(in) :: beta
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: sep_dp

    real(dp) :: iloc, isca

    iloc = 0.0_dp
    if (present(loc)) iloc = loc
    isca = 1.0_dp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig

    sep_dp = sep01((x-iloc)/isca, xi, beta)/isca

  end function sep_dp

  elemental pure function sep_sp(x, loc, sca, xi, beta, sig)

    use mo_kind, only: sp

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: loc
    real(sp), optional, intent(in) :: sca
    real(sp), optional, intent(in) :: xi
    real(sp), optional, intent(in) :: beta
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: sep_sp

    real(sp) :: iloc, isca

    iloc = 0.0_sp
    if (present(loc)) iloc = loc
    isca = 1.0_sp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig

    sep_sp = sep01((x-iloc)/isca, xi, beta)/isca

  end function sep_sp

  ! ------------------------------------------------------------------

  elemental pure function sep01_dp(x, xi, beta)

    use mo_kind,      only: dp

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: xi
    real(dp), optional, intent(in) :: beta
    real(dp)                       :: sep01_dp

    real(dp) :: ixi, ibeta
    real(dp) :: mu, sig, z

    ixi = 1.0_dp
    if (present(xi)) ixi = xi
    ibeta = 0.0_dp
    if (present(beta)) ibeta = beta

    mu  = sep01_fs_mean_dp(ixi, ibeta)
    sig = sep01_fs_std_dp(ixi, ibeta)
    z   = mu + sig*x

    sep01_dp = sig * sep01_fs(z, ixi, ibeta)

  end function sep01_dp

  elemental pure function sep01_sp(x, xi, beta)

    use mo_kind,      only: sp

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: xi
    real(sp), optional, intent(in) :: beta
    real(sp)                       :: sep01_sp

    real(sp) :: ixi, ibeta
    real(sp) :: mu, sig, z

    ixi = 1.0_sp
    if (present(xi)) ixi = xi
    ibeta = 0.0_sp
    if (present(beta)) ibeta = beta

    mu  = sep01_fs_mean_sp(ixi, ibeta)
    sig = sep01_fs_std_sp(ixi, ibeta)
    z   = mu + sig*x

    sep01_sp = sig * sep01_fs(z, ixi, ibeta)

  end function sep01_sp

  ! ------------------------------------------------------------------

  elemental pure function sep_fs_dp(x, loc, sca, xi, beta, sig)

    use mo_kind, only: dp

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: loc
    real(dp), optional, intent(in) :: sca
    real(dp), optional, intent(in) :: xi
    real(dp), optional, intent(in) :: beta
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: sep_fs_dp

    real(dp) :: iloc, isca

    iloc = 0.0_dp
    if (present(loc)) iloc = loc
    isca = 1.0_dp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig

    sep_fs_dp = sep01_fs((x-iloc)/isca, xi, beta)/isca

  end function sep_fs_dp

  elemental pure function sep_fs_sp(x, loc, sca, xi, beta, sig)

    use mo_kind, only: sp

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: loc
    real(sp), optional, intent(in) :: sca
    real(sp), optional, intent(in) :: xi
    real(sp), optional, intent(in) :: beta
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: sep_fs_sp

    real(sp) :: iloc, isca

    iloc = 0.0_sp
    if (present(loc)) iloc = loc
    isca = 1.0_sp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig

    sep_fs_sp = sep01_fs((x-iloc)/isca, xi, beta)/isca

  end function sep_fs_sp

  ! ------------------------------------------------------------------

  elemental pure function sep01_fs_dp(x, xi, beta)

    use mo_kind, only: dp

    real(dp),           intent(in) :: x
    real(dp), optional, intent(in) :: xi
    real(dp), optional, intent(in) :: beta
    real(dp)                       :: sep01_fs_dp

    real(dp) :: ixi, ibeta
    real(dp) :: alpha

    ixi = 1.0_dp
    if (present(xi)) ixi = xi
    ibeta = 0.0_dp
    if (present(beta)) ibeta = beta

    if (x < 0.0_dp) then
       alpha = ixi
    else
       alpha = 1.0_dp/ixi
    endif

    sep01_fs_dp = 2.0_dp/(ixi+1./ixi) * ep01(alpha*x, ibeta)

  end function sep01_fs_dp

  elemental pure function sep01_fs_sp(x, xi, beta)

    use mo_kind, only: sp

    real(sp),           intent(in) :: x
    real(sp), optional, intent(in) :: xi
    real(sp), optional, intent(in) :: beta
    real(sp)                       :: sep01_fs_sp

    real(sp) :: ixi, ibeta
    real(sp) :: alpha

    ixi = 1.0_sp
    if (present(xi)) ixi = xi
    ibeta = 0.0_sp
    if (present(beta)) ibeta = beta

    if (x < 0.0_sp) then
       alpha = ixi
    else
       alpha = 1.0_sp/ixi
    endif

    sep01_fs_sp = 2.0_sp/(ixi+1./ixi) * ep01(alpha*x, ibeta)

  end function sep01_fs_sp

  ! ------------------------------------------------------------------

  elemental pure function sep_fs_mean_dp(loc, sca, xi, beta, sig)

    use mo_kind, only: dp

    real(dp), optional, intent(in) :: loc
    real(dp), optional, intent(in) :: sca
    real(dp), optional, intent(in) :: xi
    real(dp), optional, intent(in) :: beta
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: sep_fs_mean_dp

    real(dp) :: iloc, isca

    iloc = 0.0_dp
    if (present(loc)) iloc = loc
    isca = 1.0_dp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig

    sep_fs_mean_dp = sep01_fs_mean_dp(xi, beta) * isca + iloc

  end function sep_fs_mean_dp

  elemental pure function sep_fs_mean_sp(loc, sca, xi, beta, sig)

    use mo_kind, only: sp

    real(sp), optional, intent(in) :: loc
    real(sp), optional, intent(in) :: sca
    real(sp), optional, intent(in) :: xi
    real(sp), optional, intent(in) :: beta
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: sep_fs_mean_sp

    real(sp) :: iloc, isca

    iloc = 0.0_sp
    if (present(loc)) iloc = loc
    isca = 1.0_sp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig

    sep_fs_mean_sp = sep01_fs_mean_sp(xi, beta) * isca + iloc

  end function sep_fs_mean_sp

  ! ------------------------------------------------------------------

  elemental pure function sep01_fs_mean_dp(xi, beta)

    use mo_kind,      only: dp
    use mo_utils,     only: ne
    use mo_functions, only: gamm

    real(dp), optional, intent(in) :: xi
    real(dp), optional, intent(in) :: beta
    real(dp)                       :: sep01_fs_mean_dp

    real(dp) :: ixi, ibeta
    real(dp) :: b1, b2, b3
    real(dp) :: g1, g2, g3
    real(dp) :: M1

    ixi = 1.0_dp
    if (present(xi)) ixi = xi
    ibeta = 0.0_dp
    if (present(beta)) ibeta = beta

    ! helpers
    if (ne(ibeta, -1.0_dp)) then
       b1 = 0.5_dp*(1.0_dp + ibeta)
       b2 =         1.0_dp + ibeta
       b3 = 1.5_dp*(1.0_dp + ibeta)
       g1 = gamm(b1)
       g2 = gamm(b2)
       g3 = gamm(b3)
       ! -> sqrt(3/4)
       M1 = g2 / sqrt(g3*g1)
    else
       M1 = sqrt(0.75_dp)
    endif

    if (ne(ixi,1.0_dp)) then
       sep01_fs_mean_dp = M1*(ixi-1./ixi)
    else
       sep01_fs_mean_dp = 0.0_dp
    endif

  end function sep01_fs_mean_dp

  elemental pure function sep01_fs_mean_sp(xi, beta)

    use mo_kind,      only: sp
    use mo_utils,     only: ne
    use mo_functions, only: gamm

    real(sp), optional, intent(in) :: xi
    real(sp), optional, intent(in) :: beta
    real(sp)                       :: sep01_fs_mean_sp

    real(sp) :: ixi, ibeta
    real(sp) :: b1, b2, b3
    real(sp) :: g1, g2, g3
    real(sp) :: M1

    ixi = 1.0_sp
    if (present(xi)) ixi = xi
    ibeta = 0.0_sp
    if (present(beta)) ibeta = beta

    ! helpers
    if (ne(ibeta, -1.0_sp)) then
       b1 = 0.5_sp*(1.0_sp + ibeta)
       b2 =         1.0_sp + ibeta
       b3 = 1.5_sp*(1.0_sp + ibeta)
       g1 = gamm(b1)
       g2 = gamm(b2)
       g3 = gamm(b3)
       ! -> sqrt(3/4)
       M1 = g2 / sqrt(g3*g1)
    else
       M1 = sqrt(0.75_sp)
    endif

    if (ne(ixi,1.0_sp)) then
       sep01_fs_mean_sp = M1*(ixi-1./ixi)
    else
       sep01_fs_mean_sp = 0.0_sp
    endif

  end function sep01_fs_mean_sp

  ! ------------------------------------------------------------------

  elemental pure function sep_fs_std_dp(sca, xi, beta, sig)

    use mo_kind, only: dp

    real(dp), optional, intent(in) :: sca
    real(dp), optional, intent(in) :: xi
    real(dp), optional, intent(in) :: beta
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: sep_fs_std_dp

    real(dp) :: isca

    isca = 1.0_dp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig

    sep_fs_std_dp = sep01_fs_std_dp(xi, beta) * isca

  end function sep_fs_std_dp

  elemental pure function sep_fs_std_sp(sca, xi, beta, sig)

    use mo_kind, only: sp

    real(sp), optional, intent(in) :: sca
    real(sp), optional, intent(in) :: xi
    real(sp), optional, intent(in) :: beta
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: sep_fs_std_sp

    real(sp) :: isca

    isca = 1.0_sp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig

    sep_fs_std_sp = sep01_fs_std_sp(xi, beta) * isca

  end function sep_fs_std_sp

  ! ------------------------------------------------------------------

  elemental pure function sep01_fs_std_dp(xi, beta)

    use mo_kind,      only: dp
    use mo_utils,     only: ne

    real(dp), optional, intent(in) :: xi
    real(dp), optional, intent(in) :: beta
    real(dp)                       :: sep01_fs_std_dp

    real(dp) :: ixi, ibeta
    real(dp) :: mu, var
    real(dp) :: M1, M2

    ixi = 1.0_dp
    if (present(xi)) ixi = xi
    ibeta = 0.0_dp
    if (present(beta)) ibeta = beta

    mu = sep01_fs_mean_dp(ixi, ibeta)
    if (ne(ixi,1.0_dp)) then
       M1 = mu / (ixi-1.0_dp/ixi)
    else
       M1 = 0.0_dp
    endif
    M2  = 1.0_dp
    var = (M2-M1**2)*(xi**2+1.0_dp/xi**2) + 2.0_dp*M1**2 - M2
    if (var > 0.0_dp) then
       sep01_fs_std_dp = sqrt(var)
    else
       sep01_fs_std_dp = 0.0_dp
    endif

  end function sep01_fs_std_dp

  elemental pure function sep01_fs_std_sp(xi, beta)

    use mo_kind,      only: sp
    use mo_utils,     only: ne

    real(sp), optional, intent(in) :: xi
    real(sp), optional, intent(in) :: beta
    real(sp)                       :: sep01_fs_std_sp

    real(sp) :: ixi, ibeta
    real(sp) :: mu, var
    real(sp) :: M1, M2

    ixi = 1.0_sp
    if (present(xi)) ixi = xi
    ibeta = 0.0_sp
    if (present(beta)) ibeta = beta

    mu = sep01_fs_mean_sp(ixi, ibeta)
    if (ne(ixi,1.0_sp)) then
       M1 = mu / (ixi-1.0_sp/ixi)
    else
       M1 = 0.0_sp
    endif
    M2  = 1.0_sp
    var = (M2-M1**2)*(xi**2+1.0_sp/xi**2) + 2.0_sp*M1**2 - M2
    if (var > 0.0_sp) then
       sep01_fs_std_sp = sqrt(var)
    else
       sep01_fs_std_sp = 0.0_sp
    endif

  end function sep01_fs_std_sp

  ! ------------------------------------------------------------------

  elemental pure function st_dp(x, nu, loc, sca, xi, sig)

    use mo_kind, only: dp

    real(dp),           intent(in) :: x
    real(dp),           intent(in) :: nu
    real(dp), optional, intent(in) :: loc
    real(dp), optional, intent(in) :: sca
    real(dp), optional, intent(in) :: xi
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: st_dp

    real(dp) :: iloc, isca

    iloc = 0.0_dp
    if (present(loc)) iloc = loc
    isca = 1.0_dp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig * sqrt((nu-2.0_dp)/nu)

    st_dp = st01((x-iloc)/isca, nu, xi)/isca

  end function st_dp

  elemental pure function st_sp(x, nu, loc, sca, xi, sig)

    use mo_kind, only: sp

    real(sp),           intent(in) :: x
    real(sp),           intent(in) :: nu
    real(sp), optional, intent(in) :: loc
    real(sp), optional, intent(in) :: sca
    real(sp), optional, intent(in) :: xi
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: st_sp

    real(sp) :: iloc, isca

    iloc = 0.0_sp
    if (present(loc)) iloc = loc
    isca = 1.0_sp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig * sqrt((nu-2.0_sp)/nu)

    st_sp = st01((x-iloc)/isca, nu, xi)/isca

  end function st_sp

  ! ------------------------------------------------------------------

  elemental pure function st01_dp(x, nu, xi)

    use mo_kind, only: dp

    real(dp),           intent(in) :: x
    real(dp),           intent(in) :: nu
    real(dp), optional, intent(in) :: xi
    real(dp)                       :: st01_dp

    real(dp) :: ixi
    real(dp) :: mu, sca, z

    ixi = 0.0_dp
    if (present(xi)) ixi = xi

    mu  = st01_fs_mean_dp(nu, ixi)
    sca = st01_fs_std_dp(nu, ixi) * sqrt((nu-2.)/nu)
    z   = mu + sca*x

    st01_dp = sca * st01_fs(z, nu, ixi)

  end function st01_dp

  elemental pure function st01_sp(x, nu, xi)

    use mo_kind, only: sp

    real(sp),           intent(in) :: x
    real(sp),           intent(in) :: nu
    real(sp), optional, intent(in) :: xi
    real(sp)                       :: st01_sp

    real(sp) :: ixi
    real(sp) :: mu, sca, z

    ixi = 0.0_sp
    if (present(xi)) ixi = xi

    mu  = st01_fs_mean_sp(nu, ixi)
    sca = st01_fs_std_sp(nu, ixi) * sqrt((nu-2.)/nu)
    z   = mu + sca*x

    st01_sp = sca * st01_fs(z, nu, ixi)

  end function st01_sp

  ! ------------------------------------------------------------------

  elemental pure function st_fs_dp(x, nu, loc, sca, xi, sig)

    use mo_kind, only: dp

    real(dp),           intent(in) :: x
    real(dp),           intent(in) :: nu
    real(dp), optional, intent(in) :: loc
    real(dp), optional, intent(in) :: sca
    real(dp), optional, intent(in) :: xi
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: st_fs_dp

    real(dp) :: iloc, isca

    iloc = 0.0_dp
    if (present(loc)) iloc = loc
    isca = 1.0_dp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig * sqrt((nu-2.0_dp)/nu)

    st_fs_dp = st01((x-iloc)/isca, nu, xi)/isca

  end function st_fs_dp

  elemental pure function st_fs_sp(x, nu, loc, sca, xi, sig)

    use mo_kind, only: sp

    real(sp),           intent(in) :: x
    real(sp),           intent(in) :: nu
    real(sp), optional, intent(in) :: loc
    real(sp), optional, intent(in) :: sca
    real(sp), optional, intent(in) :: xi
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: st_fs_sp

    real(sp) :: iloc, isca

    iloc = 0.0_sp
    if (present(loc)) iloc = loc
    isca = 1.0_sp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig * sqrt((nu-2.0_sp)/nu)

    st_fs_sp = st01((x-iloc)/isca, nu, xi)/isca

  end function st_fs_sp

  ! ------------------------------------------------------------------

  elemental pure function st01_fs_dp(x, nu, xi)

    use mo_kind, only: dp

    real(dp),           intent(in) :: x
    real(dp),           intent(in) :: nu
    real(dp), optional, intent(in) :: xi
    real(dp)                       :: st01_fs_dp

    real(dp) :: ixi
    real(dp) :: alpha

    ixi = 0.0_dp
    if (present(xi)) ixi = xi

    if (x < 0.0_dp) then
       alpha = ixi
    else
       alpha = 1.0_dp/ixi
    endif

    st01_fs_dp = 2.0_dp/(ixi+1./ixi) * t01(alpha*x, nu)

  end function st01_fs_dp

  elemental pure function st01_fs_sp(x, nu, xi)

    use mo_kind, only: sp

    real(sp),           intent(in) :: x
    real(sp),           intent(in) :: nu
    real(sp), optional, intent(in) :: xi
    real(sp)                       :: st01_fs_sp

    real(sp) :: ixi
    real(sp) :: alpha

    ixi = 0.0_sp
    if (present(xi)) ixi = xi

    if (x < 0.0_sp) then
       alpha = ixi
    else
       alpha = 1.0_sp/ixi
    endif

    st01_fs_sp = 2.0_sp/(ixi+1./ixi) * t01(alpha*x, nu)

  end function st01_fs_sp

  ! ------------------------------------------------------------------

  elemental pure function st_fs_mean_dp(nu, loc, sca, xi, sig)

    use mo_kind, only: dp

    real(dp),           intent(in) :: nu
    real(dp), optional, intent(in) :: loc
    real(dp), optional, intent(in) :: sca
    real(dp), optional, intent(in) :: xi
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: st_fs_mean_dp

    real(dp) :: iloc, isca

    iloc = 0.0_dp
    if (present(loc)) iloc = loc
    isca = 1.0_dp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig * sqrt((nu-2.0_dp)/nu)

    st_fs_mean_dp = st01_fs_mean_dp(nu, xi) * isca + iloc

  end function st_fs_mean_dp

  elemental pure function st_fs_mean_sp(nu, loc, sca, xi, sig)

    use mo_kind, only: sp

    real(sp),           intent(in) :: nu
    real(sp), optional, intent(in) :: loc
    real(sp), optional, intent(in) :: sca
    real(sp), optional, intent(in) :: xi
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: st_fs_mean_sp

    real(sp) :: iloc, isca

    iloc = 0.0_sp
    if (present(loc)) iloc = loc
    isca = 1.0_sp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig * sqrt((nu-2.0_sp)/nu)

    st_fs_mean_sp = st01_fs_mean_sp(nu, xi) * isca + iloc

  end function st_fs_mean_sp

  ! ------------------------------------------------------------------

  elemental pure function st01_fs_mean_dp(nu, xi)

    use mo_kind,      only: dp
    use mo_constants, only: pi_dp
    use mo_functions, only: gamm
    use mo_utils,     only: special_value
    use mo_utils,     only: le

    real(dp),           intent(in) :: nu
    real(dp), optional, intent(in) :: xi
    real(dp)                       :: st01_fs_mean_dp

    real(dp) :: ixi

    ixi = 0.0_dp
    if (present(xi)) ixi = xi

    if (le(nu,1.0_dp)) then
       st01_fs_mean_dp = special_value(1.0_dp, 'ieee_positive_inf')
    else
       st01_fs_mean_dp = 2.0_dp * (ixi-1./ixi) * gamm(0.5_dp*(nu+1.0_dp)) / gamm(0.5_dp*nu)
       st01_fs_mean_dp = st01_fs_mean_dp * (nu-2.0_dp)/(nu-1.0_dp) * nu/(nu-2.0_dp) / sqrt(pi_dp*nu)
    endif

  end function st01_fs_mean_dp

  elemental pure function st01_fs_mean_sp(nu, xi)

    use mo_kind,      only: sp
    use mo_constants, only: pi_sp
    use mo_functions, only: gamm
    use mo_utils,     only: special_value
    use mo_utils,     only: le

    real(sp),           intent(in) :: nu
    real(sp), optional, intent(in) :: xi
    real(sp)                       :: st01_fs_mean_sp

    real(sp) :: ixi

    ixi = 0.0_sp
    if (present(xi)) ixi = xi

    if (le(nu,1.0_sp)) then
       st01_fs_mean_sp = special_value(1.0_sp, 'ieee_positive_inf')
    else
       st01_fs_mean_sp = 2.0_sp * (ixi-1./ixi) * gamm(0.5_sp*(nu+1.0_sp)) / gamm(0.5_sp*nu)
       st01_fs_mean_sp = st01_fs_mean_sp * (nu-2.0_sp)/(nu-1.0_sp) * nu/(nu-2.0_sp) / sqrt(pi_sp*nu)
    endif

  end function st01_fs_mean_sp

  ! ------------------------------------------------------------------

  elemental pure function st_fs_std_dp(nu, sca, xi, sig)

    use mo_kind, only: dp

    real(dp),           intent(in) :: nu
    real(dp), optional, intent(in) :: sca
    real(dp), optional, intent(in) :: xi
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: st_fs_std_dp

    real(dp) :: isca

    isca = 1.0_dp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig * sqrt((nu-2.0_dp)/nu)

    st_fs_std_dp = st01_fs_std_dp(nu, xi) * isca

  end function st_fs_std_dp

  elemental pure function st_fs_std_sp(nu, sca, xi, sig)

    use mo_kind, only: sp

    real(sp),           intent(in) :: nu
    real(sp), optional, intent(in) :: sca
    real(sp), optional, intent(in) :: xi
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: st_fs_std_sp

    real(sp) :: isca

    isca = 1.0_sp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig * sqrt((nu-2.0_sp)/nu)

    st_fs_std_sp = st01_fs_std_sp(nu, xi) * isca

  end function st_fs_std_sp

  ! ------------------------------------------------------------------

  elemental pure function st01_fs_std_dp(nu, xi)

    use mo_kind,  only: dp
    use mo_utils, only: le, ne
    use mo_utils, only: special_value

    real(dp),           intent(in) :: nu
    real(dp), optional, intent(in) :: xi
    real(dp)                       :: st01_fs_std_dp

    real(dp) :: ixi
    real(dp) :: mu, var
    real(dp) :: M1, M2

    ixi = 0.0_dp
    if (present(xi)) ixi = xi

    if (le(nu,2.0_dp)) then
       st01_fs_std_dp = special_value(1.0_dp, 'ieee_positive_inf')
    else
       mu = st01_fs_mean_dp(nu, ixi)
       if (ne(ixi,1.0_dp)) then
          M1 = mu / (ixi-1./ixi)
       else
          M1 = 0.0_dp
       endif
       M2  = nu/(nu-2.0_dp)
       var = (M2 - M1**2) * (ixi**2 + 1.0_dp/ixi**2) + 2.0_dp*M1**2 - M2
       if (var > 0.0_dp) then
          st01_fs_std_dp = sqrt(var)
       else
          st01_fs_std_dp = 0.0_dp
       endif
    endif

  end function st01_fs_std_dp

  elemental pure function st01_fs_std_sp(nu, xi)

    use mo_kind,  only: sp
    use mo_utils, only: le, ne
    use mo_utils, only: special_value

    real(sp),           intent(in) :: nu
    real(sp), optional, intent(in) :: xi
    real(sp)                       :: st01_fs_std_sp

    real(sp) :: ixi
    real(sp) :: mu, var
    real(sp) :: M1, M2

    ixi = 0.0_sp
    if (present(xi)) ixi = xi

    if (le(nu,2.0_sp)) then
       st01_fs_std_sp = special_value(1.0_sp, 'ieee_positive_inf')
    else
       mu = st01_fs_mean_sp(nu, ixi)
       if (ne(ixi,1.0_sp)) then
          M1 = mu / (ixi-1./ixi)
       else
          M1 = 0.0_sp
       endif
       M2  = nu/(nu-2.0_sp)
       var = (M2 - M1**2) * (ixi**2 + 1.0_sp/ixi**2) + 2.0_sp*M1**2 - M2
       if (var > 0.0_sp) then
          st01_fs_std_sp = sqrt(var)
       else
          st01_fs_std_sp = 0.0_sp
       endif
    endif

  end function st01_fs_std_sp

  ! ------------------------------------------------------------------

  elemental pure function t_dp(x, nu, loc, sca, sig)

    use mo_kind,      only: dp

    real(dp),           intent(in) :: x
    real(dp),           intent(in) :: nu
    real(dp), optional, intent(in) :: loc
    real(dp), optional, intent(in) :: sca
    real(dp), optional, intent(in) :: sig
    real(dp)                       :: t_dp

    real(dp) :: iloc, isca

    iloc = 0.0_dp
    if (present(loc)) iloc = loc
    isca = 1.0_dp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig * sqrt((nu-2.0_dp)/nu)

    t_dp = t01((x-iloc)/isca, nu)/isca

  end function t_dp

  elemental pure function t_sp(x, nu, loc, sca, sig)

    use mo_kind,      only: sp

    real(sp),           intent(in) :: x
    real(sp),           intent(in) :: nu
    real(sp), optional, intent(in) :: loc
    real(sp), optional, intent(in) :: sca
    real(sp), optional, intent(in) :: sig
    real(sp)                       :: t_sp

    real(sp) :: iloc, isca

    iloc = 0.0_sp
    if (present(loc)) iloc = loc
    isca = 1.0_sp
    if (present(sca)) isca = sca
    if (present(sig)) isca = sig * sqrt((nu-2.0_sp)/nu)

    t_sp = t01((x-iloc)/isca, nu)/isca

  end function t_sp

  ! ------------------------------------------------------------------

  elemental pure function t01_dp(x, nu)

    use mo_kind,      only: dp
    use mo_constants, only: pi_dp
    use mo_functions, only: gamm

    real(dp), intent(in) :: x
    real(dp), intent(in) :: nu
    real(dp)             :: t01_dp

    real(dp) :: c

    c = gamm(0.5_dp*(nu+1.0_dp)) / (gamm(0.5_dp*nu) * sqrt(pi_dp*nu))
    t01_dp = c * (1.0_dp + x*x/nu)**(-0.5_dp*(nu+1.0_dp))

  end function t01_dp

  elemental pure function t01_sp(x, nu)

    use mo_kind,      only: sp
    use mo_constants, only: pi_sp
    use mo_functions, only: gamm

    real(sp), intent(in) :: x
    real(sp), intent(in) :: nu
    real(sp)             :: t01_sp

    real(sp) :: c

    c = gamm(0.5_sp*(nu+1.0_sp)) / (gamm(0.5_sp*nu) * sqrt(pi_sp*nu))
    t01_sp = c * (1.0_sp + x*x/nu)**(-0.5_sp*(nu+1.0_sp))

  end function t01_sp

  ! ------------------------------------------------------------------

end module mo_distributions

!> \file mo_kernel.f90

!> \brief   Module for kernel regression and kernel density estimation.

!> \details This module provides routines for kernel regression of data as well as kernel density
!>          estimation of both probability density functions (PDF) and cumulative density functions (CDF).\n
!>          So far only a Gaussian kernel is implemented (Nadaraya-Watson)
!>          which can be used for one- and multidimensional data.\n
!>          Furthermore, the estimation of the bandwith needed for kernel methods is provided
!>          by either Silverman''s rule of thumb or a Cross-Vaildation approach.\n
!>          The Cross-Validation method is actually an optimization of the bandwith and
!>          might be the most costly part of the kernel smoother.
!>          Therefore, the bandwith estimation is not necessarily part of the kernel smoothing
!>          but can be determined first and given as an optional argument to the smoother.

!> \author Juliane Mai
!> \date Mar 2013

MODULE mo_kernel

  ! This module provides functions for kernel regression and kernel density estimation and
  ! is part of the UFZ CHS Fortran library.

  ! Written  Juliane Mai,    Mar 2013
  ! Modified Stephan Thober, Mar 2013
  !          Matthias Cuntz, Mar 2013
  !          Matthias Cuntz, May 2014 - sort -> qsort

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! It is NOT released under the GNU Lesser General Public License, yet.

  ! If you use this routine, please contact Juliane Mai.

  ! Copyright 2013-2014 Juliane Mai, Stephan Thober, Matthias Cuntz

  USE mo_kind,      ONLY: i4, sp, dp
  USE mo_constants, ONLY: twopi_sp, twopi_dp
  USE mo_moment,    ONLY: stddev
  USE mo_nelmin,    ONLY: nelminrange      ! ND optimization
  USE mo_nr,        ONLY: golden           ! 1D optimization
  USE mo_quicksort, ONLY: qsort_index
  USE mo_integrate, ONLY: int_regular

  IMPLICIT NONE

  PUBLIC :: kernel_cumdensity        ! Kernel smoothing of a CDF                 (only 1D)
  PUBLIC :: kernel_density           ! Kernel smoothing of a PDF                 (only 1D)
  PUBLIC :: kernel_density_h         ! Bandwith estimation for PDF and CDF       (only 1D)
  PUBLIC :: kernel_regression        ! Kernel regression                         (1D and ND)
  PUBLIC :: kernel_regression_h      ! Bandwith estimation for kernel regression (1D and ND)

  ! ------------------------------------------------------------------

  !     NAME
  !         kernel_cumdensity

  !     PURPOSE
  !         Approximates the cumulative density function (CDF) to a given 1D data set using a Gaussian kernel.
  !
  !>        \brief   Approximates the cumulative density function (CDF).
  !
  !>        \details Approximates the cumulative density function (CDF)
  !>                 to a given 1D data set using a Gaussian kernel.\n
  !
  !>        The bandwith of the kernel can be pre-determined using the function kernel_density_h and
  !>        specified by the optional argument h.
  !>        If h is not given the default method to approximate the bandwith h is Silverman''s rule-of-thumb
  !>               \f[ h = \frac{4}{3}^{0.2} n^{-0.2} \sigma_x \f]
  !>        where n is the number of given data points and \f$ \sigma \f$ is the standard deviation of the data.\n
  !>        If the optional argument silverman is set to false, the cross-validation method described
  !>        by Scott et al. (2005) is applied.
  !>        Therefore, the likelihood of a given h is maximized using the Nelder-Mead algorithm nelminrange.
  !>        For large data sets this might be time consuming and should be performed aforehand using the
  !>        function kernel_density_h.\n
  !>        The dataset x can be single or double precision. The result will have the same numerical precision.\n
  !>        If the CDF for other datapoints than x is needed the optional argument xout has to be specified.
  !>        The result will than be of the same size and precision as xout.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x(:)"        \f$ x_i \f$ 1D-array with data points
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "real(sp/dp), optional :: h"       Bandwith of kernel.\n
  !>                                                     If present, argument silverman is ignored.
  !>                                                     If not present, the bandwith will be approximated first.
  !>       \param[in] "logical, optional :: silverman"   By default Silverman''s Rule-of-thumb will be used to approximate
  !>                                                     the bandwith of the kernel (silverman=true).
  !>                                                     If silverman=false the Cross-Validation approach is used
  !>                                                     to estimate the bandwidth.
  !>       \param[in] "real(sp/dp), optional :: xout(:)" If present, the CDF will be approximated at this arguments,
  !>                                                     otherwise the CDF is approximated at x.
  !>       \param[in] "integer(i4), optional :: nintegrate" If present, number of sampling points between for
  !>                                                        integration between output points.
  !>                                                        Should 1 plus a multiple of 4. Default: 101.
  !>       \param[in] "logical, optional :: mask(:)"     mask x values at calculation.\n
  !>                                                     if not xout given, then kernel estimates will have nodata value.
  !>       \param[in] "real(sp/dp), optional :: nodata"  if mask and not xout given, then masked data will
  !>                                                     have nodata kernel estimate.
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>       \return     real(sp/dp), allocatable :: kernel_cumdensity(:) &mdash; smoothed CDF at either x or xout
  !
  !     RESTRICTIONS
  !>       \note Data need to be one-dimensional. Multi-dimensional data handling not implemented yet.
  !
  !     EXAMPLE
  !         ! given data, e.g. temperature
  !         x = (/ 26.1_dp, 24.5_dp, 24.8_dp, 24.5_dp, 24.1_dp /)
  !
  !         ! estimate bandwidth via cross-validation (time-consuming)
  !         h = kernel_density_h(x,silverman=.false.)
  !
  !         ! estimate bandwidth with Silverman''s rule of thumb (default)
  !         h = kernel_density_h(x,silverman=.true.)
  !
  !         ! calc cumulative density with the estimated bandwidth h at given output points xout
  !         cdf = kernel_cumdensity(x, h=h, xout=xout)
  !         ! gives cumulative density at xout values, if specified, or at x values, if xout is not present
  !         ! if bandwith h is given                 : silverman (true/false) is ignored
  !         ! if silverman=.true.  and h not present : bandwith will be estimated using Silerman''s rule of thumb
  !         ! if silverman=.false. and h not present : bandwith will be estimated using Cross-Validation approach
  !         -> see also example in test directory

  !     LITERATURE
  !         Scott, D. W., & Sain, S. R. (2005).
  !             Multi-dimensional Density Estimation. Handbook of Statistics, 24, 229-261.
  !             doi:10.1016/S0169-7161(04)24009-3
  !         Haerdle, W., & Mueller, M. (2000). Multivariate and semiparametric kernel regression.
  !             In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
  !             application (pp. 357-392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12s,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Mar 2013
  !         Modified, Matthias Cuntz, Mar 2013
  !                   Matthias Cuntz, May 2014 - sort -> qsort

  INTERFACE kernel_cumdensity
     MODULE PROCEDURE kernel_cumdensity_1d_dp, kernel_cumdensity_1d_sp
  END INTERFACE kernel_cumdensity

  ! ------------------------------------------------------------------

  !     NAME
  !         kernel_density

  !     PURPOSE
  !         Approximates the probability density function (PDF) to a given 1D data set using a Gaussian kernel.
  !
  !>        \brief   Approximates the probability density function (PDF).
  !
  !>        \details Approximates the probability density function (PDF)
  !>                 to a given 1D data set using a Gaussian kernel.\n
  !
  !>        The bandwith of the kernel can be pre-determined using the function kernel_density_h and specified
  !>        by the optional argument h.
  !>        If h is not given the default method to approximate the bandwith h is Silverman''s rule-of-thumb
  !>               \f[ h = \frac{4}{3}^{0.2} n^{-0.2} \sigma_x \f]
  !>        where n is the number of given data points and \f$ \sigma \f$ is the standard deviation of the data.\n
  !>        If the optional argument silverman is set to false, the cross-validation method described
  !>        by Scott et al. (2005) is applied.
  !>        Therefore, the likelihood of a given h is maximized using the Nelder-Mead algorithm nelminrange.
  !>        For large data sets this might be time consuming and should be performed aforehand using the function
  !>        kernel_density_h.\n
  !>        The dataset x can be single or double precision. The result will have the same numerical precision.\n
  !>        If the PDF for other datapoints than x is needed the optional argument xout has to be specified.
  !>        The result will than be of the same size and precision as xout.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x(:)"        \f$ x_i \f$ 1D-array with data points
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "real(sp/dp), optional :: h"       Bandwith of kernel.\n
  !>                                                     If present, argument silverman is ignored.
  !>                                                     If not present, the bandwith will be approximated first.
  !>       \param[in] "logical, optional :: silverman"   By default Silverman''s Rule-of-thumb will be used to approximate
  !>                                                     the bandwith of the kernel (silverman=true).
  !>                                                     If silverman=false the Cross-Validation approach is used
  !>                                                     to estimate the bandwidth.
  !>       \param[in] "real(sp/dp), optional :: xout(:)" If present, the PDF will be approximated at this arguments,
  !>                                                     otherwise the PDF is approximated at x.
  !>       \param[in] "logical, optional :: mask(:)"     mask x values at calculation.\n
  !>                                                     if not xout given, then kernel estimates will have nodata value.
  !>       \param[in] "real(sp/dp), optional :: nodata"  if mask and not xout given, then masked data will
  !>                                                     have nodata kernel estimate.
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>       \return     real(sp/dp), allocatable :: kernel_density(:) &mdash; smoothed PDF at either x or xout
  !
  !     RESTRICTIONS
  !>       \note Data need to be one-dimensional. Multi-dimensional data handling not implemented yet.
  !
  !     EXAMPLE
  !         ! given data, e.g. temperature
  !         x = (/ 26.1_dp, 24.5_dp, 24.8_dp, 24.5_dp, 24.1_dp /)
  !
  !         ! estimate bandwidth via cross-validation (time-consuming)
  !         h = kernel_density_h(x,silverman=.false.)
  !
  !         ! estimate bandwidth with Silverman''s rule of thumb (default)
  !         h = kernel_density_h(x,silverman=.true.)
  !
  !         ! calc cumulative density with the estimated bandwidth h at given output points xout
  !         pdf = kernel_density(x, h=h, xout=xout)
  !         ! gives (probability) density at either xout values, if specified, or at x values, if xout is not present
  !         ! if bandwith h is given                 : silverman (true/false) is ignored
  !         ! if silverman=.true.  and h not present : bandwith will be estimated using Silerman''s rule of thumb
  !         ! if silverman=.false. and h not present : bandwith will be estimated using Cross-Validation approach
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         Scott, D. W., & Sain, S. R. (2005).
  !             Multi-dimensional Density Estimation. Handbook of Statistics, 24, 229-261.
  !             doi:10.1016/S0169-7161(04)24009-3
  !         Haerdle, W., & Mueller, M. (2000). Multivariate and semiparametric kernel regression.
  !             In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
  !             application (pp. 357-392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12s,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Mar 2013
  !         Modified, Stephan Thober, Mar 2013 - mask and nodata
  !                   Matthias Cuntz, Mar 2013

  INTERFACE kernel_density
     MODULE PROCEDURE kernel_density_1d_dp,  kernel_density_1d_sp
  END INTERFACE kernel_density

  ! ------------------------------------------------------------------

  !     NAME
  !         kernel_density_h

  !     PURPOSE
  !         Approximates the bandwith h of the kernel.
  !
  !>        \brief   Approximates the bandwith h of the kernel.
  !
  !>        \details  Approximates the bandwith h of the kernel for a given dataset x
  !>                  either using Silverman''s rule-of-thumb or a cross-validation method.\n
  !
  !>        By default the bandwidth h is approximated by Silverman''s rule-of-thumb
  !>               \f[ h = \frac{4}{3}^{0.2} n^{-0.2} \sigma_x \f]
  !>        where n is the number of given data points and \f$ \sigma \f$ is the standard deviation of the data.\n
  !>        If the optional argument silverman is set to false, the cross-validation method described
  !>        by Scott et al. (2005) is applied.
  !>        Therefore, the likelihood of a given h is maximized using the Nelder-Mead algorithm nelminrange.
  !>        For large data sets this might be time consuming and should be performed aforehand using the
  !>        function kernel_density_h.\n
  !>        The dataset x can be single or double precision. The result will have the same numerical precision.\n
  !>        The result of this function can be used as the optional input for kernel_density and kernel_cumdensity.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x(:)"        \f$ x_i \f$ 1D-array with data points
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: silverman"   By default Silverman''s Rule-of-thumb will be used to approximate
  !>                                                     the bandwith of the kernel (silverman=true).
  !>                                                     If silverman=false the Cross-Validation approach is used
  !>                                                     to estimate the bandwidth.
  !>       \param[in] "logical, optional :: mask(:)"     mask x values at calculation.
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>       \return     real(sp/dp), allocatable :: kernel_density_h(:) &mdash; approximated bandwith h for kernel smoother
  !
  !     RESTRICTIONS
  !>       \note Data need to be one-dimensional. Multi-dimensional data handling not implemented yet.
  !
  !     EXAMPLE
  !         ! given data, e.g. temperature
  !         x = (/ 26.1_dp, 24.5_dp, 24.8_dp, 24.5_dp, 24.1_dp /)
  !
  !         ! estimate bandwidth via cross-validation (time-consuming)
  !         h = kernel_density_h(x,silverman=.false.)
  !
  !         ! estimate bandwidth with Silverman''s rule of thumb (default)
  !         h = kernel_density_h(x,silverman=.true.)
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         Scott, D. W., & Sain, S. R. (2005).
  !             Multi-dimensional Density Estimation. Handbook of Statistics, 24, 229-261.
  !             doi:10.1016/S0169-7161(04)24009-3
  !         Haerdle, W., & Mueller, M. (2000). Multivariate and semiparametric kernel regression.
  !             In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
  !             application (pp. 357-392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12s,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Mar 2013
  !         Modified, Matthias Cuntz, Mar 2013

  INTERFACE kernel_density_h
     MODULE PROCEDURE kernel_density_h_1d_dp, kernel_density_h_1d_sp
  END INTERFACE kernel_density_h

  ! ------------------------------------------------------------------

  !     NAME
  !         kernel_regression

  !     PURPOSE
  !         Multi-dimensional non-parametric kernel regression using a Gaussian kernel.
  !
  !>        \brief   Multi-dimensional non-parametric kernel regression.
  !
  !>        \details Multi-dimensional non-parametric kernel regression using a Gaussian kernel.\n
  !
  !>        The bandwith of the kernel can be pre-determined using the function kernel_regression_h and specified
  !>        by the optional argument h.
  !>        If h is not given the default method to approximate the bandwith h is Silverman''s rule-of-thumb
  !>               \f[ h = \frac{4}{d+2}^{1/(d+4)} n^{-1/(d+4)} \sigma_{x_i} \f]
  !>        where \f$ n \f$ is the number of given data points, \f$ d \f$ is the number of dimensions,
  !>        and \f$ \sigma_{x_i} \f$ is the standard deviation of the data of dimension \f$ i \f$.\n
  !>        If the optional argument silverman is set to false, the cross-validation method described
  !>        by Scott et al. (2005) is applied.
  !>        Therefore, the likelihood of a given h is maximized using the Nelder-Mead algorithm nelminrange.
  !>        For large data sets this might be time consuming and should be performed aforehand
  !>        using the function kernel_regression_h.\n
  !>        The dataset (x,y) can be single or double precision. The result will have the same numerical precision.\n
  !>        If the regression for other datapoints than x is needed the optional argument xout has to be specified.
  !>        The result will than be of the same size and precision as xout.\n
  !>        \n
  !>        The code is adapted from the kernel_regression.py of the UFZ CHS Python library.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x(:)/x(:,:)"  1D or ND array with independent variables
  !>        \param[in] "real(sp/dp) :: y(:)"         1D array of dependent variables y(i) = f(x(i:)) with unknown f
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "real(sp/dp), optional :: h"       Bandwith of kernel.\n
  !>                                                     If present, argument silverman is ignored.
  !>                                                     If not present, the bandwith will be approximated first.
  !>       \param[in] "logical, optional :: silverman"   By default Silverman''s Rule-of-thumb will be used to approximate the
  !>                                                     bandwith of the kernel (silverman=true).
  !>                                                     If silverman=false the Cross-Validation approach is used
  !>                                                     to estimate the bandwidth.
  !>       \param[in] "real(sp/dp), optional :: xout(:)/xout(:,:)"
  !>                                                     If present, the fitted values will be returned for
  !>                                                     this independent variables,
  !>                                                     otherwise the fitted values at x are returned.
  !>       \param[in] "logical, optional :: mask(:)"     mask y values at calculation.\n
  !>                                                     if not xout given, then kernel estimates will have nodata value.
  !>       \param[in] "real(sp/dp), optional :: nodata"  if mask and not xout given, then masked data will
  !>                                                     have nodata kernel estimate.
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     real(sp/dp), allocatable :: kernel_regression(:) &mdash; fitted values at eighter x or xout
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         ! given data, e.g. temperature
  !         x(:,1) = (/ 26.1_dp, 24.5_dp, 24.8_dp, 14.5_dp,  2.1_dp /)
  !         x(:,2) = (/  2.1_dp,  4.5_dp,  6.8_dp,  4.8_dp,  0.1_dp /)
  !         y      = (/ 28.2_dp, 29.0_dp, 31.6_dp, 19.3_dp,  2.2_dp /)
  !
  !         ! estimate bandwidth via cross-validation (time-consuming)
  !         h = kernel_regression_h(x,y,silverman=.false.)
  !
  !         ! estimate bandwidth with Silverman''s rule of thumb (default)
  !         h = kernel_regression_h(x,y,silverman=.true.)
  !
  !         fit_y = kernel_regression(x, y, h=h, silverman=.false., xout=xout)
  !         ! gives fitted values at either xout values, if specified, or at x values, if xout is not present
  !         ! if bandwith h is given                 : silverman (true/false) is ignored
  !         ! if silverman=.true.  and h not present : bandwith will be estimated using Silerman''s rule of thumb
  !         ! if silverman=.false. and h not present : bandwith will be estimated using Cross-Validation approach
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         Scott, D. W., & Sain, S. R. (2005).
  !             Multi-dimensional Density Estimation. Handbook of Statistics, 24, 229-261.
  !             doi:10.1016/S0169-7161(04)24009-3
  !         Haerdle, W., & Mueller, M. (2000). Multivariate and semiparametric kernel regression.
  !             In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
  !             application (pp. 357-392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12s,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \authors Matthias Cuntz, Juliane Mai
  !>        \date Mar 2013
  INTERFACE kernel_regression
     MODULE PROCEDURE kernel_regression_2d_dp, kernel_regression_2d_sp, &
          kernel_regression_1d_dp, kernel_regression_1d_sp
  END INTERFACE kernel_regression

  ! ------------------------------------------------------------------

  !     NAME
  !         kernel_regression_h

  !     PURPOSE
  !         Approximates the bandwith h of the kernel.
  !
  !>        \brief   Approximates the bandwith h of the kernel for regression.
  !
  !>        \details  Approximates the bandwith h of the kernel for a given dataset x
  !>                  either using Silverman''s rule-of-thumb or a cross-validation method.\n
  !
  !>        By default the bandwidth h is approximated by Silverman''s rule-of-thumb
  !>               \f[ h = \frac{4}{d+2}^{1/(d+4)} n^{-1/(d+4)} \sigma_{x_i} \f]
  !>        where \f$ n \f$ is the number of given data points, \f$ d \f$ is the number of dimensions,
  !>        and \f$ \sigma_{x_i} \f$ is the standard deviation of the data of dimension \f$ i \f$.\n
  !>        If the optional argument silverman is set to false, the cross-validation method described
  !>        by Scott et al. (2005) is applied.
  !>        Therefore, the likelihood of a given h is maximized using the Nelder-Mead algorithm nelminrange.
  !>        For large data sets this might be time consuming and should be performed aforehand using the function kernel_density_h.\n
  !>        The dataset x can be single or double precision. The result will have the same numerical precision.\n
  !>        The result of this function can be used as the optional input for kernel_regression.\n
  !>        \n
  !>        The code is adapted from the kernel_regression.py of the UFZ CHS Python library.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: x(:)/x(:,:)"  1D or ND array with independent variables
  !>        \param[in] "real(sp/dp) :: y(:)"         1D array of dependent variables y(i) = f(x(i:)) with unknown f
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: silverman"   By default Silverman''s Rule-of-thumb will be used to approximate the
  !>                                                     bandwith of the kernel (silverman=true).
  !>                                                     If silverman=false the Cross-Validation approach is used
  !>                                                     to estimate the bandwidth.
  !>       \param[in] "logical, optional :: mask(:)"     mask x values at calculation.
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     real(sp/dp), allocatable :: kernel_regression_h(:) &mdash; approximated bandwith h for kernel regression\n
  !>                                                                               number of bandwidths equals
  !>                                                                               number of independent variables
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         ! given data, e.g. temperature
  !         x(:,1) = (/ 26.1_dp, 24.5_dp, 24.8_dp, 14.5_dp,  2.1_dp /)
  !         x(:,2) = (/  2.1_dp,  4.5_dp,  6.8_dp,  4.8_dp,  0.1_dp /)
  !         y      = (/ 28.2_dp, 29.0_dp, 31.6_dp, 19.3_dp,  2.2_dp /)
  !
  !         ! estimate bandwidth via cross-validation (time-consuming)
  !         h = kernel_regression_h(x,y,silverman=.false.)
  !
  !         ! estimate bandwidth with Silverman''s rule of thumb (default)
  !         h = kernel_regression_h(x,y,silverman=.true.)
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         Scott, D. W., & Sain, S. R. (2005).
  !             Multi-dimensional Density Estimation. Handbook of Statistics, 24, 229-261.
  !             doi:10.1016/S0169-7161(04)24009-3
  !         Haerdle, W., & Mueller, M. (2000). Multivariate and semiparametric kernel regression.
  !             In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
  !             application (pp. 357-392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12s,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \authors Matthias Cuntz, Juliane Mai
  !>        \date Mar 2013
  INTERFACE kernel_regression_h
     MODULE PROCEDURE kernel_regression_h_2d_dp, kernel_regression_h_2d_sp, &
          kernel_regression_h_1d_dp, kernel_regression_h_1d_sp
  END INTERFACE kernel_regression_h

  ! ------------------------------------------------------------------

  INTERFACE nadaraya_watson
     MODULE PROCEDURE nadaraya_watson_2d_dp, nadaraya_watson_2d_sp, &
          nadaraya_watson_1d_dp, nadaraya_watson_1d_sp
  END INTERFACE nadaraya_watson

  INTERFACE allocate_globals
     MODULE PROCEDURE allocate_globals_2d_dp, allocate_globals_2d_sp, &
          allocate_globals_1d_dp, allocate_globals_1d_sp
  END INTERFACE allocate_globals

  INTERFACE cross_valid_density
     MODULE PROCEDURE cross_valid_density_1d_dp, cross_valid_density_1d_sp
  END INTERFACE cross_valid_density

  INTERFACE cross_valid_regression
     MODULE PROCEDURE cross_valid_regression_dp, cross_valid_regression_sp
  END INTERFACE cross_valid_regression

  INTERFACE mesh
     MODULE PROCEDURE mesh_dp, mesh_sp
  END INTERFACE mesh

  PRIVATE

  ! Module variables which need to be public for optimization of bandwith via cross-validation
  real(sp), dimension(:,:), allocatable :: global_x_sp
  real(sp), dimension(:,:), allocatable :: global_xout_sp
  real(sp), dimension(:),   allocatable :: global_y_sp

  real(dp), dimension(:,:), allocatable :: global_x_dp
  real(dp), dimension(:,:), allocatable :: global_xout_dp
  real(dp), dimension(:),   allocatable :: global_y_dp

  ! ------------------------------------------------------------------------------------------------

CONTAINS

  function kernel_cumdensity_1d_dp(ix, h, silverman, xout, nintegrate, mask, nodata)

    implicit none

    real(dp), dimension(:),                       intent(in) :: ix
    real(dp),                           optional, intent(in) :: h
    logical,                            optional, intent(in) :: silverman
    real(dp), dimension(:),             optional, intent(in) :: xout
    integer(i4),                        optional, intent(in) :: nintegrate
    logical,  dimension(:),             optional, intent(in) :: mask
    real(dp),                           optional, intent(in) :: nodata
    real(dp), dimension(:), allocatable                      :: kernel_cumdensity_1d_dp

    ! local variables
    integer(i4)                            :: nin, nout, nmesh
    integer(i4)                            :: ii
    real(dp)                               :: hh
    real(dp)                               :: lower_x
    real(dp),    dimension(:), allocatable :: xxout
    integer(i4), dimension(:), allocatable :: xindx
    real(dp),    dimension(:), allocatable :: kernel_pdf
    real(dp),    dimension(:), allocatable :: xmesh
    real(dp)                               :: delta
    ! real(dp)                               :: tmp
    real(dp),    dimension(:), allocatable :: z
    real(dp),    dimension(:), allocatable :: x

    ! consistency check - mask needs either nodata or xout
    if (present(mask) .and. (.not. present(xout)) .and. (.not. present(nodata)) ) then
       stop 'kernel_cumdensity_1d_dp: missing nodata value or xout with present(mask)'
    end if

    ! Pack if mask present
    if (present(mask)) then
       nin = count(mask)
       allocate(x(nin))
       x = pack(ix, mask)
    else
       nin = size(ix,1)
       allocate(x(nin))
       x = ix
    endif
    allocate(z(nin))

    ! allocate
    if (present(xout)) then
       nout = size(xout,1)
       allocate(xxout(nout))
       allocate(xindx(nout))
       xxout = xout
    else
       nout = nin
       allocate(xxout(nout))
       allocate(xindx(nout))
       xxout = x
    end if
    ! sort the x
    xindx = qsort_index(xxout)
    xxout = xxout(xindx)

    ! should be (n*4 + 1) for int_regular
    if (present(nintegrate)) then
       nmesh = nintegrate
    else
       nmesh = 101_i4
    endif

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_density_h(x, silverman=silverman)
       else
          hh = kernel_density_h(x, silverman=.true.)
       end if
    end if

    ! allocate PDF, mesh and CDF
    allocate(kernel_pdf(nmesh))
    allocate(xmesh(nmesh))
    allocate(kernel_cumdensity_1d_dp(nout))

    ! calculate standard deviation of x to determine left side starting point for integration of PDF
    lower_x = minval(x) - 3.0_dp * stddev(x)

    ! loop through all regression points
    do ii = 1, nout
       ! generate nmesh points between last x and this x
       ! integrate pdf and add to last point
       if (ii .eq. 1_i4) then
          xmesh                       = mesh(lower_x, xxout(1), nmesh, delta)
          kernel_pdf(:)               = kernel_density(x, hh, xout=xmesh)
          kernel_cumdensity_1d_dp(1)  = int_regular(kernel_pdf, delta)
       else
          xmesh                       = mesh(xxout(ii-1), xxout(ii), nmesh, delta)
          kernel_pdf(:)               = kernel_density(x, hh, xout=xmesh)
          kernel_cumdensity_1d_dp(ii) = kernel_cumdensity_1d_dp(ii-1) + int_regular(kernel_pdf, delta)
       end if
    end do

    ! ! scale to range [0,1]
    ! tmp = 1.0_dp / (kernel_cumdensity_1d_dp(nout) - kernel_cumdensity_1d_dp(1))
    ! kernel_cumdensity_1d_dp(:) = ( kernel_cumdensity_1d_dp(:) - kernel_cumdensity_1d_dp(1) ) * tmp
    kernel_cumdensity_1d_dp = min(kernel_cumdensity_1d_dp, 1.0_dp)

    ! resorting
    kernel_cumdensity_1d_dp(xindx) = kernel_cumdensity_1d_dp(:)

    ! check whether output has to be unpacked
    if (present(mask) .and. (.not. present(xout))) then
       deallocate(x)
       nin = size(ix,1)
       allocate(x(nin))
       x = unpack(kernel_cumdensity_1d_dp, mask, nodata)
       deallocate(kernel_cumdensity_1d_dp)
       allocate(kernel_cumdensity_1d_dp(nin))
       kernel_cumdensity_1d_dp = x
    end if

    ! clean up
    deallocate(xxout)
    deallocate(xindx)
    deallocate(kernel_pdf)
    deallocate(xmesh)
    deallocate(z)
    deallocate(x)

  end function kernel_cumdensity_1d_dp

  function kernel_cumdensity_1d_sp(ix, h, silverman, xout, nintegrate, mask, nodata)

    implicit none

    real(sp), dimension(:),                       intent(in) :: ix
    real(sp),                           optional, intent(in) :: h
    logical,                            optional, intent(in) :: silverman
    real(sp), dimension(:),             optional, intent(in) :: xout
    integer(i4),                        optional, intent(in) :: nintegrate
    logical,  dimension(:),             optional, intent(in) :: mask
    real(sp),                           optional, intent(in) :: nodata
    real(sp), dimension(:), allocatable                      :: kernel_cumdensity_1d_sp

    ! local variables
    integer(i4)                            :: nin, nout, nmesh
    integer(i4)                            :: ii
    real(sp)                               :: hh
    real(sp)                               :: lower_x
    real(sp),    dimension(:), allocatable :: xxout
    integer(i4), dimension(:), allocatable :: xindx
    real(sp),    dimension(:), allocatable :: kernel_pdf
    real(sp),    dimension(:), allocatable :: xmesh
    real(sp)                               :: delta
    ! real(sp)                               :: tmp
    real(sp),    dimension(:), allocatable :: z
    real(sp),    dimension(:), allocatable :: x

    ! consistency check - mask needs either nodata or xout
    if (present(mask) .and. (.not. present(xout)) .and. (.not. present(nodata)) ) then
       stop 'kernel_cumdensity_1d_sp: missing nodata value or xout with present(mask)'
    end if

    ! Pack if mask present
    if (present(mask)) then
       nin = count(mask)
       allocate(x(nin))
       x = pack(ix, mask)
    else
       nin = size(ix,1)
       allocate(x(nin))
       x = ix
    endif
    allocate(z(nin))

    ! allocate
    if (present(xout)) then
       nout = size(xout,1)
       allocate(xxout(nout))
       allocate(xindx(nout))
       xxout = xout
    else
       nout = nin
       allocate(xxout(nout))
       allocate(xindx(nout))
       xxout = x
    end if
    ! sort the x
    xindx = qsort_index(xxout)
    xxout = xxout(xindx)

    ! should be (n*4 + 1) for int_regular
    if (present(nintegrate)) then
       nmesh = nintegrate
    else
       nmesh = 101_i4
    endif

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_density_h(x, silverman=silverman)
       else
          hh = kernel_density_h(x, silverman=.true.)
       end if
    end if

    ! allocate PDF, mesh and CDF
    allocate(kernel_pdf(nmesh))
    allocate(xmesh(nmesh))
    allocate(kernel_cumdensity_1d_sp(nout))

    ! calculate standard deviation of x to determine left side starting point for integration of PDF
    lower_x = minval(x) - 3.0_sp * stddev(x)

    ! loop through all regression points
    do ii = 1, nout
       ! generate nmesh points between last x and this x
       ! integrate pdf and add to last point
       if (ii .eq. 1_i4) then
          xmesh                       = mesh(lower_x, xxout(1), nmesh, delta)
          kernel_pdf(:)               = kernel_density(x, hh, xout=xmesh)
          kernel_cumdensity_1d_sp(1)  = int_regular(kernel_pdf, delta)
       else
          xmesh                       = mesh(xxout(ii-1), xxout(ii), nmesh, delta)
          kernel_pdf(:)               = kernel_density(x, hh, xout=xmesh)
          kernel_cumdensity_1d_sp(ii) = kernel_cumdensity_1d_sp(ii-1) + int_regular(kernel_pdf, delta)
       end if
    end do

    ! ! scale to range [0,1]
    ! tmp = 1.0_sp / (kernel_cumdensity_1d_sp(nout) - kernel_cumdensity_1d_sp(1))
    ! kernel_cumdensity_1d_sp(:) = ( kernel_cumdensity_1d_sp(:) - kernel_cumdensity_1d_sp(1) ) * tmp
    kernel_cumdensity_1d_sp = min(kernel_cumdensity_1d_sp, 1.0_sp)

    ! resorting
    kernel_cumdensity_1d_sp(xindx) = kernel_cumdensity_1d_sp(:)

    ! check whether output has to be unpacked
    if (present(mask) .and. (.not. present(xout))) then
       deallocate(x)
       nin = size(ix,1)
       allocate(x(nin))
       x = unpack(kernel_cumdensity_1d_sp, mask, nodata)
       deallocate(kernel_cumdensity_1d_sp)
       allocate(kernel_cumdensity_1d_sp(nin))
       kernel_cumdensity_1d_sp = x
    end if

    ! clean up
    deallocate(xxout)
    deallocate(xindx)
    deallocate(kernel_pdf)
    deallocate(xmesh)
    deallocate(z)
    deallocate(x)

  end function kernel_cumdensity_1d_sp

  ! ------------------------------------------------------------------------------------------------

  function kernel_density_1d_dp(ix, h, silverman, xout, mask, nodata)

    implicit none

    real(dp), dimension(:),                       intent(in) :: ix
    real(dp),                           optional, intent(in) :: h
    logical,                            optional, intent(in) :: silverman
    real(dp), dimension(:),             optional, intent(in) :: xout
    logical,  dimension(:),             optional, intent(in) :: mask
    real(dp),                           optional, intent(in) :: nodata
    real(dp), dimension(:), allocatable                      :: kernel_density_1d_dp

    ! local variables
    integer(i4)                         :: nin, nout
    integer(i4)                         :: ii
    real(dp)                            :: hh
    real(dp), dimension(:), allocatable :: xxout
    real(dp), dimension(:), allocatable :: z
    real(dp)                            :: multiplier
    real(dp)                            :: thresh
    real(dp), dimension(:),        allocatable :: x

    ! consistency check - mask needs either nodata or xout
    if (present(mask) .and. (.not. present(xout)) .and. (.not. present(nodata)) ) then
       stop 'kernel_density_1d_dp: missing nodata value or xout with present(mask)'
    end if

    ! Pack if mask present
    if (present(mask)) then
       nin = count(mask)
       allocate(x(nin))
       x = pack(ix, mask)
    else
       nin = size(ix,1)
       allocate(x(nin))
       x = ix
    endif
    allocate(z(nin))

    ! output size
    if (present(xout)) then
       nout = size(xout,1)
       allocate(xxout(nout))
       xxout = xout
    else
       nout = nin
       allocate(xxout(nout))
       xxout = x
    end if
    ! allocate output
    allocate(kernel_density_1d_dp(nout))

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_density_h(x, silverman=silverman)
       else
          hh = kernel_density_h(x, silverman=.true.)
       end if
    end if

    multiplier = 1.0_dp/(real(nin,dp)*hh)
    if (multiplier <= 1.0_dp) then
       thresh = tiny(1.0_dp)/multiplier
    else
       thresh = 0.0_dp
    endif
    ! loop through all regression points
    do ii=1, nout
       ! scaled difference from regression point
       z(:) = (x(:) - xxout(ii)) / hh
       ! nadaraya-watson estimator of gaussian multivariate kernel
       kernel_density_1d_dp(ii) = nadaraya_watson(z)
       if (kernel_density_1d_dp(ii) .gt. thresh) kernel_density_1d_dp(ii) = multiplier * kernel_density_1d_dp(ii)
    end do

    ! check whether output has to be unpacked
    if (present(mask) .and. (.not. present(xout))) then
       deallocate(x)
       nin = size(ix,1)
       allocate(x(nin))
       x = unpack(kernel_density_1d_dp, mask, nodata)
       deallocate(kernel_density_1d_dp)
       allocate(kernel_density_1d_dp(nin))
       kernel_density_1d_dp = x
    end if

    ! clean up
    deallocate(xxout)
    deallocate(x)
    deallocate(z)

  end function kernel_density_1d_dp

  function kernel_density_1d_sp(ix, h, silverman, xout, mask, nodata)

    implicit none

    real(sp), dimension(:),                       intent(in) :: ix
    real(sp),                           optional, intent(in) :: h
    logical,                            optional, intent(in) :: silverman
    real(sp), dimension(:),             optional, intent(in) :: xout
    logical,  dimension(:),             optional, intent(in) :: mask
    real(sp),                           optional, intent(in) :: nodata
    real(sp), dimension(:), allocatable                      :: kernel_density_1d_sp

    ! local variables
    integer(i4)                         :: nin, nout
    integer(i4)                         :: ii
    real(sp)                            :: hh
    real(sp), dimension(:), allocatable :: xxout
    real(sp), dimension(:), allocatable :: z
    real(sp)                            :: multiplier
    real(sp)                            :: thresh
    real(sp), dimension(:),        allocatable :: x

    ! consistency check - mask needs either nodata or xout
    if (present(mask) .and. (.not. present(xout)) .and. (.not. present(nodata)) ) then
       stop 'kernel_density_1d_sp: missing nodata value or xout with present(mask)'
    end if

    ! Pack if mask present
    if (present(mask)) then
       nin = count(mask)
       allocate(x(nin))
       x = pack(ix, mask)
    else
       nin = size(ix,1)
       allocate(x(nin))
       x = ix
    endif
    allocate(z(nin))

    ! output size
    if (present(xout)) then
       nout = size(xout,1)
       allocate(xxout(nout))
       xxout = xout
    else
       nout = nin
       allocate(xxout(nout))
       xxout = x
    end if
    ! allocate output
    allocate(kernel_density_1d_sp(nout))

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_density_h(x, silverman=silverman)
       else
          hh = kernel_density_h(x, silverman=.true.)
       end if
    end if

    multiplier = 1.0_sp/(real(nin,sp)*hh)
    if (multiplier <= 1.0_sp) then
       thresh = tiny(1.0_sp)/multiplier
    else
       thresh = 0.0_sp
    endif
    ! loop through all regression points
    do ii=1, nout
       ! scaled difference from regression point
       z(:) = (x(:) - xxout(ii)) / hh
       ! nadaraya-watson estimator of gaussian multivariate kernel
       kernel_density_1d_sp(ii) = nadaraya_watson(z)
       if (kernel_density_1d_sp(ii) .gt. thresh) kernel_density_1d_sp(ii) = multiplier * kernel_density_1d_sp(ii)
    end do

    ! check whether output has to be unpacked
    if (present(mask) .and. (.not. present(xout))) then
       deallocate(x)
       nin = size(ix,1)
       allocate(x(nin))
       x = unpack(kernel_density_1d_sp, mask, nodata)
       deallocate(kernel_density_1d_sp)
       allocate(kernel_density_1d_sp(nin))
       kernel_density_1d_sp = x
    end if

    ! clean up
    deallocate(xxout)
    deallocate(x)
    deallocate(z)

  end function kernel_density_1d_sp

  ! ------------------------------------------------------------------------------------------------

  function kernel_density_h_1d_dp(ix, silverman, mask)

    implicit none

    real(dp), dimension(:),           intent(in) :: ix
    logical,                optional, intent(in) :: silverman
    logical,  dimension(:), optional, intent(in) :: mask
    real(dp)                                     :: kernel_density_h_1d_dp

    ! local variables
    integer(i4)              :: nin
    real(dp)                 :: nn
    real(dp)                 :: h, hmin, fmin ! fmin set but never referenced
    real(dp), dimension(2)   :: bounds
    real(dp), parameter      :: pre_h = 1.05922384104881_dp
    real(dp), dimension(:), allocatable :: x

    if (present(mask)) then
       nin = count(mask)
       allocate(x(nin))
       x = pack(ix, mask)
    else
       nin = size(ix,1)
       allocate(x(nin))
       x = ix
    endif
    nn = real(nin,dp)

    ! Default: Silverman's rule of thumb by
    ! Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    !h(1) = (4._dp/3._dp/real(nn,dp))**(0.2_dp) * stddev_x
    h = pre_h/(nn**0.2_dp) * stddev(x(:))

    if (present(silverman)) then
       if (.not. silverman) then
          bounds(1) = max(0.2_dp * h, (maxval(x)-minval(x))/nn)
          bounds(2) = 5.0_dp * h
          call allocate_globals(x)
          fmin = golden(bounds(1),h,bounds(2),cross_valid_density_1d_dp, 0.0001_dp,hmin)
          fmin = fmin * 1.0_dp  ! just to avoid warnings of "Variable FMIN set but never referenced"
          h = hmin
          call deallocate_globals()
       end if
    end if

    kernel_density_h_1d_dp = h

    deallocate(x)

  end function kernel_density_h_1d_dp

  function kernel_density_h_1d_sp(ix, silverman, mask)

    implicit none

    real(sp), dimension(:),           intent(in) :: ix
    logical,                optional, intent(in) :: silverman
    logical,  dimension(:), optional, intent(in) :: mask
    real(sp)                                     :: kernel_density_h_1d_sp

    ! local variables
    integer(i4)              :: nin
    real(sp)                 :: nn
    real(sp)                 :: h, hmin, fmin ! fmin set but never referenced
    real(sp), dimension(2)   :: bounds
    real(sp), parameter      :: pre_h = 1.05922384104881_sp
    real(sp), dimension(:), allocatable :: x

    if (present(mask)) then
       nin = count(mask)
       allocate(x(nin))
       x = pack(ix, mask)
    else
       nin = size(ix,1)
       allocate(x(nin))
       x = ix
    endif
    nn = real(nin,sp)

    ! Default: Silverman's rule of thumb by
    ! Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    !h(1) = (4._sp/3._sp/real(nn,sp))**(0.2_sp) * stddev_x
    h = pre_h/(nn**0.2_sp) * stddev(x(:))

    if (present(silverman)) then
       if (.not. silverman) then
          bounds(1) = max(0.2_sp * h, (maxval(x)-minval(x))/nn)
          bounds(2) = 5.0_sp * h
          call allocate_globals(x)
          fmin = golden(bounds(1),h,bounds(2),cross_valid_density_1d_sp, 0.0001_sp,hmin)
          fmin = fmin * 1.0_sp  ! just to avoid warnings of "Variable FMIN set but never referenced"
          h = hmin
          call deallocate_globals()
       end if
    end if

    kernel_density_h_1d_sp = h

    deallocate(x)

  end function kernel_density_h_1d_sp

  ! ------------------------------------------------------------------------------------------------

  function kernel_regression_1d_dp(ix, iy, h, silverman, xout, mask, nodata)

    implicit none

    real(dp), dimension(:),                       intent(in) :: ix
    real(dp), dimension(:),                       intent(in) :: iy
    real(dp),                           optional, intent(in) :: h
    logical,                            optional, intent(in) :: silverman
    real(dp), dimension(:),             optional, intent(in) :: xout
    logical,  dimension(:),             optional, intent(in) :: mask
    real(dp),                           optional, intent(in) :: nodata
    real(dp), dimension(:), allocatable                      :: kernel_regression_1d_dp

    ! local variables
    integer(i4)                         :: nin, nout
    integer(i4)                         :: ii
    real(dp)                            :: hh
    real(dp), dimension(:), allocatable :: xxout
    real(dp), dimension(:), allocatable :: z
    real(dp), dimension(:), allocatable :: x
    real(dp), dimension(:), allocatable :: y

    ! consistency check - mask needs either nodata or xout
    if (present(mask) .and. (.not. present(xout)) .and. (.not. present(nodata)) ) then
       stop 'kernel_regression_1d_dp: missing nodata value or xout with present(mask)'
    end if

    nin   = size(ix,1)
    ! consistency checks of inputs
    if (size(iy,1) .ne. nin) stop 'kernel_regression_1d_dp: size(x) /= size(y)'

    ! Pack if mask present
    if (present(mask)) then
       nin = count(mask)
       allocate(x(nin))
       allocate(y(nin))
       x = pack(ix, mask)
       y = pack(iy, mask)
    else
       nin = size(ix,1)
       allocate(x(nin))
       allocate(y(nin))
       x = ix
       y = iy
    endif
    allocate(z(nin))

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_regression_h(x, y, silverman=silverman)
       else
          hh = kernel_regression_h(x, y, silverman=.true.)
       end if
    end if

    ! calc regression
    if (present(xout)) then
       nout = size(xout,1)
       allocate(xxout(nout))
       xxout = xout
    else
       nout = nin
       allocate(xxout(nout))
       xxout = x
    end if
    ! allocate output
    allocate(kernel_regression_1d_dp(nout))

    ! loop through all regression points
    do ii = 1, nout
       ! scaled difference from regression point
       z(:) = (x(:) - xxout(ii)) / hh
       ! nadaraya-watson estimator of gaussian multivariate kernel
       kernel_regression_1d_dp(ii) = nadaraya_watson(z, y)
    end do

    ! check whether output has to be unpacked
    if (present(mask) .and. (.not. present(xout))) then
       deallocate(x)
       nin = size(ix,1)
       allocate(x(nin))
       x = unpack(kernel_regression_1d_dp, mask, nodata)
       deallocate(kernel_regression_1d_dp)
       allocate(kernel_regression_1d_dp(nin))
       kernel_regression_1d_dp = x
    end if

    ! clean up
    deallocate(xxout)
    deallocate(x)
    deallocate(y)
    deallocate(z)

  end function kernel_regression_1d_dp

  function kernel_regression_1d_sp(ix, iy, h, silverman, xout, mask, nodata)

    implicit none

    real(sp), dimension(:),                       intent(in) :: ix
    real(sp), dimension(:),                       intent(in) :: iy
    real(sp),                           optional, intent(in) :: h
    logical,                            optional, intent(in) :: silverman
    real(sp), dimension(:),             optional, intent(in) :: xout
    logical,  dimension(:),             optional, intent(in) :: mask
    real(sp),                           optional, intent(in) :: nodata
    real(sp), dimension(:), allocatable                      :: kernel_regression_1d_sp

    ! local variables
    integer(i4)                         :: nin, nout
    integer(i4)                         :: ii
    real(sp)                            :: hh
    real(sp), dimension(:), allocatable :: xxout
    real(sp), dimension(:), allocatable :: z
    real(sp), dimension(:), allocatable :: x
    real(sp), dimension(:), allocatable :: y

    ! consistency check - mask needs either nodata or xout
    if (present(mask) .and. (.not. present(xout)) .and. (.not. present(nodata)) ) then
       stop 'kernel_regression_1d_sp: missing nodata value or xout with present(mask)'
    end if

    nin   = size(ix,1)
    ! consistency checks of inputs
    if (size(iy,1) .ne. nin) stop 'kernel_regression_1d_sp: size(x) /= size(y)'

    ! Pack if mask present
    if (present(mask)) then
       nin = count(mask)
       allocate(x(nin))
       allocate(y(nin))
       x = pack(ix, mask)
       y = pack(iy, mask)
    else
       nin = size(ix,1)
       allocate(x(nin))
       allocate(y(nin))
       x = ix
       y = iy
    endif
    allocate(z(nin))

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_regression_h(x, y, silverman=silverman)
       else
          hh = kernel_regression_h(x, y, silverman=.true.)
       end if
    end if

    ! calc regression
    if (present(xout)) then
       nout = size(xout,1)
       allocate(xxout(nout))
       xxout = xout
    else
       nout = nin
       allocate(xxout(nout))
       xxout = x
    end if
    ! allocate output
    allocate(kernel_regression_1d_sp(nout))

    ! loop through all regression points
    do ii = 1, nout
       ! scaled difference from regression point
       z(:) = (x(:) - xxout(ii)) / hh
       ! nadaraya-watson estimator of gaussian multivariate kernel
       kernel_regression_1d_sp(ii) = nadaraya_watson(z, y)
    end do

    ! check whether output has to be unpacked
    if (present(mask) .and. (.not. present(xout))) then
       deallocate(x)
       nin = size(ix,1)
       allocate(x(nin))
       x = unpack(kernel_regression_1d_sp, mask, nodata)
       deallocate(kernel_regression_1d_sp)
       allocate(kernel_regression_1d_sp(nin))
       kernel_regression_1d_sp = x
    end if

    ! clean up
    deallocate(xxout)
    deallocate(x)
    deallocate(y)
    deallocate(z)

  end function kernel_regression_1d_sp

  function kernel_regression_2d_dp(ix, iy, h, silverman, xout, mask, nodata)

    implicit none

    real(dp), dimension(:,:),                       intent(in) :: ix
    real(dp), dimension(:),                         intent(in) :: iy
    real(dp), dimension(:),               optional, intent(in) :: h
    logical,                              optional, intent(in) :: silverman
    real(dp), dimension(:,:),             optional, intent(in) :: xout
    logical,  dimension(:),               optional, intent(in) :: mask
    real(dp),                             optional, intent(in) :: nodata
    real(dp), dimension(:),   allocatable                      :: kernel_regression_2d_dp

    ! local variables
    integer(i4)                           :: dims
    integer(i4)                           :: nin, nout
    integer(i4)                           :: ii, jj
    real(dp), dimension(size(ix,2))       :: hh
    real(dp), dimension(:,:), allocatable :: xxout
    real(dp), dimension(:,:), allocatable :: z
    real(dp), dimension(:,:), allocatable :: x
    real(dp), dimension(:),   allocatable :: y

    ! consistency check - mask needs either nodata or xout
    if (present(mask) .and. (.not. present(xout)) .and. (.not. present(nodata)) ) then
       stop 'kernel_regression_1d_dp: missing nodata value or xout with present(mask)'
    end if

    nin  = size(ix,1)
    dims = size(ix,2)
    ! consistency checks of inputs
    if (size(iy) .ne. nin) stop 'kernel_regression_2d_dp: size(y) /= size(x,1)'
    if (present(h)) then
       if (size(h) .ne. dims) stop 'kernel_regression_2d_dp: size(h) /= size(x,2)'
    end if
    if (present(xout)) then
       if (size(xout,2) .ne. dims) stop 'kernel_regression_2d_dp: size(xout) /= size(x,2)'
    end if

    ! Pack if mask present
    if (present(mask)) then
       nin = count(mask)
       allocate(x(nin,dims))
       allocate(y(nin))
       forall(ii=1:dims) x(:,ii) = pack(ix(:,ii), mask)
       y = pack(iy, mask)
    else
       nin = size(ix,1)
       allocate(x(nin,dims))
       allocate(y(nin))
       x = ix
       y = iy
    endif
    allocate(z(nin,dims))

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_regression_h(x, y, silverman=silverman)
       else
          hh = kernel_regression_h(x, y, silverman=.true.)
       end if
    end if

    ! calc regression
    if (present(xout)) then
       nout = size(xout,1)
       allocate(xxout(nout,dims))
       xxout = xout
    else
       nout = nin
       allocate(xxout(nout,dims))
       xxout = x
    end if
    ! allocate output
    allocate(kernel_regression_2d_dp(nout))

    ! loop through all regression points
    do ii = 1, nout
       forall(jj=1:dims) z(:,jj) = (x(:,jj) - xxout(ii,jj)) / hh(jj)
       ! nadaraya-watson estimator of gaussian multivariate kernel
       kernel_regression_2d_dp(ii) = nadaraya_watson(z, y)
    end do

    ! check whether output has to be unpacked
    if (present(mask) .and. (.not. present(xout))) then
       deallocate(y)
       nin = size(iy,1)
       allocate(y(nin))
       y = unpack(kernel_regression_2d_dp, mask, nodata)
       deallocate(kernel_regression_2d_dp)
       allocate(kernel_regression_2d_dp(nin))
       kernel_regression_2d_dp = y
    end if

    ! clean up
    deallocate(xxout)
    deallocate(x)
    deallocate(y)
    deallocate(z)

  end function kernel_regression_2d_dp

  function kernel_regression_2d_sp(ix, iy, h, silverman, xout, mask, nodata)

    implicit none

    real(sp), dimension(:,:),                       intent(in) :: ix
    real(sp), dimension(:),                         intent(in) :: iy
    real(sp), dimension(:),               optional, intent(in) :: h
    logical,                              optional, intent(in) :: silverman
    real(sp), dimension(:,:),             optional, intent(in) :: xout
    logical,  dimension(:),               optional, intent(in) :: mask
    real(sp),                             optional, intent(in) :: nodata
    real(sp), dimension(:),   allocatable                      :: kernel_regression_2d_sp

    ! local variables
    integer(i4)                           :: dims
    integer(i4)                           :: nin, nout
    integer(i4)                           :: ii, jj
    real(sp), dimension(size(ix,2))       :: hh
    real(sp), dimension(:,:), allocatable :: xxout
    real(sp), dimension(:,:), allocatable :: z
    real(sp), dimension(:,:), allocatable :: x
    real(sp), dimension(:),   allocatable :: y

    ! consistency check - mask needs either nodata or xout
    if (present(mask) .and. (.not. present(xout)) .and. (.not. present(nodata)) ) then
       stop 'kernel_regression_1d_sp: missing nodata value or xout with present(mask)'
    end if

    nin  = size(ix,1)
    dims = size(ix,2)
    ! consistency checks of inputs
    if (size(iy) .ne. nin) stop 'kernel_regression_2d_sp: size(y) /= size(x,1)'
    if (present(h)) then
       if (size(h) .ne. dims) stop 'kernel_regression_2d_sp: size(h) /= size(x,2)'
    end if
    if (present(xout)) then
       if (size(xout,2) .ne. dims) stop 'kernel_regression_2d_sp: size(xout) /= size(x,2)'
    end if

    ! Pack if mask present
    if (present(mask)) then
       nin = count(mask)
       allocate(x(nin,dims))
       allocate(y(nin))
       forall(ii=1:dims) x(:,ii) = pack(ix(:,ii), mask)
       y = pack(iy, mask)
    else
       nin = size(ix,1)
       allocate(x(nin,dims))
       allocate(y(nin))
       x = ix
       y = iy
    endif
    allocate(z(nin,dims))

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_regression_h(x, y, silverman=silverman)
       else
          hh = kernel_regression_h(x, y, silverman=.true.)
       end if
    end if

    ! calc regression
    if (present(xout)) then
       nout = size(xout,1)
       allocate(xxout(nout,dims))
       xxout = xout
    else
       nout = nin
       allocate(xxout(nout,dims))
       xxout = x
    end if
    ! allocate output
    allocate(kernel_regression_2d_sp(nout))

    ! loop through all regression points
    do ii = 1, nout
       forall(jj=1:dims) z(:,jj) = (x(:,jj) - xxout(ii,jj)) / hh(jj)
       ! nadaraya-watson estimator of gaussian multivariate kernel
       kernel_regression_2d_sp(ii) = nadaraya_watson(z, y)
    end do

    ! check whether output has to be unpacked
    if (present(mask) .and. (.not. present(xout))) then
       deallocate(y)
       nin = size(iy,1)
       allocate(y(nin))
       y = unpack(kernel_regression_2d_sp, mask, nodata)
       deallocate(kernel_regression_2d_sp)
       allocate(kernel_regression_2d_sp(nin))
       kernel_regression_2d_sp = y
    end if

    ! clean up
    deallocate(xxout)
    deallocate(x)
    deallocate(y)
    deallocate(z)

  end function kernel_regression_2d_sp

  ! ------------------------------------------------------------------------------------------------

  function kernel_regression_h_1d_dp(ix, iy, silverman, mask)

    implicit none

    real(dp), dimension(:),           intent(in) :: ix
    real(dp), dimension(:),           intent(in) :: iy
    logical,                optional, intent(in) :: silverman
    logical,  dimension(:), optional, intent(in) :: mask
    real(dp)                                     :: kernel_regression_h_1d_dp

    ! local variables
    integer(i4)              :: nin
    real(dp)                 :: nn
    real(dp), dimension(1)   :: h
    real(dp), dimension(1,2) :: bounds
    real(dp), parameter      :: pre_h = 1.05922384104881_dp
    real(dp), dimension(:), allocatable :: x, y

    if (present(mask)) then
       nin = count(mask)
       allocate(x(nin))
       allocate(y(nin))
       x = pack(ix, mask)
       y = pack(iy, mask)
    else
       nin = size(ix,1)
       allocate(x(nin))
       allocate(y(nin))
       x = ix
       y = iy
    endif
    nn = real(nin,dp)

    ! Silverman's rule of thumb by
    ! Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    !h(1) = (4._dp/3._dp/real(nn,dp))**(0.2_dp) * stddev_x
    h(1) = pre_h/(nn**0.2_dp) * stddev(x(:))

    if (present(silverman)) then
       if (.not. silverman) then
          bounds(1,1) = max(0.2_dp * h(1), (maxval(x)-minval(x))/nn)
          bounds(1,2) = 5.0_dp * h(1)
          call allocate_globals(x,y)
          h = nelminrange(cross_valid_regression_dp, h, bounds)
          call deallocate_globals()
       end if
    end if

    kernel_regression_h_1d_dp = h(1)

    deallocate(x)
    deallocate(y)

  end function kernel_regression_h_1d_dp

  function kernel_regression_h_1d_sp(ix, iy, silverman, mask)

    implicit none

    real(sp), dimension(:),           intent(in) :: ix
    real(sp), dimension(:),           intent(in) :: iy
    logical,                optional, intent(in) :: silverman
    logical,  dimension(:), optional, intent(in) :: mask
    real(sp)                                     :: kernel_regression_h_1d_sp

    ! local variables
    integer(i4)              :: nin
    real(sp)                 :: nn
    real(sp), dimension(1)   :: h
    real(sp), dimension(1,2) :: bounds
    real(sp), parameter      :: pre_h = 1.05922384104881_sp
    real(sp), dimension(:), allocatable :: x, y

    if (present(mask)) then
       nin = count(mask)
       allocate(x(nin))
       allocate(y(nin))
       x = pack(ix, mask)
       y = pack(iy, mask)
    else
       nin = size(ix,1)
       allocate(x(nin))
       allocate(y(nin))
       x = ix
       y = iy
    endif
    nn = real(nin,sp)

    ! Silverman's rule of thumb by
    ! Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    !h(1) = (4._sp/3._sp/real(nn,sp))**(0.2_sp) * stddev_x
    h(1) = pre_h/(nn**0.2_sp) * stddev(x(:))

    if (present(silverman)) then
       if (.not. silverman) then
          bounds(1,1) = max(0.2_sp * h(1), (maxval(x)-minval(x))/nn)
          bounds(1,2) = 5.0_sp * h(1)
          call allocate_globals(x,y)
          h = nelminrange(cross_valid_regression_sp, h, bounds)
          call deallocate_globals()
       end if
    end if

    kernel_regression_h_1d_sp = h(1)

    deallocate(x)
    deallocate(y)

  end function kernel_regression_h_1d_sp

  function kernel_regression_h_2d_dp(ix, iy, silverman, mask)

    implicit none

    real(dp), dimension(:,:),                 intent(in) :: ix
    real(dp), dimension(:),                   intent(in) :: iy
    logical,                        optional, intent(in) :: silverman
    logical, dimension(:),          optional, intent(in) :: mask
    real(dp), dimension(size(ix,2))                       :: kernel_regression_h_2d_dp

    ! local variables
    integer(i4)                      :: dims, nn, ii
    real(dp), dimension(size(ix,2))   :: h
    real(dp), dimension(size(ix,2))   :: stddev_x
    real(dp), dimension(size(ix,2),2) :: bounds
    real(dp), dimension(:,:), allocatable :: x
    real(dp), dimension(:),   allocatable :: y

    dims = size(ix,2)
    if (present(mask)) then
       nn = count(mask)
       allocate(x(nn,dims))
       allocate(y(nn))
       forall(ii=1:dims) x(:,ii) = pack(ix(:,ii), mask)
       y = pack(iy, mask)
    else
       nn = size(ix,1)
       allocate(x(nn,dims))
       allocate(y(nn))
       x = ix
       y = iy
    endif
    ! Silverman's rule of thumb by
    ! Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    do ii=1,dims
       stddev_x(ii) = stddev(x(:,ii))
    end do
    h(:) = (4.0_dp/real(dims+2,dp)/real(nn,dp))**(1.0_dp/real(dims+4,dp)) * stddev_x(:)

    if (present(silverman)) then
       if (.not. silverman) then
          bounds(:,1) = max(0.2_dp * h(:), (maxval(x,dim=1)-minval(x,dim=1))/real(nn,dp))
          bounds(:,2) = 5.0_dp * h(:)
          call allocate_globals(x,y)
          h = nelminrange(cross_valid_regression_dp, h, bounds)
          call deallocate_globals()
       end if
    end if

    kernel_regression_h_2d_dp = h

    deallocate(x)
    deallocate(y)

  end function kernel_regression_h_2d_dp

  function kernel_regression_h_2d_sp(ix, iy, silverman, mask)

    implicit none

    real(sp), dimension(:,:),                 intent(in) :: ix
    real(sp), dimension(:),                   intent(in) :: iy
    logical,                        optional, intent(in) :: silverman
    logical, dimension(:),          optional, intent(in) :: mask
    real(sp), dimension(size(ix,2))                       :: kernel_regression_h_2d_sp

    ! local variables
    integer(i4)                      :: dims, nn, ii
    real(sp), dimension(size(ix,2))   :: h
    real(sp), dimension(size(ix,2))   :: stddev_x
    real(sp), dimension(size(ix,2),2) :: bounds
    real(sp), dimension(:,:), allocatable :: x
    real(sp), dimension(:),   allocatable :: y

    dims = size(ix,2)
    if (present(mask)) then
       nn = count(mask)
       allocate(x(nn,dims))
       allocate(y(nn))
       forall(ii=1:dims) x(:,ii) = pack(ix(:,ii), mask)
       y = pack(iy, mask)
    else
       nn = size(ix,1)
       allocate(x(nn,dims))
       allocate(y(nn))
       x = ix
       y = iy
    endif
    ! Silverman's rule of thumb by
    ! Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    do ii=1,dims
       stddev_x(ii) = stddev(x(:,ii))
    end do
    h(:) = (4.0_sp/real(dims+2,sp)/real(nn,sp))**(1.0_sp/real(dims+4,sp)) * stddev_x(:)

    if (present(silverman)) then
       if (.not. silverman) then
          bounds(:,1) = max(0.2_sp * h(:), (maxval(x,dim=1)-minval(x,dim=1))/real(nn,sp))
          bounds(:,2) = 5.0_sp * h(:)
          call allocate_globals(x,y)
          h = nelminrange(cross_valid_regression_sp, h, bounds)
          call deallocate_globals()
       end if
    end if

    kernel_regression_h_2d_sp = h

    deallocate(x)
    deallocate(y)

  end function kernel_regression_h_2d_sp

  ! ------------------------------------------------------------------------------------------------
  !
  !                PRIVATE ROUTINES
  !
  ! ------------------------------------------------------------------------------------------------

  function nadaraya_watson_1d_dp(z, y, mask, valid)

    implicit none

    real(dp), dimension(:),             intent(in)  :: z
    real(dp), dimension(:),   optional, intent(in)  :: y
    logical,  dimension(:),   optional, intent(in)  :: mask
    logical,                  optional, intent(out) :: valid
    real(dp)                                        :: nadaraya_watson_1d_dp

    ! local variables
    real(dp), dimension(size(z,1)) :: w
    real(dp)                       :: sum_w
    logical,  dimension(size(z,1)) :: mask1d
    real(dp)                       :: large_z
    real(dp), dimension(size(z,1)) :: ztmp

    if (present(mask)) then
       mask1d = mask
    else
       mask1d = .true.
    end if

    large_z = sqrt(-2.0_dp*log(tiny(1.0_dp)*sqrt(twopi_dp)))
    where (mask1d) ! temporary store z in w in case of mask
       ztmp = z
    elsewhere
       ztmp = 0.0_dp
    end where
    where (mask1d .and. (abs(ztmp) .lt. large_z))
       w = (1.0_dp/sqrt(twopi_dp)) * exp(-0.5_dp*z*z)
    elsewhere
       w = 0.0_dp
    end where
    sum_w = sum(w, mask=mask1d)

    if (present(valid)) valid = .true.
    if (present(y)) then       ! kernel regression
       if (sum_w .lt. tiny(1.0_dp)) then
          nadaraya_watson_1d_dp = huge(1.0_dp)
          if (present(valid)) valid = .false.
       else
          nadaraya_watson_1d_dp = sum(w*y,mask=mask1d) / sum_w
       end if
    else                        ! kernel density
       nadaraya_watson_1d_dp = sum_w
    end if

  end function nadaraya_watson_1d_dp

  function nadaraya_watson_1d_sp(z, y, mask, valid)

    implicit none

    real(sp), dimension(:),             intent(in)  :: z
    real(sp), dimension(:),   optional, intent(in)  :: y
    logical,  dimension(:),   optional, intent(in)  :: mask
    logical,                  optional, intent(out) :: valid
    real(sp)                                        :: nadaraya_watson_1d_sp

    ! local variables
    real(sp), dimension(size(z,1)) :: w
    real(sp)                       :: sum_w
    logical,  dimension(size(z,1)) :: mask1d
    real(sp)                       :: large_z
    real(sp), dimension(size(z,1)) :: ztmp

    if (present(mask)) then
       mask1d = mask
    else
       mask1d = .true.
    end if

    large_z = sqrt(-2.0_sp*log(tiny(1.0_sp)*sqrt(twopi_sp)))
    where (mask1d) ! temporary store z in w in case of mask
       ztmp = z
    elsewhere
       ztmp = 0.0_sp
    end where
    where (mask1d .and. (abs(ztmp) .lt. large_z))
       w = (1.0_sp/sqrt(twopi_sp)) * exp(-0.5_sp*z*z)
    elsewhere
       w = 0.0_sp
    end where
    sum_w = sum(w, mask=mask1d)

    if (present(valid)) valid = .true.
    if (present(y)) then       ! kernel regression
       if (sum_w .lt. tiny(1.0_sp)) then
          nadaraya_watson_1d_sp = huge(1.0_sp)
          if (present(valid)) valid = .false.
       else
          nadaraya_watson_1d_sp = sum(w*y,mask=mask1d) / sum_w
       end if
    else                        ! kernel density
       nadaraya_watson_1d_sp = sum_w
    end if

  end function nadaraya_watson_1d_sp

  function nadaraya_watson_2d_dp(z, y, mask, valid)

    implicit none

    real(dp), dimension(:,:),           intent(in)  :: z
    real(dp), dimension(:),   optional, intent(in)  :: y
    logical,  dimension(:),   optional, intent(in)  :: mask
    logical,                  optional, intent(out) :: valid
    real(dp)                                        :: nadaraya_watson_2d_dp

    ! local variables
    real(dp), dimension(size(z,1), size(z,2)) :: kerf
    real(dp), dimension(size(z,1))            :: w
    real(dp)                                  :: sum_w
    logical,  dimension(size(z,1))            :: mask1d
    logical,  dimension(size(z,1), size(z,2)) :: mask2d
    real(dp)                                  :: large_z
    real(dp), dimension(size(z,1), size(z,2)) :: ztmp

    if (present(mask)) then
       mask1d = mask
       mask2d = spread(mask1d, dim=2, ncopies=size(z,2))
    else
       mask1d = .true.
       mask2d = .true.
    end if

    large_z = sqrt(-2.0_dp*log(tiny(1.0_dp)*sqrt(twopi_dp)))
    where (mask2d) ! temporary store z in kerf in case of mask
       ztmp = z
    elsewhere
       ztmp = 0.0_dp
    end where
    where (mask2d .and. (abs(ztmp) .lt. large_z))
       kerf = (1.0_dp/sqrt(twopi_dp)) * exp(-0.5_dp*z*z)
    elsewhere
       kerf = 0.0_dp
    end where
    w     = product(kerf, dim=2, mask=mask2d)
    sum_w = sum(w, mask=mask1d)

    if (present(valid)) valid = .true.
    if (present(y)) then       ! kernel regression
       if (sum_w .lt. tiny(1.0_dp)) then
          nadaraya_watson_2d_dp = huge(1.0_dp)
          if (present(valid)) valid = .false.
       else
          nadaraya_watson_2d_dp = sum(w*y,mask=mask1d) / sum_w
       end if
    else                       ! kernel density
       nadaraya_watson_2d_dp = sum_w
    end if

  end function nadaraya_watson_2d_dp

  function nadaraya_watson_2d_sp(z, y, mask, valid)

    implicit none

    real(sp), dimension(:,:),           intent(in)  :: z
    real(sp), dimension(:),   optional, intent(in)  :: y
    logical,  dimension(:),   optional, intent(in)  :: mask
    logical,                  optional, intent(out) :: valid
    real(sp)                                        :: nadaraya_watson_2d_sp

    ! local variables
    real(sp), dimension(size(z,1), size(z,2)) :: kerf
    real(sp), dimension(size(z,1))            :: w
    real(sp)                                  :: sum_w
    logical,  dimension(size(z,1))            :: mask1d
    logical,  dimension(size(z,1), size(z,2)) :: mask2d
    real(sp)                                  :: large_z
    real(sp), dimension(size(z,1), size(z,2)) :: ztmp

    if (present(mask)) then
       mask1d = mask
       mask2d = spread(mask1d, dim=2, ncopies=size(z,2))
    else
       mask1d = .true.
       mask2d = .true.
    end if

    large_z = sqrt(-2.0_sp*log(tiny(1.0_sp)*sqrt(twopi_sp)))
    where (mask2d) ! temporary store z in kerf in case of mask
       ztmp = z
    elsewhere
       ztmp = 0.0_sp
    end where
    where (mask2d .and. (abs(ztmp) .lt. large_z))
       kerf = (1.0_sp/sqrt(twopi_sp)) * exp(-0.5_sp*z*z)
    elsewhere
       kerf = 0.0_sp
    end where
    w     = product(kerf, dim=2, mask=mask2d)
    sum_w = sum(w, mask=mask1d)

    if (present(valid)) valid = .true.
    if (present(y)) then       ! kernel regression
       if (sum_w .lt. tiny(1.0_sp)) then
          nadaraya_watson_2d_sp = huge(1.0_sp)
          if (present(valid)) valid = .false.
       else
          nadaraya_watson_2d_sp = sum(w*y,mask=mask1d) / sum_w
       end if
    else                       ! kernel density
       nadaraya_watson_2d_sp = sum_w
    end if

  end function nadaraya_watson_2d_sp

  ! ------------------------------------------------------------------------------------------------

  function cross_valid_regression_dp(h)

    implicit none

    ! Helper function that calculates cross-validation function for the
    ! Nadaraya-Watson estimator, which is basically the mean square error
    ! between the observations and the model estimates without the specific point,
    ! i.e. the jackknife estimate (Haerdle et al. 2000).

    ! Function is always 2D because allocate_globals makes always 2D arrays.

    real(dp), dimension(:), intent(in) :: h
    real(dp)                           :: cross_valid_regression_dp

    ! local variables
    integer(i4)                                                    :: ii, kk, nn
    logical,  dimension(size(global_x_dp,1))                       :: mask
    real(dp), dimension(size(global_x_dp,1))                       :: out
    real(dp), dimension(size(global_x_dp,1),size(global_x_dp,2))   :: zz
    logical                                                        :: valid, valid_tmp

    nn = size(global_x_dp,1)
    ! Loop through each regression point and calc kernel estimate without that point (Jackknife)
    valid = .true.
    mask  = .true.
    do ii=1, nn
       mask(ii) = .false.
       zz(:,:) = 0.0_dp
       forall(kk=1:nn, mask(kk)) zz(kk,:) = (global_x_dp(kk,:) - global_x_dp(ii,:)) / h(:)
       out(ii) = nadaraya_watson(zz, y=global_y_dp, mask=mask, valid=valid_tmp)
       valid = valid .and. valid_tmp
       mask(ii) = .true.
    end do

    ! Mean square deviation
    if (valid) then
       cross_valid_regression_dp = sum((global_y_dp-out)**2) / real(nn,dp)
    else
       cross_valid_regression_dp = huge(1.0_dp)
    end if

    print*, 'cross_valid_regression_dp = ', cross_valid_regression_dp, '  ( h = ', h, ' )'

  end function cross_valid_regression_dp

  function cross_valid_regression_sp(h)

    implicit none

    ! Helper function that calculates cross-validation function for the
    ! Nadaraya-Watson estimator, which is basically the mean square error
    ! between the observations and the model estimates without the specific point,
    ! i.e. the jackknife estimate (Haerdle et al. 2000).

    ! Function is always 2D because set_global_for_opti makes always 2D arrays.

    real(sp), dimension(:), intent(in) :: h
    real(sp)                           :: cross_valid_regression_sp

    ! local variables
    integer(i4)                                                    :: ii, kk, nn
    logical,  dimension(size(global_x_sp,1))                       :: mask
    real(sp), dimension(size(global_x_sp,1))                       :: out
    real(sp), dimension(size(global_x_sp,1),size(global_x_sp,2))   :: zz
    logical                                                        :: valid, valid_tmp

    nn = size(global_x_sp,1)
    ! Loop through each regression point and calc kernel estimate without that point (Jackknife)
    valid = .true.
    mask  = .true.
    do ii=1, nn
       mask(ii) = .false.
       forall(kk=1:nn, mask(kk)) zz(kk,:) = (global_x_sp(kk,:) - global_x_sp(ii,:)) / h(:)
       out(ii) = nadaraya_watson(zz, y=global_y_sp, mask=mask, valid=valid_tmp)
       valid = valid .and. valid_tmp
       mask(ii) = .true.
    end do

    ! Mean square deviation
    if (valid) then
       cross_valid_regression_sp = sum((global_y_sp-out)**2) / real(nn,sp)
    else
       cross_valid_regression_sp = huge(1.0_sp)
    end if

  end function cross_valid_regression_sp

  ! ------------------------------------------------------------------------------------------------

  function cross_valid_density_1d_dp(h)

    implicit none

    ! Helper function that calculates cross-validation function for the
    ! Nadaraya-Watson estimator, which is basically the mean square error
    ! where model estimate is replaced by the jackknife estimate (Haerdle et al. 2000).

    real(dp), intent(in)               :: h
    real(dp)                           :: cross_valid_density_1d_dp

    ! local variables
    integer(i4)                              :: ii, nn
    logical,  dimension(size(global_x_dp,1)) :: mask
    real(dp), dimension(size(global_x_dp,1)) :: out
    real(dp), dimension(size(global_x_dp,1)) :: zz
    real(dp)                                 :: delta
    integer(i4)                              :: mesh_n
    real(dp), dimension(:), allocatable      :: xMeshed
    real(dp), dimension(:), allocatable      :: outIntegral
    real(dp), dimension(size(global_x_dp,1)) :: zzIntegral
    real(dp)                                 :: stddev_x
    real(dp)                                 :: summ, multiplier, thresh

    nn   = size(global_x_dp,1)
    if (nn .le. 100_i4) then       ! if few number of data points given, mesh consists of 100*n points
       mesh_n = 100_i4*nn
    else                           ! mesh_n such that mesh consists of not more than 10000 points
       mesh_n = Max(2_i4, 10000_i4/nn) * nn
    end if
    allocate(xMeshed(mesh_n))
    allocate(outIntegral(mesh_n))

    ! integral of squared density function, i.e. L2-norm of kernel^2
    stddev_x = stddev(global_x_dp(:,1))
    xMeshed  = mesh(minval(global_x_dp) - 3.0_dp*stddev_x, maxval(global_x_dp) + 3.0_dp*stddev_x, mesh_n, delta)

    multiplier = 1.0_dp/(real(nn,dp)*h)
    if (multiplier <= 1.0_dp) then
       thresh = tiny(1.0_dp)/multiplier
    else
       thresh = 0.0_dp
    endif
    !$OMP parallel default(shared) &
    !$OMP private(zzIntegral)
    !$OMP do
    do ii=1, mesh_n
       zzIntegral = (global_x_dp(:,1) - xMeshed(ii)) / h 
       outIntegral(ii) = nadaraya_watson(zzIntegral)
       if (outIntegral(ii) .gt. thresh) outIntegral(ii) = multiplier * outIntegral(ii)
    end do
    !$OMP end do
    !$OMP end parallel
    where (outIntegral .lt. sqrt(tiny(1.0_dp)) )
       outIntegral = 0.0_dp
    end where
    summ = int_regular(outIntegral*outIntegral, delta)

    ! leave-one-out estimate
    mask = .true.
    !$OMP parallel default(shared) &
    !$OMP private(mask, zz)
    !$OMP do
    do ii=1, nn
       mask(ii) = .false.
       zz       = (global_x_dp(:,1) - global_x_dp(ii,1)) / h 
       out(ii)  = nadaraya_watson(zz, mask=mask)
       if (out(ii) .gt. thresh) out(ii) = multiplier * out(ii)
       mask(ii) = .true.
    end do
    !$OMP end do
    !$OMP end parallel

    cross_valid_density_1d_dp = summ - 2.0_dp / (real(nn,dp)) * sum(out)
    write(*,*) 'h = ',h, '   cross_valid = ',cross_valid_density_1d_dp

    ! clean up
    deallocate(xMeshed)
    deallocate(outIntegral)

  end function cross_valid_density_1d_dp

  function cross_valid_density_1d_sp(h)

    implicit none

    ! Helper function that calculates cross-validation function for the
    ! Nadaraya-Watson estimator, which is basically the mean square error
    ! where model estimate is replaced by the jackknife estimate (Haerdle et al. 2000).

    real(sp),               intent(in) :: h
    real(sp)                           :: cross_valid_density_1d_sp

    ! local variables
    integer(i4)                              :: ii, nn
    logical,  dimension(size(global_x_sp,1)) :: mask
    real(sp), dimension(size(global_x_sp,1)) :: out
    real(sp), dimension(size(global_x_sp,1)) :: zz
    real(sp)                                 :: delta
    integer(i4)                              :: mesh_n
    real(sp), dimension(:), allocatable      :: xMeshed
    real(sp), dimension(:), allocatable      :: outIntegral
    real(sp), dimension(size(global_x_sp,1)) :: zzIntegral
    real(sp)                                 :: stddev_x
    real(sp)                                 :: summ, multiplier, thresh

    nn   = size(global_x_sp,1)
    if (nn .le. 100_i4) then       ! if few number of data points given, mesh consists of 100*n points
       mesh_n = 100_i4*nn
    else                           ! mesh_n such that mesh consists of not more than 10000 points
       mesh_n = Max(2_i4, 10000_i4/nn) * nn
    end if
    allocate(xMeshed(mesh_n))
    allocate(outIntegral(mesh_n))

    ! integral of squared density function, i.e. L2-norm of kernel^2
    stddev_x = stddev(global_x_sp(:,1))
    xMeshed  = mesh(minval(global_x_sp) - 3.0_sp*stddev_x, maxval(global_x_sp) + 3.0_sp*stddev_x, mesh_n, delta)

    multiplier = 1.0_sp/(real(nn,sp)*h)
    if (multiplier <= 1.0_sp) then
       thresh = tiny(1.0_sp)/multiplier
    else
       thresh = 0.0_sp
    endif
    !$OMP parallel default(shared) &
    !$OMP private(zzIntegral)
    !$OMP do
    do ii=1, mesh_n
       zzIntegral = (global_x_sp(:,1) - xMeshed(ii)) / h
       outIntegral(ii) = nadaraya_watson(zzIntegral)
       if (outIntegral(ii) .gt. thresh) outIntegral(ii) = multiplier * outIntegral(ii)
    end do
    !$OMP end do
    !$OMP end parallel
    where (outIntegral .lt. sqrt(tiny(1.0_sp)) )
       outIntegral = 0.0_sp
    end where
    summ = int_regular(outIntegral*outIntegral, delta)

    ! leave-one-out estimate
    mask = .true.
    !$OMP parallel default(shared) &
    !$OMP private(mask, zz)
    !$OMP do
    do ii=1, nn
       mask(ii) = .false.
       zz       = (global_x_sp(:,1) - global_x_sp(ii,1)) / h
       out(ii)  = nadaraya_watson(zz, mask=mask)
       if (out(ii) .gt. thresh) out(ii) = multiplier * out(ii)
       mask(ii) = .true.
    end do
    !$OMP end do
    !$OMP end parallel

    cross_valid_density_1d_sp = summ - 2.0_sp / (real(nn,sp)) * sum(out)

    ! clean up
    deallocate(xMeshed)
    deallocate(outIntegral)

  end function cross_valid_density_1d_sp

  ! ------------------------------------------------------------------------------------------------

  subroutine allocate_globals_1d_dp(x,y,xout)

    implicit none

    real(dp), dimension(:),           intent(in) :: x
    real(dp), dimension(:), optional, intent(in) :: y
    real(dp), dimension(:), optional, intent(in) :: xout

    allocate( global_x_dp(size(x,1),1) )
    global_x_dp(:,1) = x

    if (present(y)) then
       allocate( global_y_dp(size(y,1)) )
       global_y_dp = y
    end if

    if (present(xout)) then
       allocate( global_xout_dp(size(xout,1),1) )
       global_xout_dp(:,1) = xout
    end if

  end subroutine allocate_globals_1d_dp

  subroutine allocate_globals_1d_sp(x,y,xout)

    implicit none

    real(sp), dimension(:),           intent(in) :: x
    real(sp), dimension(:), optional, intent(in) :: y
    real(sp), dimension(:), optional, intent(in) :: xout

    allocate( global_x_sp(size(x,1),1) )
    global_x_sp(:,1) = x

    if (present(y)) then
       allocate( global_y_sp(size(y,1)) )
       global_y_sp = y
    end if

    if (present(xout)) then
       allocate( global_xout_sp(size(xout,1),1) )
       global_xout_sp(:,1) = xout
    end if

  end subroutine allocate_globals_1d_sp

  subroutine allocate_globals_2d_dp(x,y,xout)

    implicit none

    real(dp), dimension(:,:),           intent(in) :: x
    real(dp), dimension(:),   optional, intent(in) :: y
    real(dp), dimension(:,:), optional, intent(in) :: xout

    allocate( global_x_dp(size(x,1),size(x,2)) )
    global_x_dp = x

    if (present(y)) then
       allocate( global_y_dp(size(y,1)) )
       global_y_dp = y
    end if

    if (present(xout)) then
       allocate( global_xout_dp(size(xout,1),size(xout,2)) )
       global_xout_dp = xout
    end if

  end subroutine allocate_globals_2d_dp

  subroutine allocate_globals_2d_sp(x,y,xout)

    implicit none

    real(sp), dimension(:,:),           intent(in) :: x
    real(sp), dimension(:),   optional, intent(in) :: y
    real(sp), dimension(:,:), optional, intent(in) :: xout

    allocate( global_x_sp(size(x,1),size(x,2)) )
    global_x_sp = x

    if (present(y)) then
       allocate( global_y_sp(size(y,1)) )
       global_y_sp = y
    end if

    if (present(xout)) then
       allocate( global_xout_sp(size(xout,1),size(xout,2)) )
       global_xout_sp = xout
    end if

  end subroutine allocate_globals_2d_sp

  ! ------------------------------------------------------------------------------------------------

  subroutine deallocate_globals()

    implicit none

    if (allocated(global_x_dp))    deallocate(global_x_dp)
    if (allocated(global_y_dp))    deallocate(global_y_dp)
    if (allocated(global_xout_dp)) deallocate(global_xout_dp)
    if (allocated(global_x_sp))    deallocate(global_x_sp)
    if (allocated(global_y_sp))    deallocate(global_y_sp)
    if (allocated(global_xout_sp)) deallocate(global_xout_sp)

  end subroutine deallocate_globals

  ! ------------------------------------------------------------------------------------------------

  function mesh_dp(start, end, n, delta)

    implicit none

    real(dp),    intent(in)  :: start
    real(dp),    intent(in)  :: end
    integer(i4), intent(in)  :: n
    real(dp),    intent(out) :: delta
    real(dp), dimension(n)   :: mesh_dp

    ! local variables
    integer(i4) :: ii

    delta = (end-start) / real(n-1,dp)
    forall(ii=1:n) mesh_dp(ii) = start + (ii-1) * delta

  end function mesh_dp

  function mesh_sp(start, end, n, delta)

    implicit none

    real(sp),    intent(in)  :: start
    real(sp),    intent(in)  :: end
    integer(i4), intent(in)  :: n
    real(sp),    intent(out) :: delta
    real(sp), dimension(n)   :: mesh_sp

    ! local variables
    integer(i4) :: ii

    delta = (end-start) / real(n-1,sp)
    forall(ii=1:n) mesh_sp(ii) = start + (ii-1) * delta

  end function mesh_sp

END MODULE mo_kernel

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
!>          Therefore, the bandwith estimation is not necessarily a part of the kernel smoothing
!>          but can be determined first and given as an optional argument to the smoother.

!> \authors Juliane Mai
!> \date Mar 2013

MODULE mo_kernel

  ! This module provides functions for kernel regression and kernel density estimation and
  ! is part of the UFZ CHS Fortran library.

  ! Written  Juliane Mai, Mar 2013
  ! Modified

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

  ! Copyright 2013 Juliane Mai

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

  USE mo_kind,      ONLY: i4, sp, dp
  USE mo_constants, ONLY: pi_sp, pi_dp
  USE mo_moment,    ONLY: stddev
  USE mo_nelmin,    ONLY: nelminrange
  USE mo_sort,      ONLY: sort_index

  IMPLICIT NONE

  PUBLIC :: kernel_cumdensity        ! Kernel smoothing of a CDF                (only 1D)
  PUBLIC :: kernel_density           ! Kernel smoothing of a PDF                (only 1D)
  PUBLIC :: kernel_density_h         ! Bandwith estimation for PDF and CDF      (only 1D)
  PUBLIC :: kernel_regression        ! Kernel regression                        (1D and ND)
  PUBLIC :: kernel_regression_h      ! Bandwith estimatio for kernel regression (1D and ND)

  ! ------------------------------------------------------------------

  !     NAME
  !         kernel_cumdensity

  !     PURPOSE
  !         Approximates the cumulative density function (CDF) to a given 1D data set using a Gaussian kernel.
  !
  !>        \brief   Approximates the cumulative density function (CDF).
  !
  !>        \details Approximates the cummulative density function (CDF)
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
  !>        For large data sets this might be time consuming and should be performed aforehand using the function kernel_density_h.\n
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
  !>       \param[in] "logical, optional :: silverman"   By default Silverman''s Rule-of-thumb will be used to approximate the
  !>                                                     bandwith of the kernel (silverman=true).
  !>                                                     If silverman=false the Cross-Validation approach is used
  !>                                                     to estimate the bandwidth.
  !>       \param[in] "real(sp/dp), optional :: xout(:)" If present, the CDF will be approximated at this arguments,
  !>                                                     otherwise the CDF is approximated at x.
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
  !         ! pre-estimate bandwidth via cross-validation (time-consuming)
  !         h = kernel_density_h(x,silverman=.false.)
  !
  !         ! pre-estimate bandwidth with Silverman''s rule of thumb (default)
  !         h = kernel_density_h(x,silverman=.true.)
  !
  !         cdf = kernel_cumdensity(x, h=h, silverman=.false., xout=xout)
  !         ! gives cumulative density at either xout values, if specified, or at x values, if xout is not present
  !         ! if bandwith h is given                 : silverman (true/false) is ignored
  !         ! if silverman=.true.  and h not present : bandwith will be estimated using Silerman''s rule of thumb
  !         ! if silverman=.false. and h not present : bandwith will be estimated using Cross-Validation approach
  !         -> see also example in test directory

  !     LITERATURE
  !         Scott, D. W., & Sain, S. R. (2005).
  !             Multi-dimensional Density Estimation. Handbook of Statistics, 24, 229–261.
  !             doi:10.1016/S0169-7161(04)24009-3
  !         Haerdle, W., & Mueller, M. (2000). Multivariate and semiparametric kernel regression.
  !             In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
  !             application (pp. 357-392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12s,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Mar 2013
  !         Modified,

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
  !>        For large data sets this might be time consuming and should be performed aforehand using the function kernel_density_h.\n
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
  !>       \param[in] "logical, optional :: silverman"   By default Silverman''s Rule-of-thumb will be used to approximate the
  !>                                                     bandwith of the kernel (silverman=true).
  !>                                                     If silverman=false the Cross-Validation approach is used
  !>                                                     to estimate the bandwidth.
  !>       \param[in] "real(sp/dp), optional :: xout(:)" If present, the PDF will be approximated at this arguments,
  !>                                                     otherwise the PDF is approximated at x.
  !>       \param[in] "logical,     optional :: mask(:)" If present, only consider these data points for kernel density
  !>       \param[in] "real(sp/dp), optional :: nodata"  If mask is present and xout not, then output is unpacked
  !>                                                     using this value at non masked positions
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
  !         ! pre-estimate bandwidth via cross-validation (time-consuming)
  !         h = kernel_density_h(x,silverman=.false.)
  !
  !         ! pre-estimate bandwidth with Silverman''s rule of thumb (default)
  !         h = kernel_density_h(x,silverman=.true.)
  !
  !         pdf = kernel_density(x, h=h, silverman=.false., xout=xout)
  !         ! gives (probability) density at either xout values, if specified, or at x values, if xout is not present
  !         ! if bandwith h is given                 : silverman (true/false) is ignored
  !         ! if silverman=.true.  and h not present : bandwith will be estimated using Silerman''s rule of thumb
  !         ! if silverman=.false. and h not present : bandwith will be estimated using Cross-Validation approach
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         Scott, D. W., & Sain, S. R. (2005).
  !             Multi-dimensional Density Estimation. Handbook of Statistics, 24, 229–261.
  !             doi:10.1016/S0169-7161(04)24009-3
  !         Haerdle, W., & Mueller, M. (2000). Multivariate and semiparametric kernel regression.
  !             In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
  !             application (pp. 357-392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12s,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Mar 2013
  !         Modified, Stephan Thober Mar 2013 - added mask and nodata value

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
  !>        For large data sets this might be time consuming and should be performed aforehand using the function kernel_density_h.\n
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
  !>       \param[in] "logical, optional :: silverman"   By default Silverman''s Rule-of-thumb will be used to approximate the
  !>                                                     bandwith of the kernel (silverman=true).
  !>                                                     If silverman=false the Cross-Validation approach is used
  !>                                                     to estimate the bandwidth.
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
  !         ! pre-estimate bandwidth via cross-validation (time-consuming)
  !         h = kernel_density_h(x,silverman=.false.)
  !
  !         ! pre-estimate bandwidth with Silverman''s rule of thumb (default)
  !         h = kernel_density_h(x,silverman=.true.)
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         Scott, D. W., & Sain, S. R. (2005).
  !             Multi-dimensional Density Estimation. Handbook of Statistics, 24, 229–261.
  !             doi:10.1016/S0169-7161(04)24009-3
  !         Haerdle, W., & Mueller, M. (2000). Multivariate and semiparametric kernel regression.
  !             In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
  !             application (pp. 357-392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12s,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Mar 2013
  !         Modified,

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
  !         ! pre-estimate bandwidth via cross-validation (time-consuming)
  !         h = kernel_regression(x,y,silverman=.false.)
  !
  !         ! pre-estimate bandwidth with Silverman''s rule of thumb (default)
  !         h = kernel_regression(x,y,silverman=.true.)
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
  !             Multi-dimensional Density Estimation. Handbook of Statistics, 24, 229–261.
  !             doi:10.1016/S0169-7161(04)24009-3
  !         Haerdle, W., & Mueller, M. (2000). Multivariate and semiparametric kernel regression.
  !             In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
  !             application (pp. 357-392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12s,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Matthias Cuntz, Juliane Mai
  !>        \date Mar 2013
  !         Modified,

  INTERFACE kernel_regression
     MODULE PROCEDURE kernel_regression_2d_dp, kernel_regression_2d_sp, kernel_regression_1d_dp, kernel_regression_1d_sp
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
  !         ! pre-estimate bandwidth via cross-validation (time-consuming)
  !         h = kernel_regression_h(x,y,silverman=.false.)
  !
  !         ! pre-estimate bandwidth with Silverman''s rule of thumb (default)
  !         h = kernel_regression_h(x,y,silverman=.true.)
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         Scott, D. W., & Sain, S. R. (2005).
  !             Multi-dimensional Density Estimation. Handbook of Statistics, 24, 229–261.
  !             doi:10.1016/S0169-7161(04)24009-3
  !         Haerdle, W., & Mueller, M. (2000). Multivariate and semiparametric kernel regression.
  !             In M. G. Schimek (Ed.), Smoothing and regression: Approaches, computation, and
  !             application (pp. 357-392). Hoboken, NJ, USA: John Wiley & Sons, Inc. doi:10.1002/9781118150658.ch12s,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Matthias Cuntz, Juliane Mai
  !>        \date Mar 2013
  !         Modified,

  INTERFACE kernel_regression_h
     MODULE PROCEDURE kernel_regression_h_2d_dp, kernel_regression_h_2d_sp, kernel_regression_h_1d_dp, kernel_regression_h_1d_sp
  END INTERFACE kernel_regression_h

  ! ------------------------------------------------------------------

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

  function kernel_cumdensity_1d_dp(x,h,silverman,xout)

    implicit none

    real(dp), dimension(:),                       intent(in) :: x
    real(dp),                           optional, intent(in) :: h
    logical,                            optional, intent(in) :: silverman
    real(dp), dimension(:),             optional, intent(in) :: xout
    real(dp), dimension(:), allocatable                      :: kernel_cumdensity_1d_dp

    ! local variables
    integer(i4)                                   :: nn, nout, nmesh
    integer(i4)                                   :: ii, jj
    real(dp)                                      :: hh
    real(dp)                                      :: lower_x
    real(dp),    dimension(:),        allocatable :: xxout
    integer(i4), dimension(:),        allocatable :: xindx
    real(dp),    dimension(:),        allocatable :: kernel_density_1d_dp
    real(dp),    dimension(:),        allocatable :: xmesh
    real(dp),    dimension(size(x,1))             :: z
    real(dp)                                      :: multiplier
    real(dp)                                      :: delta

    nn   = size(x,1)
    nmesh = 100_i4

    ! allocate
    if (present(xout)) then
       allocate(xxout(size(xout,1)))
       allocate(xindx(size(xout,1)))
       xxout = xout
    else
       allocate(xxout(size(x,1)))
       allocate(xindx(size(x,1)))
       xxout = x
    end if
    nout   = size(xxout,1)

    ! sort the x
    xindx = sort_index(xxout)
    xxout = xxout(xindx)

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_density_h_1d_dp(x,silverman=silverman)
       else
          hh = kernel_density_h_1d_dp(x,silverman=.true.)
       end if
    end if

    ! allocate PDF array
    allocate(kernel_density_1d_dp(nmesh))
    ! allocate mesh array
    allocate(xmesh(nmesh))
    ! allocate CDF array
    allocate(kernel_cumdensity_1d_dp(nout))

    ! calculate standard devaition of x to determine left side starting point for integration of PDF
    lower_x = minval(x) - 3.0_dp * stddev(x)

    multiplier = 1.0_dp/(real(nn,dp)*hh)
    ! loop through all regression points
    do ii = 1, nout

       ! generate mesh
       if (ii .gt. 1_i4) then
          xmesh = mesh_dp(xxout(ii-1),xxout(ii),nmesh,delta)
       else
          xmesh = mesh_dp(lower_x,xxout(1),nmesh,delta)
       end if

       do jj = 1, nmesh
          ! scaled difference from regression point
          z(:) = (x(:) - xmesh(jj)) / hh
          ! nadaraya-watson estimator of gaussian multivariate kernel
          kernel_density_1d_dp(jj) = nadaraya_watson_1d_dp(z)
          if (kernel_density_1d_dp(jj) .gt. real(nn,dp)*hh*epsilon(1.0_dp)) then
             kernel_density_1d_dp(jj) = multiplier * nadaraya_watson_1d_dp(z)
          end if
       end do
       kernel_cumdensity_1d_dp(ii) = sum(kernel_density_1d_dp) * delta

       ! generate mesh
       if (ii .gt. 1_i4) then
          kernel_cumdensity_1d_dp(ii) = kernel_cumdensity_1d_dp(ii-1) + kernel_cumdensity_1d_dp(ii)
       end if

       ! print*, xxout(ii),'   ',xindx(ii),'   kernel_cumdensity_1d_dp(',ii,') = ',kernel_cumdensity_1d_dp(ii)

    end do

    ! scale to range [0,1]
    forall(ii=1:nout) kernel_cumdensity_1d_dp(ii) = ( kernel_cumdensity_1d_dp(ii) - kernel_cumdensity_1d_dp(1) ) /  &
         ( kernel_cumdensity_1d_dp(nout) - kernel_cumdensity_1d_dp(1) )

    ! resorting
    forall(ii=1:nout) kernel_cumdensity_1d_dp(xindx(ii)) = kernel_cumdensity_1d_dp(ii)

  end function kernel_cumdensity_1d_dp

  function kernel_cumdensity_1d_sp(x,h,silverman,xout)

    implicit none

    real(sp), dimension(:),                       intent(in) :: x
    real(sp),                           optional, intent(in) :: h
    logical,                            optional, intent(in) :: silverman
    real(sp), dimension(:),             optional, intent(in) :: xout
    real(sp), dimension(:), allocatable                      :: kernel_cumdensity_1d_sp

    ! local variables
    integer(i4)                                   :: nn, nout, nmesh
    integer(i4)                                   :: ii, jj
    real(sp)                                      :: hh
    real(sp)                                      :: lower_x
    real(sp),    dimension(:),        allocatable :: xxout
    integer(i4), dimension(:),        allocatable :: xindx
    real(sp),    dimension(:),        allocatable :: kernel_density_1d_sp
    real(sp),    dimension(:),        allocatable :: xmesh
    real(sp),    dimension(size(x,1))             :: z
    real(sp)                                      :: multiplier
    real(sp)                                      :: delta

    nn   = size(x,1)
    nmesh = 100_i4

    ! allocate
    if (present(xout)) then
       allocate(xxout(size(xout,1)))
       allocate(xindx(size(xout,1)))
       xxout = xout
    else
       allocate(xxout(size(x,1)))
       allocate(xindx(size(x,1)))
       xxout = x
    end if
    nout   = size(xxout,1)

    ! sort the x
    xindx = sort_index(xxout)
    xxout = xxout(xindx)

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_density_h_1d_sp(x,silverman=silverman)
       else
          hh = kernel_density_h_1d_sp(x,silverman=.true.)
       end if
    end if

    ! allocate PDF array
    allocate(kernel_density_1d_sp(nmesh))
    ! allocate mesh array
    allocate(xmesh(nmesh))
    ! allocate CDF array
    allocate(kernel_cumdensity_1d_sp(nout))

    ! calculate standard devaition of x to determine left side starting point for integration of PDF
    lower_x = minval(x) - 3.0_sp * stddev(x)

    multiplier = 1.0_sp/(real(nn,sp)*hh)
    ! loop through all regression points
    do ii = 1, nout

       ! generate mesh
       if (ii .gt. 1_i4) then
          xmesh = mesh_sp(xxout(ii-1),xxout(ii),nmesh,delta)
       else
          xmesh = mesh_sp(lower_x,xxout(1),nmesh,delta)
       end if

       do jj = 1, nmesh
          ! scaled difference from regression point
          z(:) = (x(:) - xmesh(jj)) / hh
          ! nadaraya-watson estimator of gaussian multivariate kernel
          kernel_density_1d_sp(jj) = nadaraya_watson_1d_sp(z)
          if (kernel_density_1d_sp(jj) .gt. real(nn,sp)*hh*epsilon(1.0_sp)) then
             kernel_density_1d_sp(jj) = multiplier * nadaraya_watson_1d_sp(z)
          end if
       end do
       kernel_cumdensity_1d_sp(ii) = sum(kernel_density_1d_sp) * delta

       ! generate mesh
       if (ii .gt. 1_i4) then
          kernel_cumdensity_1d_sp(ii) = kernel_cumdensity_1d_sp(ii-1) + kernel_cumdensity_1d_sp(ii)
       end if

    end do

    ! scale to range [0,1]
    forall(ii=1:nout) kernel_cumdensity_1d_sp(ii) = ( kernel_cumdensity_1d_sp(ii) - kernel_cumdensity_1d_sp(1) ) /  &
         ( kernel_cumdensity_1d_sp(nout) - kernel_cumdensity_1d_sp(1) )

    ! resorting
    forall(ii=1:nout) kernel_cumdensity_1d_sp(xindx(ii)) = kernel_cumdensity_1d_sp(ii)

  end function kernel_cumdensity_1d_sp

  ! ------------------------------------------------------------------------------------------------

  function kernel_density_1d_dp(ix,h,silverman,xout,mask,nodata)

    implicit none

    real(dp), dimension(:),                       intent(in) :: ix
    real(dp),                           optional, intent(in) :: h
    logical,                            optional, intent(in) :: silverman
    real(dp), dimension(:),             optional, intent(in) :: xout
    logical,  dimension(:),             optional, intent(in) :: mask
    real(dp),                           optional, intent(in) :: nodata
    real(dp), dimension(:), allocatable                      :: kernel_density_1d_dp

    ! local variables
    logical,  dimension(size(ix,1))            :: mm
    integer(i4)                                :: nn, nout
    integer(i4)                                :: ii
    real(dp)                                   :: hh
    real(dp), dimension(:),        allocatable :: xxout
    real(dp), dimension(:),        allocatable :: z
    real(dp), dimension(:),        allocatable :: x
    real(dp)                                   :: multiplier

    ! consitency check
    if ( .not. present(xout) .and. present(mask) .and. .not. present(nodata) ) then
       print *, ' ERROR *** missing nodata value or xout. function kernel_density'
       stop ' ERROR*** see StdOut for details'
    end if

    ! initialize x
    mm = .true.
    if ( present( mask ) ) mm = mask

    allocate( z( count(mm) ) )
    allocate( x( count(mm) ) )
    x = pack( ix, mm)

    nn   = size(x,1)

    ! calc regression
    if (present(xout)) then
       allocate(xxout(size(xout,1)))
       xxout = xout
    else
       allocate(xxout(size(x,1)))
       xxout = x
    end if
    nout   = size(xxout,1)

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_density_h_1d_dp(x,silverman=silverman)
       else
          hh = kernel_density_h_1d_dp(x,silverman=.true.)
       end if
    end if

    ! allocate output
    allocate(kernel_density_1d_dp(nout))

    multiplier = 1.0_dp/(real(nn,dp)*hh)
    ! loop through all regression points
    do ii = 1, nout
       ! scaled difference from regression point
       z(:) = (x(:) - xxout(ii)) / hh
       ! nadaraya-watson estimator of gaussian multivariate kernel
       kernel_density_1d_dp(ii) = nadaraya_watson_1d_dp(z)
       if (kernel_density_1d_dp(ii) .gt. real(nn,dp)*hh*epsilon(1.0_dp)) then
          kernel_density_1d_dp(ii) = multiplier * nadaraya_watson_1d_dp(z)
       end if
    end do
    
    ! check whether output has to be unpacked
    if ( .not. present(xout) .and. present(mask) ) then
       deallocate( x )
       allocate( x( size(ix,1)))
       x = unpack( kernel_density_1d_dp, mm, nodata )
       deallocate( kernel_density_1d_dp )
       allocate( kernel_density_1d_dp( size(ix,1) ) )
       kernel_density_1d_dp = x
    end if

  end function kernel_density_1d_dp

  function kernel_density_1d_sp(ix,h,silverman,xout,mask,nodata)

    implicit none

    real(sp), dimension(:),                       intent(in) :: ix
    real(sp),                           optional, intent(in) :: h
    logical,                            optional, intent(in) :: silverman
    real(sp), dimension(:),             optional, intent(in) :: xout
    logical,  dimension(:),             optional, intent(in) :: mask
    real(sp),                           optional, intent(in) :: nodata
    real(sp), dimension(:), allocatable                      :: kernel_density_1d_sp

    ! local variables
    logical,  dimension(size(ix,1))            :: mm
    integer(i4)                                :: nn, nout
    integer(i4)                                :: ii
    real(sp)                                   :: hh
    real(sp), dimension(:),        allocatable :: xxout
    real(sp), dimension(:),        allocatable :: z
    real(sp), dimension(:),        allocatable :: x
    real(sp)                                   :: multiplier

    ! consitency check
    if ( .not. present(xout) .and. present(mask) .and. .not. present(nodata) ) then
       print *, ' ERROR *** missing nodata value or xout. function kernel_density'
       stop ' ERROR*** see StdOut for details'
    end if
       
    ! initialize x
    mm = .true.
    if ( present( mask ) ) mm = mask

    allocate( z( count(mm) ) )
    allocate( x( count(mm) ) )
    x = pack( ix, mm)

    nn   = size(x,1)

    ! calc regression
    if (present(xout)) then
       allocate(xxout(size(xout,1)))
       xxout = xout
    else
       allocate(xxout(size(x,1)))
       xxout = x
    end if
    nout   = size(xxout,1)

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_density_h_1d_sp(x,silverman=silverman)
       else
          hh = kernel_density_h_1d_sp(x,silverman=.true.)
       end if
    end if

    ! allocate output
    allocate(kernel_density_1d_sp(nout))

    multiplier = 1.0_sp/(real(nn,sp)*hh)
    ! loop through all regression points
    do ii = 1, nout
       ! scaled difference from regression point
       z(:) = (x(:) - xxout(ii)) / hh
       ! nadaraya-watson estimator of gaussian multivariate kernel
       kernel_density_1d_sp(ii) = nadaraya_watson_1d_sp(z)
       if (kernel_density_1d_sp(ii) .gt. real(nn,sp)*hh*epsilon(1.0_sp)) then
          kernel_density_1d_sp(ii) = multiplier * nadaraya_watson_1d_sp(z)
       end if
    end do

    ! check whether output has to be unpacked
    if ( .not. present(xout) .and. present(mask) ) then
       deallocate( x )
       allocate( x( size(ix,1)))
       x = unpack( kernel_density_1d_sp, mm, nodata )
       deallocate( kernel_density_1d_sp )
       allocate( kernel_density_1d_sp( size(ix,1) ) )
       kernel_density_1d_sp = x
    end if

  end function kernel_density_1d_sp

  ! ------------------------------------------------------------------------------------------------

  function kernel_density_h_1d_dp(x,silverman)

    real(dp), dimension(:),           intent(in) :: x
    logical,                optional, intent(in) :: silverman
    real(dp)                                     :: kernel_density_h_1d_dp

    ! local variables
    integer(i4)              :: nn
    real(dp), dimension(1)   :: h
    real(dp)                 :: stddev_x
    real(dp), dimension(1,2) :: bounds

    nn   = size(x,1)

    ! Silverman's rule of thumb by
    ! Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    stddev_x = stddev(x(:))

    h(1) = (4._dp/3._dp/real(nn,dp))**(0.2_dp) * stddev_x

    if (present(silverman)) then
       if (.not. silverman) then
          bounds(1,1) = max(0.2_dp * h(1),(maxval(x)-minval(x))/real(nn,dp))
          bounds(1,2) = 5.0_dp * h(1)
          call set_globals_for_opti_1d_dp(x) !,xout=xout)
          h = nelminrange(cross_valid_density_h_2d_dp, h, bounds,varmin=0.1_dp)
       end if
    end if

    kernel_density_h_1d_dp = h(1)

  end function kernel_density_h_1d_dp

  function kernel_density_h_1d_sp(x,silverman)

    real(sp), dimension(:),           intent(in) :: x
    logical,                optional, intent(in) :: silverman
    real(sp)                                     :: kernel_density_h_1d_sp

    ! local variables
    integer(i4)              :: nn
    real(sp), dimension(1)   :: h
    real(sp)                 :: stddev_x
    real(sp), dimension(1,2) :: bounds

    nn   = size(x,1)

    ! Silverman's rule of thumb by
    ! Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    stddev_x = stddev(x(:))

    h(1) = (4._sp/3._sp/real(nn,sp))**(0.2_sp) * stddev_x

    if (present(silverman)) then
       if (.not. silverman) then
          bounds(1,1) = max(0.2_sp * h(1),(maxval(x)-minval(x))/real(nn,sp))
          bounds(1,2) = 5.0_sp * h(1)
          call set_globals_for_opti_1d_sp(x) !,xout=xout)
          h = nelminrange(cross_valid_density_h_2d_sp, h, bounds,varmin=0.1_sp)
       end if
    end if

    kernel_density_h_1d_sp = h(1)

  end function kernel_density_h_1d_sp

  ! ------------------------------------------------------------------------------------------------

  function kernel_regression_1d_dp(x,y,h,silverman,xout)

    implicit none

    real(dp), dimension(:),                       intent(in) :: x
    real(dp), dimension(:),                       intent(in) :: y
    real(dp),                           optional, intent(in) :: h
    logical,                            optional, intent(in) :: silverman
    real(dp), dimension(:),             optional, intent(in) :: xout
    real(dp), dimension(:), allocatable                      :: kernel_regression_1d_dp

    ! local variables
    integer(i4)                                :: nn, nout
    integer(i4)                                :: ii
    real(dp)                                   :: hh
    real(dp), dimension(:),        allocatable :: xxout
    real(dp), dimension(size(x,1))             :: z

    nn   = size(x,1)

    ! consistency checks of inputs
    if (size(y) .ne. nn)   stop 'kernel_regression_1d_dp : size of 2nd argument not matching'

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_regression_h_1d_dp(x,y,silverman=silverman)
       else
          hh = kernel_regression_h_1d_dp(x,y,silverman=.true.)
       end if
    end if

    ! calc regression
    if (present(xout)) then
       allocate(xxout(size(xout,1)))
       xxout = xout
    else
       allocate(xxout(size(x,1)))
       xxout = x
    end if
    nout   = size(xxout,1)

    ! allocate output
    allocate(kernel_regression_1d_dp(nout))

    ! loop through all regression points
    do ii = 1, nout
       ! scaled difference from regression point
       z(:) = (x(:) - xxout(ii)) / hh
       ! nadaraya-watson estimator of gaussian multivariate kernel
       kernel_regression_1d_dp(ii) = nadaraya_watson_1d_dp(z,y=y)
    end do

  end function kernel_regression_1d_dp

  function kernel_regression_2d_dp(x,y,h,silverman,xout)

    implicit none

    real(dp), dimension(:,:),                       intent(in) :: x
    real(dp), dimension(:),                         intent(in) :: y
    real(dp), dimension(:),               optional, intent(in) :: h
    logical,                              optional, intent(in) :: silverman
    real(dp), dimension(:,:),             optional, intent(in) :: xout
    real(dp), dimension(:),   allocatable                      :: kernel_regression_2d_dp

    ! local variables
    integer(i4)                                :: dims, dimout
    integer(i4)                                :: nn, nout
    integer(i4)                                :: ii, jj
    real(dp), dimension(size(x,2))             :: hh
    real(dp), dimension(:,:),      allocatable :: xxout
    real(dp), dimension(size(x,1), size(x,2))  :: z


    dims = size(x,2)
    nn   = size(x,1)

    ! consistency checks of inputs
    if (size(y) .ne. nn)   stop 'kernel_regression_2d_dp : size of 2nd argument not matching'
    if (present(h)) then
       if (size(h) .ne. dims) stop 'kernel_regression_2d_dp : size of 3rd argument not matching'
    end if
    if (present(xout)) then
       if (size(xout,2) .ne. dims) stop 'kernel_regression_2d_dp : size of optional argument not matching'
    end if

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_regression_h_2d_dp(x,y,silverman=silverman)
       else
          hh = kernel_regression_h_2d_dp(x,y,silverman=.true.)
       end if
    end if

    ! calc regression
    if (present(xout)) then
       allocate(xxout(size(xout,1),size(xout,2)))
       xxout = xout
    else
       allocate(xxout(size(x,1),size(x,2)))
       xxout = x
    end if
    nout   = size(xxout,1)
    dimout = size(xxout,2)

    ! allocate output
    allocate(kernel_regression_2d_dp(nout))

    ! loop through all regression points
    do ii = 1, nout
       ! scaled difference from regression point
       forall(jj=1:dimout) z(:,jj) = (x(:,jj) - xxout(ii,jj)) / hh(jj)
       ! nadaraya-watson estimator of gaussian multivariate kernel
       kernel_regression_2d_dp(ii) = nadaraya_watson_2d_dp(z,y=y)
    end do

  end function kernel_regression_2d_dp

  function kernel_regression_1d_sp(x,y,h,silverman,xout)

    implicit none

    real(sp), dimension(:),                       intent(in) :: x
    real(sp), dimension(:),                       intent(in) :: y
    real(sp),                           optional, intent(in) :: h
    logical,                            optional, intent(in) :: silverman
    real(sp), dimension(:),             optional, intent(in) :: xout
    real(sp), dimension(:), allocatable                      :: kernel_regression_1d_sp

    ! local variables
    integer(i4)                                :: nn, nout
    integer(i4)                                :: ii
    real(sp)                                   :: hh
    real(sp), dimension(:),        allocatable :: xxout
    real(sp), dimension(size(x,1))             :: z

    nn   = size(x,1)

    ! consistency checks of inputs
    if (size(y) .ne. nn)   stop 'kernel_regression_1d_sp : size of 2nd argument not matching'

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_regression_h_1d_sp(x,y,silverman=silverman)
       else
          hh = kernel_regression_h_1d_sp(x,y,silverman=.true.)
       end if
    end if

    ! calc regression
    if (present(xout)) then
       allocate(xxout(size(xout,1)))
       xxout = xout
    else
       allocate(xxout(size(x,1)))
       xxout = x
    end if
    nout   = size(xxout,1)

    ! allocate output
    allocate(kernel_regression_1d_sp(nout))

    ! loop through all regression points
    do ii = 1, nout
       ! scaled difference from regression point
       z(:) = (x(:) - xxout(ii)) / hh
       ! nadaraya-watson estimator of gaussian multivariate kernel
       kernel_regression_1d_sp(ii) = nadaraya_watson_1d_sp(z,y=y)
    end do

  end function kernel_regression_1d_sp

  function kernel_regression_2d_sp(x,y,h,silverman,xout)

    implicit none

    real(sp), dimension(:,:),                       intent(in) :: x
    real(sp), dimension(:),                         intent(in) :: y
    real(sp), dimension(:),               optional, intent(in) :: h
    logical,                              optional, intent(in) :: silverman
    real(sp), dimension(:,:),             optional, intent(in) :: xout
    real(sp), dimension(:),   allocatable                      :: kernel_regression_2d_sp

    ! local variables
    integer(i4)                                :: dims, dimout
    integer(i4)                                :: nn, nout
    integer(i4)                                :: ii, jj
    real(sp), dimension(size(x,2))             :: hh
    real(sp), dimension(:,:),      allocatable :: xxout
    real(sp), dimension(size(x,1), size(x,2))  :: z


    dims = size(x,2)
    nn   = size(x,1)

    ! consistency checks of inputs
    if (size(y) .ne. nn)   stop 'kernel_regression_2d_sp : size of 2nd argument not matching'
    if (present(h)) then
       if (size(h) .ne. dims) stop 'kernel_regression_2d_sp : size of 3rd argument not matching'
    end if
    if (present(xout)) then
       if (size(xout,2) .ne. dims) stop 'kernel_regression_2d_sp : size of optional argument not matching'
    end if

    ! determine h
    if (present(h)) then
       hh = h
    else
       if (present(silverman)) then
          hh = kernel_regression_h_2d_sp(x,y,silverman=silverman)
       else
          hh = kernel_regression_h_2d_sp(x,y,silverman=.true.)
       end if
    end if

    ! calc regression
    if (present(xout)) then
       allocate(xxout(size(xout,1),size(xout,2)))
       xxout = xout
    else
       allocate(xxout(size(x,1),size(x,2)))
       xxout = x
    end if
    nout   = size(xxout,1)
    dimout = size(xxout,2)

    ! allocate output
    allocate(kernel_regression_2d_sp(nout))

    ! loop through all regression points
    do ii = 1, nout
       ! scaled difference from regression point
       forall(jj=1:dimout) z(:,jj) = (x(:,jj) - xxout(ii,jj)) / hh(jj)
       ! nadaraya-watson estimator of gaussian multivariate kernel
       kernel_regression_2d_sp(ii) = nadaraya_watson_2d_sp(z,y=y)
    end do

  end function kernel_regression_2d_sp

  ! ------------------------------------------------------------------------------------------------

  function kernel_regression_h_1d_dp(x,y,silverman)

    real(dp), dimension(:),           intent(in) :: x
    real(dp), dimension(:),           intent(in) :: y
    logical,                optional, intent(in) :: silverman
    real(dp)                                     :: kernel_regression_h_1d_dp

    ! local variables
    integer(i4)              :: nn
    real(dp), dimension(1)   :: h
    real(dp)                 :: stddev_x
    real(dp), dimension(1,2) :: bounds

    nn   = size(x,1)

    ! Silverman's rule of thumb by
    ! Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    stddev_x = stddev(x(:))

    h(1) = (4._dp/3._dp/real(nn,dp))**(0.2_dp) * stddev_x

    if (present(silverman)) then
       if (.not. silverman) then
          bounds(1,1) = 0.2_dp * h(1)
          bounds(1,2) = 5.0_dp * h(1)
          call set_globals_for_opti_1d_dp(x,y)
          h = nelminrange(cross_valid_regression_h_2d_dp, h, bounds)
       end if
    end if

    kernel_regression_h_1d_dp = h(1)

  end function kernel_regression_h_1d_dp

  function kernel_regression_h_2d_dp(x,y,silverman)

    real(dp), dimension(:,:),                       intent(in) :: x
    real(dp), dimension(:),                         intent(in) :: y
    logical,                              optional, intent(in) :: silverman
    real(dp), dimension(size(x,2))                             :: kernel_regression_h_2d_dp

    ! local variables
    integer(i4)                      :: dims, nn, ii
    real(dp), dimension(size(x,2))   :: h
    real(dp), dimension(size(x,2))   :: stddev_x
    real(dp), dimension(size(x,2),2) :: bounds

    dims = size(x,2)
    nn   = size(x,1)

    ! Silverman's rule of thumb by
    ! Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    do ii=1,dims
       stddev_x(ii) = stddev(x(:,ii))
    end do
    h(:) = (4._dp/real(dims+2,dp)/real(nn,dp))**(1._dp/real(dims+4,dp)) * stddev_x(:)

    if (present(silverman)) then
       if (.not. silverman) then
          bounds(:,1) = 0.2_dp * h(:)
          bounds(:,2) = 5.0_dp * h(:)
          call set_globals_for_opti_2d_dp(x,y)
          h = nelminrange(cross_valid_regression_h_2d_dp, h, bounds)
       end if
    end if

    kernel_regression_h_2d_dp = h

  end function kernel_regression_h_2d_dp

  function kernel_regression_h_1d_sp(x,y,silverman)

    real(sp), dimension(:),           intent(in) :: x
    real(sp), dimension(:),           intent(in) :: y
    logical,                optional, intent(in) :: silverman
    real(sp)                                     :: kernel_regression_h_1d_sp

    ! local variables
    integer(i4)              :: nn
    real(sp), dimension(1)   :: h
    real(sp)                 :: stddev_x
    real(sp), dimension(1,2) :: bounds

    nn   = size(x,1)

    ! Silverman's rule of thumb by
    ! Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    stddev_x = stddev(x(:))

    h(1) = (4._sp/3._sp/real(nn,sp))**(0.2_sp) * stddev_x

    if (present(silverman)) then
       if (.not. silverman) then
          bounds(1,1) = 0.2_sp * h(1)
          bounds(1,2) = 5.0_sp * h(1)
          call set_globals_for_opti_1d_sp(x,y)
          h = nelminrange(cross_valid_regression_h_2d_sp, h, bounds)
       end if
    end if

    kernel_regression_h_1d_sp = h(1)

  end function kernel_regression_h_1d_sp

  function kernel_regression_h_2d_sp(x,y,silverman)

    real(sp), dimension(:,:),                       intent(in) :: x
    real(sp), dimension(:),                         intent(in) :: y
    logical,                              optional, intent(in) :: silverman
    real(sp), dimension(size(x,2))                             :: kernel_regression_h_2d_sp

    ! local variables
    integer(i4)                      :: dims, nn, ii
    real(sp), dimension(size(x,2))   :: h
    real(sp), dimension(size(x,2))   :: stddev_x
    real(sp), dimension(size(x,2),2) :: bounds

    dims = size(x,2)
    nn   = size(x,1)

    ! Silverman's rule of thumb by
    ! Silvermann (1986), Scott (1992), Bowman and Azzalini (1997)
    do ii=1,dims
       stddev_x(ii) = stddev(x(:,ii))
    end do
    h(:) = (4._sp/real(dims+2,sp)/real(nn,sp))**(1._sp/real(dims+4,sp)) * stddev_x(:)
    ! print*, h

    if (present(silverman)) then
       if (.not. silverman) then
          bounds(:,1) = 0.2_sp * h(:)
          bounds(:,2) = 5.0_sp * h(:)
          call set_globals_for_opti_2d_sp(x,y)
          h = nelminrange(cross_valid_regression_h_2d_sp, h, bounds)
       end if
    end if

    kernel_regression_h_2d_sp = h

  end function kernel_regression_h_2d_sp

  ! ------------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------------------------

  !                PRIVATE ROUTINES

  ! ------------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------------------------

  function nadaraya_watson_1d_dp(z,y,mask,valid)

    real(dp), dimension(:),             intent(in)  :: z
    real(dp), dimension(:),   optional, intent(in)  :: y
    logical,  dimension(:),   optional, intent(in)  :: mask
    logical,                  optional, intent(out) :: valid
    real(dp)                                        :: nadaraya_watson_1d_dp

    ! local variables
    real(dp), dimension(size(z,1))            :: w
    real(dp)                                  :: sum_w
    logical,  dimension(size(z,1))            :: mask1d

    if (present(mask)) then
       mask1d = mask
    else
       mask1d = .true.
    end if

    where (mask1d) w = (1.0_dp/sqrt(2.0_dp*Pi_dp)) * exp(-0.5_dp*z*z)
    sum_w = sum(w, mask=mask1d)

    if (present(y)) then
       ! kernel regression
       if (sum_w .lt. epsilon(1.0_dp)) then
          nadaraya_watson_1d_dp = huge(1.0_dp)
          if (present(valid)) valid = .false.
       else
          nadaraya_watson_1d_dp = sum(w*y,mask=mask1d) / sum_w
          if (present(valid)) valid = .true.
       end if
    else
       ! kernel density
       nadaraya_watson_1d_dp = sum_w
       if (present(valid)) valid = .true.
    end if

  end function nadaraya_watson_1d_dp

  function nadaraya_watson_1d_sp(z,y,mask,valid)

    real(sp), dimension(:),             intent(in)  :: z
    real(sp), dimension(:),   optional, intent(in)  :: y
    logical,  dimension(:),   optional, intent(in)  :: mask
    logical,                  optional, intent(out) :: valid
    real(sp)                                        :: nadaraya_watson_1d_sp

    ! local variables
    real(sp), dimension(size(z,1))            :: w
    real(sp)                                  :: sum_w
    logical,  dimension(size(z,1))            :: mask1d

    if (present(mask)) then
       mask1d = mask
    else
       mask1d = .true.
    end if

    where (mask1d) w = (1.0_sp/sqrt(2.0_sp*Pi_sp)) * exp(-0.5_sp*z*z)
    sum_w = sum(w, mask=mask1d)

    if (present(y)) then
       ! kernel regression
       if (sum_w .lt. epsilon(1.0_sp)) then
          nadaraya_watson_1d_sp = huge(1.0_sp)
          if (present(valid)) valid = .false.
       else
          nadaraya_watson_1d_sp = sum(w*y,mask=mask1d) / sum_w
          if (present(valid)) valid = .true.
       end if
    else
       ! kernel density
       nadaraya_watson_1d_sp = sum_w
       if (present(valid)) valid = .true.
    end if

  end function nadaraya_watson_1d_sp

  function nadaraya_watson_2d_dp(z,y,mask,valid)

    real(dp), dimension(:,:),           intent(in)  :: z
    real(dp), dimension(:),   optional, intent(in)  :: y
    logical,  dimension(:),   optional, intent(in)  :: mask
    logical,                  optional, intent(out) :: valid
    real(dp)                                        :: nadaraya_watson_2d_dp

    ! local variables
    integer(i4)                               :: ii
    real(dp), dimension(size(z,1), size(z,2)) :: kerf
    real(dp), dimension(size(z,1))            :: w
    real(dp)                                  :: sum_w
    logical,  dimension(size(z,1))            :: mask1d
    logical,  dimension(size(z,1), size(z,2)) :: mask2d

    if (present(mask)) then
       forall(ii=1:size(z,1)) mask2d(ii,:) = mask(ii)
       mask1d = mask
    else
       mask1d = .true.
       mask2d = .true.
    end if

    where (mask2d) kerf = (1.0_dp/sqrt(2.0_dp*Pi_dp)) * exp(-0.5_dp*z*z)
    w    = product(kerf,dim=2,mask=mask2d)
    sum_w = sum(w, mask=mask1d)

    if (present(y)) then
       ! kernel regression
       if (sum_w .lt. epsilon(1.0_dp)) then
          nadaraya_watson_2d_dp = huge(1.0_dp)
          if (present(valid)) valid = .false.
       else
          nadaraya_watson_2d_dp = sum(w*y,mask=mask1d) / sum_w
          if (present(valid)) valid = .true.
       end if
    else
       ! kernel density
       nadaraya_watson_2d_dp = sum_w
       if (present(valid)) valid = .true.
    end if

  end function nadaraya_watson_2d_dp

  function nadaraya_watson_2d_sp(z,y,mask,valid)

    real(sp), dimension(:,:),           intent(in)  :: z
    real(sp), dimension(:),   optional, intent(in)  :: y
    logical,  dimension(:),   optional, intent(in)  :: mask
    logical,                  optional, intent(out) :: valid
    real(sp)                                        :: nadaraya_watson_2d_sp

    ! local variables
    integer(i4)                               :: ii
    real(sp), dimension(size(z,1), size(z,2)) :: kerf
    real(sp), dimension(size(z,1))            :: w
    real(sp)                                  :: sum_w
    logical,  dimension(size(z,1))            :: mask1d
    logical,  dimension(size(z,1), size(z,2)) :: mask2d

    if (present(mask)) then
       forall(ii=1:size(z,1)) mask2d(ii,:) = mask(ii)
       mask1d = mask
    else
       mask1d = .true.
       mask2d = .true.
    end if

    where (mask2d) kerf = (1.0_sp/sqrt(2.0_sp*Pi_sp)) * exp(-0.5_sp*z*z)
    w    = product(kerf,dim=2,mask=mask2d)
    sum_w = sum(w, mask=mask1d)

    if (present(y)) then
       ! kernel regression
       if (sum_w .lt. epsilon(1.0_sp)) then
          nadaraya_watson_2d_sp = huge(1.0_sp)
          if (present(valid)) valid = .false.
       else
          nadaraya_watson_2d_sp = sum(w*y,mask=mask1d) / sum_w
          if (present(valid)) valid = .true.
       end if
    else
       ! kernel density
       nadaraya_watson_2d_sp = sum_w
       if (present(valid)) valid = .true.
    end if

  end function nadaraya_watson_2d_sp

  ! ------------------------------------------------------------------------------------------------

  function cross_valid_regression_h_2d_dp(h)

    ! Helper function that calculates cross-validation function for the
    ! Nadaraya-Watson estimator, which is basically the mean square error
    ! where model estimate is replaced by the jackknife estimate (Haerdle et al. 2000).

    real(dp), dimension(:), intent(in) :: h
    real(dp)                           :: cross_valid_regression_h_2d_dp

    ! local variables
    integer(i4)                                                    :: ii, jj, kk, nn, dims
    logical,  dimension(size(global_x_dp,1))                       :: mask
    real(dp), dimension(size(global_x_dp,1))                       :: out
    real(dp), dimension(size(global_x_dp,1),size(global_x_dp,2))   :: zz
    logical                                                        :: valid, valid_tmp

    nn   = size(global_x_dp,1)
    dims = size(global_x_dp,2)

    ! Loop through each regression point
    valid = .true.
    do ii=1, nn
       mask = .true.
       mask(ii) = .false.
       forall(jj=1:dims, kk=1:nn, mask(kk)) zz(kk,jj) = (global_x_dp(kk,jj) - global_x_dp(ii,jj)) / h(jj)
       out(ii) = nadaraya_watson_2d_dp(zz, y=global_y_dp, mask=mask, valid=valid_tmp)
       valid = valid .and. valid_tmp
    end do

    if ( valid ) then
       cross_valid_regression_h_2d_dp = sum((global_y_dp-out)**2) / real(nn,dp)
    else
       cross_valid_regression_h_2d_dp = huge(1.0_dp)
    end if

  end function cross_valid_regression_h_2d_dp

  function cross_valid_regression_h_2d_sp(h)

    ! Helper function that calculates cross-validation function for the
    ! Nadaraya-Watson estimator, which is basically the mean square error
    ! where model estimate is replaced by the jackknife estimate (Haerdle et al. 2000).

    real(sp), dimension(:), intent(in) :: h
    real(sp)                           :: cross_valid_regression_h_2d_sp

    ! local variables
    integer(i4)                                                    :: ii, jj, kk, nn, dims
    logical,  dimension(size(global_x_sp,1))                       :: mask
    real(sp), dimension(size(global_x_sp,1))                       :: out
    real(sp), dimension(size(global_x_sp,1),size(global_x_sp,2))   :: zz
    logical                                                        :: valid, valid_tmp

    nn   = size(global_x_sp,1)
    dims = size(global_x_sp,2)

    ! Loop through each regression point
    valid = .true.
    do ii=1, nn
       mask = .true.
       mask(ii) = .false.
       forall(jj=1:dims, kk=1:nn, mask(kk)) zz(kk,jj) = (global_x_sp(kk,jj) - global_x_sp(ii,jj)) / h(jj)
       out(ii) = nadaraya_watson_2d_sp(zz, y=global_y_sp, mask=mask, valid=valid_tmp)
       valid = valid .and. valid_tmp
    end do

    if ( valid ) then
       cross_valid_regression_h_2d_sp = sum((global_y_sp-out)**2) / real(nn,sp)
    else
       cross_valid_regression_h_2d_sp = huge(1.0_sp)
    end if

  end function cross_valid_regression_h_2d_sp

  ! ------------------------------------------------------------------------------------------------

  function cross_valid_density_h_2d_dp(h)

    ! Helper function that calculates cross-validation function for the
    ! Nadaraya-Watson estimator, which is basically the mean square error
    ! where model estimate is replaced by the jackknife estimate (Haerdle et al. 2000).

    real(dp), dimension(:), intent(in) :: h
    real(dp)                           :: cross_valid_density_h_2d_dp

    ! local variables
    integer(i4)                                                     :: ii, jj, kk, nn, dims
    logical,  dimension(size(global_x_dp,1))                        :: mask
    real(dp), dimension(size(global_x_dp,1))                        :: out
    real(dp), dimension(size(global_x_dp,1),size(global_x_dp,2))    :: zz
    real(dp), dimension(size(global_x_dp,2),2)                      :: xMinMax
    real(dp), dimension(size(global_x_dp,2))                        :: delta
    integer(i4)                                                     :: mesh_n
    real(dp), dimension(:,:), allocatable                           :: xMeshed
    real(dp), dimension(:),   allocatable                           :: outIntegral
    real(dp), dimension(size(global_x_dp,1),size(global_x_dp,2))    :: zzIntegral
    real(dp), dimension(size(global_x_dp,2))                        :: stddev_x
    real(dp)                                                        :: summ, multiplier

    nn   = size(global_x_dp,1)
    dims = size(global_x_dp,2)

    if (nn .le. 100_i4) then
       ! if few number of data points given, mesh consists of 100*n points
       mesh_n = 100_i4
    else
       ! mesh_n such that mesh consists of not more than 10000 points
       mesh_n = Max(2_i4, 10000_i4/nn)
    end if
    allocate(xMeshed(mesh_n*size(global_x_dp,1),size(global_x_dp,2)))
    allocate(outIntegral(mesh_n*size(global_x_dp,1)))

    ! integral of squared density function
    do ii=1,dims
       stddev_x(ii) = stddev(global_x_dp(:,ii))
    end do
    forall(jj=1:dims)                 xMinMax(jj,1)  = minval(global_x_dp(:,jj)) - 3.0_dp*stddev_x(jj)
    forall(jj=1:dims)                 xMinMax(jj,2)  = maxval(global_x_dp(:,jj)) + 3.0_dp*stddev_x(jj)
    forall(jj=1:dims)                 delta(jj)      = (xMinMax(jj,2)-xMinMax(jj,1)) / real(nn*(mesh_n-1),dp)
    forall(jj=1:dims, ii=1:nn*mesh_n) xMeshed(ii,jj) = xMinMax(jj,1) + delta(jj) * real(ii-1,dp)

    multiplier = 1.0_dp/(real(nn,dp)*product(h))

    !$OMP parallel default(shared) &
    !$OMP private(zzIntegral, jj)
    !$OMP do
    do ii=1,nn*mesh_n
       forall(jj=1:dims) zzIntegral(:,jj) = (global_x_dp(:,jj) - xMeshed(ii,jj)) / h(jj)
       outIntegral(ii) = nadaraya_watson_2d_dp(zzIntegral) * multiplier
    end do
    !$OMP end do
    !$OMP end parallel

    summ = sum( outIntegral * product(delta) )
    ! print*, 'Integral1: ',summ
    ! scaling to one
    outIntegral = outIntegral / summ
    summ = sum( outIntegral**2 * product(delta) )
    ! print*, 'Integral2: ',summ

    ! Loop through each density point
    !$OMP parallel default(shared) &
    !$OMP private(zzIntegral, jj, kk, mask, zz)
    !$OMP do
    do ii=1, nn
       mask = .true.
       mask(ii) = .false.
       forall(jj=1:dims, kk=1:nn, mask(kk)) zz(kk,jj) = (global_x_dp(kk,jj) - global_x_dp(ii,jj)) / h(jj)
       out(ii) = nadaraya_watson_2d_dp(zz, mask=mask) * multiplier
    end do
    !$OMP end do
    !$OMP end parallel

    cross_valid_density_h_2d_dp = summ - 2.0_dp / (real(nn,dp)) * sum(out)
    ! print*, 'cross_valid_density_h_2d_dp ',h, cross_valid_density_h_2d_dp

  end function cross_valid_density_h_2d_dp

  function cross_valid_density_h_2d_sp(h)

    ! Helper function that calculates cross-validation function for the
    ! Nadaraya-Watson estimator, which is basically the mean square error
    ! where model estimate is replaced by the jackknife estimate (Haerdle et al. 2000).

    real(sp), dimension(:), intent(in) :: h
    real(sp)                           :: cross_valid_density_h_2d_sp

    ! local variables
    integer(i4)                                                     :: ii, jj, kk, nn, dims
    logical,  dimension(size(global_x_sp,1))                        :: mask
    real(sp), dimension(size(global_x_sp,1))                        :: out
    real(sp), dimension(size(global_x_sp,1),size(global_x_sp,2))    :: zz
    real(sp), dimension(size(global_x_sp,2),2)                      :: xMinMax
    real(sp), dimension(size(global_x_sp,2))                        :: delta
    integer(i4)                                                     :: mesh_n
    real(sp), dimension(:,:), allocatable                           :: xMeshed
    real(sp), dimension(:),   allocatable                           :: outIntegral
    real(sp), dimension(size(global_x_sp,1),size(global_x_sp,2))    :: zzIntegral
    real(sp), dimension(size(global_x_sp,2))                        :: stddev_x
    real(sp)                                                        :: summ, multiplier

    nn   = size(global_x_sp,1)
    dims = size(global_x_sp,2)

    if (nn .le. 100_i4) then
       ! if few number of data points given, mesh consists of 100*n points
       mesh_n = 100_i4
    else
       ! mesh_n such that mesh consists of not more than 10000 points
       mesh_n = Max(2_i4, 10000_i4/nn)
    end if
    allocate(xMeshed(mesh_n*size(global_x_sp,1),size(global_x_sp,2)))
    allocate(outIntegral(mesh_n*size(global_x_sp,1)))

    ! integral of squared density function
    do ii=1,dims
       stddev_x(ii) = stddev(global_x_sp(:,ii))
    end do
    forall(jj=1:dims)                 xMinMax(jj,1)  = minval(global_x_sp(:,jj)) - 3.0_sp*stddev_x(jj)
    forall(jj=1:dims)                 xMinMax(jj,2)  = maxval(global_x_sp(:,jj)) + 3.0_sp*stddev_x(jj)
    forall(jj=1:dims)                 delta(jj)      = (xMinMax(jj,2)-xMinMax(jj,1)) / real(nn*(mesh_n-1),sp)
    forall(jj=1:dims, ii=1:nn*mesh_n) xMeshed(ii,jj) = xMinMax(jj,1) + delta(jj) * real(ii-1,sp)

    multiplier = 1.0_sp/(real(nn,sp)*product(h))

    !$OMP parallel default(shared) &
    !$OMP private(zzIntegral, jj)
    !$OMP do
    do ii=1,nn*mesh_n
       forall(jj=1:dims) zzIntegral(:,jj) = (global_x_sp(:,jj) - xMeshed(ii,jj)) / h(jj)
       outIntegral(ii) = nadaraya_watson_2d_sp(zzIntegral) * multiplier
    end do
    !$OMP end do
    !$OMP end parallel

    summ = sum( outIntegral * product(delta) )
    ! print*, 'Integral1: ',summ
    ! scaling to one
    outIntegral = outIntegral / summ
    summ = sum( outIntegral**2 * product(delta) )
    ! print*, 'Integral2: ',summ

    ! Loop through each density point
    !$OMP parallel default(shared) &
    !$OMP private(zzIntegral, jj, kk, mask, zz)
    !$OMP do
    do ii=1, nn
       mask = .true.
       mask(ii) = .false.
       forall(jj=1:dims, kk=1:nn, mask(kk)) zz(kk,jj) = (global_x_sp(kk,jj) - global_x_sp(ii,jj)) / h(jj)
       out(ii) = nadaraya_watson_2d_sp(zz, mask=mask) * multiplier
    end do
    !$OMP end do
    !$OMP end parallel

    cross_valid_density_h_2d_sp = summ - 2.0_sp / (real(nn,sp)) * sum(out)
    ! print*, 'cross_valid_density_h_2d_sp ',h, cross_valid_density_h_2d_sp

  end function cross_valid_density_h_2d_sp

  ! ------------------------------------------------------------------------------------------------

  subroutine set_globals_for_opti_1d_dp(x,y,xout)

    real(dp), dimension(:),           intent(in) :: x
    real(dp), dimension(:), optional, intent(in) :: y
    real(dp), dimension(:), optional, intent(in) :: xout

    if (allocated(global_x_dp)) deallocate(global_x_dp)
    allocate( global_x_dp(size(x,1),1) )
    global_x_dp(:,1) = x

    if (present(y)) then
       if (allocated(global_y_dp)) deallocate(global_y_dp)
       allocate( global_y_dp(size(y,1)) )
       global_y_dp = y
    end if

    if (present(xout)) then
       if (allocated(global_xout_dp)) deallocate(global_xout_dp)
       allocate( global_xout_dp(size(xout,1),1) )
       global_xout_dp(:,1) = xout
    end if

  end subroutine set_globals_for_opti_1d_dp

  subroutine set_globals_for_opti_1d_sp(x,y,xout)

    real(sp), dimension(:),           intent(in) :: x
    real(sp), dimension(:), optional, intent(in) :: y
    real(sp), dimension(:), optional, intent(in) :: xout

    if (allocated(global_x_sp)) deallocate(global_x_sp)
    allocate( global_x_sp(size(x,1),1) )
    global_x_sp(:,1) = x

    if (present(y)) then
       if (allocated(global_y_sp)) deallocate(global_y_sp)
       allocate( global_y_sp(size(y,1)) )
       global_y_sp = y
    end if

    if (present(xout)) then
       if (allocated(global_xout_sp)) deallocate(global_xout_sp)
       allocate( global_xout_sp(size(xout,1),1) )
       global_xout_sp(:,1) = xout
    end if

  end subroutine set_globals_for_opti_1d_sp

  subroutine set_globals_for_opti_2d_dp(x,y,xout)

    real(dp), dimension(:,:),           intent(in) :: x
    real(dp), dimension(:),   optional, intent(in) :: y
    real(dp), dimension(:,:), optional, intent(in) :: xout

    if (allocated(global_x_dp)) deallocate(global_x_dp)
    allocate( global_x_dp(size(x,1),size(x,2)) )
    global_x_dp = x

    if (present(y)) then
       if (allocated(global_y_dp)) deallocate(global_y_dp)
       allocate( global_y_dp(size(y,1)) )
       global_y_dp = y
    end if

    if (present(xout)) then
       if (allocated(global_xout_dp)) deallocate(global_xout_dp)
       allocate( global_xout_dp(size(xout,1),size(xout,2)) )
       global_xout_dp = xout
    end if

  end subroutine set_globals_for_opti_2d_dp

  subroutine set_globals_for_opti_2d_sp(x,y,xout)

    real(sp), dimension(:,:),           intent(in) :: x
    real(sp), dimension(:),   optional, intent(in) :: y
    real(sp), dimension(:,:), optional, intent(in) :: xout

    if (allocated(global_x_sp)) deallocate(global_x_sp)
    allocate( global_x_sp(size(x,1),size(x,2)) )
    global_x_sp = x

    if (present(y)) then
       if (allocated(global_y_sp)) deallocate(global_y_sp)
       allocate( global_y_sp(size(y,1)) )
       global_y_sp = y
    end if

    if (present(xout)) then
       if (allocated(global_xout_sp)) deallocate(global_xout_sp)
       allocate( global_xout_sp(size(xout,1),size(xout,2)) )
       global_xout_sp = xout
    end if

  end subroutine set_globals_for_opti_2d_sp

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

MODULE mo_errormeasures

  ! This module contains routines for the masked calculation of
  ! error measures like MSE, RMSE, BIAS, SSE, NSE, KGE, ...

  ! Note: all except variance and standard deviation are population and not sample moments,
  !       i.e. they are normally divided by n and not (n-1)

  ! Written  Aug 2012,  Matthias Zink
  ! Modified 2012-2018, Juliane Mai, Stephan Thober, Matthias Cuntz

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2012 Matthias Zink, Juliane Mai, Stephan Thober, Matthias Cuntz - mc (at) macu (dot) de
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

  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE

  PUBLIC :: BIAS                         ! Bias
  PUBLIC :: KGE                          ! Kling-Gupta efficiency measure
  PUBLIC :: LNNSE                        ! Logarithmic Nash Sutcliffe efficiency
  PUBLIC :: MAE                          ! Mean of absolute errors
  PUBLIC :: MAE_PROB_ONE ! Mean absolute error of occurence probability with ONE cdf for obs and mod
  ! comment until mo_empcdf is commited
  ! PUBLIC :: MAE_PROB_TWO ! Mean absolute error of occurence probability with TWO separate cdfs for obs and mod
  PUBLIC :: MSE                          ! Mean of squared errors
  PUBLIC :: NSE                          ! Nash Sutcliffe efficiency
  PUBLIC :: SSE                          ! Sum of squared errors
  PUBLIC :: SAE                          ! Sum of absolute errors
  PUBLIC :: RMSE                         ! Root mean squared error

  ! ------------------------------------------------------------------

  !     NAME
  !         BIAS

  !     PURPOSE
  !         Calculates the bias
  !             BIAS = mean(y) - mean(x)
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = BIAS(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: BIAS        bias

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = BIAS(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE BIAS
     MODULE PROCEDURE BIAS_1d_sp, BIAS_1d_dp, BIAS_2d_sp, BIAS_2d_dp, BIAS_3d_sp, BIAS_3d_dp
  END INTERFACE BIAS

  ! ------------------------------------------------------------------

  !      NAME
  !          KGE

  !>        \brief Kling-Gupta-Efficiency measure.

  !>        \details The Kling-Gupta model efficiency coefficient \f$ KGE \f$ is
  !>                     \f[ KGE = 1 - \sqrt{( (1-r)^2 + (1-\alpha)^2 + (1-\beta)^2 )} \f]
  !>                 where \n
  !>                     \f$ r \f$      = Pearson product-moment correlation coefficient \n
  !>                     \f$ \alpha \f$ = ratio of simulated mean to observed mean  \n
  !>                     \f$ \beta  \f$ = ratio of simulated standard deviation to
  !>                                      observed standard deviation \n
  !>                 This three measures are calculated between two arrays (1d, 2d, or 3d).
  !>                 Usually, one is an observation and the second is a modelled variable.\n
  !>
  !>                 The higher the KGE the better the observation and simulation are matching.
  !>                 The upper limit of KGE is 1.\n
  !>
  !>                 Therefore, if you apply a minimization algorithm to calibrate regarding
  !>                 KGE you have to use the objective function
  !>                     \f[ obj\_value = 1.0 - KGE \f]
  !>                 which has then the optimum at 0.0.
  !>                 (Like for the NSE where you always optimize 1-NSE.)\n
  !>

  !     INTENT(IN)
  !>        real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !>        real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !>        real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !>        logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !>        logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return  kge &mdash; Kling-Gupta-Efficiency (value less equal 1.0)

  !     RESTRICTIONS
  !>       \note Input values must be floating points. \n

  !     EXAMPLE
  !         para = (/ 1., 2, 3., -999., 5., 6. /)
  !         kge = kge(x,y,mask=mask)

  !     LITERATURE
  !>        Gupta, Hoshin V., et al.
  !>           "Decomposition of the mean squared error and NSE performance criteria:
  !>           Implications for improving hydrological modelling."
  !>           Journal of Hydrology 377.1 (2009): 80-91.


  !     HISTORY
  !>        \author Rohini Kumar
  !>        \date August 2014
  !         Modified, R. Kumar & O. Rakovec - Sep. 2014
  !                   J. Mai                - remove double packing of input data (bug)
  !                                         - KGE instead of 1.0-KGE
  !                                         - 1d, 2d, 3d, version in sp and dp

  INTERFACE KGE
     MODULE PROCEDURE KGE_1d_dp, KGE_2d_dp, KGE_3d_dp, KGE_1d_sp, KGE_2d_sp, KGE_3d_sp
  END INTERFACE KGE

  ! ------------------------------------------------------------------

  !     NAME
  !         LNNSE

  !     PURPOSE
  !         Calculates the Logarithmic Nash Sutcliffe Efficiency
  !             LNNSE = sum((ln(y) - ln(x))**2) / sum( (ln(x) - ln(mean(x)))**2 )
  !         where x is the observation and y is the modelled data.
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         Note that the mask is intent inout, since values which are less or equal zero will be masked additionally.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = LNNSE(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: LNNSE         Logarithmic Nash Sutcliffe Efficiency

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.
  !         The mask will be updated if non-masked values are less equal zero.

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = LNNSE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Juliane Mai, May 2013
  !         updated,  Rohin Kumar, May 2013  ! for mean of logQ
  INTERFACE LNNSE
     MODULE PROCEDURE LNNSE_1d_sp, LNNSE_1d_dp, LNNSE_2d_dp, LNNSE_2d_sp, LNNSE_3d_sp, LNNSE_3d_dp
  END INTERFACE LNNSE

  ! ------------------------------------------------------------------

  !     NAME
  !         MAE

  !     PURPOSE
  !         Calculates the mean absolute error
  !             MAE = sum(abs(y - x)) / count(mask)
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = MAE(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: MAE         Mean Absolute Error

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = MAE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE MAE
     MODULE PROCEDURE MAE_1d_sp, MAE_1d_dp, MAE_2d_sp, MAE_2d_dp, MAE_3d_sp, MAE_3d_dp
  END INTERFACE MAE

  ! ! ------------------------------------------------------------------

  ! !     NAME
  ! !         MAE_PROB_TWO

  ! !     PURPOSE
  ! !         Calculate mean absolute error of occurence probabilities
  ! !             MAE_PROB_TWO = mean( |P_sim(Q_obs(t)) - P_obs(Q_obs(t))| )
  ! !
  ! !         If an optional mask is given, the calculations are over those locations that correspond to true values in the mask.
  ! !         x and y have to be double precision. The result will have the same numerical precision.

  ! !     CALLING SEQUENCE
  ! !         out = MAE_PROB_TWO(dat, mask=mask)

  ! !     INTENT(IN)
  ! !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers

  ! !     INTENT(INOUT)
  ! !         None

  ! !     INTENT(OUT)
  ! !         real(sp/dp) :: BIAS        bias

  ! !     INTENT(IN), OPTIONAL
  ! !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  ! !
  ! !         If present, only those locations in vec corresponding to the true values in mask are used.

  ! !     INTENT(INOUT), OPTIONAL
  ! !         None

  ! !     INTENT(OUT), OPTIONAL
  ! !         None

  ! !     RESTRICTIONS
  ! !         Input values must be floating points.

  ! !     EXAMPLE
  ! !         vec1 = (/ 1., 2, 3., -9999., 5., 6. /)
  ! !         vec2 = (/ 1., 2, 3., -9999., 5., 6. /)
  ! !         m   = MAE_PROB_TWO(vec1, vec2, )
  ! !         -> see also example in test directory

  ! !     LITERATURE
  ! !         None

  ! !     HISTORY
  ! !         Written,  Stephan Thober, May 2016
  ! INTERFACE MAE_PROB_TWO
  !    MODULE PROCEDURE MAE_PROB_TWO_1D_DP
  ! END INTERFACE MAE_PROB_TWO

  ! ------------------------------------------------------------------

  !     NAME
  !         MAE_PROB_ONE

  !     PURPOSE
  !         Calculate mean absolute error of occurence probabilities
  !             MAE_PROB_ONE = mean( |P_obs(Q_obs(t)) - P_obs(Q_sim(t))| )
  !
  !         If an optional mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y have to be double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = MAE_PROB_ONE(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: BIAS        bias

  !     INTENT(IN), OPTIONAL
  !         logical                   :: mask(:)     1D-array of logical values with size(x/y).
  !                                                  If present, only those locations in vec corresponding to
  !                                                  the true values in mask are used.
  !         real(sp/dp), dimension(:) :: cdfx        1D-array with kernel_cumdensity of x
  !         real(sp/dp)               :: h           Silverman estimate of kernel_density bandwidth for x
  !

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !           m   = MAE_PROB(x, y, mask=mask)
  !         same as
  !           h    = kernel_density_h(x, silverman=.true., mask=mask)
  !           cdfx = kernel_cumdensity(x, xout=x, h=h, mask=mask)
  !           m    = MAE_PROB(vec1, vec2, mask=vec1, cdfx=cdf, h=h)
  !         but cumdensity of x is calculated outside and h is reused.
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Stephan Thober, May 2016
  !         Modified, Matthias Cuntz, Jun 2016 - rm dummy=pack(x), cdfx, h
  INTERFACE MAE_PROB_ONE
     MODULE PROCEDURE MAE_PROB_ONE_1D_DP
  END INTERFACE MAE_PROB_ONE

  ! ------------------------------------------------------------------

  !     NAME
  !         MSE

  !     PURPOSE
  !         Calculates the mean squared error
  !             MSE = sum((y - x)**2) / count(mask)
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = MSE(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: MSE         Mean squared error

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = MSE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE MSE
     MODULE PROCEDURE MSE_1d_sp, MSE_1d_dp, MSE_2d_sp, MSE_2d_dp, MSE_3d_sp, MSE_3d_dp
  END INTERFACE MSE

  ! ------------------------------------------------------------------

  !     NAME
  !         NSE

  !     PURPOSE
  !         Calculates the Nash Sutcliffe Efficiency
  !             NSE = sum((y - x)**2) / sum( (x - mean(x))**2)
  !         where x is the observation and y is the modelled data.
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = NSE(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: NSE         Nash Sutcliffe Efficiency

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = NSE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         NASH, J., & SUTCLIFFE, J. (1970). River flow forecasting through conceptual models part I: A discussion of
  !                  principles. Journal of Hydrology, 10(3), 282-290. doi:10.1016/0022-1694(70)90255-6

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE NSE
     MODULE PROCEDURE NSE_1d_sp, NSE_1d_dp, NSE_2d_dp, NSE_2d_sp, NSE_3d_sp, NSE_3d_dp
  END INTERFACE NSE

  ! ------------------------------------------------------------------

  !     NAME
  !         SAE

  !     PURPOSE
  !         Calculates the sum of absolute errors
  !             SAE = sum(abs(y - x))
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = SAE(x, y, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: SAE         sum of absolute errors

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = SAE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         none

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE SAE
     MODULE PROCEDURE SAE_1d_sp, SAE_1d_dp, SAE_2d_sp, SAE_2d_dp, SAE_3d_sp, SAE_3d_dp
  END INTERFACE SAE

  ! ------------------------------------------------------------------

  !     NAME
  !         SSE

  !     PURPOSE
  !         Calculates the sum of squared errors
  !             SSE = sum((y - x)**2)
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = SSE(x, y, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: SSE         sum of squared errors

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = SSE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         none

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE SSE
     MODULE PROCEDURE SSE_1d_sp, SSE_1d_dp, SSE_2d_sp, SSE_2d_dp, SSE_3d_sp, SSE_3d_dp
  END INTERFACE SSE

  ! ------------------------------------------------------------------

  !     NAME
  !         RMSE

  !     PURPOSE
  !         Calculates the root-mean-square error
  !             RMSE = sqrt(sum((y - x)**2) / count(mask))
  !
  !         If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !         x and y can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = RMSE(dat, mask=mask)

  !     INTENT(IN)
  !         real(sp/dp), dimension(:)     :: x, y    1D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:)   :: x, y    2D-array with input numbers
  !             OR
  !         real(sp/dp), dimension(:,:,:) :: x, y    3D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: RMSE        Root-mean-square error

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)     1D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:)   2D-array of logical values with size(x/y).
  !             OR
  !         logical     :: mask(:,:,:) 3D-array of logical values with size(x/y).
  !
  !         If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec1 = (/ 1., 2, 3., -999., 5., 6. /)
  !         vec2 = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = RMSE(vec1, vec2, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Zink, Sept 2012
  INTERFACE RMSE
     MODULE PROCEDURE RMSE_1d_sp, RMSE_1d_dp, RMSE_2d_sp, RMSE_2d_dp, RMSE_3d_sp, RMSE_3d_dp
  END INTERFACE RMSE

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION BIAS_1d_sp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                      :: BIAS_1d_sp

    INTEGER(i4)                                   :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )       :: shapemask
    LOGICAL,  DIMENSION(size(x))                  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'BIAS_1d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    !
    if (n .LE. 1_i4) stop 'BIAS_1d_sp: number of arguments must be at least 2'
    !
    BIAS_1d_sp = average(y, mask=maske) - average(x, mask=maske)

  END FUNCTION BIAS_1d_sp

  FUNCTION BIAS_1d_dp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                      :: BIAS_1d_dp

    INTEGER(i4)                                   :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )       :: shapemask
    LOGICAL,  DIMENSION(size(x))                  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'BIAS_1d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'BIAS_1d_dp: number of arguments must be at least 2'
    !
    BIAS_1d_dp = average(y, mask=maske) - average(x, mask=maske)

  END FUNCTION BIAS_1d_dp

  FUNCTION BIAS_2d_sp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: BIAS_2d_sp

    INTEGER(i4)                                       :: n

    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL, DIMENSION(size(x, dim=1), size(x, dim=2)):: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'BIAS_2d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2)
    endif
    !
    if (n .LE. 1_i4) stop 'BIAS_2d_sp: number of arguments must be at least 2'
    !
    BIAS_2d_sp = average(reshape(y,    (/size(y,dim=1) * size(y,dim=2)/)),      &
         mask=reshape(maske,(/size(y,dim=1) * size(y,dim=2)/)))   -  &
         average(reshape(x,    (/size(x,dim=1) * size(x,dim=2)/)),      &
         mask=reshape(maske,(/size(x,dim=1) * size(x,dim=2)/)))
    !
  END FUNCTION BIAS_2d_sp

  FUNCTION BIAS_2d_dp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: BIAS_2d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL, DIMENSION(size(x, dim=1), size(x, dim=2)):: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'BIAS_2d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2)
    endif
    !
    if (n .LE. 1_i4) stop 'BIAS_2d_dp: number of arguments must be at least 2'
    !
    BIAS_2d_dp = average(reshape(y,    (/size(y,dim=1) * size(y,dim=2)/)),      &
         mask=reshape(maske,(/size(y,dim=1) * size(y,dim=2)/)))   -  &
         average(reshape(x,    (/size(x,dim=1) * size(x,dim=2)/)),      &
         mask=reshape(maske,(/size(x,dim=1) * size(x,dim=2)/)))
    !
  END FUNCTION BIAS_2d_dp

  FUNCTION BIAS_3d_sp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: BIAS_3d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1), &
         size(x, dim=2), size(x, dim=3))    :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'BIAS_3d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2) * size(x, dim=3)
    endif
    !
    ! not really sopisticated, it has to be checked if the 3 numbers of x and y are matching in arry position
    if (n .LE. 1_i4) stop 'BIAS_3d_sp: number of arguments must be at least 2'
    !
    BIAS_3d_sp = average(reshape(y,    (/size(y,dim=1) * size(y,dim=2) * size(y,dim=3)/)),      &
         mask=reshape(maske,(/size(y,dim=1) * size(y,dim=2) * size(y,dim=3)/)))   -  &
         average(reshape(x,    (/size(x,dim=1) * size(x,dim=2) * size(x,dim=3)/)),      &
         mask=reshape(maske,(/size(x,dim=1) * size(x,dim=2) * size(x,dim=3)/)))
    !
  END FUNCTION BIAS_3d_sp

  FUNCTION BIAS_3d_dp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: BIAS_3d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1), &
         size(x, dim=2), size(x, dim=3))    :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'BIAS_3d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2) * size(x, dim=3)
    endif
    !
    ! not really sopisticated, it has to be checked if the 3 numbers of x and y are matching in arry position
    if (n .LE. 1_i4) stop 'BIAS_3d_dp: number of arguments must be at least 2'
    !
    BIAS_3d_dp = average(reshape(y,    (/size(y,dim=1) * size(y,dim=2) * size(y,dim=3)/)),      &
         mask=reshape(maske,(/size(y,dim=1) * size(y,dim=2) * size(y,dim=3)/)))   -  &
         average(reshape(x,    (/size(x,dim=1) * size(x,dim=2) * size(x,dim=3)/)),      &
         mask=reshape(maske,(/size(x,dim=1) * size(x,dim=2) * size(x,dim=3)/)))
    !
  END FUNCTION BIAS_3d_dp

  ! ------------------------------------------------------------------

  FUNCTION KGE_1d_sp(x, y, mask)

    USE mo_moment, ONLY: average, stddev, correlation

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                          :: KGE_1d_sp

    ! local variables
    INTEGER(i4)                             :: n
    INTEGER(i4), DIMENSION(size(shape(x)) ) :: shapemask
    LOGICAL,     DIMENSION(size(x))         :: maske

    REAL(sp)                                :: mu_Obs, mu_Sim       ! Mean          of x and y
    REAL(sp)                                :: sigma_Obs, sigma_Sim ! Standard dev. of x and y
    REAL(sp)                                :: pearson_coor         ! Pearson Corr. of x and y

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'KGE_1d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'KGE_1d_sp: sample size must be at least 2'

    ! Mean
    mu_Obs = average(x, mask=maske)
    mu_Sim = average(y, mask=maske)
    ! Standard Deviation
    sigma_Obs = stddev(x, mask=maske, ddof=1_i4)
    sigma_Sim = stddev(y, mask=maske, ddof=1_i4)
    ! Pearson product-moment correlation coefficient
    pearson_coor = correlation(x, y, mask=maske, ddof=1_i4)
    !
    KGE_1d_sp = 1.0 - SQRT( &
         ( 1.0_sp - (mu_Sim/mu_Obs)       )**2 + &
         ( 1.0_sp - (sigma_Sim/sigma_Obs) )**2 + &
         ( 1.0_sp - pearson_coor)**2             &
         )

  END FUNCTION KGE_1d_sp

  FUNCTION KGE_2d_sp(x, y, mask)

    USE mo_moment, ONLY: average, stddev, correlation

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                            :: KGE_2d_sp

    ! local variables
    INTEGER(i4)                                            :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )                :: shapemask
    LOGICAL,     DIMENSION(size(x, dim=1), size(x, dim=2)) :: maske
    REAL(sp)                                               :: mu_Obs, mu_Sim       ! Mean          of x and y
    REAL(sp)                                               :: sigma_Obs, sigma_Sim ! Standard dev. of x and y
    REAL(sp)                                               :: pearson_coor         ! Pearson Corr. of x and y

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'KGE_2d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'KGE_2d_sp: sample size must be at least 2'

    ! Mean
    mu_Obs = average( &
         reshape(x(:,:), (/size(x, dim=1)*size(x, dim=2)/)), &
         mask=reshape(maske(:,:),  (/size(x, dim=1)*size(x, dim=2)/)))
    mu_Sim = average( &
         reshape(y(:,:), (/size(y, dim=1)*size(y, dim=2)/)), &
         mask=reshape(maske(:,:),  (/size(y, dim=1)*size(y, dim=2)/)))
    ! Standard Deviation
    sigma_Obs = stddev( &
         reshape(x(:,:), (/size(x, dim=1)*size(x, dim=2)/)), &
         mask=reshape(maske(:,:),  (/size(x, dim=1)*size(x, dim=2)/)), ddof=1_i4)
    sigma_Sim = stddev( &
         reshape(y(:,:), (/size(y, dim=1)*size(y, dim=2)/)), &
         mask=reshape(maske(:,:),  (/size(y, dim=1)*size(y, dim=2)/)), ddof=1_i4)
    ! Pearson product-moment correlation coefficient
    pearson_coor = correlation(&
         reshape(x(:,:), (/size(x, dim=1)*size(x, dim=2)/)), &
         reshape(y(:,:), (/size(y, dim=1)*size(y, dim=2)/)), &
         mask=reshape(maske(:,:),  (/size(y, dim=1)*size(y, dim=2)/)), ddof=1_i4)
    !
    KGE_2d_sp = 1.0 - SQRT( &
         ( 1.0_sp - (mu_Sim/mu_Obs)       )**2 + &
         ( 1.0_sp - (sigma_Sim/sigma_Obs) )**2 + &
         ( 1.0_sp - pearson_coor)**2             &
         )

  END FUNCTION KGE_2d_sp

  FUNCTION KGE_3d_sp(x, y, mask)

    USE mo_moment, ONLY: average, stddev, correlation

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                            :: KGE_3d_sp

    ! local variables
    INTEGER(i4)                                                            :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )                                :: shapemask
    LOGICAL,     DIMENSION(size(x, dim=1), size(x, dim=2), size(x, dim=3)) :: maske
    REAL(sp)                                                               :: mu_Obs, mu_Sim       ! Mean          of x and y
    REAL(sp)                                                               :: sigma_Obs, sigma_Sim ! Standard dev. of x and y
    REAL(sp)                                                               :: pearson_coor         ! Pearson Corr. of x and y

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'KGE_3d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'KGE_3d_sp: sample size must be at least 2'

    ! Mean
    mu_Obs = average( &
         reshape(x(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), &
         mask=reshape(maske(:,:,:),  (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)))
    mu_Sim = average( &
         reshape(y(:,:,:), (/size(y, dim=1)*size(y, dim=2)*size(y, dim=3)/)), &
         mask=reshape(maske(:,:,:),  (/size(y, dim=1)*size(y, dim=2)*size(y, dim=3)/)))
    ! Standard Deviation
    sigma_Obs = stddev( &
         reshape(x(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), &
         mask=reshape(maske(:,:,:),  (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), ddof=1_i4)
    sigma_Sim = stddev( &
         reshape(y(:,:,:), (/size(y, dim=1)*size(y, dim=2)*size(y, dim=3)/)), &
         mask=reshape(maske(:,:,:),  (/size(y, dim=1)*size(y, dim=2)*size(y, dim=3)/)), ddof=1_i4)
    ! Pearson product-moment correlation coefficient
    pearson_coor = correlation(&
         reshape(x(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), &
         reshape(y(:,:,:), (/size(y, dim=1)*size(y, dim=2)*size(y, dim=3)/)), &
         mask=reshape(maske(:,:,:),  (/size(y, dim=1)*size(y, dim=2)*size(y, dim=3)/)), ddof=1_i4)
    !
    KGE_3d_sp = 1.0 - SQRT( &
         ( 1.0_sp - (mu_Sim/mu_Obs)       )**2 + &
         ( 1.0_sp - (sigma_Sim/sigma_Obs) )**2 + &
         ( 1.0_sp - pearson_coor)**2             &
         )

  END FUNCTION KGE_3d_sp

  FUNCTION KGE_1d_dp(x, y, mask)

    USE mo_moment, ONLY: average, stddev, correlation

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                          :: KGE_1d_dp

    ! local variables
    INTEGER(i4)                             :: n
    INTEGER(i4), DIMENSION(size(shape(x)) ) :: shapemask
    LOGICAL,     DIMENSION(size(x))         :: maske

    REAL(dp)                                :: mu_Obs, mu_Sim       ! Mean          of x and y
    REAL(dp)                                :: sigma_Obs, sigma_Sim ! Standard dev. of x and y
    REAL(dp)                                :: pearson_coor         ! Pearson Corr. of x and y

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'KGE_1d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'KGE_1d_dp: sample size must be at least 2'

    ! Mean
    mu_Obs = average(x, mask=maske)
    mu_Sim = average(y, mask=maske)
    ! Standard Deviation
    sigma_Obs = stddev(x, mask=maske, ddof=1_i4)
    sigma_Sim = stddev(y, mask=maske, ddof=1_i4)
    ! Pearson product-moment correlation coefficient
    pearson_coor = correlation(x, y, mask=maske, ddof=1_i4)
    !
    KGE_1d_dp = 1.0 - SQRT( &
         ( 1.0_dp - (mu_Sim/mu_Obs)       )**2 + &
         ( 1.0_dp - (sigma_Sim/sigma_Obs) )**2 + &
         ( 1.0_dp - pearson_coor)**2             &
         )

  END FUNCTION KGE_1d_dp

  FUNCTION KGE_2d_dp(x, y, mask)

    USE mo_moment, ONLY: average, stddev, correlation

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                            :: KGE_2d_dp

    ! local variables
    INTEGER(i4)                                            :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )                :: shapemask
    LOGICAL,     DIMENSION(size(x, dim=1), size(x, dim=2)) :: maske
    REAL(dp)                                               :: mu_Obs, mu_Sim       ! Mean          of x and y
    REAL(dp)                                               :: sigma_Obs, sigma_Sim ! Standard dev. of x and y
    REAL(dp)                                               :: pearson_coor         ! Pearson Corr. of x and y

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'KGE_2d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'KGE_2d_dp: sample size must be at least 2'

    ! Mean
    mu_Obs = average( &
         reshape(x(:,:), (/size(x, dim=1)*size(x, dim=2)/)), &
         mask=reshape(maske(:,:),  (/size(x, dim=1)*size(x, dim=2)/)))
    mu_Sim = average( &
         reshape(y(:,:), (/size(y, dim=1)*size(y, dim=2)/)), &
         mask=reshape(maske(:,:),  (/size(y, dim=1)*size(y, dim=2)/)))
    ! Standard Deviation
    sigma_Obs = stddev( &
         reshape(x(:,:), (/size(x, dim=1)*size(x, dim=2)/)), &
         mask=reshape(maske(:,:),  (/size(x, dim=1)*size(x, dim=2)/)), ddof=1_i4)
    sigma_Sim = stddev( &
         reshape(y(:,:), (/size(y, dim=1)*size(y, dim=2)/)), &
         mask=reshape(maske(:,:),  (/size(y, dim=1)*size(y, dim=2)/)), ddof=1_i4)
    ! Pearson product-moment correlation coefficient
    pearson_coor = correlation(&
         reshape(x(:,:), (/size(x, dim=1)*size(x, dim=2)/)), &
         reshape(y(:,:), (/size(y, dim=1)*size(y, dim=2)/)), &
         mask=reshape(maske(:,:),  (/size(y, dim=1)*size(y, dim=2)/)), ddof=1_i4)
    !
    KGE_2d_dp = 1.0 - SQRT( &
         ( 1.0_dp - (mu_Sim/mu_Obs)       )**2 + &
         ( 1.0_dp - (sigma_Sim/sigma_Obs) )**2 + &
         ( 1.0_dp - pearson_coor)**2             &
         )

  END FUNCTION KGE_2d_dp

  FUNCTION KGE_3d_dp(x, y, mask)

    USE mo_moment, ONLY: average, stddev, correlation

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                            :: KGE_3d_dp

    ! local variables
    INTEGER(i4)                                                            :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )                                :: shapemask
    LOGICAL,     DIMENSION(size(x, dim=1), size(x, dim=2), size(x, dim=3)) :: maske
    REAL(dp)                                                               :: mu_Obs, mu_Sim       ! Mean          of x and y
    REAL(dp)                                                               :: sigma_Obs, sigma_Sim ! Standard dev. of x and y
    REAL(dp)                                                               :: pearson_coor         ! Pearson Corr. of x and y

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'KGE_3d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'KGE_3d_dp: sample size must be at least 2'

    ! Mean
    mu_Obs = average( &
         reshape(x(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), &
         mask=reshape(maske(:,:,:),  (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)))
    mu_Sim = average( &
         reshape(y(:,:,:), (/size(y, dim=1)*size(y, dim=2)*size(y, dim=3)/)), &
         mask=reshape(maske(:,:,:),  (/size(y, dim=1)*size(y, dim=2)*size(y, dim=3)/)))
    ! Standard Deviation
    sigma_Obs = stddev( &
         reshape(x(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), &
         mask=reshape(maske(:,:,:),  (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), ddof=1_i4)
    sigma_Sim = stddev( &
         reshape(y(:,:,:), (/size(y, dim=1)*size(y, dim=2)*size(y, dim=3)/)), &
         mask=reshape(maske(:,:,:),  (/size(y, dim=1)*size(y, dim=2)*size(y, dim=3)/)), ddof=1_i4)
    ! Pearson product-moment correlation coefficient
    pearson_coor = correlation(&
         reshape(x(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), &
         reshape(y(:,:,:), (/size(y, dim=1)*size(y, dim=2)*size(y, dim=3)/)), &
         mask=reshape(maske(:,:,:),  (/size(y, dim=1)*size(y, dim=2)*size(y, dim=3)/)), ddof=1_i4)
    !
    KGE_3d_dp = 1.0 - SQRT( &
         ( 1.0_dp - (mu_Sim/mu_Obs)       )**2 + &
         ( 1.0_dp - (sigma_Sim/sigma_Obs) )**2 + &
         ( 1.0_dp - pearson_coor)**2             &
         )

  END FUNCTION KGE_3d_dp

  ! ------------------------------------------------------------------

  FUNCTION LNNSE_1d_sp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(INOUT)   :: mask
    REAL(sp)                                          :: LNNSE_1d_sp

    INTEGER(i4)                            :: n
    INTEGER(i4), DIMENSION(size(shape(x))) :: shapemask
    REAL(sp)                               :: xmean
    REAL(sp), DIMENSION(size(x))           :: logx, logy, v1, v2
    LOGICAL,  DIMENSION(size(x))           :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'LNNSE_1d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif

    ! mask all negative and zero entries
    where (x .lt. tiny(1.0_sp) .or. y .lt. tiny(1.0_sp))
       maske = .false.
    end where
    n = count(maske)
    if (n .LE. 1_i4) stop 'LNNSE_1d_sp: number of arguments must be at least 2'

    ! logarithms
    logx = 0.0_sp
    logy = 0.0_sp
    where (maske)
       logx = log(x)
       logy = log(y)
    end where

    ! mean of x
    xmean = average(logx, mask=maske)

    ! NSE
    v1 = merge(logy - logx,  0.0_sp, maske)
    v2 = merge(logx - xmean, 0.0_sp, maske)
    LNNSE_1d_sp = 1.0_sp - dot_product(v1,v1) / dot_product(v2,v2)

  END FUNCTION LNNSE_1d_sp

  ! ------------------------------------------------------------------

  FUNCTION LNNSE_1d_dp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(INOUT)   :: mask
    REAL(dp)                                          :: LNNSE_1d_dp

    INTEGER(i4)                            :: n
    INTEGER(i4), DIMENSION(size(shape(x))) :: shapemask
    REAL(dp)                               :: xmean
    REAL(dp), DIMENSION(size(x))           :: logx, logy, v1, v2
    LOGICAL,  DIMENSION(size(x))           :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'LNNSE_1d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif

    ! mask all negative and zero entries
    where (x .lt. tiny(1.0_dp) .or. y .lt. tiny(1.0_dp))
       maske = .false.
    end where
    n = count(maske)
    if (n .LE. 1_i4) stop 'LNNSE_1d_dp: number of arguments must be at least 2'

    ! logarithms
    logx = 0.0_dp
    logy = 0.0_dp
    where (maske)
       logx = log(x)
       logy = log(y)
    end where

    ! mean of x
    xmean = average(logx, mask=maske)

    ! NSE
    v1 = merge(logy - logx,  0.0_dp, maske)
    v2 = merge(logx - xmean, 0.0_dp, maske)
    LNNSE_1d_dp = 1.0_dp - sum(v1*v1, mask=maske) / sum(v2*v2, mask=maske)

  END FUNCTION LNNSE_1d_dp

  ! ------------------------------------------------------------------

  FUNCTION LNNSE_2d_sp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: mask
    REAL(sp)                                          :: LNNSE_2d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)))            :: shapemask
    REAL(sp)                                          :: xmean
    REAL(sp), DIMENSION(size(x,dim=1),size(x,dim=2))  :: logx, logy, v1, v2
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'LNNSE_2d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif

    ! mask all negative and zero entries
    where (x .lt. tiny(1.0_sp) .or. y .lt. tiny(1.0_sp))
       maske = .false.
    end where
    n = count(maske)
    if (n .LE. 1_i4) stop 'LNNSE_2d_sp: number of arguments must be at least 2'

    ! logarithms
    logx = 0.0_sp
    logy = 0.0_sp
    where (maske)
       logx = log(x)
       logy = log(y)
    end where

    ! mean of x
    xmean = average(pack(logx,maske))

    ! NSE
    v1 = merge(logy - logx,  0.0_sp, maske)
    v2 = merge(logx - xmean, 0.0_sp, maske)
    LNNSE_2d_sp = 1.0_sp - sum(v1*v1, mask=maske) / sum(v2*v2, mask=maske)

  END FUNCTION LNNSE_2d_sp

  ! ------------------------------------------------------------------

  FUNCTION LNNSE_2d_dp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: mask
    REAL(dp)                                          :: LNNSE_2d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)))            :: shapemask
    REAL(dp)                                          :: xmean
    REAL(dp), DIMENSION(size(x,dim=1),size(x,dim=2))  :: logx, logy, v1, v2
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'LNNSE_2d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif

    ! mask all negative and zero entries
    where (x .lt. tiny(1.0_dp) .or. y .lt. tiny(1.0_dp))
       maske = .false.
    end where
    n = count(maske)
    if (n .LE. 1_i4) stop 'LNNSE_2d_dp: number of arguments must be at least 2'

    ! logarithms
    logx = 0.0_dp
    logy = 0.0_dp
    where (maske)
       logx = log(x)
       logy = log(y)
    end where

    ! mean of x
    xmean = average(pack(logx,maske))

    ! NSE
    v1 = merge(logy - logx,  0.0_dp, maske)
    v2 = merge(logx - xmean, 0.0_dp, maske)
    LNNSE_2d_dp = 1.0_dp - sum(v1*v1, mask=maske) / sum(v2*v2, mask=maske)

  END FUNCTION LNNSE_2d_dp

  ! ------------------------------------------------------------------

  FUNCTION LNNSE_3d_sp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: mask
    REAL(sp)                                            :: LNNSE_3d_sp

    INTEGER(i4)                                         :: n
    INTEGER(i4), DIMENSION(size(shape(x)))              :: shapemask
    REAL(sp)                                            :: xmean
    REAL(sp), DIMENSION(size(x,dim=1),size(x,dim=2),size(x,dim=3)) :: logx, logy, v1, v2
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2),size(x,dim=3)) :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'LNNSE_3d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif

    ! mask all negative and zero entries
    where (x .lt. tiny(1.0_sp) .or. y .lt. tiny(1.0_sp))
       maske = .false.
    end where
    n = count(maske)
    if (n .LE. 1_i4) stop 'LNNSE_3d_sp: number of arguments must be at least 2'

    ! logarithms
    logx = 0.0_sp
    logy = 0.0_sp
    where (maske)
       logx = log(x)
       logy = log(y)
    end where

    ! mean of x
    xmean = average(pack(logx,maske))

    ! NSE
    v1 = merge(logy - logx,  0.0_sp, maske)
    v2 = merge(logx - xmean, 0.0_sp, maske)
    LNNSE_3d_sp = 1.0_sp - sum(v1*v1, mask=maske) / sum(v2*v2, mask=maske)

  END FUNCTION LNNSE_3d_sp

  ! ------------------------------------------------------------------

  FUNCTION LNNSE_3d_dp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: mask
    REAL(dp)                                            :: LNNSE_3d_dp

    INTEGER(i4)                                         :: n
    INTEGER(i4), DIMENSION(size(shape(x)))              :: shapemask
    REAL(dp)                                            :: xmean
    REAL(dp), DIMENSION(size(x,dim=1),size(x,dim=2),size(x,dim=3)) :: logx, logy, v1, v2
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2),size(x,dim=3)) :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'LNNSE_3d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif

    ! mask all negative and zero entries
    where (x .lt. tiny(1.0_dp) .or. y .lt. tiny(1.0_dp))
       maske = .false.
    end where
    n = count(maske)
    if (n .LE. 1_i4) stop 'LNNSE_3d_dp: number of arguments must be at least 2'

    ! logarithms
    logx = 0.0_dp
    logy = 0.0_dp
    where (maske)
       logx = log(x)
       logy = log(y)
    end where

    ! mean of x
    xmean = average(pack(logx,maske))

    ! NSE
    v1 = merge(logy - logx,  0.0_dp, maske)
    v2 = merge(logx - xmean, 0.0_dp, maske)
    LNNSE_3d_dp = 1.0_dp - sum(v1*v1, mask=maske) / sum(v2*v2, mask=maske)

  END FUNCTION LNNSE_3d_dp

  ! ------------------------------------------------------------------

  FUNCTION MAE_1d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                          :: MAE_1d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_1d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1)
    endif
    if (n .LE. 1_i4) stop 'MAE_1d_sp: number of arguments must be at least 2'
    !
    MAE_1d_sp = SAE_1d_sp(x,y,mask=maske) / real(n, sp)

  END FUNCTION MAE_1d_sp

  FUNCTION MAE_1d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                          :: MAE_1d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_1d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1)
    endif
    if (n .LE. 1_i4) stop 'MAE_1d_dp: number of arguments must be at least 2'
    !
    MAE_1d_dp = SAE_1d_dp(x,y,mask=maske) / real(n, dp)

  END FUNCTION MAE_1d_dp

  FUNCTION MAE_2d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: MAE_2d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_2d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'MAE_2d_sp: number of arguments must be at least 2'
    !
    MAE_2d_sp = SAE_2d_sp(x,y,mask=maske) / real(n, sp)

  END FUNCTION MAE_2d_sp

  FUNCTION MAE_2d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: MAE_2d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_2d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'MAE_2d_dp: number of arguments must be at least 2'
    !
    MAE_2d_dp = SAE_2d_dp(x,y,mask=maske) / real(n, dp)

  END FUNCTION MAE_2d_dp

  FUNCTION MAE_3d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: MAE_3d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_3d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'MAE_3d_sp: number of arguments must be at least 2'
    !
    MAE_3d_sp = SAE_3d_sp(x,y,mask=maske) / real(n, sp)

  END FUNCTION MAE_3d_sp

  FUNCTION MAE_3d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: MAE_3d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_3d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'MAE_3d_dp: number of arguments must be at least 2'
    !
    MAE_3d_dp = SAE_3d_dp(x,y,mask=maske) / real(n, dp)

  END FUNCTION MAE_3d_dp

  ! ------------------------------------------------------------------

  FUNCTION MAE_PROB_ONE_1d_dp(x, y, mask, cdfx, h)

    USE mo_moment, ONLY: average
    USE mo_kernel, ONLY: kernel_cumdensity

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp), DIMENSION(:), OPTIONAL, INTENT(IN)  :: cdfx
    REAL(dp),               OPTIONAL, INTENT(IN)  :: h
    REAL(dp)                                      :: MAE_PROB_ONE_1d_dp

    INTEGER(i4)                                   :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )       :: shapemask
    LOGICAL,     DIMENSION(size(x))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if

    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MAE_PROB_1d_dp: shapes of inputs(x,y) or mask are not matching'

    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'MAE_PROB_1d_dp: number of arguments must be at least 2'

    if (present(cdfx)) then
       if (present(h)) then
          MAE_PROB_ONE_1d_dp = average(abs(cdfx - &
               kernel_cumdensity(x, xout=y, h=h, mask=maske, romberg=.true., epsint=1.0e-4_dp)))
       else
          MAE_PROB_ONE_1d_dp = average(abs(cdfx - &
               kernel_cumdensity(x, xout=y, silverman=.true., mask=maske, romberg=.true., epsint=1.0e-4_dp)))
       endif
    else
       if (present(h)) then
          MAE_PROB_ONE_1d_dp = average(abs(kernel_cumdensity(x, xout=x, h=h, mask=maske, romberg=.true., epsint=1.0e-4_dp) &
               - kernel_cumdensity(x, xout=y, h=h, mask=maske, romberg=.true., epsint=1.0e-4_dp)))
       else
          MAE_PROB_ONE_1d_dp = average(abs( &
               kernel_cumdensity(x, xout=x, silverman=.true., mask=maske, romberg=.true., epsint=1.0e-4_dp) &
               - kernel_cumdensity(x, xout=y, silverman=.true., mask=maske, romberg=.true., epsint=1.0e-4_dp)))
       endif
    endif
    MAE_PROB_ONE_1d_dp = MAE_PROB_ONE_1d_dp * 100._dp ! unit is [%]

  END FUNCTION MAE_PROB_ONE_1d_dp

  ! ------------------------------------------------------------------

  ! FUNCTION MAE_PROB_TWO_1d_dp(x, y, mask)

  !   USE mo_moment, ONLY: average
  !   USE mo_empcdf, ONLY: empcdf

  !   IMPLICIT NONE

  !   REAL(dp), DIMENSION(:),           INTENT(IN)  :: x, y
  !   LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
  !   REAL(dp)                                      :: MAE_PROB_TWO_1d_dp

  !   INTEGER(i4)                                   :: n
  !   INTEGER(i4), DIMENSION(size(shape(x)) )       :: shapemask
  !   LOGICAL,  DIMENSION(size(x))                  :: maske

  !   if (present(mask)) then
  !      shapemask = shape(mask)
  !   else
  !      shapemask = shape(x)
  !   end if
  !   !
  !   if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
  !        stop 'MAE_PROB_1d_dp: shapes of inputs(x,y) or mask are not matching'
  !   !
  !   if (present(mask)) then
  !      maske = mask
  !      n = count(maske)
  !   else
  !      maske = .true.
  !      n = size(x)
  !   endif
  !   if (n .LE. 1_i4) stop 'MAE_PROB_1d_dp: number of arguments must be at least 2'
  !   !
  !   MAE_PROB_TWO_1d_dp = average(abs(EMPCDF(y, maske=maske, print_info=.False.) - EMPCDF(x, maske=maske, print_info=.False.)))

  ! END FUNCTION MAE_PROB_TWO_1d_dp

  ! ------------------------------------------------------------------

  FUNCTION MSE_1d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                          :: MSE_1d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MSE_1d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1)
    endif
    if (n .LE. 1_i4) stop 'MSE_1d_sp: number of arguments must be at least 2'
    !
    MSE_1d_sp = SSE_1d_sp(x,y,mask=maske) / real(n, sp)

  END FUNCTION MSE_1d_sp

  FUNCTION MSE_1d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                          :: MSE_1d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MSE_1d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1)
    endif
    if (n .LE. 1_i4) stop 'MSE_1d_dp: number of arguments must be at least 2'
    !
    MSE_1d_dp = SSE_1d_dp(x,y,mask=maske) / real(n, dp)

  END FUNCTION MSE_1d_dp

  FUNCTION MSE_2d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: MSE_2d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MSE_2d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'MSE_2d_sp: number of arguments must be at least 2'
    !
    MSE_2d_sp = SSE_2d_sp(x,y,mask=maske) / real(n, sp)

  END FUNCTION MSE_2d_sp

  FUNCTION MSE_2d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: MSE_2d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MSE_2d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'MSE_2d_dp: number of arguments must be at least 2'
    !
    MSE_2d_dp = SSE_2d_dp(x,y,mask=maske) / real(n, dp)

  END FUNCTION MSE_2d_dp

  FUNCTION MSE_3d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: MSE_3d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MSE_3d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'MSE_3d_sp: number of arguments must be at least 2'
    !
    MSE_3d_sp = SSE_3d_sp(x,y,mask=maske) / real(n, sp)

  END FUNCTION MSE_3d_sp

  FUNCTION MSE_3d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: MSE_3d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'MSE_3d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'MSE_3d_dp: number of arguments must be at least 2'
    !
    MSE_3d_dp = SSE_3d_dp(x,y,mask=maske) / real(n, dp)

  END FUNCTION MSE_3d_dp

  ! ------------------------------------------------------------------

  FUNCTION NSE_1d_sp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                          :: NSE_1d_sp

    INTEGER(i4)                            :: n
    INTEGER(i4), DIMENSION(size(shape(x))) :: shapemask
    REAL(sp)                               :: xmean
    REAL(sp), DIMENSION(size(x))           :: v1, v2
    LOGICAL,  DIMENSION(size(x))           :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'NSE_1d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'NSE_1d_sp: number of arguments must be at least 2'
    ! mean of x
    xmean = average(x, mask=maske)
    !
    v1 = merge(y - x    , 0.0_sp, maske)
    v2 = merge(x - xmean, 0.0_sp, maske)
    !
    NSE_1d_sp = 1.0_sp - dot_product(v1,v1) / dot_product(v2,v2)

  END FUNCTION NSE_1d_sp

  FUNCTION NSE_1d_dp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                          :: NSE_1d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(dp)                                          :: xmean
    REAL(dp), DIMENSION(size(x))                      :: v1, v2
    LOGICAL,  DIMENSION(size(x))                      :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'NSE_1d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'NSE_1d_dp: number of arguments must be at least 2'
    ! mean of x
    xmean = average(x, mask=maske)
    !
    v1 = merge(y - x    , 0.0_dp, maske)
    v2 = merge(x - xmean, 0.0_dp, maske)
    !
    NSE_1d_dp = 1.0_dp - dot_product(v1,v1) / dot_product(v2,v2)

  END FUNCTION NSE_1d_dp

  FUNCTION NSE_2d_sp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: NSE_2d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(sp)                                          :: xmean
    LOGICAL, DIMENSION(size(x, dim=1), size(x, dim=2)):: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'NSE_2d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2)
    endif
    !
    if (n .LE. 1_i4) stop 'NSE_2d_sp: number of arguments must be at least 2'
    ! mean of x
    xmean = average(reshape(x(:,:), (/size(x, dim=1)*size(x, dim=2)/)), &
         mask=reshape(maske(:,:), (/size(x, dim=1)*size(x, dim=2)/)))
    !
    NSE_2d_sp = 1.0_sp - sum((y-x)*(y-x), mask=maske) / sum((x-xmean)*(x-xmean), mask=maske)
    !
  END FUNCTION NSE_2d_sp

  FUNCTION NSE_2d_dp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: NSE_2d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(dp)                                          :: xmean
    LOGICAL, DIMENSION(size(x, dim=1), size(x, dim=2)):: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'NSE_2d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2)
    endif
    !
    if (n .LE. 1_i4) stop 'NSE_2d_dp: number of arguments must be at least 2'
    ! mean of x
    xmean = average(reshape(x(:,:), (/size(x, dim=1)*size(x, dim=2)/)), &
         mask=reshape(maske(:,:), (/size(x, dim=1)*size(x, dim=2)/)))
    !
    NSE_2d_dp = 1.0_dp - sum((y-x)*(y-x), mask=maske) / sum((x-xmean)*(x-xmean), mask=maske)
    !
  END FUNCTION NSE_2d_dp

  FUNCTION NSE_3d_sp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: NSE_3d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(sp)                                          :: xmean
    LOGICAL,  DIMENSION(size(x, dim=1), &
         size(x, dim=2), size(x, dim=3))    :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'NSE_3d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2) * size(x, dim=3)
    endif
    !
    if (n .LE. 1_i4) stop 'NSE_3d_sp: number of arguments must be at least 2'
    ! mean of x
    xmean = average(reshape(x(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), &
         mask=reshape(maske(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)))
    !
    NSE_3d_sp = 1.0_sp - sum((y-x)*(y-x), mask=maske) / sum((x-xmean)*(x-xmean), mask=maske)
    !
  END FUNCTION NSE_3d_sp

  FUNCTION NSE_3d_dp(x, y, mask)

    USE mo_moment, ONLY: average

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: NSE_3d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    REAL(dp)                                          :: xmean
    LOGICAL,  DIMENSION(size(x, dim=1), &
         size(x, dim=2), size(x, dim=3))    :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'NSE_3d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n     = count(maske)
    else
       maske = .true.
       n     = size(x, dim=1) * size(x, dim=2) * size(x, dim=3)
    endif
    !
    if (n .LE. 1_i4) stop 'NSE_3d_dp: number of arguments must be at least 2'
    ! Average of x
    xmean = average(reshape(x(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)), &
         mask=reshape(maske(:,:,:), (/size(x, dim=1)*size(x, dim=2)*size(x, dim=3)/)))
    !
    NSE_3d_dp = 1.0_dp - sum((y-x)*(y-x), mask=maske) / sum((x-xmean)*(x-xmean), mask=maske)
    !
  END FUNCTION NSE_3d_dp


  ! ------------------------------------------------------------------

  FUNCTION SAE_1d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                          :: SAE_1d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x))                      :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SAE_1d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SAE_1d_sp: number of arguments must be at least 2'
    !
    SAE_1d_sp = sum(abs(y - x) ,mask = maske)

  END FUNCTION SAE_1d_sp

  FUNCTION SAE_1d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                          :: SAE_1d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x))                      :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SAE_1d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SAE_1d_dp: number of arguments must be at least 2'
    !
    SAE_1d_dp = sum(abs(y - x) ,mask = maske)

  END FUNCTION SAE_1d_dp

  FUNCTION SAE_2d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: SAE_2d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SAE_2d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'SAE_2d_sp: number of arguments must be at least 2'
    !
    SAE_2d_sp = SAE_1d_sp(reshape(x, (/size(x, dim=1) * size(x, dim=2)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2)/)) )

  END FUNCTION SAE_2d_sp

  FUNCTION SAE_2d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: SAE_2d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SAE_2d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'SAE_2d_dp: number of arguments must be at least 2'
    !
    SAE_2d_dp = SAE_1d_dp(reshape(x, (/size(x, dim=1) * size(x, dim=2)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2)/)) )

  END FUNCTION SAE_2d_dp

  FUNCTION SAE_3d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: SAE_3d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SAE_3d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'SAE_3d_sp: number of arguments must be at least 2'
    !
    SAE_3d_sp = SAE_1d_sp(reshape(x, (/size(x, dim=1) * size(x, dim=2) * size(x, dim=3)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2) * size(x, dim=3)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2) &
         * size(maske, dim=3)/)) )

  END FUNCTION SAE_3d_sp

  FUNCTION SAE_3d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: SAE_3d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SAE_3d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'SAE_3d_dp: number of arguments must be at least 2'
    !
    SAE_3d_dp = SAE_1d_dp(reshape(x, (/size(x, dim=1) * size(x, dim=2) * size(x, dim=3)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2) * size(x, dim=3)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2) &
         * size(maske, dim=3)/)) )

  END FUNCTION SAE_3d_dp

  ! ------------------------------------------------------------------

  FUNCTION SSE_1d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                      :: SSE_1d_sp

    INTEGER(i4)                                   :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )       :: shapemask
    LOGICAL,  DIMENSION(size(x))                  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    !
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SSE_1d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SSE_1d_sp: number of arguments must be at least 2'
    !
    SSE_1d_sp = sum((y - x)**2_i4 ,mask = maske)

  END FUNCTION SSE_1d_sp

  FUNCTION SSE_1d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                      :: SSE_1d_dp

    INTEGER(i4)                                   :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )       :: shapemask
    LOGICAL,  DIMENSION(size(x))                  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SSE_1d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SSE_1d_dp: number of arguments must be at least 2'
    !
    SSE_1d_dp = sum((y - x)**2_i4 ,mask = maske)

  END FUNCTION SSE_1d_dp

  FUNCTION SSE_2d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: SSE_2d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)))            :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1), size(x,dim=2)) :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SSE_2d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SSE_2d_sp: number of arguments must be at least 2'
    !
    SSE_2d_sp = SSE_1d_sp(reshape(x, (/size(x, dim=1) * size(x, dim=2)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2)/)) )

  END FUNCTION SSE_2d_sp

  FUNCTION SSE_2d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: SSE_2d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)))            :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1), size(x,dim=2)) :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SSE_2d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SSE_2d_dp: number of arguments must be at least 2'
    !
    SSE_2d_dp = SSE_1d_dp(reshape(x, (/size(x, dim=1) * size(x, dim=2)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2)/)) )

  END FUNCTION SSE_2d_dp

  FUNCTION SSE_3d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: SSE_3d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)))            :: shapemask
    LOGICAL, DIMENSION(size(x,dim=1), size(x,dim=2),&
         size(x,dim=3))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SSE_3d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SSE_3d_sp: number of arguments must be at least 2'
    !
    SSE_3d_sp = SSE_1d_sp(reshape(x, (/size(x, dim=1) * size(x, dim=2) * size(x, dim=3)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2) * size(x, dim=3)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2)   &
         * size(maske, dim=3)/)))

  END FUNCTION SSE_3d_sp

  FUNCTION SSE_3d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: SSE_3d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)))            :: shapemask
    LOGICAL, DIMENSION(size(x,dim=1), size(x,dim=2),&
         size(x,dim=3))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'SSE_3d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x)
    endif
    if (n .LE. 1_i4) stop 'SSE_3d_dp: number of arguments must be at least 2'
    !
    SSE_3d_dp = SSE_1d_dp(reshape(x, (/size(x, dim=1) * size(x, dim=2) * size(x, dim=3)/)),                 &
         reshape(y, (/size(y, dim=1) * size(y, dim=2) * size(x, dim=3)/)),                 &
         mask=reshape(maske, (/size(maske, dim=1) * size(maske, dim=2)   &
         * size(maske, dim=3)/)))

  END FUNCTION SSE_3d_dp

  ! ------------------------------------------------------------------

  FUNCTION RMSE_1d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(sp)                                          :: RMSE_1d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'RMSE_1d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1)
    endif
    if (n .LE. 1_i4) stop 'RMSE_1d_sp: number of arguments must be at least 2'
    !
    RMSE_1d_sp = sqrt(MSE_1d_sp(x,y,mask=maske))

  END FUNCTION RMSE_1d_sp

  FUNCTION RMSE_1d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),           INTENT(IN)      :: x, y
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN)      :: mask
    REAL(dp)                                          :: RMSE_1d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x, dim=1))               :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'RMSE_1d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1)
    endif
    if (n .LE. 1_i4) stop 'RMSE_1d_dp: number of arguments must be at least 2'
    !
    RMSE_1d_dp = sqrt(MSE_1d_dp(x,y,mask=maske))

  END FUNCTION RMSE_1d_dp

  FUNCTION RMSE_2d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(sp)                                          :: RMSE_2d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'RMSE_2d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'RMSE_2d_sp: number of arguments must be at least 2'
    !
    RMSE_2d_sp = sqrt(MSE_2d_sp(x,y,mask=maske))

  END FUNCTION RMSE_2d_sp

  FUNCTION RMSE_2d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)    :: x, y
    LOGICAL,  DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: mask
    REAL(dp)                                          :: RMSE_2d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2))  :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'RMSE_2d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2)
    endif
    if (n .LE. 1_i4) stop 'RMSE_2d_dp: number of arguments must be at least 2'
    !
    RMSE_2d_dp = sqrt(MSE_2d_dp(x,y,mask=maske))

  END FUNCTION RMSE_2d_dp

  FUNCTION RMSE_3d_sp(x, y, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: RMSE_3d_sp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'RMSE_3d_sp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'RMSE_3d_sp: number of arguments must be at least 2'
    !
    RMSE_3d_sp = sqrt(MSE_3d_sp(x,y,mask=maske))

  END FUNCTION RMSE_3d_sp

  FUNCTION RMSE_3d_dp(x, y, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)  :: x, y
    LOGICAL,  DIMENSION(:,:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: RMSE_3d_dp

    INTEGER(i4)                                       :: n
    INTEGER(i4), DIMENSION(size(shape(x)) )           :: shapemask
    LOGICAL,  DIMENSION(size(x,dim=1),size(x,dim=2), &
         size(x,dim=3))                :: maske

    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask =  shape(x)
    end if
    if ( (any(shape(x) .NE. shape(y))) .OR. (any(shape(x) .NE. shapemask)) ) &
         stop 'RMSE_3d_dp: shapes of inputs(x,y) or mask are not matching'
    !
    if (present(mask)) then
       maske = mask
       n = count(maske)
    else
       maske = .true.
       n = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)
    endif
    if (n .LE. 1_i4) stop 'RMSE_3d_dp: number of arguments must be at least 2'
    !
    RMSE_3d_dp = sqrt(MSE_3d_dp(x,y,mask=maske))

  END FUNCTION RMSE_3d_dp

END MODULE mo_errormeasures

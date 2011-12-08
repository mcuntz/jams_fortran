MODULE mo_fit

  ! This module contains the numerical recipes routines to fit arbitrary functions to discrete data

  ! Usage:
  !   USE mo_fit, ONLY: fitfun, &                    ! Fit function: wrapper of svdfit & svdvar
  !                     fpoly, fpoly_sp, fpoly_dp, & ! Polynomial base functions
  !                     linfit, &                    ! Linear fit (without errors), Model I and Model II
  !                     polyfit, &                   ! Polynomial fitting routine using Lapack
  !                     svdfit, svdvar               ! fitting routines + errors

  ! Literature
  !   LinFit, SVDFit
  !       WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !           Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !           Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996
  !   For Model II (geometric mean regression) in LinFit
  !       Robert R. Sokal and F. James Rohlf, Biometry: the principle and practice of statistics
  !           in biological research, Freeman & Co., ISBN 0-7167-2411-1
  !   Polyfit
  !       http://rosettacode.org/wiki/Polynomial_regression#Fortran

  ! Attention: polyfit uses LAPACK

  ! Written March 2011, Matthias Cuntz - copied and adapted numerical recipes subroutines
  !                                    - fitfun for easy use with mask
  !                                    - linfit Model II: geometric mean regression
  ! Modified Matthias Cuntz,Nov 2011   - sp, dp
  !                                    - documentation

  USE mo_kind, ONLY: i4, sp, dp

  Implicit NONE

  PRIVATE

  PUBLIC :: fitfun                        ! Wrapper of svdfit and svdvar
  PUBLIC :: fpoly, fpoly_dp, fpoly_sp     ! Routine to fit polynomial with fitfun or svdfit
  PUBLIC :: linfit                        ! Fitting straight line (without error bars on input), Model I or Model II
  PUBLIC :: polyfit                       ! Fit a polynomial (with Lapack)
  PUBLIC :: svdfit                        ! Parameter fitting with singular value decomposition
  PUBLIC :: svdvar                        ! Variance of fitted parameters with svdfit

  ! Public interfaces
  INTERFACE fitfun
     MODULE PROCEDURE fitfun_sp, fitfun_dp
  END INTERFACE fitfun
  INTERFACE fpoly
     MODULE PROCEDURE fpoly_sp, fpoly_dp
  END INTERFACE fpoly
  INTERFACE linfit
     MODULE PROCEDURE linfit_sp, linfit_dp
  END INTERFACE linfit
  INTERFACE polyfit
     MODULE PROCEDURE polyfit_sp, polyfit_dp
  END INTERFACE polyfit
  INTERFACE svdfit
     MODULE PROCEDURE svdfit_sp, svdfit_dp
  END INTERFACE svdfit
  INTERFACE svdvar
     MODULE PROCEDURE svdvar_sp, svdvar_dp
  END INTERFACE svdvar

  ! Private interfaces, mostly from numerical recipes
  INTERFACE geop
     MODULE PROCEDURE geop_dp, geop_i4, geop_sp, geop_v_dp, geop_v_sp
  END INTERFACE geop
  INTERFACE outerprod
     MODULE PROCEDURE outerprod_dp, outerprod_sp
  END INTERFACE outerprod
  INTERFACE pythag
     MODULE PROCEDURE pythag_dp, pythag_sp
  END INTERFACE pythag
  INTERFACE svbksb
     MODULE PROCEDURE svbksb_dp, svbksb_sp
  END INTERFACE svbksb
  INTERFACE svdcmp
     MODULE PROCEDURE svdcmp_dp, svdcmp_sp
  END INTERFACE svdcmp
  INTERFACE vabs
     MODULE PROCEDURE vabs_dp, vabs_sp
  END INTERFACE vabs

  ! for geop
  INTEGER(i4), PARAMETER :: NPAR_GEOP=4, NPAR2_GEOP=2

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         fitfun

  !     PURPOSE
  !         Fits a given function func with the parameters a to the input data by using singular value decomposition.

  !         The function func must give a vector of the base functions for a given x.
  !         For example func for a polynomial should return the vector [1,x,x**2,x**3,x**4,...].

  !         From the numerical recipes documentation for svdfit:
  !         Given a set of N data points x, y with individual standard deviations sig, all arrays of length
  !         N, use chi^2 minimization to determine the M coefficients a of a function that depends linearly
  !         on a, y = sum_i=1^M(a x func_i(x)), with x meaning the dot-product.
  !         Here we solve the fitting equations using singular value decomposition of the NxM matrix, as in paragraph 2.6.
  !         The program returns the M basis functions evaluated at x=X in fit, the M fit parameters a,
  !         their sigma-errors siga, and chi^2 chisq.

  !         If an optinal mask is given, the function is only fit at those locations that correspond to true values in the mask.

  !     CALLING SEQUENCE
  !         call fitfun(xin, yin, sigin, a, func, mask=mask, fit=fit, chisq=chisq, siga=siga)
  
  !     INDENT(IN)
  !         real(sp/dp) :: xin(:)                1D-array with input x
  !         real(sp/dp) :: yin(:)                1D-array with input y
  !         real(sp/dp) :: sigin(:)              1D-array with input sigma on y
  !         real(sp/dp), dimension(M) :: FUNCTION func(x,n)    Function that outputs the n basis functions  

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(sp/dp), dimension(M) :: a       fitted M parameters

  !     INDENT(IN), OPTIONAL
  !         logical :: mask(:)        1D-array of logical values with size(vec).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         real(sp/dp), dimension(N) :: fit     Function evaluated at x=X
  !         real(sp/dp), dimension(M) :: siga    Error of fitted parameters
  !         real(sp/dp)               :: chisq   Minimum chi^2

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! fit polynomial with degree 5, then detrend original
  !         integer(i4), parameter    :: npt=100, npol=4 ! npol=degrees+1
  !         real(dp), dimension(npt)  :: x, y, sig
  !         real(dp), dimension(npol) :: a
  !         logical,  dimension(npt)  :: maske
  !         maske = .true.
  !         maske(10) = .false.
  !         ! fit
  !         call fitfun(x, y, sig, a, fpoly_dp, fit=fit, mask=maske)
  !         ! detrend
  !         do i=1, npt
  !            y(i) = y(i) - dot_product(a, fpoly(x(i),npol))
  !         end do
  !         ! or
  !         y = y - fit

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written, Matthias Cuntz, Mar 2011

  SUBROUTINE fitfun_dp(xin, yin, sigin, a, func, mask, fit, chisq, siga)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN)  :: xin
    REAL(dp), DIMENSION(:), INTENT(IN)  :: yin
    REAL(dp), DIMENSION(:), INTENT(IN)  :: sigin
    REAL(dp), DIMENSION(:), INTENT(OUT) :: a
    INTERFACE
       FUNCTION func(x,n)
         USE mo_kind, ONLY: i4, dp
         IMPLICIT NONE
         REAL(dp),    INTENT(IN)   :: x
         INTEGER(i4), INTENT(IN)   :: n
         REAL(dp),    DIMENSION(n) :: func
       END FUNCTION func
    END INTERFACE
    LOGICAL,  DIMENSION(:),         OPTIONAL, INTENT(IN)  :: mask
    REAL(dp), DIMENSION(size(xin)), OPTIONAL, INTENT(OUT) :: fit
    REAL(dp),                       OPTIONAL, INTENT(OUT) :: chisq
    REAL(dp), DIMENSION(size(a)),   OPTIONAL, INTENT(OUT) :: siga

    LOGICAL,  DIMENSION(size(xin))       :: maske
    REAL(dp), DIMENSION(size(a))         :: w
    REAL(dp), DIMENSION(size(a),size(a)) :: cvm, v
    REAL(dp), DIMENSION(:), ALLOCATABLE  :: x, y, sig
    REAL(dp), DIMENSION(:), ALLOCATABLE  :: out
    REAL(dp) :: ichisq
    INTEGER :: i, nmask, na
    
    if (size(xin) /= size(yin))   stop 'Error fitfun_dp: size(x) /= size(y)'
    if (size(sigin) /= size(yin)) stop 'Error fitfun_dp: size(sig) /= size(y)'
    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= size(xin)) stop 'Error fitfun_dp: size(mask) /= size(x)'
       maske = mask
    endif
    nmask = count(maske)
    allocate(x(nmask))
    allocate(y(nmask))
    allocate(sig(nmask))
    x(1:nmask)   = pack(xin,maske)
    y(1:nmask)   = pack(yin,maske)
    sig(1:nmask) = pack(sigin,maske)

    call svdfit(x, y, sig, a, v, w, ichisq, func)

    if (present(chisq)) chisq = ichisq
    na = size(a)
    if (present(siga)) then
       call svdvar(v, w, cvm)
       do i=1, na
          siga(i) = sqrt(cvm(i,i))
       end do
    end if
    if (present(fit)) then
       allocate(out(nmask))
       do i=1, nmask
          out(i) = dot_product(a, real(func(real(x(i),dp),na),dp))
       end do
       fit = unpack(out,maske,0.0_dp)
       deallocate(out)
    endif

    deallocate(x)
    deallocate(y)
    deallocate(sig)

  END SUBROUTINE fitfun_dp


  SUBROUTINE fitfun_sp(xin, yin, sigin, a, func, mask, fit, chisq, siga)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN)  :: xin
    REAL(sp), DIMENSION(:), INTENT(IN)  :: yin
    REAL(sp), DIMENSION(:), INTENT(IN)  :: sigin
    REAL(sp), DIMENSION(:), INTENT(OUT) :: a
    INTERFACE
       FUNCTION func(x,n)
         USE mo_kind, ONLY: i4, sp
         IMPLICIT NONE
         REAL(sp),    INTENT(IN)   :: x
         INTEGER(i4), INTENT(IN)   :: n
         REAL(sp),    DIMENSION(n) :: func
       END FUNCTION func
    END INTERFACE
    LOGICAL,  DIMENSION(:),         OPTIONAL, INTENT(IN)  :: mask
    REAL(sp), DIMENSION(size(xin)), OPTIONAL, INTENT(OUT) :: fit
    REAL(sp),                       OPTIONAL, INTENT(OUT) :: chisq
    REAL(sp), DIMENSION(size(a)),   OPTIONAL, INTENT(OUT) :: siga

    LOGICAL,  DIMENSION(size(xin))       :: maske
    REAL(sp), DIMENSION(size(a))         :: w
    REAL(sp), DIMENSION(size(a),size(a)) :: cvm, v
    REAL(sp), DIMENSION(:), ALLOCATABLE  :: x, y, sig
    REAL(sp), DIMENSION(:), ALLOCATABLE  :: out
    REAL(sp) :: ichisq
    INTEGER :: i, nmask, na
    
    if (size(xin) /= size(yin))   stop 'Error fitfun_sp: size(x) /= size(y)'
    if (size(sigin) /= size(yin)) stop 'Error fitfun_sp: size(sig) /= size(y)'
    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= size(xin)) stop 'Error fitfun_sp: size(mask) /= size(x)'
       maske = mask
    endif
    nmask = count(maske)
    allocate(x(nmask))
    allocate(y(nmask))
    allocate(sig(nmask))
    x(1:nmask)   = pack(xin,maske)
    y(1:nmask)   = pack(yin,maske)
    sig(1:nmask) = pack(sigin,maske)

    call svdfit(x, y, sig, a, v, w, ichisq, func)

    if (present(chisq)) chisq = ichisq
    na = size(a)
    if (present(siga)) then
       call svdvar(v, w, cvm)
       do i=1, na
          siga(i) = sqrt(cvm(i,i))
       end do
    end if
    if (present(fit)) then
       allocate(out(nmask))
       do i=1, nmask
          out(i) = dot_product(a, real(func(real(x(i),sp),na),sp))
       end do
       fit = unpack(out,maske,0.0_sp)
       deallocate(out)
    endif

    deallocate(x)
    deallocate(y)
    deallocate(sig)

  END SUBROUTINE fitfun_sp

  ! ------------------------------------------------------------------

  !     NAME
  !         fpoly

  !     PURPOSE
  !         Fitting routine for a polynomial of degree n-1, i.e. with n coefficients.
  !         Returns vector with [1,x,x**2,x**3,...,x**(n-1)].
  !         To use with fitfun or svdfit.

  !     CALLING SEQUENCE
  !         vec = func(x,n)
  
  !     INDENT(IN)
  !         real(sp/dp) :: x          x
  !         real(sp/dp) :: n          n-1 powers of x

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(sp/dp), dimension(n) :: func       vec with [1,x,x**2,x**3,...,x**(n-1)]

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         -> see example at fitfun

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written, Matthias Cuntz, Mar 2011

  FUNCTION fpoly_dp(x,n)

    IMPLICIT NONE

    REAL(dp),    INTENT(IN)   :: x
    INTEGER(i4), INTENT(IN)   :: n
    REAL(dp),    DIMENSION(n) :: fpoly_dp

    fpoly_dp = geop(1.0_dp,x,n)

  END FUNCTION fpoly_dp


  FUNCTION fpoly_sp(x,n)

    IMPLICIT NONE

    REAL(sp),    INTENT(IN)   :: x
    INTEGER(i4), INTENT(IN)   :: n
    REAL(sp),    DIMENSION(n) :: fpoly_sp

    fpoly_sp = geop(1.0_sp,x,n)

  END FUNCTION fpoly_sp

  ! ------------------------------------------------------------------

  ! Return a geometrical progression as an array.
  FUNCTION geop_dp(first,factor,n)

    REAL(dp),    INTENT(IN)   :: first,factor
    INTEGER(i4), INTENT(IN)   :: n
    REAL(dp),    DIMENSION(n) :: geop_dp

    INTEGER(i4) :: k,k2
    REAL(dp)     :: temp

    if (n > 0) geop_dp(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_dp(k)=geop_dp(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_dp(k)=geop_dp(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_dp(k+1:min(k2,n))=temp*geop_dp(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if

  END FUNCTION geop_dp


  FUNCTION geop_i4(first,factor,n)

    INTEGER(i4), INTENT(IN)   :: first, factor, n
    INTEGER(i4), DIMENSION(n) :: geop_i4

    INTEGER(i4) :: k, k2, temp

    if (n > 0) geop_i4(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_i4(k)=geop_i4(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_i4(k)=geop_i4(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_i4(k+1:min(k2,n))=temp*geop_i4(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if

  END FUNCTION geop_i4


  FUNCTION geop_sp(first,factor,n)

    REAL(sp),    INTENT(IN)   :: first,factor
    INTEGER(i4), INTENT(IN)   :: n
    REAL(sp),    DIMENSION(n) :: geop_sp

    INTEGER(i4) :: k,k2
    REAL(sp)     :: temp

    if (n > 0) geop_sp(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_sp(k)=geop_sp(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_sp(k)=geop_sp(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_sp(k+1:min(k2,n))=temp*geop_sp(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if

  END FUNCTION geop_sp


  FUNCTION geop_v_dp(first,factor,n)

    REAL(dp),     DIMENSION(:), INTENT(IN) :: first, factor
    INTEGER(i4),               INTENT(IN) :: n
    REAL(dp),     DIMENSION(size(first),n) :: geop_v_dp

    INTEGER(i4) :: k,k2
    REAL(dp), DIMENSION(size(first)) :: temp

    if (n > 0) geop_v_dp(:,1)=first(:)
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_v_dp(:,k)=geop_v_dp(:,k-1)*factor(:)
       end do
    else
       do k=2,NPAR2_GEOP
          geop_v_dp(:,k)=geop_v_dp(:,k-1)*factor(:)
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_v_dp(:,k+1:min(k2,n))=geop_v_dp(:,1:min(k,n-k))*&
               spread(temp,2,size(geop_v_dp(:,1:min(k,n-k)),2))
          temp=temp*temp
          k=k2
       end do
    end if

  END FUNCTION geop_v_dp


  FUNCTION geop_v_sp(first,factor,n)

    REAL(sp),     DIMENSION(:), INTENT(IN) :: first, factor
    INTEGER(i4),               INTENT(IN) :: n
    REAL(sp),     DIMENSION(size(first),n) :: geop_v_sp

    INTEGER(i4) :: k,k2
    REAL(sp), DIMENSION(size(first)) :: temp

    if (n > 0) geop_v_sp(:,1)=first(:)
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_v_sp(:,k)=geop_v_sp(:,k-1)*factor(:)
       end do
    else
       do k=2,NPAR2_GEOP
          geop_v_sp(:,k)=geop_v_sp(:,k-1)*factor(:)
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_v_sp(:,k+1:min(k2,n))=geop_v_sp(:,1:min(k,n-k))*&
               spread(temp,2,size(geop_v_sp(:,1:min(k,n-k)),2))
          temp=temp*temp
          k=k2
       end do
    end if

  END FUNCTION geop_v_sp

  ! ------------------------------------------------------------------

  !     NAME
  !         linfit

  !     PURPOSE
  !         Fits a straight line to input data by minimizing chi^2.

  !         Given a set of data points x(1:ndata), y(1:ndata), fit them to a straight liney = a+bx 
  !         by minimizing chi2.
  !         Model I minimizes y vs. x while Model II takes the geometric mean of y vs. x and x vs. y.
  !         Returned is the fitted line at x.
  !         Optional returns are a, b and their respective probable uncertainties siga and sigb,
  !         and the chi-square chi2.

  !     CALLING SEQUENCE
  !         out = linfit(x, y, a=a, b=b, siga=siga, sigb=sigb, chi2=chi2, model2=model2)
  
  !     INDENT(IN)
  !         real(sp/dp) :: x(:)                1D-array with input x
  !         real(sp/dp) :: y(:)                1D-array with input y

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(sp/dp), dimension(size(x)) :: out     fitted values

  !     INDENT(IN), OPTIONAL
  !         logical :: model2                  If present, use geometric mean regression instead of ordinary least square

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         real(sp/dp), dimension(M) :: a      intercept
  !         real(sp/dp), dimension(M) :: b      slope
  !         real(sp/dp), dimension(M) :: siga   error on intercept
  !         real(sp/dp), dimension(M) :: sigb   error on slope
  !         real(sp/dp)               :: chisq  Minimum chi^2

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ytmp = linfit(x,y, a=inter, b=slope, model2=.true.)

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written, Matthias Cuntz, Mar 2011

  FUNCTION linfit_dp(x, y, a, b, siga, sigb, chi2, model2)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN)  :: x, y
    REAL(dp), OPTIONAL,     INTENT(OUT) :: a, b, siga, sigb, chi2
    LOGICAL , OPTIONAL,     INTENT(IN)  :: model2
    REAL(dp), DIMENSION(size(x))        :: linfit_dp

    REAL(dp) :: sigdat, nx, sx, sxoss, sy, st2
    REAL(dp), DIMENSION(size(x)), TARGET :: t
    REAL(dp) :: aa, bb, cchi2
    LOGICAL :: mod2
    REAL(dp) :: mx, my, sx2, sy2, sxy, sxy2, syx2, r, ssigb

    if (size(x) /= size(y))   stop 'linfit_dp: size(x) /= size(y)'
    if (present(model2)) then
       mod2 = model2
    else
       mod2 = .false.
    endif
    
    if (mod2) then
       nx = real(size(x),dp)
       mx = sum(x)/nx
       my = sum(y)/nx
       sx2 = sum((x-mx)*(x-mx))
       sy2 = sum((y-my)*(y-my))
       sxy = sum((x-mx)*(y-my))
       r = sxy / sqrt(sx2) / sqrt(sy2)
       bb = sign(sqrt(sy2/sx2),sxy)
       aa = my - bb*mx
       linfit_dp(:) = aa + bb*x(:)
       if (present(a)) a = aa
       if (present(b)) b = bb       
       if (present(chi2)) then
          t(:) = y(:)-linfit_dp(:)
          chi2 = dot_product(t,t)
       endif
       if (present(siga) .or. present(sigb)) then
          syx2 = (sy2 - sxy*sxy/sx2) / (nx-2.0_dp)
          sxy2 = (sx2 - sxy*sxy/sy2) / (nx-2.0_dp)
          ssigb = sqrt(syx2/sx2)
          if (present(sigb)) sigb = ssigb
          if (present(siga)) then
             siga = sqrt(syx2*(1.0_dp/nX+mx*mx/sx2))
             ! Add Extra Term for Error in xmean which is not in Sokal & Rohlf.
             ! They take the error estimate of the chi-squared error for a.
             siga = sqrt(siga*siga + bb*bb*sxy2/nx)
          endif
       endif
    else
       nx = real(size(x),dp)
       sx = sum(x)
       sy = sum(y)
       sxoss = sx/nx
       t(:) = x(:)-sxoss
       bb = dot_product(t,y)
       st2 = dot_product(t,t)
       bb = bb/st2
       aa = (sy-sx*bb)/nx
       linfit_dp(:) = aa + bb*x(:)
       if (present(a)) a = aa
       if (present(b)) b = bb
       if (present(chi2) .or. present(siga) .or. present(sigb)) then
          t(:) = y(:)-linfit_dp(:)
          cchi2 = dot_product(t,t)
          if (present(chi2)) chi2 = cchi2
          if (present(siga) .or. present(sigb)) then
             sigdat = sqrt(cchi2/(size(x)-2))
             if (present(siga)) then
                siga = sqrt((1.0_dp+sx*sx/(nx*st2))/nx)
                siga = siga*sigdat
             endif
             if (present(sigb)) then
                sigb = sqrt(1.0_dp/st2)
                sigb = sigb*sigdat
             endif
          endif
       endif
    endif

  END FUNCTION linfit_dp


  FUNCTION linfit_sp(x, y, a, b, siga, sigb, chi2, model2)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN)  :: x, y
    REAL(sp), OPTIONAL,     INTENT(OUT) :: a, b, siga, sigb, chi2
    LOGICAL , OPTIONAL,     INTENT(IN)  :: model2
    REAL(sp), DIMENSION(size(x))        :: linfit_sp

    REAL(sp) :: sigdat, nx, sx, sxoss, sy, st2
    REAL(sp), DIMENSION(size(x)), TARGET :: t
    REAL(sp) :: aa, bb, cchi2
    LOGICAL :: mod2
    REAL(sp) :: mx, my, sx2, sy2, sxy, sxy2, syx2, r, ssigb

    if (size(x) /= size(y))   stop 'linfit_sp: size(x) /= size(y)'
    if (present(model2)) then
       mod2 = model2
    else
       mod2 = .false.
    endif
    
    if (mod2) then
       nx = real(size(x),sp)
       mx = sum(x)/nx
       my = sum(y)/nx
       sx2 = sum((x-mx)*(x-mx))
       sy2 = sum((y-my)*(y-my))
       sxy = sum((x-mx)*(y-my))
       r = sxy / sqrt(sx2) / sqrt(sy2)
       bb = sign(sqrt(sy2/sx2),sxy)
       aa = my - bb*mx
       linfit_sp(:) = aa + bb*x(:)
       if (present(a)) a = aa
       if (present(b)) b = bb       
       if (present(chi2)) then
          t(:) = y(:)-linfit_sp(:)
          chi2 = dot_product(t,t)
       endif
       if (present(siga) .or. present(sigb)) then
          syx2 = (sy2 - sxy*sxy/sx2) / (nx-2.0_sp)
          sxy2 = (sx2 - sxy*sxy/sy2) / (nx-2.0_sp)
          ssigb = sqrt(syx2/sx2)
          if (present(sigb)) sigb = ssigb
          if (present(siga)) then
             siga = sqrt(syx2*(1.0_sp/nX+mx*mx/sx2))
             ! Add Extra Term for Error in xmean which is not in Sokal & Rohlf.
             ! They take the error estimate of the chi-squared error for a.
             siga = sqrt(siga*siga + bb*bb*sxy2/nx)
          endif
       endif
    else
       nx = real(size(x),sp)
       sx = sum(x)
       sy = sum(y)
       sxoss = sx/nx
       t(:) = x(:)-sxoss
       bb = dot_product(t,y)
       st2 = dot_product(t,t)
       bb = bb/st2
       aa = (sy-sx*bb)/nx
       linfit_sp(:) = aa + bb*x(:)
       if (present(a)) a = aa
       if (present(b)) b = bb
       if (present(chi2) .or. present(siga) .or. present(sigb)) then
          t(:) = y(:)-linfit_sp(:)
          cchi2 = dot_product(t,t)
          if (present(chi2)) chi2 = cchi2
          if (present(siga) .or. present(sigb)) then
             sigdat = sqrt(cchi2/(size(x)-2))
             if (present(siga)) then
                siga = sqrt((1.0_sp+sx*sx/(nx*st2))/nx)
                siga = siga*sigdat
             endif
             if (present(sigb)) then
                sigb = sqrt(1.0_sp/st2)
                sigb = sigb*sigdat
             endif
          endif
       endif
    endif

  END FUNCTION linfit_sp

  ! ------------------------------------------------------------------

  ! a*b with size(a)/=size(b)
  FUNCTION outerprod_dp(a,b)

    REAL(dp), DIMENSION(:), INTENT(IN)   :: a, b
    REAL(dp), DIMENSION(size(a),size(b)) :: outerprod_dp

    outerprod_dp = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))

  END FUNCTION outerprod_dp


  FUNCTION outerprod_sp(a,b)

    REAL(sp), DIMENSION(:), INTENT(IN)   :: a, b
    REAL(sp), DIMENSION(size(a),size(b)) :: outerprod_sp

    outerprod_sp = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))

  END FUNCTION outerprod_sp

  ! ------------------------------------------------------------------

  ! sqrt(a**2+b**2)
  FUNCTION pythag_dp(a,b)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: a,b
    REAL(dp)             :: pythag_dp

    REAL(dp) :: absa, absb

    absa = abs(a)
    absb = abs(b)
    if (absa > absb) then
       pythag_dp = absa*sqrt(1.0_dp+(absb/absa)**2)
    else
       if (absb == 0.0_dp) then
          pythag_dp = 0.0_dp
       else
          pythag_dp = absb*sqrt(1.0_dp+(absa/absb)**2)
       end if
    end if

  END FUNCTION pythag_dp


  FUNCTION pythag_sp(a,b)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: a,b
    REAL(sp)             :: pythag_sp

    REAL(sp) :: absa, absb

    absa = abs(a)
    absb = abs(b)
    if (absa > absb) then
       pythag_sp = absa*sqrt(1.0_sp+(absb/absa)**2)
    else
       if (absb == 0.0_sp) then
          pythag_sp = 0.0_sp
       else
          pythag_sp = absb*sqrt(1.0_sp+(absa/absb)**2)
       end if
    end if

  END FUNCTION pythag_sp

  ! ------------------------------------------------------------------

  !     NAME
  !         polyfit

  !     PURPOSE
  !         Find an approximating polynom of known degree for a given data.

  !         The function polyfit implements least squares approximation of a function
  !         defined in the points as specified by the arrays xi and yi.
  !         The basis dj is xj, j=0,1,..,N. The implementation is straightforward.
  !         First the plane matrix A is created. Aji=dj(xi). Then the linear problem
  !         AATc=Ay is solved. The result cj are the coefficients. Constraint_Error
  !         is propagated when dimensions of X and Y differ or else when the problem is ill-defined.

  !         If an optinal mask is given, the function is only fit at those locations that correspond to true values in the mask.

  !     CALLING SEQUENCE
  !         out = polyfit(x, y, d, mask=mask)
  
  !     INDENT(IN)
  !         real(sp/dp) :: x(:)                x
  !         real(sp/dp) :: y(:)                y
  !         integer     :: d                   d+1 parameters = polynom of degree d

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(sp/dp), dimension(d+1) :: polyfit       fitted d+1 parameters for x^0, x^1, ..., x^d

  !     INDENT(IN), OPTIONAL
  !         logical :: mask(:)        1D-array of logical values with size(x)
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         integer, parameter      :: degree = 2
  !         integer                 :: i
  !         real(dp), dimension(11) :: x = (/ (i,i=0,10) /)
  !         real(dp), dimension(11) :: y = (/ 1,   6,  17,  34, 57,  86, 121, 162, 209, 262, 321 /)
  !         real(dp), dimension(degree+1) :: a
  !         a = polyfit(x, y, degree)            
  !         write (*, '(F9.4)') a
  !         -> gives 1.0, 2.0, 3.0

  !     LITERATURE
  !         http://rosettacode.org/wiki/Polynomial_regression#Fortran

  !     HISTORY
  !         Written, Matthias Cuntz, Mar 2011

  FUNCTION polyfit_dp(vxin, vyin, d, mask)

    IMPLICIT NONE

    real(dp), dimension(:), intent(in)    :: vxin, vyin
    integer,                intent(in)    :: d
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    real(dp), dimension(d+1)              :: polyfit_dp

    real(dp), dimension(:,:), allocatable :: X
    real(dp), dimension(:,:), allocatable :: XT
    real(dp), dimension(:,:), allocatable :: XTX
    integer :: i
    integer :: n, lda, lwork
    integer :: info
    integer, dimension(:), allocatable :: ipiv
    real(dp), dimension(:), allocatable :: work
    real(dp), dimension(:), allocatable :: vx, vy
    LOGICAL,   DIMENSION(size(vxin)) :: maske
    integer :: nmask

    EXTERNAL :: DGETRF, DGETRI

    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= size(vxin)) stop 'Error polyfit_dp: size(mask) /= size(vxin)'
       maske = mask
    endif
    nmask = count(maske)
    allocate(vx(nmask))
    allocate(vy(nmask))
    vx(1:nmask) = pack(vxin,maske)
    vy(1:nmask) = pack(vyin,maske)

    n     = d+1
    lda   = n
    lwork = n
    allocate(ipiv(n))
    allocate(work(lwork))
    allocate(XT(n,size(vx)))
    allocate(X(size(vx),n))
    allocate(XTX(n,n))

    ! prepare the matrix
    X(1:nmask,1) = 1.0_dp
    if (n > 1) X(1:nmask,2) = vx
    do i=3, n
       X(1:nmask,i) = vx(1:nmask)**(i-1)
    end do

    XT  = transpose(X)
    XTX = matmul(XT,X)

    ! calls to LAPACK subs DGETRF and DGETRI
    call DGETRF(n, n, XTX, lda, ipiv, info)
    if ( info /= 0 ) stop 'Problem polyfit_dp 1'
    call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
    if ( info /= 0 ) stop 'Problem polyfit_dp 2'

    polyfit_dp = matmul(matmul(XTX,XT), vy)

    deallocate(ipiv)
    deallocate(work)
    deallocate(X)
    deallocate(XT)
    deallocate(XTX)

  END FUNCTION polyfit_dp


  FUNCTION polyfit_sp(vxin, vyin, d, mask)

    IMPLICIT NONE

    real(sp), dimension(:), intent(in)    :: vxin, vyin
    integer,                intent(in)    :: d
    LOGICAL,  DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    real(sp), dimension(d+1)              :: polyfit_sp

    real(sp), dimension(:,:), allocatable :: X
    real(sp), dimension(:,:), allocatable :: XT
    real(sp), dimension(:,:), allocatable :: XTX
    integer :: i
    integer :: n, lda, lwork
    integer :: info
    integer, dimension(:), allocatable :: ipiv
    real(sp), dimension(:), allocatable :: work
    real(sp), dimension(:), allocatable :: vx, vy
    LOGICAL,   DIMENSION(size(vxin)) :: maske
    integer :: nmask

    EXTERNAL :: SGETRF, SGETRI

    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= size(vxin)) stop 'Error polyfit_sp: size(mask) /= size(vxin)'
       maske = mask
    endif
    nmask = count(maske)
    allocate(vx(nmask))
    allocate(vy(nmask))
    vx(1:nmask) = pack(vxin,maske)
    vy(1:nmask) = pack(vyin,maske)

    n     = d+1
    lda   = n
    lwork = n
    allocate(ipiv(n))
    allocate(work(lwork))
    allocate(XT(n,size(vx)))
    allocate(X(size(vx),n))
    allocate(XTX(n,n))

    ! prepare the matrix
    X(1:nmask,1) = 1.0_sp
    if (n > 1) X(1:nmask,2) = vx
    do i=3, n
       X(1:nmask,i) = vx(1:nmask)**(i-1)
    end do

    XT  = transpose(X)
    XTX = matmul(XT,X)

    ! calls to LAPACK subs SGETRF and SGETRI
    call SGETRF(n, n, XTX, lda, ipiv, info)
    if ( info /= 0 ) stop 'Problem polyfit_sp 1'
    call SGETRI(n, XTX, lda, ipiv, work, lwork, info)
    if ( info /= 0 ) stop 'Problem polyfit_sp 2'

    polyfit_sp = matmul(matmul(XTX,XT), vy)

    deallocate(ipiv)
    deallocate(work)
    deallocate(X)
    deallocate(XT)
    deallocate(XTX)

  END FUNCTION polyfit_sp

  ! ------------------------------------------------------------------

  ! From Numerical Recipes in F77 book
  ! Solves A*X=B for a vector X, where A is specified by the arrays u, w, v as returned by 
  ! svdcmp. m and n are the logical dimensions of a, and will be equal for square matrices. mp 
  ! and np are the physical dimensions of a. b(1:m) is the input right-handside. x(1:n) is 
  ! the output solution vector. No input quantities are destroyed, so the routine may be called 
  ! sequentially with different bs.
  SUBROUTINE svbksb_dp(u,w,v,b,x)

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: u,v
    REAL(dp), DIMENSION(:), INTENT(IN) :: w,b
    REAL(dp), DIMENSION(:), INTENT(OUT) :: x
    INTEGER(i4) :: mdum,ndum
    REAL(dp), DIMENSION(size(x)) :: tmp
    REAL(dp), DIMENSION(size(w)) :: wtmp

    if (size(u,1) /= size(b))   stop 'svbksb_dp: size(u,1) /= size(b)'
    if (size(u,2) /= size(v,1)) stop 'svbksb_dp: size(u,2) /= size(v,1)'
    if (size(u,2) /= size(v,2)) stop 'svbksb_dp: size(u,2) /= size(v,2)'
    if (size(u,2) /= size(w))   stop 'svbksb_dp: size(u,2) /= size(w)'
    if (size(u,2) /= size(x))   stop 'svbksb_dp: size(u,2) /= size(x)'
    mdum = size(u,1)
    ndum = size(u,2)
    wtmp = w
    where (abs(w) < tiny(1.0_dp)) wtmp = 1.0_dp
    tmp = matmul(b,u) / wtmp
    where (abs(w) < tiny(1.0_dp)) tmp = 0.0_dp
    x = matmul(v,tmp)

  END SUBROUTINE svbksb_dp


  SUBROUTINE svbksb_sp(u,w,v,b,x)

    REAL(sp), DIMENSION(:,:), INTENT(IN) :: u,v
    REAL(sp), DIMENSION(:), INTENT(IN) :: w,b
    REAL(sp), DIMENSION(:), INTENT(OUT) :: x
    INTEGER(i4) :: mdum,ndum
    REAL(sp), DIMENSION(size(x)) :: tmp
    REAL(sp), DIMENSION(size(w)) :: wtmp

    if (size(u,1) /= size(b))   stop 'svbksb_sp: size(u,1) /= size(b)'
    if (size(u,2) /= size(v,1)) stop 'svbksb_sp: size(u,2) /= size(v,1)'
    if (size(u,2) /= size(v,2)) stop 'svbksb_sp: size(u,2) /= size(v,2)'
    if (size(u,2) /= size(w))   stop 'svbksb_sp: size(u,2) /= size(w)'
    if (size(u,2) /= size(x))   stop 'svbksb_sp: size(u,2) /= size(x)'
    mdum = size(u,1)
    ndum = size(u,2)
    wtmp = w
    where (abs(w) < tiny(1.0_sp)) wtmp = 1.0_sp
    tmp = matmul(b,u) / wtmp
    where (abs(w) < tiny(1.0_sp)) tmp = 0.0_sp
    x = matmul(v,tmp)

  END SUBROUTINE svbksb_sp

  ! ------------------------------------------------------------------

  ! From Numerical Recipes in F77 book
  ! Given a matrix a(1:m,1:n), with physical dimensions mp by np, this routine computes its 
  ! singular value decomposition, A=U*W*VT. The matrix U replaces a on output. The 
  ! diagonal matrix of singular values W is output as a vector w(1:n). The matrix V (not the 
  ! transpose VT) is output as v(1:n,1:n).
  SUBROUTINE svdcmp_dp(a,w,v)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: a
    REAL(dp), DIMENSION(:),   INTENT(OUT)   :: w
    REAL(dp), DIMENSION(:,:), INTENT(OUT)   :: v

    INTEGER(i4) :: i, its, j, k, l, m, n, nm
    REAL(dp)    :: anorm, c, f, g, h, s, scale, x, y, z
    REAL(dp), DIMENSION(size(a,1)) :: tempm
    REAL(dp), DIMENSION(size(a,2)) :: rv1, tempn

    m = size(a,1)
    if (size(a,2) /= size(v,1)) stop 'Error svdcmp_dp: size(a,2) /= size(v,1)'
    if (size(a,2) /= size(v,2)) stop 'Error svdcmp_dp: size(a,2) /= size(v,2)'
    if (size(a,2) /= size(w))   stop 'Error svdcmp_dp: size(a,2) /= size(w)'
    n = size(a,2)
    g     = 0.0_dp
    scale = 0.0_dp
    do i=1, n
       l = i+1
       rv1(i) = scale*g
       g     = 0.0_dp
       scale = 0.0_dp
       if (i <= m) then
          scale = sum(abs(a(i:m,i)))
          if (scale /= 0.0_dp) then
             a(i:m,i) = a(i:m,i)/scale
             s = dot_product(a(i:m,i),a(i:m,i))
             f = a(i,i)
             g = -sign(sqrt(s),f)
             h = f*g-s
             a(i,i) = f-g
             tempn(l:n) = matmul(a(i:m,i),a(i:m,l:n))/h
             a(i:m,l:n) = a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
             a(i:m,i)   = scale*a(i:m,i)
          end if
       end if
       w(i) = scale*g
       g     = 0.0_dp
       scale = 0.0_dp
       if ((i <= m) .and. (i /= n)) then
          scale = sum(abs(a(i,l:n)))
          if (scale /= 0.0_dp) then
             a(i,l:n) = a(i,l:n)/scale
             s = dot_product(a(i,l:n),a(i,l:n))
             f = a(i,l)
             g = -sign(sqrt(s),f)
             h = f*g-s
             a(i,l) = f-g
             rv1(l:n) = a(i,l:n)/h
             tempm(l:m) = matmul(a(l:m,l:n),a(i,l:n))
             a(l:m,l:n) = a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
             a(i,l:n)   = scale*a(i,l:n)
          end if
       end if
    end do
    anorm = maxval(abs(w)+abs(rv1))
    do i=n, 1, -1
       if (i < n) then
          if (g /= 0.0_dp) then
             v(l:n,i)   = (a(i,l:n)/a(i,l))/g
             tempn(l:n) = matmul(a(i,l:n),v(l:n,l:n))
             v(l:n,l:n) = v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
          end if
          v(i,l:n) = 0.0_dp
          v(l:n,i) = 0.0_dp
       end if
       v(i,i) = 1.0_dp
       g = rv1(i)
       l = i
    end do
    do i=min(m,n), 1, -1
       l = i+1
       g = w(i)
       a(i,l:n) = 0.0_dp
       if (g /= 0.0_dp) then
          g = 1.0_dp/g
          tempn(l:n) = (matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
          a(i:m,l:n) = a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
          a(i:m,i)   = a(i:m,i)*g
       else
          a(i:m,i) = 0.0_dp
       end if
       a(i,i) = a(i,i)+1.0_dp
    end do
    do k=n, 1, -1
       do its=1, 30
          do l=k, 1, -1
             nm = l-1
             if ((abs(rv1(l))+anorm) == anorm) exit
             if ((abs(w(nm))+anorm) == anorm) then
                c = 0.0_dp
                s = 1.0_dp
                do i=l, k
                   f = s*rv1(i)
                   rv1(i) = c*rv1(i)
                   if ((abs(f)+anorm) == anorm) exit
                   g = w(i)
                   h = pythag(f,g)
                   w(i) = h
                   h = 1.0_dp/h
                   c =  (g*h)
                   s = -(f*h)
                   tempm(1:m) = a(1:m,nm)
                   a(1:m,nm)  = a(1:m,nm)*c+a(1:m,i)*s
                   a(1:m,i)   = -tempm(1:m)*s+a(1:m,i)*c
                end do
                exit
             end if
          end do
          z = w(k)
          if (l == k) then
             if (z < 0.0_dp) then
                w(k) = -z
                v(1:n,k) = -v(1:n,k)
             end if
             exit
          end if
          if (its == 30) stop 'svdcmp_dp: no convergence in svdcmp_dp'
          x = w(l)
          nm = k-1
          y = w(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_dp*h*y)
          g = pythag(f,1.0_dp)
          f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c = 1.0_dp
          s = 1.0_dp
          do j=l, nm
             i = j+1
             g = rv1(i)
             y = w(i)
             h = s*g
             g = c*g
             z = pythag(f,h)
             rv1(j) = z
             c = f/z
             s = h/z
             f =  (x*c)+(g*s)
             g = -(x*s)+(g*c)
             h = y*s
             y = y*c
             tempn(1:n) = v(1:n,j)
             v(1:n,j)   = v(1:n,j)*c+v(1:n,i)*s
             v(1:n,i)   = -tempn(1:n)*s+v(1:n,i)*c
             z = pythag(f,h)
             w(j) = z
             if (z /= 0.0_dp) then
                z = 1.0_dp/z
                c = f*z
                s = h*z
             end if
             f =  (c*g)+(s*y)
             x = -(s*g)+(c*y)
             tempm(1:m) = a(1:m,j)
             a(1:m,j)   = a(1:m,j)*c+a(1:m,i)*s
             a(1:m,i)   = -tempm(1:m)*s+a(1:m,i)*c
          end do
          rv1(l) = 0.0_dp
          rv1(k) = f
          w(k)   = x
       end do
    end do

  END SUBROUTINE svdcmp_dp


  SUBROUTINE svdcmp_sp(a,w,v)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(INOUT) :: a
    REAL(sp), DIMENSION(:),   INTENT(OUT)   :: w
    REAL(sp), DIMENSION(:,:), INTENT(OUT)   :: v

    INTEGER(i4) :: i, its, j, k, l, m, n, nm
    REAL(sp)    :: anorm, c, f, g, h, s, scale, x, y, z
    REAL(sp), DIMENSION(size(a,1)) :: tempm
    REAL(sp), DIMENSION(size(a,2)) :: rv1, tempn

    m = size(a,1)
    if (size(a,2) /= size(v,1)) stop 'Error svdcmp_sp: size(a,2) /= size(v,1)'
    if (size(a,2) /= size(v,2)) stop 'Error svdcmp_sp: size(a,2) /= size(v,2)'
    if (size(a,2) /= size(w))   stop 'Error svdcmp_sp: size(a,2) /= size(w)'
    n = size(a,2)
    g     = 0.0_sp
    scale = 0.0_sp
    do i=1, n
       l = i+1
       rv1(i) = scale*g
       g     = 0.0_sp
       scale = 0.0_sp
       if (i <= m) then
          scale = sum(abs(a(i:m,i)))
          if (scale /= 0.0_sp) then
             a(i:m,i) = a(i:m,i)/scale
             s = dot_product(a(i:m,i),a(i:m,i))
             f = a(i,i)
             g = -sign(sqrt(s),f)
             h = f*g-s
             a(i,i) = f-g
             tempn(l:n) = matmul(a(i:m,i),a(i:m,l:n))/h
             a(i:m,l:n) = a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
             a(i:m,i)   = scale*a(i:m,i)
          end if
       end if
       w(i) = scale*g
       g     = 0.0_sp
       scale = 0.0_sp
       if ((i <= m) .and. (i /= n)) then
          scale = sum(abs(a(i,l:n)))
          if (scale /= 0.0_sp) then
             a(i,l:n) = a(i,l:n)/scale
             s = dot_product(a(i,l:n),a(i,l:n))
             f = a(i,l)
             g = -sign(sqrt(s),f)
             h = f*g-s
             a(i,l) = f-g
             rv1(l:n) = a(i,l:n)/h
             tempm(l:m) = matmul(a(l:m,l:n),a(i,l:n))
             a(l:m,l:n) = a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
             a(i,l:n)   = scale*a(i,l:n)
          end if
       end if
    end do
    anorm = maxval(abs(w)+abs(rv1))
    do i=n, 1, -1
       if (i < n) then
          if (g /= 0.0_sp) then
             v(l:n,i)   = (a(i,l:n)/a(i,l))/g
             tempn(l:n) = matmul(a(i,l:n),v(l:n,l:n))
             v(l:n,l:n) = v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
          end if
          v(i,l:n) = 0.0_sp
          v(l:n,i) = 0.0_sp
       end if
       v(i,i) = 1.0_sp
       g = rv1(i)
       l = i
    end do
    do i=min(m,n), 1, -1
       l = i+1
       g = w(i)
       a(i,l:n) = 0.0_sp
       if (g /= 0.0_sp) then
          g = 1.0_sp/g
          tempn(l:n) = (matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
          a(i:m,l:n) = a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
          a(i:m,i)   = a(i:m,i)*g
       else
          a(i:m,i) = 0.0_sp
       end if
       a(i,i) = a(i,i)+1.0_sp
    end do
    do k=n, 1, -1
       do its=1, 30
          do l=k, 1, -1
             nm = l-1
             if ((abs(rv1(l))+anorm) == anorm) exit
             if ((abs(w(nm))+anorm) == anorm) then
                c = 0.0_sp
                s = 1.0_sp
                do i=l, k
                   f = s*rv1(i)
                   rv1(i) = c*rv1(i)
                   if ((abs(f)+anorm) == anorm) exit
                   g = w(i)
                   h = pythag(f,g)
                   w(i) = h
                   h = 1.0_sp/h
                   c =  (g*h)
                   s = -(f*h)
                   tempm(1:m) = a(1:m,nm)
                   a(1:m,nm)  = a(1:m,nm)*c+a(1:m,i)*s
                   a(1:m,i)   = -tempm(1:m)*s+a(1:m,i)*c
                end do
                exit
             end if
          end do
          z = w(k)
          if (l == k) then
             if (z < 0.0_sp) then
                w(k) = -z
                v(1:n,k) = -v(1:n,k)
             end if
             exit
          end if
          if (its == 30) stop 'svdcmp_sp: no convergence in svdcmp_sp'
          x = w(l)
          nm = k-1
          y = w(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_sp*h*y)
          g = pythag(f,1.0_sp)
          f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c = 1.0_sp
          s = 1.0_sp
          do j=l, nm
             i = j+1
             g = rv1(i)
             y = w(i)
             h = s*g
             g = c*g
             z = pythag(f,h)
             rv1(j) = z
             c = f/z
             s = h/z
             f =  (x*c)+(g*s)
             g = -(x*s)+(g*c)
             h = y*s
             y = y*c
             tempn(1:n) = v(1:n,j)
             v(1:n,j)   = v(1:n,j)*c+v(1:n,i)*s
             v(1:n,i)   = -tempn(1:n)*s+v(1:n,i)*c
             z = pythag(f,h)
             w(j) = z
             if (z /= 0.0_sp) then
                z = 1.0_sp/z
                c = f*z
                s = h*z
             end if
             f =  (c*g)+(s*y)
             x = -(s*g)+(c*y)
             tempm(1:m) = a(1:m,j)
             a(1:m,j)   = a(1:m,j)*c+a(1:m,i)*s
             a(1:m,i)   = -tempm(1:m)*s+a(1:m,i)*c
          end do
          rv1(l) = 0.0_sp
          rv1(k) = f
          w(k)   = x
       end do
    end do

  END SUBROUTINE svdcmp_sp

  ! ------------------------------------------------------------------

  !     NAME
  !         svdfit

  !     PURPOSE
  !         Fits a given function func with the parameters a to the input data by using singular value decomposition.

  !         The function func must give a vector of the base functions for a given x.
  !         For example func for a polynomial should return the vector [1,x,x**2,x**3,x**4,...].

  !         From the numerical recipes documentation:
  !         Given a set of datapoints x(1:ndata), y(1:ndata) with individual standard deviations
  !         sig(1:ndata), use chi2 minimization to determine the ma coefficients a of the fitting function
  !         y=sum(ai*afunci(x))). Here we solve the fitting equations using singular value decomposition
  !         of the n data by ma matrix, as in paragraph 2.6. Arrays u(1:mp,1:np), v(1:np,1:np),
  !         w(1:np) provide workspace on input; on output they define the singular value decomposition,
  !         and can be used to obtain the covariance matrix. mp,np are the physical dimensions
  !         of the matrices u,v,w, as indicated above. It is necessary that mp>=ndata, np>=ma. The
  !         program returns values for the ma fit parameters a, and chi2, chisq. The user supplies a
  !         subroutine funcs(x,afunc,ma) that returns them a basis functions evaluated at x=x
  !         in the array afunc.

  !     CALLING SEQUENCE
  !         call svdfit(x, y, sig, a, v, w, chisq, func)
  
  !     INDENT(IN)
  !         real(sp/dp) :: x(:)                1D-array with input x
  !         real(sp/dp) :: y(:)                1D-array with input y
  !         real(sp/dp) :: sig(:)              1D-array with input sigma on y
  !         real(sp/dp), dimension(M) :: FUNCTION func(x,n)    Function that outputs the n basis functions  

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(sp/dp), dimension(M)   :: a       fitted M parameters
  !         real(sp/dp)                 :: chisq   Minimum chi^2
  !         real(sp/dp), dimension(M,M) :: v       output for svdvar
  !         real(sp/dp), dimension(M)   :: w       output for svdvar

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! fit polynomial with degree 5
  !         integer(i4b), parameter :: npt=100, npol=5
  !         integer(i4b) :: i
  !         real(dp) :: chisq
  !         real(dp), dimension(npol)      :: a, w
  !         real(dp), dimension(npt)       :: x, y, sig
  !         real(dp), dimension(npol,npol) :: cvm, v
  !         call svdfit(x, y, sig, a, v, w, chisq, fpoly_dp)
  !         call svdvar(v, w, cvm)
  !         do i=1, npol
  !            write(*,'(1x,f12.6,a,f10.6)') a(i), '  +-', sqrt(cvm(i,i))
  !         end do

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written, Matthias Cuntz, Mar 2011

  SUBROUTINE svdfit_dp(ix, iy, isig, oa, ov, ow, ochisq, func)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),   INTENT(IN)  :: ix, iy, isig
    REAL(dp), DIMENSION(:),   INTENT(OUT) :: oa, ow
    REAL(dp), DIMENSION(:,:), INTENT(OUT) :: ov
    REAL(dp),                 INTENT(OUT) :: ochisq
    INTERFACE
       FUNCTION func(x,n)
         USE mo_kind, ONLY: i4, dp
         IMPLICIT NONE
         REAL(dp),    INTENT(IN)   :: x
         INTEGER(i4), INTENT(IN)   :: n
         REAL(dp),    DIMENSION(n) :: func
       END FUNCTION func
    END INTERFACE

    REAL(dp), DIMENSION(size(ix))   :: x
    REAL(dp), DIMENSION(size(iy))   :: y
    REAL(dp), DIMENSION(size(isig)) :: sig
    REAL(dp), DIMENSION(size(oa))   :: a
    REAL(dp), DIMENSION(size(ow))   :: w
    REAL(dp), DIMENSION(size(ov,1),size(ov,2)) :: v
    REAL(dp) :: chisq

    REAL(dp), PARAMETER :: TOL = 1.0e-5_dp
    INTEGER(i4) :: i, ma, n
    REAL(dp), DIMENSION(size(x))         :: b, sigi
    REAL(dp), DIMENSION(size(x),size(a)) :: u, usav

    if (size(ix) /= size(iy))   stop 'Error svdfit_dp: size(x) /= size(y)'
    if (size(iy) /= size(isig)) stop 'Error svdfit_dp: size(y) /= size(sig)'
    if (size(oa) /= size(ov,1)) stop 'Error svdfit_dp: size(a) /= size(v,1)'
    if (size(oa) /= size(ov,2)) stop 'Error svdfit_dp: size(a) /= size(v,2)'
    if (size(oa) /= size(ow))   stop 'Error svdfit_dp: size(a) /= size(w)'

    x   = real(ix,dp)
    y   = real(iy,dp)
    sig = real(isig,dp)

    n  = size(x)
    ma = size(a)
    sigi = 1.0_dp/sig
    b    = y*sigi
    do i=1, n
       usav(i,:) = func(x(i),ma)
    end do
    u    = usav*spread(sigi,dim=2,ncopies=ma)
    usav = u
    call svdcmp(u,w,v)
    where (w < TOL*maxval(w)) w = 0.0_dp
    call svbksb(u,w,v,b,a)
    chisq = vabs(matmul(usav,a)-b)**2

    oa = real(a,dp)
    ow = real(w,dp)
    ov = real(v,dp)
    ochisq = real(chisq,dp)

  END SUBROUTINE svdfit_dp


  SUBROUTINE svdfit_sp(ix, iy, isig, oa, ov, ow, ochisq, func)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),   INTENT(IN)  :: ix, iy, isig
    REAL(sp), DIMENSION(:),   INTENT(OUT) :: oa, ow
    REAL(sp), DIMENSION(:,:), INTENT(OUT) :: ov
    REAL(sp),                 INTENT(OUT) :: ochisq
    INTERFACE
       FUNCTION func(x,n)
         USE mo_kind, ONLY: i4, sp
         IMPLICIT NONE
         REAL(sp),    INTENT(IN)   :: x
         INTEGER(i4), INTENT(IN)   :: n
         REAL(sp),    DIMENSION(n) :: func
       END FUNCTION func
    END INTERFACE

    REAL(sp), DIMENSION(size(ix))   :: x
    REAL(sp), DIMENSION(size(iy))   :: y
    REAL(sp), DIMENSION(size(isig)) :: sig
    REAL(sp), DIMENSION(size(oa))   :: a
    REAL(sp), DIMENSION(size(ow))   :: w
    REAL(sp), DIMENSION(size(ov,1),size(ov,2)) :: v
    REAL(sp) :: chisq

    REAL(sp), PARAMETER :: TOL = 1.0e-5_sp
    INTEGER(i4) :: i, ma, n
    REAL(sp), DIMENSION(size(x))         :: b, sigi
    REAL(sp), DIMENSION(size(x),size(a)) :: u, usav

    if (size(ix) /= size(iy))   stop 'Error svdfit_sp: size(x) /= size(y)'
    if (size(iy) /= size(isig)) stop 'Error svdfit_sp: size(y) /= size(sig)'
    if (size(oa) /= size(ov,1)) stop 'Error svdfit_sp: size(a) /= size(v,1)'
    if (size(oa) /= size(ov,2)) stop 'Error svdfit_sp: size(a) /= size(v,2)'
    if (size(oa) /= size(ow))   stop 'Error svdfit_sp: size(a) /= size(w)'

    x   = real(ix,sp)
    y   = real(iy,sp)
    sig = real(isig,sp)

    n  = size(x)
    ma = size(a)
    sigi = 1.0_sp/sig
    b    = y*sigi
    do i=1, n
       usav(i,:) = func(x(i),ma)
    end do
    u    = usav*spread(sigi,dim=2,ncopies=ma)
    usav = u
    call svdcmp(u,w,v)
    where (w < TOL*maxval(w)) w = 0.0_sp
    call svbksb(u,w,v,b,a)
    chisq = vabs(matmul(usav,a)-b)**2

    oa = real(a,sp)
    ow = real(w,sp)
    ov = real(v,sp)
    ochisq = real(chisq,sp)

  END SUBROUTINE svdfit_sp

  ! ------------------------------------------------------------------

  !     NAME
  !         svdfit

  !     PURPOSE
  !         Calculates the variance/covariance for the singular value decomposition estimates.

  !         From the numerical recipes documentation:
  !         To evaluate the covariance matrix cvm of the fit for ma parameters obtained by svdfit,
  !         call this routine with matrices v,w as returned from svdfit. np,ncvm give the physical
  !         dimensions of v,w,cvm as indicated.

  !     CALLING SEQUENCE
  !         call svdvar(v, w, cvm)
  
  !     INDENT(IN)
  !         real(sp/dp), dimension(M,M) :: v       output from svdfit
  !         real(sp/dp), dimension(M)   :: w       output from svdfit

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(sp/dp), dimension(M,M)   :: cvm     covariance matrix for M fitted parameters

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         ! fit polynomial with degree 5
  !         integer(i4b), parameter :: npt=100, npol=5
  !         integer(i4b) :: i
  !         real(dp) :: chisq
  !         real(dp), dimension(npol)      :: a, w
  !         real(dp), dimension(npt)       :: x, y, sig
  !         real(dp), dimension(npol,npol) :: cvm, v
  !         call svdfit(x, y, sig, a, v, w, chisq, fpoly_dp)
  !         call svdvar(v, w, cvm)
  !         do i=1, npol
  !            write(*,'(1x,f12.6,a,f10.6)') a(i), '  +-', sqrt(cvm(i,i))
  !         end do

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written, Matthias Cuntz, Mar 2011

  SUBROUTINE svdvar_dp(v, w, cvm)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN)  :: v
    REAL(dp), DIMENSION(:),   INTENT(IN)  :: w
    REAL(dp), DIMENSION(:,:), INTENT(OUT) :: cvm

    INTEGER(i4) :: ma
    REAL(dp), DIMENSION(size(w)) :: wti
    REAL(dp), DIMENSION(size(w)) :: wtmp

    if (size(v,1) /= size(v,2))   stop 'Error svdvar_dp: size(v,1) /= size(v,2)'
    if (size(v,1) /= size(w))     stop 'Error svdvar_dp: size(v,1) /= size(w)'
    if (size(v,1) /= size(cvm,1)) stop 'Error svdvar_dp: size(v,1) /= size(cvm,1)'
    if (size(v,1) /= size(cvm,2)) stop 'Error svdvar_dp: size(v,1) /= size(cvm,2)'

    ma = size(v,1)
    wtmp = w
    where (abs(w) < tiny(1.0_dp)) wtmp = 1.0_dp
    wti = 1.0_dp/(wtmp*wtmp)
    where (abs(w) < tiny(1.0_dp)) wti = 0.0_dp
    cvm = v*spread(wti,dim=1,ncopies=ma)
    cvm = matmul(cvm,transpose(v))

  END SUBROUTINE svdvar_dp


  SUBROUTINE svdvar_sp(v, w, cvm)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN)  :: v
    REAL(sp), DIMENSION(:),   INTENT(IN)  :: w
    REAL(sp), DIMENSION(:,:), INTENT(OUT) :: cvm

    INTEGER(i4) :: ma
    REAL(sp), DIMENSION(size(w)) :: wti
    REAL(sp), DIMENSION(size(w)) :: wtmp

    if (size(v,1) /= size(v,2))   stop 'Error svdvar_sp: size(v,1) /= size(v,2)'
    if (size(v,1) /= size(w))     stop 'Error svdvar_sp: size(v,1) /= size(w)'
    if (size(v,1) /= size(cvm,1)) stop 'Error svdvar_sp: size(v,1) /= size(cvm,1)'
    if (size(v,1) /= size(cvm,2)) stop 'Error svdvar_sp: size(v,1) /= size(cvm,2)'

    ma = size(v,1)
    wtmp = w
    where (abs(w) < tiny(1.0_sp)) wtmp = 1.0_sp
    wti = 1.0_sp/(wtmp*wtmp)
    where (abs(w) < tiny(1.0_sp)) wti = 0.0_sp
    cvm = v*spread(wti,dim=1,ncopies=ma)
    cvm = matmul(cvm,transpose(v))

  END SUBROUTINE svdvar_sp

  ! ------------------------------------------------------------------

  ! Vector length
  FUNCTION vabs_dp(v)
    ! vector length
    REAL(dp), DIMENSION(:), INTENT(IN) :: v
    REAL(dp)                           :: vabs_dp

    vabs_dp = sqrt(dot_product(v,v))

  END FUNCTION vabs_dp


  FUNCTION vabs_sp(v)
    ! vector length
    REAL(sp), DIMENSION(:), INTENT(IN) :: v
    REAL(sp)                           :: vabs_sp

    vabs_sp = sqrt(dot_product(v,v))

  END FUNCTION vabs_sp

  ! ------------------------------------------------------------------

END MODULE mo_fit

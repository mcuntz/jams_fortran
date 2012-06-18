MODULE mo_boxcox

  ! This module contains routines to calculate the Box-Cox transformation
  ! as well as estimating the best exponent for the Box-Cox transformation

  ! Usage:
  ! USE mo_boxcox, ONLY: boxcox, get_boxcox
  ! lmbda    = get_boxcox(data)
  ! new_data = boxcox(data, lmbda)
  ! data     = invboxcox(new_data, lmbda)

  ! Written March 2011, Matthias Cuntz
  !   - modified Python code of Travis Oliphant (2002): boxcox, llf_boxcox, get_boxcox
  !   - modified numerical recipes: brent, mnbrak, swap, shft

  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: boxcox     ! Calculate Box-Cox power transformed values given the original values and the exponent lambda
  PUBLIC :: get_boxcox ! Find the exponent that maximises the log-likelihood function
  PUBLIC :: invboxcox  ! Calculate the inverse Box-Cox given the transformed values and the exponent lambda

  INTERFACE boxcox
     MODULE PROCEDURE boxcox_sp, boxcox_dp
  END INTERFACE boxcox
  INTERFACE get_boxcox
     MODULE PROCEDURE get_boxcox_sp, get_boxcox_dp
  END INTERFACE get_boxcox
  INTERFACE invboxcox
     MODULE PROCEDURE invboxcox_sp, invboxcox_dp
  END INTERFACE invboxcox

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         boxcox

  !     PURPOSE
  !         Transform a positive dataset with a Box-Cox power transformation.
  !             boxcox(x) = (x**lambda - 1)/lambda    if lambda > 0
  !                         ln(x)                     if lambda = 0
  !
  !         If an optinal mask is given, then the Box-Cox transformation is only performed on
  !         those locations that correspond to true values in the mask.
  !         x can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = boxcox(x, lmbda, mask=mask)
  
  !     INDENT(IN)
  !         real(sp/dp) :: x(:)       1D-array with input numbers (>0.)
  !         real(sp/dp) :: lmbda      Exponent power of Box-Cox transform (>= 0.)

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(sp/dp) :: boxcox     Power transformed values (at mask=.True.)

  !     INDENT(IN), OPTIONAL
  !         logical :: mask(:)        1D-array of logical values with size(x).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be positive, i.e. > 0.

  !     EXAMPLE
  !         lmbda    = get_boxcox(data, mask=(data>0.))
  !         new_data = boxcox(data, lmbda, mask=(data>0.))
  !         idata    = invboxcox(new_data, lmbda, mask=(data>0.))
  !         -> see also example in test directory

  !     HISTORY
  !         Written,  Matthias Cuntz, March 2011
  !            - Modified Python code of Travis Oliphant (2002): boxcox, llf_boxcox, get_boxcox
  !            - Modified numerical recipes: brent, mnbrak, swap, shft

  FUNCTION boxcox_sp(x, lmbda, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(in)           :: x
    REAL(sp),               INTENT(in)           :: lmbda
    LOGICAL,  DIMENSION(:), INTENT(in), OPTIONAL :: mask
    REAL(sp), DIMENSION(size(x))                 :: boxcox_sp

    LOGICAL, DIMENSION(size(x)) :: maske
    REAL(sp)                    :: lmbda1

    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error boxcox_sp: size(mask) /= size(x)'
       maske = mask
    endif
    if (any((x <= 0.0_sp) .and. maske)) stop 'Error boxcox_sp: x <= 0'
    if (lmbda == 0.0_sp) then
       where (maske)
          boxcox_sp = log(x)
       elsewhere
          boxcox_sp = x
       end where
    else
       lmbda1    = 1.0_sp / lmbda
       boxcox_sp = merge((exp(lmbda*log(x)) - 1.0_sp) * lmbda1, x, maske)
    endif

  END FUNCTION boxcox_sp


  FUNCTION boxcox_dp(x, lmbda, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(in)           :: x
    REAL(dp),               INTENT(in)           :: lmbda
    LOGICAL,  DIMENSION(:), INTENT(in), OPTIONAL :: mask
    REAL(dp), DIMENSION(size(x))                 :: boxcox_dp

    LOGICAL, DIMENSION(size(x)) :: maske
    REAL(dp)                    :: lmbda1

    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error boxcox_dp: size(mask) /= size(x)'
       maske = mask
    endif
    if (any((x <= 0.0_dp) .and. maske)) stop 'Error boxcox_dp: x <= 0'
    if (lmbda == 0.0_dp) then
       where (maske)
          boxcox_dp = log(x)
       elsewhere
          boxcox_dp = x
       end where
    else
       lmbda1    = 1.0_dp / lmbda
       boxcox_dp = merge((exp(lmbda*log(x)) - 1.0_dp) * lmbda1, x, maske)
    endif

  END FUNCTION boxcox_dp

  ! ------------------------------------------------------------------

  ! Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is 
  ! between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates 
  ! the minimum to a fractional precision of about tol using Brent''s method. The abscissa of 
  ! the minimum is returned.
  ! Parameters: Maximum allowed number of iterations; goldenratio; and a small number that 
  ! protects against trying to achieve fractional accuracy for a minimum that happens to be 
  ! exactly zero.
  Function brent_sp(ax, bx, cx, func_sp, tol, dat)

    IMPLICIT NONE

    REAL(sp),               INTENT(in) :: ax, bx, cx
    INTERFACE
       FUNCTION func_sp(x, da)
         USE mo_kind, ONLY: sp
         IMPLICIT NONE
         REAL(sp),               INTENT(in) :: x
         REAL(sp), DIMENSION(:), INTENT(in) :: da
         REAL(sp)                           :: func_sp
       END FUNCTION func_sp
    END INTERFACE
    REAL(sp),               INTENT(in) :: tol
    REAL(sp), DIMENSION(:), INTENT(in) :: dat
    REAL(sp)                           :: brent_sp

    INTEGER(i4), PARAMETER :: ITMAX = 500
    REAL(sp),    PARAMETER :: CGOLD = 0.3819660_sp
    REAL(sp),    PARAMETER :: ZEPS  = 1.0e-3_sp*epsilon(ax)
    INTEGER(i4) :: iter
    REAL(sp)    :: a, b, d, e, etemp, fu, fv, fw, fx
    REAL(sp)    :: p, q, r, tol1, tol2, u, v, w, x, xm

    d = 0.0_sp

    a  = min(ax,cx)
    b  = max(ax,cx)
    v  = bx
    w  = v
    x  = v
    e  = 0.0_sp
    fx = func_sp(x, dat)
    fv = fx
    fw = fx
    do iter=1, ITMAX
       xm   = 0.5_sp*(a+b)
       tol1 = tol*abs(x) + ZEPS
       tol2 = 2.0_sp*tol1
       if (abs(x-xm) <= (tol2-0.5_sp*(b-a))) then
          brent_sp = x
          return
       end if
       if (abs(e) > tol1) then
          r = (x-w)*(fx-fv)
          q = (x-v)*(fx-fw)
          p = (x-v)*q-(x-w)*r
          q = 2.0_sp*(q-r)
          if (q > 0.0_sp) p = -p
          q     = abs(q)
          etemp = e
          e     = d
          if (abs(p) >= abs(0.5_sp*q*etemp) .or. &
               p <= q*(a-x) .or. p >= q*(b-x)) then
             e = merge(a-x, b-x, x>=xm)
             d = CGOLD*e
          else
             d = p/q
             u = x+d
             if (u-a < tol2 .or. b-u < tol2) d = sign(tol1,xm-x)
          end if
       else
          e = merge(a-x, b-x, x>=xm)
          d = CGOLD*e
       end if
       u  = merge(x+d, x+sign(tol1,d), abs(d)>=tol1)
       fu = func_sp(u, dat)
       if (fu <= fx) then
          if (u >= x) then
             a = x
          else
             b = x
          end if
          call shft_sp(v,w,x,u)
          call shft_sp(fv,fw,fx,fu)
       else
          if (u < x) then
             a = u
          else
             b = u
          end if
          if (fu <= fw .or. w == x) then
             v  = w
             fv = fw
             w  = u
             fw = fu
          else if (fu <= fv .or. v == x .or. v == w) then
             v  = u
             fv = fu
          end if
       end if
    end do
    ! If it comes to here, there were too much iterations
    STOP 'brent_sp: exceed maximum iterations'

  END FUNCTION brent_sp


  FUNCTION brent_dp(ax, bx, cx, func_dp, tol, dat)

    IMPLICIT NONE

    REAL(dp),               INTENT(in) :: ax, bx, cx
    INTERFACE
       FUNCTION func_dp(x, da)
         USE mo_kind, ONLY: dp
         IMPLICIT NONE
         REAL(dp),               INTENT(in) :: x
         REAL(dp), DIMENSION(:), INTENT(in) :: da
         REAL(dp)                           :: func_dp
       END FUNCTION func_dp
    END INTERFACE
    REAL(dp),               INTENT(in) :: tol
    REAL(dp), DIMENSION(:), INTENT(in) :: dat
    REAL(dp)                           :: brent_dp

    INTEGER(i4), PARAMETER :: ITMAX=500
    REAL(dp),    PARAMETER :: CGOLD=0.3819660_dp
    REAL(dp),    PARAMETER :: ZEPS=1.0e-3_dp*epsilon(ax)
    INTEGER(i4) :: iter
    REAL(dp)    :: a, b, d, e, etemp, fu, fv, fw, fx
    REAL(dp)    :: p, q, r, tol1, tol2, u, v, w, x, xm

    d = 0.0_dp

    a  = min(ax,cx)
    b  = max(ax,cx)
    v  = bx
    w  = v
    x  = v
    e  = 0.0_dp
    fx = func_dp(x, dat)
    fv = fx
    fw = fx
    do iter=1, ITMAX
       xm   = 0.5_dp*(a+b)
       tol1 = tol*abs(x) + ZEPS
       tol2 = 2.0_dp*tol1
       if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
          brent_dp = x
          return
       end if
       if (abs(e) > tol1) then
          r = (x-w)*(fx-fv)
          q = (x-v)*(fx-fw)
          p = (x-v)*q-(x-w)*r
          q = 2.0_dp*(q-r)
          if (q > 0.0_dp) p = -p
          q     = abs(q)
          etemp = e
          e     = d
          if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
               p <= q*(a-x) .or. p >= q*(b-x)) then
             e = merge(a-x, b-x, x>=xm)
             d = CGOLD*e
          else
             d = p/q
             u = x+d
             if (u-a < tol2 .or. b-u < tol2) d = sign(tol1,xm-x)
          end if
       else
          e = merge(a-x, b-x, x>=xm)
          d = CGOLD*e
       end if
       u  = merge(x+d, x+sign(tol1,d), abs(d)>=tol1)
       fu = func_dp(u, dat)
       if (fu <= fx) then
          if (u >= x) then
             a = x
          else
             b = x
          end if
          call shft_dp(v,w,x,u)
          call shft_dp(fv,fw,fx,fu)
       else
          if (u < x) then
             a = u
          else
             b = u
          end if
          if (fu <= fw .or. w == x) then
             v  = w
             fv = fw
             w  = u
             fw = fu
          else if (fu <= fv .or. v == x .or. v == w) then
             v  = u
             fv = fu
          end if
       end if
    end do
    ! If it comes to here, there were too much iterations
    STOP 'brent_dp: exceed maximum iterations'

  END FUNCTION brent_dp

  ! ------------------------------------------------------------------

  ! Find the lambda that maximizes the log-likelihood function
  FUNCTION get_boxcox_sp(x, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(in)           :: x
    LOGICAL,  DIMENSION(:), INTENT(in), OPTIONAL :: mask
    REAL(sp)                                     :: get_boxcox_sp

    REAL(sp), PARAMETER :: tol = 1.48e-08
    REAL(sp) :: ax, bx, cx
    REAL(sp) :: fa, fb, fc
    LOGICAL,  DIMENSION(size(x))        :: maske
    REAL(sp), DIMENSION(:), ALLOCATABLE :: xx
    INTEGER(i4) :: nn
    

    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error get_boxcox_sp: size(mask) /= size(x)'
       maske = mask
    endif
    nn = count(maske)
    if (allocated(xx)) deallocate(xx)
    allocate(xx(nn))
    xx = pack(x, mask=maske)
    ax = -2.0_sp
    bx = 2.0_sp
    call mnbrak_sp(ax, bx, cx, fa, fb, fc, llf_boxcox_sp, xx)
    get_boxcox_sp = brent_sp(ax, bx, cx, llf_boxcox_sp, tol, xx)

  END FUNCTION get_boxcox_sp


  FUNCTION get_boxcox_dp(x, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(in)           :: x
    LOGICAL,  DIMENSION(:), INTENT(in), OPTIONAL :: mask
    REAL(dp)                                     :: get_boxcox_dp

    REAL(dp), PARAMETER :: tol = 1.48e-08
    REAL(dp) :: ax, bx, cx
    REAL(dp) :: fa, fb, fc
    LOGICAL,  DIMENSION(size(x))        :: maske
    REAL(dp), DIMENSION(:), ALLOCATABLE :: xx
    INTEGER(i4) :: nn
    

    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error get_boxcox_dp: size(mask) /= size(x)'
       maske = mask
    endif
    nn = count(maske)
    if (allocated(xx)) deallocate(xx)
    allocate(xx(nn))
    xx = pack(x, mask=maske)
    ax = -2.0_dp
    bx = 2.0_dp
    call mnbrak_dp(ax, bx, cx, fa, fb, fc, llf_boxcox_dp, xx)
    get_boxcox_dp = brent_dp(ax, bx, cx, llf_boxcox_dp, tol, xx)

  END FUNCTION get_boxcox_dp

  ! ------------------------------------------------------------------

  !     NAME
  !         invboxcox

  !     PURPOSE
  !         Back-transformation of Box-Cox-transformed data.
  !             boxcox(x)    = (x**lambda - 1)/lambda        if lambda > 0
  !                            ln(x)                         if lambda = 0
  !             invboxcox(y) = (lambda*y + 1)**(1/lambda)    if lambda > 0
  !                            exp(y)                        if lambda = 0
  !
  !         If an optinal mask is given, then the inverse Box-Cox transformation is only performed on
  !         those locations that correspond to true values in the mask.
  !         x can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = invboxcox(x, lmbda, mask=mask)
  
  !     INDENT(IN)
  !         real(sp/dp) :: x(:)       1D-array with input numbers (>0.)
  !         real(sp/dp) :: lmbda      Exponent power of Box-Cox transform (>= 0.)

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(sp/dp) :: boxcox     Back transformed values (at mask=.True.)

  !     INDENT(IN), OPTIONAL
  !         logical :: mask(:)        1D-array of logical values with size(x).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         lmbda    = get_boxcox(data, mask=(data>0.))
  !         new_data = boxcox(data, lmbda, mask=(data>0.))
  !         idata    = invboxcox(new_data, lmbda, mask=(data>0.))
  !         -> see also example in test directory

  !     HISTORY
  !         Written,  Matthias Cuntz, March 2011
  !            - Modified Python code of Travis Oliphant (2002): boxcox, llf_boxcox, get_boxcox
  !            - Modified numerical recipes: brent, mnbrak, swap, shft

  FUNCTION invboxcox_sp(x, lmbda, mask)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(in)           :: x
    REAL(sp),               INTENT(in)           :: lmbda
    LOGICAL,  DIMENSION(:), INTENT(in), OPTIONAL :: mask
    REAL(sp), DIMENSION(size(x))                 :: invboxcox_sp

    LOGICAL, DIMENSION(size(x)) :: maske
    REAL(sp)                    :: lmbda1

    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error invboxcox_sp: size(mask) /= size(x)'
       maske = mask
    endif
    if (lmbda == 0.0_sp) then
       where (maske)
          invboxcox_sp = exp(x)
       elsewhere
          invboxcox_sp = x
       end where
    else
       lmbda1    = 1.0_sp / lmbda
       invboxcox_sp = merge(exp(lmbda1*log(lmbda*x+1.0_sp)), x, maske)
    endif

  END FUNCTION invboxcox_sp

  FUNCTION invboxcox_dp(x, lmbda, mask)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(in)           :: x
    REAL(dp),               INTENT(in)           :: lmbda
    LOGICAL,  DIMENSION(:), INTENT(in), OPTIONAL :: mask
    REAL(dp), DIMENSION(size(x))                 :: invboxcox_dp

    LOGICAL, DIMENSION(size(x)) :: maske
    REAL(dp)                    :: lmbda1

    maske(:) = .true.
    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error invboxcox_dp: size(mask) /= size(x)'
       maske = mask
    endif
    if (lmbda == 0.0_dp) then
       where (maske)
          invboxcox_dp = exp(x)
       elsewhere
          invboxcox_dp = x
       end where
    else
       lmbda1    = 1.0_dp / lmbda
       invboxcox_dp = merge(exp(lmbda1*log(lmbda*x+1.0_dp)), x, maske)
    endif

  END FUNCTION invboxcox_dp

  ! ------------------------------------------------------------------

  ! The boxcox (negative) log-likelihood function.
  FUNCTION llf_boxcox_sp(lmbda, x)

    IMPLICIT NONE

    REAL(sp),               INTENT(in) :: lmbda
    REAL(sp), DIMENSION(:), INTENT(in) :: x
    REAL(sp)                           :: llf_boxcox_sp

    REAL(sp), DIMENSION(size(x)) :: y
    REAL(sp)                     :: N, my, f

    N             = real(size(x),sp)
    y             = boxcox_sp(x, lmbda)
    my            = sum(y) / N
    f             = (lmbda-1.0_sp) * sum(log(x))
    f             = f - 0.5_sp*N * log(sum((y-my)*(y-my))/N)
    llf_boxcox_sp = -f

  END FUNCTION llf_boxcox_sp


  FUNCTION llf_boxcox_dp(lmbda, x)

    IMPLICIT NONE

    REAL(dp),               INTENT(in) :: lmbda
    REAL(dp), DIMENSION(:), INTENT(in) :: x
    REAL(dp)                           :: llf_boxcox_dp

    REAL(dp), DIMENSION(size(x)) :: y
    REAL(dp)                     :: N, my, f

    N             = real(size(x),dp)
    y             = boxcox_dp(x, lmbda)
    my            = sum(y) / N
    f             = (lmbda-1.0_dp) * sum(log(x))
    f             = f - 0.5_dp*N * log(sum((y-my)*(y-my))/N)
    llf_boxcox_dp = -f

  END FUNCTION llf_boxcox_dp

  ! ------------------------------------------------------------------

  ! Given a function func, and given distinct initial points ax and bx, this routine searches 
  ! in the downhill direction (defined by the functionas evaluated at the initial points) and 
  ! returns new points ax, bx, cx that bracket a minimum of the function.
  ! Parameters: GOLD is the default ratio by which successive intervals are magnified; GLIMIT 
  ! is the maximum magnification allowed for a parabolic-fit step.
  SUBROUTINE mnbrak_sp(ax, bx, cx, fa, fb, fc, func_sp, dat)

    IMPLICIT NONE

    REAL(sp),               INTENT(inout) :: ax, bx
    REAL(sp),               INTENT(out)   :: cx
    REAL(sp),               INTENT(out)   :: fa, fb, fc
    INTERFACE
       FUNCTION func_sp(x, da)
         USE mo_kind, ONLY: sp
         IMPLICIT NONE
         REAL(sp),               INTENT(in) :: x
         REAL(sp), DIMENSION(:), INTENT(in) :: da
         REAL(sp)                           :: func_sp
       END FUNCTION func_sp
    END INTERFACE
    REAL(sp), DIMENSION(:), INTENT(in) :: dat

    REAL(sp), PARAMETER :: GOLD   = 1.618034_sp
    REAL(sp), PARAMETER :: GLIMIT = 100.0_sp
    REAL(sp), PARAMETER :: TINY   = 1.0e-20_sp
    REAL(sp) :: fu, q, r, u, ulim

    fa = func_sp(ax, dat)
    fb = func_sp(bx, dat)
    if (fb > fa) then
       call swap_sp(ax,bx)
       call swap_sp(fa,fb)
    end if
    cx = bx + GOLD*(bx-ax)
    fc = func_sp(cx, dat)
    do
       if (fb < fc) return
       r    = (bx-ax)*(fb-fc)
       q    = (bx-cx)*(fb-fa)
       u    = bx-((bx-cx)*q-(bx-ax)*r)/(2.0_sp*sign(max(abs(q-r),TINY),q-r))
       ulim = bx + GLIMIT*(cx-bx)
       if ((bx-u)*(u-cx) > 0.0_sp) then
          fu = func_sp(u, dat)
          if (fu < fc) then
             ax = bx
             fa = fb
             bx = u
             fb = fu
             return
          else if (fu > fb) then
             cx = u
             fc = fu
             return
          end if
          u  = cx + GOLD*(cx-bx)
          fu = func_sp(u, dat)
       else if ((cx-u)*(u-ulim) > 0.0_sp) then
          fu = func_sp(u, dat)
          if (fu < fc) then
             bx = cx
             cx = u
             u  = cx + GOLD*(cx-bx)
             call shft_sp(fb,fc,fu,func_sp(u, dat))
          end if
       else if ((u-ulim)*(ulim-cx) >= 0.0_sp) then
          u  = ulim
          fu = func_sp(u, dat)
       else
          u  = cx + GOLD*(cx-bx)
          fu = func_sp(u, dat)
       end if
       call shft_sp(ax,bx,cx,u)
       call shft_sp(fa,fb,fc,fu)
    end do

  END SUBROUTINE mnbrak_sp


  SUBROUTINE mnbrak_dp(ax, bx, cx, fa, fb, fc, func_dp, dat)

    IMPLICIT NONE

    REAL(dp),               INTENT(inout) :: ax, bx
    REAL(dp),               INTENT(out)   :: cx
    REAL(dp),               INTENT(out)   :: fa, fb, fc
    INTERFACE
       FUNCTION func_dp(x, da)
         USE mo_kind, ONLY: dp
         IMPLICIT NONE
         REAL(dp),               INTENT(in) :: x
         REAL(dp), DIMENSION(:), INTENT(in) :: da
         REAL(dp)                           :: func_dp
       END FUNCTION func_dp
    END INTERFACE
    REAL(dp), DIMENSION(:), INTENT(in) :: dat

    REAL(dp), PARAMETER :: GOLD   = 1.618034_dp
    REAL(dp), PARAMETER :: GLIMIT = 100.0_dp
    REAL(dp), PARAMETER :: TINY   = 1.0e-20_dp
    REAL(dp) :: fu, q, r, u, ulim

    fa = func_dp(ax, dat)
    fb = func_dp(bx, dat)
    if (fb > fa) then
       call swap_dp(ax,bx)
       call swap_dp(fa,fb)
    end if
    cx = bx + GOLD*(bx-ax)
    fc = func_dp(cx, dat)
    do
       if (fb < fc) return
       r    = (bx-ax)*(fb-fc)
       q    = (bx-cx)*(fb-fa)
       u    = bx-((bx-cx)*q-(bx-ax)*r)/(2.0_dp*sign(max(abs(q-r),TINY),q-r))
       ulim = bx + GLIMIT*(cx-bx)
       if ((bx-u)*(u-cx) > 0.0_dp) then
          fu = func_dp(u, dat)
          if (fu < fc) then
             ax = bx
             fa = fb
             bx = u
             fb = fu
             return
          else if (fu > fb) then
             cx = u
             fc = fu
             return
          end if
          u  = cx + GOLD*(cx-bx)
          fu = func_dp(u, dat)
       else if ((cx-u)*(u-ulim) > 0.0_dp) then
          fu = func_dp(u, dat)
          if (fu < fc) then
             bx = cx
             cx = u
             u  = cx + GOLD*(cx-bx)
             call shft_dp(fb,fc,fu,func_dp(u, dat))
          end if
       else if ((u-ulim)*(ulim-cx) >= 0.0_dp) then
          u  = ulim
          fu = func_dp(u, dat)
       else
          u  = cx + GOLD*(cx-bx)
          fu = func_dp(u, dat)
       end if
       call shft_dp(ax,bx,cx,u)
       call shft_dp(fa,fb,fc,fu)
    end do

  END SUBROUTINE mnbrak_dp

  ! ------------------------------------------------------------------

  ! Helper for brent and mnbrak
  SUBROUTINE shft_sp(a,b,c,d)

    IMPLICIT NONE

    REAL(sp), INTENT(out)   :: a
    REAL(sp), INTENT(inout) :: b,c
    REAL(sp), INTENT(in)    :: d

    a = b
    b = c
    c = d

  END SUBROUTINE shft_sp


  SUBROUTINE shft_dp(a,b,c,d)

    IMPLICIT NONE

    REAL(dp), INTENT(out)   :: a
    REAL(dp), INTENT(inout) :: b,c
    REAL(dp), INTENT(in)    :: d

    a = b
    b = c
    c = d

  END SUBROUTINE shft_dp

  ! ------------------------------------------------------------------

  ! Helper for mnbrak
  SUBROUTINE swap_sp(a,b)

    IMPLICIT NONE

    REAL(sp), INTENT(inout) :: a, b
    REAL(sp) :: dum
    dum = a
    a   = b
    b   = dum

  END SUBROUTINE swap_sp


  SUBROUTINE swap_dp(a,b)

    IMPLICIT NONE

    REAL(dp), INTENT(inout) :: a, b
    REAL(dp) :: dum
    dum = a
    a   = b
    b   = dum

  END SUBROUTINE swap_dp

END MODULE mo_boxcox

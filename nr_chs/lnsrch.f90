SUBROUTINE lnsrch_sp(xold,fold,g,p,x,f,stpmax,check,func)

  use mo_kind,   only: sp, i4, lgt
  use mo_nrutil, only: assert_eq,nrerror,vabs

  IMPLICIT NONE

  REAL(SP), DIMENSION(:), INTENT(IN)    :: xold,g
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
  REAL(SP),               INTENT(IN)    :: fold,stpmax
  REAL(SP), DIMENSION(:), INTENT(OUT)   :: x
  REAL(SP),               INTENT(OUT)   :: f
  LOGICAL(lgt),           INTENT(OUT)   :: check

  INTERFACE
     FUNCTION func(x)
       use mo_kind, only: sp
       IMPLICIT NONE 
       REAL(SP)                           :: func
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
     END FUNCTION func
  END INTERFACE

  REAL(SP), PARAMETER :: ALF=1.0e-4_sp,TOLX=epsilon(x)
  INTEGER(I4)         :: ndum
  REAL(SP)            :: a,alam,alam2,alamin,b,disc,f2,fold2
  REAL(SP)            :: pabs,rhs1,rhs2,slope,tmplam

  ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
  check=.false.
  pabs=vabs(p(:))
  if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
  slope=dot_product(g,p)
  alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_sp))
  alam=1.0_sp
  do
     x(:)=xold(:)+alam*p(:)
     f=func(x)
     if (alam < alamin) then
        x(:)=xold(:)
        check=.true.
        RETURN
     else if (f <= fold+ALF*alam*slope) then
        RETURN
     else
        if (alam == 1.0_sp) then
           tmplam=-slope/(2.0_sp*(f-fold-slope))
        else
           rhs1=f-fold-alam*slope
           rhs2=f2-fold2-alam2*slope
           a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
           b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
                (alam-alam2)
           if (a == 0.0_sp) then
              tmplam=-slope/(2.0_sp*b)
           else
              disc=b*b-3.0_sp*a*slope
              if (disc < 0.0_sp) call nrerror('roundoff problem in lnsrch')
              tmplam=(-b+sqrt(disc))/(3.0_sp*a)
           end if
           if (tmplam > 0.5_sp*alam) tmplam=0.5_sp*alam
        end if
     end if
     alam2=alam
     f2=f
     fold2=fold
     alam=max(tmplam,0.1_sp*alam)
  end do
END SUBROUTINE lnsrch_sp

SUBROUTINE lnsrch_dp(xold,fold,g,p,x,f,stpmax,check,func)

  use mo_kind,   only: dp, i4, lgt
  use mo_nrutil, only: assert_eq,nrerror,vabs

  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(IN)    :: xold,g
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
  REAL(DP),               INTENT(IN)    :: fold,stpmax
  REAL(DP), DIMENSION(:), INTENT(OUT)   :: x
  REAL(DP),               INTENT(OUT)   :: f
  LOGICAL(lgt),           INTENT(OUT)   :: check

  INTERFACE
     FUNCTION func(x)
       use mo_kind, only: dp
       IMPLICIT NONE
       REAL(DP)                           :: func
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
     END FUNCTION func
  END INTERFACE

  REAL(DP), PARAMETER :: ALF=1.0e-4_dp,TOLX=epsilon(x)
  INTEGER(I4)         :: ndum
  REAL(DP)            :: a,alam,alam2,alamin,b,disc,f2,fold2
  REAL(DP)            :: pabs,rhs1,rhs2,slope,tmplam

  ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
  check=.false.
  pabs=vabs(p(:))
  if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
  slope=dot_product(g,p)
  alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_dp))
  alam=1.0_dp
  do
     x(:)=xold(:)+alam*p(:)
     f=func(x)
     if (alam < alamin) then
        x(:)=xold(:)
        check=.true.
        RETURN
     else if (f <= fold+ALF*alam*slope) then
        RETURN
     else
        if (alam == 1.0_dp) then
           tmplam=-slope/(2.0_dp*(f-fold-slope))
        else
           rhs1=f-fold-alam*slope
           rhs2=f2-fold2-alam2*slope
           a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
           b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
                (alam-alam2)
           if (a == 0.0_dp) then
              tmplam=-slope/(2.0_dp*b)
           else
              disc=b*b-3.0_dp*a*slope
              if (disc < 0.0_dp) call nrerror('roundoff problem in lnsrch')
              tmplam=(-b+sqrt(disc))/(3.0_dp*a)
           end if
           if (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
        end if
     end if
     alam2=alam
     f2=f
     fold2=fold
     alam=max(tmplam,0.1_dp*alam)
  end do
END SUBROUTINE lnsrch_dp


SUBROUTINE dfpmin_sp(p,gtol,iter,fret,func,dfunc)
  use mo_kind,   only: i4, sp
  use mo_nrutil, only: nrerror,outerprod,unit_matrix,vabs
  use mo_nr,     only: lnsrch

  IMPLICIT NONE

  REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
  REAL(SP),               INTENT(IN)    :: gtol
  INTEGER(I4),            INTENT(OUT)   :: iter
  REAL(SP),               INTENT(OUT)   :: fret

  INTERFACE
     FUNCTION func(p)
       use mo_kind, only: sp
       IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: p
       REAL(SP) :: func
     END FUNCTION func
     !BL
     FUNCTION dfunc(p)
       use mo_kind, only: sp
       IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: p
       REAL(SP), DIMENSION(size(p)) :: dfunc
     END FUNCTION dfunc
  END INTERFACE

  ! local variables
  INTEGER(I4), PARAMETER               :: ITMAX=200
  REAL(SP),    PARAMETER               :: STPMX=100.0_sp,EPS=epsilon(p),TOLX=4.0_sp*EPS
  INTEGER(I4)                          :: its
  LOGICAL                              :: check
  REAL(SP)                             :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
  REAL(SP), DIMENSION(size(p))         :: dg,g,hdg,pnew,xi
  REAL(SP), DIMENSION(size(p),size(p)) :: hessin

  fp=func(p)
  g=dfunc(p)
  call unit_matrix(hessin)
  xi=-g
  stpmax=STPMX*max(vabs(p),real(size(p),sp))
  do its=1,ITMAX
     iter=its
     call lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,func)
     fp=fret
     xi=pnew-p
     p=pnew
     if (maxval(abs(xi)/max(abs(p),1.0_sp)) < TOLX) RETURN
     dg=g
     g=dfunc(p)
     den=max(fret,1.0_sp)
     if (maxval(abs(g)*max(abs(p),1.0_sp)/den) < gtol) RETURN
     dg    = g-dg
     hdg   = matmul(hessin,dg)
     fac   = dot_product(dg,xi)
     fae   = dot_product(dg,hdg)
     sumdg = dot_product(dg,dg)
     sumxi = dot_product(xi,xi)
     if (fac**2 > EPS*sumdg*sumxi) then
        fac    = 1.0_sp/fac
        fad    = 1.0_sp/fae
        dg     = fac*xi-fad*hdg
        hessin = hessin+fac*outerprod(xi,xi)-&
             fad*outerprod(hdg,hdg)+fae*outerprod(dg,dg)
     end if
     xi = -matmul(hessin,g)
  end do
  call nrerror('dfpmin: too many iterations')
END SUBROUTINE dfpmin_sp

SUBROUTINE dfpmin_dp(p,gtol,iter,fret,func,dfunc)
  use mo_kind,   only: i4, dp
  use mo_nrutil, only: nrerror,outerprod,unit_matrix,vabs
  use mo_nr,     only: lnsrch

  IMPLICIT NONE

  REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
  REAL(DP),               INTENT(IN)    :: gtol
  INTEGER(I4),            INTENT(OUT)   :: iter
  REAL(DP),               INTENT(OUT)   :: fret

  INTERFACE
     FUNCTION func(p)
       use mo_kind, only: dp
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: p
       REAL(DP) :: func
     END FUNCTION func
     !BL
     FUNCTION dfunc(p)
       use mo_kind, only: dp
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(IN) :: p
       REAL(DP), DIMENSION(size(p)) :: dfunc
     END FUNCTION dfunc
  END INTERFACE

  ! local variables
  INTEGER(I4), PARAMETER               :: ITMAX=200
  REAL(DP),    PARAMETER               :: STPMX=100.0_dp,EPS=epsilon(p),TOLX=4.0_dp*EPS
  INTEGER(I4)                          :: its
  LOGICAL                              :: check
  REAL(DP)                             :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
  REAL(DP), DIMENSION(size(p))         :: dg,g,hdg,pnew,xi
  REAL(DP), DIMENSION(size(p),size(p)) :: hessin

  fp=func(p)
  g=dfunc(p)
  call unit_matrix(hessin)
  xi=-g
  stpmax=STPMX*max(vabs(p),real(size(p),dp))
  do its=1,ITMAX
     iter=its
     call lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,func)
     fp=fret
     xi=pnew-p 
     p=pnew
     if (maxval(abs(xi)/max(abs(p),1.0_dp)) < TOLX) RETURN
     dg=g
     g=dfunc(p)
     den=max(fret,1.0_dp)
     if (maxval(abs(g)*max(abs(p),1.0_dp)/den) < gtol) RETURN
     dg    = g-dg
     hdg   = matmul(hessin,dg)
     fac   = dot_product(dg,xi)
     fae   = dot_product(dg,hdg)
     sumdg = dot_product(dg,dg)
     sumxi = dot_product(xi,xi)
     if (fac**2 > EPS*sumdg*sumxi) then
        fac    = 1.0_dp/fac
        fad    = 1.0_dp/fae
        dg     = fac*xi-fad*hdg
        hessin = hessin+fac*outerprod(xi,xi)-&
             fad*outerprod(hdg,hdg)+fae*outerprod(dg,dg)
     end if
     xi = -matmul(hessin,g)
  end do
  call nrerror('dfpmin: too many iterations')
END SUBROUTINE dfpmin_dp


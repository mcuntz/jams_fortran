FUNCTION gser_s(a,x,gln)
  USE mo_kind, only: sp, i4
  USE mo_nrutil, ONLY : nrerror
  use mo_nr, ONLY : gammln
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,x
  REAL(SP), OPTIONAL, INTENT(OUT) :: gln
  REAL(SP) :: gser_s
  INTEGER(I4), PARAMETER :: ITMAX=100
  REAL(SP), PARAMETER :: EPS=epsilon(x)
  INTEGER(I4) :: n
  REAL(SP) :: ap,del,summ
  if (x == 0.0) then
     gser_s=0.0
     RETURN
  end if
  ap=a
  summ=1.0_sp/a
  del=summ
  do n=1,ITMAX
     ap=ap+1.0_sp
     del=del*x/ap
     summ=summ+del
     if (abs(del) < abs(summ)*EPS) exit
  end do
  if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_s')
  if (present(gln)) then
     gln=gammln(a)
     gser_s=summ*exp(-x+a*log(x)-gln)
  else
     gser_s=summ*exp(-x+a*log(x)-gammln(a))
  end if
END FUNCTION gser_s


FUNCTION gser_v(a,x,gln)
  USE mo_kind, only: sp, i4, lgt
  USE mo_nrutil, ONLY : assert_eq,nrerror
  use mo_nr, ONLY : gammln
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
  REAL(SP), DIMENSION(size(a)) :: gser_v
  INTEGER(I4), PARAMETER :: ITMAX=100
  REAL(SP), PARAMETER :: EPS=epsilon(x)
  INTEGER(I4) :: n
  REAL(SP), DIMENSION(size(a)) :: ap,del,summ
  LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
  n=assert_eq(size(a),size(x),'gser_v')
  zero=(x == 0.0)
  where (zero) gser_v=0.0
  ap=a
  summ=1.0_sp/a
  del=summ
  converged=zero
  do n=1,ITMAX
     where (.not. converged)
        ap=ap+1.0_sp
        del=del*x/ap
        summ=summ+del
        converged = (abs(del) < abs(summ)*EPS)
     end where
     if (all(converged)) exit
  end do
  if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_v')
  if (present(gln)) then
     if (size(gln) < size(a)) call &
          nrerror('gser_v: Not enough space for gln')
     gln=gammln(a)
     where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gln)
  else
     where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gammln(a))
  end if
END FUNCTION gser_v


FUNCTION dgser_s(a,x,gln)
  USE mo_kind, only: dp, i4
  USE mo_nrutil, ONLY : nrerror
  use mo_nr, ONLY : gammln
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,x
  REAL(DP), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP) :: dgser_s
  INTEGER(I4), PARAMETER :: ITMAX=100
  REAL(DP), PARAMETER :: EPS=epsilon(x)
  INTEGER(I4) :: n
  REAL(DP) :: ap,del,summ
  if (x == 0.0_dp) then
     dgser_s=0.0_dp
     RETURN
  end if
  ap=a
  summ=1.0_dp/a
  del=summ
  do n=1,ITMAX
     ap=ap+1.0_dp
     del=del*x/ap
     summ=summ+del
     if (abs(del) < abs(summ)*EPS) exit
  end do
  if (n > ITMAX) call nrerror('a too large, ITMAX too small in dgser_s')
  if (present(gln)) then
     gln=gammln(a)
     dgser_s=summ*exp(-x+a*log(x)-gln)
  else
     dgser_s=summ*exp(-x+a*log(x)-gammln(a))
  end if
END FUNCTION dgser_s


FUNCTION dgser_v(a,x,gln)
  USE mo_kind, only: dp, i4, lgt
  USE mo_nrutil, ONLY : assert_eq,nrerror
  use mo_nr, ONLY : gammln
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP), DIMENSION(size(a)) :: dgser_v
  INTEGER(I4), PARAMETER :: ITMAX=100
  REAL(DP), PARAMETER :: EPS=epsilon(x)
  INTEGER(I4) :: n
  REAL(DP), DIMENSION(size(a)) :: ap,del,summ
  LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
  n=assert_eq(size(a),size(x),'dgser_v')
  zero=(x == 0.0_dp)
  where (zero) dgser_v=0.0_dp
  ap=a
  summ=1.0_dp/a
  del=summ
  converged=zero
  do n=1,ITMAX
     where (.not. converged)
        ap=ap+1.0_dp
        del=del*x/ap
        summ=summ+del
        converged = (abs(del) < abs(summ)*EPS)
     end where
     if (all(converged)) exit
  end do
  if (n > ITMAX) call nrerror('a too large, ITMAX too small in dgser_v')
  if (present(gln)) then
     if (size(gln) < size(a)) call &
          nrerror('dgser_v: Not enough space for gln')
     gln=gammln(a)
     where (.not. zero) dgser_v=summ*exp(-x+a*log(x)-gln)
  else
     where (.not. zero) dgser_v=summ*exp(-x+a*log(x)-gammln(a))
  end if
END FUNCTION dgser_v

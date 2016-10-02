FUNCTION gcf_s(a,x,gln)
  USE mo_kind, only: sp, i4
  USE mo_nrutil, ONLY : nrerror
  use mo_nr, ONLY : gammln
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,x
  REAL(SP), OPTIONAL, INTENT(OUT) :: gln
  REAL(SP) :: gcf_s
  INTEGER(I4), PARAMETER :: ITMAX=100
  REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  INTEGER(I4) :: i
  REAL(SP) :: an,b,c,d,del,h
  if (x == 0.0) then
     gcf_s=1.0
     RETURN
  end if
  b=x+1.0_sp-a
  c=1.0_sp/FPMIN
  d=1.0_sp/b
  h=d
  do i=1,ITMAX
     an=-i*(i-a)
     b=b+2.0_sp
     d=an*d+b
     if (abs(d) < FPMIN) d=FPMIN
     c=b+an/c
     if (abs(c) < FPMIN) c=FPMIN
     d=1.0_sp/d
     del=d*c
     h=h*del
     if (abs(del-1.0_sp) <= EPS) exit
  end do
  if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_s')
  if (present(gln)) then
     gln=gammln(a)
     gcf_s=exp(-x+a*log(x)-gln)*h
  else
     gcf_s=exp(-x+a*log(x)-gammln(a))*h
  end if
END FUNCTION gcf_s


FUNCTION gcf_v(a,x,gln)
  USE mo_kind, only: sp, i4, lgt
  USE mo_nrutil, ONLY : assert_eq,nrerror
  use mo_nr, ONLY : gammln
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
  REAL(SP), DIMENSION(size(a)) :: gcf_v
  INTEGER(I4), PARAMETER :: ITMAX=100
  REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  INTEGER(I4) :: i
  REAL(SP), DIMENSION(size(a)) :: an,b,c,d,del,h
  LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
  i=assert_eq(size(a),size(x),'gcf_v')
  zero=(x == 0.0)
  where (zero)
     gcf_v=1.0
  elsewhere
     b=x+1.0_sp-a
     c=1.0_sp/FPMIN
     d=1.0_sp/b
     h=d
  end where
  converged=zero
  do i=1,ITMAX
     where (.not. converged)
        an=-i*(i-a)
        b=b+2.0_sp
        d=an*d+b
        d=merge(FPMIN,d, abs(d)<FPMIN )
        c=b+an/c
        c=merge(FPMIN,c, abs(c)<FPMIN )
        d=1.0_sp/d
        del=d*c
        h=h*del
        converged = (abs(del-1.0_sp)<=EPS)
     end where
     if (all(converged)) exit
  end do
  if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_v')
  if (present(gln)) then
     if (size(gln) < size(a)) call &
          nrerror('gcf_v: Not enough space for gln')
     gln=gammln(a)
     where (.not. zero) gcf_v=exp(-x+a*log(x)-gln)*h
  else
     where (.not. zero) gcf_v=exp(-x+a*log(x)-gammln(a))*h
  end if
END FUNCTION gcf_v


FUNCTION dgcf_s(a,x,gln)
  USE mo_kind, only: dp, i4
  USE mo_nrutil, ONLY : nrerror
  use mo_nr, ONLY : gammln
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,x
  REAL(DP), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP) :: dgcf_s
  INTEGER(I4), PARAMETER :: ITMAX=100
  REAL(DP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  INTEGER(I4) :: i
  REAL(DP) :: an,b,c,d,del,h
  if (x == 0.0_dp) then
     dgcf_s=1.0_dp
     RETURN
  end if
  b=x+1.0_dp-a
  c=1.0_dp/FPMIN
  d=1.0_dp/b
  h=d
  do i=1,ITMAX
     an=-i*(i-a)
     b=b+2.0_dp
     d=an*d+b
     if (abs(d) < FPMIN) d=FPMIN
     c=b+an/c
     if (abs(c) < FPMIN) c=FPMIN
     d=1.0_dp/d
     del=d*c
     h=h*del
     if (abs(del-1.0_dp) <= EPS) exit
  end do
  if (i > ITMAX) call nrerror('a too large, ITMAX too small in dgcf_s')
  if (present(gln)) then
     gln=gammln(a)
     dgcf_s=exp(-x+a*log(x)-gln)*h
  else
     dgcf_s=exp(-x+a*log(x)-gammln(a))*h
  end if
END FUNCTION dgcf_s


FUNCTION dgcf_v(a,x,gln)
  USE mo_kind, only: dp, i4, lgt
  USE mo_nrutil, ONLY : assert_eq,nrerror
  use mo_nr, ONLY : gammln
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP), DIMENSION(size(a)) :: dgcf_v
  INTEGER(I4), PARAMETER :: ITMAX=100
  REAL(DP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  INTEGER(I4) :: i
  REAL(DP), DIMENSION(size(a)) :: an,b,c,d,del,h
  LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
  i=assert_eq(size(a),size(x),'dgcf_v')
  zero=(x == 0.0_dp)
  where (zero)
     dgcf_v=1.0_dp
  elsewhere
     b=x+1.0_dp-a
     c=1.0_dp/FPMIN
     d=1.0_dp/b
     h=d
  end where
  converged=zero
  do i=1,ITMAX
     where (.not. converged)
        an=-i*(i-a)
        b=b+2.0_dp
        d=an*d+b
        d=merge(FPMIN,d, abs(d)<FPMIN )
        c=b+an/c
        c=merge(FPMIN,c, abs(c)<FPMIN )
        d=1.0_dp/d
        del=d*c
        h=h*del
        converged = (abs(del-1.0_dp)<=EPS)
     end where
     if (all(converged)) exit
  end do
  if (i > ITMAX) call nrerror('a too large, ITMAX too small in dgcf_v')
  if (present(gln)) then
     if (size(gln) < size(a)) call &
          nrerror('dgcf_v: Not enough space for gln')
     gln=gammln(a)
     where (.not. zero) dgcf_v=exp(-x+a*log(x)-gln)*h
  else
     where (.not. zero) dgcf_v=exp(-x+a*log(x)-gammln(a))*h
  end if
END FUNCTION dgcf_v

SUBROUTINE asolve_sp(b,x,itrnsp)
  USE mo_kind, only: sp, i4
  USE mo_nrutil, ONLY : assert_eq,nrerror
  use mo_nr, ONLY : sprsdiag
  USE xlinbcg_data_sp
  REAL(SP), DIMENSION(:), INTENT(IN) :: b
  REAL(SP), DIMENSION(:), INTENT(OUT) :: x
  INTEGER(I4), INTENT(IN) :: itrnsp
  INTEGER(I4) :: ndum
  ndum=assert_eq(size(b),size(x),'asolve_sp')
  call sprsdiag(sa,x)
  if (any(x == 0.0_sp)) call nrerror('asolve_sp: singular diagonal matrix')
  x=b/x
END SUBROUTINE asolve_sp


SUBROUTINE asolve_dp(b,x,itrnsp)
  USE mo_kind, only: dp, i4
  USE mo_nrutil, ONLY : assert_eq,nrerror
  use mo_nr, ONLY : sprsdiag
  USE xlinbcg_data_dp
  REAL(DP), DIMENSION(:), INTENT(IN) :: b
  REAL(DP), DIMENSION(:), INTENT(OUT) :: x
  INTEGER(I4), INTENT(IN) :: itrnsp
  INTEGER(I4) :: ndum
  ndum=assert_eq(size(b),size(x),'asolve_dp')
  call sprsdiag(sa,x)
  if (any(x == 0.0_dp)) call nrerror('asolve_dp: singular diagonal matrix')
  x=b/x
END SUBROUTINE asolve_dp

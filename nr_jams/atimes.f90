SUBROUTINE atimes_sp(x,r,itrnsp)
  USE mo_kind, only: sp, i4
  USE mo_nrutil, ONLY : assert_eq
  use mo_nr, ONLY : sprsax,sprstx
  USE xlinbcg_data_sp
  REAL(SP), DIMENSION(:), INTENT(IN) :: x
  REAL(SP), DIMENSION(:), INTENT(OUT) :: r
  INTEGER(I4), INTENT(IN) :: itrnsp
  INTEGER(I4) :: n
  n=assert_eq(size(x),size(r),'atimes_sp')
  if (itrnsp == 0) then
     call sprsax(sa,x,r)
  else
     call sprstx(sa,x,r)
  end if
END SUBROUTINE atimes_sp


SUBROUTINE atimes_dp(x,r,itrnsp)
  USE mo_kind, only: dp, i4
  USE mo_nrutil, ONLY : assert_eq
  use mo_nr, ONLY : sprsax,sprstx
  USE xlinbcg_data_dp
  REAL(DP), DIMENSION(:), INTENT(IN) :: x
  REAL(DP), DIMENSION(:), INTENT(OUT) :: r
  INTEGER(I4), INTENT(IN) :: itrnsp
  INTEGER(I4) :: n
  n=assert_eq(size(x),size(r),'atimes_dp')
  if (itrnsp == 0) then
     call sprsax(sa,x,r)
  else
     call sprstx(sa,x,r)
  end if
END SUBROUTINE atimes_dp

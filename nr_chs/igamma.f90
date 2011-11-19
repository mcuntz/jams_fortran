FUNCTION igamma_s(a,x)
  USE mo_kind
  USE mo_nrutil, ONLY : assert
  use mo_nr, ONLY : gcf,gser
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,x
  REAL(SP) :: gln
  REAL(SP) :: igamma_s
  call assert( x >= 0.0,  a > 0.0, 'igamma_s args')
  if (x<a+1.0_sp) then
     igamma_s=1.0_sp-gser(a,x,gln)
  else
     igamma_s=gcf(a,x,gln)
  end if
  igamma_s = igamma_s * exp(gln)
END FUNCTION igamma_s


FUNCTION igamma_v(a,x)
  USE mo_kind
  USE mo_nrutil, ONLY : assert, assert_eq
  use mo_nr, ONLY : gcf,gser
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(SP), DIMENSION(size(a)) :: gln
  REAL(SP), DIMENSION(size(a)) :: igamma_v
  LOGICAL(LGT), DIMENSION(size(x)) :: mask
  INTEGER(I4) :: ndum
  ndum=assert_eq(size(a),size(x),'igamma_v')
  call assert( all(x >= 0.0),  all(a > 0.0), 'igamma_v args')
  mask = (x<a+1.0_sp)
  igamma_v=merge(1.0_sp-gser(a,merge(x,0.0_sp,mask),gln), &
       gcf(a,merge(x,0.0_sp,.not. mask),gln),mask)
  igamma_v = igamma_v * exp(gln)
END FUNCTION igamma_v


FUNCTION digamma_s(a,x)
  USE mo_kind
  USE mo_nrutil, ONLY : assert
  use mo_nr, ONLY : gcf,gser
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,x
  REAL(DP) :: gln
  REAL(DP) :: digamma_s
  call assert( x >= 0.0_dp,  a > 0.0_dp, 'digamma_s args')
  if (x<a+1.0_dp) then
     digamma_s=1.0_dp-gser(a,x,gln)
  else
     digamma_s=gcf(a,x,gln)
  end if
  digamma_s = digamma_s * exp(gln)
END FUNCTION digamma_s


FUNCTION digamma_v(a,x)
  USE mo_kind
  USE mo_nrutil, ONLY : assert, assert_eq
  use mo_nr, ONLY : gcf, gser
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(DP), DIMENSION(size(a)) :: gln
  REAL(DP), DIMENSION(size(a)) :: digamma_v
  LOGICAL(LGT), DIMENSION(size(x)) :: mask
  INTEGER(I4) :: ndum
  ndum=assert_eq(size(a),size(x),'digamma_v')
  call assert( all(x >= 0.0_dp),  all(a > 0.0_dp), 'digamma_v args')
  mask = (x<a+1.0_dp)
  digamma_v=merge(1.0_dp-gser(a,merge(x,0.0_dp,mask),gln), &
       gcf(a,merge(x,0.0_dp,.not. mask),gln),mask)
  digamma_v = digamma_v * exp(gln)
END FUNCTION digamma_v

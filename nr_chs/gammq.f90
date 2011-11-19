FUNCTION gammq_s(a,x,gln)
  USE mo_kind
  USE mo_nrutil, ONLY : assert
  use mo_nr, ONLY : gcf,gser
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,x
  REAL(SP), OPTIONAL, INTENT(OUT) :: gln
  REAL(SP) :: gammq_s
  call assert( x >= 0.0,  a > 0.0, 'gammq_s args')
  if (x<a+1.0_sp) then
     if (present(gln)) then
        gammq_s=1.0_sp-gser(a,x,gln)
     else
        gammq_s=1.0_sp-gser(a,x)
     endif
  else
     if (present(gln)) then
        gammq_s=gcf(a,x,gln)
     else
        gammq_s=gcf(a,x)
     endif
  end if
END FUNCTION gammq_s


FUNCTION gammq_v(a,x,gln)
  USE mo_kind
  USE mo_nrutil, ONLY : assert, assert_eq, nrerror
  use mo_nr, ONLY : gcf,gser
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
  REAL(SP), DIMENSION(size(a)) :: gammq_v
  LOGICAL(LGT), DIMENSION(size(x)) :: mask
  INTEGER(I4) :: ndum
  ndum=assert_eq(size(a),size(x),'gammq_v')
  call assert( all(x >= 0.0),  all(a > 0.0), 'gammq_v args')
  mask = (x<a+1.0_sp)
  if (present(gln)) then
     if (size(gln) < size(a)) call nrerror('gammq_v: Not enough space for gln')
     gammq_v=merge(1.0_sp-gser(a,merge(x,0.0_sp,mask),gln), &
          gcf(a,merge(x,0.0_sp,.not. mask),gln),mask)
  else
     gammq_v=merge(1.0_sp-gser(a,merge(x,0.0_sp,mask),gln), &
          gcf(a,merge(x,0.0_sp,.not. mask),gln),mask)
  endif
END FUNCTION gammq_v


FUNCTION dgammq_s(a,x,gln)
  USE mo_kind
  USE mo_nrutil, ONLY : assert
  use mo_nr, ONLY : gcf,gser
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,x
  REAL(DP), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP) :: dgammq_s
  call assert( x >= 0.0_dp,  a > 0.0_dp, 'dgammq_s args')
  if (x<a+1.0_dp) then
     if (present(gln)) then
        dgammq_s=1.0_dp-gser(a,x,gln)
     else
        dgammq_s=1.0_dp-gser(a,x)
     endif
  else
     if (present(gln)) then
        dgammq_s=gcf(a,x,gln)
     else
        dgammq_s=gcf(a,x)
     endif
  end if
END FUNCTION dgammq_s


FUNCTION dgammq_v(a,x,gln)
  USE mo_kind
  USE mo_nrutil, ONLY : assert, assert_eq, nrerror
  use mo_nr, ONLY : gcf, gser
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP), DIMENSION(size(a)) :: dgammq_v
  LOGICAL(LGT), DIMENSION(size(x)) :: mask
  INTEGER(I4) :: ndum
  ndum=assert_eq(size(a),size(x),'dgammq_v')
  call assert( all(x >= 0.0_dp),  all(a > 0.0_dp), 'dgammq_v args')
  mask = (x<a+1.0_dp)
  if (present(gln)) then
     if (size(gln) < size(a)) call nrerror('dgammq_v: Not enough space for gln')
     dgammq_v=merge(1.0_dp-gser(a,merge(x,0.0_dp,mask),gln), &
          gcf(a,merge(x,0.0_dp,.not. mask),gln),mask)
  else
     dgammq_v=merge(1.0_dp-gser(a,merge(x,0.0_dp,mask),gln), &
          gcf(a,merge(x,0.0_dp,.not. mask),gln),mask)
  endif
END FUNCTION dgammq_v

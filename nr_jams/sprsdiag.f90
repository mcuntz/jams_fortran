SUBROUTINE sprsdiag_sp(sa,b)
  USE mo_kind, only: sp, i4, lgt, sprs2_sp
  USE mo_nrutil, ONLY : array_copy,assert_eq
  IMPLICIT NONE
  TYPE(sprs2_sp), INTENT(IN) :: sa
  REAL(SP), DIMENSION(:), INTENT(OUT) :: b
  REAL(SP), DIMENSION(size(b)) :: val
  INTEGER(I4) :: k,l,ndum,nerr
  INTEGER(I4), DIMENSION(size(b)) :: i
  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE :: mask
  ndum=assert_eq(sa%n,size(b),'sprsdiag_sp')
  l=sa%len
  allocate(mask(l))
  mask = (sa%irow(1:l) == sa%jcol(1:l))
  call array_copy(pack(sa%val(1:l),mask),val,k,nerr)
  i(1:k)=pack(sa%irow(1:l),mask)
  deallocate(mask)
  b=0.0
  b(i(1:k))=val(1:k)
END SUBROUTINE sprsdiag_sp

SUBROUTINE sprsdiag_dp(sa,b)
  USE mo_kind, only: dp, i4, lgt, sprs2_dp
  USE mo_nrutil, ONLY : array_copy,assert_eq
  IMPLICIT NONE
  TYPE(sprs2_dp), INTENT(IN) :: sa
  REAL(DP), DIMENSION(:), INTENT(OUT) :: b
  REAL(DP), DIMENSION(size(b)) :: val
  INTEGER(I4) :: k,l,ndum,nerr
  INTEGER(I4), DIMENSION(size(b)) :: i
  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE :: mask
  ndum=assert_eq(sa%n,size(b),'sprsdiag_dp')
  l=sa%len
  allocate(mask(l))
  mask = (sa%irow(1:l) == sa%jcol(1:l))
  call array_copy(pack(sa%val(1:l),mask),val,k,nerr)
  i(1:k)=pack(sa%irow(1:l),mask)
  deallocate(mask)
  b=0.0
  b(i(1:k))=val(1:k)
END SUBROUTINE sprsdiag_dp

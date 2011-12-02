SUBROUTINE choldc_sp(a,p)
  use mo_kind, only: i4, sp
  use mo_nrutil, ONLY : assert_eq,nrerror
  IMPLICIT NONE
  !
  REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
  REAL(SP), DIMENSION(:), INTENT(OUT)     :: p
  INTEGER(I4)                             :: i,n
  REAL(SP)                                :: summ
  !
  n=assert_eq(size(a,1),size(a,2),size(p),'choldc_sp')
  do i=1,n
     summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
     if (summ <= 0.0) call nrerror('choldc_sp failed')
     p(i)=sqrt(summ)
     a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
  end do
END SUBROUTINE choldc_sp

SUBROUTINE choldc_dp(a,p)
  use mo_kind, only: i4, dp
  use mo_nrutil, ONLY : assert_eq,nrerror
  IMPLICIT NONE
  !
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
  REAL(DP), DIMENSION(:), INTENT(OUT)     :: p
  INTEGER(I4)                             :: i,n
  REAL(DP)                                :: summ
  !
  n=assert_eq(size(a,1),size(a,2),size(p),'choldc_dp')
  do i=1,n
     summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
     if (summ <= 0.0) call nrerror('choldc_dp failed')
     p(i)=sqrt(summ)
     a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
  end do
END SUBROUTINE choldc_dp

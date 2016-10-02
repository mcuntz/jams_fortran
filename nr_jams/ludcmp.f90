SUBROUTINE ludcmp_sp(a,indx,d)
  use mo_kind,   ONLY : i4, sp
  use mo_nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
  IMPLICIT NONE
  REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
  INTEGER(I4), DIMENSION(:), INTENT(OUT) :: indx
  REAL(SP), INTENT(OUT) :: d
  REAL(SP), DIMENSION(size(a,1)) :: vv
  REAL(SP), PARAMETER :: TINY=1.0e-20_sp
  INTEGER(I4) :: j,n,imax
  n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
  d=1.0_sp
  vv=maxval(abs(a),dim=2)
  if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
  vv=1.0_sp/vv
  do j=1,n
     imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
     if (j /= imax) then
        call swap(a(imax,:),a(j,:))
        d=-d
        vv(imax)=vv(j)
     end if
     indx(j)=imax
     if (a(j,j) == 0.0) a(j,j)=TINY
     a(j+1:n,j)=a(j+1:n,j)/a(j,j)
     a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
  end do
END SUBROUTINE ludcmp_sp

SUBROUTINE ludcmp_dp(a,indx,d)
  use mo_kind,   ONLY : i4, dp
  use mo_nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
  INTEGER(I4), DIMENSION(:), INTENT(OUT) :: indx
  REAL(DP), INTENT(OUT) :: d
  REAL(DP), DIMENSION(size(a,1)) :: vv
  REAL(DP), PARAMETER :: TINY=1.0e-20_dp
  INTEGER(I4) :: j,n,imax
  n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
  d=1.0_dp
  vv=maxval(abs(a),dim=2)
  if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
  vv=1.0_dp/vv
  do j=1,n
     imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
     if (j /= imax) then
        call swap(a(imax,:),a(j,:))
        d=-d
        vv(imax)=vv(j)
     end if
     indx(j)=imax
     if (a(j,j) == 0.0) a(j,j)=TINY
     a(j+1:n,j)=a(j+1:n,j)/a(j,j)
     a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
  end do
END SUBROUTINE ludcmp_dp

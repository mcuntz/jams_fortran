SUBROUTINE tridag_ser(a,b,c,r,u)
  USE mo_kind, only: sp, i4
  USE mo_nrutil, ONLY : assert_eq,nrerror
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
  REAL(SP), DIMENSION(:), INTENT(OUT) :: u
  REAL(SP), DIMENSION(size(b)) :: gam
  INTEGER(I4) :: n,j
  REAL(SP) :: bet
  n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
  bet=b(1)
  if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
  u(1)=r(1)/bet
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j-1)*gam(j)
     if (bet == 0.0) &
          call nrerror('tridag_ser: Error at code stage 2')
     u(j)=(r(j)-a(j-1)*u(j-1))/bet
  end do
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  end do
END SUBROUTINE tridag_ser

RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
  USE mo_kind, only: sp, i4
  USE mo_nrutil, ONLY : assert_eq,nrerror
  use mo_nr, ONLY : tridag_ser
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
  REAL(SP), DIMENSION(:), INTENT(OUT) :: u
  INTEGER(I4), PARAMETER :: NPAR_TRIDAG=4
  INTEGER(I4) :: n,n2,nm,nx
  REAL(SP), DIMENSION(size(b)/2) :: y,q,piva
  REAL(SP), DIMENSION(size(b)/2-1) :: x,z
  REAL(SP), DIMENSION(size(a)/2) :: pivc
  n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_par')
  if (n < NPAR_TRIDAG) then
     call tridag_ser(a,b,c,r,u)
  else
     if (maxval(abs(b(1:n))) == 0.0) &
          call nrerror('tridag_par: possible singular matrix')
     n2=size(y)
     nm=size(pivc)
     nx=size(x)
     piva = a(1:n-1:2)/b(1:n-1:2)
     pivc = c(2:n-1:2)/b(3:n:2)
     y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
     q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
     if (nm < n2) then
        y(n2) = b(n)-piva(n2)*c(n-1)
        q(n2) = r(n)-piva(n2)*r(n-1)
     end if
     x = -piva(2:n2)*a(2:n-2:2)
     z = -pivc(1:nx)*c(3:n-1:2)
     call tridag_par(x,y,z,q,u(2:n:2))
     u(1) = (r(1)-c(1)*u(2))/b(1)
     u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
          -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
     if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
  end if
END SUBROUTINE tridag_par


SUBROUTINE dtridag_ser(a,b,c,r,u)
  USE mo_kind, only: dp, i4
  USE mo_nrutil, ONLY : assert_eq,nrerror
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
  REAL(DP), DIMENSION(:), INTENT(OUT) :: u
  REAL(DP), DIMENSION(size(b)) :: gam
  INTEGER(I4) :: n,j
  REAL(DP) :: bet
  n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'dtridag_ser')
  bet=b(1)
  if (bet == 0.0_dp) call nrerror('dtridag_ser: Error at code stage 1')
  u(1)=r(1)/bet
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j-1)*gam(j)
     if (bet == 0.0_dp) &
          call nrerror('tridag_ser: Error at code stage 2')
     u(j)=(r(j)-a(j-1)*u(j-1))/bet
  end do
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  end do
END SUBROUTINE dtridag_ser


RECURSIVE SUBROUTINE dtridag_par(a,b,c,r,u)
  USE mo_kind, only: dp, i4
  USE mo_nrutil, ONLY : assert_eq,nrerror
  use mo_nr, ONLY : dtridag_ser
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
  REAL(DP), DIMENSION(:), INTENT(OUT) :: u
  INTEGER(I4), PARAMETER :: NPAR_TRIDAG=4
  INTEGER(I4) :: n,n2,nm,nx
  REAL(DP), DIMENSION(size(b)/2) :: y,q,piva
  REAL(DP), DIMENSION(size(b)/2-1) :: x,z
  REAL(DP), DIMENSION(size(a)/2) :: pivc
  n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'dtridag_par')
  if (n < NPAR_TRIDAG) then
     call dtridag_ser(a,b,c,r,u)
  else
     if (maxval(abs(b(1:n))) == 0.0_dp) &
          call nrerror('dtridag_par: possible singular matrix')
     n2=size(y)
     nm=size(pivc)
     nx=size(x)
     piva = a(1:n-1:2)/b(1:n-1:2)
     pivc = c(2:n-1:2)/b(3:n:2)
     y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
     q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
     if (nm < n2) then
        y(n2) = b(n)-piva(n2)*c(n-1)
        q(n2) = r(n)-piva(n2)*r(n-1)
     end if
     x = -piva(2:n2)*a(2:n-2:2)
     z = -pivc(1:nx)*c(3:n-1:2)
     call dtridag_par(x,y,z,q,u(2:n:2))
     u(1) = (r(1)-c(1)*u(2))/b(1)
     u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
          -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
     if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
  end if
END SUBROUTINE dtridag_par

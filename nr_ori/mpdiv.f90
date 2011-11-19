	SUBROUTINE mpdiv(q,r,u,v,n,m)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : mpinv,mpmul
	USE mpops, ONLY : mpmov,mpsub
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: q,r
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
!	The logical dimensions are: CHARACTER(1) :: q(n-m+1),r(m),u(n),v(m)
	INTEGER(I4B), INTENT(IN) :: n,m
	INTEGER(I4B), PARAMETER :: MACC=3
	INTEGER(I4B) :: is
	CHARACTER(1), DIMENSION(:), ALLOCATABLE :: s
	CHARACTER(1), DIMENSION(:), ALLOCATABLE, TARGET :: rr
	CHARACTER(1), DIMENSION(:), POINTER :: rr2
	allocate(rr(2*(n+MACC)),s(n+MACC))
	rr2=>rr(2:)
	call mpinv(s,v,n-m+MACC,m)
	call mpmul(rr,s,u,n-m+MACC,n)
	call mpmov(q,rr2,n-m+1)
	call mpmul(rr,q,v,n-m+1,m)
	call mpsub(is,rr2,u,rr2,n)
	if (is /= 0) call nrerror('MACC too small in mpdiv')
	call mpmov(r,rr(n-m+2:),m)
	deallocate(rr,s)
	END SUBROUTINE mpdiv

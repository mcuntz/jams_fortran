MODULE f1dim_mod
	use mo_nrtype
	INTEGER(I4B) :: ncom
	REAL(SP), DIMENSION(:), POINTER :: pcom,xicom
CONTAINS
!BL
	FUNCTION f1dim(x)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: f1dim
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), DIMENSION(:), ALLOCATABLE :: xt
	allocate(xt(ncom))
	xt(:)=pcom(:)+x*xicom(:)
	f1dim=func(xt)
	deallocate(xt)
	END FUNCTION f1dim
END MODULE f1dim_mod

	SUBROUTINE linmin(p,xi,fret)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : mnbrak,brent
	USE f1dim_mod
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: fret
	REAL(SP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
	REAL(SP), PARAMETER :: TOL=1.0e-4_sp
	REAL(SP) :: ax,bx,fa,fb,fx,xmin,xx
	ncom=assert_eq(size(p),size(xi),'linmin')
	pcom=>p
	xicom=>xi
	ax=0.0
	xx=1.0
	call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
	fret=brent(ax,xx,bx,f1dim,TOL,xmin)
	xi=xmin*xi
	p=p+xi
	END SUBROUTINE linmin

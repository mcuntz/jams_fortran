!	FUNCTION shoot(v) is named "funcv" for use with "newt"
	FUNCTION funcv(v)
	use mo_nrtype
	use mo_nr, ONLY : odeint,rkqs
	USE sphoot_caller, ONLY: nvar,x1,x2
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: v
	REAL(SP), DIMENSION(size(v)) :: funcv
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	REAL(SP) :: h1,hmin
	REAL(SP), DIMENSION(nvar) :: y
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
!BL
		SUBROUTINE load(x1,v,y)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x1
		REAL(SP), DIMENSION(:), INTENT(IN) :: v
		REAL(SP), DIMENSION(:), INTENT(OUT) :: y
		END SUBROUTINE load
!BL
		SUBROUTINE score(x2,y,f)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x2
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: f
		END SUBROUTINE score
	END INTERFACE
	h1=(x2-x1)/100.0_sp
	hmin=0.0
	call load(x1,v,y)
	call odeint(y,x1,x2,EPS,h1,hmin,derivs,rkqs)
	call score(x2,y,funcv)
	END FUNCTION funcv

	FUNCTION qsimp(func,a,b)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : trapzd
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qsimp
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp, UNLIKELY=-1.0e30_sp
	INTEGER(I4B) :: j
	REAL(SP) :: os,ost,st
	ost=UNLIKELY
	os= UNLIKELY
	do j=1,JMAX
		call trapzd(func,a,b,st,j)
		qsimp=(4.0_sp*st-ost)/3.0_sp
		if (abs(qsimp-os) < EPS*abs(os)) RETURN
		if (qsimp == 0.0 .and. os == 0.0 .and. j > 6) RETURN
		os=qsimp
		ost=st
	end do
	call nrerror('qsimp: too many steps')
	END FUNCTION qsimp

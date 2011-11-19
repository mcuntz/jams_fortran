	SUBROUTINE gaucof(a,b,amu0,x,w)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,unit_matrix
	use mo_nr, ONLY : eigsrt,tqli
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
	REAL(SP), INTENT(IN) :: amu0
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(SP), DIMENSION(size(a),size(a)) :: z
	INTEGER(I4B) :: n
	n=assert_eq(size(a),size(b),size(x),size(w),'gaucof')
	b(2:n)=sqrt(b(2:n))
	call unit_matrix(z)
	call tqli(a,b,z)
	call eigsrt(a,z)
	x=a
	w=amu0*z(1,:)**2
	END SUBROUTINE gaucof

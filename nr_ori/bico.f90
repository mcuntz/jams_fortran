	FUNCTION bico_s(n,k)
	use mo_nrtype
	use mo_nr, ONLY : factln
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n,k
	REAL(SP) :: bico_s
	bico_s=nint(exp(factln(n)-factln(k)-factln(n-k)))
	END FUNCTION bico_s


	FUNCTION bico_v(n,k)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : factln
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n,k
	REAL(SP), DIMENSION(size(n)) :: bico_v
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(n),size(k),'bico_v')
	bico_v=nint(exp(factln(n)-factln(k)-factln(n-k)))
	END FUNCTION bico_v

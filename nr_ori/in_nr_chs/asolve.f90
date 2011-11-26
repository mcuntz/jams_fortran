	SUBROUTINE asolve(b,x,itrnsp)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : sprsdiag
	USE xlinbcg_data
	REAL(DP), DIMENSION(:), INTENT(IN) :: b
	REAL(DP), DIMENSION(:), INTENT(OUT) :: x
	INTEGER(I4B), INTENT(IN) :: itrnsp
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(b),size(x),'asolve')
	call sprsdiag(sa,x)
	if (any(x == 0.0)) call nrerror('asolve: singular diagonal matrix')
	x=b/x
	END SUBROUTINE asolve

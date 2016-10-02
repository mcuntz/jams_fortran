	SUBROUTINE strapzd(func,a,b,s,n)
	use mo_kind
        use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), INTENT(INOUT) :: s
	INTEGER(I4), INTENT(IN) :: n
	INTERFACE
		FUNCTION func(x)
		use mo_kind
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP) :: del,fsum
	INTEGER(I4) :: it
	if (n == 1) then
		s=0.5_sp*(b-a)*sum(func( (/ a,b /) ))
	else
		it=2**(n-2)
		del=(b-a)/it
		fsum=sum(func(arth(a+0.5_sp*del,del,it)))
		s=0.5_sp*(s+del*fsum)
	end if
	END SUBROUTINE strapzd

	SUBROUTINE dtrapzd(func,a,b,s,n)
	use mo_kind
        use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,b
	REAL(DP), INTENT(INOUT) :: s
	INTEGER(I4), INTENT(IN) :: n
	INTERFACE
		FUNCTION func(x)
		use mo_kind
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(DP) :: del,fsum
	INTEGER(I4) :: it
	if (n == 1) then
		s=0.5_dp*(b-a)*sum(func( (/ a,b /) ))
	else
		it=2**(n-2)
		del=(b-a)/it
		fsum=sum(func(arth(a+0.5_dp*del,del,it)))
		s=0.5_dp*(s+del*fsum)
	end if
	END SUBROUTINE dtrapzd


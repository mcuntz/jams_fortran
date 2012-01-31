	SUBROUTINE trapzd(func,a,b,s,n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP) :: del,fsum
	INTEGER(I4B) :: it
	if (n == 1) then
		s=0.5_sp*(b-a)*sum(func( (/ a,b /) ))
	else
		it=2**(n-2)
		del=(b-a)/it
		fsum=sum(func(arth(a+0.5_sp*del,del,it)))
		s=0.5_sp*(s+del*fsum)
	end if
	END SUBROUTINE trapzd

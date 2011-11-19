	SUBROUTINE sort2(arr,slave)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : indexx
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave
	INTEGER(I4B) :: ndum
	INTEGER(I4B), DIMENSION(size(arr)) :: index
	ndum=assert_eq(size(arr),size(slave),'sort2')
	call indexx(arr,index)
	arr=arr(index)
	slave=slave(index)
	END SUBROUTINE sort2

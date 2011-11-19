	SUBROUTINE four2(data,isign)
	use mo_nrtype
	use mo_nr, ONLY : fourrow
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(size(data,2),size(data,1)) :: temp
	call fourrow(data,isign)
	temp=transpose(data)
	call fourrow(temp,isign)
	data=transpose(temp)
	END SUBROUTINE four2

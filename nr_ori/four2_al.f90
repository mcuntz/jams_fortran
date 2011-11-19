	SUBROUTINE four2_alt(data,isign)
	use mo_nrtype
	use mo_nr, ONLY : fourcol
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(size(data,2),size(data,1)) :: temp
	temp=transpose(data)
	call fourcol(temp,isign)
	data=transpose(temp)
	call fourcol(data,isign)
	END SUBROUTINE four2_alt

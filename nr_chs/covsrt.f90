	SUBROUTINE covsrt_sp(covar,maska)
	use mo_kind
        use mo_nrutil, ONLY : assert_eq,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
	INTEGER(I4) :: ma,mfit,j,k
	ma=assert_eq(size(covar,1),size(covar,2),size(maska),'covsrt')
	mfit=count(maska)
	covar(mfit+1:ma,1:ma)=0.0_sp
	covar(1:ma,mfit+1:ma)=0.0_sp
	k=mfit
	do j=ma,1,-1
		if (maska(j)) then
			call swap(covar(1:ma,k),covar(1:ma,j))
			call swap(covar(k,1:ma),covar(j,1:ma))
			k=k-1
		end if
	end do
	END SUBROUTINE covsrt_sp

	SUBROUTINE covsrt_dp(covar,maska)
	use mo_kind
        use mo_nrutil, ONLY : assert_eq,swap
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: covar
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
	INTEGER(I4) :: ma,mfit,j,k
	ma=assert_eq(size(covar,1),size(covar,2),size(maska),'covsrt')
	mfit=count(maska)
	covar(mfit+1:ma,1:ma)=0.0_dp
	covar(1:ma,mfit+1:ma)=0.0_dp
	k=mfit
	do j=ma,1,-1
		if (maska(j)) then
			call swap(covar(1:ma,k),covar(1:ma,j))
			call swap(covar(k,1:ma),covar(j,1:ma))
			k=k-1
		end if
	end do
	END SUBROUTINE covsrt_dp

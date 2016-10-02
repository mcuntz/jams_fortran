	FUNCTION bessi_s(n,x)
	use mo_kind
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessi0
	IMPLICIT NONE
	INTEGER(I4), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessi_s
	INTEGER(I4), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
	INTEGER(I4) :: j,m
	REAL(SP) :: bi,bim,bip,tox
	call assert(n >= 2, 'bessi_s args')
	bessi_s=0.0
	if (x*x <= 8.0_sp*tiny(x)) RETURN
	tox=2.0_sp/abs(x)
	bip=0.0
	bi=1.0
	m=2*((n+int(sqrt(real(IACC*n,sp)))))
	do j=m,1,-1
		bim=bip+j*tox*bi
		bip=bi
		bi=bim
		if (exponent(bi) > IEXP) then
			bessi_s=scale(bessi_s,-IEXP)
			bi=scale(bi,-IEXP)
			bip=scale(bip,-IEXP)
		end if
		if (j == n) bessi_s=bip
	end do
	bessi_s=bessi_s*bessi0(x)/bi
	if (x < 0.0 .and. mod(n,2) == 1) bessi_s=-bessi_s
	END FUNCTION bessi_s


	FUNCTION bessi_v(n,x)
	use mo_kind
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessi0
	IMPLICIT NONE
	INTEGER(I4), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessi_v
	INTEGER(I4), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
	INTEGER(I4) :: j,m
	REAL(SP), DIMENSION(size(x)) :: bi,bim,bip,tox
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	call assert(n >= 2, 'bessi_v args')
	bessi_v=0.0
	mask = (x <= 8.0_sp*tiny(x))
	tox=2.0_sp/merge(2.0_sp,abs(x),mask)
	bip=0.0
	bi=1.0_sp
	m=2*((n+int(sqrt(real(IACC*n,sp)))))
	do j=m,1,-1
		bim=bip+j*tox*bi
		bip=bi
		bi=bim
		where (exponent(bi) > IEXP)
			bessi_v=scale(bessi_v,-IEXP)
			bi=scale(bi,-IEXP)
			bip=scale(bip,-IEXP)
		end where
		if (j == n) bessi_v=bip
	end do
	bessi_v=bessi_v*bessi0(x)/bi
	where (mask) bessi_v=0.0_sp
	where (x < 0.0 .and. mod(n,2) == 1) bessi_v=-bessi_v
	END FUNCTION bessi_v


	FUNCTION dbessi_s(n,x)
	use mo_kind
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessi0
	IMPLICIT NONE
	INTEGER(I4), INTENT(IN) :: n
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: dbessi_s
	INTEGER(I4), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
	INTEGER(I4) :: j,m
	REAL(DP) :: bi,bim,bip,tox
	call assert(n >= 2, 'dbessi_s args')
	dbessi_s=0.0
	if (x*x <= 8.0_dp*tiny(x)) RETURN
	tox=2.0_dp/abs(x)
	bip=0.0
	bi=1.0
	m=2*((n+int(sqrt(real(IACC*n,dp)))))
	do j=m,1,-1
		bim=bip+j*tox*bi
		bip=bi
		bi=bim
		if (exponent(bi) > IEXP) then
			dbessi_s=scale(dbessi_s,-IEXP)
			bi=scale(bi,-IEXP)
			bip=scale(bip,-IEXP)
		end if
		if (j == n) dbessi_s=bip
	end do
	dbessi_s=dbessi_s*bessi0(x)/bi
	if (x < 0.0 .and. mod(n,2) == 1) dbessi_s=-dbessi_s
	END FUNCTION dbessi_s


	FUNCTION dbessi_v(n,x)
	use mo_kind
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessi0
	IMPLICIT NONE
	INTEGER(I4), INTENT(IN) :: n
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(size(x)) :: dbessi_v
	INTEGER(I4), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
	INTEGER(I4) :: j,m
	REAL(DP), DIMENSION(size(x)) :: bi,bim,bip,tox
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	call assert(n >= 2, 'dbessi_v args')
	dbessi_v=0.0
	mask = (x <= 8.0_dp*tiny(x))
	tox=2.0_dp/merge(2.0_dp,abs(x),mask)
	bip=0.0
	bi=1.0_dp
	m=2*((n+int(sqrt(real(IACC*n,dp)))))
	do j=m,1,-1
		bim=bip+j*tox*bi
		bip=bi
		bi=bim
		where (exponent(bi) > IEXP)
			dbessi_v=scale(dbessi_v,-IEXP)
			bi=scale(bi,-IEXP)
			bip=scale(bip,-IEXP)
		end where
		if (j == n) dbessi_v=bip
	end do
	dbessi_v=dbessi_v*bessi0(x)/bi
	where (mask) dbessi_v=0.0_dp
	where (x < 0.0 .and. mod(n,2) == 1) dbessi_v=-dbessi_v
	END FUNCTION dbessi_v


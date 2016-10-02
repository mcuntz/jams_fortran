	FUNCTION bessk_s(n,x)
	use mo_kind
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessk0,bessk1
	IMPLICIT NONE
	INTEGER(I4), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessk_s
	INTEGER(I4) :: j
	REAL(SP) :: bk,bkm,bkp,tox
	call assert(n >= 2, x > 0.0, 'bessk_s args')
	tox=2.0_sp/x
	bkm=bessk0(x)
	bk=bessk1(x)
	do j=1,n-1
		bkp=bkm+j*tox*bk
		bkm=bk
		bk=bkp
	end do
	bessk_s=bk
	END FUNCTION bessk_s


	FUNCTION bessk_v(n,x)
	use mo_kind
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessk0,bessk1
	IMPLICIT NONE
	INTEGER(I4), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessk_v
	INTEGER(I4) :: j
	REAL(SP), DIMENSION(size(x)) :: bk,bkm,bkp,tox
	call assert(n >= 2, all(x > 0.0), 'bessk_v args')
	tox=2.0_sp/x
	bkm=bessk0(x)
	bk=bessk1(x)
	do j=1,n-1
		bkp=bkm+j*tox*bk
		bkm=bk
		bk=bkp
	end do
	bessk_v=bk
	END FUNCTION bessk_v

	FUNCTION dbessk_s(n,x)
	use mo_kind
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessk0,bessk1
	IMPLICIT NONE
	INTEGER(I4), INTENT(IN) :: n
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: dbessk_s
	INTEGER(I4) :: j
	REAL(DP) :: bk,bkm,bkp,tox
	call assert(n >= 2, x > 0.0, 'dbessk_s args')
	tox=2.0_dp/x
	bkm=bessk0(x)
	bk=bessk1(x)
	do j=1,n-1
		bkp=bkm+j*tox*bk
		bkm=bk
		bk=bkp
	end do
	dbessk_s=bk
	END FUNCTION dbessk_s


	FUNCTION dbessk_v(n,x)
	use mo_kind
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessk0,bessk1
	IMPLICIT NONE
	INTEGER(I4), INTENT(IN) :: n
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(size(x)) :: dbessk_v
	INTEGER(I4) :: j
	REAL(DP), DIMENSION(size(x)) :: bk,bkm,bkp,tox
	call assert(n >= 2, all(x > 0.0), 'dbessk_v args')
	tox=2.0_dp/x
	bkm=bessk0(x)
	bk=bessk1(x)
	do j=1,n-1
		bkp=bkm+j*tox*bk
		bkm=bk
		bk=bkp
	end do
	dbessk_v=bk
	END FUNCTION dbessk_v


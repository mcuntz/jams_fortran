	FUNCTION ellf_s(phi,ak)
	use mo_nrtype
	use mo_nr, ONLY : rf
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: phi,ak
	REAL(SP) :: ellf_s
	REAL(SP) :: s
	s=sin(phi)
	ellf_s=s*rf(cos(phi)**2,(1.0_sp-s*ak)*(1.0_sp+s*ak),1.0_sp)
	END FUNCTION ellf_s


	FUNCTION ellf_v(phi,ak)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : rf
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: phi,ak
	REAL(SP), DIMENSION(size(phi)) :: ellf_v
	REAL(SP), DIMENSION(size(phi)) :: s
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(phi),size(ak),'ellf_v')
	s=sin(phi)
	ellf_v=s*rf(cos(phi)**2,(1.0_sp-s*ak)*(1.0_sp+s*ak),&
		spread(1.0_sp,1,size(phi)))
	END FUNCTION ellf_v

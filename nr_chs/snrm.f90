FUNCTION snrm_sp(sx,itol)
  USE mo_kind
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: sx
  INTEGER(I4), INTENT(IN) :: itol
  REAL(SP) :: snrm_sp
  if (itol <= 3) then
     snrm_sp=sqrt(dot_product(sx,sx))
  else
     snrm_sp=maxval(abs(sx))
  end if
END FUNCTION snrm_sp


FUNCTION snrm_dp(sx,itol)
  USE mo_kind
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: sx
  INTEGER(I4), INTENT(IN) :: itol
  REAL(DP) :: snrm_dp
  if (itol <= 3) then
     snrm_dp=sqrt(dot_product(sx,sx))
  else
     snrm_dp=maxval(abs(sx))
  end if
END FUNCTION snrm_dp

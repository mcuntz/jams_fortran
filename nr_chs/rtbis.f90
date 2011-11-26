FUNCTION rtbis_s(func,x1,x2,xacc)
  USE mo_kind, only: sp, i4
  USE mo_nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: x1,x2,xacc
  REAL(SP) :: rtbis_s
  INTERFACE
     FUNCTION func(x)
       USE mo_kind
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: func
     END FUNCTION func
  END INTERFACE
  INTEGER(I4), PARAMETER :: MAXIT=40
  INTEGER(I4) :: j
  REAL(SP) :: dx,f,fmid,xmid
  fmid=func(x2)
  f=func(x1)
  if (f*fmid >= 0.0_sp) call nrerror('rtbis_s: root must be bracketed')
  if (f < 0.0_sp) then
     rtbis_s=x1
     dx=x2-x1
  else
     rtbis_s=x2
     dx=x1-x2
  end if
  do j=1,MAXIT
     dx=dx*0.5_sp
     xmid=rtbis_s+dx
     fmid=func(xmid)
     if (fmid <= 0.0_sp) rtbis_s=xmid
     if (abs(dx) < xacc .or. fmid == 0.0_sp) RETURN
  end do
  call nrerror('rtbis_s: too many bisections')
END FUNCTION rtbis_s

FUNCTION rtbis_d(func,x1,x2,xacc)
  USE mo_kind, only: dp, i4
  USE mo_nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x1,x2,xacc
  REAL(DP) :: rtbis_d
  INTERFACE
     FUNCTION func(x)
       USE mo_kind
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: func
     END FUNCTION func
  END INTERFACE
  INTEGER(I4), PARAMETER :: MAXIT=40
  INTEGER(I4) :: j
  REAL(DP) :: dx,f,fmid,xmid
  fmid=func(x2)
  f=func(x1)
  if (f*fmid >= 0.0_dp) call nrerror('rtbis_d: root must be bracketed')
  if (f < 0.0_dp) then
     rtbis_d=x1
     dx=x2-x1
  else
     rtbis_d=x2
     dx=x1-x2
  end if
  do j=1, MAXIT
     dx=dx*0.5_dp
     xmid=rtbis_d+dx
     fmid=func(xmid)
     if (fmid <= 0.0_dp) rtbis_d=xmid
     if (abs(dx) < xacc .or. fmid == 0.0_dp) RETURN
  end do
  call nrerror('rtbis_d: too many bisections')
END FUNCTION rtbis_d

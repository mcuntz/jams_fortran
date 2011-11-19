FUNCTION gammln_s(xx)
  USE mo_kind
  USE mo_nrutil, ONLY : arth,assert
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: xx
  REAL(SP) :: gammln_s
  REAL(DP) :: tmp,x
  REAL(DP) :: stp = 2.5066282746310005_dp
  REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
       -86.50532032941677_dp,24.01409824083091_dp,&
       -1.231739572450155_dp,0.1208650973866179e-2_dp,&
       -0.5395239384953e-5_dp/)
  call assert(xx > 0.0, 'gammln_s arg')
  x=xx
  tmp=x+5.5_dp
  tmp=(x+0.5_dp)*log(tmp)-tmp
  gammln_s=tmp+log(stp*(1.000000000190015_dp+&
       sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
END FUNCTION gammln_s


FUNCTION gammln_v(xx)
  USE mo_kind
  USE mo_nrutil, ONLY: assert
  IMPLICIT NONE
  INTEGER(I4) :: i
  REAL(SP), DIMENSION(:), INTENT(IN) :: xx
  REAL(SP), DIMENSION(size(xx)) :: gammln_v
  REAL(DP), DIMENSION(size(xx)) :: ser,tmp,x,y
  REAL(DP) :: stp = 2.5066282746310005_dp
  REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
       -86.50532032941677_dp,24.01409824083091_dp,&
       -1.231739572450155_dp,0.1208650973866179e-2_dp,&
       -0.5395239384953e-5_dp/)
  if (size(xx) == 0) RETURN
  call assert(all(xx > 0.0), 'gammln_v arg')
  x=xx
  tmp=x+5.5_dp
  tmp=(x+0.5_dp)*log(tmp)-tmp
  ser=1.000000000190015_dp
  y=x
  do i=1,size(coef)
     y=y+1.0_dp
     ser=ser+coef(i)/y
  end do
  gammln_v=tmp+log(stp*ser/x)
END FUNCTION gammln_v


FUNCTION dgammln_s(z)
  !  Uses Lanczos-type approximation to ln(gamma) for z > 0.
  !  Reference:
  !       Lanczos, C. ''A precision approximation of the gamma
  !               function'', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
  !  Accuracy: About 14 significant digits except for small regions
  !            in the vicinity of 1 and 2.
  !  Programmer: Alan Miller
  !              1 Creswick Street, Brighton, Vic. 3187, Australia
  !  e-mail: amiller @ bigpond.net.au
  !  Latest revision - 14 October 1996
  USE mo_kind
  USE mo_nrutil, ONLY : assert
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: z
  REAL(DP) :: dgammln_s
  ! Local variables
  REAL(DP)  :: a(9) = (/ 0.9999999999995183_dp, 676.5203681218835_dp, &
       -1259.139216722289_dp, 771.3234287757674_dp, &
       -176.6150291498386_dp, 12.50734324009056_dp, &
       -0.1385710331296526_dp, 0.9934937113930748D-05, &
       0.1659470187408462D-06 /)
  REAL(DP)  :: zero = 0.0_dp,   &
       one = 1.0_dp, &
       lnsqrt2pi =  0.9189385332046727_dp, &
       half = 0.5_dp, &
       sixpt5 = 6.5_dp, &
       seven = 7.0_dp
  REAL(DP)  :: tmp
  INTEGER   :: j
  call assert(z > zero, 'dgammln_s arg')
  dgammln_s = zero
  tmp = z + seven
  DO j = 9, 2, -1
     dgammln_s = dgammln_s + a(j)/tmp
     tmp = tmp - one
  END DO
  dgammln_s = dgammln_s + a(1)
  dgammln_s = LOG(dgammln_s) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
  RETURN
END FUNCTION dgammln_s


FUNCTION dgammln_v(z)
  !  Uses Lanczos-type approximation to ln(gamma) for z > 0.
  !  Reference:
  !       Lanczos, C. ''A precision approximation of the gamma
  !               function'', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
  !  Accuracy: About 14 significant digits except for small regions
  !            in the vicinity of 1 and 2.
  !  Programmer: Alan Miller
  !              1 Creswick Street, Brighton, Vic. 3187, Australia
  !  e-mail: amiller @ bigpond.net.au
  !  Latest revision - 14 October 1996
  USE mo_kind
  USE mo_nrutil, ONLY : assert
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: z
  REAL(DP), DIMENSION(size(z)) :: dgammln_v
  ! Local variables
  REAL(DP)  :: a(9) = (/ 0.9999999999995183_dp, 676.5203681218835_dp, &
       -1259.139216722289_dp, 771.3234287757674_dp, &
       -176.6150291498386_dp, 12.50734324009056_dp, &
       -0.1385710331296526_dp, 0.9934937113930748D-05, &
       0.1659470187408462D-06 /)
  REAL(DP)  :: zero = 0.0_dp,   &
       one = 1.0_dp, &
       lnsqrt2pi =  0.9189385332046727_dp, &
       half = 0.5_dp, &
       sixpt5 = 6.5_dp, &
       seven = 7.0_dp
  REAL(DP), DIMENSION(size(z)) :: tmp
  INTEGER   :: j
  if (size(z) == 0) RETURN
  call assert(all(z > 0.0_dp), 'dgammln_v arg')
  dgammln_v = zero
  tmp = z + seven
  DO j = 9, 2, -1
     dgammln_v = dgammln_v + a(j)/tmp
     tmp = tmp - one
  END DO
  dgammln_v = dgammln_v  + a(1)
  dgammln_v = LOG(dgammln_v) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
END FUNCTION dgammln_v

FUNCTION bessk0_s(x)
  use mo_kind
  use mo_nrutil, ONLY : assert,poly
  use mo_nr, ONLY : bessi0
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: x
  REAL(SP) :: bessk0_s
  REAL(sp) :: y
  REAL(sp), DIMENSION(7) :: p = (/-0.57721566_sp,0.42278420_sp,&
       0.23069756_sp,0.3488590e-1_sp,0.262698e-2_sp,0.10750e-3_sp,&
       0.74e-5_sp/)
  REAL(sp), DIMENSION(7) :: q = (/1.25331414_sp,-0.7832358e-1_sp,&
       0.2189568e-1_sp,-0.1062446e-1_sp,0.587872e-2_sp,&
       -0.251540e-2_sp,0.53208e-3_sp/)
  call assert(x > 0.0, 'bessk0_s arg')
  if (x <= 2.0) then
     y=x*x/4.0_sp
     bessk0_s=(-log(x/2.0_sp)*bessi0(x))+poly(y,p)
  else
     y=(2.0_sp/x)
     bessk0_s=(exp(-x)/sqrt(x))*poly(y,q)
  end if
END FUNCTION bessk0_s


FUNCTION bessk0_v(x)
  use mo_kind
  use mo_nrutil, ONLY : assert,poly
  use mo_nr, ONLY : bessi0
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: x
  REAL(SP), DIMENSION(size(x)) :: bessk0_v
  REAL(sp), DIMENSION(size(x)) :: y
  LOGICAL(LGT), DIMENSION(size(x)) :: mask
  REAL(sp), DIMENSION(7) :: p = (/-0.57721566_sp,0.42278420_sp,&
       0.23069756_sp,0.3488590e-1_sp,0.262698e-2_sp,0.10750e-3_sp,&
       0.74e-5_sp/)
  REAL(sp), DIMENSION(7) :: q = (/1.25331414_sp,-0.7832358e-1_sp,&
       0.2189568e-1_sp,-0.1062446e-1_sp,0.587872e-2_sp,&
       -0.251540e-2_sp,0.53208e-3_sp/)
  call assert(all(x > 0.0), 'bessk0_v arg')
  mask = (x <= 2.0)
  where (mask)
     y=x*x/4.0_sp
     bessk0_v=(-log(x/2.0_sp)*bessi0(x))+poly(y,p,mask)
  elsewhere
     y=(2.0_sp/x)
     bessk0_v=(exp(-x)/sqrt(x))*poly(y,q,.not. mask)
  end where
END FUNCTION bessk0_v

FUNCTION dbessk0_s(x)
  use mo_kind
  use mo_nrutil, ONLY : assert,poly
  use mo_nr, ONLY : bessi0
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x
  REAL(DP) :: dbessk0_s
  REAL(DP) :: y
  REAL(DP), DIMENSION(7) :: p = (/-0.57721566_dp,0.42278420_dp,&
       0.23069756_dp,0.3488590e-1_dp,0.262698e-2_dp,0.10750e-3_dp,&
       0.74e-5_dp/)
  REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,-0.7832358e-1_dp,&
       0.2189568e-1_dp,-0.1062446e-1_dp,0.587872e-2_dp,&
       -0.251540e-2_dp,0.53208e-3_dp/)
  call assert(x > 0.0, 'dbessk0_s arg')
  if (x <= 2.0) then
     y=x*x/4.0_dp
     dbessk0_s=(-log(x/2.0_dp)*bessi0(x))+poly(y,p)
  else
     y=(2.0_dp/x)
     dbessk0_s=(exp(-x)/sqrt(x))*poly(y,q)
  end if
END FUNCTION dbessk0_s


FUNCTION dbessk0_v(x)
  use mo_kind
  use mo_nrutil, ONLY : assert,poly
  use mo_nr, ONLY : bessi0
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: x
  REAL(DP), DIMENSION(size(x)) :: dbessk0_v
  REAL(DP), DIMENSION(size(x)) :: y
  LOGICAL(LGT), DIMENSION(size(x)) :: mask
  REAL(DP), DIMENSION(7) :: p = (/-0.57721566_dp,0.42278420_dp,&
       0.23069756_dp,0.3488590e-1_dp,0.262698e-2_dp,0.10750e-3_dp,&
       0.74e-5_dp/)
  REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,-0.7832358e-1_dp,&
       0.2189568e-1_dp,-0.1062446e-1_dp,0.587872e-2_dp,&
       -0.251540e-2_dp,0.53208e-3_dp/)
  call assert(all(x > 0.0), 'dbessk0_v arg')
  mask = (x <= 2.0)
  where (mask)
     y=x*x/4.0_dp
     dbessk0_v=(-log(x/2.0_dp)*bessi0(x))+poly(y,p,mask)
  elsewhere
     y=(2.0_dp/x)
     dbessk0_v=(exp(-x)/sqrt(x))*poly(y,q,.not. mask)
  end where
END FUNCTION dbessk0_v


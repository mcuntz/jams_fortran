function erfcc_s(x)

  use mo_kind, only: sp
  use mo_nrutil, only: poly
  
  implicit none

  real(sp), intent(in) :: x
  real(sp)             :: erfcc_s

  real(sp) :: t,z
  real(sp), dimension(10) :: coef = (/ -1.26551223_sp, 1.00002368_sp, &
       0.37409196_sp, 0.09678418_sp, -0.18628806_sp, 0.27886807_sp, &
       -1.13520398_sp, 1.48851587_sp, -0.82215223_sp, 0.17087277_sp /)

  z = abs(x)
  t = 1.0_sp / (1.0_sp+0.5_sp*z)
  
  erfcc_s = t * exp(-z*z + poly(t,coef))
  if (x < 0.0_sp) erfcc_s = 2.0_sp - erfcc_s
  
end function erfcc_s


function erfcc_v(x)

  use mo_kind,   only: sp
  use mo_nrutil, only: poly
  
  implicit none
  
  real(sp), dimension(:), intent(in) :: x
  real(sp), dimension(size(x))       :: erfcc_v

  real(sp), dimension(size(x)) :: t,z
  real(sp), dimension(10) :: coef = (/ -1.26551223_sp, 1.00002368_sp, &
       0.37409196_sp, 0.09678418_sp, -0.18628806_sp, 0.27886807_sp, &
       -1.13520398_sp, 1.48851587_sp, -0.82215223_sp, 0.17087277_sp /)

  z = abs(x)
  t = 1.0_sp / (1.0_sp+0.5_sp*z)

  erfcc_v = t * exp(-z*z + poly(t,coef))
  where (x < 0.0_sp) erfcc_v = 2.0_sp - erfcc_v
  
end function erfcc_v



function derfcc_s(x)

  use mo_kind, only: dp
  use mo_nrutil, only: poly
  
  implicit none

  real(dp), intent(in) :: x
  real(dp)             :: derfcc_s

  real(dp) :: t,z
  real(dp), dimension(10) :: coef = (/ -1.26551223_dp, 1.00002368_dp, &
       0.37409196_dp, 0.09678418_dp, -0.18628806_dp, 0.27886807_dp, &
       -1.13520398_dp, 1.48851587_dp, -0.82215223_dp, 0.17087277_dp /)

  z = abs(x)
  t = 1.0_dp / (1.0_dp+0.5_dp*z)
  
  derfcc_s = t * exp(-z*z + poly(t,coef))
  if (x < 0.0_dp) derfcc_s = 2.0_dp - derfcc_s
  
end function derfcc_s


function derfcc_v(x)

  use mo_kind,   only: dp
  use mo_nrutil, only: poly
  
  implicit none
  
  real(dp), dimension(:), intent(in) :: x
  real(dp), dimension(size(x))       :: derfcc_v

  real(dp), dimension(size(x)) :: t,z
  real(dp), dimension(10) :: coef = (/ -1.26551223_dp, 1.00002368_dp, &
       0.37409196_dp, 0.09678418_dp, -0.18628806_dp, 0.27886807_dp, &
       -1.13520398_dp, 1.48851587_dp, -0.82215223_dp, 0.17087277_dp /)

  z = abs(x)
  t = 1.0_dp / (1.0_dp+0.5_dp*z)

  derfcc_v = t * exp(-z*z + poly(t,coef))
  where (x < 0.0_dp) derfcc_v = 2.0_dp - derfcc_v
  
end function derfcc_v

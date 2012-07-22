Program main

  use mo_kind, only: dp
  use ModLine, only: line

  implicit none
  
  integer, parameter :: nmax = 14
  integer  :: n
  real(dp) :: x(nmax), y(nmax), t(2), dt(2), s

  n = 14
  x = (/ 3, 4, 6, 6, 6, 7, 8, 9, 11, 11, 12, 12, 14, 16 /)
  y = (/ 24.82, 23.26, 14.77, 19.06, 14.79, 17.66, 11.83, 14.82, 5.16, 11.12, 8.04, 2.72, 0.74, 1.21 /)

  call line(n,x,y,t,dt,s)

end Program main

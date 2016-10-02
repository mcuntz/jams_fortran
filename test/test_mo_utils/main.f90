program test_utils

  use mo_kind,  only: sp, dp, i4
  use mo_utils, only: eq, ge, le, ne
  use mo_utils, only: swap, locate, cumsum, arange, linspace
  use mo_utils, only: is_finite, is_nan, is_normal, special_value

  implicit none

  real(dp) :: a_dp
  real(dp) :: b_dp
  real(sp) :: a_sp
  real(sp) :: b_sp

  integer(i4), parameter :: nn = 100
  real(dp), dimension(nn) :: dat1, dat2, dat3
  real(sp), dimension(nn) :: sat1, sat2, sat3
  integer(i4), dimension(nn) :: iat1, iat2, iat3
  complex(dp), dimension(nn) :: cat1, cat2
  real(dp), dimension(:), allocatable :: adat1

  real(dp) :: d1
  real(dp), dimension(5) :: d5
  real(sp) :: s1
  real(sp), dimension(5) :: s5
  integer(i4) :: ii1
  integer(i4), dimension(5) :: ii5

  integer(i4) :: i
  logical  :: isgood
  logical  :: compare

  isgood = .true.

  ! -----------------------------------------------------
  ! DOUBLE PRECISON
  ! -----------------------------------------------------

  write(*,*) ''
  write(*,*) 'Test: eq/ equal: dp'
  ! 0.1 == 0.1+eps --> .False.
  a_dp = 0.1_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = eq(a_dp, b_dp)
  write(*,'(E24.17,A4,E24.17,A5,L2)') a_dp,' == ',b_dp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! 1.0 == 1.0+eps --> .True.
  a_dp = 1.0_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = eq(a_dp, b_dp)
  write(*,'(E24.17,A4,E24.17,A5,L2)') a_dp,' == ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  write(*,*) ''
  write(*,*) 'Test: ne/ notequal: dp'
  ! 0.1 /= 0.1+eps --> .True.
  a_dp = 0.1_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = ne(a_dp, b_dp)
  write(*,'(E24.17,A4,E24.17,A5,L2)') a_dp,' /= ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 1.0 /= 1.0+eps --> .False.
  a_dp = 1.0_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = ne(a_dp, b_dp)
  write(*,'(E24.17,A4,E24.17,A5,L2)') a_dp,' /= ',b_dp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  write(*,*) ''
  write(*,*) 'Test: le/ lesserequal: dp'
  ! 0.1 <= 0.1+eps  --> .True.
  a_dp = 0.1_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = le(a_dp, b_dp)
  write(*,'(E24.17,A4,E24.17,A5,L2)') a_dp,' <= ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 1.0 <= 1.0+eps  --> .True.
  a_dp = 1.0_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = le(a_dp, b_dp)
  write(*,'(E24.17,A4,E24.17,A5,L2)') a_dp,' <= ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 0.1 <= 0.2  --> .True.
  a_dp = 0.1_dp
  b_dp = 0.2_dp
  compare = le(a_dp, b_dp)
  write(*,'(E24.17,A4,E24.17,A5,L2)') a_dp,' <= ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  ! tiny <= 2*tiny  --> .True.
  a_dp = tiny(1.0_dp)
  b_dp = 2.0_dp*a_dp
  compare = le(a_dp, b_dp)
  write(*,'(E24.17,A4,E24.17,A5,L2)') a_dp,' <= ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  write(*,*) ''
  write(*,*) 'Test: ge/ greaterequal: dp'
  ! 0.1 >= 0.1+eps  --> .False.
  a_dp = 0.1_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = ge(a_dp, b_dp)
  write(*,'(E24.17,A4,E24.17,A5,L2)') a_dp,' >= ',b_dp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! 1.0 >= 1.0+eps  --> .True.
  a_dp = 1.0_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = ge(a_dp, b_dp)
  write(*,'(E24.17,A4,E24.17,A5,L2)') a_dp,' >= ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 0.1 >= 0.2  --> .False.
  a_dp = 0.1_dp
  b_dp = 0.2_dp
  compare = ge(a_dp, b_dp)
  write(*,'(E24.17,A4,E24.17,A5,L2)') a_dp,' >= ',b_dp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! tiny >= 2*tiny  --> .False.
  a_dp = tiny(1.0_dp)
  b_dp = 2.0_dp*a_dp
  compare = ge(a_dp, b_dp)
  write(*,'(E24.17,A4,E24.17,A5,L2)') a_dp,' >= ',b_dp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! -----------------------------------------------------
  ! SINGLE PRECISON
  ! -----------------------------------------------------

  write(*,*) ''
  write(*,*) 'Test: eq/ equal: sp'
  ! 0.1 == 0.1+eps --> .False.
  a_sp = 0.1_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = eq(a_sp, b_sp)
  write(*,'(E15.8,A4,E15.8,A5,L2)') a_sp,' == ',b_sp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! 1.0 == 1.0+eps --> .True.
  a_sp = 1.0_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = eq(a_sp, b_sp)
  write(*,'(E15.8,A4,E15.8,A5,L2)') a_sp,' == ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  write(*,*) ''
  write(*,*) 'Test: ne/ notequal: sp'
  ! 0.1 /= 0.1+eps --> .True.
  a_sp = 0.1_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = ne(a_sp, b_sp)
  write(*,'(E15.8,A4,E15.8,A5,L2)') a_sp,' /= ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 1.0 /= 1.0+eps --> .False.
  a_sp = 1.0_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = ne(a_sp, b_sp)
  write(*,'(E15.8,A4,E15.8,A5,L2)') a_sp,' /= ',b_sp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  write(*,*) ''
  write(*,*) 'Test: le/ lesserequal: sp'
  ! 0.1 <= 0.1+eps  --> .True.
  a_sp = 0.1_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = le(a_sp, b_sp)
  write(*,'(E15.8,A4,E15.8,A5,L2)') a_sp,' <= ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 1.0 <= 1.0+eps  --> .True.
  a_sp = 1.0_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = le(a_sp, b_sp)
  write(*,'(E15.8,A4,E15.8,A5,L2)') a_sp,' <= ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 0.1 <= 0.2  --> .True.
  a_sp = 0.1_sp
  b_sp = 0.2_sp
  compare = le(a_sp, b_sp)
  write(*,'(E15.8,A4,E15.8,A5,L2)') a_sp,' <= ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  ! tiny <= 2*tiny  --> .True.
  a_sp = tiny(1.0_sp)
  b_sp = 2.0_sp*a_sp
  compare = le(a_sp, b_sp)
  write(*,'(E15.8,A4,E15.8,A5,L2)') a_sp,' <= ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  write(*,*) ''
  write(*,*) 'Test: ge/ greaterequal: sp'
  ! 0.1 >= 0.1+eps  --> .False.
  a_sp = 0.1_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = ge(a_sp, b_sp)
  write(*,'(E15.8,A4,E15.8,A5,L2)') a_sp,' >= ',b_sp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! 1.0 >= 1.0+eps  --> .True.
  a_sp = 1.0_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = ge(a_sp, b_sp)
  write(*,'(E15.8,A4,E15.8,A5,L2)') a_sp,' >= ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 0.1 >= 0.2  --> .False.
  a_sp = 0.1_sp
  b_sp = 0.2_sp
  compare = ge(a_sp, b_sp)
  write(*,'(E15.8,A4,E15.8,A5,L2)') a_sp,' >= ',b_sp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! tiny >= 2*tiny  --> .False.
  a_sp = tiny(1.0_sp)
  b_sp = 2.0_sp*a_sp
  compare = ge(a_sp, b_sp)
  write(*,'(E15.8,A4,E15.8,A5,L2)') a_sp,' >= ',b_sp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! -----------------------------------------------------
  ! Swap

  call random_number(dat1)
  call random_number(dat2)
  dat3 = dat1
  call swap(dat1, dat2)
  isgood = isgood .and. all(eq(dat2,dat3))

  call swap(dat2, 1, nn)
  call swap(dat2, nn, 1)
  isgood = isgood .and. all(eq(dat2,dat3))

  call random_number(sat1)
  call random_number(sat2)
  sat3 = sat1
  call swap(sat1, sat2)
  isgood = isgood .and. all(eq(sat2,sat3))

  call swap(sat2, 1, nn)
  call swap(sat2, nn, 1)
  isgood = isgood .and. all(eq(sat2,sat3))

  call random_number(dat1)
  call random_number(dat2)
  iat1 = int(dat1, i4)
  iat2 = int(dat2, i4)
  iat3 = iat1
  call swap(iat1, iat2)
  isgood = isgood .and. all(iat2 == iat3)

  call swap(iat2, 1, nn)
  call swap(iat2, nn, 1)
  isgood = isgood .and. all(iat2 == iat3)

  ! -----------------------------------------------------
  ! Locate

  ! double precision
  forall(i=1:nn) dat1(i) = real(i,dp)
  ! 0d
  d1 = 1.1_dp
  ii1 = locate(dat1, d1)
  isgood = isgood .and. (ii1 == 1)
  ! 1d
  d5 = (/ 0.1_dp, 5.5_dp, 10.1_dp, 50.5_dp, 200.1_dp /)
  ii5 = locate(dat1, d5)
  isgood = isgood .and. all(ii5 == (/0, 5, 10, 50, 100/))

  ! single precision
  forall(i=1:nn) sat1(i) = real(i,sp)
  ! 0d
  s1 = 1.1_sp
  ii1 = locate(sat1, s1)
  isgood = isgood .and. (ii1 == 1)
  ! 1d
  s5 = (/ 0.1_sp, 5.5_sp, 10.1_sp, 50.5_sp, 200.1_sp /)
  ii5 = locate(sat1, s5)
  isgood = isgood .and. all(ii5 == (/0, 5, 10, 50, 100/))

  ! -----------------------------------------------------
  ! is_finite, is_nan, is_normal

  call random_number(dat1)
  isgood = isgood .and. all(is_finite(dat1))
  isgood = isgood .and. all(.not. is_nan(dat1))
  isgood = isgood .and. all(is_normal(dat1))
  ! NaN
  dat1(1) = 0.0_dp
  dat1(1) = dat1(1)/dat1(1)
  isgood = isgood .and. any(is_nan(dat1))
  isgood = isgood .and. is_nan(dat1(1))
  ! Inf
  dat1(2) = huge(1.0_dp)
  dat1(2) = dat1(2)*dat1(2)
  isgood = isgood .and. any(.not. is_finite(dat1))
  isgood = isgood .and. (.not. is_finite(dat1(2)))
  ! Both
  isgood = isgood .and. (.not. all(is_normal(dat1)))
  isgood = isgood .and. (.not. any(is_normal(dat1(1:2))))

  call random_number(sat1)
  isgood = isgood .and. all(is_finite(sat1))
  isgood = isgood .and. all(.not. is_nan(sat1))
  isgood = isgood .and. all(is_normal(sat1))
  ! NaN
  sat1(1) = 0.0_sp
  sat1(1) = sat1(1)/sat1(1)
  isgood = isgood .and. any(is_nan(sat1))
  isgood = isgood .and. is_nan(sat1(1))
  ! Inf
  sat1(2) = huge(1.0_sp)
  sat1(2) = sat1(2)*sat1(2)
  isgood = isgood .and. any(.not. is_finite(sat1))
  isgood = isgood .and. (.not. is_finite(sat1(2)))
  ! Both
  isgood = isgood .and. (.not. all(is_normal(sat1)))
  isgood = isgood .and. (.not. any(is_normal(sat1(1:2))))

  ! -----------------------------------------------------
  ! special_value

  isgood = isgood .and. is_nan(special_value(1.0_dp, 'IEEE_SIGNALING_NAN'))
  isgood = isgood .and. is_nan(special_value(1.0_dp, 'IEEE_QUIET_NAN'))
  isgood = isgood .and. (.not. is_finite(special_value(1.0_dp, 'IEEE_POSITIVE_INF')))
  isgood = isgood .and. (.not. is_finite(special_value(1.0_dp, 'IEEE_NEGATIVE_INF')))
  isgood = isgood .and. is_finite(special_value(1.0_dp, 'IEEE_POSITIVE_DENORMAL'))
  isgood = isgood .and. is_finite(special_value(1.0_dp, 'IEEE_NEGATIVE_DENORMAL'))
  isgood = isgood .and. is_finite(special_value(1.0_dp, 'IEEE_POSITIVE_NORMAL'))
  isgood = isgood .and. is_finite(special_value(1.0_dp, 'IEEE_NEGATIVE_NORMAL'))
  isgood = isgood .and. is_finite(special_value(1.0_dp, 'IEEE_POSITIVE_ZERO'))
  isgood = isgood .and. is_finite(special_value(1.0_dp, 'IEEE_NEGATIVE_ZERO'))
  isgood = isgood .and. eq(special_value(1.0_dp, 'IEEE_NEGATIVE_INF'), -special_value(1.0_dp, 'IEEE_POSITIVE_INF'))
  isgood = isgood .and. eq(special_value(1.0_dp, 'IEEE_NEGATIVE_DENORMAL'), -special_value(1.0_dp, 'IEEE_POSITIVE_DENORMAL'))
  isgood = isgood .and. eq(special_value(1.0_dp, 'IEEE_NEGATIVE_NORMAL'), -special_value(1.0_dp, 'IEEE_POSITIVE_NORMAL'))
  isgood = isgood .and. eq(special_value(1.0_dp, 'IEEE_NEGATIVE_ZERO'), -special_value(1.0_dp, 'IEEE_POSITIVE_ZERO'))
  isgood = isgood .and. eq(abs(special_value(1.0_dp, 'IEEE_POSITIVE_ZERO')), 0.0_dp)
  isgood = isgood .and. eq(abs(special_value(1.0_dp, 'IEEE_NEGATIVE_ZERO')), 0.0_dp)

  isgood = isgood .and. is_nan(special_value(1.0_sp, 'IEEE_SIGNALING_NAN'))
  isgood = isgood .and. is_nan(special_value(1.0_sp, 'IEEE_QUIET_NAN'))
  isgood = isgood .and. (.not. is_finite(special_value(1.0_sp, 'IEEE_POSITIVE_INF')))
  isgood = isgood .and. (.not. is_finite(special_value(1.0_sp, 'IEEE_NEGATIVE_INF')))
  isgood = isgood .and. is_finite(special_value(1.0_sp, 'IEEE_POSITIVE_DENORMAL'))
  isgood = isgood .and. is_finite(special_value(1.0_sp, 'IEEE_NEGATIVE_DENORMAL'))
  isgood = isgood .and. is_finite(special_value(1.0_sp, 'IEEE_POSITIVE_NORMAL'))
  isgood = isgood .and. is_finite(special_value(1.0_sp, 'IEEE_NEGATIVE_NORMAL'))
  isgood = isgood .and. is_finite(special_value(1.0_sp, 'IEEE_POSITIVE_ZERO'))
  isgood = isgood .and. is_finite(special_value(1.0_sp, 'IEEE_NEGATIVE_ZERO'))
  isgood = isgood .and. eq(special_value(1.0_sp, 'IEEE_NEGATIVE_INF'), -special_value(1.0_sp, 'IEEE_POSITIVE_INF'))
  isgood = isgood .and. eq(special_value(1.0_sp, 'IEEE_NEGATIVE_DENORMAL'), -special_value(1.0_sp, 'IEEE_POSITIVE_DENORMAL'))
  isgood = isgood .and. eq(special_value(1.0_sp, 'IEEE_NEGATIVE_NORMAL'), -special_value(1.0_sp, 'IEEE_POSITIVE_NORMAL'))
  isgood = isgood .and. eq(special_value(1.0_sp, 'IEEE_NEGATIVE_ZERO'), -special_value(1.0_sp, 'IEEE_POSITIVE_ZERO'))
  isgood = isgood .and. eq(abs(special_value(1.0_sp, 'IEEE_POSITIVE_ZERO')), 0.0_sp)
  isgood = isgood .and. eq(abs(special_value(1.0_sp, 'IEEE_NEGATIVE_ZERO')), 0.0_sp)

  ! -----------------------------------------------------
  ! Cumsum
  ! double precision
  dat1(:) = 1.0_dp
  dat2    = cumsum(dat1)
  isgood  = isgood .and. (eq(dat2(1),1.0_dp))
  isgood  = isgood .and. (eq(dat2(nn),real(nn,dp)))
  ! single precision
  sat1(:) = 2.0_sp
  sat2    = cumsum(sat1)
  isgood  = isgood .and. (eq(sat2(1),2.0_sp))
  isgood  = isgood .and. (eq(sat2(nn),real(2*nn,sp)))
  ! integer
  iat1(:) = 3
  iat2    = cumsum(iat1)
  isgood  = isgood .and. (iat2(1) == 3)
  isgood  = isgood .and. (iat2(nn) == 3*nn)
  ! complex
  cat1(:) = (1.0_dp,1.0_dp)
  cat2    = cumsum(cat1)
  isgood  = isgood .and. (eq(real(cat2(1)),1.0_dp))
  isgood  = isgood .and. (eq(aimag(cat2(nn)),real(nn,dp)))

  ! -----------------------------------------------------
  ! Range
  ! double precision
  dat1(:) = 1.0_dp
  dat2    = cumsum(dat1)
  dat1    = arange(real(nn,dp))
  isgood  = isgood .and. all(eq(dat1,dat2))
  ! single precision
  sat1(:) = 2.0_dp
  sat2    = cumsum(sat1)
  sat1    = arange(real(nn,sp)) * 2.0_sp
  isgood  = isgood .and. all(eq(sat1,sat2))
  ! integer
  iat1(:) = 1
  iat2    = cumsum(iat1) - 2
  iat1    = arange(-1,nn-2)
  isgood  = isgood .and. all(iat1==iat2)
  ! allocatable out
  adat1 = arange(real(nn,dp))

  ! -----------------------------------------------------
  ! Linspace
  ! double precision
  dat1   = arange(real(nn,dp))
  dat2   = linspace(1.0_dp,real(nn,dp),nn)
  isgood = isgood .and. all(eq(dat1,dat2))
  ! single precision
  sat1   = arange(real(nn,sp))/real(nn,sp)
  sat2   = linspace(0.01_sp,1.0_sp,nn)
  isgood = isgood .and. all(eq(sat1,sat2))
  ! integer
  iat1(:) = 3
  iat2    = cumsum(iat1)
  iat1    = linspace(3,3*nn,nn)
  isgood  = isgood .and. all(iat1==iat2)


  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_utils o.k.'
  else
     write(*,*) 'mo_utils failed!'
  endif

end program test_utils

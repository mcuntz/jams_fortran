program test_utils

  use mo_kind,  only: sp, dp
  use mo_utils, only: eq, ge, le, ne

  implicit none

  real(dp) :: a_dp
  real(dp) :: b_dp
  real(sp) :: a_sp
  real(sp) :: b_sp

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
  write(*,'(E23.17,A4,E23.17,A5,L2)') a_dp,' == ',b_dp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! 1.0 == 1.0+eps --> .True.
  a_dp = 1.0_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = eq(a_dp, b_dp)
  write(*,'(E23.17,A4,E23.17,A5,L2)') a_dp,' == ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  write(*,*) ''
  write(*,*) 'Test: ne/ notequal: dp'
  ! 0.1 /= 0.1+eps --> .True.
  a_dp = 0.1_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = ne(a_dp, b_dp)
  write(*,'(E23.17,A4,E23.17,A5,L2)') a_dp,' /= ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 1.0 /= 1.0+eps --> .False.
  a_dp = 1.0_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = ne(a_dp, b_dp)
  write(*,'(E23.17,A4,E23.17,A5,L2)') a_dp,' /= ',b_dp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  write(*,*) ''
  write(*,*) 'Test: le/ lesserequal: dp'
  ! 0.1 <= 0.1+eps  --> .True.
  a_dp = 0.1_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = le(a_dp, b_dp)
  write(*,'(E23.17,A4,E23.17,A5,L2)') a_dp,' <= ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 1.0 <= 1.0+eps  --> .True.
  a_dp = 1.0_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = le(a_dp, b_dp)
  write(*,'(E23.17,A4,E23.17,A5,L2)') a_dp,' <= ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 0.1 <= 0.2  --> .True.
  a_dp = 0.1_dp
  b_dp = 0.2_dp
  compare = le(a_dp, b_dp)
  write(*,'(E23.17,A4,E23.17,A5,L2)') a_dp,' <= ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  ! tiny <= 2*tiny  --> .True.
  a_dp = tiny(1.0_dp)
  b_dp = 2.0_dp*a_dp
  compare = le(a_dp, b_dp)
  write(*,'(E23.17,A4,E23.17,A5,L2)') a_dp,' <= ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  write(*,*) ''
  write(*,*) 'Test: ge/ greaterequal: dp'
  ! 0.1 >= 0.1+eps  --> .False.
  a_dp = 0.1_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = ge(a_dp, b_dp)
  write(*,'(E23.17,A4,E23.17,A5,L2)') a_dp,' >= ',b_dp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! 1.0 >= 1.0+eps  --> .True.
  a_dp = 1.0_dp
  b_dp = a_dp + epsilon(1.0_dp)
  compare = ge(a_dp, b_dp)
  write(*,'(E23.17,A4,E23.17,A5,L2)') a_dp,' >= ',b_dp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 0.1 >= 0.2  --> .False.
  a_dp = 0.1_dp
  b_dp = 0.2_dp
  compare = ge(a_dp, b_dp)
  write(*,'(E23.17,A4,E23.17,A5,L2)') a_dp,' >= ',b_dp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! tiny >= 2*tiny  --> .False.
  a_dp = tiny(1.0_dp)
  b_dp = 2.0_dp*a_dp
  compare = ge(a_dp, b_dp)
  write(*,'(E23.17,A4,E23.17,A5,L2)') a_dp,' >= ',b_dp,' --> ',compare
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
  write(*,'(E14.8,A4,E14.8,A5,L2)') a_sp,' == ',b_sp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! 1.0 == 1.0+eps --> .True.
  a_sp = 1.0_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = eq(a_sp, b_sp)
  write(*,'(E14.8,A4,E14.8,A5,L2)') a_sp,' == ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  write(*,*) ''
  write(*,*) 'Test: ne/ notequal: sp'
  ! 0.1 /= 0.1+eps --> .True.
  a_sp = 0.1_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = ne(a_sp, b_sp)
  write(*,'(E14.8,A4,E14.8,A5,L2)') a_sp,' /= ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 1.0 /= 1.0+eps --> .False.
  a_sp = 1.0_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = ne(a_sp, b_sp)
  write(*,'(E14.8,A4,E14.8,A5,L2)') a_sp,' /= ',b_sp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  write(*,*) ''
  write(*,*) 'Test: le/ lesserequal: sp'
  ! 0.1 <= 0.1+eps  --> .True.
  a_sp = 0.1_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = le(a_sp, b_sp)
  write(*,'(E14.8,A4,E14.8,A5,L2)') a_sp,' <= ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 1.0 <= 1.0+eps  --> .True.
  a_sp = 1.0_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = le(a_sp, b_sp)
  write(*,'(E14.8,A4,E14.8,A5,L2)') a_sp,' <= ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 0.1 <= 0.2  --> .True.
  a_sp = 0.1_sp
  b_sp = 0.2_sp
  compare = le(a_sp, b_sp)
  write(*,'(E14.8,A4,E14.8,A5,L2)') a_sp,' <= ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  ! tiny <= 2*tiny  --> .True.
  a_sp = tiny(1.0_sp)
  b_sp = 2.0_sp*a_sp
  compare = le(a_sp, b_sp)
  write(*,'(E14.8,A4,E14.8,A5,L2)') a_sp,' <= ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  write(*,*) ''
  write(*,*) 'Test: ge/ greaterequal: sp'
  ! 0.1 >= 0.1+eps  --> .False.
  a_sp = 0.1_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = ge(a_sp, b_sp)
  write(*,'(E14.8,A4,E14.8,A5,L2)') a_sp,' >= ',b_sp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! 1.0 >= 1.0+eps  --> .True.
  a_sp = 1.0_sp
  b_sp = a_sp + epsilon(1.0_sp)
  compare = ge(a_sp, b_sp)
  write(*,'(E14.8,A4,E14.8,A5,L2)') a_sp,' >= ',b_sp,' --> ',compare
  isgood = isgood .and. (compare)

  ! 0.1 >= 0.2  --> .False.
  a_sp = 0.1_sp
  b_sp = 0.2_sp
  compare = ge(a_sp, b_sp)
  write(*,'(E14.8,A4,E14.8,A5,L2)') a_sp,' >= ',b_sp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  ! tiny >= 2*tiny  --> .False.
  a_sp = tiny(1.0_sp)
  b_sp = 2.0_sp*a_sp
  compare = ge(a_sp, b_sp)
  write(*,'(E14.8,A4,E14.8,A5,L2)') a_sp,' >= ',b_sp,' --> ',compare
  isgood = isgood .and. (.not. compare)

  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_utils o.k.'
  else
     write(*,*) 'mo_utils failed!'
  endif

end program test_utils



!*******************************************************
!
!   TO TEST xor4096_apps
!
!*******************************************************
program combinatorics

  use mo_kind,          only: i4
  use mo_combinatorics, only: binomcoeffi

  implicit none

  integer(i4) :: n_i4,k_i4
  integer(i4) :: c_i4

  logical :: isgood


  isgood = .true.
  write(*,*) ''

  !--------------------------------------------------
  ! Test Binomial Coefficient
  !--------------------------------------------------

  n_i4 = 5_i4
  k_i4 = 3_i4
  c_i4 = binomcoeffi(n_i4,k_i4)
  write(*,*) '5 choose 3 = ', c_i4

  if ( c_i4 .ne. 10_i4 ) isgood = .false.

  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_combinatorics o.k.'
  else
     write(*,*) 'mo_combinatorics failed'
  end if

end program combinatorics

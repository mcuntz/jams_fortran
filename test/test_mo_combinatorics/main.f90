!*******************************************************
!
!   TO TEST xor4096_apps
!
!*******************************************************
program combinatorics

  use mo_kind,          only: i4, sp
  use mo_combinatorics, only: binomcoeffi, next_kofn, all_kofn, random_kofn
  use mo_xor4096,       only: xor4096, n_save_state

  implicit none

  integer(i4) :: c_i4
  integer(i4) :: i
  integer(i4) :: next_i4(3)
  integer(i4) :: iseed
  real(sp)    :: rn
  integer(i4), dimension(n_save_state) :: save_state_i4
  integer(i4) :: random_i4(3)
  integer(i4), allocatable :: all_i4(:,:)

  logical :: isgood


  isgood = .true.
  write(*,*) ''

  !--------------------------------------------------
  ! Test: Binomial Coefficient
  !--------------------------------------------------

  c_i4 = binomcoeffi(5, 3)

  write(*,*) 'BINOMCOEFFI (5,3)                : ', c_i4
  if ( c_i4 .ne. 10_i4 ) isgood = .false.

  !--------------------------------------------------
  ! Test: Next k of n subset
  !--------------------------------------------------

  next_i4 = next_kofn(5,3,(/1,3,5/))

  write(*,*) 'NEXT_KOFN subset after (/1,3,5/) : ', next_i4
  if ( any(next_i4 .ne. (/ 1_i4, 4_i4, 5_i4 /)) ) isgood = .false.

  !--------------------------------------------------
  ! Test: All k of n subset
  !--------------------------------------------------

  all_i4 = all_kofn(5,3)

  write(*,*) 'ALL_KOFN (5,3)                  : '
  do i=1,size(all_i4,1)
     write(*,*) '                              subset #',i,'   ',all_i4(i,:)
  end do
  if ( any(all_i4(5,:) .ne. (/ 1_i4, 3_i4, 5_i4 /)) ) isgood = .false.

  !--------------------------------------------------
  ! Test: Random k of n subset
  !--------------------------------------------------

  iseed=1045_i4
  call xor4096(iSeed, rn, save_state=save_state_i4)

  write(*,*) 'RANDOM_KOFN (5,3)               : '
  do i=1,6
     random_i4 = random_kofn(5,3,save_state=save_state_i4)
     write(*,*) '                              subset #',i,'   ',random_i4
  end do
  if ( any(random_i4(:) .ne. (/ 3_i4, 4_i4, 5_i4 /)) ) isgood = .false.


  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_combinatorics o.k.'
  else
     write(*,*) 'mo_combinatorics failed'
  end if

end program combinatorics

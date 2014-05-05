!*******************************************************
!
!   TO TEST mo_combinatorics
!
!*******************************************************
program combinatorics

  use mo_kind,          only: i4, sp
  use mo_combinatorics, only: binomcoeffi!, factorial
  use mo_functions,     only: factorial
  use mo_combinatorics, only: next_kofn,   all_kofn,   random_kofn
  use mo_combinatorics, only: next_index_permut, all_index_permut, random_index_permut, random_permut
  use mo_xor4096,       only: xor4096, n_save_state

  implicit none

  integer(i4)                          :: c_i4
  integer(i4)                          :: i
  integer(i4), dimension(3)            :: next_i4
  integer(i4), dimension(4)            :: next_4_i4
  integer(i4)                          :: iseed
  real(sp)                             :: rn
  integer(i4), dimension(n_save_state) :: save_state_i4
  integer(i4), dimension(3)            :: random_i4
  integer(i4), allocatable             :: all_i4(:,:)
  integer(i4), dimension(4)            :: permut_i4
  logical                              :: isgood
  ! integer(i4)                          :: j, nn
  ! integer(i4), dimension(4,4)          :: count_permut_i4


  isgood = .true.
  write(*,*) ''

  !--------------------------------------------------
  ! Test: Binomial Coefficient
  !--------------------------------------------------

  c_i4 = binomcoeffi(5, 3)
  write(*,*) 'BINOMCOEFFI (5,3)                : ', c_i4
  if ( c_i4 .ne. 10_i4 ) isgood = .false.

  !--------------------------------------------------
  ! Test: Factorial
  !--------------------------------------------------

  c_i4 = factorial(5)
  write(*,*) 'FACTORIAL 5!                     : ', c_i4
  if ( c_i4 .ne. 120_i4 ) isgood = .false.
  write(*,*) ' '

  !--------------------------------------------------
  ! Test: Random k of n subset
  !--------------------------------------------------

  iseed=1045_i4
  call xor4096(iSeed, rn, save_state=save_state_i4)

  random_i4 = random_kofn(5,3,save_state=save_state_i4)
  write(*,*) 'RANDOM_KOFN (5,3)                : ', random_i4
  if ( any(random_i4(:) .ne. (/ 1_i4, 3_i4, 4_i4 /)) ) isgood = .false.

  !--------------------------------------------------
  ! Test: Next k of n subset
  !--------------------------------------------------

  next_i4 = next_kofn(5,3,(/1,4,5/))
  write(*,*) 'NEXT_KOFN subset after (/1,4,5/) : ', next_i4
  if ( any(next_i4 .ne. (/ 2_i4, 3_i4, 4_i4 /)) ) isgood = .false.

  !--------------------------------------------------
  ! Test: All k of n subset
  !--------------------------------------------------

  allocate(all_i4(binomcoeffi(5,3),3))

  all_i4 = all_kofn(5,3)
  write(*,*) 'ALL_KOFN (5,3)                   : '
  do i=1,size(all_i4,1)
     write(*,*) '         subset #',i,'   ',all_i4(i,:)
  end do
  if ( any(all_i4(5,:) .ne. (/ 1_i4, 3_i4, 5_i4 /)) ) isgood = .false.
  write(*,*) ' '

  deallocate(all_i4)

  !--------------------------------------------------
  ! Test: Random index permutation
  !--------------------------------------------------

  iseed=1206_i4
  call xor4096(iSeed, rn, save_state=save_state_i4)

  ! count_permut_i4 = 0
  ! nn = 100000
  ! do i=1, nn
  !    permut_i4 = random_index_permut(4_i4)
  !    forall(j=1:4) count_permut_i4(permut_i4(j),j) = count_permut_i4(permut_i4(j),j) + 1
  !    !write(*,*) 'RANDOM_INDEX_PERMUT of 4                  : ',permut_i4
  ! end do
  ! print*, nn/4
  ! print*, count_permut_i4
  permut_i4 = random_index_permut(4_i4)
  write(*,*) 'RANDOM_INDEX_PERMUT of 4                  : ',permut_i4
  if ( any(permut_i4(:) .ne. (/ 4_i4, 1_i4, 2_i4, 3_i4 /)) ) isgood = .false.

  !--------------------------------------------------
  ! Test: Next index permutation
  !--------------------------------------------------

  next_4_i4 = next_index_permut(4,(/4,1,2,3/))
  write(*,*) 'NEXT_INDEX_PERMUT permutation after (/4,1,2,3/) : ', next_4_i4
  if ( any(next_4_i4(:) .ne. (/ 2_i4, 1_i4, 4_i4, 3_i4 /)) ) isgood = .false.

  !--------------------------------------------------
  ! Test: All index permutations
  !--------------------------------------------------
  allocate(all_i4(factorial(3),3))

  all_i4 = all_index_permut(3)
  write(*,*) 'ALL_INDEX_PERMUT of 3                     : '
  do i=1,size(all_i4,1)
     write(*,*) '    permutation #',i,'   ',all_i4(i,:)
  end do
  if ( any(all_i4(6,:) .ne. (/ 3_i4, 2_i4, 1_i4 /)) ) isgood = .false.

  deallocate(all_i4)

  !--------------------------------------------------
  ! Test: Random permutation
  !--------------------------------------------------

  iseed=1206_i4
  call xor4096(iSeed, rn, save_state=save_state_i4)

  forall(i=1:4) permut_i4(i) = i+6
  call random_permut(permut_i4)
  write(*,*) 'RANDOM_PERMUT of 4                  : ',permut_i4
  if ( any(permut_i4(:) .ne. (/ 10_i4, 7_i4, 8_i4, 9_i4 /)) ) isgood = .false.


  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_combinatorics o.k.'
  else
     write(*,*) 'mo_combinatorics failed'
  end if

end program combinatorics

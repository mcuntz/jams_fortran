!*******************************************************
!
!   TO TEST mo_combinatorics
!
!*******************************************************
program combinatorics

  use mo_kind,          only: i4, i8, sp, dp
  use mo_combinatorics, only: binomcoeffi!, factorial
  use mo_combinatorics, only: nextpart, prevpart, isgoodpart
  use mo_functions,     only: factorial
  use mo_combinatorics, only: next_kofn,   all_kofn,   random_kofn
  use mo_combinatorics, only: next_index_permut, all_index_permut, random_index_permut, random_permut
  use mo_xor4096,       only: xor4096, n_save_state

  implicit none

  integer(i4)                          :: c_i4
  integer(i8)                          :: c_i8
  real(sp)                             :: c_sp
  real(dp)                             :: c_dp
  integer(i4)                          :: i
  integer(i4), dimension(3)            :: next_i4
  integer(i4), dimension(4)            :: next_4_i4
  integer(i4)                          :: iseed
  real(sp)                             :: rn
  integer(i4), dimension(n_save_state) :: save_state_i4
  integer(i4), dimension(3)            :: random_i4
  integer(i4), allocatable             :: all_i4(:,:)
  integer(i4), dimension(4)            :: permut_i4

  integer(i4), dimension(4)            :: part_i4
  integer(i4), dimension(4)            :: part_i8

!  integer(i4), dimension(5)            :: part5_i4

  logical                              :: isgood

  ! integer(i4)                          :: j, nn
  ! integer(i4), dimension(4,4)          :: count_permut_i4


  isgood = .true.
  write(*,*) ''

  !--------------------------------------------------
  ! Test: Binomial Coefficient
  !--------------------------------------------------

  c_i4 = binomcoeffi(5_i4, 3_i4)
  write(*,*) 'BINOMCOEFFI (5,3)                : ', c_i4
  if ( c_i4 .ne. 10_i4 ) isgood = .false.

  c_i8 = binomcoeffi(5_i8, 3_i8)
  write(*,*) 'BINOMCOEFFI (5,3)                : ', c_i8
  if ( c_i8 .ne. 10_i8 ) isgood = .false.

  c_sp = binomcoeffi(5.5_sp, 3_i4)
  write(*,*) 'BINOMCOEFFI (11/2,3)             : ', c_sp
  if ( int(c_sp*10000) .ne. 144375 ) isgood = .false.

  c_dp = binomcoeffi(5.5_dp, 3_i8)
  write(*,*) 'BINOMCOEFFI (11/2,3)             : ', c_dp
  if ( int(c_dp*10000) .ne. 144375 ) isgood = .false.

  
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

  !--------------------------------------------------
  ! Test: Partitions
  !--------------------------------------------------

  part_i4       = 0_i4
  part_i8       = 0_i8
  
  part_i4(4)    = 1_i4
  part_i8(4)    = 1_i8

  write (*,*) " "
  write (*,*) "Calculation of the Partition:"
  write (*,*) "1. Partition_i4 of 4: ", part_i4
  if(.not. isgoodpart(part_i4)) isgood  = .false.
  part_i4       =   nextpart(part_i4)
  write (*,*) "2. Partition_i4 of 4: ", part_i4
  if(.not. (isgoodpart(part_i4) .and. all(part_i4 == (/1_i4,0_i4,1_i4,0_i4/) ))) isgood  = .false.
  part_i4       =   nextpart(part_i4)
  write (*,*) "3. Partition_i4 of 4: ", part_i4
  if(.not. (isgoodpart(part_i4) .and. all(part_i4 == (/0_i4,2_i4,0_i4,0_i4/) ))) isgood  = .false.
  part_i4       =   nextpart(part_i4)
  write (*,*) "4. Partition_i4 of 4: ", part_i4
  if(.not. (isgoodpart(part_i4) .and. all(part_i4 == (/2_i4,1_i4,0_i4,0_i4/) ))) isgood  = .false.
  part_i4       =   nextpart(part_i4)
  write (*,*) "5. Partition_i4 of 4: ", part_i4
  if(.not. (isgoodpart(part_i4) .and. all(part_i4 == (/4_i4,0_i4,0_i4,0_i4/) ))) isgood  = .false.

  write (*,*) "1. Partition_i8 of 4: ", part_i8
  if(.not. isgoodpart(part_i8)) isgood  = .false.
  part_i8       =   nextpart(part_i8)
  write (*,*) "2. Partition_i8 of 4: ", part_i8
  if(.not. (isgoodpart(part_i8) .and. all(part_i8 == (/1_i8,0_i8,1_i8,0_i8/) ))) isgood  = .false.
  part_i8       =   nextpart(part_i8)
  write (*,*) "3. Partition_i8 of 4: ", part_i8
  if(.not. (isgoodpart(part_i8) .and. all(part_i8 == (/0_i8,2_i8,0_i8,0_i8/) ))) isgood  = .false.
  part_i8       =   nextpart(part_i8)
  write (*,*) "4. Partition_i8 of 4: ", part_i8
  if(.not. (isgoodpart(part_i8) .and. all(part_i8 == (/2_i8,1_i8,0_i8,0_i8/) ))) isgood  = .false.
  part_i8       =   nextpart(part_i8)
  write (*,*) "5. Partition_i8 of 4: ", part_i8
  if(.not. (isgoodpart(part_i8) .and. all(part_i8 == (/4_i8,0_i8,0_i8,0_i8/) ))) isgood  = .false.

  write (*,*) " "
  write (*,*) "Reverse calculation of the Partition:"
  write (*,*) "5. Partition_i4 of 4: ", part_i4
  if(.not. (isgoodpart(part_i4) .and. all(part_i4 == (/4_i4,0_i4,0_i4,0_i4/) ))) isgood  = .false.
  part_i4       =   prevpart(part_i4)
  write (*,*) "4. Partition_i4 of 4: ", part_i4
  if(.not. (isgoodpart(part_i4) .and. all(part_i4 == (/2_i4,1_i4,0_i4,0_i4/) ))) isgood  = .false.
  part_i4       =   prevpart(part_i4)
  write (*,*) "3. Partition_i4 of 4: ", part_i4
  if(.not. (isgoodpart(part_i4) .and. all(part_i4 == (/0_i4,2_i4,0_i4,0_i4/) ))) isgood  = .false.
  part_i4       =   prevpart(part_i4)
  write (*,*) "2. Partition_i4 of 4: ", part_i4
  if(.not. (isgoodpart(part_i4) .and. all(part_i4 == (/1_i4,0_i4,1_i4,0_i4/) ))) isgood  = .false.
  part_i4       =   prevpart(part_i4)
  write (*,*) "1. Partition_i4 of 4: ", part_i4
  if(.not. (isgoodpart(part_i4) .and. all(part_i4 == (/0_i4,0_i4,0_i4,1_i4/) ))) isgood  = .false.

  write (*,*) "5. Partition_i8 of 4: ", part_i8
  if(.not. (isgoodpart(part_i8) .and. all(part_i8 == (/4_i8,0_i8,0_i8,0_i8/) ))) isgood  = .false.
  part_i8       =   prevpart(part_i8)
  write (*,*) "4. Partition_i8 of 4: ", part_i8
  if(.not. (isgoodpart(part_i8) .and. all(part_i8 == (/2_i8,1_i8,0_i8,0_i8/) ))) isgood  = .false.
  part_i8       =   prevpart(part_i8)
  write (*,*) "3. Partition_i8 of 4: ", part_i8
  if(.not. (isgoodpart(part_i8) .and. all(part_i8 == (/0_i8,2_i8,0_i8,0_i8/) ))) isgood  = .false.
  part_i8       =   prevpart(part_i8)
  write (*,*) "2. Partition_i8 of 4: ", part_i8
  if(.not. (isgoodpart(part_i8) .and. all(part_i8 == (/1_i8,0_i8,1_i8,0_i8/) ))) isgood  = .false.
  part_i8       =   prevpart(part_i8)
  write (*,*) "1. Partition_i8 of 4: ", part_i8
  if(.not. (isgoodpart(part_i8) .and. all(part_i8 == (/0_i8,0_i8,0_i8,1_i8/) ))) isgood  = .false.

!  part5_i4      = 0_i4
!  part5_i4(5)   = 1_i4
!  
!  write(*,*) ""
!  write(*,*) "Partitions of 5:"
!  do i=1_i4, 6
!  write (*,*) part5_i4
!  part5_i4      =   nextpart(part5_i4)
!  end do
!  write (*,*) part5_i4

!  write(*,*) ""
!  write(*,*) "reverse Partitions of 5:"
!  do i=1_i4, 6
!  write (*,*) part5_i4
!  part5_i4      =   prevpart(part5_i4)
!  end do
!  write (*,*) part5_i4

  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_combinatorics o.k.'
  else
     write(*,*) 'mo_combinatorics failed'
  end if

end program combinatorics

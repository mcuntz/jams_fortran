PROGRAM main
  
  USE mo_kind,   ONLY: i4, i8
  USE mo_ranper, ONLY: ranper

  IMPLICIT NONE
  
  
  LOGICAL :: isgood
  integer(i4), dimension(10) :: A1
  integer(i4), dimension(10) :: A2
  integer(i8), dimension(10) :: B

  Write(*,*) ''
  Write(*,*) 'Test mo_ranper.f90'

  isgood = .true.
  
  call ranper( A1, .true., .true.)
  call ranper( A2, .true., .false.)
  
  isgood = any(A1 /= A2) .and. (sum( A1 ) == (size(A1,1) * (size(A1,1)+1)) / 2)

  B = (/ 2_i8, 4_i8, 5_i8, 3_i8, 9_i8, 7_i8, 1_i8, 10_i8, 8_i8, 6_i8 /)
  call ranper( B, .false.)
  isgood = isgood .and. ( sum( B ) == size(B,1) * (size(B,1)+1) / 2 )

  Write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_ranper o.k.'
  else
     write(*,*) 'mo_ranper failed!'
  endif

END PROGRAM main

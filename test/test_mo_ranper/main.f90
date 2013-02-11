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
  
  call ranper( A1, .true.)
  call ranper( A2, .true., .false.)
  
  isgood = any(A1 /= A2) .and. (sum( A1 ) == (size(A1,1) * (size(A1,1)+1)) / 2)

  B = (/ 2_i4, 4_i4, 5_i4, 3_i4, 9_i4, 7_i4, 1_i4, 10_i4, 8_i4, 6_i4 /)
  call ranper( B, .false.)
  isgood = isgood .and. ( sum( B ) == size(B,1) * (size(B,1)+1) / 2 )

  Write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_ranper o.k.'
  else
     write(*,*) 'mo_ranper failed!'
  endif

END PROGRAM main

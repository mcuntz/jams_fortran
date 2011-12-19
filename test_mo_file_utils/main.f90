PROGRAM main
  
  USE mo_file_utils, ONLY: find_next_unit

  IMPLICIT NONE
  
  INTEGER :: inew
  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_file_utils.f90'

  isgood = .true.
  open(unit=20, file="main.f90", action="read", status="old", form="formatted")
  open(unit=21, file="mo_file_utils.f90", action="read", status="old", form="formatted")
  inew = find_next_unit(20,30)
  if (inew /= 22) isgood =.false.

  if (isgood) then
     write(*,*) 'mo_file_utils o.k.'
  else
     write(*,*) 'mo_file_utils failed!'
  endif

END PROGRAM main

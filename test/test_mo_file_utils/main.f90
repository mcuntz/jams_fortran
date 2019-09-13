PROGRAM main

  use mo_kind,       only: i4
  use mo_file_utils, only: find_next_unit, lines_in_file, lif, next_unit

  implicit none

  integer(i4) :: inew, ilines, islines
  logical     :: isgood

  write(*,*) ''
  write(*,*) 'Test mo_file_utils.f90'

  ! find_next_unit
  isgood = .true.
  open(unit=20, file="main.f90", action="read", status="old", form="formatted")
  open(unit=21, file="mo_file_utils.f90", action="read", status="old", form="formatted")
  inew = find_next_unit(20,30)
  if (inew /= 22) isgood =.false.
  inew = find_next_unit(20)
  if (inew /= 22) isgood =.false.
  inew = find_next_unit()
  if (inew /= 22) isgood =.false.

  close(20)
  close(21)


  ! next_unit
  isgood = .true.
  open(unit=20, file="main.f90", action="read", status="old", form="formatted")
  open(unit=21, file="mo_file_utils.f90", action="read", status="old", form="formatted")
  inew = next_unit(20,30)
  if (inew /= 22) isgood =.false.
  inew = next_unit(20)
  if (inew /= 22) isgood =.false.
  inew = next_unit()
  if (inew /= 22) isgood =.false.

  close(20)
  close(21)


  ! lines_in_file
  islines = 70
  ilines = lines_in_file("main.f90")
  if (ilines /= islines) isgood =.false.
  ilines = lines_in_file("main.f90", comment='!')
  if (ilines /= islines-4) isgood =.false.
  ilines = lines_in_file("main.f90", comment='!', noblank=.true.)
  if (ilines /= islines-4-16) isgood =.false.


  ! lif
  islines = 70
  ilines = lif("main.f90")
  if (ilines /= islines) isgood =.false.
  ilines = lif("main.f90", comment='!')
  if (ilines /= islines-4) isgood =.false.
  ilines = lif("main.f90", comment='!', noblank=.true.)
  if (ilines /= islines-4-16) isgood =.false.


  if (isgood) then
     write(*,*) 'mo_file_utils o.k.'
  else
     write(*,*) 'mo_file_utils failed!'
  endif

end program main

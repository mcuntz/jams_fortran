PROGRAM main
  
  USE mo_kind, ONLY: i4
  use mo_ansi_colors, only: color, c_red, c_green
  USE mo_nml,  ONLY: open_nml, close_nml, position_nml

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: nnml = 100
  CHARACTER(len=*), PARAMETER :: nfile = 'namelist.txt'
  INTEGER(i4) :: stat

  LOGICAL :: isgood
  
  Write(*,*) ''
  Write(*,*) 'Test mo_nml.f90'

  isgood = .true.
  call open_nml(nfile, nnml)
  call position_nml('name1', unit=nnml, first=.true., status=stat)
  if (stat /= 0) isgood =.false.
  call position_nml('name2', unit=nnml, status=stat)
  !call position_nml('kkkk', unit=nnml)
  if (stat /= 0) isgood =.false.
  call position_nml('name3', unit=nnml, first=.true., status=stat)
  call close_nml()
  
  Write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_nml ', color('o.k.', c_green)
  else
     write(*,*) 'mo_nml ', color('failed!', c_red)
  endif

END PROGRAM main

PROGRAM main
  
  USE mo_finish, ONLY: finish
  use mo_ansi_colors, only: color, c_green

  IMPLICIT NONE
  
  Write(*,*) ''
  Write(*,*) 'Test mo_finish.f90'

  call finish(' mo_finish',color(' o.k.', c_green))

END PROGRAM main

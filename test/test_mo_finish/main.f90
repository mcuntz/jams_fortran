PROGRAM main
  
  USE mo_finish, ONLY: finish

  IMPLICIT NONE
  
  Write(*,*) ''
  Write(*,*) 'Test mo_finish.f90'

  call finish(' mo_finish',' o.k.')

END PROGRAM main

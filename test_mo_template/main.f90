PROGRAM main
  
  USE mo_kind,     ONLY: dp, sp
  USE mo_template, ONLY: mean, PI_sp, PI_dp
  
  IMPLICIT NONE
  
  REAL(dp), DIMENSION(10) :: dat
  REAL(sp), DIMENSION(10) :: sat
  
  Write(*,*) ''
  Write(*,*) 'Test mo_template.f90'

  ! Double precision
  Write(*,*) ''
  Write(*,*) '  Test double precision'

  ! Test the parameters
  Write(*,'(A,f17.15)') '     PI_dp is (target=3.141592653589793): ', PI_dp

  ! Test the function mean
  dat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  Write(*,'(A,10f5.1)') '     dat is: ', dat
  Write(*,'(A,f5.2)') '     Mean of dat (target=5.5):                ', mean(dat)
  ! Test the mask
  Write(*,'(A,f5.2)') '     Mean of dat without last (target=5.0):   ', mean(dat(1:9))
  Write(*,'(A,f5.2)') '     Mean of dat with 10 masked (target=5.0): ', mean(dat, mask=(dat /= 10.))
  

  ! Single precision
  Write(*,*) ''
  Write(*,*) '  Test single precision'

  ! Test the parameters
  Write(*,'(A,f8.6)') '     PI_sp is (target=3.141593): ', PI_sp

  ! Test the function mean
  sat = (/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)
  Write(*,'(A,10f5.1)') '     sat is: ', sat
  Write(*,'(A,f5.2)') '     Mean of sat (target=5.5):                ', mean(sat)
  ! Test the mask
  Write(*,'(A,f5.2)') '     Mean of sat without last (target=5.0):   ', mean(sat(1:9))
  Write(*,'(A,f5.2)') '     Mean of sat with 10 masked (target=5.0): ', mean(sat, mask=(sat /= 10.))

END PROGRAM main

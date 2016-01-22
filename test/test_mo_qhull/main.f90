PROGRAM main

  USE mo_kind,  ONLY: dp, i4
  USE mo_qhull, ONLY: qhull ! sums between two elements in first dimension (multiply with constant)

  IMPLICIT NONE

  INTEGER(i4), PARAMETER :: dim = 2
  INTEGER(i4), PARAMETER :: npoints = 3

  REAL(dp), DIMENSION(dim,npoints) :: points

  INTEGER(i4) :: out, i
  LOGICAL :: isgood
  CHARACTER(len=250) :: str1, str2, ofile

  isgood = .true.
  
  points(:,1) = (/ 0.0_dp, 1.0_dp /)
  points(:,2) = (/ 0.38782323776762651_dp, 0.49585832137731267_dp /)
  points(:,3) = (/ 1.0_dp, 0.0_dp /)

  ofile = 'qhull_make_check_test_file'

  out = qhull(points, outfile=ofile)

  open(unit=101, file='../FORTRAN_chs_lib/test/test_mo_qhull/qhull.res', action='read')
  open(unit=102, file=ofile, action='read')
  do i=1, 11
     read(101,*) str1
     read(102,*) str2
     if (trim(str1) /= trim(str2)) isgood = .false.
  end do

  if (isgood) then
     write(*,*) 'mo_qhull o.k.'
  else
     write(*,*) 'mo_qhull failed!'
  endif

END PROGRAM main

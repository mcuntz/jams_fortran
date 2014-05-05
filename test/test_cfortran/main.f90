PROGRAM main

  USE mo_kind,     ONLY: dp, i4
  USE mo_cfortran, ONLY: fctest ! sums between two elements in first dimension (multiply with constant)

  IMPLICIT NONE

  INTEGER(i4), PARAMETER :: d1 = 5
  INTEGER(i4), PARAMETER :: d2 = 4
  INTEGER(i4), PARAMETER :: n1 = 2
  INTEGER(i4), PARAMETER :: n2 = 3

  REAL(dp), DIMENSION(d1,d2) :: A

  REAL(dp), parameter :: c=2.0

  INTEGER(i4) :: i, j
  CHARACTER(LEN=30) :: form1

  forall(i=1:d1, j=1:d2) A(i,j) = real(i**j,dp)

  write(form1,"(A,I3,A)") "(", d1, "f5.0)"
  do j=1, d2
     write(*,form1) A(:,j)
  end do

  ! Beware C-indexes start with 0
  ! Also C is colunm-major so that one wants to transpose A, perhaps
  write(*,'(A,2f6.1)') 'sum(A(2+1:3+1,1)*2 = (3+4)*2 = 14: ', fctest(A, n1, n2, c), sum(A(n1+1:n2+1,1)*c)
  write(*,'(A,2f6.1)') 'sum(A(1,2+1:3+1)*2 = (1+1)*2 =  4: ', fctest(transpose(A), n1, n2, c), sum(A(1,n1+1:n2+1)*c)

END PROGRAM main

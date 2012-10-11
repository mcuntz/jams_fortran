PROGRAM append_test

use mo_kind,   only: i4
#ifndef ABSOFT
use mo_append, only: append, paste

implicit none

integer(i4), dimension(:),   allocatable  :: vector1_i4, vector2_i4
integer(i4), dimension(2)                 :: vector3_i4

integer(i4), dimension(:,:), allocatable  :: matrix1_i4, matrix2_i4, matrix3_i4, matrix4_i4
#endif

logical                                   :: isgood

isgood = .true.

#ifndef ABSOFT
! ---------------------------------------------------------------------------
!                           TEST APPEND
! ---------------------------------------------------------------------------

allocate(vector2_i4(4))
allocate(matrix1_i4(2,2), matrix2_i4(3,2))

vector2_i4 = 2_i4
vector3_i4 = 3_i4

matrix1_i4 = 5_i4
matrix2_i4 = 3_i4

print*, '(A) non-allocated vector will be appended by  vector2'
print*, '    vector2 = (/ ',vector2_i4,' /) ^ T'
call append(vector1_i4, vector2_i4)
print*, '    result  = (/ ',vector1_i4,' /) ^ T'
print*, ' '
if ( any(vector1_i4 .ne. (/ 2_i4,2_i4,2_i4,2_i4 /)) ) isgood=.false.

print*, '(B) vector1  will be appended by  vector3'
print*, '    vector1 = (/ ',vector1_i4,' /) ^ T'
print*, '    vector3 = (/ ',vector3_i4,' /) ^ T'
call append(vector1_i4, vector3_i4)
print*, '    result  = (/ ',vector1_i4,' /) ^ T'
print*, ' '
if ( any(vector1_i4 .ne. (/ 2_i4,2_i4,2_i4,2_i4, 3_i4,3_i4 /)) ) isgood=.false.

print*, '(C) vector1  will be appended by scalar 5_i4'
print*, '    vector1 = (/ ',vector1_i4,' /) ^ T'
call append(vector1_i4, 5_i4)
print*, '    result  = (/ ',vector1_i4,' /)'
print*, ' '
if ( any(vector1_i4 .ne. (/ 2_i4,2_i4,2_i4,2_i4, 3_i4,3_i4, 5_i4 /)) ) isgood=.false.

print*, '(D) matrix1 will be appended by matrix 2'
print*, '    matrix1 = (/', matrix1_i4(1,:)
print*, '                ', matrix1_i4(2,:),' /)'
print*, '    matrix2 = (/', matrix2_i4(1,:)
print*, '                ', matrix2_i4(2,:)
print*, '                ', matrix2_i4(3,:),' /)'
call append(matrix1_i4, matrix2_i4)
print*, '    result  = (/ ',matrix1_i4(1,:)
print*, '                 ',matrix1_i4(2,:)
print*, '                 ',matrix1_i4(3,:)
print*, '                 ',matrix1_i4(4,:)
print*, '                 ',matrix1_i4(5,:),' /)'
print*, ' '
if ( any(matrix1_i4(1,:) .ne. (/ 5_i4,5_i4 /)) ) isgood=.false.
if ( any(matrix1_i4(2,:) .ne. (/ 5_i4,5_i4 /)) ) isgood=.false.
if ( any(matrix1_i4(3,:) .ne. (/ 3_i4,3_i4 /)) ) isgood=.false.
if ( any(matrix1_i4(4,:) .ne. (/ 3_i4,3_i4 /)) ) isgood=.false.
if ( any(matrix1_i4(5,:) .ne. (/ 3_i4,3_i4 /)) ) isgood=.false.

deallocate(vector2_i4)
deallocate(matrix1_i4, matrix2_i4)

! ---------------------------------------------------------------------------
!                           TEST PASTE
! ---------------------------------------------------------------------------

allocate(matrix2_i4(2,2), matrix3_i4(2,3), matrix4_i4(1,3))

matrix2_i4 = 5_i4
matrix3_i4 = 3_i4
matrix4_i4 = 2_i4

print*, '(A) non-allocated matrix will be pasted by  matrix2'
print*, '    matrix2 = (/', matrix2_i4(1,:)
print*, '                 ',matrix2_i4(2,:),' /)'
call paste(matrix1_i4, matrix2_i4)
print*, '    result  = (/ ',matrix1_i4(1,:)
print*, '                 ',matrix1_i4(2,:),' /)'
print*, ' '
if ( any(matrix1_i4(1,:) .ne. (/ 5_i4,5_i4 /)) ) isgood=.false.
if ( any(matrix1_i4(2,:) .ne. (/ 5_i4,5_i4 /)) ) isgood=.false.

print*, '(B) matrix1  will be pasted by matrix3'
print*, '    matrix1 = (/ ', matrix1_i4(1,:)
print*, '                 ', matrix1_i4(2,:),' /)'
print*, '    matrix3 = (/ ', matrix3_i4(1,:)
print*, '                 ', matrix3_i4(2,:),' /)'
call paste(matrix1_i4, matrix3_i4)
print*, '    result  = (/ ', matrix1_i4(1,:)
print*, '                 ', matrix1_i4(2,:),' /)'
print*, ' '
if ( any(matrix1_i4(1,:) .ne. (/ 5_i4,5_i4, 3_i4,3_i4,3_i4 /)) ) isgood=.false.
if ( any(matrix1_i4(2,:) .ne. (/ 5_i4,5_i4, 3_i4,3_i4,3_i4 /)) ) isgood=.false.

print*, '(C) matrix4  will be pasted by scalar 5_i4'
print*, '    matrix4 = (/ ', matrix4_i4(1,:),' /)'
call paste(matrix4_i4, 5_i4)
print*, '    result  = (/ ', matrix4_i4(1,:),' /)'
print*, ' '
if ( any(matrix4_i4(1,:) .ne. (/ 2_i4,2_i4,2_i4, 5_i4 /)) ) isgood=.false.

deallocate(matrix1_i4, matrix2_i4, matrix3_i4, matrix4_i4)
#endif

if (isgood) then
   write(*,*) 'mo_append o.k.'
else
   write(*,*) 'mo_append failed!'
endif

END PROGRAM append_test

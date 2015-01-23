PROGRAM append_test

use mo_kind,   only: i4
#ifndef ABSOFT
use mo_append, only: append, paste

implicit none

integer(i4), dimension(:),   allocatable  :: vector1_i4, vector2_i4
integer(i4), dimension(2)                 :: vector3_i4

integer(i4), dimension(:,:), allocatable  :: matrix1_i4, matrix2_i4, matrix3_i4, matrix4_i4

character(256), dimension(:),   allocatable :: vector1_c, vector2_c
character(256), dimension(2)                :: vector3_c
character(256), dimension(:,:), allocatable :: matrix1_c, matrix2_c, matrix3_c, matrix4_c
#endif

logical :: isgood

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

deallocate(vector1_i4)
deallocate(vector2_i4)
deallocate(matrix1_i4, matrix2_i4)

! Character
allocate(vector2_c(4))
allocate(matrix1_c(2,2), matrix2_c(3,2))

vector2_c = 'B'
vector3_c = 'C'

matrix1_c = 'E'
matrix2_c = 'C'

print*, '(A) non-allocated character vector will be appended by character vector2'
print*, '    vector2(1) = ',trim(vector2_c(1))
call append(vector1_c, vector2_c)
print*, '    result(1)  = ',trim(vector1_c(1))
print*, ' '
if (trim(vector1_c(1)) /= 'B') isgood=.false.

print*, '(B) character vector1  will be appended by character vector3'
print*, '    vector3(1) = ',trim(vector3_c(1))
call append(vector1_c, vector3_c)
print*, '    result(4)  = ',trim(vector1_c(4))
print*, '    result(5)  = ',trim(vector1_c(5))
print*, ' '
if (trim(vector1_c(5)) /= 'C') isgood=.false.

print*, '(C) character vector1 will be appended by scalar E'
call append(vector1_c, 'E')
print*, '    result(6)  = ',trim(vector1_c(6))
print*, '    result(7)  = ',trim(vector1_c(7))
print*, ' '
if (trim(vector1_c(7)) /= 'E') isgood=.false.

print*, '(D) character matrix1 will be appended by character matrix 2'
print*, '    matrix1(1,1) = ', trim(matrix1_c(2,2))
print*, '    matrix2(1,1) = ', trim(matrix2_c(1,1))
call append(matrix1_c, matrix2_c)
print*, '    matrix1(2,2) = ', trim(matrix1_c(2,2))
print*, '    matrix1(3,2) = ', trim(matrix1_c(3,2))
print*, ' '
if (trim(matrix1_c(3,2)) /= 'C') isgood=.false.


print*, '(E) character matrix1 will be appended by first column only of character matrix 2'
call append(matrix1_c, matrix2_c(:,:1), fill_value = 'X' )
print*, '    matrix1(7,2) = ', trim(matrix1_c(7,2))

if (trim( matrix1_c(7,2 ) ) .ne. 'X' ) isgood=.false.

deallocate(vector1_c)
deallocate(vector2_c)
deallocate(matrix1_c, matrix2_c)

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

deallocate( matrix2_i4, matrix3_i4 )
allocate(matrix2_i4(2,2), matrix3_i4(2,3) )
matrix2_i4 = 5_i4
matrix3_i4 = 3_i4

print*, '(D) matrix2 will be pasted by the first row of matrix3'
print*, '    matrix2 = (/', matrix2_i4(1,:)
print*, '                ',matrix2_i4(2,:),' /)'
print*, '    matrix3 = (/', matrix3_i4(1,:)
print*, '                ',matrix3_i4(2,:),' /)'
call paste(matrix2_i4, matrix3_i4(:1,:), fill_value = -9999_i4 )
print*, '    result  = (/', matrix2_i4(1,:)
print*, '                ',matrix2_i4(2,:),'/)'
if ( matrix2_i4(2,3) .ne. -9999_i4 ) isgood=.false.

deallocate(matrix1_i4, matrix2_i4, matrix3_i4, matrix4_i4)

! Character
allocate(matrix2_c(2,2), matrix3_c(2,3), matrix4_c(1,3))

matrix2_c = 'E'
matrix3_c = 'C'
matrix4_c = 'B'

print*, '(A) non-allocated character matrix will be pasted by character matrix2'
print*, '    matrix2(1,1) = ', trim(matrix2_c(1,1))
call paste(matrix1_c, matrix2_c)
print*, '    matrix1(1,1) = ', trim(matrix1_c(1,1))
print*, ' '
if (trim(matrix1_c(1,1)) /= 'E') isgood=.false.

print*, '(B) character matrix1 will be pasted by character matrix3'
print*, '    matrix1(2,2) = ', trim(matrix1_c(2,2))
print*, '    matrix3(1,1) = ', trim(matrix3_c(1,1))
call paste(matrix1_c, matrix3_c)
print*, '    matrix1(2,2) = ', trim(matrix1_c(2,2))
print*, '    matrix1(2,3) = ', trim(matrix1_c(2,3))
print*, ' '
if (trim(matrix1_c(2,3)) /= 'C') isgood=.false.

print*, '(C) character matrix4 will be pasted by E'
print*, '    matrix4(1,3) = ', trim(matrix4_c(1,3))
call paste(matrix4_c, 'E')
print*, '    matrix4(1,4) = ', trim(matrix4_c(1,4))
print*, ' '
if (trim(matrix4_c(1,4)) /= 'E') isgood=.false.

deallocate(matrix1_c, matrix2_c, matrix3_c, matrix4_c)
#endif

if (isgood) then
   write(*,*) 'mo_append o.k.'
else
   write(*,*) 'mo_append failed!'
endif

END PROGRAM append_test

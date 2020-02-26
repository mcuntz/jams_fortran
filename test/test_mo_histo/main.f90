PROGRAM histo_test

use mo_kind, only: dp, i4
use mo_ansi_colors, only: color, c_red, c_green
#ifndef __ABSOFT__
use mo_histo, only: histo

implicit none

real(dp), dimension(10)                :: x_dp
real(dp), dimension(5,3)               :: x_dp_2d
real(dp), allocatable, dimension(:)    :: binx_dp
real(dp), allocatable, dimension(:,:)  :: binx_dp_2d
real(dp)                               :: w_dp

logical , dimension(10)                :: maske
integer(i4), allocatable, dimension(:) :: bincount

real(dp),    dimension(:),   allocatable :: Has2Be_binx_dp
real(dp),    dimension(:,:), allocatable :: Has2Be_binx_dp_2d
integer(i4), dimension(:),   allocatable :: Has2Be_bincount

integer(i4)                            :: i,n

print*, '  '
print*, 'YOUR DATA:'
x_dp = (/ 1.2_dp, 1.5_dp, 3.4_dp, 5.3_dp , 6.3_dp, 7.2_dp, 1.8_dp, 5.4_dp, 4.3_dp , 10.3_dp /)
print*, '# Values: ',size(x_dp)
print*, x_dp
print*, '  '
print*, '----------------------------'
print*, '  '

print*, 'Histogram (default)'
call histo(x_dp,binx_dp,bincount,w_dp)

n = size(binx_dp,1)
do i=1,n
   print '(A4,I2,A8,F8.3,A13,I3)', &
         'Bin(',i,'): x = (',binx_dp(i),')   Bincount: ',bincount(i)
end do
print*, '  '

allocate(Has2Be_bincount(3), Has2Be_binx_dp(3))
Has2Be_binx_dp  = (/2.717_dp, 5.750_dp, 8.783_dp/)
Has2Be_bincount = (/4, 5, 1/)
if ( all(bincount .eq. Has2Be_bincount) .and. &
     all(nint(binx_dp*1000.0_dp) .eq. nint(Has2Be_binx_dp*1000.0_dp))) then
   print*, 'mo_histo: default ', color('o.k.', c_green)
else
   print*, 'mo_histo: default ', color('failed', c_red)
end if
deallocate(Has2Be_bincount, Has2Be_binx_dp)

print*, '  '
print*, '----------------------------'
print*, '  '

print*, 'Histogram (bins=7)'
call histo(x_dp,binx_dp,bincount,w_dp,bins=7)

n = size(binx_dp,1)
do i=1,n
   print '(A4,I2,A8,F8.3,A13,I3)', &
         'Bin(',i,'): x = (',binx_dp(i),')   Bincount: ',bincount(i)
end do
print*, '  '

allocate(Has2Be_bincount(n), Has2Be_binx_dp(n))
Has2Be_binx_dp  = (/1.850_dp, 3.150_dp, 4.450_dp, 5.750_dp, 7.050_dp,9.650_dp/)
Has2Be_bincount = (/3,1,1,3,1,1/)
if ( all(bincount .eq. Has2Be_bincount) .and. &
     all(nint(binx_dp*1000.0_dp) .eq. nint(Has2Be_binx_dp*1000.0_dp))) then
   print*, 'mo_histo: fixed number of bins ', color('o.k.', c_green)
else
   print*, 'mo_histo: fixed number of bins ', color('failed', c_red)
end if
deallocate(Has2Be_bincount, Has2Be_binx_dp)

print*, '  '
print*, '----------------------------'
print*, '  '

print*, 'Histogram (binwidth=1.5)'
call histo(x_dp,binx_dp,bincount,w_dp,binwidth=1.5_dp)

n = size(binx_dp,1)
do i=1,n
   print '(A4,I2,A8,F8.3,A13,I3)', &
         'Bin(',i,'): x = (',binx_dp(i),')   Bincount: ',bincount(i)
end do
print*, '  '

allocate(Has2Be_bincount(n), Has2Be_binx_dp(n))
Has2Be_binx_dp  = (/1.950_dp, 3.450_dp, 4.950_dp, 6.450_dp, 7.950_dp,10.950_dp/)
Has2Be_bincount = (/3,1,3,1,1,1/)
if ( all(bincount .eq. Has2Be_bincount) .and. &
     all(nint(binx_dp*1000.0_dp) .eq. nint(Has2Be_binx_dp*1000.0_dp))) then
   print*, 'mo_histo: fixed width of bins ', color('o.k.', c_green)
else
   print*, 'mo_histo: fixed width of bins ', color('failed', c_red)
end if
deallocate(Has2Be_bincount, Has2Be_binx_dp)

print*, '  '
print*, '----------------------------'
print*, '  '

print*, 'Histogram (mask= TTTTT FFFFF)'
maske = .false.
do i=1,5
   maske(i) = .true.
end do

call histo(x_dp,binx_dp,bincount,w_dp,mask=maske)

n = size(binx_dp,1)
do i=1,n
   print '(A4,I2,A8,F8.3,A13,I3)', &
         'Bin(',i,'): x = (',binx_dp(i),')   Bincount: ',bincount(i)
end do
print*, '  '

allocate(Has2Be_bincount(n), Has2Be_binx_dp(n))
Has2Be_binx_dp  = (/2.475_dp, 5.025_dp/)
Has2Be_bincount = (/3,2/)
if ( all(bincount .eq. Has2Be_bincount) .and. &
     all(nint(binx_dp*1000.0_dp) .eq. nint(Has2Be_binx_dp*1000.0_dp))) then
   print*, 'mo_histo: with mask ', color('o.k.', c_green)
else
   print*, 'mo_histo: with mask ', color('failed', c_red)
end if
deallocate(Has2Be_bincount, Has2Be_binx_dp)

print*, '  '
print*, '----------------------------'
print*, '  '

print*, '  '
print*, 'YOUR DATA:'
x_dp_2d(:,1) = (/ 1.2_dp, 1.5_dp, 3.4_dp, 5.3_dp , 6.3_dp /)
x_dp_2d(:,2) = (/ 10.0_dp, 8.0_dp, 4.0_dp, 3.0_dp , 3.0_dp /)
x_dp_2d(:,3) = (/ 7.0_dp, 5.0_dp, 3.0_dp, 3.0_dp , 4.0_dp /)
print*, '# Values:  ( ',size(x_dp_2d,1),' , ',size(x_dp_2d,2),' )'
print*, 'x(:,1) = ',x_dp_2d(:,1)
print*, 'x(:,2) = ',x_dp_2d(:,2)
print*, 'x(:,3) = ',x_dp_2d(:,3)
print*, '  '
print*, '----------------------------'
print*, '  '

print*, 'Histogram (default)'
call histo(x_dp_2d(:,1),binx_dp,bincount,w_dp)

n = size(binx_dp,1)
do i=1,n
   print '(A4,I2,A8,F8.3,A13,I3)', &
         'Bin(',i,'): x = (',binx_dp(i),')   Bincount: ',bincount(i)
end do
print*, '  '

print*, 'Histogram (with averaging)'
call histo(x_dp_2d,binx_dp_2d,bincount,w_dp)

n = size(binx_dp,1)
do i=1,n
   print '(A4,I2,A8,F8.3,A3,F8.3,A3,F8.3,A13,I3)', &
         'Bin(',i,'): x = (',binx_dp_2d(i,1),',',binx_dp_2d(i,2),',',binx_dp_2d(i,3),')   Bincount: ',bincount(i)
end do
print*, '  '

allocate(Has2Be_bincount(2), Has2Be_binx_dp_2d(2,3))
Has2Be_binx_dp_2d(1,:)  = (/2.475_dp, 7.333_dp, 5.000_dp/)
Has2Be_binx_dp_2d(2,:)  = (/5.025_dp, 3.000_dp, 3.500_dp/)
Has2Be_bincount = (/3,2/)
if ( all(bincount .eq. Has2Be_bincount) .and. &
     all(nint(binx_dp_2d*1000.0_dp) .eq. nint(Has2Be_binx_dp_2d*1000.0_dp))) then
#endif
   print*, 'mo_histo: with averaging (default) ', color('o.k.', c_green)
#ifndef __ABSOFT__
else
   print*, 'mo_histo: with averaging (default) ', color('failed', c_red)
end if
deallocate(Has2Be_bincount, Has2Be_binx_dp_2d)

print*, '  '
print*, '----------------------------'
print*, '  '
#endif


END PROGRAM histo_test

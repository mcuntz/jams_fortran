PROGRAM inpoly_test

use mo_kind,   only: dp, i4
use mo_poly, only: inpoly, areapoly, center_of_mass

implicit none

integer(i4)                             :: N
real(dp), dimension(:,:), allocatable   :: coord
integer(i4)                             :: inside
real(dp)                                :: area
real(dp),dimension(2)                   :: com

logical                                 :: isgood

! USAGE: VERY EASY TESTING POLYGON

N=4

allocate (coord(N,2))
coord(:,1) = (/ 1.0_dp,2.0_dp,2.0_dp,1.0_dp /)
coord(:,2) = (/ 1.0_dp,1.0_dp,2.0_dp,2.0_dp /)
isgood = .true.

call inpoly( (/1.5_dp,1.5_dp/) , coord, inside)
if (inside .ne. 1_i4) isgood=.false.

select case (inside)
case(-1)    
   print*, 'THE POINT IS OUTSIDE THE POLYGON'
case(1)     
   print*, 'THE POINT IS INSIDE THE POLYGON'
case(0)     
   print*, 'THE POINT IS ON AN EDGE OR AT A VERTEX'
end select

call inpoly( (/0.5_dp,1.5_dp/) , coord, inside)
if (inside .ne. -1_i4) isgood=.false.

select case (inside)
case(-1)    
   print*, 'THE POINT IS OUTSIDE THE POLYGON'
case(1)     
   print*, 'THE POINT IS INSIDE THE POLYGON'
case(0)     
   print*, 'THE POINT IS ON AN EDGE OR AT A VERTEX'
end select

call inpoly( (/1.5_dp,1.0_dp/) , coord, inside)
if (inside .ne. 0_i4) isgood=.false.

select case (inside)
case(-1)    
   print*, 'THE POINT IS OUTSIDE THE POLYGON'
case(1)     
   print*, 'THE POINT IS INSIDE THE POLYGON'
case(0)     
   print*, 'THE POINT IS ON AN EDGE OR AT A VERTEX'
end select


area = areapoly(coord)
if ( abs(area - 1.0_dp) .gt. tiny(1.0_dp) ) isgood = .false.
print*, 'AREA OF POLYGON IS ', area

com = center_of_mass(coord)
if ( ( abs(com(1)- 1.5_dp) .gt. tiny(1.0_dp) ) .or. ( abs(com(2)- 1.5_dp) .gt. tiny(1.0_dp) ) )  isgood = .false.
print*, 'CENTER OF MASS FOR POLYGON IS ', com




if (isgood) then
   write(*,*) 'mo_poly o.k.'
else
   write(*,*) 'mo_poly failed!'
endif

END PROGRAM inpoly_test

MODULE mo_poly

  ! This module determines wether a 2D point lies inside, outside, or 
  ! on the vertice/edge of a 2D polygon
  ! and is part of the UFZ CHS Fortran library.
  !
  ! Written  Juliane Mai, July 2012

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! The UFZ Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The UFZ Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the UFZ Fortran library. If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2012 Juliane Mai

  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: inpoly          ! Test if 2D point is inside, outside or on vertice/edge of a 2D polygon

  INTERFACE inpoly
     MODULE PROCEDURE inpoly_dp, inpoly_sp
  END INTERFACE inpoly

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         inpoly

  !     PURPOSE
  !         Determines whether a 2D point falls in a polygon or is located outside or 
  !         lies on a vertex or an edge of the polygon.
  !         The polygon can be convex or not.
  !         
  !         The method is only applicable for 2D polygons and points.
  !
  !         The original version of the source code (pnpoly) was implemented by 
  !         W. Randolph Franklin. It was insufficiently assigning vertex/edge points.

  !     CALLING SEQUENCE
  !         call inpoly( p, coord, inside )
  
  !     INDENT(IN)
  !         real(sp/dp) :: p(2)           coordinates of the point in question
  !         real(sp/dp) :: coord(:,2)     (x,y) coordinates of edges of the polygon    

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         integer(i4) :: inside   whether point is inside (=1), outside (=-1) or
  !                                 on a vertex/edge of the polygon

  !     INDENT(IN), OPTIONAL
  !         None  

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Only available in 2D version 

  !     EXAMPLE
  !         polygon(:,1) = (/ 1.0_dp,2.0_dp,2.0_dp,1.0_dp /)
  !         polygon(:,2) = (/ 1.0_dp,1.0_dp,2.0_dp,2.0_dp /)
  !
  !         point = (/ 1.5, 1.5 /)
  !
  !         call inpoly( point, polygon, inside )
  !
  !         --> inside = 1 ... point is inside the polygon
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html

  !     HISTORY
  !         Written,  Juliane Mai, July 2012

SUBROUTINE inpoly_dp(P,coord,erg) 
!
  REAL(dp),   dimension(2),   INTENT(IN)  :: P           ! Point in question
  REAL(dp),   dimension(:,:), INTENT(IN)  :: coord       ! Coordinates of the polygon
  INTEGER(i4),                INTENT(OUT) :: erg         ! Result: 
                                                         !     inside:         erg =  1
                                                         !     outside:        erg = -1
                                                         !     on vertex/edge: erg =  0
  
  ! Local variables
  REAL(dp), dimension(:), allocatable   :: X, Y 
  REAL(dp)                              :: lx, ly
  LOGICAL                               :: MX,MY,NX,NY 
  INTEGER(i4)                           :: N, I, J

  N  = size(coord,1)
  allocate(x(N))
  allocate(y(N))

  do i=1,N
     X(I)=coord(I,1)-P(1)
     Y(I)=coord(I,2)-P(2)
     ! Edge
     if ( (X(I) .EQ. 0.0_dp) .and. (Y(I) .EQ. 0.0_dp) ) then
        erg=0_i4
        return
     end if
  end do

  erg=-1_i4
  DO 2 I=1,N 
     J=1+MOD(I,N) 
     ! vertical Vertex
     If ( ( coord(I,1) .eq. coord(J,1) ) .and. ( coord(I,1) .eq. P(1) ) ) then
        ly = (P(2)-coord(J,2)) / (coord(I,2)-coord(J,2))
        If ( ( ly .ge. 0.0_dp ) .and. ( ly .le. 1.0_dp) ) then
           erg=0_i4
           return
        end if
     end if
     ! horizontal Vertex
     If ( ( coord(I,2) .eq. coord(J,2) ) .and. ( coord(I,2) .eq. P(2) ) ) then
        lx = (P(1)-coord(J,1)) / (coord(I,1)-coord(J,1))
        If ( ( lx .ge. 0.0_dp ) .and. ( lx .le. 1.0_dp) ) then
           erg=0_i4
           return
        end if
     end if
     ! 
     MX=X(I).GE.0.0_dp 
     NX=X(J).GE.0.0_dp 
     MY=Y(I).GE.0.0_dp 
     NY=Y(J).GE.0.0_dp 
     
     IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2 
     IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3 
     erg=-erg 
     GO TO 2 
3    IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5 
4    erg=0_i4 
     RETURN 
5    erg=-erg 
2    CONTINUE

     deallocate(x)
     deallocate(y)

     RETURN
 
END SUBROUTINE inpoly_dp

SUBROUTINE inpoly_sp(P,coord,erg) 
!
  REAL(sp),   dimension(2),   INTENT(IN)  :: P           ! Point in question
  REAL(sp),   dimension(:,:), INTENT(IN)  :: coord       ! Coordinates of the polygon
  INTEGER(i4),                INTENT(OUT) :: erg         ! Result: 
                                                         !     inside:         erg =  1
                                                         !     outside:        erg = -1
                                                         !     on vertex/edge: erg =  0
  
  ! Local variables
  REAL(sp), dimension(:), allocatable   :: X, Y 
  REAL(sp)                              :: lx, ly
  LOGICAL                               :: MX,MY,NX,NY 
  INTEGER(i4)                           :: N, I, J

  N  = size(coord,1)
  allocate(x(N))
  allocate(y(N))

  do i=1,N
     X(I)=coord(I,1)-P(1)
     Y(I)=coord(I,2)-P(2)
     ! Edge
     if ( (X(I) .EQ. 0.0_sp) .and. (Y(I) .EQ. 0.0_sp) ) then
        erg=0_i4
        return
     end if
  end do

  erg=-1_i4
  DO 2 I=1,N 
     J=1+MOD(I,N) 
     ! vertical Vertex
     If ( ( coord(I,1) .eq. coord(J,1) ) .and. ( coord(I,1) .eq. P(1) ) ) then
        ly = (P(2)-coord(J,2)) / (coord(I,2)-coord(J,2))
        If ( ( ly .ge. 0.0_sp ) .and. ( ly .le. 1.0_sp) ) then
           erg=0_i4
           return
        end if
     end if
     ! horizontal Vertex
     If ( ( coord(I,2) .eq. coord(J,2) ) .and. ( coord(I,2) .eq. P(2) ) ) then
        lx = (P(1)-coord(J,1)) / (coord(I,1)-coord(J,1))
        If ( ( lx .ge. 0.0_sp ) .and. ( lx .le. 1.0_sp) ) then
           erg=0_i4
           return
        end if
     end if
     ! 
     MX=X(I).GE.0.0_sp 
     NX=X(J).GE.0.0_sp 
     MY=Y(I).GE.0.0_sp 
     NY=Y(J).GE.0.0_sp 
     
     IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2 
     IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3 
     erg=-erg 
     GO TO 2 
3    IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5 
4    erg=0_i4 
     RETURN 
5    erg=-erg 
2    CONTINUE

     deallocate(x)
     deallocate(y)

     RETURN
 
END SUBROUTINE inpoly_sp

END MODULE mo_poly

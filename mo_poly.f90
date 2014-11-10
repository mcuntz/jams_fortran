MODULE mo_poly

  ! This module determines wether a 2D point lies inside, outside, or
  ! on the vertice/edge of a 2D polygon
  ! and is part of the UFZ CHS Fortran library.
  !
  ! Written  Juliane Mai, July 2012
  ! Modified Maren Goehler, July 2012 - area & center of mass

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
  ! along with the UFZ Fortran library (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2012 Juliane Mai


  USE mo_kind, ONLY: i4, sp, dp
  USE mo_utils, only: eq, ge, le

  IMPLICIT NONE

  PUBLIC :: areapoly        ! compute the area of a polygon
  PUBLIC :: center_of_mass  ! compute the center of mass of a polygon
  PUBLIC :: inpoly          ! Test if 2D point is inside, outside or&
  ! on vertice/edge of a 2D polygon

  ! ------------------------------------------------------------------

  !     NAME
  !         areapoly

  !     PURPOSE
  !         Function for computing the area of a polygon
  !         The polygon can be convex or not.
  !
  !         The method is only applicable for 2D polygons and points.
  !
  !         http://de.wikipedia.org/wiki/Geometrischer_Schwerpunkt
  !
  !     CALLING SEQUENCE
  !         area =  areapoly(coord)

  !     INTENT(IN)
  !         real(sp/dp) :: coord(:,2)     (x,y) coordinates of edges of the polygon

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Only available in 2D version

  !     EXAMPLE
  !         polygon(:,1) = (/ 1.0_dp,2.0_dp,2.0_dp,1.0_dp /)
  !         polygon(:,2) = (/ 1.0_dp,1.0_dp,2.0_dp,2.0_dp /)
  !
  !         area = areapoly(polygon )
  !
  !         --> area = 1
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         http://de.wikipedia.org/wiki/Geometrischer_Schwerpunkt

  !     HISTORY
  !         Written,  Maren Goehler, July 2012
  INTERFACE areapoly
     MODULE PROCEDURE areapoly_sp, areapoly_dp
  END INTERFACE areapoly

  ! ------------------------------------------------------------------

  !     NAME
  !         center_of_mass

  !     PURPOSE
  !         Function for computing the center of mass of a polygon
  !         Computation of polygon area needed for center of mass.
  !         The polygon can be convex or not.
  !
  !         The method is only applicable for 2D polygons and points.
  !
  !         http://de.wikipedia.org/wiki/Geometrischer_Schwerpunkt
  !
  !         A = sum(xiyi+1-xi+1yi)
  !         xs = 1/(6*A) * sum(xi+xi+1)(xiyi+1-xi+1yi)
  !         ys = 1/(6*A) * sum(yi+yi+1)(xiyi+1-xi+1yi)

  !     CALLING SEQUENCE
  !         com =  center_of_mass(coord)

  !     INTENT(IN)
  !         real(sp/dp) :: coord(:,2)     (x,y) coordinates of edges of the polygon

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Only available in 2D version

  !     EXAMPLE
  !         polygon(:,1) = (/ 1.0_dp,2.0_dp,2.0_dp,1.0_dp /)
  !         polygon(:,2) = (/ 1.0_dp,1.0_dp,2.0_dp,2.0_dp /)
  !
  !         com = center_of_mass(polygon)
  !
  !         --> com = (/1.5_dp, 1.5_dp/)
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         http://de.wikipedia.org/wiki/Geometrischer_Schwerpunkt

  !     HISTORY
  !         Written,  Maren Goehler, July 2012
  INTERFACE center_of_mass
     MODULE PROCEDURE  center_of_mass_sp,center_of_mass_dp
  END INTERFACE center_of_mass

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

  !     INTENT(IN)
  !         real(sp/dp) :: p(2)           coordinates of the point in question
  !         real(sp/dp) :: coord(:,2)     (x,y) coordinates of edges of the polygon

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         integer(i4) :: inside   whether point is inside (=1), outside (=-1) or
  !                                 on a vertex/edge of the polygon

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
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
  !         Modified, Juliane Mai, July 2012 - removing GOTO statements
  INTERFACE inpoly
     MODULE PROCEDURE inpoly_dp, inpoly_sp
  END INTERFACE inpoly

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------


  FUNCTION areapoly_dp(coord)
    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),     INTENT(IN)     :: coord        ! Coordinates of the polygon in question
    REAL(dp)                                     :: areapoly_dp

    ! Local variables
    INTEGER(i4)                                  :: i,k          ! loop
    INTEGER(i4)                                  :: Nedges       ! number of coordinates
    REAL(dp)                                     :: xsum         ! for summing up
    REAL(dp)                                     :: ysum         ! for summing up

    xsum   = 0.0_dp
    ysum   = 0.0_dp
    Nedges = size(coord,1)

    do i = 1,  Nedges
       if (i .eq. Nedges) then
          k = 1_i4
       else
          k = i + 1_i4
       end if
       xsum = xsum + ( coord(i,1) * coord(k,2) )
       ysum = ysum + ( coord(i,2) * coord(k,1) )
    end do

    areapoly_dp = 0.5_sp * (xsum - ysum)

    RETURN
  END FUNCTION areapoly_dp


  FUNCTION areapoly_sp(coord)
    IMPLICIT NONE
    REAL(sp),   dimension(:,:),   INTENT(IN)    :: coord       ! Coordinates of the polygon in question
    REAL(sp)                                    :: areapoly_sp

    ! Local variables
    INTEGER(i4)                                 :: i,k         ! loop
    INTEGER(i4)                                 :: Nedges      ! number of coordinates
    REAL(sp)                                    :: xsum        ! for summing up
    REAL(sp)                                    :: ysum        ! for summing up

    xsum   = 0.0_sp
    ysum   = 0.0_sp
    Nedges = size(coord,1)

    do i = 1,  Nedges
       if (i .eq. Nedges) then
          k = 1_i4
       else
          k = i + 1_i4
       end if
       xsum = xsum + ( coord(i,1) * coord(k,2) )
       ysum = ysum + ( coord(i,2) * coord(k,1) )
    end do


    areapoly_sp = 0.5_sp * (xsum - ysum)

  END FUNCTION areapoly_sp

  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------

  FUNCTION center_of_mass_dp(coord)
    IMPLICIT NONE

    REAL(dp),   dimension(:,:),   INTENT(IN)    :: coord              ! Coordinates of polygon in question
    REAL(dp),dimension(2)                       :: center_of_mass_dp

    ! Local variables
    INTEGER(i4)                                 :: i,k       ! loop
    INTEGER(i4)                                 :: Nedges    ! number of coordinates
    REAL(dp)                                    :: area      ! area of the polygon
    REAL(dp)                                    :: xsum      ! for summing up
    REAL(dp)                                    :: ysum      ! for summing up
    REAL(dp)                                    :: hotspot_x ! for summing up
    REAL(dp)                                    :: hotspot_y ! for summing up

    xsum = 0.0_dp
    ysum = 0.0_dp
    Nedges = size(coord,1)

    area = areapoly_dp(coord)

    do i = 1, Nedges
       if (i .eq. Nedges ) then
          k = 1_i4
       else
          k = i + 1_i4
       end if
       ! multiply x coord by the y coord of next vertex
       xsum = xsum + ((coord(i,1) + coord(k,1)) * &
            ((coord(i,1) * coord(k,2) - coord(k,1) * coord(i,2))))

       ysum = ysum + ((coord(i,2) + coord(k,2)) * &
            ((coord(i,1) * coord(k,2) - coord(k,1) * coord(i,2))))
    end do

    hotspot_x = 1.0_dp / (6.0_dp * area) * xsum
    hotspot_y = 1.0_dp / (6.0_dp * area) * ysum

    center_of_mass_dp(1) = hotspot_x
    center_of_mass_dp(2) = hotspot_y

    RETURN
  END FUNCTION  center_of_mass_dp


  FUNCTION center_of_mass_sp(coord)
    IMPLICIT NONE

    REAL(sp),dimension(:,:),   INTENT(IN)       :: coord               ! Coordinates of polygon in question
    REAL(sp),dimension(2)                       :: center_of_mass_sp

    ! Local variables
    INTEGER(i4)                                 :: i,k       ! loop
    INTEGER(i4)                                 :: Nedges    ! number of coordinates
    REAL(sp)                                    :: area      ! area of the polygon
    REAL(sp)                                    :: xsum      ! for summing up
    REAL(sp)                                    :: ysum      ! for summing up
    REAL(sp)                                    :: hotspot_x ! for summing up
    REAL(sp)                                    :: hotspot_y ! for summing up

    xsum   = 0.0_sp
    ysum   = 0.0_sp
    Nedges = size(coord,1)

    area = areapoly_sp(coord)

    do i = 1, Nedges
       if (i .eq. Nedges ) then
          k = 1_i4
       else
          k = i + 1_i4
       end if
       ! multiply x coord by the y coord of next vertex
       xsum = xsum + ((coord(i,1) + coord(k,1)) * &
            ((coord(i,1) * coord(k,2) - coord(k,1) * coord(i,2))))

       ysum = ysum + ((coord(i,2) + coord(k,2)) * &
            ((coord(i,1) * coord(k,2) - coord(k,1) * coord(i,2))))
    end do

    hotspot_x = 1.0_sp / (6.0_sp * area) * xsum
    hotspot_y = 1.0_sp / (6.0_sp * area) * ysum

    center_of_mass_sp(1) = hotspot_x
    center_of_mass_sp(2) = hotspot_y

    RETURN
  END FUNCTION  center_of_mass_sp

  ! ------------------------------------------------------------------

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
    LOGICAL                               :: MX,MY,NX,NY, test1, test2
    INTEGER(i4)                           :: N, I, J

    N  = size(coord,1)
    allocate(x(N))
    allocate(y(N))

    do i=1,N
       X(I)=coord(I,1)-P(1)
       Y(I)=coord(I,2)-P(2)
       ! Edge
       if ( eq(X(I),0.0_dp) .and. eq(Y(I),0.0_dp) ) then
          erg=0_i4
          return
       end if
    end do

    erg=-1_i4

    DO I=1,N
       J=1+MOD(I,N)
       ! vertical Vertex
       If ( eq(coord(I,1),coord(J,1)) .and. eq(coord(I,1),P(1)) ) then
          ly = (P(2)-coord(J,2)) / (coord(I,2)-coord(J,2))
          If ( ge(ly,0.0_dp) .and. le(ly,1.0_dp) ) then
             erg=0_i4
             return
          end if
       end if
       ! horizontal Vertex
       If ( eq(coord(I,2),coord(J,2)) .and. eq(coord(I,2),P(2)) ) then
          lx = (P(1)-coord(J,1)) / (coord(I,1)-coord(J,1))
          If ( ge(lx,0.0_dp ) .and. le(lx,1.0_dp) ) then
             erg=0_i4
             return
          end if
       end if
       !
       MX = ge(X(I),0.0_dp)
       NX = ge(X(J),0.0_dp)
       MY = ge(Y(I),0.0_dp)
       NY = ge(Y(J),0.0_dp)

       test1 = .NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)
       test2 = .NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))

       if (.not. test1) then
          if (test2) then
             IF ((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I)) .lt. 0.0_dp) then
                cycle
             ELSE
                IF ((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I)) .gt. 0.0_dp) then
                   erg = -erg
                   cycle
                ELSE
                   erg = 0_i4
                   return
                END IF
             END IF
          else
             erg=-erg
          end if
       end if

    END DO

    deallocate(x)
    deallocate(y)

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
    LOGICAL                               :: MX,MY,NX,NY, test1, test2
    INTEGER(i4)                           :: N, I, J

    N  = size(coord,1)
    allocate(x(N))
    allocate(y(N))

    do i=1,N
       X(I)=coord(I,1)-P(1)
       Y(I)=coord(I,2)-P(2)
       ! Edge
       if ( eq(X(I),0.0_sp) .and. eq(Y(I),0.0_sp) ) then
          erg=0_i4
          return
       end if
    end do

    erg=-1_i4

    DO I=1,N
       J=1+MOD(I,N)
       ! vertical Vertex
       If ( eq(coord(I,1),coord(J,1)) .and. eq(coord(I,1),P(1)) ) then
          ly = (P(2)-coord(J,2)) / (coord(I,2)-coord(J,2))
          If ( ge(ly,0.0_sp) .and. le(ly,1.0_sp) ) then
             erg=0_i4
             return
          end if
       end if
       ! horizontal Vertex
       If ( eq(coord(I,2),coord(J,2)) .and. eq(coord(I,2),P(2)) ) then
          lx = (P(1)-coord(J,1)) / (coord(I,1)-coord(J,1))
          If ( ge(lx,0.0_sp ) .and. le(lx,1.0_sp) ) then
             erg=0_i4
             return
          end if
       end if
       !
       MX = ge(X(I),0.0_sp)
       NX = ge(X(J),0.0_sp)
       MY = ge(Y(I),0.0_sp)
       NY = ge(Y(J),0.0_sp)

       test1 = .NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)
       test2 = .NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))

       if (.not. test1) then
          if (test2) then
             IF ((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I)) .lt. 0.0_sp) then
                cycle
             ELSE
                IF ((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I)) .gt. 0.0_sp) then
                   erg = -erg
                   cycle
                ELSE
                   erg = 0_i4
                   return
                END IF
             END IF
          else
             erg=-erg
          end if
       end if

    END DO

    deallocate(x)
    deallocate(y)

  END SUBROUTINE inpoly_sp

END MODULE mo_poly

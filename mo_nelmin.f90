MODULE mo_nelmin


  ! This module provides NELMIN, which minimizes a function using the Nelder-Mead algorithm
  ! with the Applied Statistics algorithms No. 047.

  ! Written  Matthias Cuntz, Jul 2012 - extended asa047 of John Burkardt

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

  ! Copyright 2012 Matthias Cuntz

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nelmin          ! minimizes a function using the Nelder-Mead algorithm
  PUBLIC :: nelminxy        ! same as nelmin but pass x and y to function

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         nelmin
  !         nelminxy

  !     PURPOSE
  !         Minimizes a user-specified function using the Nelder-Mead algorithm.
  !
  !         Simplex function minimisation procedure due to Nelder and Mead (1965),
  !         as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
  !         subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
  !         25, 97) and Hill(1978, 27, 380-2)
  !
  !         The function to be minimized is the first argument of nelmin and
  !         must be defined as
  !           FUNCTION func(p)
  !             USE mo_kind, ONLY: dp
  !             IMPLICIT NONE
  !             REAL(dp), DIMENSION(:),           INTENT(IN) :: p
  !             REAL(dp) :: func
  !           END FUNCTION func
  !
  !         and for nelminxy
  !           FUNCTION func(p,x,y)
  !             USE mo_kind, ONLY: dp
  !             IMPLICIT NONE
  !             REAL(dp), DIMENSION(:), INTENT(IN) :: p
  !             REAL(dp), DIMENSION(:), INTENT(IN) :: x
  !             REAL(dp), DIMENSION(:), INTENT(IN) :: y
  !             REAL(dp) :: func
  !           END FUNCTION func
  !
  !         This routine does not include a termination test using the
  !         fitting of a quadratic surface.

  !     CALLING SEQUENCE
  !         pmin = nelmin(func, pstart, funcmin, varmin, step, konvge, maxeval, &
  !                       neval, numrestart, ierror)
  !         pmin = nelminxy(func, pstart, xx, yy, funcmin, varmin, step, konvge, maxeval, &
  !                         neval, numrestart, ierror)

  !     INDENT(IN)
  !         real(dp) :: func(p,xx,yy)    Function on which to search the minimum
  !         real(dp) :: pstart(:)        Starting point for the iteration.
  !         real(dp) :: xx               First values to pass as function arguments
  !         real(dp) :: yy               Second values to pass as function arguments

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(dp) :: nelmin(size(start))    the coordinates of the point which is estimated to minimize the function.

  !     INDENT(IN), OPTIONAL
  !         real(dp) :: varmin           the terminating limit for the variance
  !                                      of the function values. varmin>0 is required.
  !         real(dp) :: step(:)          determines the size and shape of the initial simplex. 
  !                                      The relative magnitudes of its elements should reflect
  !                                      the units of the variables. size(step)=size(start)
  !         integer(i4) :: konvge        the convergence check is carried out every konvge iterations.
  !                                      konvge>0 is required.
  !         integer(i4) :: maxeval       the maximum number of function evaluations.

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         real(dp)    :: funcmin             the minimum value of the function.
  !         integer(i4) :: neval               the number of function evaluations used.
  !         integer(i4) :: numrestart          the number of restarts.
  !         integer(i4) :: ierror              error indicator.
  !                                            0: no errors detected.
  !                                            1: varmin or konvge have an illegal value.
  !                                            2: iteration terminated because maxeval was exceeded without convergence.

  !     RESTRICTIONS
  !         None.

  !     EXAMPLE
  !         pstart(1:n) = (/ -1.2_dp, 1.0_dp /)
  !         step(1:n) = (/ 1.0_dp, 1.0_dp /)
  !         pmin = nelmin(rosenbrock, pstart, step=step, ierror=ierror)
  !         -> see also example in test directory

  !     LITERATURE
  !         Simplex function minimisation procedure due to Nelder and Mead (1965),
  !         as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
  !         subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
  !         25, 97) and Hill(1978, 27, 380-2)

  !         John Nelder, Roger Mead, A simplex method for function minimization,
  !         Computer Journal, Volume 7, 1965, pages 308-313.

  !         R ONeill, Algorithm AS 47: Function Minimization Using a Simplex Procedure,
  !         Applied Statistics, Volume 20, Number 3, 1971, pages 338-345.

  !     HISTORY
  !         Original FORTRAN77 version by R ONeill.
  !         FORTRAN90 version by John Burkardt.
  !         Modified,  Matthias Cuntz, Jul 2012 - i4, dp, intent
  !                                             - function, optional
  !                                             - nelimxy

  function nelmin(func, pstart, funcmin, varmin, step, konvge, maxeval, &
       neval, numrestart, ierror)

    implicit none

    INTERFACE
       FUNCTION func(pp)
         USE mo_kind, ONLY: dp
         IMPLICIT NONE
         REAL(dp), DIMENSION(:), INTENT(IN) :: pp
         REAL(dp) :: func
       END FUNCTION func
    END INTERFACE
    real(dp),              intent(IN)  :: pstart(:)
    real(dp),    optional, intent(OUT) :: funcmin
    real(dp),    optional, intent(IN)  :: varmin
    real(dp),    optional, intent(IN)  :: step(:)
    integer(i4), optional, intent(IN)  :: konvge
    integer(i4), optional, intent(IN)  :: maxeval
    integer(i4), optional, intent(OUT) :: neval
    integer(i4), optional, intent(OUT) :: numrestart
    integer(i4), optional, intent(OUT) :: ierror
    real(dp) :: nelmin(size(pstart))

    real(dp), parameter :: ccoeff = 0.5_dp
    real(dp), parameter :: ecoeff = 2.0_dp
    real(dp), parameter :: rcoeff = 1.0_dp
    real(dp), parameter :: eps    = 0.001_dp
    real(dp) :: ipstart(size(pstart))
    real(dp) :: ifuncmin
    real(dp) :: ivarmin
    real(dp) :: istep(size(pstart))
    integer(i4) :: ikonvge
    integer(i4) :: imaxeval
    integer(i4) :: ineval
    integer(i4) :: inumrestart
    integer(i4) :: iierror
    integer(i4) :: n, nn
    integer(i4) :: i
    integer(i4) :: ihi
    integer(i4) :: ilo
    integer(i4) :: j
    integer(i4) :: jcount
    integer(i4) :: l
    real(dp) :: del
    real(dp) :: p(size(pstart),size(pstart)+1)
    real(dp) :: p2star(size(pstart))
    real(dp) :: pbar(size(pstart))
    real(dp) :: pstar(size(pstart))
    real(dp) :: rq
    real(dp) :: x
    real(dp) :: y(size(pstart)+1)
    real(dp) :: y2star
    real(dp) :: ylo
    real(dp) :: ystar
    real(dp) :: z
    real(dp) :: dn, dnn
    real(dp) :: p0(size(pstart)), y0
    !
    ! Defaults
    !
    nelmin(:) = 0.
    if (present(varmin)) then
       if (varmin <= 0.0_dp) stop 'Error nelim: varmin<0'
       ivarmin = varmin
    else
       ivarmin = 1.0e-9_dp
    endif
    if (present(maxeval)) then
       if (maxeval <= 1) stop 'Error nelim: maxeval<=1'
       imaxeval = maxeval
    else
       imaxeval = 1000
    endif
    if (present(konvge)) then
       ikonvge = konvge
    else
       ikonvge = imaxeval / 10
    endif
    if (ikonvge < 1) stop 'Error nelim: konvg<1'
    !
    if (present(step)) then
       istep = step
    else ! if not given, deviate initial by 1%
       y0 = func(pstart)
       do i=1, size(pstart)
          p0 = pstart
          if (p0(i) /= 0.0_dp) then
             p0(i) = 1.01_dp*p0(i)
          else
             p0(i) = 0.01_dp
          endif
          istep(i) = abs(func(p0) - y0)
       enddo
    endif
    !
    !  Check the input parameters.
    !
    !
    !  Initialization.
    !
    ipstart = pstart
    n       = size(ipstart)
    ineval  = 0_i4
    inumrestart = 0_i4
    jcount = ikonvge
    dn     = real(n, dp)
    nn     = n + 1
    dnn    = real(nn, dp)
    del    = 1.0_dp
    rq     = ivarmin * dn
    !
    !  Initial or restarted loop.
    !
    do
       p(1:n,nn) = ipstart(1:n)
       y(nn)     = func(ipstart)
       ineval     = ineval + 1
       !
       !  Define the initial simplex.
       !
       do j = 1, n
          x        = ipstart(j)
          ipstart(j) = ipstart(j) + istep(j) * del
          p(1:n,j) = ipstart(1:n)
          y(j)     = func(ipstart)
          ineval    = ineval + 1
          ipstart(j) = x
       end do
       !
       !  Find highest and lowest Y values.  FUNCMIN = Y(IHI) indicates
       !  the vertex of the simplex to be replaced.
       !
       ilo = minloc(y(1:n+1), 1)
       ylo = y(ilo)
       !
       !  Inner loop.
       !
       do while ( ineval < imaxeval )
          !
          !  FUNCMIN is, of course, the HIGHEST value???
          !
          ihi     = maxloc(y(1:n+1), 1)
          ifuncmin = y(ihi)
          !
          !  Calculate PBAR, the centroid of the simplex vertices
          !  excepting the vertex with Y value FUNCMIN.
          !
          do i = 1, n
             pbar(i) = ( sum(p(i,1:n+1)) - p(i,ihi) ) / dn
          end do
          !
          !  Reflection through the centroid.
          !
          pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
          ystar      = func(pstar)
          ineval      = ineval + 1
          !
          !  Successful reflection, so extension.
          !
          if ( ystar < ylo ) then

             p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
             y2star      = func(p2star)
             ineval       = ineval + 1
             !
             !  Retain extension or contraction.
             !
             if ( ystar < y2star ) then
                p(1:n,ihi) = pstar(1:n)
                y(ihi)     = ystar
             else
                p(1:n,ihi) = p2star(1:n)
                y(ihi)     = y2star
             end if
             !
             !  No extension.
             !
          else
             l = 0
             do i = 1, nn
                if ( ystar < y(i) ) l = l + 1
             end do

             if ( 1 < l ) then
                p(1:n,ihi) = pstar(1:n)
                y(ihi)     = ystar
                !
                !  Contraction on the Y(IHI) side of the centroid.
                !
             else if ( l == 0 ) then
                p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
                y2star      = func(p2star)
                ineval       = ineval + 1
                !
                !  Contract the whole simplex.
                !
                if ( y(ihi) < y2star ) then
                   do j = 1, n + 1
                      p(1:n,j)  = ( p(1:n,j) + p(1:n,ilo) ) * 0.5_dp
                      nelmin(1:n) = p(1:n,j)
                      y(j)      = func(nelmin)
                      ineval     = ineval + 1
                   end do
                   ilo = minloc(y(1:n+1), 1)
                   ylo = y(ilo)

                   cycle
                   !
                   !  Retain contraction.
                   !
                else
                   p(1:n,ihi) = p2star(1:n)
                   y(ihi)     = y2star
                end if
                !
                !  Contraction on the reflection side of the centroid.
                !
             else if ( l == 1 ) then

                p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
                y2star      = func(p2star)
                ineval       = ineval + 1
                !
                !  Retain reflection?
                !
                if ( y2star <= ystar ) then
                   p(1:n,ihi) = p2star(1:n)
                   y(ihi)     = y2star
                else
                   p(1:n,ihi) = pstar(1:n)
                   y(ihi)     = ystar
                end if

             end if

          end if
          !
          !  Check if YLO improved.
          !
          if ( y(ihi) < ylo ) then
             ylo = y(ihi)
             ilo = ihi
          end if

          jcount = jcount - 1

          if ( 0 < jcount ) cycle
          !
          !  Check to see if minimum reached.
          !
          if ( ineval <= imaxeval ) then
             jcount = ikonvge
             x = sum ( y(1:n+1) ) / dnn
             z = sum ( ( y(1:n+1) - x )**2 )
             if ( z <= rq ) exit
          end if

       end do
       !
       !  Factorial tests to check that FUNCMIN is a local minimum.
       !
       nelmin(1:n) = p(1:n,ilo)
       ifuncmin   = y(ilo)

       if ( imaxeval < ineval ) then
          iierror = 2
          exit
       end if

       iierror = 0

       do i = 1, n
          del     = istep(i) * eps
          nelmin(i) = nelmin(i) + del
          z       = func(nelmin)
          ineval  = ineval + 1
          if ( z < ifuncmin ) then
             iierror = 2
             exit
          end if
          nelmin(i) = nelmin(i) - del - del
          z       = func(nelmin)
          ineval  = ineval + 1
          if ( z < ifuncmin ) then
             iierror = 2
             exit
          end if
          nelmin(i) = nelmin(i) + del
       end do

       if ( iierror == 0 ) exit
       !
       !  Restart the procedure.
       !
       ipstart(1:n) = nelmin(1:n)
       del          = eps
       inumrestart  = inumrestart + 1

    end do

    if (present(funcmin)) then
       funcmin = ifuncmin
    endif
    if (present(neval)) then
       neval = ineval
    endif
    if (present(numrestart)) then
       numrestart = inumrestart
    endif
    if (present(ierror)) then
       ierror = iierror
    endif

  end function nelmin


  function nelminxy(func, pstart, xx, yy, funcmin, varmin, step, konvge, maxeval, &
       neval, numrestart, ierror)

    implicit none

    INTERFACE
       FUNCTION func(pp, xxx, yyy)
         USE mo_kind, ONLY: dp
         IMPLICIT NONE
         REAL(dp), DIMENSION(:), INTENT(IN) :: pp
         REAL(dp), DIMENSION(:), INTENT(IN) :: xxx
         REAL(dp), DIMENSION(:), INTENT(IN) :: yyy
         REAL(dp) :: func
       END FUNCTION func
    END INTERFACE
    real(dp),              intent(IN)  :: pstart(:)
    real(dp),              intent(IN)  :: xx(:)
    real(dp),              intent(IN)  :: yy(:)
    real(dp),    optional, intent(OUT) :: funcmin
    real(dp),    optional, intent(IN)  :: varmin
    real(dp),    optional, intent(IN)  :: step(:)
    integer(i4), optional, intent(IN)  :: konvge
    integer(i4), optional, intent(IN)  :: maxeval
    integer(i4), optional, intent(OUT) :: neval
    integer(i4), optional, intent(OUT) :: numrestart
    integer(i4), optional, intent(OUT) :: ierror
    real(dp) :: nelminxy(size(pstart))

    real(dp), parameter :: ccoeff = 0.5_dp
    real(dp), parameter :: ecoeff = 2.0_dp
    real(dp), parameter :: rcoeff = 1.0_dp
    real(dp), parameter :: eps    = 0.001_dp
    real(dp) :: ipstart(size(pstart))
    real(dp) :: ifuncmin
    real(dp) :: ivarmin
    real(dp) :: istep(size(pstart))
    integer(i4) :: ikonvge
    integer(i4) :: imaxeval
    integer(i4) :: ineval
    integer(i4) :: inumrestart
    integer(i4) :: iierror
    integer(i4) :: n, nn
    integer(i4) :: i
    integer(i4) :: ihi
    integer(i4) :: ilo
    integer(i4) :: j
    integer(i4) :: jcount
    integer(i4) :: l
    real(dp) :: del
    real(dp) :: p(size(pstart),size(pstart)+1)
    real(dp) :: p2star(size(pstart))
    real(dp) :: pbar(size(pstart))
    real(dp) :: pstar(size(pstart))
    real(dp) :: rq
    real(dp) :: x
    real(dp) :: y(size(pstart)+1)
    real(dp) :: y2star
    real(dp) :: ylo
    real(dp) :: ystar
    real(dp) :: z
    real(dp) :: dn, dnn
    real(dp) :: p0(size(pstart)), y0
    !
    ! Defaults
    !
    nelminxy(:) = 0.
    if (present(varmin)) then
       if (varmin <= 0.0_dp) stop 'Error nelimxy: varmin<0'
       ivarmin = varmin
    else
       ivarmin = 1.0e-9_dp
    endif
    if (present(maxeval)) then
       if (maxeval <= 1) stop 'Error nelimxy: maxeval<=1'
       imaxeval = maxeval
    else
       imaxeval = 1000
    endif
    if (present(konvge)) then
       ikonvge = konvge
    else
       ikonvge = imaxeval / 10
    endif
    if (ikonvge < 1) stop 'Error nelimxy: konvg<1'
    if (size(xx) /= size(yy)) stop 'Error nelimxy: size(xx) /= size(yy)'
    !
    if (present(step)) then
       istep = step
    else ! if not given, deviate initial by 1%
       y0 = func(pstart, xx, yy)
       do i=1, size(pstart)
          p0 = pstart
          if (p0(i) /= 0.0_dp) then
             p0(i) = 1.01_dp*p0(i)
          else
             p0(i) = 0.01_dp
          endif
          istep(i) = abs(func(p0, xx, yy) - y0)
       enddo
    endif
    !
    !  Check the input parameters.
    !
    !
    !  Initialization.
    !
    ipstart = pstart
    n       = size(ipstart)
    ineval  = 0_i4
    inumrestart = 0_i4
    jcount = ikonvge
    dn     = real(n, dp)
    nn     = n + 1
    dnn    = real(nn, dp)
    del    = 1.0_dp
    rq     = ivarmin * dn
    !
    !  Initial or restarted loop.
    !
    do
       p(1:n,nn) = ipstart(1:n)
       y(nn)     = func(ipstart, xx, yy)
       ineval     = ineval + 1
       !
       !  Define the initial simplex.
       !
       do j = 1, n
          x        = ipstart(j)
          ipstart(j) = ipstart(j) + istep(j) * del
          p(1:n,j) = ipstart(1:n)
          y(j)     = func(ipstart, xx, yy)
          ineval    = ineval + 1
          ipstart(j) = x
       end do
       !
       !  Find highest and lowest Y values.  FUNCMIN = Y(IHI) indicates
       !  the vertex of the simplex to be replaced.
       !
       ilo = minloc(y(1:n+1), 1)
       ylo = y(ilo)
       !
       !  Inner loop.
       !
       do while ( ineval < imaxeval )
          !
          !  FUNCMIN is, of course, the HIGHEST value???
          !
          ihi     = maxloc(y(1:n+1), 1)
          ifuncmin = y(ihi)
          !
          !  Calculate PBAR, the centroid of the simplex vertices
          !  excepting the vertex with Y value FUNCMIN.
          !
          do i = 1, n
             pbar(i) = ( sum(p(i,1:n+1)) - p(i,ihi) ) / dn
          end do
          !
          !  Reflection through the centroid.
          !
          pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
          ystar      = func(pstar, xx, yy)
          ineval      = ineval + 1
          !
          !  Successful reflection, so extension.
          !
          if ( ystar < ylo ) then

             p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
             y2star      = func(p2star, xx, yy)
             ineval       = ineval + 1
             !
             !  Retain extension or contraction.
             !
             if ( ystar < y2star ) then
                p(1:n,ihi) = pstar(1:n)
                y(ihi)     = ystar
             else
                p(1:n,ihi) = p2star(1:n)
                y(ihi)     = y2star
             end if
             !
             !  No extension.
             !
          else
             l = 0
             do i = 1, nn
                if ( ystar < y(i) ) l = l + 1
             end do

             if ( 1 < l ) then
                p(1:n,ihi) = pstar(1:n)
                y(ihi)     = ystar
                !
                !  Contraction on the Y(IHI) side of the centroid.
                !
             else if ( l == 0 ) then
                p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
                y2star      = func(p2star, xx, yy)
                ineval       = ineval + 1
                !
                !  Contract the whole simplex.
                !
                if ( y(ihi) < y2star ) then
                   do j = 1, n + 1
                      p(1:n,j)  = ( p(1:n,j) + p(1:n,ilo) ) * 0.5_dp
                      nelminxy(1:n) = p(1:n,j)
                      y(j)      = func(nelminxy, xx, yy)
                      ineval     = ineval + 1
                   end do
                   ilo = minloc(y(1:n+1), 1)
                   ylo = y(ilo)

                   cycle
                   !
                   !  Retain contraction.
                   !
                else
                   p(1:n,ihi) = p2star(1:n)
                   y(ihi)     = y2star
                end if
                !
                !  Contraction on the reflection side of the centroid.
                !
             else if ( l == 1 ) then

                p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
                y2star      = func(p2star, xx, yy)
                ineval       = ineval + 1
                !
                !  Retain reflection?
                !
                if ( y2star <= ystar ) then
                   p(1:n,ihi) = p2star(1:n)
                   y(ihi)     = y2star
                else
                   p(1:n,ihi) = pstar(1:n)
                   y(ihi)     = ystar
                end if

             end if

          end if
          !
          !  Check if YLO improved.
          !
          if ( y(ihi) < ylo ) then
             ylo = y(ihi)
             ilo = ihi
          end if

          jcount = jcount - 1

          if ( 0 < jcount ) cycle
          !
          !  Check to see if minimum reached.
          !
          if ( ineval <= imaxeval ) then
             jcount = ikonvge
             x = sum ( y(1:n+1) ) / dnn
             z = sum ( ( y(1:n+1) - x )**2 )
             if ( z <= rq ) exit
          end if

       end do
       !
       !  Factorial tests to check that FUNCMIN is a local minimum.
       !
       nelminxy(1:n) = p(1:n,ilo)
       ifuncmin   = y(ilo)

       if ( imaxeval < ineval ) then
          iierror = 2
          exit
       end if

       iierror = 0

       do i = 1, n
          del     = istep(i) * eps
          nelminxy(i) = nelminxy(i) + del
          z       = func(nelminxy, xx, yy)
          ineval  = ineval + 1
          if ( z < ifuncmin ) then
             iierror = 2
             exit
          end if
          nelminxy(i) = nelminxy(i) - del - del
          z       = func(nelminxy, xx, yy)
          ineval  = ineval + 1
          if ( z < ifuncmin ) then
             iierror = 2
             exit
          end if
          nelminxy(i) = nelminxy(i) + del
       end do

       if ( iierror == 0 ) exit
       !
       !  Restart the procedure.
       !
       ipstart(1:n) = nelminxy(1:n)
       del          = eps
       inumrestart  = inumrestart + 1

    end do

    if (present(funcmin)) then
       funcmin = ifuncmin
    endif
    if (present(neval)) then
       neval = ineval
    endif
    if (present(numrestart)) then
       numrestart = inumrestart
    endif
    if (present(ierror)) then
       ierror = iierror
    endif

  end function nelminxy

  ! ------------------------------------------------------------------

END MODULE mo_nelmin

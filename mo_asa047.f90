MODULE mo_asa047


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
  PUBLIC :: nelmin_opt      ! same as nelmin but with optional inputs to function

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         nelmin(_opt)

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
  !           FUNCTION func(x)
  !             USE mo_kind, ONLY: dp
  !             IMPLICIT NONE
  !             REAL(dp), DIMENSION(:),           INTENT(IN) :: x
  !             REAL(dp) :: func
  !           END FUNCTION func
  !
  !         and for nelmin_opt
  !           FUNCTION func(x,a)
  !             USE mo_kind, ONLY: dp
  !             IMPLICIT NONE
  !             REAL(dp), DIMENSION(:),           INTENT(IN) :: x
  !             REAL(dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: a
  !             REAL(dp) :: func
  !           END FUNCTION func
  !
  !         This routine does not include a termination test using the
  !         fitting of a quadratic surface.

  !     CALLING SEQUENCE
  !         call nelmin(func, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
  !                     icount, numres, ifault)
  !         call nelmin_opt(func, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
  !                         icount, numres, ifault, par=par)

  !     INDENT(IN)
  !         real(dp) :: reqmin       the terminating limit for the variance
  !                                  of the function values. reqmin>0 is required.
  !
  !         real(dp) :: step(:)      determines the size and shape of the initial simplex. 
  !                                  The relative magnitudes of its elements should reflect
  !                                  the units of the variables. size(step)=size(start)
  !
  !         integer(i4) :: konvge    the convergence check is carried out every konvge iterations.
  !                                  konvge>0 is required.
  !
  !         integer(i4) :: kcount    the maximum number of function evaluations.

  !     INDENT(INOUT)
  !         real(dp) :: start(:)     on input, a starting point for the iteration.
  !                                  on output, this data may have been overwritten.

  !     INDENT(OUT)
  !         real(dp) :: xmin(size(start))  the coordinates of the point which is estimated to minimize the function.
  !
  !         real(dp) :: ynewlo             the minimum value of the function.
  !
  !         integer(i4) :: icount          the number of function evaluations used.
  !
  !         integer(i4) :: numres          the number of restarts.
  !
  !         integer(i4) :: ifault          error indicator.
  !                                        0: no errors detected.
  !                                        1: reqmin or konvge have an illegal value.
  !                                        2: iteration terminated because kcount was exceeded without convergence.

  !     INDENT(IN), OPTIONAL
  !         real(dp) :: par(:)        Parameters to pass to function.

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None.

  !     EXAMPLE
  !         start(1:n) = (/ -1.2_dp, 1.0_dp /)
  !         reqmin = 1.0E-08_dp
  !         step(1:n) = (/ 1.0_dp, 1.0_dp /)
  !         konvge = 10
  !         kcount = 500
  !         ynewlo = rosenbrock(start)
  !         call nelmin(rosenbrock, start, xmin, ynewlo, reqmin, step, &
  !                     konvge, kcount, icount, numres, ifault)
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
  !         Modified,  Matthias Cuntz, Jul 2012 - optional parameters to function
  !                                             - i4, dp, intent

  subroutine nelmin(func, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
       icount, numres, ifault)

    implicit none

    INTERFACE
       FUNCTION func(xx)
         USE mo_kind, ONLY: dp
         IMPLICIT NONE
         REAL(dp), DIMENSION(:),           INTENT(IN) :: xx
         REAL(dp) :: func
       END FUNCTION func
    END INTERFACE
    real(dp),              intent(INOUT) :: start(:)
    real(dp),              intent(OUT)   :: xmin(size(start))
    real(dp),              intent(OUT)   :: ynewlo
    real(dp),              intent(IN)    :: reqmin
    real(dp),              intent(IN)    :: step(:)
    integer(i4),           intent(IN)    :: konvge
    integer(i4),           intent(IN)    :: kcount
    integer(i4),           intent(OUT)   :: icount
    integer(i4),           intent(OUT)   :: numres
    integer(i4),           intent(OUT)   :: ifault

    real(dp), parameter :: ccoeff = 0.5_dp
    real(dp), parameter :: ecoeff = 2.0_dp
    real(dp), parameter :: eps = 0.001_dp
    real(dp), parameter :: rcoeff = 1.0_dp
    integer(i4) :: n
    integer(i4) :: i
    integer(i4) :: ihi
    integer(i4) :: ilo
    integer(i4) :: j
    integer(i4) :: jcount
    integer(i4) :: l
    real(dp) :: del
    real(dp) :: p(size(start),size(start)+1)
    real(dp) :: p2star(size(start))
    real(dp) :: pbar(size(start))
    real(dp) :: pstar(size(start))
    real(dp) :: rq
    real(dp) :: x
    real(dp) :: y(size(start)+1)
    real(dp) :: y2star
    real(dp) :: ylo
    real(dp) :: ystar
    real(dp) :: z
    !
    !  Check the input parameters.
    !
    if ( reqmin <= 0.0_dp ) then
       ifault = 1
       return
    end if

    if ( konvge < 1 ) then
       ifault = 1
       return
    end if
    !
    !  Initialization.
    !
    n      = size(start)
    icount = 0_i4
    numres = 0_i4
    jcount = konvge
    del    = 1.0_dp
    rq     = reqmin * real(n,dp)
    !
    !  Initial or restarted loop.
    !
    do
       p(1:n,n+1) = start(1:n)
       y(n+1)     = func(start)
       icount     = icount + 1
       !
       !  Define the initial simplex.
       !
       do j = 1, n
          x        = start(j)
          start(j) = start(j) + step(j) * del
          p(1:n,j) = start(1:n)
          y(j)     = func(start)
          icount   = icount + 1
          start(j) = x
       end do
       !
       !  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
       !  the vertex of the simplex to be replaced.
       !
       ilo = minloc(y(1:n+1), 1)
       ylo = y(ilo)
       !
       !  Inner loop.
       !
       do while ( icount < kcount )
          !
          !  YNEWLO is, of course, the HIGHEST value???
          !
          ihi    = maxloc(y(1:n+1), 1)
          ynewlo = y(ihi)
          !
          !  Calculate PBAR, the centroid of the simplex vertices
          !  excepting the vertex with Y value YNEWLO.
          !
          do i = 1, n
             pbar(i) = ( sum(p(i,1:n+1)) - p(i,ihi) ) / real(n,dp)
          end do
          !
          !  Reflection through the centroid.
          !
          pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
          ystar      = func(pstar)
          icount     = icount + 1
          !
          !  Successful reflection, so extension.
          !
          if ( ystar < ylo ) then

             p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
             y2star      = func(p2star)
             icount      = icount + 1
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
             do i = 1, n + 1
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
                icount      = icount + 1
                !
                !  Contract the whole simplex.
                !
                if ( y(ihi) < y2star ) then
                   do j = 1, n + 1
                      p(1:n,j)  = ( p(1:n,j) + p(1:n,ilo) ) * 0.5_dp
                      xmin(1:n) = p(1:n,j)
                      y(j)      = func(xmin)
                      icount    = icount + 1
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
                icount      = icount + 1
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
          if ( icount <= kcount ) then
             jcount = konvge
             x = sum ( y(1:n+1) ) / real(n + 1,dp)
             z = sum ( ( y(1:n+1) - x )**2 )
             if ( z <= rq ) exit
          end if

       end do
       !
       !  Factorial tests to check that YNEWLO is a local minimum.
       !
       xmin(1:n) = p(1:n,ilo)
       ynewlo    = y(ilo)

       if ( kcount < icount ) then
          ifault = 2
          exit
       end if

       ifault = 0

       do i = 1, n
          del     = step(i) * eps
          xmin(i) = xmin(i) + del
          z       = func(xmin)
          icount  = icount + 1
          if ( z < ynewlo ) then
             ifault = 2
             exit
          end if
          xmin(i) = xmin(i) - del - del
          z       = func(xmin)
          icount  = icount + 1
          if ( z < ynewlo ) then
             ifault = 2
             exit
          end if
          xmin(i) = xmin(i) + del
       end do

       if ( ifault == 0 ) exit
       !
       !  Restart the procedure.
       !
       start(1:n) = xmin(1:n)
       del        = eps
       numres     = numres + 1

    end do

    return

  end subroutine nelmin


  subroutine nelmin_opt(func, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
       icount, numres, ifault, par)

    implicit none

    INTERFACE
       FUNCTION func(xx,aa)
         USE mo_kind, ONLY: dp
         IMPLICIT NONE
         REAL(dp), DIMENSION(:),           INTENT(IN) :: xx
         REAL(dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: aa
         REAL(dp) :: func
       END FUNCTION func
    END INTERFACE
    real(dp),              intent(INOUT) :: start(:)
    real(dp),              intent(OUT)   :: xmin(size(start))
    real(dp),              intent(OUT)   :: ynewlo
    real(dp),              intent(IN)    :: reqmin
    real(dp),              intent(IN)    :: step(:)
    integer(i4),           intent(IN)    :: konvge
    integer(i4),           intent(IN)    :: kcount
    integer(i4),           intent(OUT)   :: icount
    integer(i4),           intent(OUT)   :: numres
    integer(i4),           intent(OUT)   :: ifault
    real(dp),    optional, intent(IN)    :: par(:)

    real(dp), parameter :: ccoeff = 0.5_dp
    real(dp), parameter :: ecoeff = 2.0_dp
    real(dp), parameter :: eps = 0.001_dp
    real(dp), parameter :: rcoeff = 1.0_dp
    integer(i4) :: n
    integer(i4) :: i
    integer(i4) :: ihi
    integer(i4) :: ilo
    integer(i4) :: j
    integer(i4) :: jcount
    integer(i4) :: l
    real(dp) :: del
    real(dp) :: p(size(start),size(start)+1)
    real(dp) :: p2star(size(start))
    real(dp) :: pbar(size(start))
    real(dp) :: pstar(size(start))
    real(dp) :: rq
    real(dp) :: x
    real(dp) :: y(size(start)+1)
    real(dp) :: y2star
    real(dp) :: ylo
    real(dp) :: ystar
    real(dp) :: z
    !
    !  Check the input parameters.
    !
    if ( reqmin <= 0.0_dp ) then
       ifault = 1
       return
    end if

    if ( konvge < 1 ) then
       ifault = 1
       return
    end if
    !
    !  Initialization.
    !
    n      = size(start)
    icount = 0_i4
    numres = 0_i4
    jcount = konvge
    del    = 1.0_dp
    rq     = reqmin * real(n,dp)
    !
    !  Initial or restarted loop.
    !
    do
       p(1:n,n+1) = start(1:n)
       if (present(par)) then
          y(n+1)     = func(start, par)
       else
          y(n+1)     = func(start)
       endif
       icount     = icount + 1
       !
       !  Define the initial simplex.
       !
       do j = 1, n
          x        = start(j)
          start(j) = start(j) + step(j) * del
          p(1:n,j) = start(1:n)
          if (present(par)) then
             y(j)     = func(start, par)
          else
             y(j)     = func(start)
          endif
          icount   = icount + 1
          start(j) = x
       end do
       !
       !  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
       !  the vertex of the simplex to be replaced.
       !
       ilo = minloc(y(1:n+1), 1)
       ylo = y(ilo)
       !
       !  Inner loop.
       !
       do while ( icount < kcount )
          !
          !  YNEWLO is, of course, the HIGHEST value???
          !
          ihi    = maxloc(y(1:n+1), 1)
          ynewlo = y(ihi)
          !
          !  Calculate PBAR, the centroid of the simplex vertices
          !  excepting the vertex with Y value YNEWLO.
          !
          do i = 1, n
             pbar(i) = ( sum(p(i,1:n+1)) - p(i,ihi) ) / real(n,dp)
          end do
          !
          !  Reflection through the centroid.
          !
          pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
          if (present(par)) then
             ystar      = func(pstar, par)
          else
             ystar      = func(pstar)
          endif
          icount     = icount + 1
          !
          !  Successful reflection, so extension.
          !
          if ( ystar < ylo ) then

             p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
             if (present(par)) then
                y2star      = func(p2star, par)
             else
                y2star      = func(p2star)
             endif
             icount      = icount + 1
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
             do i = 1, n + 1
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
                if (present(par)) then
                   y2star      = func(p2star, par)
                else
                   y2star      = func(p2star)
                endif
                icount      = icount + 1
                !
                !  Contract the whole simplex.
                !
                if ( y(ihi) < y2star ) then
                   do j = 1, n + 1
                      p(1:n,j)  = ( p(1:n,j) + p(1:n,ilo) ) * 0.5_dp
                      xmin(1:n) = p(1:n,j)
                      if (present(par)) then
                         y(j)      = func(xmin, par)
                      else
                         y(j)      = func(xmin)
                      endif
                      icount    = icount + 1
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
                if (present(par)) then
                   y2star      = func(p2star, par)
                else
                   y2star      = func(p2star)
                endif
                icount      = icount + 1
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
          if ( icount <= kcount ) then
             jcount = konvge
             x = sum ( y(1:n+1) ) / real(n + 1,dp)
             z = sum ( ( y(1:n+1) - x )**2 )
             if ( z <= rq ) exit
          end if

       end do
       !
       !  Factorial tests to check that YNEWLO is a local minimum.
       !
       xmin(1:n) = p(1:n,ilo)
       ynewlo    = y(ilo)

       if ( kcount < icount ) then
          ifault = 2
          exit
       end if

       ifault = 0

       do i = 1, n
          del     = step(i) * eps
          xmin(i) = xmin(i) + del
          if (present(par)) then
             z       = func(xmin, par)
          else
             z       = func(xmin)
          endif
          icount  = icount + 1
          if ( z < ynewlo ) then
             ifault = 2
             exit
          end if
          xmin(i) = xmin(i) - del - del
          if (present(par)) then
             z       = func(xmin, par)
          else
             z       = func(xmin)
          endif
          icount  = icount + 1
          if ( z < ynewlo ) then
             ifault = 2
             exit
          end if
          xmin(i) = xmin(i) + del
       end do

       if ( ifault == 0 ) exit
       !
       !  Restart the procedure.
       !
       start(1:n) = xmin(1:n)
       del        = eps
       numres     = numres + 1

    end do

    return

  end subroutine nelmin_opt

  ! ------------------------------------------------------------------

END MODULE mo_asa047

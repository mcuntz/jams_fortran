MODULE mo_nelmin


  ! This module provides NELMIN, which minimizes a function using the Nelder-Mead algorithm
  ! with the Applied Statistics algorithms No. 047.

  ! Written  Matthias Cuntz, Jul 2012 - extended asa047 of John Burkardt

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! It is NOT released under the GNU Lesser General Public License, yet.

  ! If you use this routine, please contact Matthias Cuntz.

  ! Copyright 2012 Matthias Cuntz

  USE mo_kind,  ONLY: i4, sp, dp
  USE mo_utils, ONLY: ne

  IMPLICIT NONE

  PUBLIC :: nelmin          ! minimizes a function using the Nelder-Mead algorithm
  PUBLIC :: nelminxy        ! same as nelmin but pass x and y to function
  PUBLIC :: nelminrange     ! same as nelmin but with given range of parameter

  ! ------------------------------------------------------------------

  !     NAME
  !         nelmin
  !         nelminxy
  !         nelminrange

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
  !         pmin = nelmin(func, pstart, prange, funcmin, varmin, step, konvge, maxeval, &
  !                       neval, numrestart, ierror)

  !     INTENT(IN)
  !         real(dp) :: func(p,xx,yy)    Function on which to search the minimum
  !         real(dp) :: pstart(:)        Starting point for the iteration.
  !         real(dp) :: xx               First values to pass as function arguments
  !         real(dp) :: yy               Second values to pass as function arguments
  !         real(dp) :: prange(:,2)      Range of parameters (upper and lower bound).

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(dp) :: nelmin(size(start))    the coordinates of the point which is estimated to minimize the function.

  !     INTENT(IN), OPTIONAL
  !         real(dp) :: varmin           the terminating limit for the variance
  !                                      of the function values. varmin>0 is required.
  !         real(dp) :: step(:)          determines the size and shape of the initial simplex. 
  !                                      The relative magnitudes of its elements should reflect
  !                                      the units of the variables. size(step)=size(start)
  !         integer(i4) :: konvge        the convergence check is carried out every konvge iterations.
  !                                      konvge>0 is required.
  !         integer(i4) :: maxeval       the maximum number of function evaluations.
  !                                      default: 1000

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         real(dp)    :: funcmin             the minimum value of the function.
  !         integer(i4) :: neval               the number of function evaluations used.
  !         integer(i4) :: numrestart          the number of restarts.
  !         integer(i4) :: ierror              error indicator.
  !                                            0: no errors detected.
  !                                            1: varmin or konvge have an illegal value.
  !                                            2: iteration terminated because maxeval was exceeded without convergence.
  !         real(dp), allocatable 
  !                     :: history(:)          the history of best function values, history(neval)=funcmin

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
  !                                             - nelminxy
  !                    Juliane Mai,    Aug 2012 - nelminrange
  !                    Juliane Mai,    Dec 2012 - history output

  INTERFACE nelminrange
     MODULE PROCEDURE nelminrange_dp, nelminrange_sp
  END INTERFACE nelminrange

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  function nelmin(func, pstart, varmin, step, konvge, maxeval, &
       funcmin, neval, numrestart, ierror, history)

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
    real(dp),    optional, intent(IN)  :: varmin
    real(dp),    optional, intent(IN)  :: step(:)
    integer(i4), optional, intent(IN)  :: konvge
    integer(i4), optional, intent(IN)  :: maxeval
    real(dp),    optional, intent(OUT) :: funcmin
    integer(i4), optional, intent(OUT) :: neval
    integer(i4), optional, intent(OUT) :: numrestart
    integer(i4), optional, intent(OUT) :: ierror
    real(dp),    optional, intent(OUT), allocatable :: history(:)  ! History of objective function values
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
    real(dp), allocatable :: history_tmp(:)
    !
    ! Defaults
    !
    nelmin(:) = 0.
    if (present(varmin)) then
       if (varmin <= 0.0_dp) stop 'Error nelmin: varmin<0'
       ivarmin = varmin
    else
       ivarmin = 1.0e-9_dp
    endif
    ! maximal number of function evaluations
    if (present(maxeval)) then
       if (maxeval <= 1) stop 'Error nelmin: maxeval<=1'
       imaxeval = maxeval
    else
       imaxeval = 1000
    endif
    ! history output
    if (present(history)) then
       ! worst case length
       allocate(history_tmp(imaxeval+3*size(ipstart)+1))
    end if
    if (present(konvge)) then
       ikonvge = konvge
    else
       ikonvge = imaxeval / 10
    endif
    if (ikonvge < 1) stop 'Error nelmin: konvg<1'
    !
    if (present(step)) then
       istep = step
    else ! if not given, deviate initial by 1%
       y0 = func(pstart)
       do i=1, size(pstart)
          p0 = pstart
          if (ne(p0(i),0.0_dp)) then
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
       if (present(history)) history_tmp(ineval) = y(nn)
       !
       !  Define the initial simplex.
       !
       do j = 1, n
          x          = ipstart(j)
          ipstart(j) = ipstart(j) + istep(j) * del
          p(1:n,j)   = ipstart(1:n)
          y(j)       = func(ipstart)
          ineval     = ineval + 1
          ipstart(j) = x
          if (present(history)) history_tmp(ineval) = Min( y(j),  history_tmp(ineval-1) )
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
          if (present(history)) history_tmp(ineval) = Min( ystar, history_tmp(ineval-1) )
          !
          !  Successful reflection, so extension.
          !
          if ( ystar < ylo ) then

             p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
             y2star      = func(p2star)
             ineval       = ineval + 1
             if (present(history)) history_tmp(ineval) = Min( y2star, history_tmp(ineval-1) )
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
                if (present(history)) history_tmp(ineval) = Min( y2star, history_tmp(ineval-1) )
                !
                !  Contract the whole simplex.
                !
                if ( y(ihi) < y2star ) then
                   do j = 1, n + 1
                      p(1:n,j)  = ( p(1:n,j) + p(1:n,ilo) ) * 0.5_dp
                      nelmin(1:n) = p(1:n,j)
                      y(j)      = func(nelmin)
                      ineval     = ineval + 1
                      if (present(history)) history_tmp(ineval) = Min( y(j), history_tmp(ineval-1) )
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
                if (present(history)) history_tmp(ineval) = Min( y2star, history_tmp(ineval-1) )
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
          if (present(history)) history_tmp(ineval) = Min( z, history_tmp(ineval-1) )
          if ( z < ifuncmin ) then
             iierror = 2
             exit
          end if
          nelmin(i) = nelmin(i) - del - del
          z       = func(nelmin)
          ineval  = ineval + 1
          if (present(history)) history_tmp(ineval) = Min( z, history_tmp(ineval-1) )
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
    if (present(history)) then
       allocate(history(ineval))
       history(:) = history_tmp(1:ineval)
    end if

  end function nelmin


  function nelminxy(func, pstart, xx, yy, varmin, step, konvge, maxeval, &
       funcmin, neval, numrestart, ierror, history)

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
    real(dp),    optional, intent(OUT), allocatable :: history(:)  ! History of objective function values
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
    real(dp), allocatable :: history_tmp(:)
    !
    ! Defaults
    !
    nelminxy(:) = 0.
    if (present(varmin)) then
       if (varmin <= 0.0_dp) stop 'Error nelminxy: varmin<0'
       ivarmin = varmin
    else
       ivarmin = 1.0e-9_dp
    endif
    ! maximal number of function evaluations
    if (present(maxeval)) then
       if (maxeval <= 1) stop 'Error nelminxy: maxeval<=1'
       imaxeval = maxeval
    else
       imaxeval = 1000
    endif
    ! history output
    if (present(history)) then
       ! worst case length
       allocate(history_tmp(imaxeval+3*size(ipstart)+1))
    end if
    if (present(konvge)) then
       ikonvge = konvge
    else
       ikonvge = imaxeval / 10
    endif
    if (ikonvge < 1) stop 'Error nelminxy: konvg<1'
    if (size(xx) /= size(yy)) stop 'Error nelminxy: size(xx) /= size(yy)'
    !
    if (present(step)) then
       istep = step
    else ! if not given, deviate initial by 1%
       y0 = func(pstart, xx, yy)
       do i=1, size(pstart)
          p0 = pstart
          if (ne(p0(i),0.0_dp)) then
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
       if (present(history)) history_tmp(ineval) = y(nn)
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
          if (present(history)) history_tmp(ineval) = Min( y(j),  history_tmp(ineval-1) )
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
          if (present(history)) history_tmp(ineval) = Min( ystar,  history_tmp(ineval-1) )
          !
          !  Successful reflection, so extension.
          !
          if ( ystar < ylo ) then

             p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
             y2star      = func(p2star, xx, yy)
             ineval       = ineval + 1
             if (present(history)) history_tmp(ineval) = Min( y2star,  history_tmp(ineval-1) )
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
                if (present(history)) history_tmp(ineval) = Min( y2star,  history_tmp(ineval-1) )
                !
                !  Contract the whole simplex.
                !
                if ( y(ihi) < y2star ) then
                   do j = 1, n + 1
                      p(1:n,j)  = ( p(1:n,j) + p(1:n,ilo) ) * 0.5_dp
                      nelminxy(1:n) = p(1:n,j)
                      y(j)      = func(nelminxy, xx, yy)
                      ineval     = ineval + 1
                      if (present(history)) history_tmp(ineval) = Min( y(j),  history_tmp(ineval-1) )
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
                if (present(history)) history_tmp(ineval) = Min( y2star,  history_tmp(ineval-1) )
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
          if (present(history)) history_tmp(ineval) = Min( z,  history_tmp(ineval-1) )
          if ( z < ifuncmin ) then
             iierror = 2
             exit
          end if
          nelminxy(i) = nelminxy(i) - del - del
          z       = func(nelminxy, xx, yy)
          ineval  = ineval + 1
          if (present(history)) history_tmp(ineval) = Min( z,  history_tmp(ineval-1) )
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
    if (present(history)) then
       allocate(history(ineval))
       history(:) = history_tmp(1:ineval)
    end if

  end function nelminxy

  function nelminrange_dp(func, pstart, prange, varmin, step, konvge, maxeval, &
       funcmin, neval, numrestart, ierror, history)

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
    real(dp),              intent(IN)  :: prange(:,:)
    real(dp),    optional, intent(IN)  :: varmin
    real(dp),    optional, intent(IN)  :: step(:)
    integer(i4), optional, intent(IN)  :: konvge
    integer(i4), optional, intent(IN)  :: maxeval
    real(dp),    optional, intent(OUT) :: funcmin
    integer(i4), optional, intent(OUT) :: neval
    integer(i4), optional, intent(OUT) :: numrestart
    integer(i4), optional, intent(OUT) :: ierror
    real(dp),    optional, intent(OUT), allocatable :: history(:)  ! History of objective function values
    real(dp)                           :: nelminrange_dp(size(pstart))

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
    real(dp), allocatable :: history_tmp(:)
    !
    ! Defaults
    !
    nelminrange_dp(:) = 0.5_dp * ( prange(:,1) + prange(:,2) )

    if (present(varmin)) then
       if (varmin <= 0.0_dp) stop 'Error nelmin: varmin<0'
       ivarmin = varmin
    else
       ivarmin = 1.0e-9_dp
    endif
    ! maximal number of function evaluations
    if (present(maxeval)) then
       if (maxeval <= 1) stop 'Error nelmin: maxeval<=1'
       imaxeval = maxeval
    else
       imaxeval = 1000
    endif
    ! history output
    if (present(history)) then
       ! worst case length
       allocate(history_tmp(imaxeval+3*size(ipstart)+1))
    end if
    if (present(konvge)) then
       ikonvge = konvge
    else
       ikonvge = imaxeval / 10
    endif
    if (ikonvge < 1) stop 'Error nelmin: konvg<1'
    !
    if (present(step)) then
       istep = step
    else ! if not given, deviate initial by 1%
       y0 = func(pstart)
       do i=1, size(pstart)
          p0 = pstart
          if (ne(p0(i),0.0_dp)) then
             ! bound to range
             p0(i) = min( prange(i,2) , max( prange(i,1) , 1.01_dp*p0(i) ) )
          else
             ! bound to range
             p0(i) = min( prange(i,2) , max( prange(i,1) , 0.01_dp ) )
          endif
          istep(i) = abs(func(p0) - y0)
       enddo
    endif
    !
    !  Check the input parameters.
    !
    if (size(prange,2) .ne. 2_i4)           stop 'nelminrange_dp: range has to be array of size (:,2)'
    if (size(prange,1) .ne.  size(ipstart)) stop 'nelminrange_dp: range has to be given for each dimension'
    if (.not. (all(prange(:,1) .le. pstart(:)) .and. all(prange(:,2) .ge. pstart(:)) )) then
       stop 'nelminrange_dp: starting point is not in range'
    end if
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
       if (present(history)) history_tmp(ineval) = y(nn)
       !
       !  Define the initial simplex.
       !
       do j = 1, n
          x          = ipstart(j)
          ! bound to range
          ipstart(j) =  min( prange(j,2) , max( prange(j,1) , ipstart(j) + istep(j) * del ) )
          p(1:n,j)   = ipstart(1:n)
          y(j)       = func(ipstart)
          ineval     = ineval + 1
          ipstart(j) = x
          if (present(history)) history_tmp(ineval) = Min( y(j),  history_tmp(ineval-1) )
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
          ! bound to range
          pstar(1:n) = min( prange(1:n,2) , max( prange(1:n,1) , pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) ) ) )  
          ystar      = func(pstar)
          ineval      = ineval + 1
          if (present(history)) history_tmp(ineval) = Min( ystar,  history_tmp(ineval-1) )
          !
          !  Successful reflection, so extension.
          !
          if ( ystar < ylo ) then

             ! bound to range
             p2star(1:n) = min( prange(1:n,2) , max( prange(1:n,1) , pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) ) ) )
             y2star      = func(p2star)
             ineval       = ineval + 1
             if (present(history)) history_tmp(ineval) = Min( y2star,  history_tmp(ineval-1) )
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
                ! bound to range
                p2star(1:n)  = min( prange(1:n,2) , max( prange(1:n,1) , pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) ) ) )  
                y2star       = func(p2star)
                ineval       = ineval + 1
                if (present(history)) history_tmp(ineval) = Min( y2star,  history_tmp(ineval-1) )
                !
                !  Contract the whole simplex.
                !
                if ( y(ihi) < y2star ) then
                   do j = 1, n + 1
                      ! bound to range
                      p(1:n,j)  = min( prange(1:n,2) , max( prange(1:n,1) , ( p(1:n,j) + p(1:n,ilo) ) * 0.5_dp ) ) 
                      nelminrange_dp(1:n) = p(1:n,j)
                      y(j)      = func(nelminrange_dp)
                      ineval     = ineval + 1
                      if (present(history)) history_tmp(ineval) = Min( y(j),  history_tmp(ineval-1) )
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

                ! bound to range
                p2star(1:n) =  min( prange(1:n,2) , max( prange(1:n,1) , pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) ) ) ) 
                y2star      = func(p2star)
                ineval      = ineval + 1
                if (present(history)) history_tmp(ineval) = Min( y2star,  history_tmp(ineval-1) )
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
       ! bound to range
       nelminrange_dp(1:n) = min( prange(1:n,2) , max( prange(1:n,1) , p(1:n,ilo) ) ) 
       ifuncmin   = y(ilo)

       if ( imaxeval < ineval ) then
          iierror = 2
          exit
       end if

       iierror = 0

       do i = 1, n
          del     = istep(i) * eps
          ! bound to range: probably not necessary
          nelminrange_dp(i) = min( prange(i,2) , max( prange(i,1) , nelminrange_dp(i) + del ) ) 
          z       = func(nelminrange_dp)
          ineval  = ineval + 1
          if (present(history)) history_tmp(ineval) = Min( z,  history_tmp(ineval-1) )
          if ( z < ifuncmin ) then
             iierror = 2
             exit
          end if
          ! bound to range: probably not necessary
          nelminrange_dp(i) = min( prange(i,2) , max( prange(i,1) , nelminrange_dp(i) - del - del ) ) 
          z       = func(nelminrange_dp)
          ineval  = ineval + 1
          if (present(history)) history_tmp(ineval) = Min( z,  history_tmp(ineval-1) )
          if ( z < ifuncmin ) then
             iierror = 2
             exit
          end if
          ! bound to range: probably not necessary
          nelminrange_dp(i) = min( prange(i,2) , max( prange(i,1) , nelminrange_dp(i) + del ) ) 
       end do

       if ( iierror == 0 ) exit
       !
       !  Restart the procedure.
       !
       ! bound to range: probably not necessary
       ipstart(1:n) = min( prange(1:n,2) , max( prange(1:n,1) , nelminrange_dp(1:n) ) ) 
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
    if (present(history)) then
       allocate(history(ineval))
       history(:) = history_tmp(1:ineval)
    end if

  end function nelminrange_dp

  function nelminrange_sp(func, pstart, prange, varmin, step, konvge, maxeval, &
       funcmin, neval, numrestart, ierror, history)

    implicit none

    INTERFACE
       FUNCTION func(pp)
         USE mo_kind, ONLY: sp
         IMPLICIT NONE
         REAL(sp), DIMENSION(:), INTENT(IN) :: pp
         REAL(sp) :: func
       END FUNCTION func
    END INTERFACE
    real(sp),              intent(IN)  :: pstart(:)
    real(sp),              intent(IN)  :: prange(:,:)
    real(sp),    optional, intent(IN)  :: varmin
    real(sp),    optional, intent(IN)  :: step(:)
    integer(i4), optional, intent(IN)  :: konvge
    integer(i4), optional, intent(IN)  :: maxeval
    real(sp),    optional, intent(OUT) :: funcmin
    integer(i4), optional, intent(OUT) :: neval
    integer(i4), optional, intent(OUT) :: numrestart
    integer(i4), optional, intent(OUT) :: ierror
    real(sp),    optional, intent(OUT), allocatable :: history(:)  ! History of objective function values
    real(sp)                           :: nelminrange_sp(size(pstart))

    real(sp), parameter :: ccoeff = 0.5_sp
    real(sp), parameter :: ecoeff = 2.0_sp
    real(sp), parameter :: rcoeff = 1.0_sp
    real(sp), parameter :: eps    = 0.001_sp
    real(sp) :: ipstart(size(pstart))
    real(sp) :: ifuncmin
    real(sp) :: ivarmin
    real(sp) :: istep(size(pstart))
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
    real(sp) :: del
    real(sp) :: p(size(pstart),size(pstart)+1)
    real(sp) :: p2star(size(pstart))
    real(sp) :: pbar(size(pstart))
    real(sp) :: pstar(size(pstart))
    real(sp) :: rq
    real(sp) :: x
    real(sp) :: y(size(pstart)+1)
    real(sp) :: y2star
    real(sp) :: ylo
    real(sp) :: ystar
    real(sp) :: z
    real(sp) :: dn, dnn
    real(sp) :: p0(size(pstart)), y0
    real(sp), allocatable :: history_tmp(:)
    !
    ! Defaults
    !
    nelminrange_sp(:) = 0.5_sp * ( prange(:,1) + prange(:,2) )

    if (present(varmin)) then
       if (varmin <= 0.0_sp) stop 'Error nelmin: varmin<0'
       ivarmin = varmin
    else
       ivarmin = 1.0e-9_sp
    endif
    ! maximal number of function evaluations
    if (present(maxeval)) then
       if (maxeval <= 1) stop 'Error nelmin: maxeval<=1'
       imaxeval = maxeval
    else
       imaxeval = 1000
    endif
    ! history output
    if (present(history)) then
       ! worst case length
       allocate(history_tmp(imaxeval+3*size(ipstart)+1))
    end if
    if (present(konvge)) then
       ikonvge = konvge
    else
       ikonvge = imaxeval / 10
    endif
    if (ikonvge < 1) stop 'Error nelmin: konvg<1'
    !
    if (present(step)) then
       istep = step
    else ! if not given, deviate initial by 1%
       y0 = func(pstart)
       do i=1, size(pstart)
          p0 = pstart
          if (ne(p0(i),0.0_sp)) then
             ! bound to range
             p0(i) = min( prange(i,2) , max( prange(i,1) , 1.01_sp*p0(i) ) )
          else
             ! bound to range
             p0(i) = min( prange(i,2) , max( prange(i,1) , 0.01_sp ) )
          endif
          istep(i) = abs(func(p0) - y0)
       enddo
    endif
    !
    !  Check the input parameters.
    !
    if (size(prange,2) .ne. 2_i4)           stop 'nelminrange_sp: range has to be array of size (:,2)'
    if (size(prange,1) .ne.  size(ipstart)) stop 'nelminrange_sp: range has to be given for each dimension'
    if (.not. (all(prange(:,1) .le. pstart(:)) .and. all(prange(:,2) .ge. pstart(:)) )) then
       stop 'nelminrange_sp: starting point is not in range'
    end if
    !
    !  Initialization.
    !
    ipstart = pstart
    n       = size(ipstart)
    ineval  = 0_i4
    inumrestart = 0_i4
    jcount = ikonvge
    dn     = real(n, sp)
    nn     = n + 1
    dnn    = real(nn, sp)
    del    = 1.0_sp
    rq     = ivarmin * dn
    !
    !  Initial or restarted loop.
    !
    do
       p(1:n,nn) = ipstart(1:n)
       y(nn)     = func(ipstart)
       ineval     = ineval + 1
       if (present(history)) history_tmp(ineval) = y(nn)
       !
       !  Define the initial simplex.
       !
       do j = 1, n
          x          = ipstart(j)
          ! bound to range
          ipstart(j) =  min( prange(j,2) , max( prange(j,1) , ipstart(j) + istep(j) * del ) )
          p(1:n,j)   = ipstart(1:n)
          y(j)       = func(ipstart)
          ineval     = ineval + 1
          ipstart(j) = x
          if (present(history)) history_tmp(ineval) = Min( y(j),  history_tmp(ineval-1) )
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
          ! bound to range
          pstar(1:n) = min( prange(1:n,2) , max( prange(1:n,1) , pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) ) ) )  
          ystar      = func(pstar)
          ineval      = ineval + 1
          if (present(history)) history_tmp(ineval) = Min( ystar,  history_tmp(ineval-1) )
          !
          !  Successful reflection, so extension.
          !
          if ( ystar < ylo ) then

             ! bound to range
             p2star(1:n) = min( prange(1:n,2) , max( prange(1:n,1) , pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) ) ) )
             y2star      = func(p2star)
             ineval       = ineval + 1
             if (present(history)) history_tmp(ineval) = Min( y2star,  history_tmp(ineval-1) )
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
                ! bound to range
                p2star(1:n)  = min( prange(1:n,2) , max( prange(1:n,1) , pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) ) ) )  
                y2star       = func(p2star)
                ineval       = ineval + 1
                if (present(history)) history_tmp(ineval) = Min( y2star,  history_tmp(ineval-1) )
                !
                !  Contract the whole simplex.
                !
                if ( y(ihi) < y2star ) then
                   do j = 1, n + 1
                      ! bound to range
                      p(1:n,j)  = min( prange(1:n,2) , max( prange(1:n,1) , ( p(1:n,j) + p(1:n,ilo) ) * 0.5_sp ) ) 
                      nelminrange_sp(1:n) = p(1:n,j)
                      y(j)      = func(nelminrange_sp)
                      ineval     = ineval + 1
                      if (present(history)) history_tmp(ineval) = Min( y(j),  history_tmp(ineval-1) )
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

                ! bound to range
                p2star(1:n) =  min( prange(1:n,2) , max( prange(1:n,1) , pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) ) ) ) 
                y2star      = func(p2star)
                ineval      = ineval + 1
                if (present(history)) history_tmp(ineval) = Min( y2star,  history_tmp(ineval-1) )
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
       ! bound to range
       nelminrange_sp(1:n) = min( prange(1:n,2) , max( prange(1:n,1) , p(1:n,ilo) ) ) 
       ifuncmin   = y(ilo)

       if ( imaxeval < ineval ) then
          iierror = 2
          exit
       end if

       iierror = 0

       do i = 1, n
          del     = istep(i) * eps
          ! bound to range: probably not necessary
          nelminrange_sp(i) = min( prange(i,2) , max( prange(i,1) , nelminrange_sp(i) + del ) ) 
          z       = func(nelminrange_sp)
          ineval  = ineval + 1
          if (present(history)) history_tmp(ineval) = Min( z,  history_tmp(ineval-1) )
          if ( z < ifuncmin ) then
             iierror = 2
             exit
          end if
          ! bound to range: probably not necessary
          nelminrange_sp(i) = min( prange(i,2) , max( prange(i,1) , nelminrange_sp(i) - del - del ) ) 
          z       = func(nelminrange_sp)
          ineval  = ineval + 1
          if (present(history)) history_tmp(ineval) = Min( z,  history_tmp(ineval-1) )
          if ( z < ifuncmin ) then
             iierror = 2
             exit
          end if
          ! bound to range: probably not necessary
          nelminrange_sp(i) = min( prange(i,2) , max( prange(i,1) , nelminrange_sp(i) + del ) ) 
       end do

       if ( iierror == 0 ) exit
       !
       !  Restart the procedure.
       !
       ! bound to range: probably not necessary
       ipstart(1:n) = min( prange(1:n,2) , max( prange(1:n,1) , nelminrange_sp(1:n) ) ) 
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
    if (present(history)) then
       allocate(history(ineval))
       history(:) = history_tmp(1:ineval)
    end if

  end function nelminrange_sp


  ! ------------------------------------------------------------------

END MODULE mo_nelmin

module mo_dds

  ! This module contains routines for Dynamically Dimensioned Search (DDS)

  ! Written Jul 2012, Matthias Cuntz - module version of modified DDS v1.1 of R. Kumar to
  ! original DDS of Bryan Tolson


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

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: DDS    ! Dynamically Dimensioned Search (DDS)
  PUBLIC :: DDSxy  ! same as DDS but pass x and y to function


  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         DDS
  !         DDSxy

  !     PURPOSE
  !         Searches Minimum or Maximum of a user-specified function using
  !         Dynamically Dimensioned Search (DDS)
  !
  !         DDS is an n-dimensional continuous global optimization algorithm.
  !         It is coded as a minimizer but one can give maxit=1.0 in a maximization problem,
  !         so that the algorithm minimizes the negative of the objective function F (-1*F).
  !
  !         The function to be minimized is the first argument of DDS and
  !         must be defined as
  !           FUNCTION func(p)
  !             USE mo_kind, ONLY: dp
  !             IMPLICIT NONE
  !             REAL(dp), DIMENSION(:),           INTENT(IN) :: p
  !             REAL(dp) :: func
  !           END FUNCTION func
  !
  !         and for DDSxy
  !           FUNCTION func(p,x,y)
  !             USE mo_kind, ONLY: dp
  !             IMPLICIT NONE
  !             REAL(dp), DIMENSION(:), INTENT(IN) :: p
  !             REAL(dp), DIMENSION(:), INTENT(IN) :: x
  !             REAL(dp), DIMENSION(:), INTENT(IN) :: y
  !             REAL(dp) :: func
  !           END FUNCTION func

  !     CALLING SEQUENCE
  !         popt = DDS(obj_func, pini, pmin, pmax, r, seed, maxiter, maxit, nobj)
  !         popt = DDSxy(obj_func, pini, pmin, pmax, r, seed, maxiter, maxit, nobj)

  !     INDENT(IN)
  !         real(dp) :: obj_func(p,xx,yy)    Function on which to search the minimum
  !         real(dp) :: pini(:)              inital value of decision variables 
  !         real(dp) :: pmin(:)              Min value of decision variables 
  !         real(dp) :: pmax(:)              Max value of decision variables 
  !         real(dp) :: pstart(:)            Starting point for the iteration.
  !         real(dp) :: xx                   First values to pass as function arguments
  !         real(dp) :: yy                   Second values to pass as function arguments

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(dp) :: DDS(size(pini))      The parameters of the point which is estimated to minimize the function.

  !     INDENT(IN), OPTIONAL
  !         real(dp)    :: r                  DDS perturbation parameter (default: 0.2)
  !         integer(i4) :: seed               User seed to initialise the random number generator (default: None)
  !         integer(i8) :: maxiter            Maximum number of iteration or function evaluation (default: 1000)
  !         real(dp)    :: maxit              Maximization or minimization of function
  !                                           1.0 if minimization, -1.0 if maximization (default: 1.0)

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         integer(i8) :: nobj               Number of obj. function evaluations

  !     RESTRICTIONS
  !         None.

  !     EXAMPLE
  !         dv_min = (/ -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0 /)
  !         dv_max = (/ 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0 /)
  !         dv_ini = (/ -.226265E+01, -.130187E+01, -.151219E+01, 0.133983E+00, 0.988159E+00, &
  !                     -.495074E+01, -.126574E+02, 0.572684E+00, 0.303864E+01, 0.343031E+01 /)
  !         dv_opt = DDS(griewank, dv_ini, dv_min, dv_max, r_val, seed, nIterMax, to_max, nobj)
  !         -> see also example in test directory

  !     LITERATURE
  !         Tolson, B. A., and C. A. Shoemaker (2007)
  !             Dynamically dimensioned search algorithm for computationally efficient watershed
  !             model calibration, Water Resour. Res., 43, W01413, doi:10.1029/2005WR004723.

  !     HISTORY
  !         Written original Bryan Tolson, Feb 2007 - DDS v1.1
  !         Modified, Rohini Kumar, Feb 2008
  !         Modified, Matthias Cuntz, Jul 2012 - module

  function DDS(obj_func, pini, pmin, pmax, r, seed, maxiter, maxit, nobj)

    use mo_kind, only: i4, i8, dp
#ifdef IMSL
    use RNSET_INT
    use RNUN_INT
#endif

    implicit none

    INTERFACE
       FUNCTION obj_func(pp)
         USE mo_kind, ONLY: dp
         IMPLICIT NONE
         REAL(dp), DIMENSION(:), INTENT(IN) :: pp
         REAL(dp) :: obj_func
       END FUNCTION obj_func
    END INTERFACE
    real(dp),    dimension(:), intent(in)  :: pini    ! inital value of decision variables 
    real(dp),    dimension(:), intent(in)  :: pmin    ! Min value of decision variables 
    real(dp),    dimension(:), intent(in)  :: pmax    ! Max value of decision variables 
    real(dp),        optional, intent(in)  :: r       ! DDS perturbation parameter (-> 0.2 by default)
    integer(i4),     optional, intent(in)  :: seed    ! User seed to initialise the random number generator
    integer(i8),     optional, intent(in)  :: maxiter ! Maximum number of iteration or function evaluation
    real(dp),        optional, intent(in)  :: maxit   ! Maximization or minimization of function
    integer(i8),     optional, intent(out) :: nobj ! Counter for Obj. function evaluations
    real(dp), dimension(size(pini)) :: DDS ! Best value of decision variables

    ! Local variables
    integer(i4)    :: pnum     ! Total number of decision variables  
    integer(i4)    :: iseed    ! User given seed
    integer(i8)    :: imaxiter ! Maximum number of iteration or function evaluation
    integer(i8)    :: iobj     ! Counter for Obj. function evaluations
    real(dp)       :: ir       ! DDS perturbation parameter
    real(dp)       :: imaxit   ! Maximization or minimization of function
    real(dp)       :: of_value, of_new, of_best, Pn, new_value
    real(dp), dimension(size(pini)) :: pnew   ! Test value of decision variables (dummy) 
    real(dp), dimension(1)          :: ranval
    integer(i8)    :: i, k, dvn_count, dv       
    integer(i4)    :: j
#ifndef IMSL
    real(dp)       :: harvest  ! uniform random number variable
    integer(i4)    :: ssize
    integer(i4), dimension(:), allocatable :: seeds
#endif

    pnum = size(pini)
    !
    ! Check input
    if (size(pmin) /= pnum) stop 'Error DDS: size(pmin) /= size(pini)'
    if (size(pmax) /= pnum) stop 'Error DDS: size(pmax) /= size(pini)'
    if (present(r)) then
       ! r limit (Pertrubation parameter)
       if (r <= 0.0_dp .or. r > 1.0_dp) stop 'Error DDS: Enter DDS perturbation parameter (0.0, 1.0]'
       ir = r
    else
       ir = 0.2_dp
    endif
    if (present(maxiter)) then
       ! max. iteration limit
       if (maxiter < 6 .or. maxIter > 1000000_i8) &
            stop 'Error DDS: Enter max no. of function evals in range 6 to 1000000'
       imaxiter = maxiter
    else
       imaxiter = 1000
    endif
    if (present(maxit)) then
       ! maxit -1 or 1
       if (abs(maxit) /= 1.0_dp) stop 'Error DDS: Enter maxit "-1.0" to maximise or "1.0" to minimise'
       imaxit = maxit
    else
       imaxit = 1.0_dp
    endif
    ! Seed random number
    if (present(seed)) then
       ! user seed limit
       if (seed < 0 .or. seed > 2147483562_i4) then
          stop 'Error DDS: Enter random seed input in range 0 to 2147483562'
       endif
       iseed = seed
    else
       iseed = 0
    endif
#ifdef IMSL
    call RNSET(iseed)
#else
    ! INITIALIZATION OF RANDOM SEEDS	
    ! set the seed so that results easily replicated
    call random_seed(SIZE=ssize)  ! size is number elements in array, 2 usually
    ALLOCATE(seeds(ssize))
    seeds(1)=iseed
    do j=1,ssize ! define initial set of seeds
       seeds(j)=max(1,iseed-j)
    end do
    call random_seed(PUT=seeds(1:ssize)) ! set generator initially
    ! define seeds randomly based on user seed
    seeds(1)=iseed
    do j=2,ssize
       call random_number(harvest) ! get one uniform random number
       seeds(j)=INT(harvest*2147483398.0) ! limit of seed 2 is constant
    end do
    ! Now have replicable, random set of seed arrays, set the generator:
    call random_seed(PUT=seeds(1:ssize))
    deallocate(seeds)
    ! Note - there are probably easier ways to do this but above works...
    ! DONE random seed setting
#endif

    ! Evaluate initial solution and return objective function value (of_value)
    !  and Initialise the other variables (e.g. of_best)
    of_value = obj_func(pini)

    ! maxit is 1.0 for MIN problems, -1 for MAX problems
    of_new = imaxit * of_value

    ! accumulate DDS initialization outputs
    iobj    = 1
    DDS     = pini
    of_best = of_new
    !
    ! Code below is now the DDS algorithm as presented in Figure 1 of Tolson and Shoemaker (2007)
    !
    do i=1, imaxiter-1
       ! Determine Decision Variable (DV) selected for perturbation:
       Pn = 1.0_dp - log(real(i,dp))/log(real(imaxiter-1,dp)) ! probability each DV selected
       dvn_count = 0                                           ! counter for how many DVs selected for perturbation
       pnew = DDS                                                ! define pnew initially as best current solution
       !
       do j=1, pnum
#ifdef IMSL
          call D_RNUN(ranval)                                    ! selects next uniform random number in sequence
#else
          call random_number(ranval)
#endif
          if(ranval(1) < Pn) then                                    ! jth DV selected for perturbation
             dvn_count = dvn_count + 1
             ! call 1-D perturbation function to get new DV value (new_value)
             call neigh_value(DDS(j), pmin(j), pmax(j), ir, new_value) 
             pnew(j) = new_value 
          end if
       end do
       !
       if (dvn_count == 0) then                                 ! no DVs selected at random, so select one
#ifdef IMSL
          call D_RNUN(ranval)                                    ! selects next uniform random number in sequence
#else
          call random_number(ranval)
#endif
          dv = ceiling(real(pnum,dp) * ranval(1))                     ! index for one DV   
          ! call 1-D perturbation function to get new DV value (new_value):
          call neigh_value(DDS(dv), pmin(dv), pmax(dv), ir, new_value)
          pnew(dv) = new_value                                    ! change relevant DV value in stest
       end if
       !
       ! Evaluate obj function value (of_value) for stest
       of_value = obj_func(pnew)
       of_new = imaxit * of_value                                  ! maxit handles min(= 1) and max(= -1) problems 
       ! update current best solution
       if (of_new <= of_best) then                                               
          of_best = of_new
          DDS = pnew
       end if
       ! accumulate DDS search history
       iobj = i + 1
       !
    end do
    close(89)
    !
    if (present(nobj)) nobj = iobj
    !
  end function DDS


  function DDSxy(obj_func, pini, pmin, pmax, xx, yy, r, seed, maxiter, maxit, nobj)

    use mo_kind, only: i4, i8, dp
#ifdef IMSL
    use RNSET_INT
    use RNUN_INT
#endif

    implicit none

    INTERFACE
       FUNCTION obj_func(pp, xx, yy)
         USE mo_kind, ONLY: dp
         IMPLICIT NONE
         REAL(dp), DIMENSION(:), INTENT(IN) :: pp
         REAL(dp), DIMENSION(:), INTENT(IN) :: xx
         REAL(dp), DIMENSION(:), INTENT(IN) :: yy
         REAL(dp) :: obj_func
       END FUNCTION obj_func
    END INTERFACE
    real(dp),    dimension(:), intent(in)  :: pini    ! inital value of decision variables 
    real(dp),    dimension(:), intent(in)  :: pmin    ! Min value of decision variables 
    real(dp),    dimension(:), intent(in)  :: pmax    ! Max value of decision variables 
    real(dp),    dimension(:), intent(in)  :: xx      ! First argument to obj_func
    real(dp),    dimension(:), intent(in)  :: yy      ! Second argument fo obj_func
    real(dp),        optional, intent(in)  :: r       ! DDS perturbation parameter (-> 0.2 by default)
    integer(i4),     optional, intent(in)  :: seed    ! User seed to initialise the random number generator
    integer(i8),     optional, intent(in)  :: maxiter ! Maximum number of iteration or function evaluation
    real(dp),        optional, intent(in)  :: maxit   ! Maximization or minimization of function
    integer(i8),     optional, intent(out) :: nobj ! Counter for Obj. function evaluations
    real(dp), dimension(size(pini)) :: DDSxy ! Best value of decision variables

    ! Local variables
    integer(i4)    :: pnum     ! Total number of decision variables  
    integer(i4)    :: iseed    ! User given seed
    integer(i8)    :: imaxiter ! Maximum number of iteration or function evaluation
    integer(i8)    :: iobj     ! Counter for Obj. function evaluations
    real(dp)       :: ir       ! DDS perturbation parameter
    real(dp)       :: imaxit   ! Maximization or minimization of function
    real(dp)       :: of_value, of_new, of_best, Pn, new_value
    real(dp), dimension(size(pini)) :: pnew   ! Test value of decision variables (dummy) 
    real(dp), dimension(1)          :: ranval
    integer(i8)    :: i, k, dvn_count, dv       
    integer(i4)    :: j
#ifndef IMSL
    real(dp)       :: harvest  ! uniform random number variable
    integer(i4)    :: ssize
    integer(i4), dimension(:), allocatable :: seeds
#endif

    pnum = size(pini)
    !
    ! Check input
    if (size(pmin) /= pnum) stop 'Error DDSxy: size(pmin) /= size(pini)'
    if (size(pmax) /= pnum) stop 'Error DDSxy: size(pmax) /= size(pini)'
    if (size(xx) /= size(yy)) stop 'Error DDSxy: size(xx) /= size(yy)'
    if (present(r)) then
       ! r limit (Pertrubation parameter)
       if (r <= 0.0_dp .or. r > 1.0_dp) stop 'Error DDSxy: Enter DDS perturbation parameter (0.0, 1.0]'
       ir = r
    else
       ir = 0.2_dp
    endif
    if (present(maxiter)) then
       ! max. iteration limit
       if (maxiter < 6 .or. maxIter > 1000000_i8) &
            stop 'Error DDSxy: Enter max no. of function evals in range 6 to 1000000'
       imaxiter = maxiter
    else
       imaxiter = 1000
    endif
    if (present(maxit)) then
       ! maxit -1 or 1
       if (abs(maxit) /= 1.0_dp) stop 'Error DDSxy: Enter maxit "-1.0" to maximise or "1.0" to minimise'
       imaxit = maxit
    else
       imaxit = 1.0_dp
    endif
    ! Seed random number
    if (present(seed)) then
       ! user seed limit
       if (seed < 0 .or. seed > 2147483562_i4) then
          stop 'Error DDSxy: Enter random seed input in range 0 to 2147483562'
       endif
       iseed = seed
    else
       iseed = 0
    endif
#ifdef IMSL
    call RNSET(iseed)
#else
    ! INITIALIZATION OF RANDOM SEEDS	
    ! set the seed so that results easily replicated
    call random_seed(SIZE=ssize)  ! size is number elements in array, 2 usually
    ALLOCATE(seeds(ssize))
    seeds(1)=iseed
    do j=1,ssize ! define initial set of seeds
       seeds(j)=max(1,iseed-j)
    end do
    call random_seed(PUT=seeds(1:ssize)) ! set generator initially
    ! define seeds randomly based on user seed
    seeds(1)=iseed
    do j=2,ssize
       call random_number(harvest) ! get one uniform random number
       seeds(j)=INT(harvest*2147483398.0) ! limit of seed 2 is constant
    end do
    ! Now have replicable, random set of seed arrays, set the generator:
    call random_seed(PUT=seeds(1:ssize))
    deallocate(seeds)
    ! Note - there are probably easier ways to do this but above works...
    ! DONE random seed setting
#endif

    ! Evaluate initial solution and return objective function value (of_value)
    !  and Initialise the other variables (e.g. of_best)
    of_value = obj_func(pini, xx, yy)

    ! maxit is 1.0 for MIN problems, -1 for MAX problems
    of_new = imaxit * of_value

    ! accumulate DDS initialization outputs
    iobj    = 1
    DDSxy     = pini
    of_best = of_new
    !
    ! Code below is now the DDS algorithm as presented in Figure 1 of Tolson and Shoemaker (2007)
    !
    do i=1, imaxiter-1
       ! Determine Decision Variable (DV) selected for perturbation:
       Pn = 1.0_dp - log(real(i,dp))/log(real(imaxiter-1,dp)) ! probability each DV selected
       dvn_count = 0                                           ! counter for how many DVs selected for perturbation
       pnew = DDSxy                                                ! define pnew initially as best current solution
       !
       do j=1, pnum
#ifdef IMSL
          call D_RNUN(ranval)                                    ! selects next uniform random number in sequence
#else
          call random_number(ranval)
#endif
          if(ranval(1) < Pn) then                                    ! jth DV selected for perturbation
             dvn_count = dvn_count + 1
             ! call 1-D perturbation function to get new DV value (new_value)
             call neigh_value(DDSxy(j), pmin(j), pmax(j), ir, new_value) 
             pnew(j) = new_value 
          end if
       end do
       !
       if (dvn_count == 0) then                                 ! no DVs selected at random, so select one
#ifdef IMSL
          call D_RNUN(ranval)                                    ! selects next uniform random number in sequence
#else
          call random_number(ranval)
#endif
          dv = ceiling(real(pnum,dp) * ranval(1))                     ! index for one DV   
          ! call 1-D perturbation function to get new DV value (new_value):
          call neigh_value(DDSxy(dv), pmin(dv), pmax(dv), ir, new_value)
          pnew(dv) = new_value                                    ! change relevant DV value in stest
       end if
       !
       ! Evaluate obj function value (of_value) for stest
       of_value = obj_func(pnew, xx, yy)
       of_new = imaxit * of_value                                  ! maxit handles min(= 1) and max(= -1) problems 
       ! update current best solution
       if (of_new <= of_best) then                                               
          of_best = of_new
          DDSxy = pnew
       end if
       ! accumulate DDS search history
       iobj = i + 1
       !
    end do
    close(89)
    !
    if (present(nobj)) nobj = iobj
    !
  end function DDSxy

  ! ------------------------------------------------------------------

  ! Purpose is to generate a neighboring decision variable value for a single
  !  decision variable value being perturbed by the DDS optimization algorithm.
  !  New DV value respects the upper and lower DV bounds.
  !  Coded by Bryan Tolson, Nov 2005.
  ! 
  ! I/O variable definitions:
  !  x_cur     - current decision variable (DV) value
  !  x_min     - min DV value
  !  x_max     - max DV value
  !  r         - the neighborhood perturbation factor
  !  new_value - new DV variable value (within specified min and max)

  subroutine neigh_value(x_cur, x_min, x_max, r, new_value)

    use mo_kind, only: dp
#ifdef IMSL
    USE RNNOA_INT
#endif

    implicit none
    !
    real(dp), intent(in)   :: x_cur, x_min, x_max, r
    real(dp), intent(out)  :: new_value
    real(dp)               :: x_range
    real(dp), dimension(1) :: zvalue
#ifndef IMSL
    real(dp)               :: work1, work2, work3, ranval
#endif

    x_range = x_max - x_min

    !
    ! generate a standard normal random variate (zvalue)
    ! perturb current value with normal random variable
    ! generates a standard normal random deviate
    !
#ifdef IMSL
    call D_RNNOA(zvalue) ! ISML Stat Library 2 routine - Acceptance/rejection 
#else
    !	Below returns a standard Gaussian random number based upon Numerical recipes gasdev and 
    !	Marsagalia-Bray Algorithm
    Work3=2.0 
    Do While (Work3.ge.1.0.or.Work3.eq.0.0)
       call random_number(ranval) ! get one uniform random number
       Work1 = 2.0 * real(ranval,dp) - 1.0
       call random_number(ranval) ! get one uniform random number
       Work2 = 2.0 * real(ranval,dp) - 1.0
       Work3 = Work1 * Work1 + Work2 * Work2
    enddo
    Work3 = ((-2.0 * log(Work3)) / Work3)**0.5  ! natural log	
    ! pick one of two deviates at random (don't worry about trying to use both):
    call random_number(ranval) ! get one uniform random number
    IF (ranval.LT.0.5) THEN
       zvalue(1) = Work1 * Work3
    else
       zvalue(1) = Work2 * Work3
    endif
#endif

    ! calculate new decision variable value:			   	
    new_value = x_cur + zvalue(1)*r*x_range

    !  check new value is within DV bounds.  If not, bounds are reflecting.
    if(new_value < x_min) then
       new_value = x_min + (x_min - new_value)
       if (new_value > x_max) then
          ! if reflection goes past x_max then value should be x_min since 
          ! without reflection the approach goes way past lower bound.  
          ! This keeps x close to lower bound when x_cur is close to lower bound
          !  Practically speaking, this should never happen with r values <0.3
          new_value = x_min
       end if
    else if(new_value > x_max) then
       new_value = x_max - (new_value - x_max)
       if (new_value < x_min) then
          ! if reflection goes past x_min then value should be x_max for same reasons as above
          new_value = x_max
       endif
    endif
    !
  end subroutine neigh_value

  ! ------------------------------------------------------------------

end module mo_dds

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


  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         DDS

  !     PURPOSE
  !         Searches Minimum or Maximum of a user-specified function using
  !         Dynamically Dimensioned Search (DDS)
  !
  !         DDS is an n-dimensional continuous global optimization algorithm.
  !         It is coded as a minimizer but one can give maxit=True in a maximization problem,
  !         so that the algorithm minimizes the negative of the objective function F=(-1*F).
  !
  !         The function to be minimized is the first argument of DDS and must be defined as
  !             function func(p)
  !               use mo_kind, only: dp
  !               implicit none
  !               real(dp), dimension(:), intent(in) :: p
  !               real(dp) :: func
  !             end function func
  !
  !     CALLING SEQUENCE
  !         popt = DDS(obj_func, pini, prange, r=r, seed=seed, maxiter=maxiter, maxit=maxit)

  !     INDENT(IN)
  !         real(dp) :: obj_func(p,xx,yy)       Function on which to search the minimum
  !         real(dp) :: pini(:)                 inital value of decision variables
  !         real(dp) :: prange(size(pini),2)    Min/max range of decision variables

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(dp) :: DDS(size(pini))         The parameters of the point which is estimated to minimize the function.

  !     INDENT(IN), OPTIONAL
  !         real(dp)    :: r                    DDS perturbation parameter (default: 0.2)
  !         integer(i8) :: seed                 User seed to initialise the random number generator (default: None)
  !         integer(i8) :: maxiter              Maximum number of iteration or function evaluation (default: 1000)
  !         logical     :: maxit                Maximization (.True.) or minimization (.False.) of function
  !                                             (default: .False.)

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None.

  !     EXAMPLE
  !         dv_range(:,1) = (/ -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0 /)
  !         dv_range(:,2) = (/ 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0 /)
  !         dv_ini        = (/ -.226265E+01, -.130187E+01, -.151219E+01, 0.133983E+00, 0.988159E+00, &
  !                            -.495074E+01, -.126574E+02, 0.572684E+00, 0.303864E+01, 0.343031E+01 /)
  !         dv_opt = DDS(griewank, dv_ini, dv_range)
  !         -> see also example in test directory

  !     LITERATURE
  !         Tolson, B. A., and C. A. Shoemaker (2007)
  !             Dynamically dimensioned search algorithm for computationally efficient watershed
  !             model calibration, Water Resour. Res., 43, W01413, doi:10.1029/2005WR004723.

  !     HISTORY
  !         Written original Bryan Tolson, Feb 2007 - DDS v1.1
  !         Modified, Rohini Kumar, Feb 2008
  !         Modified, Matthias Cuntz, Jul 2012 - module

  function DDS(obj_func, pini, prange, r, seed, maxiter, maxit)

    use mo_kind,    only: i4, i8, dp
    use mo_xor4096, only: xor4096, xor4096g

    implicit none

    INTERFACE
       function obj_func(pp)
         use mo_kind, only: dp
         implicit none
         real(dp), dimension(:), intent(in) :: pp
         real(dp) :: obj_func
       end function obj_func
    END INTERFACE
    real(dp),    dimension(:),   intent(in)  :: pini    ! inital value of decision variables 
    real(dp),    dimension(:,:), intent(in)  :: prange  ! Min/max values of decision variables 
    real(dp),          optional, intent(in)  :: r       ! DDS perturbation parameter (-> 0.2 by default)
    integer(i8),       optional, intent(in)  :: seed    ! User seed to initialise the random number generator
    integer(i8),       optional, intent(in)  :: maxiter ! Maximum number of iteration or function evaluation
    logical,           optional, intent(in)  :: maxit   ! Maximization or minimization of function
    real(dp),          dimension(size(pini)) :: DDS     ! Best value of decision variables

    ! Local variables
    integer(i4)    :: pnum     ! Total number of decision variables  
    integer(i8)    :: iseed    ! User given seed
    integer(i8)    :: imaxiter ! Maximum number of iteration or function evaluation
    real(dp)       :: ir       ! DDS perturbation parameter
    real(dp)       :: imaxit   ! Maximization or minimization of function
    real(dp)       :: of_new, of_best, Pn, new_value ! intermediate results
    real(dp), dimension(size(pini)) :: pnew                    ! Test value of decision variables
    real(dp)       :: ranval                                   ! random value
    integer(i8)    :: i                                        ! maxiter=i8
    integer(i4)    :: j, dvn_count, dv                         ! pnum=i4
    integer, dimension(8) :: sdate                             ! date_and_time return

    !
    ! Check input
    pnum = size(pini)
    if (size(prange,1) /= pnum) stop 'Error DDS: size(prange,1) /= size(pini)'
    if (size(prange,2) /= 2)    stop 'Error DDS: size(prange,2) /= 2'
    ! r Perturbation parameter
    ir = 0.2_dp
    if (present(r)) ir = r
    if (ir <= 0.0_dp .or. ir > 1.0_dp) stop 'Error DDS: DDS perturbation parameter (0.0, 1.0]'
    ! max. iteration
    imaxiter = 1000
    if (present(maxiter)) imaxiter = maxiter
    if (imaxiter < 6) stop 'Error DDS: max function evals must be minimum 6'
    ! Min or max objective function
    imaxit = 1.0_dp
    if (present(maxit)) then
       if (maxit) imaxit = -1.0_dp
    endif
    ! Given seed
    iseed = 0
    if (present(seed)) iseed = seed
    iseed = max(iseed, 0_i8)
    !
    ! Seed random numbers
    if (iseed == 0) then
       call date_and_time(values=sdate)
       iseed = sdate(1)*31536000000_i8 + sdate(2)*2592000000_i8 + sdate(3)*86400000_i8 + &
            sdate(5)*3600000_i8 + sdate(6)*60000_i8 + sdate(7)*1000_i8 + sdate(8)
       call xor4096(iseed,ranval)
       call xor4096g(iseed,ranval)
    else
       call xor4096(iseed,ranval)
       call xor4096g(iseed,ranval)
    endif

    ! Evaluate initial solution and return objective function value
    ! and Initialise the other variables (e.g. of_best)
    ! imaxit is 1.0 for MIN problems, -1 for MAX problems
    DDS     = pini
    of_new  = imaxit * obj_func(pini)
    of_best = of_new

    !
    ! Code below is now the DDS algorithm as presented in Figure 1 of Tolson and Shoemaker (2007)
    !
    do i=1, imaxiter-1
       ! Determine Decision Variable (DV) selected for perturbation:
       Pn        = 1.0_dp - log(real(i,dp))/log(real(imaxiter-1,dp)) ! probability each DV selected
       dvn_count = 0                                                 ! counter for how many DVs selected for perturbation
       pnew      = DDS                                               ! define pnew initially as best current solution
       !
       ! Step 3 of Fig 1 of Tolson and Shoemaker (2007)
       do j=1, pnum
          call xor4096(0_i8,ranval)                           ! selects next uniform random number in sequence
          ! Step 4 of Fig 1 of Tolson and Shoemaker (2007)
          if (ranval < Pn) then                               ! jth DV selected for perturbation
             dvn_count = dvn_count + 1
             ! call 1-D perturbation function to get new DV value (new_value)
             call neigh_value(DDS(j), prange(j,1), prange(j,2), ir, new_value) 
             pnew(j) = new_value 
          end if
       end do
       !
       ! Step 3 of Fig 1 of Tolson and Shoemaker (2007) in case {N} empty
       if (dvn_count == 0) then                               ! no DVs selected at random, so select one
          call xor4096(0_i8,ranval)                           ! selects next uniform random number in sequence
          dv = ceiling(real(pnum,dp) * ranval)                     ! index for one DV   
          ! call 1-D perturbation function to get new DV value (new_value):
          call neigh_value(DDS(dv), prange(dv,1), prange(dv,2), ir, new_value)
          pnew(dv) = new_value                                ! change relevant DV value in stest
       end if
       !
       ! Step 5 of Fig 1 of Tolson and Shoemaker (2007)
       ! Evaluate obj function value for test
       of_new = imaxit * obj_func(pnew)                       ! imaxit handles min(=1) and max(=-1) problems 
       ! update current best solution
       if (of_new <= of_best) then                                               
          of_best = of_new
          DDS     = pnew
       end if
    end do
    !
  end function DDS


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

    use mo_kind,    only: i8, dp
    use mo_xor4096, only: xor4096g

    implicit none
    !
    real(dp), intent(in)   :: x_cur, x_min, x_max, r
    real(dp), intent(out)  :: new_value
    real(dp)               :: x_range
    real(dp) :: zvalue

    x_range = x_max - x_min

    ! generate a standard normal random variate (zvalue)
    call xor4096g(0_i8,zvalue)

    ! calculate new decision variable value:			   	
    new_value = x_cur + zvalue*r*x_range

    !  check new value is within bounds. If not, bounds are reflecting.
    if (new_value < x_min) then
       new_value = x_min + (x_min - new_value)
       if (new_value > x_max) then
          ! if reflection goes past x_max then value should be x_min since 
          ! without reflection the approach goes way past lower bound.
          ! This keeps x close to lower bound when x_cur is close to lower bound
          ! Practically speaking, this should never happen with r values <0.3.
          new_value = x_min
       end if
    else if (new_value > x_max) then
       new_value = x_max - (new_value - x_max)
       if (new_value < x_min) then
          ! if reflection goes past x_min then value should be x_max for same reasons as above.
          new_value = x_max
       endif
    endif

  end subroutine neigh_value

  ! ------------------------------------------------------------------

end module mo_dds

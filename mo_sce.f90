!> \file mo_sce.f90

!> \brief Shuffled Complex Evolution optimization algorithm.

!> \details Optimization algorithm using Shuffled Complex Evolution strategy. 
!>          Original version 2.1 of Qingyun Duan (1992) rewritten in Fortran 90.

!> \authors Juliane Mai
!> \date Feb 2013

MODULE mo_sce

  ! This module is a template for the UFZ CHS Fortran library.

  ! Written  Juliane Mai, Feb 2013

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

  ! Copyright 2011-2012 Matthias Cuntz, Juliane Mai


  ! Note on Numerical Recipes License
  ! ---------------------------------
  ! Be aware that some code is under the Numerical Recipes License 3rd
  ! edition <http://www.nr.com/aboutNR3license.html>

  ! The Numerical Recipes Personal Single-User License lets you personally
  ! use Numerical Recipes code ("the code") on any number of computers,
  ! but only one computer at a time. You are not permitted to allow anyone
  ! else to access or use the code. You may, under this license, transfer
  ! precompiled, executable applications incorporating the code to other,
  ! unlicensed, persons, providing that (i) the application is
  ! noncommercial (i.e., does not involve the selling or licensing of the
  ! application for a fee), and (ii) the application was first developed,
  ! compiled, and successfully run by you, and (iii) the code is bound
  ! into the application in such a manner that it cannot be accessed as
  ! individual routines and cannot practicably be unbound and used in
  ! other programs. That is, under this license, your application user
  ! must not be able to use Numerical Recipes code as part of a program
  ! library or "mix and match" workbench.

  ! Businesses and organizations that purchase the disk or code download,
  ! and that thus acquire one or more Numerical Recipes Personal
  ! Single-User Licenses, may permanently assign those licenses, in the
  ! number acquired, to individual employees. Such an assignment must be
  ! made before the code is first used and, once made, it is irrevocable
  ! and can not be transferred.

  ! If you do not hold a Numerical Recipes License, this code is only for
  ! informational and educational purposes but cannot be used.

  IMPLICIT NONE

  PUBLIC :: sce        ! sce optimization

  ! ------------------------------------------------------------------

  !     NAME
  !         sce

  !     PURPOSE
  !>        \brief Shuffled Complex Evolution algorithm for global optimization.
  !
  !>        \details  Shuffled Complex Evolution method for global optimization/n
  !>                  -- version 2.1/n
  !>                  /n
  !>                  by Qingyun Duan/n
  !>                  Department of Hydrology & Water Resources/n
  !>                  University of Arizona, Tucson, AZ 85721/n
  !>                  (602) 621-9360, email: duan@hwr.arizona.edu/n
  !>                  /n
  !>                  Written by Qingyun Duan,   Oct 1990./n
  !>                  Revised by Qingyun Duan,   Aug 1991./n
  !>                  Revised by Qingyun Duan,   Apr 1992./n
  !>                  /n
  !>                  Re-written by Juliane Mai, Feb 2013./n
  !>                  /n
  !>                  Statement by Qingyun Duan:/n
  !>                  ----------------------------------/n
  !>                  /n
  !>                     This general purpose global optimization program is developed at
  !>                     the Department of Hydrology & Water Resources of the University
  !>                     of Arizona.  Further information regarding the SCE-UA method can
  !>                     be obtained from Dr. Q. Duan, Dr. S. Sorooshian or Dr. V.K. Gupta
  !>                     at the address and phone number listed above.  We request all
  !>                     users of this program make proper reference to the paper entitled
  !>                     'Effective and Efficient Global Optimization for Conceptual
  !>                     Rainfall-runoff Models' by Duan, Q., S. Sorooshian, and V.K. Gupta,
  !>                     Water Resources Research, Vol 28(4), pp.1015-1031, 1992.
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>       \return     real(dp) :: bestx &mdash; best parameter set found
  !
  !     RESTRICTIONS
  !>       \note Function will be minimized./n
  !>             No mask of parameters implemented yet./n
  !>             No maximization implemented yet. /n
  !
  !     EXAMPLE
  !         None
  !         -> see also example in test directory

  !     LITERATURE  
  !         Duan, Q., S. Sorooshian, and V.K. Gupta - 
  !             Effective and Efficient Global Optimization for Conceptual Rainfall-runoff Models
  !             Water Resources Research, Vol 28(4), pp.1015-1031, 1992.

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Feb 2013
  !         Modified, 

  ! ------------------------------------------------------------------

  PRIVATE

  ! chkcst   ! check if a trial point satisfies all constraints.
  ! ran1     ! returns uniform distributed random numbers
  ! gasdev   ! return gaussian distributed random numbers
  ! comp     ! reduces the number of complexes and points in a population
  ! parstt   ! checks for parameter convergence
  ! getpnt   ! generated a new point within its range
  ! cce      ! generates new point(s) from sub-complex

  ! ------------------------------------------------------------------

CONTAINS

  function sce(functn,pini,prange,                        & ! IN
       mymaxn,mymaxit,mykstop,mypcento,myseed,            & ! Optional IN
       myngs,mynpg,mynps,mynspl,mymings,myiniflg,myprint, & ! Optional IN
       myalpha, mybeta,                                   & ! Optional IN
       bestf,neval,history                                & ! Optional OUT
    ) result(bestx)

    use mo_kind,    only: i4, i8, dp
    use mo_sort,    only: sort
    use mo_xor4096, only: get_timeseed

    implicit none
    ! 
    interface
       function functn(paraset)
         ! calculates the cost function at a certain parameter set paraset
         use mo_kind, only: dp
         real(dp), dimension(:),    intent(in)       :: paraset
         real(dp)                                    :: functn
       end function functn
    end interface
    real(dp), dimension(:),   intent(in)             :: pini        ! initial parameter set
    real(dp), dimension(:,:), intent(in)             :: prange      ! lower and upper bound on parameters
    !
    ! algorithmic control & convergence parameter 
    integer(i8), optional,               intent(in)  :: mymaxn      ! max no. of trials allowed before optimization is terminated
    !                                                               !     DEFAULT: 1000_i8
    logical,     optional,               intent(in)  :: mymaxit     ! maximization (.true.) or minimization (.false.) of function
    !                                                               !     DEFAULT: false
    integer(i4), optional,               intent(in)  :: mykstop     ! number of shuffling loops in which the criterion value must
    !                                                               !     change by given percentage before optimiz. is terminated
    !                                                               !     DEFAULT: 10_i4
    real(dp),    optional,               intent(in)  :: mypcento    ! percentage by which the criterion value must change in
    !                                                               !     given number of shuffling loops
    !                                                               !     DEFAULT: 0.0001_dp
    integer(i8), optional,               intent(in)  :: myseed      ! initial random seed
    !                                                               !     DEFAULT: get_timeseed
    integer(i4), optional,               intent(in)  :: myngs       ! number of complexes in the initial population
    !                                                               !     DEFAULT: 2_i4
    integer(i4), optional,               intent(in)  :: mynpg       ! number of points in each complex
    !                                                               !     DEFAULT: 2*n+1
    integer(i4), optional,               intent(in)  :: mynps       ! number of points in a sub-complex
    !                                                               !     DEFAULT: n+1
    integer(i4), optional,               intent(in)  :: mynspl      ! number of evolution steps allowed for each complex before
    !                                                               !     complex shuffling
    !                                                               !     DEFAULT: 2*n+1
    integer(i4), optional,               intent(in)  :: mymings     ! minimum number of complexes required, if the number of 
    !                                                               !     complexes is allowed to reduce as the 
    !                                                               !     optimization proceeds
    !                                                               !     DEFAULT: ngs = number of complexes in initial population
    integer(i4), optional,               intent(in)  :: myiniflg    ! flag on whether to include the initial point in population
    !                                                               !     0, not included
    !                                                               !     1, included (DEFAULT)
    integer(i4), optional,               intent(in)  :: myprint     ! flag for controlling print-out after each shuffling loop
    !                                                               !     0, print information on the best point of the population
    !                                                               !     1, print information on every point of the population
    !                                                               !     2, no printing (DEFAULT)
    real(dp),    optional,               intent(in)  :: myalpha     ! parameter for reflection  of points in complex
    !                                                               !     DEFAULT: 0.8_dp
    real(dp),    optional,               intent(in)  :: mybeta      ! parameter for contraction of points in complex
    !                                                               !     DEFAULT: 0.45_dp
    real(dp),    optional,               intent(out) :: bestf       ! function value of bestx(.)
    integer(i4), optional,               intent(out) :: neval       ! number of function evaluations
    real(dp),    optional, &
         dimension(:), allocatable,      intent(out) :: history     ! history of best function values after each iteration
    real(dp), dimension(size(pini,1))                :: bestx       ! best point at current shuffling loop (is RETURN)
    !
    ! optionals transfer
    real(dp), dimension(size(pini,1))                :: bl          ! lower bound on parameters
    real(dp), dimension(size(pini,1))                :: bu          ! upper bound on parameters
    integer(i8)                                      :: maxn        ! max no. of trials allowed before optimization is terminated
    logical                                          :: maxit       ! minimization (false) or maximization (true)
    integer(i4)                                      :: kstop       ! number of shuffling loops in which the criterion value
    !                                                               !     must change
    real(dp)                                         :: pcento      ! percentage by which the criterion value must change
    integer(i8)                                      :: iseed       ! initial random seed  
    integer(i4)                                      :: ngs         ! number of complexes in the initial population
    integer(i4)                                      :: npg         ! number of points in each complex
    integer(i4)                                      :: nps         ! number of points in a sub-complex
    integer(i4)                                      :: nspl        ! number of evolution steps allowed for each complex before 
    !                                                               !     complex shuffling
    integer(i4)                                      :: mings       ! minimum number of complexes required
    integer(i4)                                      :: iniflg      ! flag on whether to include the initial point in population
    integer(i4)                                      :: iprint      ! flag for controlling print-out after each shuffling loop
    real(dp)                                         :: alpha       ! parameter for reflection  of points in complex
    real(dp)                                         :: beta        ! parameter for contraction of points in complex
    real(dp), dimension(:), allocatable              :: history_tmp ! history of best function values after each iteration
    real(dp)                                         :: bestf_tmp   ! function value of bestx(.)
    !
    ! local variables
    integer(i4)                                      :: nopt        ! number of parameters to be optimized
    integer(i4)                                      :: npt         ! total number of points in initial population (npt=ngs*npg)
    real(dp)                                         :: fpini       ! function value at initial point
    real(dp),    dimension(:,:), allocatable         :: x           ! coordinates of points in the population 
    real(dp),    dimension(:),   allocatable         :: xf          ! function values of x                    
    real(dp),    dimension(:),   allocatable         :: xx          ! coordinates of a single point in x
    real(dp),    dimension(:,:), allocatable         :: cx          ! coordinates of points in a complex      
    real(dp),    dimension(:),   allocatable         :: cf          ! function values of cx
    real(dp),    dimension(:,:), allocatable         :: s           ! coordinates of points in the current sub-complex !simplex
    real(dp),    dimension(:),   allocatable         :: sf          ! function values of s
    real(dp),    dimension(:),   allocatable         :: worstx      ! worst point at current shuffling loop
    real(dp)                                         :: worstf      ! function value of worstx(.)
    real(dp),    dimension(:),   allocatable         :: xnstd       ! standard deviation of parameters in the population
    real(dp)                                         :: gnrng       ! normalized geometric mean of parameter ranges
    integer(i4), dimension(:),   allocatable         :: lcs         ! indices locating position of s(.,.) in x(.,.)  
    real(dp),    dimension(:),   allocatable         :: bound       ! length of range of ith variable being optimized
    real(dp),    dimension(:),   allocatable         :: unit        ! ?????
    integer(i4)                                      :: ngs1        ! number of complexes in current population
    integer(i4)                                      :: ngs2        ! number of complexes in last population
    integer(i8)                                      :: iseed1      ! current random seed
    real(dp),    dimension(:),   allocatable         :: criter      ! vector containing the best criterion values of the last
    !                                                               !     (kstop+1) shuffling loops
    integer(i4)                                      :: ipcnvg      ! flag indicating whether parameter convergence is reached
    !                                                               !      (i.e., check if gnrng is less than 0.001)
    !                                                               !      0, parameter convergence not satisfied
    !                                                               !      1, parameter convergence satisfied
    integer(i4)                                      :: nloop
    integer(i4)                                      :: loop
    integer(i4)                                      :: i, j, k, l
    integer(i4)                                      :: lpos
    logical                                          :: lpos_ok     ! for selction of points based on triangular 
    !                                                               ! probability distribution
    integer(i4)                                      :: npt1
    integer(i8)                                      :: icall       ! counter for function evaluations
    integer(i4)                                      :: igs, ipr
    integer(i4)                                      :: k1, k2
    real(dp)                                         :: denomi      ! for checking improvement of last steps
    real(dp)                                         :: timeou      ! for checking improvement of last steps
    real(dp)                                         :: rtiny       ! for checking improvement of last steps
    integer(i4)                                      :: nopt1       ! only for printing of parameter sets
    integer(i4)                                      :: nopt2       ! only for printing of parameter sets
    character(4), dimension(:),  allocatable         :: xname       ! only for printing of parameter sets
    real(dp)                                         :: rand        ! random number
    !
    nopt = size(pini,dim=1)
    bl(:)   = prange(:,1)
    bu(:)   = prange(:,2)
    !
    ! input checking
    if (size(prange,dim=1) .ne. nopt) then
       stop 'sceua: prange has not matching rows'
    end if
    if (size(prange,dim=2) .ne. 2) then
       stop 'sceua: two colums expected for prange'
    end if
    !
    ! optionals checking
    if (present(mymaxn)) then
       if (mymaxn .lt. 2_i4) stop 'sceua: maxn has to be at least 2'
       maxn = mymaxn
    else
       maxn = 1000_i8
    end if
    if (present(mymaxit)) then
       maxit = mymaxit
    else
       maxit = .false.
    end if
    if(present(mykstop)) then
       if (mykstop .lt. 1_i4) stop 'sceua: kstop has to be at least 1'
       kstop = mykstop
    else
       kstop = 10_i4
    end if
    if(present(mypcento)) then
       if (mypcento .lt. 0_dp) stop 'sceua: pcento should be positive'
       pcento = mypcento
    else
       pcento = 0.0001_dp
    end if
    if(present(myseed)) then
       if (myseed .lt. 1_i8) stop 'sceua: seed should be non-negative'
       iseed = myseed
    else
       call get_timeseed(iseed)
    end if
    if(present(myngs)) then
       if (myngs .lt. 1_i4) stop 'sceua: ngs has to be at least 1'
       ngs = myngs
    else
       ngs = 2_i4
    end if
    if(present(mynpg)) then
       if (mynpg .lt. 3_i4) stop 'sceua: npg has to be at least 3'
       npg = mynpg
    else
       npg = 2*nopt + 1
    end if
    if(present(mynps)) then
       if (mynps .lt. 2_i4) stop 'sceua: nps has to be at least 2'
       nps = mynps
    else
       nps = nopt + 1_i4
    end if
    if(present(mynspl)) then
       if (mynspl .lt. 3_i4) stop 'sceua: nspl has to be at least 3'
       nspl = mynspl
    else
       nspl = 2*nopt + 1
    end if
    if(present(mymings)) then
       if (mymings .lt. 1_i4) stop 'sceua: mings has to be at least 1'
       mings = mymings
    else
       mings = ngs  ! no reduction of complexes
    end if
    if(present(myiniflg)) then
       if ( (myiniflg .ne. 1_i4) .and. (myiniflg .ne. 0_i4)) stop 'sceua: iniflg has to be 0 or 1'
       iniflg = myiniflg
    else
       iniflg = 1_i4
    end if
    if(present(myprint)) then
       if ( (myprint .lt. 0_i4) .or. (myprint .gt. 2_i4)) stop 'sceua: iprint has to be 0, 1 or 2'
       iprint = myprint
    else
       iprint = 2_i4  ! no printing
       iprint = 0_i4
    end if
    if(present(myalpha)) then
       alpha = myalpha
    else
       alpha = 0.8_dp
    end if
    if(present(mybeta)) then
       beta = mybeta
    else
       beta = 0.45_dp
    end if
    
    ! allocation of arrays
    allocate(x(ngs*npg,nopt))
    allocate(xf(ngs*npg))
    allocate(xx(nopt))
    allocate(cx(npg,nopt)) 
    allocate(cf(npg))
    allocate(s(nps,nopt))
    allocate(sf(nps))
    allocate(worstx(nopt))
    allocate(xnstd(nopt))
    allocate(lcs(nps))
    allocate(bound(nopt))
    allocate(unit(nopt))
    allocate(criter(kstop+1))   ! kstop+1
    allocate(xname(nopt))
    allocate(history_tmp(maxn))

    if (iprint .lt. 2) then
       write (*,*) ' enter the sceua subroutine --- '   
    end if
    !
    !  initialize variables
    do i=1,nopt
       xname(i) = 'XXXX'
    end do

    nloop = 0_i4
    loop  = 0_i4
    igs   = 0_i4
    nopt1 = 8_i4
    if (nopt .lt. 8)  nopt1 = nopt
    nopt2 = 12_i4
    if (nopt .lt. 12) nopt2 = nopt
    !
    !  initialize random seed to a negative integer
    iseed1 = -abs(iseed)
    !
    !  compute the total number of points in initial popualtion
    npt = ngs * npg
    ngs1 = ngs
    npt1 = npt
    !
    if (iprint .lt. 2) then
       write(ipr,400)
       write (*,*) ' ***  Evolution Loop Number ',nloop
    end if
    !
    !  compute the bound for parameters being optimized
    do j = 1, nopt
       bound(j) = bu(j) - bl(j)
       unit(j) = 1.0_dp
    end do
    !--------------------------------------------------
    !  compute the function value of the initial point
    !--------------------------------------------------
    ! function evaluation will be counted later...
    if (.not. maxit) then
       fpini = functn(pini)
       history_tmp(1) = fpini
    else
       fpini = -functn(pini)
       history_tmp(1) = -fpini
    end if

    !  print the initial point and its criterion value
    if (iprint .lt. 2) then
       write(ipr,500)
       write(ipr,510) (xname(j),j=1,nopt2)
       if (.not. maxit) then
          write(ipr,520) fpini,(pini(j),j=1,nopt2)
       else
          write(ipr,520) -fpini,(pini(j),j=1,nopt2)
       end if
       if (nopt.gt.12) then
          write(ipr,530) (xname(j),j=13,nopt)
          write(ipr,540) (pini(j),j=13,nopt)
       end if
    end if
    !
    !  generate an initial set of npt1 points in the parameter space
    !  if iniflg is equal to 1, set x(1,.) to initial point pini(.)
    if (iniflg .eq. 1) then
       do j = 1, nopt
          x(1,j) = pini(j)
       end do
       xf(1) = fpini
       !
       !  else, generate a point randomly and set it equal to x(1,.)
    else
       call getpnt(1,bl(1:nopt),bu(1:nopt),unit(1:nopt),bl(1:nopt),iseed1,xx)
       do j=1, nopt
          x(1,j) = xx(j)
       end do
       xf(1) = functn(xx)
       if (.not. maxit) then
          xf(1) = functn(pini)
       else
          xf(1) = -functn(pini)
       end if
    end if
    !
    ! count function evaluation of the first point
    icall = 1_i8
    ! if (icall .ge. maxn) return 
    !
    !  generate npt1-1 random points distributed uniformly in the parameter
    !  space, and compute the corresponding function values
    do i = 2, npt1
       call getpnt(1,bl(1:nopt),bu(1:nopt),unit(1:nopt),bl(1:nopt),iseed1,xx)
       do j = 1, nopt
          x(i,j) = xx(j)
       end do
       if (.not. maxit) then
          xf(i) = functn(xx)
          icall = icall + 1_i8
          history_tmp(icall) = min(history_tmp(icall-1),xf(i))
       else
          xf(i) = -functn(xx)
          icall = icall + 1_i8
          history_tmp(icall) = max(history_tmp(icall-1),-xf(i))
       end if
       
       if (icall .ge. maxn) then
          npt1 = i
          exit
       end if
    end do
    !
    !  arrange the points in order of increasing function value
    call sort_matrix(x(1:npt1,1:nopt),xf(1:npt1))
    !
    !  record the best and worst points
    do j = 1, nopt
       bestx(j) = x(1,j)
       worstx(j) = x(npt1,j)
    end do
    bestf_tmp = xf(1)
    worstf = xf(npt1)
    ! if (.not. maxit) then
    !    bestf_tmp = xf(1)
    !    worstf = xf(npt1)
    ! else
    !    bestf_tmp = -xf(1)
    !    worstf    = -xf(npt1)
    ! end if
    !
    !  compute the parameter range for the initial population
    call parstt(x(1:npt1,1:nopt),bound,xnstd,gnrng,ipcnvg)
    !
    !  print the results for the initial population
    print: if (iprint .lt. 2) then
       write(ipr,600)
       write(ipr,610) (xname(j),j=1,nopt1)
       if (nopt .gt. 8) then 
          write(ipr,620) (xname(j),j=9,nopt)
       end if
       if (.not. maxit) then
          write(ipr,630) nloop,icall,ngs1,bestf_tmp,worstf,gnrng, (bestx(j),j=1,nopt1)
       else
          write(ipr,630) nloop,icall,ngs1,-bestf_tmp,-worstf,gnrng, (bestx(j),j=1,nopt1)
       end if
       if (nopt .gt. 8) then
          write(ipr,640) (bestx(j),j=9,nopt)
       end if
       if (iprint .eq. 1) then
          write(ipr,650) nloop
          do i = 1, npt1
             write(ipr,660) xf(i),(x(i,j),j=1,nopt1)
             if (nopt .gt. 8) then
                write(ipr,640) (x(i,j),j=9,nopt)
             end if
          end do
       end if
    end if print
    !
    ! Maximum number of function evaluations reached?
    if (icall .ge. maxn) then 
       if (iprint .lt. 2) then
          ! maximum trials reached
          write(ipr,800) maxn,loop,igs,nloop
          ! parameter set
          write(ipr,830)
          write(ipr,510) (xname(j),j=1,nopt2)
          if (.not. maxit) then
             write(ipr,520) bestf_tmp,(bestx(j),j=1,nopt2)
          else
             write(ipr,520) -bestf_tmp,(bestx(j),j=1,nopt2)
          end if
          if (nopt .gt. 12) then 
             write(ipr,530) (xname(j),j=13,nopt)
             write(ipr,540) (bestx(j),j=13,nopt)
          end if
       end if
       if (present(neval)) neval = icall
       if (present(bestf) .and. .not. maxit) bestf = bestf_tmp
       if (present(bestf) .and. maxit) bestf = -bestf_tmp
       if (present(history)) then
          allocate(history(icall))
          history(:) = history_tmp(1:icall)
       end if
       ! -----------------------
       ! Abort subroutine
       ! -----------------------
       return
    end if
    !
    ! Feasible parameter space converged?
    if (ipcnvg .eq. 1) then 
       if (iprint .lt. 2) then
          ! converged because feasible parameter space small
          write(ipr,820) gnrng*100.
          ! parameter set
          write(ipr,830)
          write(ipr,510) (xname(j),j=1,nopt2)
          if (.not. maxit) then
             write(ipr,520) bestf_tmp,(bestx(j),j=1,nopt2)
          else
             write(ipr,520) -bestf_tmp,(bestx(j),j=1,nopt2)
          end if
          if (nopt .gt. 12) then
             write(ipr,530) (xname(j),j=13,nopt)
             write(ipr,540) (bestx(j),j=13,nopt)
          end if
       end if
       if (present(neval)) neval = icall
       if (present(bestf) .and. .not. maxit) bestf = bestf_tmp
       if (present(bestf) .and. maxit)       bestf = -bestf_tmp
       if (present(history)) then
          allocate(history(icall))
          history(:) = history_tmp(1:icall)
       end if
       ! -----------------------
       ! Abort subroutine
       ! -----------------------
       return
    end if
    !
    !  begin the main loop 
    mainloop:do while (icall .lt. maxn)
       nloop = nloop + 1
       !
       if (iprint .lt. 2) then
          write (*,*) ' ***  Evolution Loop Number ',nloop
       end if
       !
       !  begin loop on complexes
       !     <beta> loop from duan(1993)
       comploop: do igs = 1, ngs1
          !
          !  assign points into complexes
          do k1 = 1, npg
             k2 = (k1-1) * ngs1 + igs
             do j = 1, nopt
                cx(k1,j) = x(k2,j)
             end do
             cf(k1) = xf(k2)
          end do
          !
          !  begin inner loop - random selection of sub-complexes 
          !     <alpha> loop from duan(1993)
          subcomploop: do loop = 1, nspl
             !
             !  choose a sub-complex (nps points) according to a linear
             !  probability distribution
             !
             ! if number of points in subcomplex (nps) = number of points in complex (npg)
             ! --> select all
             if (nps .eq. npg) then
                do k = 1, nps
                   lcs(k) = k
                end do
             else
                rand = ran1(iseed1)
                lcs(1) = 1 + int(npg + 0.5_dp - sqrt( (npg+.5)**2 - npg * (npg+1) * rand ))
                do k = 2, nps
                   lpos_ok = .false.
                   do while (.not. lpos_ok)
                      lpos_ok = .true.
                      rand = ran1(iseed1)
                      lpos = 1 + int(npg + 0.5_dp - sqrt((npg+.5)**2 -  npg * (npg+1) * rand ))
                      ! test if point already chosen: returns lpos_ok=false if any point already exists
                      do k1 = 1, k-1
                         if (lpos .eq. lcs(k1)) lpos_ok = (lpos_ok .and. .false.)
                      end do
                   end do
                   lcs(k) = lpos
                end do
                !
                !  arrange the sub-complex in order of inceasing function value
                call sort(lcs(1:nps))
             end if
             !
             !  create the sub-complex arrays
             do k = 1, nps
                do j = 1, nopt
                   s(k,j) = cx(lcs(k),j)
                end do
                sf(k) = cf(lcs(k))
             end do
             !
             !  use the sub-complex to generate new point(s)
             call cce(s(1:nps,1:nopt),sf(1:nps),bl(1:nopt),bu(1:nopt),xnstd(1:nopt),  &
                      icall,maxn,maxit,iseed1,functn, alpha,beta,history_tmp)
             !
             !  if the sub-complex is accepted, replace the new sub-complex
             !  into the complex
             do k = 1, nps
                do j = 1, nopt
                   cx(lcs(k),j) = s(k,j)
                end do
                cf(lcs(k)) = sf(k)
             end do
             !
             !  sort the points
             call sort_matrix(cx(1:npg,1:nopt),cf(1:npg))
             !
             !  if maximum number of runs exceeded, break out of the loop
             if (icall .ge. maxn) exit
             !
          end do subcomploop ! <alpha loop>
          !
          !  replace the new complex into original array x(.,.)
          do k1 = 1, npg
             k2 = (k1-1) * ngs1 + igs
             do j = 1, nopt
                x(k2,j) = cx(k1,j)
             end do
             xf(k2) = cf(k1)
          end do
          if (icall .ge. maxn) exit
          !
          !  end loop on complexes
       end do comploop  ! <beta loop> 
       !
       !  re-sort the points
       call sort_matrix(x(1:npt1,1:nopt),xf(1:npt1))
       !
       !  record the best and worst points
       do j = 1, nopt
          bestx(j) = x(1,j)
          worstx(j) = x(npt1,j)
       end do
       bestf_tmp = xf(1)
       worstf    = xf(npt1)
       ! if (.not. maxit) then
       !    bestf_tmp = xf(1)
       !    worstf    = xf(npt1)
       ! else
       !    bestf_tmp = -xf(1)
       !    worstf    = -xf(npt1)
       ! end if
       !
       !  test the population for parameter convergence
       call parstt(x(1:npt1,1:nopt),bound,xnstd,gnrng,ipcnvg)
       !
       !  print the results for current population
       if (iprint .lt. 2) then
          if (mod(nloop,5) .eq. 0) then
             write(ipr,610) (xname(j),j=1,nopt1)
             if (nopt .gt. 8) then
                write(ipr,620) (xname(j),j=9,nopt)
             end if
          end if
          if (.not. maxit) then
             write(ipr,630) nloop,icall,ngs1,bestf_tmp,worstf,gnrng, (bestx(j),j=1,nopt1)
          else
             write(ipr,630) nloop,icall,ngs1,-bestf_tmp,-worstf,gnrng, (bestx(j),j=1,nopt1)
          end if
          if (nopt .gt. 8) then 
             write(ipr,640) (bestx(j),j=9,nopt)
          end if
          if (iprint .eq. 1) then
             write(ipr,650) nloop
             do i = 1, npt1
                write(ipr,660) xf(i),(x(i,j),j=1,nopt1)
                if (nopt .gt. 8) then 
                   write(ipr,640) (x(i,j),j=9,nopt)
                end if
             end do
          end if
       end if
       !
       !  test if maximum number of function evaluations exceeded
       ! Maximum number of function evaluations reached?
       if (icall .ge. maxn) then
          if (iprint .lt. 2) then
             ! maximum trials reached
             write(ipr,800) maxn,loop,igs,nloop
             ! parameter set
             write(ipr,830)
             write(ipr,510) (xname(j),j=1,nopt2)
             if (.not. maxit) then
                write(ipr,520) bestf_tmp,(bestx(j),j=1,nopt2)
             else
                write(ipr,520) -bestf_tmp,(bestx(j),j=1,nopt2)
              end if
             if (nopt .gt. 12) then
                write(ipr,530) (xname(j),j=13,nopt)
                write(ipr,540) (bestx(j),j=13,nopt)
             end if
          end if
          if (present(neval)) neval = icall
          if (present(bestf) .and. .not. maxit) bestf = bestf_tmp
          if (present(bestf) .and. maxit)       bestf = -bestf_tmp
          if (present(history)) then
             allocate(history(icall))
             history(:) = history_tmp(1:icall)
          end if
          ! -----------------------
          ! Abort subroutine
          ! -----------------------
          return
       end if
       !
       !  compute the count on successive loops w/o function improvement
       criter(kstop+1) = bestf_tmp
       if (nloop .gt. kstop) then 
          denomi = dabs(criter((kstop+1)-kstop) + criter(kstop+1)) / 2.
          rtiny = 1.0/(1.0**15)
          denomi = max(denomi, rtiny)
          timeou = dabs(criter((kstop+1)-kstop) - criter(kstop+1)) / denomi
          if (timeou .lt. pcento) then 
             if (iprint .lt. 2) then
                ! criterion value has not changed during last loops
                write(ipr,810) pcento*100.,kstop
                ! parameter set
                write(ipr,830)
                write(ipr,510) (xname(j),j=1,nopt2)
                if (.not. maxit) then
                   write(ipr,520) bestf_tmp,(bestx(j),j=1,nopt2)
                else
                   write(ipr,520) -bestf_tmp,(bestx(j),j=1,nopt2)
                end if
                if (nopt .gt. 12) then 
                   write(ipr,530) (xname(j),j=13,nopt)
                   write(ipr,540) (bestx(j),j=13,nopt)
                end if
             end if
             if (present(neval)) neval = icall
             if (present(bestf) .and. .not. maxit) bestf = bestf_tmp
             if (present(bestf) .and. maxit)       bestf = -bestf_tmp
             if (present(history)) then
                allocate(history(icall))
                history(:) = history_tmp(1:icall)
             end if
             ! -----------------------
             ! Abort subroutine
             ! -----------------------
             return
          end if
       end if
       do l = 1, kstop
          criter(l) = criter(l+1)
       end do
       !
       !  if population is converged into a sufficiently small space
       if (ipcnvg .eq. 1) then 
          if (iprint .lt. 2) then
             ! converged because feasible parameter space small
             write(ipr,820) gnrng*100.
             ! parameter set
             write(ipr,830)
             write(ipr,510) (xname(j),j=1,nopt2)
             if (.not. maxit) then
                write(ipr,520) bestf_tmp,(bestx(j),j=1,nopt2)
             else
                write(ipr,520) -bestf_tmp,(bestx(j),j=1,nopt2)
             end if
             if (nopt .gt. 12) then
                write(ipr,530) (xname(j),j=13,nopt)
                write(ipr,540) (bestx(j),j=13,nopt)
             end if
          end if
          if (present(neval)) neval = icall
          if (present(bestf) .and. .not. maxit) bestf = bestf_tmp
          if (present(bestf) .and. maxit)       bestf = -bestf_tmp
          if (present(history)) then
             allocate(history(icall))
             history(:) = history_tmp(1:icall)
          end if
          ! -----------------------
          ! Abort subroutine
          ! -----------------------
          return
       end if
       !
       !  none of the stopping criteria is satisfied, continue search
       !
       !  check for complex number reduction
       if (ngs1 .gt. mings) then
          ngs2 = ngs1
          ngs1 = ngs1 - 1
          npt1 = ngs1 * npg
          call comp(ngs2,npg,x(1:ngs2*npg,1:nopt),xf(1:ngs2*npg),cx(1:ngs2*npg,1:nopt),cf(1:ngs2*npg))
       end if
    end do mainloop
    
    if (present(neval)) neval = icall
    if (present(bestf) .and. .not. maxit) bestf = bestf_tmp
    if (present(bestf) .and. maxit)       bestf = -bestf_tmp
    if (present(history)) then
       allocate(history(icall))
       history(:) = history_tmp(1:icall)
    end if

400 format('==================================================',/, &
         'ENTER THE SHUFFLED COMPLEX EVOLUTION GLOBAL SEARCH',/, &
         '==================================================',/   )
500 format(//,'*** PRINT THE INITIAL POINT AND ITS CRITERION VALUE ***')
510 format(/,' CRITERION',12(6x,a4),/, &
              '------------------------------------------------------------')
520 format(g10.3,12f10.3)
530 format(10x,12(6x,a4))
540 format(10x,12f10.3)
600 format(//,1x,'*** PRINT THE RESULTS OF THE SCE SEARCH ***')
610 format(/,1x,'LOOP',1x,'TRIALS',1x,'COMPLXS',2x,'BEST F',3x,'WORST F',3x,'PAR RNG',1x,8(6x,a4))
620 format(49x,8(6x,a4))
630 format(i5,1x,i5,3x,i5,3g10.3,8(f10.3))
640 format(49x,8(f10.3))
650 format(/,1x,'POPULATION AT LOOP ',i3,/,1x,'---------------------------')
660 format(15x,g10.3,20x,8(f10.3))
800 format(//,1x,'*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE',  &
         ' LIMIT ON THE MAXIMUM',/,5x,'NUMBER OF TRIALS ',i5,     &
         ' EXCEEDED.  SEARCH WAS STOPPED AT',/,5x,'SUB-COMPLEX ', &
         i3,' OF COMPLEX ',i3,' IN SHUFFLING LOOP ',i3,' ***')
810 format(//,1x,'*** OPTIMIZATION TERMINATED BECAUSE THE CRITERION', &
         ' VALUE HAS NOT CHANGED ',/,5x,f5.2,' PERCENT IN',i3,        &
         ' SHUFFLING LOOPS ***')
820 format(//,1x,'*** OPTIMIZATION TERMINATED BECAUSE THE POPULATION', &
         ' HAS CONVERGED INTO ',/,4x,f5.2,' PERCENT OF THE', &
         ' FEASIBLE SPACE ***')
830 format(//,'*** PRINT THE FINAL PARAMETER ESTIMATE AND ITS', &
         ' CRITERION VALUE ***')

  end function sce

  subroutine parstt(x,bound,xnstd,gnrng,ipcnvg)   
    !
    !  subroutine checking for parameter convergence
    !
    use mo_kind, only: i4, dp

    implicit none  

    real(dp),    dimension(:,:),           intent(in)  :: x      ! points in population, cols=nopt, rows=npt1 (!)
    real(dp),    dimension(:),             intent(in)  :: bound  ! difference of upper and lower limit per parameter
    real(dp),    dimension(size(bound,1)), intent(out) :: xnstd  ! std. deviation of points in population per parameter
    real(dp),                              intent(out) :: gnrng  ! fraction of feasible space covered by complexes
    integer(i4),                           intent(out) :: ipcnvg ! 1 : population converged into sufficiently small space
    !                                                            ! 0 : not converged
    !
    ! local variables
    integer(i4)                         :: nopt    ! number of parameters
    integer(i4)                         :: npt1    ! number of points in current population
    integer(i4)                         :: i, k
    real(dp)                            :: xsum1   ! sum           of all values per parameter
    real(dp)                            :: xsum2   ! sum of square of all values per parameter
    real(dp)                            :: gsum    ! sum of all (range scaled) currently covered parameter ranges: 
    !                                              !    (max-min)/bound
    real(dp), dimension(size(bound,1))  :: xmin    ! minimal value per parameter
    real(dp), dimension(size(bound,1))  :: xmax    ! maximal value per parameter
    real(dp), dimension(size(bound,1))  :: xmean   ! mean    value per parameter
    real(dp),                 parameter :: delta = tiny(1.0_dp)
    real(dp),                 parameter :: peps  = 0.001_dp
    !
    nopt = size(bound,1)
    npt1 = size(x,1)

    !  compute maximum, minimum and standard deviation of parameter values
    gsum = 0._dp
    do k = 1, nopt
       xmax(k) = -huge(1.0_dp)
       xmin(k) =  huge(1.0_dp)
       xsum1 = 0._dp
       xsum2 = 0._dp
       do i = 1, npt1
          xmax(k) = dmax1(x(i,k), xmax(k))
          xmin(k) = dmin1(x(i,k), xmin(k))
          xsum1 = xsum1 + x(i,k)
          xsum2 = xsum2 + x(i,k)*x(i,k)
       end do
       xmean(k) = xsum1 / real(npt1,dp)
       xnstd(k) = (xsum2 / real(npt1,dp) - xmean(k)*xmean(k))
       if (xnstd(k) .le. delta) then
          xnstd(k) = delta
       end if
       xnstd(k) = sqrt(xnstd(k))
       xnstd(k) = xnstd(k) / bound(k)
       gsum = gsum + log( delta + (xmax(k)-xmin(k))/bound(k) )
    end do
    gnrng = exp(gsum/real(nopt,dp))
    !
    !  check if normalized standard deviation of parameter is <= eps
    ipcnvg = 0_i4
    if (gnrng .le. peps) then
       ipcnvg = 1_i4
    end if

  end subroutine parstt

  subroutine comp(ngs2,npg,a,af,b,bf)   
    !
    !  This subroutine reduce 
    !  input matrix a(n,ngs2*npg) to matrix b(n,ngs1*npg) and 
    !  input vector af(ngs2*npg)  to vector bf(ngs1*npg)
    !
    !  Needed for complex reduction: 
    !     Number of complexes decreases from ngs2 to ngs1
    !     Therefore, number of points in population decreases (npt+npg) to npt
    !
    use mo_kind, only: i4, dp

    implicit none  

    integer(i4),                              intent(in)    :: ngs2  ! OLD number of complexes = ngs1+1
    integer(i4),                              intent(in)    :: npg   ! number of points in each complex
    real(dp), dimension(:,:),                 intent(inout) :: a     ! OLD points,          ncols=n, nrows=ngs2*npg
    real(dp), dimension(:),                   intent(inout) :: af    ! OLD function values,          nrows=ngs2*npg
    real(dp), dimension(size(a,1),size(a,2)), intent(out)   :: b     ! NEW points,          ncols=n, nrows=ngs1*npg=npt 
    real(dp), dimension(size(af)),            intent(out)   :: bf    ! NEW function values,          nrows=ngs1*npg=npt
    !
    ! local variables
    integer(i4) :: n      ! cols: number of parameters: nopt
    integer(i4) :: npt    ! rows: number of points in NEW population: npt1
    integer(i4) :: ngs1   ! NEW number of complexes
    integer(i4) :: i, j
    integer(i4) :: igs    ! current new complex
    integer(i4) :: ipg    ! current new point in complex igs
    integer(i4) :: k1, k2

    n = size(a,dim=2)
    ! NEW number of complexes
    ngs1 = ngs2 - 1
    ! number of points in NEW population
    npt  = ngs1 * npg

    do igs=1, ngs1
       do ipg=1, npg
          k1=(ipg-1)*ngs2 + igs
          k2=(ipg-1)*ngs1 + igs
          do i=1, n
             b(k2,i) = a(k1,i)
          end do
          bf(k2) = af(k1)
       end do
    end do
    !
    do j=1, npt
       do i=1, n
          a(j,i) = b(j,i)
       end do
       af(j) = bf(j)
    end do
    !
  end subroutine comp

  subroutine sort_matrix(rb,ra)
    !
    ! This subroutine is adapted from "Numerical Recipes" by Press et al., pp. 233-234
    !
    use mo_sort, only: sort_index
    use mo_kind, only: i4, dp

    implicit none    

    real(dp), dimension(:),   intent(inout) :: ra    ! rows: number of points in population --> npt1
    real(dp), dimension(:,:), intent(inout) :: rb    ! rows: number of points in population --> npt1
    !                                                ! cols: number of parameters --> nopt
    ! local variables
    integer(i4)                        :: n          ! number of points in population --> npt1
    integer(i4)                        :: m          ! number of parameters --> nopt
    integer(i4)                        :: i
    integer(i4), dimension(size(rb,1)) :: iwk
    !
    n = size(rb,1)
    m = size(rb,2)

    ! indexes of sorted reference vector
    iwk(:) = sort_index(ra(1:n))

    ! sort reference vector
    ra(1:n) = ra(iwk)

    ! sort each column of second array
    do i=1, m
       rb(1:n,i) = rb(iwk,i)
    end do

  end subroutine sort_matrix

  function ran1(seed)

    !  This subroutine is from "Numerical Recipes" by Press et al.

    use mo_kind, only: i4, i8, dp
    implicit none 

    integer(i8), intent(inout)    :: seed
    real(dp)                      :: ran1

    ! local variables to be saved for next call
    integer(i8), parameter        :: m1 = 259200, m2 = 134456, m3 = 243000
    integer(i8), parameter        :: ia1 = 7141,  ia2 = 8121,  ia3 = 4561
    integer(i8), parameter        :: ic1 = 54773, ic2 = 28411, ic3 = 51349
    real(dp),    parameter        :: rm1 = 3.8580247e-6, rm2 = 7.4373773e-6
    real(dp), dimension(97), save :: r
    integer(i8),             save :: iff = 0       ! for first call 0, afterwards 1
    integer(i8),             save :: ix1, ix2, ix3

    ! local variables
    integer(i4) :: j

    if ((seed .lt. 0) .or. (iff .eq. 0)) then
       iff = 1_i8
       ix1 = mod(ic1 - seed,m1)
       ix1 = mod((ia1 * ix1) + ic1,m1)
       ix2 = mod(ix1,m2)
       ix1 = mod((ia1 * ix1) + ic1,m1)
       ix3 = mod(ix1,m3)
       do j = 1, 97
          ix1 = mod((ia1 * ix1) + ic1,m1)
          ix2 = mod((ia2 * ix2) + ic2,m2)
          r(j) = (real(ix1,dp) + (real(ix2,dp) * rm2)) * rm1
       end do
       seed = 1_i8
    end if
    ix1 = mod((ia1 * ix1) + ic1,m1)
    ix2 = mod((ia2 * ix2) + ic2,m2)
    ix3 = mod((ia3 * ix3) + ic3,m3)
    j = 1 + ((97 * ix3) / m3)
    if ((j .gt. 97) .or. (j .lt. 1)) stop 'Error ran1: j out of range'
    ran1 = r(j)
    r(j) = (real(ix1,dp) + (real(ix2,dp) * rm2)) * rm1

  end function ran1

  function gasdev(seed)

    !  THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.

    use mo_kind, only: i4, i8, dp
    implicit none

    integer(i8), intent(inout) :: seed 
    real(dp)                   :: gasdev

    ! local variables to be saved for next call
    integer(i4) :: iset = 0
    ! local variables
    real(dp) :: v1, v2
    real(dp) :: r
    real(dp) :: fac
    real(dp) :: gset

    iset = 0
    r    = 2.0_dp
    if (iset .eq. 0) then
       do while (r .gt. 1._dp)
          v1 = (2._dp * ran1(seed)) - 1._dp
          v2 = (2._dp * ran1(seed)) - 1._dp
          r = (v1 ** 2) + (v2 ** 2)
       end do
       fac    = sqrt(- ((2._dp * log(r)) / r))
       gset   = v1 * fac
       gasdev = v2 * fac
       iset   = 1
    else
       gasdev = gset
       iset   = 0
    end if

  end function gasdev

  subroutine chkcst(x,bl,bu,ibound)
    !
    !     This subroutine check if the trial point satisfies all
    !     constraints.
    !
    !     ibound - violation indicator
    !            = 0  no violation
    !            = 1  violation
    !     nopt = number of optimizing variables
    !     ii = the ii'th variable of the arrays x, bl, and bu
    !
    !     Note: removed checking for implicit constraints (was anyway empty)
    !
    use mo_kind, only: i4, dp
    implicit none

    real(dp), dimension(:),    intent(in)  :: x      ! trial point dim=(nopt)
    real(dp), dimension(:),    intent(in)  :: bl     ! lower bound dim=(nopt)
    real(dp), dimension(:),    intent(in)  :: bu     ! upper bound dim=(nopt)
    integer(i4),               intent(out) :: ibound ! violation indicator
    ! local variables
    integer(i4) :: ii
    integer(i4) :: nopt
    !
    ibound = 0_i4
    nopt   = size(x,1)

    do ii=1, nopt
       if ( (x(ii) .lt. bl(ii)) .or. (x(ii) .gt. bu(ii)) ) then
          ! This constraint is violated   
          ibound = 1_i4
          return
       end if
    end do

  end subroutine chkcst

  subroutine getpnt(idist,bl,bu,std,xi,iseed,x)
    !
    !     This subroutine generates a new point within feasible region
    !
    !     Note: checking of implicit constraints removed
    !
    use mo_kind, only: i4, i8, dp

    implicit none 

    integer(i4),                     intent(in)    :: idist   ! idist = probability flag
    !                                                         !       = 1 - uniform distribution
    !                                                         !       = 2 - Gaussian distribution
    real(dp), dimension(:),          intent(in)    :: bl      ! lower bound
    real(dp), dimension(:),          intent(in)    :: bu      ! upper bound
    real(dp), dimension(:),          intent(in)    :: std     ! standard deviation of probability distribution
    real(dp), dimension(:),          intent(in)    :: xi      ! focal point
    integer(i8),                     intent(inout) :: iseed   ! seed for random number generator
    real(dp), dimension(size(xi,1)), intent(out)   :: x       ! new point
    ! 
    ! local variables
    integer(i4) :: nopt    ! number of parameters
    integer(i4) :: j
    integer(i4) :: ibound   ! 0=point in bound, 1=point out of bound
    real(dp)    :: rand
    !
    nopt = size(xi,1)
    do j=1, nopt
       ibound = 1
       do while (ibound .eq. 1)
          if (idist .eq. 1) rand = ran1(iseed)
          if (idist .eq. 2) rand = gasdev(iseed)
          x(j) = xi(j) + std(j) * rand * (bu(j) - bl(j))
          !
          !     Check explicit constraints
          !        
          call chkcst((/x(j)/),(/bl(j)/),(/bu(j)/),ibound)
       end do
    end do

  end subroutine getpnt

  subroutine cce(s,sf,bl,bu,xnstd,icall,maxn,maxit,iseed, functn, &
       alpha,beta,history)
    !
    !  ALGORITHM GENERATE A NEW POINT(S) FROM A SUB-COMPLEX
    !
    !  Note: new intent IN    variables for flexible reflection & contraction: alpha, beta
    !        new intent INOUT variable for history of objective function values
    !
    use mo_kind, only: i4, i8, dp

    implicit none    

    real(dp),    dimension(:,:), intent(inout) :: s       ! points in sub-complex, cols=nopt, rows=nps
    real(dp),    dimension(:),   intent(inout) :: sf      ! objective function value of points in 
    !                                                     ! sub-complex, rows=nps
    real(dp),    dimension(:),   intent(in)    :: bl      ! lower bound per parameter
    real(dp),    dimension(:),   intent(in)    :: bu      ! upper bound per parameter
    real(dp),    dimension(:),   intent(in)    :: xnstd   ! standard deviation of points in sub-complex per parameter
    integer(i8),                 intent(inout) :: icall   ! number of function evaluations
    integer(i8),                 intent(in)    :: maxn    ! maximal number of function evaluations allowed
    logical,                     intent(in)    :: maxit   ! minimization (true) or maximization (false)
    integer(i8),                 intent(inout) :: iseed   ! seed for random number generation
    interface
       function functn(paraset)
         ! calculates the cost function at a certain parameter set paraset
         use mo_kind
         real(dp), dimension(:), intent(in)    :: paraset
         real(dp)                              :: functn
       end function functn
    end interface
    real(dp),                    intent(in)    :: alpha   ! parameter for reflection steps
    real(dp),                    intent(in)    :: beta    ! parameter for contraction steps
    real(dp), dimension(maxn),   intent(inout) :: history ! history of best function value
    !
    ! local variables
    integer(i4)                                :: nopt    ! number of parameters
    integer(i4)                                :: nps     ! number of points in a sub-complex
    integer(i4)                                :: i, j
    integer(i4)                                :: n       ! nps  = number of points in sub-complex 
    integer(i4)                                :: m       ! nopt = number of parameters
    integer(i4)                                :: ibound  ! ibound = flag indicating if constraints are violated
    !                                                     !        = 1   yes
    !                                                     !        = 0   no
    real(dp), dimension(size(s,2))             :: sw      ! the worst point of the simplex
    real(dp), dimension(size(s,2))             :: sb      ! the best point of the simplex
    real(dp), dimension(size(s,2))             :: ce      ! the centroid of the simplex excluding wo
    real(dp), dimension(size(s,2))             :: snew    ! new point generated from the simplex
    real(dp)                                   :: fw      ! function value of the worst point
    real(dp)                                   :: fnew    ! function value of the new point

    nps  = size(s,1)
    nopt = size(s,2)
    ! equivalence of variables for readabilty of code
    n = nps
    m = nopt
    !
    ! identify the worst point wo of the sub-complex s
    ! compute the centroid ce of the remaining points
    ! compute step, the vector between wo and ce
    ! identify the worst function value fw
    do j = 1, m
       sb(j) = s(1,j)
       sw(j) = s(n,j)
       ce(j) = 0.0_dp
       do i = 1, n-1
          ce(j) = ce(j) + s(i,j)
       end do
       ce(j) = ce(j)/real(n-1,dp)
    end do
    fw = sf(n)
    !
    ! compute the new point snew
    !
    ! first try a reflection step
    do j = 1, m
       snew(j) = ce(j) + alpha * (ce(j) - sw(j))
    end do
    !
    ! check if snew satisfies all constraints
    call chkcst(snew(1:nopt),bl(1:nopt),bu(1:nopt),ibound)
    !
    ! snew is outside the bound,
    ! choose a point at random within feasible region according to
    ! a normal distribution with best point of the sub-complex
    ! as mean and standard deviation of the population as std
    if (ibound .eq. 1) then
       call getpnt(2,bl(1:nopt),bu(1:nopt),xnstd(1:nopt),sb(1:nopt),iseed,snew)
    end if
    !
    ! compute the function value at snew
    if (.not. maxit) then
       fnew = functn(snew)
       icall = icall + 1_i8
       history(icall) = min(history(icall-1),fnew)
    else
       fnew = -functn(snew)
       icall = icall + 1_i8
       history(icall) = max(history(icall-1),-fnew)
    end if
    !
    ! maximum numbers of function evaluations reached
    if (icall .ge. maxn) return
    !
    ! compare fnew with the worst function value fw
    ! fnew is greater than fw, so try a contraction step
    if (fnew .gt. fw) then
       do j = 1, m
          snew(j) = ce(j) - beta * (ce(j) - sw(j))
       end do
       !
       ! compute the function value of the contracted point
       if (.not. maxit) then
          fnew = functn(snew)
          icall = icall + 1_i8
          history(icall) = min(history(icall-1),fnew)
       else
          fnew = -functn(snew)
          icall = icall + 1_i8
          history(icall) = max(history(icall-1),-fnew)
       end if
       !
       ! maximum numbers of function evaluations reached
       if (icall .ge. maxn) return
       !
       ! compare fnew to the worst value fw
       ! if fnew is less than or equal to fw, then accept the point and return
       if (fnew .gt. fw) then
          !
          ! if both reflection and contraction fail, choose another point
          ! according to a normal distribution with best point of the sub-complex
          ! as mean and standard deviation of the population as std
          call getpnt(2,bl(1:nopt),bu(1:nopt),xnstd(1:nopt),sb(1:nopt),iseed,snew)
          !
          ! compute the function value at the random point
          if (.not. maxit) then
             fnew = functn(snew)
             icall = icall + 1_i8
             history(icall) = min(history(icall-1),fnew)
          else
             fnew = -functn(snew)
             icall = icall + 1_i8
             history(icall) = max(history(icall-1),-fnew)
          end if
          !
          ! maximum numbers of function evaluations reached
          if (icall .ge. maxn) return
          !
          ! successful mutation
          ! replace the worst point by the new point
          do j = 1, m
             s(n,j) = snew(j)
          end do
          sf(n) = fnew
       end if
    end if
    !
    ! replace the worst point by the new point
    do j = 1, m
       s(n,j) = snew(j)
    end do
    sf(n) = fnew

  end subroutine cce


END MODULE mo_sce

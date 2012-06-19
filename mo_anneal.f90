MODULE mo_anneal

  ! This module is minimizing a cost function via Simulated Annealing
  ! and is part of the UFZ CHS Fortran library.

  ! 

  ! Written  Juliane Mai, March 2012
  ! Modified Juliane Mai, May   2012 : anneal: sp version
  !          Juliane Mai, May   2012 : anneal: documentation
  !          Juliane Mai, May   2012 : GetTemperature: sp and dp version
  !          Juliane Mai, June  2012 : weighted parameter selection

  USE mo_kind,    ONLY: i4, i8, sp, dp
  USE mo_xor4096, ONLY: xor4096

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: anneal            ! Minimization of a cost function via Simaulated Annealing
  PUBLIC :: GetTemperature    ! Function which returns the optimal initial temperature for 
                              ! a given acceptance ratio and initial parameter set

  ! Interfaces for single and double precision routines; sort alphabetically
  INTERFACE anneal
     MODULE PROCEDURE anneal_dp, anneal_sp
  END INTERFACE anneal

  INTERFACE GetTemperature
     MODULE PROCEDURE GetTemperature_dp,  GetTemperature_sp
  END INTERFACE GetTemperature


  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         anneal

  !     PURPOSE
  !         Minimizes a user provided cost function using the Simulated Annealing strategy
  !    

  !     CALLING SEQUENCE
  !         para = (/ 1.0_dp , 2.0_dp /)
  !         temperature = 10.0_dp
  !
  !         user defined function cost_dp which calculates the cost function value for a 
  !         parameter set (the interface given below has to be used for this function!)
  !
  !         user defined subroutine range_dp which returns the range for a parameter  
  !         (the subroutine has to have a specific interface givene below!)
  !
  !         call anneal(cost_dp, para, range_dp, temperature, costbest, parabest)
  !
  !         see also test folder for a detailed example

  
  !     INDENT(IN)
  !         REAL(SP/DP),    DIMENSION(:)  :: para  ... initial parameter set
  !
  !         REAL(SP/DP)                   :: temp  ... initial temperature

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         REAL(SP/DP),                              :: costbest ... minimized value of cost 
  !                                                                   function
  !
  !         REAL(SP/DP),    DIMENSION(size(para,1))   :: parabest ... parameter set minimizing the
  !                                                                   cost function

  !     INDENT(IN), OPTIONAL
  !         REAL(SP/DP)    :: DT_in          ... geometrical decreement of temperature
  !                                              0.7<DT<0.999
  !                                              DEFAULT: 0.9
  !
  !         INTEGER(I4)    :: nITERmax_in    ... maximal number of iterations
  !                                              will be increased by 10% if stopping criteria of
  !                                              acceptance ratio or epsilon decreement of cost 
  !                                              function is not fullfilled
  !                                              DEFAULT: 1000
  !
  !         INTEGER(I4)    :: LEN_in         ... Length of Markov Chain
  !                                              DEFAULT: MAX(250, size(para,1))
  !
  !         INTEGER(I4)    :: nST_in         ... Number of consecutive LEN steps
  !                                              DEFAULT: 5
  !
  !         REAL(SP/DP)    :: eps_in         ... Stopping criteria of epsilon decreement of 
  !                                              cost function
  !                                              DEFAULT: 0.01
  !
  !         REAL(SP/DP)    :: acc_in         ... Stopping criteria for Acceptance Ratio
  !                                              acc_in <= 0.1
  !                                              DEFAULT: 0.1
  !
  !         INTEGER(I4/I8) :: seeds_in(3)    ... Seeds of random numbers used for random parameter
  !                                              set generation
  !                                              DEFAULT: dependent on current time
  !
  !         LOGICAL        :: printflag_in   ... If .true. detailed command line output is written
  !                                              DEFAULT: .false.
  !
  !         LOGICAL        :: coststatus_in  ... If .true. the status of cost function value is  
  !                                              checked, i.e. rejects parameter sets which are not
  !                                              valid therefore the calculation of the user 
  !                                              provided FUNCTION cost has to set the OPTIONAL 
  !                                              coststatus variable to .false. in case of an 
  !                                              invalid parameter set
  !                                              DEFAULT = .false. (assuming that each parameter set
  !                                                                 is valid)
  !
  !         LOGICAL, DIMENSION(size(para,1)) :: 
  !                           maskpara_in    ... vector of logicals
  !                                              maskpara(i) = .true.  --> parameter is optimized
  !                                              maskpara(i) = .false. --> parameter is discarded 
  !                                                                        from optimiztaion
  !                                              DEFAULT = .true.
  !         REAL(DP), DIMENSION(size(para,1)) ::
  !                           weight_in      ... vector of weights per parameter
  !                                              gives the frequency of parameter to be chosen for
  !                                              optimization (will be scaled to a CDF internally)
  !                                              eg. [1,2,1] --> parameter 2 is chosen twice as 
  !                                                              often as parameter 1 and 2
  !                                              DEFAULT: weight_in = 1.0_dp

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Needs a user defined:
  !         (1) FUNCTION to evaluate the cost function value for a given parameter set.
  !             this function needs to have the following interface:
  !                     FUNCTION cost(paraset,status_in)
  !                     ! calculates the cost function at a certain parameter set paraset and 
  !                     ! returns optionally if it is a valid parameter set (status_in)
  !                         USE mo_kind
  !                         REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset
  !                         REAL(DP)                            :: cost
  !                         LOGICAL,  OPTIONAL,     INTENT(OUT) :: status_in
  !                     END FUNCTION cost
  !         (2) SUBROUTINE to specify the valid range of each parameter 
  !             (at a certain parameter set).
  !             this routine needs to have the following interface:
  !                     SUBROUTINE range(paraset,iPar,rangePar)
  !                     ! gives the range (min,max) of the parameter iPar at a 
  !                     ! certain parameter set paraset
  !                         USE mo_kind
  !                         REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset
  !                         INTEGER(I4),            INTENT(IN)  :: iPar
  !                         REAL(DP), DIMENSION(2), INTENT(OUT) :: rangePar
  !                     END SUBROUTINE range
  !     EXAMPLE
  !         see test/test_mo_anneal/

  !     LITERATURE
  !         (1) S. Kirkpatrick, C. D. Gelatt, and M. P. Vecchi. 
  !             Optimization by simulated annealing. 
  !             Science, 220:671–680, 1983.
  !
  !         (2) N. Metropolis, A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller, and E. Teller. 
  !             Equation of state calculations by fast computing machines. 
  !             J. Chem. Phys., 21:1087–1092, June 1953.

  !     HISTORY
  !        Written  Juliane Mai, March 2012
  !        Modified Juliane Mai, May   2012 : sp version
  !                 Juliane Mai, May   2012 : documentation

  SUBROUTINE anneal_dp(cost, para, range, temp, costbest, parabest, & 
                       DT_in, nITERmax_in, LEN_in, nST_in, & 
                       eps_in, acc_in, seeds_in, printflag_in, coststatus_in, maskpara_in, weight_in)

    IMPLICIT NONE

    INTERFACE
       FUNCTION cost(paraset,status_in)
       ! calculates the cost function at a certain parameter set paraset and 
       ! returns optionally if it is a valid parameter set (status_in)
           use mo_kind
           REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset
           REAL(DP)                            :: cost
           LOGICAL,  OPTIONAL,     INTENT(OUT) :: status_in
       END FUNCTION cost
    END INTERFACE
    
    INTERFACE
       SUBROUTINE range(paraset,iPar,rangePar)
       ! gives the range (min,max) of the parameter iPar at a certain parameter set paraset
           use mo_kind
           REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset
           INTEGER(I4),            INTENT(IN)  :: iPar
           REAL(DP), DIMENSION(2), INTENT(OUT) :: rangePar
       END SUBROUTINE range
    END INTERFACE

    REAL(DP),    DIMENSION(:),            INTENT(IN)    :: para  ! initial parameter i
    REAL(DP),                             INTENT(IN)    :: temp  ! starting temperature
    REAL(DP),                             INTENT(OUT)   :: costbest  ! minimized value of cost function
    REAL(DP),    DIMENSION(size(para,1)), INTENT(OUT)   :: parabest  ! parameter set minimizing the cost function
    REAL(DP),    OPTIONAL,       INTENT(IN)  :: DT_in    ! geometrical decreement, 0.7<DT<0.999
    INTEGER(I4), OPTIONAL,       INTENT(IN)  :: nITERmax_in ! maximal number of iterations
    INTEGER(I4), OPTIONAL,       INTENT(IN)  :: LEN_in   ! Length of Markov Chain, MAX(250, size(para,1))
    INTEGER(I4), OPTIONAL,       INTENT(IN)  :: nST_in   ! Number of consecutive LEN steps
    REAL(DP),    OPTIONAL,       INTENT(IN)  :: eps_in   ! epsilon decreement of cost function
    REAL(DP),    OPTIONAL,       INTENT(IN)  :: acc_in   ! Acceptance Ratio, <0.1 stopping criteria
    INTEGER(I8), OPTIONAL,       INTENT(IN)  :: seeds_in(3)   ! Seeds of random numbers
    LOGICAL,     OPTIONAL,       INTENT(IN)  :: printflag_in ! If command line output is written (.true.) 
    LOGICAL,     OPTIONAL,       INTENT(IN)  :: coststatus_in  ! Checks status of cost function value, 
                                                               ! i.e. is parameter set is feasible (.true.)
                                                               ! default = .false.
    LOGICAL,     OPTIONAL, DIMENSION(size(para,1)), INTENT(IN)  :: maskpara_in  
                                                         ! true if parameter will be optimized
                                                         ! false if parameter is discarded in optimization
                                                         ! default = .true.
    REAL(DP),    OPTIONAL, DIMENSION(size(para,1)), INTENT(IN)  :: weight_in  
                                                         ! vector of weights per parameter
                                                         ! gives the frequency of parameter to be chosen for optimization 

    INTEGER(I4) :: n    ! Number of parameters
    REAL(DP)    :: T
    REAL(DP)    :: DT  
    INTEGER(I4) :: nITERmax  
    INTEGER(I4) :: LEN   
    INTEGER(I4) :: nST   
    REAL(DP)    :: eps   
    REAL(DP)    :: acc 
    INTEGER(I8) :: seeds(3)
    LOGICAL     :: printflag
    LOGICAL     :: coststatus
    LOGICAL,     DIMENSION(size(para,1))   :: maskpara    ! true if parameter will be optimized
    INTEGER(I4), DIMENSION(:), ALLOCATABLE :: truepara    ! indexes of parameters to be optimized
    REAL(DP),    DIMENSION(size(para,1))   :: weight      ! CDF of parameter to chose for optimization
    
    type paramLim
       real(DP)                                 :: min                 !            minimum value
       real(DP)                                 :: max                 !            maximum value
       real(DP)                                 :: new                 !            new state value
       real(DP)                                 :: old                 !            old state value
       real(DP)                                 :: best                !            best value found
       real(DP)                                 :: dMult               !            sencitivity multiplier for parameter search
    end type paramLim
    type (paramLim), dimension (size(para,1)), target   :: gamma      !            Transfer function parameters

    ! for random numbers
    real(DP)                     :: RN1, RN2, RN3     ! Random numbers
    integer(I8)                  :: iin1, iin2, iin3  ! optional arguments for restarting RN streams
    integer(I8)                  :: win1, win2, win3  ! optional arguments for restarting RN streams
    integer(I8), dimension(0:63) :: xin1, xin2,xin3   ! optional arguments for restarting RN streams

    ! for SA
    integer(I4)            :: idummy,i
    character(10)          :: timeseed
    real(DP)               :: NormPhi
    real(DP)               :: ac_ratio, pa
    real(DP)               :: fo, fn, df, fBest, rho, fInc, fbb, dr
    real(DP)               :: T0, DT0
    real(DP),  parameter   :: small = -700._DP
    !integer(I4)            :: nITER
    integer(I4)            :: j, iter, kk
    integer(I4)            :: Ipos, Ineg
    integer(I4)            :: iConL, iConR, iConF
    integer(I4)            :: iTotalCounter                        ! includes reheating for final conditions
    integer(I4)            :: iTotalCounterR                       ! counter of interations in one reheating
    logical                :: iStop 
    
    integer(I4)            :: iPar
    real(DP), DIMENSION(2) :: iParRange

    n = size(para,1)

    if (present(DT_in)) then
       if ( (DT_in .lt. 0.7_dp) .or. (DT_in .gt. 0.999_dp) ) then
          stop 'Input argument DT must lie between 0.7 and 0.999'
       else
          DT = DT_in
       end if
    else
       DT = 0.9_dp
    endif

    if (present(nITERmax_in)) then
       if (nITERmax_in .lt. 1_I4) then
          stop 'Input argument nITERmax must be greater 0'
       else
          nITERmax = nITERmax_in
       end if
    else
       nITERmax = 1000_i4
    endif

    if (present(LEN_in)) then
       if (LEN_in .lt. Max(20_i4*n,250_i4)) then
          !stop 'Input argument LEN must be greater than Max(250,20*N), N=number of parameters'
          print*, 'WARNING: Input argument LEN should be greater than Max(250,20*N), N=number of parameters'
          LEN = LEN_in
       else
          LEN = LEN_in
       end if
    else
       LEN = Max(20_i4*n,250_i4)
    endif

    if (present(nST_in)) then
       if (nST_in .lt. 0_i4) then
          stop 'Input argument nST must be greater than 0'
       else
          nST = nST_in
       end if
    else
       nST = 5_i4
    endif

    if (present(eps_in)) then
       if (eps_in .le. 0.0_dp) then
          stop 'Input argument eps must be greater than 0'
       else
          eps = eps_in
       end if
    else
       eps = 0.01_dp
    endif

    if (present(acc_in)) then
       if ((acc_in .le. 0.0_dp) .or. (acc_in .ge. 1.0_dp)) then
          stop 'Input argument acc must lie between 0.0 and 1.0'
       else
          acc = acc_in
       end if
    else
       acc = 0.1_dp
    endif

    if (present(seeds_in)) then
       seeds = seeds_in
    else
       ! Seeds depend on actual time
       call date_and_time(time=timeseed)
       read(timeseed,'(i6,1x,i3)') idummy, seeds(1)
       print*,'seeds(1)=', seeds(1)
       seeds(2) = seeds(1)+1000_i8
       seeds(3) = seeds(2)+1000_i8
    endif

    if (present(printflag_in)) then
       printflag = printflag_in
    else
       printflag = .false.
    endif

    if (present(maskpara_in)) then
       if (count(maskpara_in) .eq. 0_i4) then
          stop 'Input argument maskpara: At least one element has to be true'
       else
          maskpara = maskpara_in
       end if
    else
       maskpara = .true.
    endif
    
    allocate ( truepara(count(maskpara)) )
    idummy = 0_i4
    do i=1,N
       if ( maskpara(i) ) then
          idummy = idummy+1_i4
          truepara(idummy) = i
       end if
    end do
    
    if (printflag) then
       print*, 'Following parameters will be optimized: ',truepara
    end if

    weight    = 0.0_dp
    if (present(weight_in)) then
       where ( maskpara(:) ) 
          weight(:) = weight_in(:)
       end where
    else
       where ( maskpara(:) ) 
          weight(:) = 1.0_dp
       end where
    endif
    ! scaling the weights
    weight = weight/sum(weight)
    ! cummulating the weights
    do i=2,n
       weight(i) = weight(i) + weight(i-1)
    end do

    T = temp

    call xor4096(seeds(1), RN1, iin=iin1, win=win1, xin=xin1)
    call xor4096(seeds(2), RN2, iin=iin2, win=win2, xin=xin2)
    call xor4096(seeds(3), RN3, iin=iin3, win=win3, xin=xin3)
    seeds = 0_i8

    ! Start Simulated Annealing routine
    !gamma(:)%min = para(:,2)
    !gamma(:)%max = para(:,3)
    gamma(:)%dmult = 1.0_dp
    gamma(:)%new = para(:)
    gamma(:)%old = para(:)
    gamma(:)%best = para(:)
    NormPhi = -9999.9_dp
    T0=      T
    DT0=     DT
    !nITER =  2_i4*N
    !
    ! Generate and evaluate the initial solution state
    fo = cost(gamma(:)%old)   
    !
    ! initialize counters /var for new SA    
    iConL=           0_i4
    iConR=           1_i4
    iConF=           0_i4
    iTotalCounter=   0_i4
    iTotalCounterR = 0_i4   
    fn =             0.0_DP
    ac_ratio =       1.0_DP
    iStop=           .TRUE.
    iter =           0_i4
    ! restoring initial T param.
    !T=      T0
    DT=     DT0
    ! Storing the best solution so far
    NormPhi = fo
    fo      = fo/NormPhi
    fBest   = 1.0_dp   != fo/NormPhi

    if (printflag) then
       print '(A15,E15.7,A4,E15.7,A4)',     ' start NSe   = ', fBest , '  ( ',fBest*normPhi,' )  '
       print '(I8, 2I5, 4E15.7)',   1_i4, 0_i4, 0_i4, 1._dp, T, fo*normPhi, fBest*normPhi
    end if
    !
    ! ****************** Stop Criterium  *******************
    ! Repeat until the % reduction of nST consecutive Markov
    ! chains (LEN) of the objective function (f) <= epsilon
    ! ******************************************************
    loopTest: do while (iStop)
    ! loopTest: do iter = 1, nIter
       iter = iter + 1_i4    
       Ipos= 0_i4
       Ineg= 0_i4
       fbb=  fBest
       !LEN = int( real(iter,dp)*sqrt(real(iter,dp)) / nITER + 1.5*real(N,dp),i4) 
       ! Repeat LEN times with feasible solution
       j=1
       loopLEN: do while (j .le. LEN)
         !print*, 'iPar: ',j
         iTotalCounterR =  iTotalCounterR + 1_i4
         iTotalCounter = iTotalCounter + 1_i4
         !if (mod(iTotalCounter,1000_i4) == 0_i4) then
         !  call writeresults_o(4,fn) 
         !end if
         !
	 ! Generate a random subsequent state and evaluate its objective function
         ! (1) Generate new parameter set
         dR=(1.0_DP - real(iTotalCounterR,dp) / real(nITERmax,dp))**2.0_DP
         if (  .not. (dR >= 0.05_DP ) .and. &
               (iTotalCounterR <= int(real(nIterMax,dp)/3._dp*4_dp,i4))) then
            dR = 0.05_DP
         end if   
         ! (1a) Select parameter to be changed    
         call xor4096(seeds(2),RN2, iin=iin2, win=win2, xin=xin2)
         ! V1: select one parameter from all n
         !iPar = int(real(N,dp)*RN2+1._DP,i4)
         !
         ! V2: select only masked parameters
         !iPar = truepara(int(real(count(maskpara),dp)*RN2+1._DP,i4))
         !
         ! V3: select according to CDF based on weights and mask
         iPar=1_i4
         do while (weight(iPar) .lt. RN2)
            iPar = iPar + 1_i4
         end do
         !print*, 'iPar=',iPar !
         ! (1b) Generate new value of selected parameter
         call xor4096(seeds(3),RN3, iin=iin3, win=win3, xin=xin3)
         call range(gamma(:)%old, iPar, iParRange )
         gamma(iPar)%min = iParRange(1)
         gamma(iPar)%max = iParRange(2)
         gamma(iPar)%new = parGen_dp( gamma(iPar)%old, gamma(iPar)%dMult*dR, &
                                      gamma(iPar)%min, gamma(iPar)%max, 8_i4, 1_i4,RN3)
         ! (2) Calculate new objective function value and normalize it
         !print*, 'coststatus before=',coststatus
         if (present(coststatus_in)) then
            if (coststatus_in) then
               fn = cost(gamma(:)%new,coststatus)  !  coststatus is INTENT(OUT)
            else
               fn = cost(gamma(:)%new) 
               coststatus = .true.
            end if
         else
            fn = cost(gamma(:)%new) 
            coststatus = .true.
         end if
         !print*, 'coststatus after=',coststatus
          
         feasible: if (coststatus) then   ! feasible parameter set
          fn = fn/normPhi
          !
          ! Change in cost function value
          df = fn-fo             
          !
          ! analyze change in the objective function: df
          if (df < 0.0_DP) then
             !
             ! accept the new state
             Ipos=Ipos+1_i4
             !print*, 'Ipos=',Ipos
             fo = fn
             gamma(:)%old   = gamma(:)%new 
             !
             ! keep best solution
             if (fo < fBest) then
               fBest =  fo
               gamma(:)%best   = gamma(:)%new
             endif
          else
             if ( df >  eps ) then
               rho=-df/T
               if (rho < small) then   !small = -700._dp
                 pa=0.0_DP
               else
                 pa=EXP(rho)
               end if    
               !print*, 'rho=',rho,'  pa=',pa                 
  	           !
               call xor4096(seeds(1), RN1, iin=iin1, win=win1, xin=xin1)
               !print*, RN1
               !
               if (pa > RN1) then
                 ! accept new state with certain probability
                 Ineg=Ineg+1_i4
                 !print*, 'Ineg=',Ineg
                 !pause
                 fo = fn
                 ! save old state
                 gamma(:)%old   = gamma(:)%new 
               end if
             !else
                !print*,'NO CHANGE IN OBJ FUNCTION',gamma(:)%new
             end if
          end if
          j=j+1
         else
           iTotalCounterR =  iTotalCounterR - 1_i4
           iTotalCounter = iTotalCounter - 1_i4
         end if feasible !valid parameter set
      end do loopLEN
      !
      ! estimate acceptance ratio
      ac_ratio=(real(Ipos,dp) + real(Ineg,dp))/real(LEN,dp)
      !    
      if (printflag) then
         print '(I8, 2I5, E15.7, 3E15.7)', ITotalCounter, Ipos, Ineg, ac_ratio, T, &
                                           fo*NormPhi, fBest*NormPhi
      end if
      !
      ! Cooling schedule
      fInc= (fbb-fBest)/fbb
      if (fInc < 0.00000001_dp) then
        iConF= iConF+1_i4
      else
        iConF=0_i4
      end if
      !
      if ((ac_ratio < 0.15_dp) .and. (iConF > 5_i4) .and. (iConR <= -3_i4)) then     ! - iConR  no reheating
        ! Re-heating
        if (printflag) then
           print *, 'Re-heating: ', iConR
        end if
        iConR= iConR+1_i4
        T=T0/2._DP
        iter = 0_i4            ! for LEN
        iTotalCounterR = 0_i4  ! for dR
        ! start from current best
        gamma(:)%old   = gamma(:)%best
      else
        ! Update Temperature (geometrical decrement)
        if (ac_ratio < 0.4_dp)  then 
          DT=0.995_dp
        else
          DT=DT0
        end if
        T=T*DT
      end if
      !
      ! Stop Criteria: consecutive MC with marginal decrements and acceptance ratio
      if (fInc < eps) then
        iConL= iConL+1_i4
      else
        iConL=0_i4
      end if
      if ((iConL > nST) .and. (ac_ratio < acc) .and. (iConR > 2_i4) ) iStop=.FALSE.
      ! a way out maximum
      if ( iTotalCounter > nITERmax .and. ac_ratio > acc) then
         nITERmax = int(real(nITERmax,dp)* 1.10_dp,i4)
         if (printflag) then
            print *,                       'nITERmax changed to =', nITERmax
         end if
      end if
      ! STOP condition
      if ( iTotalCounter > nITERmax ) then
         !print*, 'iStop = .false.'
         iStop=.FALSE.          
      end if                          
      !
    end do loopTest
    !
    ! heuristic *glopbal* minimum (end solution)
    !
    ! write final results
    !
    ! calculate cost function again (only for check and return values)  
    parabest = gamma(:)%best
    costbest = cost(parabest)
    fo = costbest/NormPhi

    if (printflag) then
       print *, '   '
       !print '(A15,E15.7,A4,E15.7,A4)',       ' end cost    = ', fBest , '  ( ',fBest*normPhi,' )  '
       print '(A15,E15.7)',                   ' end cost    = ',fBest*normPhi
       print *,           'end parameter: '
       do kk = 1,N
          print '(A10,I3,A3, E15.7)' ,    '    para #',kk,' = ',gamma(kk)%best
       end do
    
       print *, 'Final check:    ', (fo - fBest)
    end if
    
    deallocate(truepara)

END SUBROUTINE anneal_dp

SUBROUTINE anneal_sp(cost, para, range, temp, costbest, parabest, & 
                       DT_in, nITERmax_in, LEN_in, nST_in, & 
                       eps_in, acc_in, seeds_in, printflag_in, coststatus_in, maskpara_in, weight_in)

    IMPLICIT NONE

    INTERFACE
       FUNCTION cost(paraset,status_in)
       ! calculates the cost function at a certain parameter set paraset and 
       ! returns optionally if it is a valid parameter set (status_in)
           use mo_kind
           REAL(SP), DIMENSION(:), INTENT(IN)  :: paraset
           REAL(SP)                            :: cost
           LOGICAL,  OPTIONAL,     INTENT(OUT) :: status_in
       END FUNCTION cost
    END INTERFACE
    
    INTERFACE
       SUBROUTINE range(paraset,iPar,rangePar)
       ! gives the range (min,max) of the parameter iPar at a certain parameter set paraset
           use mo_kind
           REAL(SP), DIMENSION(:), INTENT(IN)  :: paraset
           INTEGER(I4),            INTENT(IN)  :: iPar
           REAL(SP), DIMENSION(2), INTENT(OUT) :: rangePar
       END SUBROUTINE range
    END INTERFACE

    REAL(SP),    DIMENSION(:),            INTENT(IN)    :: para  ! initial parameter i
    REAL(SP),                             INTENT(IN)    :: temp  ! starting temperature
    REAL(SP),                             INTENT(OUT)   :: costbest  ! minimized value of cost function
    REAL(SP),    DIMENSION(size(para,1)), INTENT(OUT)   :: parabest  ! parameter set minimizing the cost function
    REAL(SP),    OPTIONAL,       INTENT(IN)  :: DT_in    ! geometrical decreement, 0.7<DT<0.999
    INTEGER(I4), OPTIONAL,       INTENT(IN)  :: nITERmax_in ! maximal number of iterations
    INTEGER(I4), OPTIONAL,       INTENT(IN)  :: LEN_in   ! Length of Markov Chain, MAX(250, size(para,1))
    INTEGER(I4), OPTIONAL,       INTENT(IN)  :: nST_in   ! Number of consecutive LEN steps
    REAL(SP),    OPTIONAL,       INTENT(IN)  :: eps_in   ! epsilon decreement of cost function
    REAL(SP),    OPTIONAL,       INTENT(IN)  :: acc_in   ! Acceptance Ratio, <0.1 stopping criteria
    INTEGER(I4), OPTIONAL,       INTENT(IN)  :: seeds_in(3)   ! Seeds of random numbers
    LOGICAL,     OPTIONAL,       INTENT(IN)  :: printflag_in ! If command line output is written (.true.) 
    LOGICAL,     OPTIONAL,       INTENT(IN)  :: coststatus_in  ! Checks status of cost function value, 
                                                               ! i.e. is parameter set is feasible (.true.)
                                                               ! default = .false.
    LOGICAL,     OPTIONAL, DIMENSION(size(para,1)), INTENT(IN)  :: maskpara_in  
                                                         ! true if parameter will be optimized
                                                         ! false if parameter is discarded in optimization
                                                         ! default = .true.
    REAL(SP),    OPTIONAL, DIMENSION(size(para,1)), INTENT(IN)  :: weight_in  
                                                         ! vector of weights per parameter
                                                         ! gives the frequency of parameter to be chosen for optimization 


    INTEGER(I4) :: n    ! Number of parameters
    REAL(SP)    :: T
    REAL(SP)    :: DT  
    INTEGER(I4) :: nITERmax  
    INTEGER(I4) :: LEN   
    INTEGER(I4) :: nST   
    REAL(SP)    :: eps   
    REAL(SP)    :: acc 
    INTEGER(I4) :: seeds(3)
    LOGICAL     :: printflag
    LOGICAL     :: coststatus
    LOGICAL,     DIMENSION(size(para,1))   :: maskpara    ! true if parameter will be optimized
    INTEGER(I4), DIMENSION(:), ALLOCATABLE :: truepara    ! indexes of parameters to be optimized
    REAL(SP),    DIMENSION(size(para,1))   :: weight      ! CDF of parameter to chose for optimization
    
    type paramLim
       real(SP)                                 :: min                 !            minimum value
       real(SP)                                 :: max                 !            maximum value
       real(SP)                                 :: new                 !            new state value
       real(SP)                                 :: old                 !            old state value
       real(SP)                                 :: best                !            best value found
       real(SP)                                 :: dMult               !            sencitivity multiplier for parameter search
    end type paramLim
    type (paramLim), dimension (size(para,1)), target   :: gamma      !            Transfer function parameters

    ! for random numbers
    real(SP)                      :: RN1, RN2, RN3     ! Random numbers
    integer(I4)                   :: iin1, iin2, iin3  ! optional arguments for restarting RN streams
    integer(I4)                   :: win1, win2, win3  ! optional arguments for restarting RN streams
    integer(I4), dimension(0:127) :: xin1, xin2,xin3   ! optional arguments for restarting RN streams

    ! for SA
    integer(I4)            :: idummy,i
    character(10)          :: timeseed
    real(SP)               :: NormPhi
    real(SP)               :: ac_ratio, pa
    real(SP)               :: fo, fn, df, fBest, rho, fInc, fbb, dr
    real(SP)               :: T0, DT0
    real(SP),  parameter   :: small = -700._DP
    !integer(I4)            :: nITER
    integer(I4)            :: j, iter, kk
    integer(I4)            :: Ipos, Ineg
    integer(I4)            :: iConL, iConR, iConF
    integer(I4)            :: iTotalCounter                        ! includes reheating for final conditions
    integer(I4)            :: iTotalCounterR                       ! counter of interations in one reheating
    logical                :: iStop 
    
    integer(I4)            :: iPar
    real(SP), DIMENSION(2) :: iParRange

    n = size(para,1)

    if (present(DT_in)) then
       if ( (DT_in .lt. 0.7_sp) .or. (DT_in .gt. 0.999_sp) ) then
          stop 'Input argument DT must lie between 0.7 and 0.999'
       else
          DT = DT_in
       end if
    else
       DT = 0.9_sp
    endif

    if (present(nITERmax_in)) then
       if (nITERmax_in .lt. 1_I4) then
          stop 'Input argument nITERmax must be greater 0'
       else
          nITERmax = nITERmax_in
       end if
    else
       nITERmax = 1000_i4
    endif

    if (present(LEN_in)) then
       if (LEN_in .lt. Max(20_i4*n,250_i4)) then
          !stop 'Input argument LEN must be greater than Max(250,20*N), N=number of parameters'
          print*, 'WARNING: Input argument LEN should be greater than Max(250,20*N), N=number of parameters'
          LEN = LEN_in
       else
          LEN = LEN_in
       end if
    else
       LEN = Max(20_i4*n,250_i4)
    endif

    if (present(nST_in)) then
       if (nST_in .lt. 0_i4) then
          stop 'Input argument nST must be greater than 0'
       else
          nST = nST_in
       end if
    else
       nST = 5_i4
    endif

    if (present(eps_in)) then
       if (eps_in .le. 0.0_sp) then
          stop 'Input argument eps must be greater than 0'
       else
          eps = eps_in
       end if
    else
       eps = 0.01_sp
    endif

    if (present(acc_in)) then
       if ((acc_in .le. 0.0_sp) .or. (acc_in .ge. 1.0_sp)) then
          stop 'Input argument acc must lie between 0.0 and 1.0'
       else
          acc = acc_in
       end if
    else
       acc = 0.1_sp
    endif

    if (present(seeds_in)) then
       seeds = seeds_in
    else
       ! Seeds depend on actual time
       call date_and_time(time=timeseed)
       read(timeseed,'(i6,1x,i3)') idummy, seeds(1)
       seeds(2) = seeds(1)+1000_i4
       seeds(3) = seeds(2)+1000_i4
    endif

    if (present(printflag_in)) then
       printflag = printflag_in
    else
       printflag = .false.
    endif

    if (present(maskpara_in)) then
       if (count(maskpara_in) .eq. 0_i4) then
          stop 'Input argument maskpara: At least one element has to be true'
       else
          maskpara = maskpara_in
       end if
    else
       maskpara = .true.
    endif
    
    allocate ( truepara(count(maskpara)) )
    idummy = 0_i4
    do i=1,N
       if ( maskpara(i) ) then
          idummy = idummy+1_i4
          truepara(idummy) = i
       end if
    end do
    
    if (printflag) then
       print*, 'Following parameters will be optimized: ',truepara
    end if

    weight    = 0.0_sp
    if (present(weight_in)) then
       where ( maskpara(:) ) 
          weight(:) = weight_in(:)
       end where
    else
       where ( maskpara(:) ) 
          weight(:) = 1.0_sp
       end where
    endif
    ! scaling the weights
    weight = weight/sum(weight)
    ! cummulating the weights
    do i=2,n
       weight(i) = weight(i) + weight(i-1)
    end do

    T = temp

    call xor4096(seeds(1), RN1, iin=iin1, win=win1, xin=xin1)
    call xor4096(seeds(2), RN2, iin=iin2, win=win2, xin=xin2)
    call xor4096(seeds(3), RN3, iin=iin3, win=win3, xin=xin3)
    seeds = 0_i4

    ! Start Simulated Annealing routine
    !gamma(:)%min = para(:,2)
    !gamma(:)%max = para(:,3)
    gamma(:)%dmult = 1.0_sp
    gamma(:)%new = para(:)
    gamma(:)%old = para(:)
    gamma(:)%best = para(:)
    NormPhi = -9999.9_sp
    T0=      T
    DT0=     DT
    !nITER =  2_i4*N
    !
    ! Generate and evaluate the initial solution state
    fo = cost(gamma(:)%old)   
    !
    ! initialize counters /var for new SA    
    iConL=           0_i4
    iConR=           1_i4
    iConF=           0_i4
    iTotalCounter=   0_i4
    iTotalCounterR = 0_i4   
    fn =             0.0_SP
    ac_ratio =       1.0_SP
    iStop=           .TRUE.
    iter =           0_i4
    ! restoring initial T param.
    !T=      T0
    DT=     DT0
    ! Storing the best solution so far
    NormPhi = fo
    fo      = fo/NormPhi
    fBest   = 1.0_sp   != fo/NormPhi

    if (printflag) then
       print '(A15,E15.7,A4,E15.7,A4)',     ' start NSe   = ', fBest , '  ( ',fBest*normPhi,' )  '
       print '(I8, 2I5, 4E15.7)',   1_i4, 0_i4, 0_i4, 1._sp, T, fo*normPhi, fBest*normPhi
    end if
    !
    ! ****************** Stop Criterium  *******************
    ! Repeat until the % reduction of nST consecutive Markov
    ! chains (LEN) of the objective function (f) <= epsilon
    ! ******************************************************
    loopTest: do while (iStop)
    ! loopTest: do iter = 1, nIter
       iter = iter + 1_i4    
       Ipos= 0_i4
       Ineg= 0_i4
       fbb=  fBest
       !LEN = int( real(iter,sp)*sqrt(real(iter,sp)) / nITER + 1.5*real(N,sp),i4) 
       ! Repeat LEN times with feasible solution
       j=1_i4
       loopLEN: do while (j .le. LEN)
         !print*, 'iPar: ',j
         iTotalCounterR =  iTotalCounterR + 1_i4
         iTotalCounter = iTotalCounter + 1_i4
         !if (mod(iTotalCounter,1000_i4) == 0_i4) then
         !  call writeresults_o(4,fn) 
         !end if
         !
	 ! Generate a random subsequent state and evaluate its objective function
         ! (1) Generate new parameter set
         dR=(1.0_SP - real(iTotalCounterR,sp) / real(nITERmax,sp))**2.0_SP
         if (  .not. (dR >= 0.05_SP ) .and. &
               (iTotalCounterR <= int(real(nIterMax,sp)/3._sp*4_sp,i4))) then
            dR = 0.05_SP
         end if   
         ! (1a) Select parameter to be changed    
         call xor4096(seeds(2),RN2, iin=iin2, win=win2, xin=xin2)
         ! V1: select one parameter from all n
         !iPar = int(real(N,sp)*RN2+1._SP,i4)
         !
         ! V2: select only masked parameters
         !iPar = truepara(int(real(count(maskpara),sp)*RN2+1._SP,i4))
         !
         ! V3: select according to CDF based on weights and mask
         iPar=1_i4
         do while (weight(iPar) .lt. RN2)
            iPar = iPar + 1_i4
         end do
         !print*, 'iPar=',iPar !
         ! (1b) Generate new value of selected parameter
         call xor4096(seeds(3),RN3, iin=iin3, win=win3, xin=xin3)
         call range(gamma(:)%old, iPar, iParRange )
         gamma(iPar)%min = iParRange(1)
         gamma(iPar)%max = iParRange(2)
         gamma(iPar)%new = parGen_sp( gamma(iPar)%old, gamma(iPar)%dMult*dR, &
                                      gamma(iPar)%min, gamma(iPar)%max, 8_i4, 1_i4,RN3)
         ! (2) Calculate new objective function value and normalize it
         !print*, 'coststatus before=',coststatus
         if (present(coststatus_in)) then
            if (coststatus_in) then
               fn = cost(gamma(:)%new,coststatus)  !  coststatus is INTENT(OUT)
            else
               fn = cost(gamma(:)%new) 
               coststatus = .true.
            end if
         else
            fn = cost(gamma(:)%new) 
            coststatus = .true.
         end if
         !print*, 'coststatus after=',coststatus
          
         feasible: if (coststatus) then   ! feasible parameter set
          fn = fn/normPhi
          !
          !
          df = fn-fo             
          !
          ! analyze change in the objective function: df
          if (df < 0.0_SP) then
             !
             ! accept the new state
             Ipos=Ipos+1_i4
             !print*, 'Ipos=',Ipos
             fo = fn
             gamma(:)%old   = gamma(:)%new 
             !
             ! keep best solution
             if (fo < fBest) then
               fBest =  fo
               gamma(:)%best   = gamma(:)%new
             endif
          else
             if ( df >  eps ) then
               rho=-df/T
               if (rho < small) then   !small = -700._dp
                 pa=0.0_SP
               else
                 pa=EXP(rho)
               end if    
               !print*, 'rho=',rho,'  pa=',pa                 
  	           !
               call xor4096(seeds(1), RN1, iin=iin1, win=win1, xin=xin1)
               !print*, RN1
               !
               if (pa > RN1) then
                 ! accept new state with certain probability
                 Ineg=Ineg+1_i4
                 !print*, 'Ineg=',Ineg
                 !pause
                 fo = fn
                 ! save old state
                 gamma(:)%old   = gamma(:)%new 
               end if
             !else
                !print*,'NO CHANGE IN OBJ FUNCTION',gamma(:)%new
             end if
          end if
           j=j+1
         else
           iTotalCounterR =  iTotalCounterR - 1_i4
           iTotalCounter = iTotalCounter - 1_i4
         end if feasible !valid parameter set
      end do loopLEN
      !
      ! estimate acceptance ratio
      ac_ratio=(real(Ipos,sp) + real(Ineg,sp))/real(LEN,sp)
      !    
      if (printflag) then
         print '(I8, 2I5, E15.7, 3E15.7)', ITotalCounter, Ipos, Ineg, ac_ratio, T, &
                                           fo*NormPhi, fBest*NormPhi
      end if
      !
      ! Cooling schedule
      fInc= (fbb-fBest)/fbb
      if (fInc < 0.00000001_sp) then
        iConF= iConF+1_i4
      else
        iConF=0_i4
      end if
      !
      if ((ac_ratio < 0.15_sp) .and. (iConF > 5_i4) .and. (iConR <= -3_i4)) then     ! - iConR  no reheating
        ! Re-heating
        if (printflag) then
           print *, 'Re-heating: ', iConR
        end if
        iConR= iConR+1_i4
        T=T0/2._SP
        iter = 0_i4            ! for LEN
        iTotalCounterR = 0_i4  ! for dR
        ! start from current best
        gamma(:)%old   = gamma(:)%best
      else
        ! Update Temperature (geometrical decrement)
        if (ac_ratio < 0.4_sp)  then 
          DT=0.995_sp
        else
          DT=DT0
        end if
        T=T*DT
      end if
      !
      ! Stop Criteria: consecutive MC with marginal decrements and acceptance ratio
      if (fInc < eps) then
        iConL= iConL+1_i4
      else
        iConL=0_i4
      end if
      if ((iConL > nST) .and. (ac_ratio < acc) .and. (iConR > 2_i4) ) iStop=.FALSE.
      ! a way out maximum
      if ( iTotalCounter > nITERmax .and. ac_ratio > acc) then
         nITERmax = int(real(nITERmax,dp)* 1.10_sp,i4)
         if (printflag) then
            print *,                       'nITERmax changed to =', nITERmax
         end if
      end if
      ! STOP condition
      if ( iTotalCounter > nITERmax ) then
         !print*, 'iStop = .false.'
         iStop=.FALSE.          
      end if                          
      !
    end do loopTest
    !
    ! heuristic *glopbal* minimum (end solution)
    !
    ! write final results
    !
    ! calculate cost function again (only for check and return values)  
    parabest = gamma(:)%best
    costbest = cost(parabest)
    fo = costbest/NormPhi

    if (printflag) then
       print *, '   '
       !print '(A15,E15.7,A4,E15.7,A4)',       ' end cost    = ', fBest , '  ( ',fBest*normPhi,' )  '
       print '(A15,E15.7)',                   ' end cost    = ',fBest*normPhi
       print *,           'end parameter: '
       do kk = 1,N
          print '(A10,I3,A3, E15.7)' ,    '    para #',kk,' = ',gamma(kk)%best
       end do
    
       print *, 'Final check:    ', (fo - fBest)
    end if
    
    deallocate(truepara)

END SUBROUTINE anneal_sp

real(DP) function GetTemperature_dp( paraset, cost, range, acc_goal, &
                                     samplesize_in, maskpara_in,seeds_in,printflag_in, weight_in)
  use mo_kind, only: dp, i4, i8
  implicit none
  !
  real(dp),               dimension(:),               intent(in)    :: paraset         ! a valid parameter set of the model
  real(dp),                                           intent(in)    :: acc_goal        ! acceptance ratio to achieve
  integer(i4), optional,                              intent(in)    :: samplesize_in   ! size of random set the 
                                                                                       ! acc_estimate is based on
                                                                                       ! default = Max(250, 20*<Number of parameters>)
  logical,     optional,  dimension(size(paraset,1)), intent(in)    :: maskpara_in     ! true if parameter will be optimized
                                                                                       ! false if parameter is discarded in 
                                                                                       ! optimization
                                                                                       ! default = .true.
  integer(i8), optional,  dimension(2),               intent(in)    :: seeds_in        ! Seeds of random numbers
                                                                                       ! default = time dependent
  logical,     optional,                              intent(in)    :: printflag_in    ! .true. if detailed temperature estimation 
                                                                                       ! is printed
                                                                                       ! default = .false.
  REAL(DP),    OPTIONAL,  DIMENSION(size(paraset,1)), INTENT(IN)     :: weight_in  
                                                         ! vector of weights per parameter
                                                         ! gives the frequency of parameter to be chosen for optimization 

  INTERFACE
     FUNCTION cost(paraset,status_in)
       ! calculates the cost function at a certain parameter set paraset and 
       ! returns optionally if it is a valid parameter set (status_in)
       use mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset
       REAL(DP)                            :: cost
       LOGICAL,  OPTIONAL,     INTENT(OUT) :: status_in
     END FUNCTION cost
  END INTERFACE
    
  INTERFACE
     SUBROUTINE range(paraset,iPar,rangePar)
       ! gives the range (min,max) of the parameter iPar at a certain parameter set paraset
       use mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset
       INTEGER(I4),            INTENT(IN)  :: iPar
       REAL(DP), DIMENSION(2), INTENT(OUT) :: rangePar
     END SUBROUTINE range
  END INTERFACE
  !
  ! Local variables
  !
  integer(i4)                           :: n
  integer(i4)                           :: samplesize
  integer(i4)                           :: idummy, i, j
  integer(i4)                           :: iPar
  real(dp),      dimension(2)           :: iParRange
  character(10)                         :: timeseed
  real(dp)                              :: NormPhi
  real(dp)                              :: fo, fn, dr
  real(dp)                              :: T

  LOGICAL                                   :: coststatus
  LOGICAL,     DIMENSION(size(paraset,1))   :: maskpara    ! true if parameter will be optimized
  INTEGER(I4), DIMENSION(:), ALLOCATABLE    :: truepara    ! indexes of parameters to be optimized
  LOGICAL                                   :: printflag   ! if detailed estimation of temperature is printed
  REAL(DP),    DIMENSION(size(paraset,1))   :: weight      ! CDF of parameter to chose for optimization
  
  type paramLim
     real(DP)                                 :: min                 !            minimum value
     real(DP)                                 :: max                 !            maximum value
     real(DP)                                 :: new                 !            new state value
     real(DP)                                 :: old                 !            old state value
     real(DP)                                 :: best                !            best value found
     real(DP)                                 :: dMult               !            sencitivity multiplier for parameter search
  end type paramLim
  type (paramLim), dimension (size(paraset,1)), target   :: gamma    !            Transfer function parameters

  ! for random numbers
  INTEGER(I8)                  :: seeds(2)
  real(DP)                     :: RN1, RN2     ! Random numbers
  integer(I8)                  :: iin1, iin2   ! optional arguments for restarting RN streams
  integer(I8)                  :: win1, win2   ! optional arguments for restarting RN streams
  integer(I8), dimension(0:63) :: xin1, xin2   ! optional arguments for restarting RN streams
  
  ! for initial temperature estimate
  real(DP)                              :: acc_estim  ! estimate of acceptance probability
                                                      ! depends on temperature
                                                      ! goal: find a temperature such that acc_estim ~ 0.9
  real(DP), dimension(:,:), allocatable :: Energy     ! dim(LEN,2) cost function values before (:,1) and after (:,2) transition
 
  n = size(paraset,1)

  if (present(samplesize_in)) then
     if (samplesize_in .lt. Max(20_i4*n,250_i4)) then
        !stop 'Input argument LEN must be greater than Max(250,20*N), N=number of parameters'
        print*, 'WARNING (GetTemperature): '
        print*, 'Input argument samplesize_in should be greater than Max(250,20*N), N=number of parameters'
        samplesize = samplesize_in
     else
        samplesize = samplesize_in
     end if
  else
     samplesize = Max(20_i4*n,250_i4)
  endif

  if (present(maskpara_in)) then
     if (count(maskpara_in) .eq. 0_i4) then
        stop 'Input argument maskpara: At least one element has to be true'
     else
        maskpara = maskpara_in
     end if
  else
     maskpara = .true.
  endif

  if (present(seeds_in)) then
     seeds = seeds_in
  else
     ! Seeds depend on actual time
     call date_and_time(time=timeseed)
     read(timeseed,'(i6,1x,i3)') idummy, seeds(1)
     print*,'temp: seeds(1)=', seeds(1)
     seeds(2) = seeds(1)+1000_i8
  endif

  if (present(printflag_in)) then
     printflag = printflag_in
  else
     printflag = .false.
  endif

  allocate(Energy(samplesize,2))

  allocate ( truepara(count(maskpara)) )
  idummy = 0_i4
  do i=1,N
     if ( maskpara(i) ) then
        idummy = idummy+1_i4
        truepara(idummy) = i
     end if
  end do

  weight    = 0.0_dp
  if (present(weight_in)) then
     where ( maskpara(:) ) 
        weight(:) = weight_in(:)
     end where
  else
     where ( maskpara(:) ) 
        weight(:) = 1.0_dp
     end where
  endif
  ! scaling the weights
  weight = weight/sum(weight)
  ! cummulating the weights
  do i=2,n
     weight(i) = weight(i) + weight(i-1)
  end do

  ! Setting up the RNG
  ! (1) Seeds depend on actual time or on input seeds
  ! (2) Initialize the streams
        call xor4096(seeds(1), RN1, iin=iin1, win=win1, xin=xin1)
        call xor4096(seeds(2), RN2, iin=iin2, win=win2, xin=xin2)
        seeds = 0_i8
  ! (3) Now ready for calling


!TODO
! status? 

  gamma(:)%dmult = 1.0_dp
  gamma(:)%new   = paraset(:)
  gamma(:)%old   = paraset(:)
  gamma(:)%best  = paraset(:)
  NormPhi        = -9999.9_dp

  fo =  cost(paraset) 
  NormPhi = fo
  fo      = fo/NormPhi

  j=1
  loopSamplesize: do while (j .le. samplesize)
     ! Generate a random subsequent state and evaluate its objective function
     ! (1)  Generate new parameter set
            dR=1.0_DP
     ! (1a) Select parameter to be changed    
            call xor4096(seeds(1),RN1, iin=iin1, win=win1, xin=xin1)
            ! V1: select one parameter from all n
            !iPar = int(real(N,dp)*RN1+1._DP,i4)
            !
            ! V2: select only masked parameters
            !iPar = truepara(int(real(count(maskpara),dp)*RN1+1._DP,i4))
            !
            ! V3: select according to CDF based on weights and mask
            iPar=1_i4
            do while (weight(iPar) .lt. RN1)
               iPar = iPar + 1_i4
            end do
     ! (1b) Generate new value of selected parameter
            call xor4096(seeds(2),RN2, iin=iin2, win=win2, xin=xin2)
            call range(gamma(:)%old, iPar, iParRange )
            gamma(iPar)%min = iParRange(1)
            gamma(iPar)%max = iParRange(2)
            gamma(iPar)%new = parGen_dp( gamma(iPar)%old, gamma(iPar)%dMult*dR, &
                 gamma(iPar)%min, gamma(iPar)%max, 8_i4, 1_i4,RN2)
     ! (2)  Calculate new objective function value and normalize it
            fn = cost(gamma(:)%new) 
            coststatus = .true.
     
     feasible: if (coststatus) then   ! feasible parameter set
        fn = fn/normPhi
        ! Save Energy states of feasible transitions
        ! for adaption of optimal initial temperature
        ! Walid Ben-Ameur: "Comput. the Initial Temperature of Sim. Annealing"
        ! Comput. Opt. and App. 2004
        Energy(j,2) = fn     ! E_max_t
        Energy(j,1) = fo     ! E_min_t
        j=j+1
      end if feasible !valid parameter set
   end do loopSamplesize

   ! estimation of the acceptance probability based on the random set ||<Samplesize>||
   ! only if actual temperature (T) equals initial temperature (temp)
   T = maxval(Energy)  !1.0_dp
   !print*, Energy(:,2)/T
   acc_estim = sum(exp(-(Energy(:,2)/T))) / sum(exp(-(Energy(:,1)/T)))
   if (printflag) then
      print*, "acc_estimate = ", acc_estim, "    ( T = ",T," )"
   end if
   Do While ( (acc_estim .lt. 1.0_dp) .and. (abs(acc_estim - acc_goal) .gt. 0.0001_dp))
      T = T * (Log(acc_estim)/Log(acc_goal))**(0.5_dp) ! **(1.0/p)  with p=1.0
      acc_estim = sum(exp(-(Energy(:,2)/T))) / sum(exp(-(Energy(:,1)/T)))
      if (printflag) then
         print*, "acc_estimate = ", acc_estim, "    ( T = ",T," )"
      end if
   end do
   GetTemperature_dp = T
   !
end function GetTemperature_dp

real(sp) function GetTemperature_sp( paraset, cost, range, acc_goal, & 
                                     samplesize_in, maskpara_in,seeds_in,printflag_in, weight_in)
  use mo_kind, only: sp, i4
  implicit none
  !
  real(sp),               dimension(:),               intent(in)    :: paraset         ! a valid parameter set of the model
  real(sp),                                           intent(in)    :: acc_goal        ! acceptance ratio to achieve
  integer(i4), optional,                              intent(in)    :: samplesize_in   ! size of random set the 
                                                                                       ! acc_estimate is based on
                                                                                       ! default = Max(250, 20*<Number of parameters>)
  logical,     optional,  dimension(size(paraset,1)), intent(in)    :: maskpara_in     ! true if parameter will be optimized
                                                                                       ! false if parameter is discarded in 
                                                                                       ! optimization
                                                                                       ! default = .true.
  integer(i4), optional,  dimension(2),               intent(in)    :: seeds_in        ! Seeds of random numbers
                                                                                       ! default = time dependent
  logical,     optional,                              intent(in)    :: printflag_in    ! .true. if detailed temperature estimation 
                                                                                       ! is printed
                                                                                       ! default = .false.
  REAL(sp),    OPTIONAL,  DIMENSION(size(paraset,1)), INTENT(IN)     :: weight_in  
                                                         ! vector of weights per parameter
                                                         ! gives the frequency of parameter to be chosen for optimization 

  INTERFACE
     FUNCTION cost(paraset,status_in)
       ! calculates the cost function at a certain parameter set paraset and 
       ! returns optionally if it is a valid parameter set (status_in)
       use mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN)  :: paraset
       REAL(SP)                            :: cost
       LOGICAL,  OPTIONAL,     INTENT(OUT) :: status_in
     END FUNCTION cost
  END INTERFACE
    
  INTERFACE
     SUBROUTINE range(paraset,iPar,rangePar)
       ! gives the range (min,max) of the parameter iPar at a certain parameter set paraset
       use mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN)  :: paraset
       INTEGER(I4),            INTENT(IN)  :: iPar
       REAL(SP), DIMENSION(2), INTENT(OUT) :: rangePar
     END SUBROUTINE range
  END INTERFACE
  !
  ! Local variables
  !
  integer(i4)                           :: n
  integer(i4)                           :: samplesize
  integer(i4)                           :: idummy, i, j
  integer(i4)                           :: iPar
  real(sp),      dimension(2)           :: iParRange
  character(10)                         :: timeseed
  real(sp)                              :: NormPhi
  real(sp)                              :: fo, fn, dr
  real(sp)                              :: T

  LOGICAL                                   :: coststatus
  LOGICAL,     DIMENSION(size(paraset,1))   :: maskpara    ! true if parameter will be optimized
  INTEGER(I4), DIMENSION(:), ALLOCATABLE    :: truepara    ! indexes of parameters to be optimized
  LOGICAL                                   :: printflag   ! if detailed estimation of temperature is printed
  REAL(SP),    DIMENSION(size(paraset,1))   :: weight      ! CDF of parameter to chose for optimization
  
  type paramLim
     real(SP)                                 :: min                 !            minimum value
     real(SP)                                 :: max                 !            maximum value
     real(SP)                                 :: new                 !            new state value
     real(SP)                                 :: old                 !            old state value
     real(SP)                                 :: best                !            best value found
     real(SP)                                 :: dMult               !            sencitivity multiplier for parameter search
  end type paramLim
  type (paramLim), dimension (size(paraset,1)), target   :: gamma    !            Transfer function parameters

  ! for random numbers
  INTEGER(I4)                   :: seeds(2)
  real(SP)                      :: RN1, RN2     ! Random numbers
  integer(I4)                   :: iin1, iin2   ! optional arguments for restarting RN streams
  integer(I4)                   :: win1, win2   ! optional arguments for restarting RN streams
  integer(I4), dimension(0:127) :: xin1, xin2   ! optional arguments for restarting RN streams
  
  ! for initial temperature estimate
  real(SP)                              :: acc_estim  ! estimate of acceptance probability
                                                      ! depends on temperature
                                                      ! goal: find a temperature such that acc_estim ~ 0.9
  real(SP), dimension(:,:), allocatable :: Energy     ! dim(LEN,2) cost function values before (:,1) and after (:,2) transition
 
  n = size(paraset,1)

  if (present(samplesize_in)) then
     if (samplesize_in .lt. Max(20_i4*n,250_i4)) then
        !stop 'Input argument LEN must be greater than Max(250,20*N), N=number of parameters'
        print*, 'WARNING (GetTemperature): '
        print*, 'Input argument samplesize_in should be greater than Max(250,20*N), N=number of parameters'
        samplesize = samplesize_in
     else
        samplesize = samplesize_in
     end if
  else
     samplesize = Max(20_i4*n,250_i4)
  endif

  if (present(maskpara_in)) then
     if (count(maskpara_in) .eq. 0_i4) then
        stop 'Input argument maskpara: At least one element has to be true'
     else
        maskpara = maskpara_in
     end if
  else
     maskpara = .true.
  endif

  if (present(seeds_in)) then
     seeds = seeds_in
  else
     ! Seeds depend on actual time
     call date_and_time(time=timeseed)
     read(timeseed,'(i6,1x,i3)') idummy, seeds(1)
     print*,'temp: seeds(1)=', seeds(1)
     seeds(2) = seeds(1)+1000_i4
  endif

  if (present(printflag_in)) then
     printflag = printflag_in
  else
     printflag = .false.
  endif

  allocate(Energy(samplesize,2))

  allocate ( truepara(count(maskpara)) )
  idummy = 0_i4
  do i=1,N
     if ( maskpara(i) ) then
        idummy = idummy+1_i4
        truepara(idummy) = i
     end if
  end do

  weight    = 0.0_sp
  if (present(weight_in)) then
     where ( maskpara(:) ) 
        weight(:) = weight_in(:)
     end where
  else
     where ( maskpara(:) ) 
        weight(:) = 1.0_sp
     end where
  endif
  ! scaling the weights
  weight = weight/sum(weight)
  ! cummulating the weights
  do i=2,n
     weight(i) = weight(i) + weight(i-1)
  end do
  
  ! Setting up the RNG
  ! (1) Seeds depend on actual time or on input seeds
  ! (2) Initialize the streams
        call xor4096(seeds(1), RN1, iin=iin1, win=win1, xin=xin1)
        call xor4096(seeds(2), RN2, iin=iin2, win=win2, xin=xin2)
        seeds = 0_i4
  ! (3) Now ready for calling


!TODO
! status? 

  gamma(:)%dmult = 1.0_sp
  gamma(:)%new   = paraset(:)
  gamma(:)%old   = paraset(:)
  gamma(:)%best  = paraset(:)
  NormPhi        = -9999.9_sp

  fo      =  cost(paraset) 
  NormPhi = fo
  fo      = fo/NormPhi

  j=1
  loopSamplesize: do while (j .le. samplesize)
     ! Generate a random subsequent state and evaluate its objective function
     ! (1)  Generate new parameter set
            dR=1.0_SP
     ! (1a) Select parameter to be changed    
            call xor4096(seeds(1),RN1, iin=iin1, win=win1, xin=xin1)
            ! V1: select one parameter from all n
            !iPar = int(real(N,sp)*RN1+1._SP,i4)
            !
            ! V2: select only masked parameters
            !iPar = truepara(int(real(count(maskpara),sp)*RN1+1._SP,i4))
            !
            ! V3: select according to CDF based on weights and mask
            iPar=1_i4
            do while (weight(iPar) .lt. RN1)
               iPar = iPar + 1_i4
            end do
     ! (1b) Generate new value of selected parameter
            call xor4096(seeds(2),RN2, iin=iin2, win=win2, xin=xin2)
            call range(gamma(:)%old, iPar, iParRange )
            gamma(iPar)%min = iParRange(1)
            gamma(iPar)%max = iParRange(2)
            gamma(iPar)%new = parGen_sp( gamma(iPar)%old, gamma(iPar)%dMult*dR, &
                 gamma(iPar)%min, gamma(iPar)%max, 8_i4, 1_i4,RN2)
     ! (2)  Calculate new objective function value and normalize it
            fn = cost(gamma(:)%new) 
            coststatus = .true.
     
     feasible: if (coststatus) then   ! feasible parameter set
        fn = fn/normPhi
        ! Save Energy states of feasible transitions
        ! for adaption of optimal initial temperature
        ! Walid Ben-Ameur: "Comput. the Initial Temperature of Sim. Annealing"
        ! Comput. Opt. and App. 2004
        Energy(j,2) = fn     ! E_max_t
        Energy(j,1) = fo     ! E_min_t
        j=j+1
      end if feasible !valid parameter set
   end do loopSamplesize

   ! estimation of the acceptance probability based on the random set ||<Samplesize>||
   ! only if actual temperature (T) equals initial temperature (temp)
   T = maxval(Energy)  !1.0_sp
   acc_estim = sum(exp(-(Energy(:,2)/T))) / sum(exp(-(Energy(:,1)/T)))
   if (printflag) then
      print*, "acc_estimate = ", acc_estim, "    ( T = ",T," )"
   end if
    Do While ( (acc_estim .lt. 1.0_sp) .and. (abs(acc_estim - acc_goal) .gt. 0.0001_sp))
      T = T * (Log(acc_estim)/Log(acc_goal))**(0.5_sp) ! **(1.0/p)  with p=1.0
      acc_estim = sum(exp(-(Energy(:,2)/T))) / sum(exp(-(Energy(:,1)/T)))
      if (printflag) then
         print*, "acc_estimate = ", acc_estim, "    ( T = ",T," )"
      end if
   end do
   GetTemperature_sp = T
   !
end function GetTemperature_sp

!***************************************************
!*               PRIVATE FUNCTIONS                 *
!***************************************************

!
real(DP) function parGen_dp(old, dMax, oMin, oMax,iDigit,isZero,RN)
  use mo_kind, only: dp,i4
  implicit none
  !
  real(DP),    intent(IN)   :: old, dMax, oMin,oMax
  integer(I4), intent(IN)   :: iDigit,iszero
  real(DP),    intent(IN)   :: RN
  real(DP)                  :: oi, ox, delta, dMaxScal
  !
  oi=oMin
  ox=oMax
  ! scaling (new)
  dMaxScal = dMax * ABS(oMax - oMin)
  
  if(oi >= ox) then
    parGen_dp = (oi+ox)/2.0_DP
  else
    if ( (old - dMaxScal) > oi)  oi = old - dMaxScal
    if ( (old + dMaxScal) < ox)  ox = old + dMaxScal
    delta = (oi + RN * (ox - oi) ) - old
    if ( delta > dMaxScal) delta =  dMaxScal
    if ( delta <-dMaxScal) delta = -dMaxScal
    parGen_dp = old + dChange_dp(delta,iDigit,isZero)
  endif
  !
  ! Parameter is bounded between Max and Min.
  ! Correction from Kumar and Luis
  !
  if (parGen_dp < oMin) then
      parGen_dp = oMin
  elseif(parGen_dp > oMax)then
      parGen_dp = oMax
  end if
end function parGen_dp

real(SP) function parGen_sp(old, dMax, oMin, oMax,iDigit,isZero,RN)
  use mo_kind, only: sp,i4
  implicit none
  !
  real(SP),    intent(IN)   :: old, dMax, oMin,oMax
  integer(I4), intent(IN)   :: iDigit,iszero
  real(SP),    intent(IN)   :: RN
  real(SP)                  :: oi, ox, delta, dMaxScal
  !
  oi=oMin
  ox=oMax
  ! scaling (new)
  dMaxScal = dMax * ABS(oMax - oMin)
  
  if(oi >= ox) then
    parGen_sp = (oi+ox)/2.0_SP
  else
    if ( (old - dMaxScal) > oi)  oi = old - dMaxScal
    if ( (old + dMaxScal) < ox)  ox = old + dMaxScal
    delta = (oi + RN * (ox - oi) ) - old
    if ( delta > dMaxScal) delta =  dMaxScal
    if ( delta <-dMaxScal) delta = -dMaxScal
    parGen_sp = old + dChange_sp(delta,iDigit,isZero)
  endif
  !
  ! Parameter is bounded between Max and Min.
  ! Correction from Kumar and Luis
  !
  if (parGen_sp < oMin) then
      parGen_sp = oMin
  elseif(parGen_sp > oMax)then
      parGen_sp = oMax
  end if
end function parGen_sp

! ************************************************************************
! delta change (steps & accuracy)
! ************************************************************************
real(DP) function  dChange_dp(delta,iDigit,isZero)
  use mo_kind, only: i4,i8, dp
  implicit none

  integer(I4), intent(IN)   :: iDigit,isZero
  real(DP),    intent(IN)   :: delta
  integer(I4)               :: ioszt
  integer(I8)               :: iDelta
  !
  ioszt = 10**iDigit
  iDelta = int( delta * real(ioszt,dp),i8 )
  if(isZero == 1_i4) then
    if(iDelta == 0_i8) then
      if(delta < 0_i4) then
        iDelta = -1_i8
      else
        iDelta = 1_i8
      endif
    endif
  endif
  dChange_dp=real(iDelta,dp)/real(ioszt,dp)
end function  dChange_dp

! ************************************************************************
! delta change (steps & accuracy)
! ************************************************************************
real(SP) function  dChange_sp(delta,iDigit,isZero)
  use mo_kind, only: i4,i8, sp
  implicit none

  integer(I4), intent(IN)   :: iDigit,isZero
  real(SP),    intent(IN)   :: delta
  integer(I4)               :: ioszt
  integer(I8)               :: iDelta
  !
  ioszt = 10**iDigit
  iDelta = int( delta * real(ioszt,sp),i8 )
  if(isZero == 1_i4) then
    if(iDelta == 0_i8) then
      if(delta < 0_i4) then
        iDelta = -1_i8
      else
        iDelta = 1_i8
      endif
    endif
  endif
  dChange_sp=real(iDelta,sp)/real(ioszt,sp)
end function  dChange_sp

  ! ------------------------------------------------------------------

END MODULE mo_anneal

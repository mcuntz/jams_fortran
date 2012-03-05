MODULE mo_anneal

  ! This module is minimizing a cost function via Simulated Annealing
  ! and is part of the UFZ CHS Fortran library.

  ! 

  ! Written  Juliane Mai, March 2012
  ! Modified 

  USE mo_kind,    ONLY: i4, i8, sp, dp
  USE mo_xor4096, ONLY: xor4096

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: anneal            ! Minimization of a cost function via Simaulated Annealing

  ! Interfaces for single and double precision routines; sort alphabetically
  INTERFACE anneal
     MODULE PROCEDURE anneal_dp
  END INTERFACE anneal

  !INTERFACE parGen
  !   MODULE PROCEDURE parGen_dp
  !END INTERFACE parGen

  !INTERFACE dChange
  !   MODULE PROCEDURE dChange_dp
  !END INTERFACE dChange

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         anneal

  !     PURPOSE
  !         Minimizes a user provided cost function using the Simulated Annealing strategy
  !
  !         

  !     CALLING SEQUENCE
  !         
  
  !     INDENT(IN)
  !         None

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         

  !     EXAMPLE
  !         

  !     LITERATURE
  !         

  !     HISTORY
  !         Written,  Juliane Mai, March 2012

  SUBROUTINE anneal_dp(cost, para, temp, costbest, parabest, & 
                       DT_in, nITERmax_in, LEN_in, nST_in, & 
                       eps_in, acc_in, seeds_in, printflag_in)

    IMPLICIT NONE

    INTERFACE
       FUNCTION cost(paraset)
           use mo_kind
           REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset
           REAL(DP)                            :: cost
       END FUNCTION cost
    END INTERFACE

    REAL(DP),    DIMENSION(:,:),          INTENT(IN)    :: para  ! parameter(i), min(parameter(i)), max(parameter(i))
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
    integer(I8)            :: idummy
    character(10)          :: timeseed
    real(DP)               :: NormPhi
    real(DP)               :: ac_ratio, pa
    real(DP)               :: fo, fn, df, fBest, rho, fInc, fbb, dr
    real(DP)               :: T0, DT0
    real(DP),  parameter   :: small = -700._DP
    integer(I4)            :: iPar
    integer(I4)            :: nITER
    integer(I4)            :: j, iter, kk
    integer(I4)            :: Ipos, Ineg
    integer(I4)            :: iConL, iConR, iConF
    integer(I4)            :: iTotalCounter                        ! includes reheating for final conditions
    integer(I4)            :: iTotalCounterR                       ! counter of interations in one reheating
    logical(1)             :: iStop 

    n = size(para,1)

    ! Input check
    if (size(para,2) .ne. 3) then
       stop 'Input argument para must have 3 columns (initial parameter value, min value, max value)'
    end if

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
          stop 'Input argument LEN must be greater than Max(250,20*N), N=number of parameters'
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
       acc = 0.1_i4
    endif

    if (present(seeds_in)) then
       seeds = seeds_in
    else
       ! Seeds depend on actual time
       call date_and_time(time=timeseed)
       read(timeseed,'(i6,1x,i3)') idummy, seeds(1)
       seeds(2) = seeds(1)+1000_i8
       seeds(3) = seeds(2)+1000_i8
    endif

    if (present(printflag_in)) then
       printflag = printflag_in
    else
       printflag = .false.
    endif

    T = temp

    call xor4096(seeds(1), RN1, iin=iin1, win=win1, xin=xin1)
    call xor4096(seeds(2), RN2, iin=iin2, win=win2, xin=xin2)
    call xor4096(seeds(3), RN3, iin=iin3, win=win3, xin=xin3)
    seeds = 0_i8

    ! Start Simulated Annealing routine
    gamma(:)%min = para(:,2)
    gamma(:)%max = para(:,3)
    gamma(:)%dmult = 1.0_dp
    gamma(:)%new = para(:,1)
    gamma(:)%old = para(:,1)
    NormPhi = -9999.9_dp
    T0=      T
    DT0=     DT
    nITER =  2_i4*N
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
    fBest   = 1.0_dp   != fo/NormPhi

    if (printflag) then
       print '(A15,E15.7,A4,E15.7,A4)',     ' start NSe   = ', fBest , '  ( ',fBest*normPhi,' )  '
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
       ! Repeat LEN times
       loopLEN: do j=1,LEN
         iTotalCounterR =  iTotalCounterR + 1_i4
         iTotalCounter = iTotalCounter + 1_i4
         if ( (iTotalCounter .eq. 1_i4) .and. printflag)  then
           print '(I8, 2I5, E15.7,3E15.7)',   ITotalCounter, Ipos, Ineg, ac_ratio, T, fo, fBest
         end if
         if (mod(iTotalCounter,1000_i4) == 0_i4) then
           !call writeresults_o(4,fn) 
         end if
         !
	 ! Generate a random subsequent state and evaluate its objective function
         ! (1) Generate new parameter set
         dR=(1.0_DP - real(iTotalCounterR,dp) / real(nITERmax,dp))**2.0_DP
         if (.not. (dR >= 0.05_DP ) .and. (iTotalCounterR <= int(real(nIterMax,dp)/3._dp*4_dp))) dR = 0.05_DP   
         ! (1a) Select parameter to be changed    
         call xor4096(seeds(2),RN2, iin=iin2, win=win2, xin=xin2)
         iPar = int(real(N,dp)*RN2+1._DP,i4)
         ! (1b) Generate new value of selected parameter
         call xor4096(seeds(3),RN3, iin=iin3, win=win3, xin=xin3)
         gamma(iPar)%new = parGen_dp(gamma(iPar)%old,  gamma(iPar)%dMult*dR, gamma(iPar)%min, gamma(iPar)%max, 8_i4, 1_i4,RN3)
         ! (2) Calculate new objective function value and normalize it
         fn = cost(gamma(:)%new) 
         fn = fn/normPhi
         !
         ! write parameters
         !call writeresults_o(5,fn)
         !
         df = fn-fo             
         !
         ! analyze change in the objective function: df
         if (df < 0.0_DP) then
           !
           ! accept the new state
	   Ipos=Ipos+1_i4
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
             if (rho < small) then
               pa=0.0_DP
             else
               pa=EXP(rho)
             end if                     
  	     !
             call xor4096(seeds(1), RN1, iin=iin1, win=win1, xin=xin1)
             !
             if (pa > RN1) then
               ! accept new state with certain probability
               Ineg=Ineg+1_i4
               fo = fn
               ! save old state
               gamma(:)%old   = gamma(:)%new 
             end if
           end if
	 end if

      end do loopLEN
      !
      ! estimate acceptance ratio
      ac_ratio=(real(Ipos,dp) + real(Ineg,dp))/real(LEN,dp)
      !    
      if (printflag) then
         print '(I8, 2I5, E15.7, 3E15.7)', ITotalCounter, Ipos, Ineg, ac_ratio, T, fo, fBest
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
      if ( iTotalCounter > nITERmax ) iStop=.FALSE.                                    
      !
    end do loopTest
    !
    ! heuristic *global* minimum (end solution)
    !
    ! write final results
    if (printflag) then
       print *, '   '
       print '(A15,E15.7,A4,E15.7,A4)',       ' end NSe     = ', fBest , '  ( ',fBest*normPhi,' )  '
       print *,           'end parameter: '
       do kk = 1,N
          print '(A10,I3,A3, E15.7)' ,    '    para #',kk,' = ',gamma(kk)%best
       end do
    
       print *, 'Final check:    ', (fo - fBest)
    end if

    parabest = gamma(:)%best
    costbest = cost(parabest)

  END SUBROUTINE anneal_dp

SUBROUTINE anneal_sp(cost, para, temp, costbest, parabest, & 
                       DT_in, nITERmax_in, LEN_in, nST_in, & 
                       eps_in, acc_in, seeds_in, printflag_in)

    IMPLICIT NONE

    INTERFACE
       FUNCTION cost(paraset)
           use mo_kind
           REAL(SP), DIMENSION(:), INTENT(IN)  :: paraset
           REAL(SP)                            :: cost
       END FUNCTION cost
    END INTERFACE

    REAL(SP),    DIMENSION(:,:),          INTENT(IN)    :: para  ! parameter(i), min(parameter(i)), max(parameter(i))
    REAL(SP),                             INTENT(IN)    :: temp  ! starting temperature
    REAL(SP),                             INTENT(OUT)   :: costbest  ! minimized value of cost function
    REAL(SP),    DIMENSION(size(para,1)), INTENT(OUT)   :: parabest  ! parameter set minimizing the cost function
    REAL(SP),    OPTIONAL,       INTENT(IN)  :: DT_in    ! geometrical decreement, 0.7<DT<0.999
    INTEGER(I4), OPTIONAL,       INTENT(IN)  :: nITERmax_in ! maximal number of iterations
    INTEGER(I4), OPTIONAL,       INTENT(IN)  :: LEN_in   ! Length of Markov Chain, MAX(250, size(para,1))
    INTEGER(I4), OPTIONAL,       INTENT(IN)  :: nST_in   ! Number of consecutive LEN steps
    REAL(SP),    OPTIONAL,       INTENT(IN)  :: eps_in   ! epsilon decreement of cost function
    REAL(SP),    OPTIONAL,       INTENT(IN)  :: acc_in   ! Acceptance Ratio, <0.1 stopping criteria
    INTEGER(I8), OPTIONAL,       INTENT(IN)  :: seeds_in(3)   ! Seeds of random numbers
    LOGICAL,     OPTIONAL,       INTENT(IN)  :: printflag_in ! If command line output is written (.true.) 

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
    real(SP)                     :: RN1, RN2, RN3     ! Random numbers
    integer(I4)                  :: iin1, iin2, iin3  ! optional arguments for restarting RN streams
    integer(I4)                  :: win1, win2, win3  ! optional arguments for restarting RN streams
    integer(I4), dimension(0:127):: xin1, xin2,xin3   ! optional arguments for restarting RN streams

    ! for SA
    integer(I8)            :: idummy
    character(10)          :: timeseed
    real(SP)               :: NormPhi
    real(SP)               :: ac_ratio, pa
    real(SP)               :: fo, fn, df, fBest, rho, fInc, fbb, dr
    real(SP)               :: T0, DT0
    real(SP),  parameter   :: small = -700._SP
    integer(I4)            :: iPar
    integer(I4)            :: nITER
    integer(I4)            :: j, iter, kk
    integer(I4)            :: Ipos, Ineg
    integer(I4)            :: iConL, iConR, iConF
    integer(I4)            :: iTotalCounter                        ! includes reheating for final conditions
    integer(I4)            :: iTotalCounterR                       ! counter of interations in one reheating
    logical(1)             :: iStop 

    n = size(para,1)

    ! Input check
    if (size(para,2) .ne. 3) then
       stop 'Input argument para must have 3 columns (initial parameter value, min value, max value)'
    end if

    if (present(DT_in)) then
       if ( (DT_in .lt. 0.7_SP) .or. (DT_in .gt. 0.999_SP) ) then
          stop 'Input argument DT must lie between 0.7 and 0.999'
       else
          DT = DT_in
       end if
    else
       DT = 0.9_SP
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
          stop 'Input argument LEN must be greater than Max(250,20*N), N=number of parameters'
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
       if (eps_in .le. 0.0_SP) then
          stop 'Input argument eps must be greater than 0'
       else
          eps = eps_in
       end if
    else
       eps = 0.01_SP
    endif

    if (present(acc_in)) then
       if ((acc_in .le. 0.0_SP) .or. (acc_in .ge. 1.0_SP)) then
          stop 'Input argument acc must lie between 0.0 and 1.0'
       else
          acc = acc_in
       end if
    else
       acc = 0.1_i4
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

    T = temp

    call xor4096(seeds(1), RN1, iin=iin1, win=win1, xin=xin1)
    call xor4096(seeds(2), RN2, iin=iin2, win=win2, xin=xin2)
    call xor4096(seeds(3), RN3, iin=iin3, win=win3, xin=xin3)
    seeds = 0_i8

    ! Start Simulated Annealing routine
    gamma(:)%min = para(:,2)
    gamma(:)%max = para(:,3)
    gamma(:)%dmult = 1.0_SP
    gamma(:)%new = para(:,1)
    gamma(:)%old = para(:,1)
    NormPhi = -9999.9_SP
    T0=      T
    DT0=     DT
    nITER =  2_i4*N
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
    fBest   = 1.0_SP   != fo/NormPhi

    if (printflag) then
       print '(A15,E15.7,A4,E15.7,A4)',     ' start NSe   = ', fBest , '  ( ',fBest*normPhi,' )  '
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
       !LEN = int( real(iter,SP)*sqrt(real(iter,SP)) / nITER + 1.5*real(N,SP),i4) 
       ! Repeat LEN times
       loopLEN: do j=1,LEN
         iTotalCounterR =  iTotalCounterR + 1_i4
         iTotalCounter = iTotalCounter + 1_i4
         if ( (iTotalCounter .eq. 1_i4) .and. printflag)  then
           print '(I8, 2I5, E15.7,3E15.7)',   ITotalCounter, Ipos, Ineg, ac_ratio, T, fo, fBest
         end if
         if (mod(iTotalCounter,1000_i4) == 0_i4) then
           !call writeresults_o(4,fn) 
         end if
         !
	 ! Generate a random subsequent state and evaluate its objective function
         ! (1) Generate new parameter set
         dR=(1.0_SP - real(iTotalCounterR,SP) / real(nITERmax,SP))**2.0_SP
         if (.not. (dR >= 0.05_SP ) .and. (iTotalCounterR <= int(real(nIterMax,SP)/3._SP*4_SP))) dR = 0.05_SP   
         ! (1a) Select parameter to be changed    
         call xor4096(seeds(2),RN2, iin=iin2, win=win2, xin=xin2)
         iPar = int(real(N,SP)*RN2+1._SP,i4)
         ! (1b) Generate new value of selected parameter
         call xor4096(seeds(3),RN3, iin=iin3, win=win3, xin=xin3)
         gamma(iPar)%new = parGen_SP(gamma(iPar)%old,  gamma(iPar)%dMult*dR, gamma(iPar)%min, gamma(iPar)%max, 8_i4, 1_i4,RN3)
         ! (2) Calculate new objective function value and normalize it
         fn = cost(gamma(:)%new) 
         fn = fn/normPhi
         !
         ! write parameters
         !call writeresults_o(5,fn)
         !
         df = fn-fo             
         !
         ! analyze change in the objective function: df
         if (df < 0.0_SP) then
           !
           ! accept the new state
	   Ipos=Ipos+1_i4
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
             if (rho < small) then
               pa=0.0_SP
             else
               pa=EXP(rho)
             end if                     
  	     !
             call xor4096(seeds(1), RN1, iin=iin1, win=win1, xin=xin1)
             !
             if (pa > RN1) then
               ! accept new state with certain probability
               Ineg=Ineg+1_i4
               fo = fn
               ! save old state
               gamma(:)%old   = gamma(:)%new 
             end if
           end if
	 end if

      end do loopLEN
      !
      ! estimate acceptance ratio
      ac_ratio=(real(Ipos,SP) + real(Ineg,SP))/real(LEN,SP)
      !    
      if (printflag) then
         print '(I8, 2I5, E15.7, 3E15.7)', ITotalCounter, Ipos, Ineg, ac_ratio, T, fo, fBest
      end if
      !
      ! Cooling schedule
      fInc= (fbb-fBest)/fbb
      if (fInc < 0.00000001_SP) then
        iConF= iConF+1_i4
      else
        iConF=0_i4
      end if
      !
      if ((ac_ratio < 0.15_SP) .and. (iConF > 5_i4) .and. (iConR <= -3_i4)) then     ! - iConR  no reheating
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
        if (ac_ratio < 0.4_SP)  then 
          DT=0.995_SP
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
         nITERmax = int(real(nITERmax,sp)* 1.10_sp,i4)
         if (printflag) then
            print *,                       'nITERmax changed to =', nITERmax
         end if
      end if
      ! STOP condition
      if ( iTotalCounter > nITERmax ) iStop=.FALSE.                                    
      !
    end do loopTest
    !
    ! heuristic *global* minimum (end solution)
    !
    ! write final results
    if (printflag) then
       print *, '   '
       print '(A15,E15.7,A4,E15.7,A4)',       ' end NSe     = ', fBest , '  ( ',fBest*normPhi,' )  '
       print *,           'end parameter: '
       do kk = 1,N
          print '(A10,I3,A3, E15.7)' ,    '    para #',kk,' = ',gamma(kk)%best
       end do
    
       print *, 'Final check:    ', (fo - fBest)
    end if

    parabest = gamma(:)%best
    costbest = cost(parabest)

  END SUBROUTINE anneal_sp

!***************************************************
!*               PRIVATE FUNCTIONS                 *
!***************************************************

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

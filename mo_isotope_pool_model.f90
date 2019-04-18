!> \file mo_isotope_pool_model.f90

!> \brief Generic isotopic pool model and generic land-use change model

!> \details Explicit solution of a generic pool model,
!> given sink, sources and all the fluxes between all pools,
!> together with associated fractionation factors.

!> Also explicit solution of a generic land-use change model,
!> given all the fluxes and area changes between land use classes.

!> \author Matthias Cuntz
!> \date Apr 2019

MODULE mo_isotope_pool_model

  implicit none

  private

  public :: isotope_luc_model  ! Solves the next time step of a generic isotopic land-use change model
  public :: isotope_pool_model ! Solves the next time step of a generic isotopic pool model

  ! ------------------------------------------------------------------

  private :: diag
  interface diag
     module procedure diag_2d, diag_3d
  end interface diag  
  
  ! ------------------------------------------------------------------

  !     NAME
  !         isotope_pool_model

  !     PURPOSE
  !>        \brief Generic isotopic pool model

  !>        \details Next explicit time step of the isotopic composition of a generic pool model,
  !>        given sink, sources and all the fluxes between all pools,
  !>        together associated fractionation factors.

  !>        The pool model for pool i is:
  !>            dC(i)/dt = sum(F(:,i)) - sum(F(i,:)) + S(i) - Si(i)
  !>        where C are the pools, F the fluxes between pools (F>0, F(i,i)=0),
  !>        S a pool-independent source term, and Si a sink term, which can be pool-independent
  !>        or pool-dependent, i.e. Si(i) = beta(i)*C(i)

  !>        The isotope model for pool i is then:
  !>            dC'(i)/dt = sum(alpha(:,i)*R(:)*F(:,i)) - R(i)*sum(alpha(i,:)*F(i,:)) + Rs(i)*S(i) - alpha(i,i)*R(i)*Si(i)
  !>        with C' the isotope pools. R are the isotope ratios of the pools at time t-dt,
  !>        and alpha are possible fractionation factors.

  !>        All pools can have either dimension (n) or (n,m) with n pools and m land points (or fractions).
  !>        Fluxes and fractionation factors have dimensions (n,n) or (n,n,m).

  !     CALLING SEQUENCE
  !         call isotope_pool_model(dt, Ci, C, F, S, Rs, Si, alpha, beta, trash)

  !     INTENT
  !>        \param[in]    "real(dp) :: dt"                       Time step
  !>        \param[inout] "real(dp) :: Ci(1:n[,1:m])"            Isotope concentrations in pools
  !>        \param[in]    "real(dp) :: C(1:n[,1:m])"             Non-isotope concentrations in pools at time step t-1
  !>        \param[in]    "real(dp) :: F(1:n,1:n[,1:m])"         Non-isotope fluxes between pools (F(i,i)=0)
  !>                                   F(i,j) positive flux from pool i to pool j
  !>                                   Isotope flux is alpha(i,j)*R(i)*F(i,j)
  !>        \param[in]    "real(dp), optional :: S(1:n[,1:m])"   Non-isotope source fluxes to pools (other than between pools)
  !>                                   Default: 0.
  !>        \param[in]    "real(dp), optional :: Rs(1:n[,1:m])"  Isotopic compositions of source fluxes
  !>                                   Any fractionations during source processes should be included in Rs.
  !>                                   Default: 1.
  !>        \param[in]    "real(dp), optional :: Si(1:n[,1:m])"  Non-isotope sinks of pools (other than between pools)
  !>                                   Isotope flux is alpha(i,i)*R(i)*Si(i)
  !>                                   Default: 0.
  !>        \param[in]    "real(dp), optional :: alpha(1:n,1:n[,1:m])" Isotopic fractionation factors associated with
  !>                                   fluxes F between pools (alpha(i,j) i/=j) and of 
  !>                                   sinks Si (other than between pools) (alpha(i,i))
  !>                                   Default: 1.
  !>        \param[in]    "real(dp), optional :: beta(1:n[,1:m])"      Either Si or beta can be given.
  !>                                   If beta is given then Si(i)=beta(i)*C(i).
  !>                                   Si supercedes beta, i.e. Si will be taken if beta and Si are given.
  !>        \param[inout] "real(dp), optional :: trash(1:n[,1:m])"     Container to store possible inconsistencies,
  !>                                   might be numeric, between non-isotope and isotope model.

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  interface isotope_pool_model
     module procedure isotope_pool_model_1d, isotope_pool_model_2d
  end interface isotope_pool_model

contains

  ! ------------------------------------------------------------------

  !     NAME
  !         isotope_luc_model

  !     PURPOSE
  !>        \brief Generic isotopic land-use change model

  !>        \details Next explicit time step of the isotopic composition of a generic land-use change model,
  !>        given all fluxes and area changes between land use classes.

  !>        The land-use change model for land-use class i is:
  !>            d(A(i)*C(i)) = -C(i)*sum(dA(i,:)) + sum(dA(:,i)*C(:))
  !>        where C are the carbon concentrations in the land-use classes, A are the areas of the land-use classes
  !>        and dA are the area changes between land-use classes for the time step (A>0, A(i,i)=0).

  !>        Note: dA is change per time step, NOT per second.

  !>        The land-use model is then solved from time step t to t+dt
  !>            A(i,t+dt)*C(i,t+dt) - A(i,t)*C(i,t) = -C(i,t)*sum(dA(i,:)) + sum(dA(:,i)*C(:,t))
  !>            C(i,t+dt) = (C(i,t)*(A(i,t)-sum(dA(i,:))) + sum(dA(:,i)*C(:,t))) / A(i,t+dt)
  
  !>        The isotopic land-use model is then:
  !>            A(i,t+dt)*C'(i,t+dt) - A(i,t)*C'(i,t) = -C'(i,t)*sum(dA(i,:)) + sum(dA(:,i)*C'(:,t))
  !>            C'(i,t+dt) = (C'(i,t)*(A(i,t)-sum(dA(i,:))) + sum(dA(:,i)*C'(:,t))) / A(i,t+dt)
  !>        with C' the isotope carbon concentrations.

  !>        The use can also give the carbon concentrations C of the land use classes
  !>        at time t
  !>            R(i,t) = C'(i,t) / C(i,t)
  !>            A(i,t+dt)*C'(i,t+dt) - A(i,t)*R(i,t)*C(i,t) = -R(i,t)*C(i,t)*sum(dA(i,:)) + sum(dA(:,i)*R(:,t)*C(:,t))
  !>            C'(i,t+dt) = (R(i,t)*C(i,t)*(A(i,t)-sum(dA(i,:))) + sum(dA(:,i)*R(:,t)*C(:,t))) / A(i,t+dt)

  !     CALLING SEQUENCE
  !         call isotope_luc_model(Ci, A, dA, C, trash)

  !     INTENT
  !>        \param[in]    "real(dp) :: dt"                       Time step
  !>        \param[inout] "real(dp) :: Ci(1:n)"
  !>                                   On entry: isotope concentrations in land-use classes at time step t-1
  !>                                   On exit:  updated isotope concentrations of land-use classes at time step t
  !>        \param[in]    "real(dp) :: A(1:n)"                   Areas of land-use classes at time step t-1
  !>        \param[in]    "real(dp) :: dA(1:n,1:n)"              Change from each land-use class to the others
  !>                                   during the time step
  !>        where C are the carbon concentrations in the land-use classes, A are the areas of the land-use classes
  !>        and dA are the area changes between land-use classes (A>0, A(i,i)=0),
  !>        \param[in]    "real(dp), optional :: C(1:n)"         Non-isotope concentrations in land-use classes
  !>                                   at time step t-1
  !>        \param[inout] "real(dp), optional :: trash(1:n)"     Container to store possible inconsistencies,
  !>                                   might be numeric, between non-isotope and isotope model.

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  subroutine isotope_luc_model(Ci, A, dA, C, trash)

    use mo_kind,  only: dp, i4
    use mo_utils, only: eq, ne

    implicit none

    real(dp), dimension(:),   intent(inout)           :: Ci     ! Iso land-use class
    real(dp), dimension(:),   intent(in)              :: A      ! Land area land-use class
    real(dp), dimension(:,:), intent(in)              :: dA     ! Land area changes between land-use classes
    real(dp), dimension(:),   intent(in),    optional :: C      ! Non-iso land-use class
    real(dp), dimension(:),   intent(inout), optional :: trash  ! garbage can for numerical inconsistencies

    ! Local variables
    integer(i4) :: i  ! counter
    integer(i4) :: nn ! number of pools
    real(dp), dimension(size(Ci,1)) :: R    ! Isotope ratio of pool
    real(dp), dimension(size(Ci,1)) :: Anew ! A at t+dt
    ! defaults for optional inputs
    real(dp), dimension(size(Ci,1)) :: iC, itrash

    ! Check sizes
    nn = size(Ci,1)
    if ( (size(A,1) /= nn) .or. (size(dA,1) /= nn) .or. (size(dA,2) /= nn) ) then
       write(*,*) 'Error isotope_luc_model: non-fitting dimensions between isotopic concentrations,'
       write(*,*) '                         land areas and land-area changes.'
       write(*,*) '    size(Ci):   ', size(Ci,1)
       write(*,*) '    size(A):    ', size(A,1)
       write(*,*) '    size(dA,1): ', size(dA,1)
       write(*,*) '    size(dA,2): ', size(dA,2)
       stop 9
    endif

    ! Check dA >= 0
    if (any(dA < 0._dp)) then
       write(*,*) 'Error isotope_luc_model: land area changes between land use classes must be >= 0.'
       write(*,*) '    dA: ', dA
       stop 9
    endif

    ! Check dA(i,:) == 0. if Ci(i) == 0.
    if (any(eq(Ci,0._dp))) then
       do i=1, nn
          if (eq(Ci(i),0._dp)) then
             if (any(ne(dA(i,:),0._dp))) then
                write(*,*) 'Error isotope_luc_model: land area changes from land-use class i must be 0'
                write(*,*) '                         if isotope carbon concentration of land-use class is 0.'
                write(*,*) '    i, Ci(i): ', i, Ci(i)
                write(*,*) '       dA(i): ', dA(i,:)
                stop 9
             endif
          endif
       end do
    endif

    if (present(C)) then
       ! Check dA(i,:) == 0. if C(i) == 0.
       if (any(eq(C,0._dp))) then
          do i=1, nn
             if (eq(C(i),0._dp)) then
                if (any(ne(dA(i,:),0._dp))) then
                   write(*,*) 'Error isotope_luc_model: land area changes from land-use class i must be 0'
                   write(*,*) '                         if carbon concentration of land-use class is 0.'
                   write(*,*) '     i, C(i): ', i, C(i)
                   write(*,*) '       dA(i): ', dA(i,:)
                   stop 9
                endif
             endif
          end do
       endif
    endif
    
    ! Set optionals
    if (present(C)) then
       iC = C
    else
       iC = 1._dp
    endif
    if (present(trash)) then
       itrash = trash
    else
       itrash = 0._dp
    endif

    ! Isotope ratio
    R(:) = 1._dp
    where (iC > 0._dp) R = Ci / iC

    ! Land areas at t+dt
    Anew = A - sum(dA, dim=2) + sum(dA, dim=1)

    ! isotopic LUC model
    Ci = R * iC * (A - sum(dA, dim=2)) + sum(dA * spread(R*iC, dim=2, ncopies=nn), dim=1)
    where (Anew > 0._dp)
       Ci = Ci / Anew
    elsewhere
       Ci     = 0._dp
       itrash = itrash + Ci
    endwhere

    ! Check final land-use classes
    ! Isotope land-use class became < 0.
    if (any(Ci < 0._dp)) then
       itrash = itrash + merge(abs(Ci), 0._dp, Ci < 0._dp)
       Ci = merge(0._dp, Ci, Ci < 0._dp)
    endif
    ! Non-isotope land-use class == 0. but isotope land-use class > 0.
    if (any(eq(iC,0._dp) .and. (Ci > 0._dp))) then
       itrash = itrash + merge(Ci, 0._dp, eq(iC,0._dp) .and. (Ci > 0._dp))
       Ci = merge(0._dp, Ci, eq(iC,0._dp) .and. (Ci > 0._dp))
    endif
    ! Non-isotope land-use class >0. but isotope land-use class == 0.
    ! ???

    if (present(trash)) trash = itrash

    return

  end subroutine isotope_luc_model

  
  ! ------------------------------------------------------------------
  
  subroutine isotope_pool_model_1d(dt, Ci, C, F, S, Rs, Si, alpha, beta, trash)

    use mo_kind,  only: dp, i4
    use mo_utils, only: eq, ne

    implicit none

    real(dp),                 intent(in)              :: dt     ! time step
    real(dp), dimension(:),   intent(inout)           :: Ci     ! Iso pool
    real(dp), dimension(:),   intent(in)              :: C      ! Non-iso pool
    real(dp), dimension(:,:), intent(in)              :: F      ! Fluxes between pools
    real(dp), dimension(:),   intent(in),    optional :: S      ! Sources not between pools
    real(dp), dimension(:),   intent(in),    optional :: Rs     ! Isotope ratio of S
    real(dp), dimension(:),   intent(in),    optional :: Si     ! Sinks not between pools
    real(dp), dimension(:,:), intent(in),    optional :: alpha  ! Fractionation factors sinks and fluxes between pools
    real(dp), dimension(:),   intent(in),    optional :: beta   ! Alternative to sinks = beta*Ct
    real(dp), dimension(:),   intent(inout), optional :: trash  ! garbage can for numerical inconsistencies

    ! Local variables
    integer(i4) :: i  ! counter
    integer(i4) :: nn ! number of pools
    real(dp), dimension(size(Ci,1)) :: R                 ! Isotope ratio of pool
    ! defaults for optional inputs
    real(dp), dimension(size(Ci,1))            :: iS, iRs, iSi, itrash
    real(dp), dimension(size(Ci,1),size(Ci,1)) :: ialpha
    real(dp), dimension(size(Ci,1),size(Ci,1)) :: alphaF ! alpha*F
    real(dp), dimension(size(Ci,1)) :: sink              ! not-between pools sink
    real(dp), dimension(size(Ci,1)) :: source            ! not-between pools source
    real(dp), dimension(size(Ci,1)) :: isink             ! between pools sink
    real(dp), dimension(size(Ci,1)) :: isource           ! between pools source

    ! Check sizes
    nn = size(Ci,1)
    if ( (size(C,1) /= nn) .or. (size(F,1) /= nn) .or. (size(F,2) /= nn) ) then
       write(*,*) 'Error isotope_pool_model_1d: non-fitting dimensions between isotopic pools,'
       write(*,*) '                             non-isotopic pools and fluxes.'
       write(*,*) '    size(Ci):  ', size(Ci,1)
       write(*,*) '    size(C):   ', size(C,1)
       write(*,*) '    size(F,1): ', size(F,1)
       write(*,*) '    size(F,2): ', size(F,2)
       stop 9
    endif

    ! Check F >= 0
    if (any(F < 0._dp)) then
       write(*,*) 'Error isotope_pool_model_1d: fluxes between pools must be >= 0.'
       write(*,*) '    F: ', F
       stop 9
    endif

    ! Check F(i,:) == 0. if C(i) == 0.
    if (any(eq(C,0._dp))) then
       do i=1, nn
          if (eq(C(i),0._dp)) then
             if (any(ne(F(i,:),0._dp))) then
                write(*,*) 'Error isotope_pool_model_1d: fluxes from pool i must be 0 if concentration in pool is 0.'
                write(*,*) '    i, C(i): ', i, C(i)
                write(*,*) '       F(i): ', F(i,:)
                stop 9
             endif
          endif
       end do
    endif

    ! Set optionals
    if (present(S)) then
       iS = S
    else
       iS = 0._dp
    endif
    if (present(Rs)) then
       iRs = Rs
    else
       iRs = 1._dp
    endif
    if (present(Si)) then
       iSi = Si
    else
       iSi = 0._dp
    endif
    if (present(alpha)) then
       ialpha = alpha
    else
       ialpha = 1._dp
    endif
    if (present(alpha)) then
       ialpha = alpha
    else
       ialpha = 1._dp
    endif
    if (present(beta) .and. (.not. present(Si))) then
       iSi = beta * C
    endif
    if (present(trash)) then
       itrash = trash
    else
       itrash = 0._dp
    endif

    ! Isotope ratio
    R(:) = 0._dp
    where (C > 0._dp) R = Ci / C

    ! alpha * F
    alphaF = ialpha * F

    ! between and not-between source and sinks
    sink    = diag(ialpha) * R * iSi * dt
    source  = iRs * iS * dt
    isink   = sum(alphaF, dim=2) * R * dt
    isource = sum(alphaF * spread(R, dim=2, ncopies=nn), dim=1) * dt

    ! Explicit solution
    Ci = Ci - sink + source - isink + isource

    ! Check final pools
    ! Isotope pool became < 0.
    if (any(Ci < 0._dp)) then
       itrash = itrash + merge(abs(Ci), 0._dp, Ci < 0._dp)
       Ci = merge(0._dp, Ci, Ci < 0._dp)
    endif
    ! Non-isotope pool == 0. but isotope pool > 0.
    if (any(eq(C,0._dp) .and. (Ci > 0._dp))) then
       itrash = itrash + merge(Ci, 0._dp, eq(C,0._dp) .and. (Ci > 0._dp))
       Ci = merge(0._dp, Ci, eq(C,0._dp) .and. (Ci > 0._dp))
    endif
    ! Non-isotope pool >0. but isotope pool == 0.
    ! ???

    if (present(trash)) trash = itrash

    return

  end subroutine isotope_pool_model_1d

  
  subroutine isotope_pool_model_2d(dt, Ci, C, F, S, Rs, Si, alpha, beta, trash)

    use mo_kind,  only: dp, i4
    use mo_utils, only: eq, ne

    implicit none

    real(dp),                   intent(in)              :: dt     ! time step
    real(dp), dimension(:,:),   intent(inout)           :: Ci     ! Iso pool
    real(dp), dimension(:,:),   intent(in)              :: C      ! Non-iso pool
    real(dp), dimension(:,:,:), intent(in)              :: F      ! Fluxes between pools
    real(dp), dimension(:,:),   intent(in),    optional :: S      ! Sources not between pools
    real(dp), dimension(:,:),   intent(in),    optional :: Rs     ! Isotope ratio of S
    real(dp), dimension(:,:),   intent(in),    optional :: Si     ! Sinks not between pools
    real(dp), dimension(:,:,:), intent(in),    optional :: alpha  ! Fractionation factors sinks and fluxes between pools
    real(dp), dimension(:,:),   intent(in),    optional :: beta   ! Alternative to sinks = beta*Ct
    real(dp), dimension(:,:),   intent(inout), optional :: trash  ! garbage can for numerical inconsistencies

    ! Local variables
    integer(i4) :: i, j  ! counter
    integer(i4) :: nland ! number of land points
    integer(i4) :: nn    ! number of pools
    real(dp), dimension(size(Ci,1),size(Ci,2)) :: R                 ! Isotope ratio of pool
    ! defaults for optional inputs
    real(dp), dimension(size(Ci,1),size(Ci,2))            :: iS, iRs, iSi, itrash
    real(dp), dimension(size(Ci,1),size(Ci,1),size(Ci,2)) :: ialpha
    real(dp), dimension(size(Ci,1),size(Ci,1),size(Ci,2)) :: alphaF ! alpha*F
    real(dp), dimension(size(Ci,1),size(Ci,2)) :: sink              ! not-between pools sink
    real(dp), dimension(size(Ci,1),size(Ci,2)) :: source            ! not-between pools source
    real(dp), dimension(size(Ci,1),size(Ci,2)) :: isink             ! between pools sink
    real(dp), dimension(size(Ci,1),size(Ci,2)) :: isource           ! between pools source

    ! Check sizes
    nn = size(Ci,1)
    if ( (size(C,1) /= nn) .or. (size(F,1) /= nn) .or. (size(F,2) /= nn) ) then
       write(*,*) 'Error isotope_pool_model_2d: non-fitting dimensions between isotopic pools,'
       write(*,*) '                             non-isotopic pools and fluxes.'
       write(*,*) '    size(Ci,1): ', size(Ci,1)
       write(*,*) '    size(C,1):  ', size(C,1)
       write(*,*) '    size(F,1):  ', size(F,1)
       write(*,*) '    size(F,2):  ', size(F,2)
       stop 9
    endif
    nland = size(Ci,2)
    if ( (size(C,2) /= nland) .or. (size(F,3) /= nland) ) then
       write(*,*) 'Error isotope_pool_model_2d: non-fitting first dimensions between isotopic pools,'
       write(*,*) '                             non-isotopic pools and fluxes.'
       write(*,*) '    size(Ci,2): ', size(Ci,2)
       write(*,*) '    size(C,2):  ', size(C,2)
       write(*,*) '    size(F,3):  ', size(F,3)
       stop 9
    endif

    ! Check F >= 0
    if (any(F < 0._dp)) then
       write(*,*) 'Error isotope_pool_model_2d: fluxes between pools must be >= 0.'
       write(*,*) '    F: ', F
       stop 9
    endif

    ! Check F(i,:) == 0. if C(i) == 0.
    if (any(eq(C,0._dp))) then
       do j=1, nland
          do i=1, nn
             if (eq(C(i,j),0._dp)) then
                if (any(ne(F(i,:,j),0._dp))) then
                   write(*,*) 'Error isotope_pool_model_2d:'
                   write(*,*) '    fluxes from pool i at land point j must be 0 if concentration in pool is 0.'
                   write(*,*) '    i, j, C(i,j):   ', i, j, C(i,j)
                   write(*,*) '          F(i,:,j): ', F(i,:,j)
                   stop 9
                endif
             endif
          end do
       end do
    endif

    ! Set optionals
    if (present(S)) then
       iS = S
    else
       iS = 0._dp
    endif
    if (present(Rs)) then
       iRs = Rs
    else
       iRs = 1._dp
    endif
    if (present(Si)) then
       iSi = Si
    else
       iSi = 0._dp
    endif
    if (present(alpha)) then
       ialpha = alpha
    else
       ialpha = 1._dp
    endif
    if (present(alpha)) then
       ialpha = alpha
    else
       ialpha = 1._dp
    endif
    if (present(beta) .and. (.not. present(Si))) then
       iSi = beta * C
    endif
    if (present(trash)) then
       itrash = trash
    else
       itrash = 0._dp
    endif

    ! Isotope ratio
    R(:,:) = 0._dp
    where (C > 0._dp) R = Ci / C

    ! alpha * F
    alphaF = ialpha * F

    ! between and not-between source and sinks
    sink    = diag(ialpha) * R * iSi * dt
    source  = iRs * iS * dt
    isink   = sum(alphaF, dim=2) * R * dt
    isource = sum(alphaF * spread(R, dim=2, ncopies=nn), dim=1) * dt
    ! Explicit solution
    Ci = Ci - sink + source - isink + isource

    ! Check final pools
    ! Isotope pool became < 0.
    if (any(Ci < 0._dp)) then
       itrash = itrash + merge(abs(Ci), 0._dp, Ci < 0._dp)
       Ci = merge(0._dp, Ci, Ci < 0._dp)
    endif
    ! Non-isotope pool == 0. but isotope pool > 0.
    if (any(eq(C,0._dp) .and. (Ci > 0._dp))) then
       itrash = itrash + merge(Ci, 0._dp, eq(C,0._dp) .and. (Ci > 0._dp))
       Ci = merge(0._dp, Ci, eq(C,0._dp) .and. (Ci > 0._dp))
    endif
    ! Non-isotope pool >0. but isotope pool == 0.
    ! ???

    if (present(trash)) trash = itrash

    return

  end subroutine isotope_pool_model_2d
  
  ! ------------------------------------------------------------------

  ! Diagonal elements of a matrix
  function diag_2d(matrix)

    use mo_kind, only: dp, i4

    implicit none

    real(dp), dimension(:,:), intent(in) :: matrix
    real(dp), dimension(:), allocatable  :: diag_2d

    integer(i4) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_2d: array must be squared matrix.'
    if (.not. allocated(diag_2d)) allocate(diag_2d(size(matrix,1)))

    forall(i=1:size(matrix,1)) diag_2d(i) = matrix(i,i)

  end function diag_2d

  ! Diagonal elements of the two first dimensions of a matrix
  function diag_3d(matrix)

    use mo_kind, only: dp, i4

    implicit none

    real(dp), dimension(:,:,:), intent(in) :: matrix
    real(dp), dimension(:,:), allocatable  :: diag_3d

    integer(i4) :: i, j

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_3d: array must be squared matrix in the first and second dimensions.'
    if (.not. allocated(diag_3d)) allocate(diag_3d(size(matrix,1),size(matrix,3)))

    do j=1, size(matrix,3)
       forall(i=1:size(matrix,1)) diag_3d(i,j) = matrix(i,i,j)
    end do

  end function diag_3d

END MODULE mo_isotope_pool_model

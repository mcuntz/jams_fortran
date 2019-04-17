!> \file mo_isotope_pool_model.f90

!> \brief Generic isotopic pool model

!> \details Explicit solution of a generic pool model,
!> given sink, sources and all the fluxes between all pools,
!> together associated fractionation factors.

!> \author Matthias Cuntz
!> \date Apr 2019

MODULE mo_isotope_pool_model

  implicit none

  private

  public :: isotope_pool_model ! Solves the next time step of the generic isotope model

contains

  ! ------------------------------------------------------------------

  !     NAME
  !         isotope_pool_model

  !     PURPOSE
  !>        \brief Generic isotopic pool model

  !>        \details Next explicit time step of the isotopic composition of a generic pool model,
  !>        given sink, sources and all the fluxes between all pools,
  !>        together associated fractionation factors.

  !>        The pool model is:
  !>            dC(i)/dt = sum(F(:,i)) - sum(F(i,:)) + S(i) - Si(i)
  !>        where C are the pools, F the fluxes between pools (F>0, F(i,i)=0),
  !>        S a pool-independent source term, and Si a sink term, which can be pool-independent
  !>        or pool-dependent, i.e. Si(i) = beta(i)*C(i)

  !>        The isotope model is then:
  !>            dC'(i)/dt = sum(alpha(:,i)*R(:)*F(:,i)) - R(i)*sum(alpha(i,:)*F(i,:)) + Rs(i)*S(i) - alpha(i,i)*R(i)*Si(i)
  !>        with C' the isotope pools. R are the isotope ratios of the pools at time t-dt,
  !>        and alpha are possible fractionation factors.

  !     CALLING SEQUENCE
  !         call isotope_pool_model(dt, Ci, C, F, S, Rs, Si, alpha, beta, trash)

  !     INTENT
  !>        \param[in]    "real(dp) :: dt"                       Time step
  !>        \param[inout] "real(dp) :: Ci(1:n)"                  Isotope concentrations in pools
  !>        \param[in]    "real(dp) :: C(1:n)"                   Non-isotope concentrations in pools at time step t-1
  !>        \param[in]    "real(dp) :: F(1:n,1:n)"               Non-isotope fluxes between pools (F(i,i)=0)
  !>                                   F(i,j) positive flux from pool i to pool j
  !>                                   Isotope flux is alpha(i,j)*R(i)*F(i,j)
  !>        \param[in]    "real(dp), optional :: S(1:n)"         Non-isotope source fluxes to pools (other than between pools)
  !>                                   Default: 0.
  !>        \param[in]    "real(dp), optional :: Rs(1:n)"        Isotopic compositions of source fluxes
  !>                                   Any fractionations during source processes should be included in Rs.
  !>                                   Default: 1.
  !>        \param[in]    "real(dp), optional :: Si(1:n)"        Non-isotope sinks of pools (other than between pools)
  !>                                   Isotope flux is alpha(i,i)*R(i)*Si(i)
  !>                                   Default: 0.
  !>        \param[in]    "real(dp), optional :: alpha(1:n,1:n)" Isotopic fractionation factors associated with
  !>                                   fluxes F between pools (alpha(i,j) i/=j) and of 
  !>                                   sinks Si (other than between pools) (alpha(i,i))
  !>                                   Default: 1.
  !>        \param[in]    "real(dp), optional :: beta(1:n)"      Either Si or beta can be given.
  !>                                   If beta is given then Si(i)=beta(i)*C(i).
  !>                                   Si supercedes beta, i.e. Si will be taken if beta and Si are given.
  !>        \param[inout] "real(dp), optional :: trash(1:n)"     Container to store possible inconsistencies,
  !>                                   might be numeric, between non-isotope and isotope model.

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  subroutine isotope_pool_model(dt, Ci, C, F, S, Rs, Si, alpha, beta, trash)

    use mo_kind,           only: dp, i4
    use mo_utils,          only: eq, ne
    use mo_linear_algebra, only: diag

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
       write(*,*) 'Error isotope_pool_model: non-fitting dimensions betwwen isotopic pools, non-isotopic pools and fluxes.'
       write(*,*) '    size(Ci):  ', size(Ci,1)
       write(*,*) '    size(C):   ', size(C,1)
       write(*,*) '    size(F,1): ', size(F,1)
       write(*,*) '    size(F,2): ', size(F,2)
       stop 9
    endif
    
    ! Check F >= 0
    if (any(F < 0._dp)) then
       write(*,*) 'Error isotope_pool_model: fluxes between pools must be >= 0.'
       write(*,*) '    F: ', F
       stop 9
    endif

    ! Check F(i,:) == 0. if C(i) == 0.
    if (any(eq(C,0._dp))) then
       do i=1, nn
          if (eq(C(i),0._dp)) then
             if (any(ne(F(i,:),0._dp))) then
                write(*,*) 'Error isotope_pool_model: fluxes from pool i must be 0 if concentration in pool is 0.'
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

    return

end subroutine isotope_pool_model

END MODULE mo_isotope_pool_model
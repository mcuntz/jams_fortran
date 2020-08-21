! ============================================================================
! Name        : mo_eq_lv.f90
! Author      : Giovanni Dalmasso
! Version     : 1.0
! Copyright   : Department of Computational Hydrosystems CHS - UZF Leipzig
! Description : Lotka Volterra equations
! ============================================================================

module mo_eq_lv

    ! ------------------------------------------------------------------
    ! The Lotka Volterra equations, also known as the predator-prey equations, are a pair of first-order, non-linear,
    ! differential equations frequently used to describe the dynamics of biological systems in which two species interact,
    ! one a predator and one its prey. They evolve in time according to the following of equations
    ! ------------------------------------------------------------------

    implicit none

    private

    public :: LV_eqn, LV_eqn_sp, LV_eqn_dp, &
         LV_eqn_para, LV_eqn_para_sp, LV_eqn_para_dp, &
         deriv_testRB_sp, jacobn_testRB_sp, &
         deriv_testRB_dp, jacobn_testRB_dp, &
         deriv2_testRB_sp, jacobn2_testRB_sp, &
         deriv2_testRB_dp, jacobn2_testRB_dp, &
         deriv2_testRB_para_sp, jacobn2_testRB_para_sp, &
         deriv2_testRB_para_dp, jacobn2_testRB_para_dp


    ! Interfaces for single and double precision routines
    interface LV_eqn
        module procedure LV_eqn_sp, LV_eqn_dp
    end interface

    interface LV_eqn_para
        module procedure LV_eqn_para_sp, LV_eqn_para_dp
    end interface

contains

    ! SINGLE PRECISION
    subroutine LV_eqn_sp( x, y, dydx )    ! returns derivatives dydx at x

        use mo_kind, only : sp

        implicit none

        real(sp), intent(in) :: x                       ! time
        real(sp), dimension(:), intent(in) :: y         ! unknowns of the equations
        real(sp), dimension(:), intent(out) :: dydx     ! derivatives of y

        real(sp), parameter :: epsilonn = 0.5_sp        ! regulates the relationship between the two species
                                                        ! keeping the optimum balance

        real(sp), parameter :: mu = 2._sp               ! represents the ability to reproduce of the first species
                                                        ! with respect to the destructive capacity of the second one

        real(sp), parameter :: eta = -1._sp             ! regulates the relationship between the influences of nutrition
                                                        ! related reproductive capacity

        ! ============================================================================
        ! System of ODE
        !
        ! y(1) --> predators
        ! y(2) --> preys
        ! does not depend on time  --> 0.0 * x (to avoid warnings of unused variables)
        ! ============================================================================

        dydx(1) = y(1) - epsilonn*mu*y(2)        + 0.0_sp*x
        dydx(2) = (mu/epsilonn)*y(1) + eta*y(2)  + 0.0_sp*x

    end subroutine LV_eqn_sp

    ! DOUBLE PRECISION
    subroutine LV_eqn_dp( x, y, dydx )   ! returns derivatives dydx at x

        use mo_kind, only : dp

        implicit none

        real(dp), intent(in) :: x                       ! time
        real(dp), dimension(:), intent(in) :: y         ! unknowns of the equations
        real(dp), dimension(:), intent(out) :: dydx     ! derivatives of y

        real(dp), parameter :: epsilonn = 0.5_dp        ! regulates the relationship between the two species
                                                        ! keeping the optimum balance

        real(dp), parameter :: mu = 2._dp               ! represents the ability to reproduce of the first species
                                                        ! with respect to the destructive capacity of the second one

        real(dp), parameter :: eta = -1._dp             ! regulates the relationship between the influences of nutrition
                                                        ! related reproductive capacity

        ! ============================================================================
        ! System of ODE
        !
        ! y(1) --> predators
        ! y(2) --> preys
        ! does not depend on time  --> 0.0 * x (to avoid warnings of unused variables)
        ! ============================================================================

        dydx(1) = y(1) - epsilonn*mu*y(2)        + 0.0_dp*x
        dydx(2) = (mu/epsilonn)*y(1) + eta*y(2)  + 0.0_dp*x

    end subroutine LV_eqn_dp

!-----------
! with parameter input
!-----------

    ! SINGLE PRECISION
    subroutine LV_eqn_para_sp( x, y, para, dydx )    ! returns derivatives dydx at x

        use mo_kind, only : sp

        implicit none

        real(sp), intent(in) :: x                       ! time
        real(sp), dimension(:), intent(in) :: y         ! unknowns of the equations
        real(sp), dimension(:), intent(out) :: dydx     ! derivatives of y

        real(sp), dimension(:), intent(in) :: para      ! parameter

!        real(sp), parameter :: epsilonn = 0.5_sp        ! regulates the relationship between the two species
!                                                        ! keeping the optimum balance

!        real(sp), parameter :: mu = 2._sp               ! represents the ability to reproduce of the first species
!                                                        ! with respect to the destructive capacity of the second one

!        real(sp), parameter :: eta = -1._sp             ! regulates the relationship between the influences of nutrition
!                                                        ! related reproductive capacity

        ! ============================================================================
        ! System of ODE
        !
        ! y(1) --> predators
        ! y(2) --> preys
        ! does not depend on time  --> 0.0 * x (to avoid warnings of unused variables)
        ! ============================================================================

        dydx(1) = y(1) - para(1)*para(2)*y(2)           + 0.0_sp*x
        dydx(2) = (para(2)/para(1))*y(1) + para(3)*y(2) + 0.0_sp*x

    end subroutine LV_eqn_para_sp

    ! DOUBLE PRECISION
    subroutine LV_eqn_para_dp( x, y, para, dydx )   ! returns derivatives dydx at x

        use mo_kind, only : dp

        implicit none

        real(dp), intent(in) :: x                       ! time
        real(dp), dimension(:), intent(in) :: y         ! unknowns of the equations
        real(dp), dimension(:), intent(out) :: dydx     ! derivatives of y

        real(dp), dimension(:), intent(in) :: para      ! parameter

!        real(dp), parameter :: epsilonn = 0.5_dp        ! regulates the relationship between the two species
!                                                        ! keeping the optimum balance

!        real(dp), parameter :: mu = 2._dp               ! represents the ability to reproduce of the first species
!                                                        ! with respect to the destructive capacity of the second one

!        real(dp), parameter :: eta = -1._dp             ! regulates the relationship between the influences of nutrition
!                                                        ! related reproductive capacity

        ! ============================================================================
        ! System of ODE
        !
        ! y(1) --> predators
        ! y(2) --> preys
        ! does not depend on time  --> 0.0 * x (to avoid warnings of unused variables)
        ! ============================================================================

        dydx(1) = y(1) - para(1)*para(2)*y(2)           + 0.0_dp*x
        dydx(2) = (para(2)/para(1))*y(1) + para(3)*y(2) + 0.0_dp*x

    end subroutine LV_eqn_para_dp


!-----------
! stiff sets of ODEs
!-----------

!SINGLE PRECISION
    subroutine deriv_testRB_sp( x, y, dydx )   ! returns derivatives dydx at x

        use mo_kind, only : sp

        implicit none

        real(sp), intent(in) :: x                       ! time
        real(sp), dimension(:), intent(in) :: y         ! unknowns of the equations
        real(sp), dimension(:), intent(out) :: dydx     ! derivatives of y

        real(sp) :: xtmp

        xtmp = x

        dydx(1) = 998.0_sp*y(1) + 1998.0_sp*y(2)
        dydx(2) = -999.0_sp*y(1)- 1999.0_sp*y(2)

    end subroutine deriv_testRB_sp

    subroutine jacobn_testRB_sp( x, y, dfdx, dfdy )   ! returns derivatives dydx at x

        use mo_kind, only : sp

        implicit none

        real(sp), intent(in) :: x                           ! time
        real(sp), dimension(:),     intent(in)  :: y        ! unknowns of the equations
        real(sp), dimension(:),     intent(out) :: dfdx     ! derivatives of f for x
        real(sp), dimension(:,:),   intent(out) :: dfdy     ! derivatives of f for y

        real(sp) :: xtmp
        real(sp), dimension(size(y,1)) :: ytmp

        xtmp = x
        ytmp = y

        dfdx(1)     = 0.0_sp
        dfdx(2)     = 0.0_sp

        dfdy(1,1)   = 998.0_sp
        dfdy(1,2)   = 1998.0_sp
        dfdy(2,1)   = -999.0_sp
        dfdy(2,2)   = 1999.0_sp

    end subroutine jacobn_testRB_sp

    subroutine deriv2_testRB_sp( x, y, dydx )   ! returns derivatives dydx at x

        use mo_kind, only : sp

        implicit none

        real(sp), intent(in) :: x                       ! time
        real(sp), dimension(:), intent(in) :: y         ! unknowns of the equations
        real(sp), dimension(:), intent(out) :: dydx     ! derivatives of y

        real(sp) :: xtmp

        xtmp = x

        dydx(1) = (-0.013_sp-1000.0_sp*y(3))*y(1)
        dydx(2) = -2500.0_sp*y(3)*y(2)
        dydx(3) = -0.013_sp*y(1) + (-1000.0_sp*y(1)-2500.0_sp*y(2))*y(3)

    end subroutine deriv2_testRB_sp

    subroutine jacobn2_testRB_sp( x, y, dfdx, dfdy )   ! returns derivatives dydx at x

        use mo_kind, only : sp

        implicit none

        real(sp), intent(in) :: x                           ! time
        real(sp), dimension(:),     intent(in)  :: y        ! unknowns of the equations
        real(sp), dimension(:),     intent(out) :: dfdx     ! derivatives of f for x
        real(sp), dimension(:,:),   intent(out) :: dfdy     ! derivatives of f for y

        real(sp) :: xtmp
        real(sp), dimension(size(y,1)) :: ytmp

        xtmp = x
        ytmp = y

        dfdx(1)     = 0.0_sp
        dfdx(2)     = 0.0_sp
        dfdx(3)     = 0.0_sp

        dfdy(1,1)   = -0.013_sp-1000.0_sp*y(3)
        dfdy(1,2)   = 0.0_sp
        dfdy(1,3)   = -1000.0_sp*y(1)
        dfdy(2,1)   = 0.0_sp
        dfdy(2,2)   = -2500.0_sp*y(3)
        dfdy(2,3)   = -2500.0_sp*y(2)
        dfdy(3,1)   = -0.013_sp-1000.0_sp*y(3)
        dfdy(3,2)   = -2500.0_sp*y(3)
        dfdy(3,3)   = -1000.0_sp*y(1)-2500.0_sp*y(2)

    end subroutine jacobn2_testRB_sp

    subroutine deriv2_testRB_para_sp( x, y, para, dydx )   ! returns derivatives dydx at x

        use mo_kind, only : sp

        implicit none

        real(sp), intent(in) :: x                       ! time
        real(sp), dimension(:), intent(in) :: y         ! unknowns of the equations
        real(sp), dimension(:), intent(out) :: dydx     ! derivatives of y

        real(sp), dimension(:), intent(in) :: para      ! parameter

        real(sp) :: xtmp

        xtmp = x


        dydx(1) = (-0.013_sp-para(1)*y(3))*y(1)
        dydx(2) = -2500.0_sp*y(3)*y(2)
        dydx(3) = -0.013_sp*y(1) + (-1000.0_sp*y(1)-2500.0_sp*y(2))*y(3)

    end subroutine deriv2_testRB_para_sp

    subroutine jacobn2_testRB_para_sp( x, y, para, dfdx, dfdy )   ! returns derivatives dydx at x

        use mo_kind, only : sp

        implicit none

        real(sp), intent(in) :: x                           ! time
        real(sp), dimension(:),     intent(in)  :: y        ! unknowns of the equations
        real(sp), dimension(:),     intent(out) :: dfdx     ! derivatives of f for x
        real(sp), dimension(:,:),   intent(out) :: dfdy     ! derivatives of f for y

        real(sp), dimension(:),     intent(in)  :: para      ! parameter

        real(sp) :: xtmp
        real(sp), dimension(size(y,1)) :: ytmp

        xtmp = x
        ytmp = y

        dfdx(1)     = 0.0_sp
        dfdx(2)     = 0.0_sp
        dfdx(3)     = 0.0_sp

        dfdy(1,1)   = -0.013_sp-para(1)*y(3)
        dfdy(1,2)   = 0.0_sp
        dfdy(1,3)   = -1000.0_sp*y(1)
        dfdy(2,1)   = 0.0_sp
        dfdy(2,2)   = -2500.0_sp*y(3)
        dfdy(2,3)   = -2500.0_sp*y(2)
        dfdy(3,1)   = -0.013_sp-1000.0_sp*y(3)
        dfdy(3,2)   = -2500.0_sp*y(3)
        dfdy(3,3)   = -1000.0_sp*y(1)-2500.0_sp*y(2)

    end subroutine jacobn2_testRB_para_sp

!DOUBLE PRECISION
    subroutine deriv_testRB_dp( x, y, dydx )   ! returns derivatives dydx at x

        use mo_kind, only : dp

        implicit none

        real(dp), intent(in) :: x                       ! time
        real(dp), dimension(:), intent(in) :: y         ! unknowns of the equations
        real(dp), dimension(:), intent(out) :: dydx     ! derivatives of y

        real(dp) :: xtmp

        xtmp = x

        dydx(1) = 998.0_dp*y(1) + 1998.0_dp*y(2)
        dydx(2) = -999.0_dp*y(1)- 1999.0_dp*y(2)

    end subroutine deriv_testRB_dp

    subroutine jacobn_testRB_dp( x, y, dfdx, dfdy )   ! returns derivatives dydx at x

        use mo_kind, only : dp

        implicit none

        real(dp), intent(in) :: x                           ! time
        real(dp), dimension(:),     intent(in)  :: y        ! unknowns of the equations
        real(dp), dimension(:),     intent(out) :: dfdx     ! derivatives of f for x
        real(dp), dimension(:,:),   intent(out) :: dfdy     ! derivatives of f for y

        real(dp) :: xtmp
        real(dp), dimension(size(y,1)) :: ytmp

        xtmp = x
        ytmp = y

        dfdx(1)     = 0.0_dp
        dfdx(2)     = 0.0_dp

        dfdy(1,1)   = 998.0_dp
        dfdy(1,2)   = 1998.0_dp
        dfdy(2,1)   = -999.0_dp
        dfdy(2,2)   = 1999.0_dp

    end subroutine jacobn_testRB_dp

    subroutine deriv2_testRB_dp( x, y, dydx )   ! returns derivatives dydx at x

        use mo_kind, only : dp

        implicit none

        real(dp), intent(in) :: x                       ! time
        real(dp), dimension(:), intent(in) :: y         ! unknowns of the equations
        real(dp), dimension(:), intent(out) :: dydx     ! derivatives of y

        real(dp) :: xtmp

        xtmp = x

        dydx(1) = (-0.013_dp-1000.0_dp*y(3))*y(1)
        dydx(2) = -2500.0_dp*y(3)*y(2)
        dydx(3) = -0.013_dp*y(1) + (-1000.0_dp*y(1)-2500.0_dp*y(2))*y(3)

    end subroutine deriv2_testRB_dp

    subroutine jacobn2_testRB_dp( x, y, dfdx, dfdy )   ! returns derivatives dydx at x

        use mo_kind, only : dp

        implicit none

        real(dp), intent(in) :: x                           ! time
        real(dp), dimension(:),     intent(in)  :: y        ! unknowns of the equations
        real(dp), dimension(:),     intent(out) :: dfdx     ! derivatives of f for x
        real(dp), dimension(:,:),   intent(out) :: dfdy     ! derivatives of f for y

        real(dp) :: xtmp
        real(dp), dimension(size(y,1)) :: ytmp

        xtmp = x
        ytmp = y

        dfdx(1)     = 0.0_dp
        dfdx(2)     = 0.0_dp
        dfdx(3)     = 0.0_dp

        dfdy(1,1)   = -0.013_dp-1000.0_dp*y(3)
        dfdy(1,2)   = 0.0_dp
        dfdy(1,3)   = -1000.0_dp*y(1)
        dfdy(2,1)   = 0.0_dp
        dfdy(2,2)   = -2500.0_dp*y(3)
        dfdy(2,3)   = -2500.0_dp*y(2)
        dfdy(3,1)   = -0.013_dp-1000.0_dp*y(3)
        dfdy(3,2)   = -2500.0_dp*y(3)
        dfdy(3,3)   = -1000.0_dp*y(1)-2500.0_dp*y(2)

    end subroutine jacobn2_testRB_dp

    subroutine deriv2_testRB_para_dp( x, y, para, dydx )   ! returns derivatives dydx at x

        use mo_kind, only : dp

        implicit none

        real(dp), intent(in) :: x                       ! time
        real(dp), dimension(:), intent(in) :: y         ! unknowns of the equations
        real(dp), dimension(:), intent(out) :: dydx     ! derivatives of y

        real(dp), dimension(:), intent(in) :: para      ! parameter

        real(dp) :: xtmp

        xtmp = x

        dydx(1) = (-0.013_dp-para(1)*y(3))*y(1)
        dydx(2) = -2500.0_dp*y(3)*y(2)
        dydx(3) = -0.013_dp*y(1) + (-1000.0_dp*y(1)-2500.0_dp*y(2))*y(3)

    end subroutine deriv2_testRB_para_dp

    subroutine jacobn2_testRB_para_dp( x, y, para, dfdx, dfdy )   ! returns derivatives dydx at x

        use mo_kind, only : dp

        implicit none

        real(dp), intent(in) :: x                           ! time
        real(dp), dimension(:),     intent(in)  :: y        ! unknowns of the equations
        real(dp), dimension(:),     intent(out) :: dfdx     ! derivatives of f for x
        real(dp), dimension(:,:),   intent(out) :: dfdy     ! derivatives of f for y

        real(dp), dimension(:),     intent(in)  :: para      ! parameter

        real(dp) :: xtmp
        real(dp), dimension(size(y,1)) :: ytmp

        xtmp = x
        ytmp = y

        dfdx(1)     = 0.0_dp
        dfdx(2)     = 0.0_dp
        dfdx(3)     = 0.0_dp

        dfdy(1,1)   = -0.013_dp-para(1)*y(3)
        dfdy(1,2)   = 0.0_dp
        dfdy(1,3)   = -1000.0_dp*y(1)
        dfdy(2,1)   = 0.0_dp
        dfdy(2,2)   = -2500.0_dp*y(3)
        dfdy(2,3)   = -2500.0_dp*y(2)
        dfdy(3,1)   = -0.013_dp-1000.0_dp*y(3)
        dfdy(3,2)   = -2500.0_dp*y(3)
        dfdy(3,3)   = -1000.0_dp*y(1)-2500.0_dp*y(2)

    end subroutine jacobn2_testRB_para_dp

end module mo_eq_lv

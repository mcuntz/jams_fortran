! ============================================================================
! Name        : mo_eq_lv.f90
! Author      : Giovanni Dalmasso
! Version     : 1.0
! Copyright   : Department of Computational Hydrosystems CHS - UZF Leipzig
! Description : Lotka Volterra equations
! ============================================================================

module mo_eq_lv

    ! ------------------------------------------------------------------
    ! The Lotka�Volterra equations, also known as the predator�prey equations, are a pair of first-order, non-linear,
    ! differential equations frequently used to describe the dynamics of biological systems in which two species interact,
    ! one a predator and one its prey. They evolve in time according to the following of equations
    ! ------------------------------------------------------------------

    implicit none

    private

    public :: LV_eqn, LV_eqn_sp, LV_eqn_dp

    ! Interfaces for single and double precision routines
    interface LV_eqn
        module procedure LV_eqn_sp, LV_eqn_dp
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
        ! ============================================================================

        dydx(1) = y(1) - epsilonn*mu*y(2)
        dydx(2) = (mu/epsilonn)*y(1) + eta*y(2)

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
        ! ============================================================================

        dydx(1) = y(1) - epsilonn*mu*y(2)
        dydx(2) = (mu/epsilonn)*y(1) + eta*y(2)

    end subroutine LV_eqn_dp

end module mo_eq_lv
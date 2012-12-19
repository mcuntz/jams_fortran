!> \file mo_constants.f90

!> \brief Provides computational, mathematical, physical, and file constants

!> \details Provides computational constants like epsilon, mathematical constants such as Pi,
!> physical constants such as the Stefan-Boltzmann constant, and file units for some standard streams
!> such as standard in.

!> \author Matthias Cuntz
!> \date Nov 2011
MODULE mo_constants

  !  This module contains basic and derived constants
  !
  !  Written Nov 2011, Matthias Cuntz

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

  ! Copyright 2011-2012 Matthias Cuntz

  USE mo_kind, ONLY: sp, dp

  IMPLICIT NONE

  ! Computational
  !> epsilon(1.0) in double precision
  REAL(dp), PARAMETER :: eps_dp = epsilon(1.0_dp)
  !> epsilon(1.0) in single precision
  REAL(sp), PARAMETER :: eps_sp = epsilon(1.0_sp)

  ! Mathematical
  !> Pi in double precision
  REAL(dp), PARAMETER :: PI_dp    = 3.141592653589793238462643383279502884197_dp    ! Pi
  !> Pi in single precision
  REAL(sp), PARAMETER :: PI_sp    = 3.141592653589793238462643383279502884197_sp
  !> Pi/2 in double precision
  REAL(dp), PARAMETER :: PIO2_dp  = 1.57079632679489661923132169163975144209858_dp  ! Pi/2
  !> Pi/2 in single precision
  REAL(sp), PARAMETER :: PIO2_sp  = 1.57079632679489661923132169163975144209858_sp
  !> 2*Pi in double precision
  REAL(dp), PARAMETER :: TWOPI_dp = 6.283185307179586476925286766559005768394_dp    ! 2*Pi
  !> 2*Pi in single precision
  REAL(sp), PARAMETER :: TWOPI_sp = 6.283185307179586476925286766559005768394_sp
  !> Square root of 2 in double precision
  REAL(dp), PARAMETER :: SQRT2_dp = 1.41421356237309504880168872420969807856967_dp  ! Sqrt(2)
  !> Square root of 2 in single precision
  REAL(sp), PARAMETER :: SQRT2_sp = 1.41421356237309504880168872420969807856967_sp
  !> 2/3 in double precision
  REAL(dp), PARAMETER :: TWOTHIRD_D = 0.6666666666666666666666666666666666667_dp      ! 2/3
  !> 2/3 in single precision
  REAL(sp), PARAMETER :: TWOTHIRD   = 0.6666666666666666666666666666666666667_sp

  ! Physical
  !> Psychrometric constant [kPa K^-1]
  REAL(sp), PARAMETER :: gamma_dp      = 0.0646_dp                   ! psychrometric constant [kPa C-1]
  !> Gravity accelaration [m^2 s^-1] in double precision
  REAL(dp), PARAMETER :: Gravity_dp   = 9.81_dp                      ! Gravity acceleration [m^2/s]
  !> Gravity accelaration [m^2 s^-1] in single precision
  REAL(sp), PARAMETER :: Gravity_sp   = 9.81_sp
  !>  Solar constant in [MJ m^-2 s^-1]
  REAL(sp), PARAMETER :: GSolarConst_dp  = 1367_dp                   ! Global solar constant in [W m-2 = kg s-3]
  !> Specific heat for vaporization of water in [J m-2 mm-1]
  REAL(sp), PARAMETER :: SpecHeatET_dp   = 2.45e06_dp                ! in [W s m-2 mm-1 = kg s-2 mm-1]
  !> Standard temperature [K] in double precision
  REAL(dp), PARAMETER :: T0_dp        = 273.15_dp                    ! Celcius <-> Kelvin [K]
  !> Standard temperature [K] in single precision
  REAL(sp), PARAMETER :: T0_sp        = 273.15_sp
  !> Standard pressure [Pa] in double precision
  REAL(dp), PARAMETER :: P0_dp        = 101325._dp                   ! Standard pressure [Pa]
  !> Standard pressure [Pa] in single precision
  REAL(sp), PARAMETER :: P0_sp        = 101325._sp
  !> Stefan-Boltzmann constant [W m^-2 K^-4] in double precision
  REAL(dp), PARAMETER :: sigma_dp     = 5.67e-08_dp                  ! Stefan-Boltzmann constant [W/m^2/K^4]
  !> Stefan-Boltzmann constant [W m^-2 K^-4] in single precision
  REAL(sp), PARAMETER :: sigma_sp     = 5.67e-08_sp

  ! Numerical Recipes
  !> Pi in double precision
  REAL(dp), PARAMETER :: PI_D    = 3.141592653589793238462643383279502884197_dp      ! Pi
  !> Pi in single precision
  REAL(sp), PARAMETER :: PI      = 3.141592653589793238462643383279502884197_sp
  !> Pi/2 in double precision
  REAL(dp), PARAMETER :: PIO2_D  = 1.57079632679489661923132169163975144209858_dp    ! Pi/2
  !> Pi/2 in single precision
  REAL(sp), PARAMETER :: PIO2    = 1.57079632679489661923132169163975144209858_sp
  !> 2*Pi in double precision
  REAL(dp), PARAMETER :: TWOPI_D = 6.283185307179586476925286766559005768394_dp      ! 2*Pi
  !> 2*Pi in single precision
  REAL(sp), PARAMETER :: TWOPI   = 6.283185307179586476925286766559005768394_sp
  !> Square root of 2 in double precision
  REAL(dp), PARAMETER :: SQRT2_D = 1.41421356237309504880168872420969807856967_dp    ! Sqrt(2)
  !> Square root of 2 in single precision
  REAL(sp), PARAMETER :: SQRT2   = 1.41421356237309504880168872420969807856967_sp
  !> Euler''s constant in double precision
  REAL(dp), PARAMETER :: EULER_D = 0.5772156649015328606065120900824024310422_dp     ! Euler
  !> Euler''s constant in single precision
  REAL(sp), PARAMETER :: EULER   = 0.5772156649015328606065120900824024310422_sp

  ! file units
  !> Standard input file unit
  INTEGER, PARAMETER :: nin  = 5   ! standard input stream
  !> Standard output file unit
  INTEGER, PARAMETER :: nout = 6   ! standard output stream
  !> Standard error file unit
  INTEGER, PARAMETER :: nerr = 0   ! error output stream
  !> Standard file unit for namelist
  INTEGER, PARAMETER :: nnml = 100 ! namelist unit

END MODULE mo_constants

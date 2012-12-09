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

  ! Mathematical
  ! Pi
  REAL(dp), PARAMETER :: PI_dp    = 3.141592653589793238462643383279502884197_dp
  REAL(sp), PARAMETER :: PI_sp    = 3.141592653589793238462643383279502884197_sp
  ! Pi/2
  REAL(dp), PARAMETER :: PIO2_dp  = 1.57079632679489661923132169163975144209858_dp
  REAL(sp), PARAMETER :: PIO2_sp  = 1.57079632679489661923132169163975144209858_sp
  ! 2*Pi
  REAL(dp), PARAMETER :: TWOPI_dp = 6.283185307179586476925286766559005768394_dp
  REAL(sp), PARAMETER :: TWOPI_sp = 6.283185307179586476925286766559005768394_sp
  ! Sqrt(2)
  REAL(dp), PARAMETER :: SQRT2_dp = 1.41421356237309504880168872420969807856967_dp
  REAL(sp), PARAMETER :: SQRT2_sp = 1.41421356237309504880168872420969807856967_sp

  ! Physical
  ! Gravity acceleration [m^2/s]
  REAL(dp), PARAMETER :: Gravity_dp   = 9.81_dp
  REAL(sp), PARAMETER :: Gravity_sp   = 9.81_sp
  ! Celcius <-> Kelvin [K]
  REAL(dp), PARAMETER :: T0_dp        = 273.15_dp
  REAL(sp), PARAMETER :: T0_sp        = 273.15_sp
  ! Standard pressure [Pa]
  REAL(dp), PARAMETER :: P0_dp        = 101325._dp
  REAL(sp), PARAMETER :: P0_sp        = 101325._sp
  ! Stefan-Boltzmann constant [W/m^2/K^4]
  REAL(dp), PARAMETER :: sigma_dp     = 5.67e-08_dp
  REAL(sp), PARAMETER :: sigma_sp     = 5.67e-08_sp

  ! Numerical Recipes
  ! Pi
  REAL(dp), PARAMETER :: PI_D    = 3.141592653589793238462643383279502884197_dp
  REAL(sp), PARAMETER :: PI      = 3.141592653589793238462643383279502884197_sp
  ! Pi/2
  REAL(dp), PARAMETER :: PIO2_D  = 1.57079632679489661923132169163975144209858_dp
  REAL(sp), PARAMETER :: PIO2    = 1.57079632679489661923132169163975144209858_sp
  ! 2*Pi
  REAL(dp), PARAMETER :: TWOPI_D = 6.283185307179586476925286766559005768394_dp
  REAL(sp), PARAMETER :: TWOPI   = 6.283185307179586476925286766559005768394_sp
  ! Sqrt(2)
  REAL(dp), PARAMETER :: SQRT2_D = 1.41421356237309504880168872420969807856967_dp
  REAL(sp), PARAMETER :: SQRT2   = 1.41421356237309504880168872420969807856967_sp
  ! ???
  REAL(dp), PARAMETER :: EULER_D = 0.5772156649015328606065120900824024310422_dp
  REAL(sp), PARAMETER :: EULER   = 0.5772156649015328606065120900824024310422_sp

  ! File units
  ! standard input stream
  INTEGER, PARAMETER :: nin  = 5
  ! standard output stream
  INTEGER, PARAMETER :: nout = 6
  ! error output stream
  INTEGER, PARAMETER :: nerr = 0
  ! namelist unit
  INTEGER, PARAMETER :: nnml = 100

END MODULE mo_constants

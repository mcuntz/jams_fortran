!> \file mo_constants.f90

!> \brief Provides computational, mathematical, physical, and file constants

!> \details Provides computational constants like epsilon, mathematical constants such as Pi,
!> physical constants such as the Stefan-Boltzmann constant, and file units for some standard streams
!> such as standard in.

!> \author Matthias Cuntz
!> \date Nov 2011

module mo_constants

  !  This module contains basic and derived constants
  !
  !  Written  Nov 2011, Matthias Cuntz
  !  Modified Mar 2014, Matthias Cuntz                   - iso_fortran_env
  !           Jun 2018, Matthias Cuntz, Johannes Brenner - constants for standard atmosphere (mtclim)

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2011-2018 Matthias Cuntz - mc (at) macu (dot) de
  !
  ! Permission is hereby granted, free of charge, to any person obtaining a copy
  ! of this software and associated documentation files (the "Software"), to deal
  ! in the Software without restriction, including without limitation the rights
  ! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  ! copies of the Software, and to permit persons to whom the Software is
  ! furnished to do so, subject to the following conditions:
  !
  ! The above copyright notice and this permission notice shall be included in all
  ! copies or substantial portions of the Software.
  !
  ! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  ! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  ! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  ! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  ! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  ! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  ! SOFTWARE.

  use mo_kind, only: sp, dp
  use, intrinsic :: iso_fortran_env, only: input_unit, output_unit, error_unit

  implicit none

  ! Computational
  !> epsilon(1.0) in double precision
  real(dp), parameter :: eps_dp = epsilon(1.0_dp)
  !> epsilon(1.0) in single precision
  real(sp), parameter :: eps_sp = epsilon(1.0_sp)

  
  ! Mathematical
  !> pi in double precision
  real(dp), parameter :: pi_dp        = 3.141592653589793238462643383279502884197_dp    ! pi
  !> pi in single precision
  real(sp), parameter :: pi_sp        = 3.141592653589793238462643383279502884197_sp
  !> pi/2 in double precision
  real(dp), parameter :: pio2_dp      = 1.57079632679489661923132169163975144209858_dp  ! pi/2
  !> pi/2 in single precision
  real(sp), parameter :: pio2_sp      = 1.57079632679489661923132169163975144209858_sp
  !> 2*pi in double precision
  real(dp), parameter :: twopi_dp     = 6.283185307179586476925286766559005768394_dp    ! 2*pi
  !> 2*pi in single precision
  real(sp), parameter :: twopi_sp     = 6.283185307179586476925286766559005768394_sp
  !> square root of 2 in double precision
  real(dp), parameter :: sqrt2_dp     = 1.41421356237309504880168872420969807856967_dp  ! sqrt(2)
  !> square root of 2 in single precision
  real(sp), parameter :: sqrt2_sp     = 1.41421356237309504880168872420969807856967_sp
  !> 1/3 in double precision
  real(dp), parameter :: onethird_dp  = 0.3333333333333333333333333333333333333_dp      ! 1/3
  !> 1/3 in single precision
  real(sp), parameter :: onethird_sp  = 0.3333333333333333333333333333333333333_sp
  !> 2/3 in double precision
  real(dp), parameter :: twothird_dp  = 1.0_dp - onethird_dp                            ! 2/3
  !> 2/3 in single precision
  real(sp), parameter :: twothird_sp  = 1.0_sp - onethird_sp
  !> degree to radian conversion (pi/180) in double precision
  real(dp), parameter :: deg2rad_dp   = pi_dp/180._dp            ! deg2rad
  !> degree to radian conversion (pi/180) in single precision
  real(sp), parameter :: deg2rad_sp   = pi_sp/180._sp
  !> radian to degree conversion (180/pi) in double precision
  real(dp), parameter :: rad2deg_dp   = 180._dp/pi_dp            ! rad2deg
  !> radian to degree conversion (180/pi) in single precision
  real(sp), parameter :: rad2deg_sp   = 180._sp/pi_sp
  !> arc to radian conversion (180/pi) in double precision
  real(dp), parameter :: arc2rad_dp   = pi_dp/12._dp             ! arc2rad
  !> arc to radian conversion (180/pi) in single precision
  real(sp), parameter :: arc2rad_sp   = pi_sp/12._sp
  ! seconds per radian of hour angle in single precision
  real(dp), parameter :: secperrad_dp = 13750.9871_dp            ! secperrad
  ! seconds per radian of hour angle in single precision
  real(sp), parameter :: secperrad_sp = 13750.9871_sp
  

  ! Physical
  !> seconds per day [s] in double precision
  real(dp), parameter :: secday_dp      = 86400._dp               ! secday [s]
  !> seconds per day [s] in single precision
  real(sp), parameter :: secday_sp      = 86400._sp
  !> psychrometric constant [kPa K^-1] in double precision
  real(dp), parameter :: psychro_dp     = 0.0646_dp               ! psychrometric constant [kPa c-1]
  !> psychrometric constant [kPa K^-1] in sibgle precision
  real(sp), parameter :: psychro_sp     = 0.0646_sp
  !> gravity accelaration [m^2 s^-1] in double precision
  real(dp), parameter :: gravity_dp     = 9.81_dp                 ! gravity acceleration [m^2/s]
  !> gravity accelaration [m^2 s^-1] in single precision
  real(sp), parameter :: gravity_sp     = 9.81_sp
  !>  solar constant in [J m^-2 s^-1] in double precision
  real(dp), parameter :: solarconst_dp  = 1367._dp                ! solar constant in [W m-2 = kg s-3]
  !>  solar constant in [J m^-2 s^-1] in single precision
  real(sp), parameter :: solarconst_sp  = 1367._sp
  !> specific heat for vaporization of water in [J m-2 mm-1] in double precision
  real(dp), parameter :: specheatet_dp  = 2.45e06_dp              ! specific heat in [W s m-2 mm-1 = kg s-2 mm-1]
  !> specific heat for vaporization of water in [J m-2 mm-1] in single precision
  real(sp), parameter :: specheatet_sp  = 2.45e06_sp
  !> standard temperature [k] in double precision
  real(dp), parameter :: T0_dp          = 273.15_dp               ! celcius <-> kelvin [K]
  !> standard temperature [k] in single precision
  real(sp), parameter :: T0_sp          = 273.15_sp
  !> stefan-boltzmann constant [W m^-2 K^-4] in double precision
  real(dp), parameter :: sigma_dp       = 5.67e-08_dp             ! stefan-boltzmann constant [w/m^2/k^4]
  !> stefan-boltzmann constant [W m^-2 K^-4] in single precision
  real(sp), parameter :: sigma_sp       = 5.67e-08_sp
  ! earth radius [m] in double precision
  real(dp), parameter :: radiusearth_dp = 6371228._dp
  ! earth radius [m] in single precision
  real(sp), parameter :: radiusearth_sp = 6371228._sp
  ! radians of earth orbit per julian day
  real(dp), parameter :: radperday_dp   = 0.017214_dp
  ! radians of earth orbit per julian day
  real(sp), parameter :: radperday_sp   = 0.017214_sp
  ! minimum declination (radians)
  real(dp), parameter :: mindecl_dp     = -0.4092797_dp
  ! minimum declination (radians)
  real(sp), parameter :: mindecl_sp     = -0.4092797_sp
  ! julian day offset of winter solstice
  real(dp), parameter :: daysoff_dp     = 11.25_dp
  ! julian day offset of winter solstice
  real(sp), parameter :: daysoff_sp     = 11.25_sp
  ! molecular weight of air (kg mol-1)
  real(dp), parameter :: Ma_dp          = 28.9644e-3_dp
  ! molecular weight of air (kg mol-1)
  real(sp), parameter :: Ma_sp          = 28.9644e-3_sp
  ! molecular weight of water (kg mol-1)
  real(dp), parameter :: Mw_dp          = 18.0148e-3_dp
  ! molecular weight of water (kg mol-1)
  real(sp), parameter :: Mw_sp          = 18.0148e-3_sp
  ! unitless ratio of molec weights (mw/ma)
  real(dp), parameter :: epsair_dp      = Mw_dp/Ma_dp
  ! unitless ratio of molec weights (mw/ma)
  real(sp), parameter :: epsair_sp      = Mw_sp/Ma_sp
  ! gas law constant (m3 pa mol-1 K-1)
  real(dp), parameter :: R_dp           = 8.3143_dp
  ! gas law constant (m3 pa mol-1 K-1)
  real(sp), parameter :: R_sp           = 8.3143_sp
  ! standard temperature lapse rate (-K m-1)
  real(dp), parameter :: lr_std_dp      = 0.0065_dp
  ! standard temperature lapse rate (-K m-1)
  real(sp), parameter :: lr_std_sp      = 0.0065_sp

  !> Standard atmosphere
  !> standard pressure [Pa] in double precision
  real(dp), parameter :: p0_dp    = 101325._dp      ! standard pressure [Pa]
  !> standard pressure [Pa] in single precision
  real(sp), parameter :: p0_sp    = 101325._sp
  !> standard density  [kg m^-3] in double precision
  real(dp), parameter :: rho0_dp  = 1.225_dp        ! standard air density
  !> standard density  [kg m^-3] in single precision
  real(sp), parameter :: rho0_sp  = 1.225_sp
  !> specific heat capacity of air [J kg^-1 K^-1] in double precision
  real(dp), parameter :: cp0_dp   = 1005.0_dp       ! standard specific heat of air
  !> specific heat capacity of air [J kg^-1 K^-1] in single precision
  real(sp), parameter :: cp0_sp   = 1005.0_sp
  !> standard temp at 0.0 m elevation (K)           ! standard temperature at 0 m asl
  real(dp), parameter :: T_std_dp = 288.15_dp
  !> standard temp at 0.0 m elevation (K)
  real(sp), parameter :: T_std_sp = 288.15_sp

  
  ! Numerical Recipes
  !> pi in double precision
  real(dp), parameter :: pi_d    = 3.141592653589793238462643383279502884197_dp      ! pi
  !> pi in single precision
  real(sp), parameter :: pi      = 3.141592653589793238462643383279502884197_sp
  !> pi/2 in double precision
  real(dp), parameter :: pio2_d  = 1.57079632679489661923132169163975144209858_dp    ! pi/2
  !> pi/2 in single precision
  real(sp), parameter :: pio2    = 1.57079632679489661923132169163975144209858_sp
  !> 2*pi in double precision
  real(dp), parameter :: twopi_d = 6.283185307179586476925286766559005768394_dp      ! 2*pi
  !> 2*pi in single precision
  real(sp), parameter :: twopi   = 6.283185307179586476925286766559005768394_sp
  !> square root of 2 in double precision
  real(dp), parameter :: sqrt2_d = 1.41421356237309504880168872420969807856967_dp    ! sqrt(2)
  !> square root of 2 in single precision
  real(sp), parameter :: sqrt2   = 1.41421356237309504880168872420969807856967_sp
  !> Euler's constant in double precision
  real(dp), parameter :: euler_d = 0.5772156649015328606065120900824024310422_dp     ! Euler
  !> Euler's constant in single precision
  real(sp), parameter :: euler   = 0.5772156649015328606065120900824024310422_sp

  ! Standard file units
  !> standard input file unit
  ! integer, parameter :: nin  = 5   ! standard input stream
  integer, parameter :: nin  = input_unit   ! standard input stream
  !> standard output file unit
  ! integer, parameter :: nout = 6   ! standard output stream
  integer, parameter :: nout = output_unit   ! standard output stream
  !> standard error file unit
  ! integer, parameter :: nerr = 0   ! error output stream
  integer, parameter :: nerr = error_unit   ! error output stream
  !> standard file unit for namelist
  integer, parameter :: nnml = 100 ! namelist unit

end module mo_constants

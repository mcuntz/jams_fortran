!> \file mo_isotope.f90

!> \brief Isotope parameters and helper routines

!> \details Provides isotope parameters such as international standards
!> and convenience routines such as delta value calculations.

!> \author Matthias Cuntz
!> \date Apr 2019

! License
! -------
! This file is part of the JAMS Fortran package, distributed under the MIT License.
!
! Copyright (c) 2019 Matthias Cuntz - mc (at) macu (dot) de
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

MODULE mo_isotope

  use mo_kind,      only: dp
  use mo_constants, only: Mw_dp, Ma_dp, twothird_dp
  
  implicit none

  private

  ! Public routines
  public :: delta     ! isotopic delta value
  public :: delta1000 ! isotopic delta value in permil
  public :: isoratio  ! ratio of rare to abundant isotope

  ! Public parameters
  ! ratio diffusion of air vs. water vapour ~ 1.6
  real(dp), parameter, public :: ratio_diff_air2vap = Ma_dp / Mw_dp
  ! ratio of "diffusion" of air vs. water vapour through boundary layer ~ 1.4
  real(dp), parameter, public :: ratio_boundary_air2vap = ratio_diff_air2vap**twothird_dp
  ! ratio diffusion of water vapour vs. air ~ 0.62
  real(dp), parameter, public :: ratio_diff_vap2air = Mw_dp / Ma_dp
  ! ratio of "diffusion" of water vapour vs. air through boundary layer ~ 0.73
  real(dp), parameter, public :: ratio_boundary_vap2air = ratio_diff_vap2air**twothird_dp
  ! 13C/12C VPDB standard
  real(dp), parameter, public :: vpdbc13 = 0.0112372_dp
  ! molecular weight of 13C (kg mol-1)
  real(dp), parameter, public :: Mc13   = 13.e-3_dp
  ! molecular weight of 13CO2 (kg mol-1)
  real(dp), parameter, public :: Mc13o2 = 45.e-3_dp

contains

  ! ------------------------------------------------------------------

  !     NAME
  !         delta

  !     PURPOSE
  !>        \brief Isotopic delta value.

  !>        \details Elemental routine to calculate the delta value
  !>        from an isotope ratio or the rare and the abundant isotope concentrations
  !>        relative to a given standard.
  
  !>        The delta value can be set to a default value if the abundant
  !>        isotopic concentration is below a give precision.

  !     CALLING SEQUENCE
  !         del = delta(rare, abundant, standard, default, precision)

  !     PARAMETER
  !>        \param[in] "real(dp) :: rare"                   Rare isotope concentration or isotope ratio
  !>        \param[in] "real(dp), optional :: abundant"     Abundant isotope concentration if rare is rare isotope concentration
  !>                                                        Default: not used.
  !>        \param[in] "real(dp), optional :: standard"     Isotop ratio of the standard material.
  !>                                                        Default: 1.
  !>        \param[in] "real(dp), optional :: default"      Set isoratio default if abundant < precision.
  !>                                                        Default: 0.
  !>        \param[in] "real(dp), optional :: precision"    Set isoratio default if abundant < precision.

  !     RETURN
  !>       \return      real(dp) :: delta &mdash; Isotopic delta value: rare / abundant / standard - 1.

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  elemental pure function delta(rare, abundant, standard, default, precision)

    implicit none

    real(dp), intent(in)           :: rare
    real(dp), intent(in), optional :: abundant
    real(dp), intent(in), optional :: standard
    real(dp), intent(in), optional :: default
    real(dp), intent(in), optional :: precision
    real(dp)                       :: delta

    real(dp) :: istandard, idefault, iprecision

    ! default optionals
    if (present(standard)) then
       istandard = standard
    else
       istandard = 0.0_dp
    endif
    if (present(default)) then
       idefault = default
    else
       idefault = 0.0_dp
    endif
    if (present(precision)) then
       iprecision = precision
    else
       iprecision = 0.0_dp
    endif

    if (present(abundant)) then
       if (abs(abundant) > iprecision) then
          delta = rare / abundant / istandard - 1.0_dp
       else
          delta = idefault
       end if
    else
       delta = rare / istandard - 1.0_dp
    end if
    
  end function delta

  
  ! ------------------------------------------------------------------

  !     NAME
  !         delta1000

  !     PURPOSE
  !>        \brief Isotopic delta value in permil.

  !>        \details Elemental routine to calculate the delta value in permil
  !>        from an isotope ratio or the rare and the abundant isotope concentrations
  !>        relative to a given standard.
  
  !>        The delta value can be set to a default value if the abundant
  !>        isotopic concentration is below a give precision.

  !     CALLING SEQUENCE
  !         del = delta1000(rare, abundant, standard, default, precision)

  !     PARAMETER
  !>        \param[in] "real(dp) :: rare"                   Rare isotope concentration or isotope ratio
  !>        \param[in] "real(dp), optional :: abundant"     Abundant isotope concentration if rare is rare isotope concentration
  !>                                                        Default: not used.
  !>        \param[in] "real(dp), optional :: standard"     Isotop ratio of the standard material.
  !>                                                        Default: 1.
  !>        \param[in] "real(dp), optional :: default"      Set isoratio default if abundant < precision.
  !>                                                        Default: 0.
  !>        \param[in] "real(dp), optional :: precision"    Set isoratio default if abundant < precision.

  !     RETURN
  !>       \return      real(dp) :: delta &mdash; Isotopic delta value in permil: (rare / abundant / standard - 1.) * 1000.

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  elemental pure function delta1000(rare, abundant, standard, default, precision)

    implicit none

    real(dp), intent(in)           :: rare
    real(dp), intent(in), optional :: abundant
    real(dp), intent(in), optional :: standard
    real(dp), intent(in), optional :: default
    real(dp), intent(in), optional :: precision
    real(dp)                       :: delta1000

    real(dp) :: istandard, idefault, iprecision

    ! default optionals
    if (present(standard)) then
       istandard = standard
    else
       istandard = 0.0_dp
    endif
    if (present(default)) then
       idefault = default
    else
       idefault = 0.0_dp
    endif
    if (present(precision)) then
       iprecision = precision
    else
       iprecision = 0.0_dp
    endif

    if (present(abundant)) then
       if (abs(abundant) > iprecision) then
          delta1000 = (rare / abundant / istandard - 1.0_dp) * 1000._dp
       else
          delta1000 = idefault
       end if
    else
       delta1000 = (rare / istandard - 1.0_dp) * 1000._dp
    end if
    
  end function delta1000

  
  ! ------------------------------------------------------------------

  !     NAME
  !         isoratio

  !     PURPOSE
  !>        \brief Ratio of rare to abundant isotopes.

  !>        \details Elemental routine to calculate the ratio of the rare
  !>        and the abundant isotope concentrations.

  !>        The ratio can be set to a default value if the abundant
  !>        isotopic concentration is below a give precision.

  !     CALLING SEQUENCE
  !         ratio = isoratio(rare, abundant, default, precision)

  !     PARAMETER
  !>        \param[in] "real(dp) :: rare"                   Rare isotope concentration
  !>        \param[in] "real(dp) :: abundant"               Abundant isotope concentration
  !>        \param[in] "real(dp), optional :: default"      Set isoratio default if abundant < precision.
  !>                                                        Default: 0.
  !>        \param[in] "real(dp), optional :: precision"    Set isoratio default if abundant < precision.

  !     RETURN
  !>       \return      real(dp) :: isoratio &mdash; Isotope ratio: rare / abundant

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  elemental pure function isoratio(rare, abundant, default, precision)

    implicit none

    real(dp), intent(in)           :: rare
    real(dp), intent(in)           :: abundant
    real(dp), intent(in), optional :: default
    real(dp), intent(in), optional :: precision
    real(dp)                       :: isoratio

    real(dp) :: idefault, iprecision

    ! default optionals
    if (present(default)) then
       idefault = default
    else
       idefault = 0.0_dp
    endif
    if (present(precision)) then
       iprecision = precision
    else
       iprecision = 0.0_dp
    endif
    
    if (abs(abundant) > iprecision) then
       isoratio = rare / abundant
    else
       isoratio = idefault
    end if
    
  end function isoratio

END MODULE mo_isotope

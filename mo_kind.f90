!> \file mo_kind.f90

!> \brief Define number representations

!> \details This module declares the desired ranges and precisions of the number representations,
!> such as single precision or double precision, 32-bit or 64-bit integer, etc.
!> It is similar to the nrtype module of Numerical Recipes in f90.

!> \authors Matthias Cuntz, Juliane Mai, Nov 2011
!> \date 2011-2017

!  Number model from which the SELECTED_REAL_KIND are requested:
!                   4 byte REAL      8 byte REAL        16 byte REAL
!          IEEE:    precision =  6   precision =   15   precision =   33
!                   exponent  = 37   exponent  =  307   exponent  = 4931
!          CRAY:        -            precision =   13
!                                    exponent  = 2465

! Written  Juliane Mai, Matthias Cuntz, Nov 2011
! Modified Matthias Cuntz, Nov 2011 - private/public
!                                   - documentation
!          Matthias Cuntz, May 2014 - iso_fortran_env and iso_c_binding
!          Matthias Cuntz, Dec 2017 - quadruple precision
!          Matthias Cuntz, Mar 2020 - i1 = int8_t, which is guaranteed to be 8 bit

! License
! -------
! This file is part of the JAMS Fortran package, distributed under the MIT License.
!
! Copyright (c) 2011-2020 Matthias Cuntz, Juliane Mai - mc (at) macu (dot) de
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

MODULE mo_kind

  ! iso_fortran_env does not work with compilers intel v11 and sun v12.2
  ! use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, real32, real64, real128
  use, intrinsic :: iso_c_binding,   only: &
       c_int8_t, c_short, c_int, c_long_long, &
       c_float, c_double, c_long_double, &
       c_float_complex, c_double_complex, c_long_double_complex !, &
       ! c_bool

  IMPLICIT NONE

  !> 1 Byte Integer Kind
  INTEGER, PARAMETER :: i1  = c_int8_t ! SELECTED_INT_KIND(2)
  ! INTEGER, PARAMETER :: i1  = int8 ! c_word does not exist; should be c_bool, probably; but see above
  !> 2 Byte Integer Kind
  ! INTEGER, PARAMETER :: i2  = SELECTED_INT_KIND(4)
  ! INTEGER, PARAMETER :: i2  = int16
  INTEGER, PARAMETER :: i2  = c_short
  !> 4 Byte Integer Kind
  ! INTEGER, PARAMETER :: i4  = SELECTED_INT_KIND(9)
  ! INTEGER, PARAMETER :: i4  = int32
  INTEGER, PARAMETER :: i4  = c_int
  !> 8 Byte Integer Kind
  ! INTEGER, PARAMETER :: i8  = SELECTED_INT_KIND(18)
  ! INTEGER, PARAMETER :: i8  = int64
  ! INTEGER, PARAMETER :: i8  = c_long
  INTEGER, PARAMETER :: i8  = c_long_long
  !> Single Precision Real Kind
  ! INTEGER, PARAMETER :: sp  = SELECTED_REAL_KIND(6,37)
  ! INTEGER, PARAMETER :: sp  = real32
  INTEGER, PARAMETER :: sp  = c_float
  !> Double Precision Real Kind
  ! INTEGER, PARAMETER :: dp  = SELECTED_REAL_KIND(15,307)
  ! INTEGER, PARAMETER :: dp  = real64
  INTEGER, PARAMETER :: dp  = c_double
  !> Quadruple Precision Real Kind
  ! INTEGER, PARAMETER :: qp  = SELECTED_REAL_KIND(33,4931)
  ! INTEGER, PARAMETER :: qp  = real128
  INTEGER, PARAMETER :: qp  = c_long_double
  !> Single Precision Complex Kind
  ! INTEGER, PARAMETER :: spc = KIND((1.0_sp,1.0_sp))
  ! INTEGER, PARAMETER :: spc = sp
  INTEGER, PARAMETER :: spc = c_float_complex
  !> Double Precision Complex Kind
  ! INTEGER, PARAMETER :: dpc = KIND((1.0_dp,1.0_dp))
  ! INTEGER, PARAMETER :: dpc = dp
  INTEGER, PARAMETER :: dpc = c_double_complex
  !> Quadruple Precision Complex Kind
  ! INTEGER, PARAMETER :: qpc = KIND((1.0_dpc,1.0_dpc))
  ! INTEGER, PARAMETER :: qpc = qp
  INTEGER, PARAMETER :: qpc  = c_long_double_complex
  !> Logical Kind
  INTEGER, PARAMETER :: lgt = KIND(.true.)
  !INTEGER, PARAMETER :: lgt = c_bool

  ! Types have to be in a public section for doxygen
  !> Single Precision Numerical Recipes types for sparse arrays
  TYPE sprs2_sp
     INTEGER(i4) :: n, len
     REAL(sp),    DIMENSION(:), POINTER :: val
     INTEGER(i4), DIMENSION(:), POINTER :: irow
     INTEGER(i4), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_sp
  !> Double Precision Numerical Recipes types for sparse arrays
  TYPE sprs2_dp
     INTEGER(i4) :: n, len
     REAL(dp),    DIMENSION(:), POINTER :: val
     INTEGER(i4), DIMENSION(:), POINTER :: irow
     INTEGER(i4), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_dp

END MODULE mo_kind

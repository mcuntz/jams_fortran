!> \file mo_ranper.f90

!> \brief calculate random permutation

!> \details The code written here is originally based on the subroutine
!> ranper presented in "Combinatorial algorithms for computers and 
!> calculators (2nd ed.)", (Nijenhuis, A., & Wilf, H. S. (1978), Academic
!> Pr.) on page 63. It was modified for compatability with our Fortran
!> Library

!> \authors Stephan Thober
!> \date Feb 2013
module mo_ranper
  
  ! This module contains subroutine for calculating random permutations.

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

  ! Copyright 2011-2012 Matthias Cuntz, Juliane Mai
  !
  implicit none
  !
  private
  !
  public :: ranper
  !
    ! ------------------------------------------------------------------

  !     NAME
  !         ranper

  !     PURPOSE
  !         calculates a random permutation of n numbers
  !
  !>        \brief random permutation.
  !
  !>        \details calculates a random permutation of n integer numbers
  !>        \f$ x_i, i =1,\ldots,n \f$
  !
  !>        If the boolean setup is set true, then the array to permute
  !>        is initialized from 1 to n, otherwise it is taken from input.
  !
  !     INTENT(IN)
  !>        \param[in]    "logical        :: setup"    When set .true. array A is 
  !>        initialized from \f$ 1,\ldots,n \f, otherwise taken from input
  !
  !     INTENT(INOUT)
  !>        \param[inout] "integer(i4,dp) :: A(:)"   \f$ x_i \f$ 1D-array with input numbers

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>        \param[in]    "logical        :: init_xor" When present and set .true., the random
  !>        number generator is reset with a new seed
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !>       \note If setup is .false., then A must contain positive numbers.
  !>       ranper uses the random number generator xor4096 for generating
  !>       i4 integer random numbers. If this random number generator is used
  !>       elsewhere in the code, then the parameters of the random number
  !>       generator have to be stored to restart the actual random number
  !>       stream, since ranper requires its own stream.
  !
  !     EXAMPLE
  !         A = (/ 5, 2, 3, 41, 4, 6 /)
  !         call ranper( A, .false.)
  !         A = (/ 2, 41, 6, 3, 5, 4 /)
  !         -> see also example in test directory

  !     LITERATURE
  !         Nijenhuis, A., & Wilf, H. S. (1978). Combinatorial algorithms for computers and calculators (2nd ed.). Academic Pr., p. 63

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Feb 2013
  !
  INTERFACE ranper
     MODULE PROCEDURE ranper_i4, ranper_i8
  END INTERFACE ranper
  !
contains
  subroutine ranper_i4( A, Setup, init_xor )
    !
    use mo_kind,    only: i4, sp
    use mo_xor4096, only: get_timeseed, xor4096
    !
    implicit none
    !
    logical,     optional,     intent(in)    :: init_xor
    logical,                   intent(in)    :: setup
    integer(i4), dimension(:), intent(inout) :: A
    !
    ! local variables
    logical     :: init
    integer(i4) :: seed
    integer(i4) :: i
    integer(i4) :: m
    integer(i4) :: l
    integer(i4) :: L1
    real(sp)    :: rn ! random number
    !
    init = .false.
    !
    if ( present( init_xor ) ) init = init_xor
    !
    if ( init ) then
       call get_timeseed( seed )
       call xor4096( seed, rn )
    end if
    !
    if ( setup ) then
       do i = 1, size(A)
          A(i) = i
       end do
    end if
    !
    do m = 1, size(A)
       call xor4096( 0_i4, rn )
       l = m + int(rn * real( size(A) + 1 - m, sp ), i4)
       L1= A(l)
       A(l) = A(m)
       A(m) = L1
    end do
    !
  end subroutine ranper_i4
  !
  subroutine ranper_i8( A, Setup, init_xor )
    !
    use mo_kind,    only: i8, dp
    use mo_xor4096, only: get_timeseed, xor4096
    !
    implicit none
    !
    logical,     optional,     intent(in)    :: init_xor
    logical,                   intent(in)    :: setup
    integer(i8), dimension(:), intent(inout) :: A
    !
    ! local variables
    logical     :: init
    integer(i8) :: seed
    integer(i8) :: i
    integer(i8) :: m
    integer(i8) :: l
    integer(i8) :: L1
    real(dp)    :: rn ! random number
    !
    init = .false.
    if ( present( init_xor ) ) init = init_xor
    !
    if ( init ) then
       call get_timeseed( seed )
       call xor4096( seed, rn )
    end if
    !
    if ( setup ) then
       do i = 1, size(A)
          A(i) = i
       end do
    end if
    !
    do m = 1, size(A)
       call xor4096( 0_i8, rn )
       l = m + int(rn * real( size(A) + 1 - m, dp ), i8)
       L1= A(l)
       A(l) = A(m)
       A(m) = L1
    end do
    !
  end subroutine ranper_i8
  !
end module mo_ranper

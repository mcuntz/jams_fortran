!> \file mo_combinatorics.f90

!> \brief Package for combinatorial calculations.

!> \details This package provides routines and functions for combinatorial calculations.

!> \authors Matthias Cuntz, Giovanni Dalmasso, Juliane Mai, Stephan Thober, Sebastian Mueller
!> \date Feb 2013

MODULE mo_combinatorics

  ! Written  Matthias Cuntz, Giovanni Dalmasso, Juliane Mai, Stephan Thober Feb 2013
  ! Modified Matthias Cuntz, May 2014 - removed numerical recipes, use mo_functions

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2011-2014 Matthias Cuntz, Giovanni Dalmasso, Juliane Mai, Stephan Thober, Sebastian Mueller
  !                         mc (at) macu (dot) de
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

  USE mo_kind,      ONLY: i4, i8, sp, dp
  USE mo_functions, only: factln, factorial

  IMPLICIT NONE

  PUBLIC :: binomcoeffi         ! Binomial coefficient (n choose k)
  PUBLIC :: random_kofn         ! Random selection of k of n
  PUBLIC :: next_kofn           ! Next selection of k of n to a given one
  PUBLIC :: all_kofn            ! All selections of k of n
  PUBLIC :: random_permut       ! Random permutation of a vector
  PUBLIC :: random_index_permut ! Random permutation of (1..n)
  PUBLIC :: next_index_permut   ! Next permutation of (1..n) to a given one
  PUBLIC :: all_index_permut    ! All permutations of (1..n)
  PUBLIC :: nextpart            ! Next partition of n to a given one
  PUBLIC :: prevpart            ! Previous partition of n to a given one
  PUBLIC :: isgoodpart          ! checks if a given array is a partition of n

  ! ------------------------------------------------------------------

  !     NAME
  !         binomcoeffi

  !     PURPOSE
  !
  !>        \brief The binomial coefficient.
  !
  !>        \details Calculates the binomial coefficient (n choose k):
  !>                     \f[ C(n,k) = \frac{n!}{k! (n-k)!} \f], 
  !>                 i.e. the number of possibilities to select from n numbers k
  !>                 for real input it also calculates the generalized binomial coefficient which is used for the binomial series
  !>                     \f[ C(x,k) = \prod_{i=1}^k \frac{x+1-i}{i} \f],
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8)/real(sp,dp) :: n/n(:)/x/x(:)"    from n numbers ...
  !>        \param[in] "integer(i4/i8)             :: k/k(:)"           ... k will be selected
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     integer(i4/i8)/real(sp,dp) :: binomcoeffi/binomcoeffi(:) &mdash; Binomial coefficient (n choose k)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         bico = binomcoeffi(5,2)
  !         bico --> 10
  !         bico = binomcoeffi(1.5,2)
  !         bico --> 0.375
  !         -> see also example in test directory

  !     LITERATURE
  !         

  !     HISTORY
  !>        \author Matthias Cuntz, Juliane Mai
  !>        \date Feb 2013
  !         Modified, Sebastian Mueller June 2014

  INTERFACE binomcoeffi
     MODULE PROCEDURE binomcoeffi_i4_d0, binomcoeffi_i8_d0, binomcoeffi_i4_d1, binomcoeffi_i8_d1,&
                      binomcoeffi_sp_d0, binomcoeffi_dp_d0, binomcoeffi_sp_d1, binomcoeffi_dp_d1,&
                      binomcoeffi_dpi4_d0, binomcoeffi_dpi4_d1
  END INTERFACE binomcoeffi
  
  ! ------------------------------------------------------------------

  !     NAME
  !         all_kofn

  !     PURPOSE
  !
  !>        \brief All possible selections of k numbers from n.
  !
  !>        \details Returns all possible selections of a k-subset from n numbers.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8) :: n"                 select from n numbers ...
  !>        \param[in] "integer(i4/i8) :: k"                 ... k numbers
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     integer(i4/i8), allocatable, dimension(:,:) :: alle &mdash; All selections. Dim_2=k.
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         alle = all_kofn(3, 2)
  !         alle --> (/ (/ 1,2 /), (/ 1,3 /), (/ 2,3 /) /) 
  !         -> see also example in test directory

  !     LITERATURE
  !         Nijenhuis, A., & Wilf, H. S. (1978). 
  !         Combinatorial algorithms for computers and calculators (2nd ed.). 
  !         Academic Pr.

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Feb 2013
  !         Modified, 

  INTERFACE all_kofn
     MODULE PROCEDURE all_kofn_i4, all_kofn_i8
  END INTERFACE all_kofn
  
  ! ------------------------------------------------------------------

  !     NAME
  !         next_kofn

  !     PURPOSE
  !
  !>        \brief The next selection of k numbers from n to a given one.
  !
  !>        \details Determines the next selection of k numbers from k to a given one.\n
  !>                 If one has n=5 numbers and want to pick k=3, the first possible selection is (1,2,3).
  !>                 The next selection will be (1,2,4). The subsequent selection to (1,3,4) will be (1,3,5).\n
  !>                 The given subset (previous) has to be sorted.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8) :: n"                 select from n numbers ...
  !>        \param[in] "integer(i4/i8) :: k"                 ... k numbers
  !>        \param[in] "integer(i4/i8) :: previous(:)"       previous selection (dim_1 = k)
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     integer(i4/i8) :: next(size(previous)) &mdash; Next selection of k from n numbers
  !
  !     RESTRICTIONS
  !>        \note The given subset (previous) has to be sorted.
  !
  !     EXAMPLE
  !         next = next_kofn(5, 3, (/1,3,4/))
  !         next --> (/ 1,3,5 /)
  !         -> see also example in test directory

  !     LITERATURE
  !         Nijenhuis, A., & Wilf, H. S. (1978). 
  !         Combinatorial algorithms for computers and calculators (2nd ed.). 
  !         Academic Pr.

  !     HISTORY
  !>        \author Giovanni Dalmasso, Juliane Mai
  !>        \date Feb 2013
  !         Modified, 

  INTERFACE next_kofn
     MODULE PROCEDURE next_kofn_i4, next_kofn_i8
  END INTERFACE next_kofn

  ! ------------------------------------------------------------------

  !     NAME
  !         random_kofn

  !     PURPOSE
  !
  !>        \brief A random selection of a k-subset of n numbers (1..n).
  !
  !>        \details Returns a random k-subset of n numbers (1..n).\n
  !>                 The returned subset will be sorted.\n
  !>                 The code is adapted from Nijenhuis (1978) - routine RANKSB.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8) :: n"                 select from n numbers ...
  !>        \param[in] "integer(i4/i8) :: k"                 ... k numbers
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !>        \param[in,out] "integer(i4/i8), optional :: save_state(n_save_state)"  an array for saving the state 
  !>                                                                               of an uniform random number stream
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     integer(i4/i8) :: set(k) &mdash; Random selection of k from n numbers
  !
  !     RESTRICTIONS
  !>        \note A random number stream for generating uniform random numbers has to be initialized before running this function. 
  !
  !     EXAMPLE
  !         set = random_kofn(5, 3)
  !         set --> (/ 2,3,5 /)
  !         -> see also example in test directory

  !     LITERATURE
  !         Nijenhuis, A., & Wilf, H. S. (1978). 
  !         Combinatorial algorithms for computers and calculators (2nd ed.). 
  !         Academic Pr.

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Feb 2013
  !         Modified, 

  INTERFACE random_kofn
     MODULE PROCEDURE random_kofn_i4 , random_kofn_i8
  END INTERFACE random_kofn

  ! ------------------------------------------------------------------

  !     NAME
  !         all_index_permut

  !     PURPOSE
  !
  !>        \brief All possible permutations of n integers 1..n.
  !
  !>        \details Determines all possible permutations of n integers 1..n. \n
  !>                 The number of possibilities is the factorial of n (n!).
  !>                 The code is adapted from Nijenhuis (1978) - routine NEXPER.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8) :: n"                 permutation of n integers 1..n
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     integer(i4/i8), allocatable :: permut(:,:) &mdash; All permutations of n integers \n
  !>                                                         size(permut,1) = n! \n
  !>                                                         size(permut,2) = n
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         permut = all_index_permut(3)
  !         permut --> (/ (/ 1,2,3 /), (/2,1,3/), (/3,1,2/), (/ 1,3,2 /), (/2,3,1/), (/3,2,1/) /)
  !         -> see also example in test directory

  !     LITERATURE
  !         Nijenhuis, A., & Wilf, H. S. (1978). 
  !         Combinatorial algorithms for computers and calculators (2nd ed.). 
  !         Academic Pr.

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Feb 2013
  !         Modified, 

  INTERFACE all_index_permut
     MODULE PROCEDURE all_index_permut_i4, all_index_permut_i8
  END INTERFACE all_index_permut

  ! ------------------------------------------------------------------

  !     NAME
  !         next_index_permut

  !     PURPOSE
  !
  !>        \brief The next permutation of n integers 1..n to a given one.
  !
  !>        \details Determines the next permutation of n integers 1..n to a given one.\n
  !>                 If one has n=5 numbers, the first possible permutation is (1,2,3,4,5).
  !>                 The next selection will be (2,1,3,4,5). \n
  !>                 The code is adapted from Nijenhuis (1978) - routine NEXPER.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8) :: n"                 permutation of n integers 1..n
  !>        \param[in] "integer(i4/i8) :: previous(:)"       previous permutation (dim_1 = n)
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     integer(i4/i8) :: next(n) &mdash; Next permutation of n integers
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         next = next_index_permut(5, (/1,2,3,4,5/))
  !         next --> (/ 2,1,3,4,5 /)
  !         -> see also example in test directory

  !     LITERATURE
  !         Nijenhuis, A., & Wilf, H. S. (1978). 
  !         Combinatorial algorithms for computers and calculators (2nd ed.). 
  !         Academic Pr.

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Feb 2013
  !         Modified, 

  INTERFACE next_index_permut
     MODULE PROCEDURE next_index_permut_i4, next_index_permut_i8
  END INTERFACE next_index_permut

  ! ------------------------------------------------------------------

  !     NAME
  !         random_index_permut

  !     PURPOSE
  !
  !>        \brief A random permutation of the integers 1..n.
  !
  !>        \details Returns a random permutation of n numbers (1..n)\n
  !>                 The code adapted from the Fortran library of A. Miller.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8) :: n"       permutation of integers 1..n
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !>        \param[in,out] "integer(i4/i8), optional :: save_state(n_save_state)"  an array for saving the state 
  !>                                                                               of an uniform random number stream
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     integer(i4/i8) :: set(n) &mdash; Random permutation of n integers 1..n
  !
  !     RESTRICTIONS
  !>        \note A random number stream for generating uniform random numbers has to be 
  !>              initialized before running this function. 
  !
  !     EXAMPLE
  !         set = random_index_permut(5, 3)
  !         set --> (/ 2,1,3,5,4 /)
  !         -> see also example in test directory

  !     LITERATURE
  !         A. Miller, CSIRO Mathematical & Information Sciences, 
  !         Clayton 3169, Victoria, Australia, Version 1.13, 2 October 2000.
  !         http://jblevins.org/mirror/amiller/

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date Feb 2013
  !         Modified, 

  INTERFACE random_index_permut
     MODULE PROCEDURE random_index_permut_i4, random_index_permut_i8
  END INTERFACE random_index_permut

  ! ------------------------------------------------------------------

  !     NAME
  !         random_permut

  !     PURPOSE
  !         Calculates a random permutation of n numbers
  !
  !>        \brief Random permutation.
  !
  !>        \details Randomly permutes the elements of a given vector.
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !>        \param[inout] "integer(i4/i8)/real(sp/dp) :: vec(:)"   1D-array to be shuffled

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !>        \param[in,out] "integer(i4/i8), optional :: save_state(n_save_state)"  an array for saving the state 
  !>                                                                               of a uniform random number stream
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !>        \note A random number stream for generating uniform random numbers has to be 
  !>              initialized before running this function. 
  !
  !     EXAMPLE
  !         a = (/5,3,2/)
  !         call random_permut(a)
  !         a --> (/ 2,5,3 /)
  !         -> see also example in test directory

  !     LITERATURE
  !         Nijenhuis, A., & Wilf, H. S. (1978).\n
  !         Combinatorial algorithms for computers and calculators (2nd ed.). Academic Pr., p. 63
  !         Algorithm RANPER

  !     HISTORY
  !>        \author Stephan Thober, Juliane Mai, Matthias Cuntz
  !>        \date Feb 2013

  INTERFACE random_permut
     MODULE PROCEDURE random_permut_i4, random_permut_i8, random_permut_sp, random_permut_dp
  END INTERFACE random_permut

  ! ------------------------------------------------------------------

  !     NAME
  !         nextpart

  !     PURPOSE
  !         Calculates the next partition of n to a given one
  !
  !>        \brief next partition of n to a given one
  !
  !>        \details next greater partition of n to a given one where the partitions are sorted by inverse lexicographical ordering
  !>                 for example: (2,0) > (0,1)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8), dimension(:), intent(in) :: part"   a given partition not equal to (n,0,..,0) [not last]
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return "integer(i4/i8), dimension(:)                :: nextpart"   next partition of n
  !
  !     RESTRICTIONS
  !>        \note the given partiton should not be the greatest one i.e. not (n,0,..,0) 
  !
  !     EXAMPLE
  !         nextpart((0,2,0,0)) = (2,1,0,0)

  !     LITERATURE
  !         none


  !     HISTORY
  !>        \author Sebastian Mueller
  !>        \date June 2014

  INTERFACE nextpart
     MODULE PROCEDURE nextpart_i4, nextpart_i8
  END INTERFACE nextpart

  ! ------------------------------------------------------------------

  !     NAME
  !         prevpart

  !     PURPOSE
  !         Calculates the previous partition of n to a given one
  !
  !>        \brief previous partition of n to a given one
  !
  !>        \details next smaller partition of n to a given one where the partitions are sorted by inverse lexicographical ordering
  !>                 for example: (2,0) > (0,1)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8), dimension(:), intent(in) :: part"   a given partition not equal to (0,..,0,1) [not first]
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return "integer(i4/i8), dimension(:)                :: prevpart"   previous partition of n
  !
  !     RESTRICTIONS
  !>        \note the given partiton should not be the smallest one i.e. not (0,..,0,1) 
  !
  !     EXAMPLE
  !         prevpart((2,1,0,0)) = (0,2,0,0)

  !     LITERATURE
  !         none


  !     HISTORY
  !>        \author Sebastian Mueller
  !>        \date July 2014

  INTERFACE prevpart
     MODULE PROCEDURE prevpart_i4, prevpart_i8
  END INTERFACE prevpart

  ! ------------------------------------------------------------------

  !     NAME
  !         isgoodpart

  !     PURPOSE
  !         Checks if a given array is a partition of n
  !
  !>        \brief Checks if a given array is a partition of n
  !
  !>        \details Checks if a given array is a partition of n, where p=(p_1,..,p_n) is a partition of n iff
  !>                 \f[ \sum_{i=1}^n p_i*i = p_1 + 2*p_2 + 3*p_3 + .. =n \f]
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8), dimension(:), intent(in) :: part"   a given partition not equal to (n,0,..,0)
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return "logical                                     :: isgoodpart"
  !
  !     RESTRICTIONS
  !>        \note the given array should be of dimension greater than one
  !
  !     EXAMPLE
  !         isgoodpart((2,1,0)) = true
  !         isgoodpart((0,2))   = false

  !     LITERATURE
  !         none


  !     HISTORY
  !>        \author Sebastian Mueller
  !>        \date June 2014

  INTERFACE isgoodpart
     MODULE PROCEDURE isgoodpart_i4, isgoodpart_i8
  END INTERFACE isgoodpart

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  FUNCTION binomcoeffi_i4_d0(n,k)
    ! Returns the binomial coefficient as an integer number.
    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: n, k
    INTEGER(i4)             :: binomcoeffi_i4_d0

    binomcoeffi_i4_d0 = nint(exp(factln(n)-factln(k)-factln(n-k)))

  END FUNCTION binomcoeffi_i4_d0

  FUNCTION binomcoeffi_i8_d0(n,k)
    ! Returns the binomial coefficient as an integer number.
    IMPLICIT NONE

    INTEGER(i8), INTENT(IN) :: n, k
    INTEGER(i8)             :: binomcoeffi_i8_d0

    binomcoeffi_i8_d0 = nint(exp(factln(n)-factln(k)-factln(n-k)))

  END FUNCTION binomcoeffi_i8_d0

  FUNCTION binomcoeffi_i4_d1(n,k)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:), INTENT(IN) :: n, k
    INTEGER(i4), DIMENSION(size(n))       :: binomcoeffi_i4_d1

    if (size(n) /= size(k)) stop 'Error binomcoeffi: size(n) /= size(k)'

    binomcoeffi_i4_d1 = nint(exp(factln(n)-factln(k)-factln(n-k)))

  END FUNCTION binomcoeffi_i4_d1

  FUNCTION binomcoeffi_i8_d1(n,k)

    IMPLICIT NONE

    INTEGER(i8), DIMENSION(:), INTENT(IN) :: n, k
    INTEGER(i8), DIMENSION(size(n))       :: binomcoeffi_i8_d1

    if (size(n) /= size(k)) stop 'Error binomcoeffi: size(n) /= size(k)'

    binomcoeffi_i8_d1 = nint(exp(factln(n)-factln(k)-factln(n-k)))

  END FUNCTION binomcoeffi_i8_d1

!generalized binomial coefficient
  FUNCTION binomcoeffi_sp_d0(x,k)
    ! Returns the generalized binomial coefficient as a real number.
    IMPLICIT NONE

    real(sp),    INTENT(IN) :: x
    INTEGER(i4), INTENT(IN) :: k
    real(sp)                :: binomcoeffi_sp_d0

    integer(i4)             :: i
    
    if (k>0_i4) then 
        binomcoeffi_sp_d0       =   1.0_sp    
        do i=1_i4,k
            binomcoeffi_sp_d0   =   binomcoeffi_sp_d0*(x+1.0_sp-real(i,sp))/real(i,sp)
        end do
    else if (k<0_i4) then
        binomcoeffi_sp_d0       =   0.0_sp
    else
        binomcoeffi_sp_d0       =   1.0_sp
    endif
    
  END FUNCTION binomcoeffi_sp_d0

  FUNCTION binomcoeffi_dp_d0(x,k)
    ! Returns the generalized binomial coefficient as a real number.
    IMPLICIT NONE

    real(dp),    INTENT(IN) :: x
    INTEGER(i8), INTENT(IN) :: k
    real(dp)                :: binomcoeffi_dp_d0

    integer(i8)             :: i
    
    if (k>0_i8) then 
        binomcoeffi_dp_d0       =   1.0_sp    
        do i=1_i8,k
            binomcoeffi_dp_d0   =   binomcoeffi_dp_d0*(x+1.0_dp-real(i,dp))/real(i,dp)
        end do
    else if (k<0_i8) then
        binomcoeffi_dp_d0       =   0.0_dp
    else
        binomcoeffi_dp_d0       =   1.0_dp
    endif

  END FUNCTION binomcoeffi_dp_d0

  FUNCTION binomcoeffi_sp_d1(x,k)

    IMPLICIT NONE

    real(sp),    DIMENSION(:), INTENT(IN) :: x
    INTEGER(i4), DIMENSION(:), INTENT(IN) :: k
    real(sp),    DIMENSION(size(x))       :: binomcoeffi_sp_d1

    integer(i4)             :: i

    if (size(x) /= size(k)) stop 'Error binomcoeffi: size(n) /= size(k)'

    do i=1_i4, size(x)
        binomcoeffi_sp_d1(i) = binomcoeffi_sp_d0(x(i),k(i))
    end do

  END FUNCTION binomcoeffi_sp_d1

  FUNCTION binomcoeffi_dp_d1(x,k)

    IMPLICIT NONE

    real(dp),    DIMENSION(:), INTENT(IN) :: x
    INTEGER(i8), DIMENSION(:), INTENT(IN) :: k
    real(dp),    DIMENSION(size(x))       :: binomcoeffi_dp_d1

    integer(i8)             :: i

    if (size(x) /= size(k)) stop 'Error binomcoeffi: size(n) /= size(k)'

    do i=1_i8, size(x)
        binomcoeffi_dp_d1(i) = binomcoeffi_dp_d0(x(i),k(i))
    end do

  END FUNCTION binomcoeffi_dp_d1

  FUNCTION binomcoeffi_dpi4_d0(x,k)
    ! Returns the generalized binomial coefficient as a real number.
    IMPLICIT NONE

    real(dp),    INTENT(IN) :: x
    INTEGER(i4), INTENT(IN) :: k
    real(dp)                :: binomcoeffi_dpi4_d0

    integer(i4)             :: i
    
    if (k>0_i4) then 
        binomcoeffi_dpi4_d0     =   1.0_sp    
        do i=1_i4,k
            binomcoeffi_dpi4_d0 =   binomcoeffi_dpi4_d0*(x+1.0_dp-real(i,dp))/real(i,dp)
        end do
    else if (k<0_i4) then
        binomcoeffi_dpi4_d0     =   0.0_dp
    else
        binomcoeffi_dpi4_d0     =   1.0_dp
    endif

  END FUNCTION binomcoeffi_dpi4_d0

  FUNCTION binomcoeffi_dpi4_d1(x,k)

    IMPLICIT NONE

    real(dp),    DIMENSION(:), INTENT(IN) :: x
    INTEGER(i4), DIMENSION(:), INTENT(IN) :: k
    real(dp),    DIMENSION(size(x))       :: binomcoeffi_dpi4_d1

    integer(i4)             :: i

    if (size(x) /= size(k)) stop 'Error binomcoeffi: size(n) /= size(k)'

    do i=1_i4, size(x)
        binomcoeffi_dpi4_d1(i) = binomcoeffi_dpi4_d0(x(i),k(i))
    end do

  END FUNCTION binomcoeffi_dpi4_d1

  ! ------------------------------------------------------------------

  FUNCTION all_kofn_i4(n, k) result(alle)

    IMPLICIT NONE

    INTEGER(I4),                         INTENT(IN)  :: n             ! from n numbers will be selected
    INTEGER(I4),                         INTENT(IN)  :: k             ! k numbers will be selected
    INTEGER(I4), DIMENSION(:,:), allocatable         :: alle           ! all subsets

    ! local variables
    integer(i4)         :: bico
    integer(i4)         :: i

    bico = binomcoeffi(n,k)
    allocate(alle(bico,k))

    forall(i=1:k) alle(1,i) = i

    do i=2,bico
       alle(i,:) = next_kofn(n, k,  alle(i-1,:))
    end do

  END FUNCTION all_kofn_i4

  FUNCTION all_kofn_i8(n, k) result(alle)

    IMPLICIT NONE

    INTEGER(I8),                         INTENT(IN)  :: k             ! k numbers will be selected
    INTEGER(I8),                         INTENT(IN)  :: n             ! from n numbers will be selected
    INTEGER(I8), DIMENSION(:,:), allocatable         :: alle          ! all subsets

    ! local variables
    integer(i8)         :: bico
    integer(i8)         :: i

    bico = binomcoeffi(n,k)
    allocate(alle(bico,k))

    forall(i=1:k) alle(1,i) = i

    do i=2,bico
       alle(i,:) = next_kofn(n, k,  alle(i-1,:))
    end do

  END FUNCTION all_kofn_i8

  ! ------------------------------------------------------------------

  FUNCTION next_kofn_i4(n, k, previous) result(next)

    IMPLICIT NONE

    INTEGER(I4),                         INTENT(IN)  :: k             ! k numbers will be selected
    INTEGER(I4),                         INTENT(IN)  :: n             ! from n numbers will be selected
    INTEGER(I4), DIMENSION(k),           INTENT(IN)  :: previous      ! previous subset
    INTEGER(I4), DIMENSION(k)                        :: next          ! next subset

    ! local variables
    integer(i4)         :: indx
    logical             :: found
    integer(i4)         :: i

    ! find index from the back to start changing the values
    indx = k
    found = .false.
    do while ( (.not. found) .and. (indx .gt. 1) ) 
       if ( previous(indx) .eq. n-k+indx ) then
          indx = indx-1
       else
          found = .true.
       end if
    end do
    if( (indx .gt. 1) .or. ( previous(1) .ne. (n-k+1)) ) then
       ! it is a subset in between
       next = previous
       next(indx) = previous(indx) + 1
       forall(i=indx+1:k) next(i) = next(indx) + i-indx
    else
       ! there do not exist a next subset and first one (/1, ..., k/) is returned
       forall(i=1:k) next(i) = i
    end if        

  END FUNCTION next_kofn_i4

  FUNCTION next_kofn_i8(n, k, previous) result(next)

    IMPLICIT NONE

    INTEGER(I8),                         INTENT(IN)  :: k             ! k numbers will be selected
    INTEGER(I8),                         INTENT(IN)  :: n             ! from n numbers will be selected
    INTEGER(I8), DIMENSION(k),           INTENT(IN)  :: previous      ! previous subset
    INTEGER(I8), DIMENSION(k)                        :: next          ! next subset

    ! local variables
    integer(i8)         :: indx
    logical             :: found
    integer(i8)         :: i

    ! find index from the back to start changing the values
    indx = k
    found = .false.
    do while ( (.not. found) .and. (indx .gt. 1) ) 
       if ( previous(indx) .eq. n-k+indx ) then
          indx = indx-1
       else
          found = .true.
       end if
    end do

     if( (indx .gt. 1) .or. ( previous(1) .ne. (n-k+1)) ) then
       ! it is a subset in between
       next       = previous
       next(indx) = previous(indx) + 1
       forall(i=indx+1:k) next(i) = next(indx) + i-indx
    else
       ! there do not exist a next subset and first one (/1, ..., k/) is returned
       forall(i=1:k) next(i) = i
    end if        

  END FUNCTION next_kofn_i8

  ! ------------------------------------------------------------------

  FUNCTION random_kofn_i4(n, k, save_state) result(set)
    
    use mo_xor4096, only: xor4096, n_save_state

    implicit none

    INTEGER(I4),                                    INTENT(IN)    :: k             ! k numbers will be selected
    INTEGER(I4),                                    INTENT(IN)    :: n             ! from n numbers will be selected
    INTEGER(I4), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state    ! for saving state of a uniform, sp 
    !                                                                              ! random number stream
    INTEGER(I4), DIMENSION(k)                                     :: set           ! random subset

    ! local variables
    logical     :: gout
    integer(i4) :: x, r, ds, p, s, c, i, l, m,m0, iseed
    real(sp)    :: rn

    iseed = 0_i4

    ! (A)
    c=k
    forall(i=1:k) set(i) = ((i-1)*n)/k
    
    ! (B)
    do while (c .gt. 0)
       if (present(save_state)) then
          call xor4096(iSeed, rn, save_state=save_state)
       else
          call xor4096(iSeed, rn)
       end if
       x = 1 + int( real(n,sp)*rn, i4 )
       l = 1 + int( real(x*k-1,sp)/real(n,sp), i4 )
       if(x.le.set(l)) cycle
       set(l) = set(l) + 1
       c = c - 1
    end do
    i = 0
    p = 0
    s = k
    m = 0

    ! (C)
    gout = .false.
    do while (.not. gout)
       i = i+1
       if (i .gt. k) then
          gout = .true.
       else
          if (set(i) .eq. int(real((i-1)*n,sp)/real(k),i4)) then
             set(i) = 0
          else
             p=p+1
             m=set(i)
             set(i) = 0
             set(p) = m
          end if
       end if
    end do

    ! (D)
    gout = .false.
    do while (.not. gout)
       l  = 1+int(real(set(p)*k-1,sp)/real(n,sp),i4)
       ds = set(p) - int( real((l-1)*n,sp)/real(k,sp) ,i4)
       set(p) = 0
       set(s) = l
       s=s-ds
       p=p-1
       if (p .le. 0) gout = .true.
    end do
    l = k
    m0 = 1+int(real((set(k)-1)*n,sp)/real(k,sp),i4)
    r = 0

    do while( l .ne. 0)
       
       ! (E)
       if (set(l) .ne. 0) then
          r = l
          m0 = 1+int(real((set(l)-1)*n,sp)/real(k,sp),i4)
          m = int(real(set(l)* n, sp) / real(k,sp),i4) - m0 + 1
       end if
       
       ! (F)
       if (present(save_state)) then
          call xor4096(iSeed, rn, save_state=save_state)
       else
          call xor4096(iSeed, rn)
       end if
       x = m0 + int(real(m,sp) * rn,i4)
       i = l
       
       ! (G)
       gout = .false.
       do while (.not. gout)
          i = i + 1
          if (i .le. r) then
             if( x .ge. set(i)) then
                set(i-1) = set(i)
                x = x+1
             else 
                gout = .true.
             end if
          else
             gout = .true.
          end if
       end do
       
       ! (H)
       set(i-1) = x
       m = m-1
       l=l-1

    end do

  END FUNCTION random_kofn_i4

  FUNCTION random_kofn_i8(n, k, save_state) result(set)
    
    use mo_xor4096, only: xor4096, n_save_state

    implicit none

    INTEGER(I8),                                    INTENT(IN)    :: k             ! k numbers will be selected
    INTEGER(I8),                                    INTENT(IN)    :: n             ! from n numbers will be selected
    INTEGER(I8), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state    ! for saving state of a uniform, dp 
    !                                                                              ! random number stream
    INTEGER(I8), DIMENSION(k)                                     :: set           ! random subset

    ! local variables
    logical     :: gout
    integer(i8) :: x, r, ds, p, s, c, i, l, m,m0, iseed
    real(dp)    :: rn

    iseed = 0_i8

    ! (A)
    c=k
    forall(i=1:k) set(i) = ((i-1)*n)/k
    
    ! (B)
    do while (c .gt. 0)
       if (present(save_state)) then
          call xor4096(iSeed, rn, save_state=save_state)
       else
          call xor4096(iSeed, rn)
       end if
       x = 1 + int( real(n,dp)*rn, i8 )
       l = 1 + int( real(x*k-1,dp)/real(n,dp), i8 )
       if(x.le.set(l)) cycle
       set(l) = set(l) + 1
       c = c - 1
    end do
    i = 0_i8
    p = 0_i8
    s = k
    m = 0_i8

    ! (C)
    gout = .false.
    do while (.not. gout)
       i = i+1
       if (i .gt. k) then
          gout = .true.
       else
          if (set(i) .eq. int(real((i-1)*n,dp)/real(k),i8)) then
             set(i) = 0
          else
             p=p+1
             m=set(i)
             set(i) = 0
             set(p) = m
          end if
       end if
    end do

    ! (D)
    gout = .false.
    do while (.not. gout)
       l  = 1+int(real(set(p)*k-1,dp)/real(n,dp),i8)
       ds = set(p) - int( real((l-1)*n,dp)/real(k,dp) ,i8)
       set(p) = 0
       set(s) = l
       s=s-ds
       p=p-1
       if (p .le. 0) gout = .true.
    end do
    l = k
    m0 = 1_i8+int(real((set(k)-1)*n,dp)/real(k,dp),i8)
    r = 0_i8

    do while( l .ne. 0)
       
       ! (E)
       if (set(l) .ne. 0) then
          r = l
          m0 = 1+int(real((set(l)-1)*n,dp)/real(k,dp),i8)
          m = int(real(set(l)* n, dp) / real(k,dp),i8) - m0 + 1
       end if
       
       ! (F)
       if (present(save_state)) then
          call xor4096(iSeed, rn, save_state=save_state)
       else
          call xor4096(iSeed, rn)
       end if
       x = m0 + int(real(m,dp) * rn,i8)
       i = l
       
       ! (G)
       gout = .false.
       do while (.not. gout)
          i = i + 1
          if (i .le. r) then
             if( x .ge. set(i)) then
                set(i-1) = set(i)
                x = x+1
             else 
                gout = .true.
             end if
          else
             gout = .true.
          end if
       end do
       
       ! (H)
       set(i-1) = x
       m = m-1
       l=l-1

    end do

  END FUNCTION random_kofn_i8

  ! ------------------------------------------------------------------

  FUNCTION all_index_permut_i4(n) result(alle)

    IMPLICIT NONE

    INTEGER(I4),                         INTENT(IN)  :: n             ! n numbers will be permuted
    INTEGER(I4), DIMENSION(:,:), allocatable         :: alle           ! all permutations

    ! local variables
    integer(i4)         :: fact
    integer(i4)         :: i

    fact = factorial(n)
    allocate(alle(fact,n))

    forall(i=1:n) alle(1,i) = i

    do i=2,fact
       alle(i,:) = next_index_permut(n, alle(i-1,:))
    end do

  END FUNCTION all_index_permut_i4

  FUNCTION all_index_permut_i8(n) result(alle)

    IMPLICIT NONE

    INTEGER(I8),                         INTENT(IN)  :: n             ! n numbers will be permuted
    INTEGER(I8), DIMENSION(:,:), allocatable         :: alle           ! all permutations

    ! local variables
    integer(i8)         :: fact
    integer(i8)         :: i

    fact = factorial(n)
    allocate(alle(fact,n))

    forall(i=1:n) alle(1,i) = i

    do i=2,fact
       alle(i,:) = next_index_permut(n, alle(i-1,:))
    end do

  END FUNCTION all_index_permut_i8

  ! ------------------------------------------------------------------

  FUNCTION next_index_permut_i4(n, previous) result(next)

    implicit none

    integer(i4),               intent(in) :: n
    integer(i4), dimension(n), intent(in) :: previous
    integer(i4), dimension(n)             :: next

    ! local variables
    integer(i4)                 :: i, j
    integer(i4), dimension(n-1) :: s, d
    integer(i4)                 :: indx, indx2, found_case
    integer(i4)                 :: largest, smallest
    logical                     :: last

    ! calculate d
    do i=1,n-1
       d(i) = 0
       s(i) = 0
       do j=1,i
          if (previous(j) .gt. previous(i+1)) then
             d(i) = d(i)+1
          end if
       end do
       do j=1, i
          s(i) = s(i) + d(j)
       end do
    end do

    ! check if last one
    last = .true.
    if( mod(n,2_i4) .eq. 0_i4 ) then
       ! (a) n odd:  last sequence is (2, 3, 4 ... n, 1)
       do i=1,n-1
          last = last .and. (previous(i) .eq. (i+1))
       end do
       last = last .and. (previous(n) .eq. 1)
    else
       ! (b) n even: last sequence is (3, 4, 5, ..., n, 2, 1)
       do i=1,n-2
          last = last .and. (previous(i) .eq. (i+2))
       end do
       last = last .and. (previous(n-1) .eq. 2)
       last = last .and. (previous(n)   .eq. 1)
    end if
    if (last) then
       do i=1,n
          next(i) = i
       end do
       ! stop already here
       return
    end if
    
    if ( mod(s(n-1),2_i4) .eq. 0_i4) then
       ! if s(n-1) is even interchange previous(1) and previous(2)
       next = previous
       next(2) = previous(1)
       next(1) = previous(2)
    else
       ! if s(n-1) is odd do some more work...
       ! (1) find index and case: 
       !         case 1 :  d(indx) < indx && s(indx) odd
       !         case 2 :  d(indx) > 0    && s(indx) even
       !         indx where case 1 or 2 fulfilled first time
       found_case = 0_i4
       indx = 0_i4
       do while (found_case .eq. 0_i4)
          indx = indx + 1_i4
          ! case 1 ?
          if ( (mod(s(indx),2_i4) .eq. 1_i4) .and. (d(indx) .lt. indx) ) then
             found_case = 1_i4
          else 
             ! case 2 ?
             if( (mod(s(indx),2_i4) .eq. 0_i4) .and. (d(indx) .gt. 0_i4) ) then
                found_case = 2_i4
             end if
          end if
       end do
       ! (2a) case 1 fulfilled: search largest number previous(i),(i=1,indx) 
       !                        less than previous(indx+1)
       ! (2b) case 2 fulfilled: search smallest number previous(i),(i=1,indx) 
       !                        greater than previous(indx+1)
       select case(found_case)
          case(1)
             indx2 = indx+1
             largest = 1
             do i=1,indx
                if( (previous(i) .lt. previous(indx+1)) .and. (largest .le. previous(i)) ) then
                   indx2 = i
                   largest = previous(i)
                end if
             end do
          case(2)
             indx2 = indx+1
             smallest = n
             do i=1,indx
                if( (previous(i) .gt. previous(indx+1)) .and. (smallest .ge. previous(i)) ) then
                   indx2 = i
                   smallest = previous(i)
                end if
             end do
          case default
             stop 'next_index_permut_i4: something went wrong here'
       end select
       ! (3) swap previous(indx2) and previous(indx+1)
       next = previous
       next(indx2) = previous(indx+1)
       next(indx+1) = previous(indx2)
    end if

  END FUNCTION next_index_permut_i4

  FUNCTION next_index_permut_i8(n, previous) result(next)

    implicit none

    integer(i8),               intent(in) :: n
    integer(i8), dimension(n), intent(in) :: previous
    integer(i8), dimension(n)             :: next

    ! local variables
    integer(i8)                 :: i, j
    integer(i8), dimension(n-1) :: s, d
    integer(i8)                 :: indx, indx2, found_case
    integer(i8)                 :: largest, smallest
    logical                     :: last

    ! calculate d
    do i=1,n-1
       d(i) = 0
       s(i) = 0
       do j=1,i
          if (previous(j) .gt. previous(i+1)) then
             d(i) = d(i)+1
          end if
       end do
       do j=1, i
          s(i) = s(i) + d(j)
       end do
    end do

    ! check if last one
    last = .true.
    if( mod(n,2_i8) .eq. 0_i8 ) then
       ! (a) n odd:  last sequence is (2, 3, 4 ... n, 1)
       do i=1,n-1
          last = last .and. (previous(i) .eq. (i+1))
       end do
       last = last .and. (previous(n) .eq. 1)
    else
       ! (b) n even: last sequence is (3, 4, 5, ..., n, 2, 1)
       do i=1,n-2
          last = last .and. (previous(i) .eq. (i+2))
       end do
       last = last .and. (previous(n-1) .eq. 2)
       last = last .and. (previous(n)   .eq. 1)
    end if
    if (last) then
       do i=1,n
          next(i) = i
       end do
       ! stop already here
       return
    end if
    
    if ( mod(s(n-1),2_i8) .eq. 0_i8) then
       ! if s(n-1) is even interchange previous(1) and previous(2)
       next = previous
       next(2) = previous(1)
       next(1) = previous(2)
    else
       ! if s(n-1) is odd do some more work...
       ! (1) find index and case: 
       !         case 1 :  d(indx) < indx && s(indx) odd
       !         case 2 :  d(indx) > 0    && s(indx) even
       !         indx where case 1 or 2 fulfilled first time
       found_case = 0_i8
       indx = 0_i8
       do while (found_case .eq. 0_i8)
          indx = indx + 1_i8
          ! case 1 ?
          if ( (mod(s(indx),2_i8) .eq. 1_i8) .and. (d(indx) .lt. indx) ) then
             found_case = 1_i8
          else 
             ! case 2 ?
             if( (mod(s(indx),2_i8) .eq. 0_i8) .and. (d(indx) .gt. 0_i8) ) then
                found_case = 2_i8
             end if
          end if
       end do
       ! (2a) case 1 fulfilled: search largest number previous(i),(i=1,indx) 
       !                        less than previous(indx+1)
       ! (2b) case 2 fulfilled: search smallest number previous(i),(i=1,indx) 
       !                        greater than previous(indx+1)
       select case(found_case)
          case(1)
             indx2 = indx+1
             largest = 1
             do i=1,indx
                if( (previous(i) .lt. previous(indx+1)) .and. (largest .le. previous(i)) ) then
                   indx2 = i
                   largest = previous(i)
                end if
             end do
          case(2)
             indx2 = indx+1
             smallest = n
             do i=1,indx
                if( (previous(i) .gt. previous(indx+1)) .and. (smallest .ge. previous(i)) ) then
                   indx2 = i
                   smallest = previous(i)
                end if
             end do
          case default
             stop 'next_index_permut_i8: something went wrong here'
       end select
       ! (3) swap previous(indx2) and previous(indx+1)
       next = previous
       next(indx2) = previous(indx+1)
       next(indx+1) = previous(indx2)
    end if

  END FUNCTION next_index_permut_i8

  ! ------------------------------------------------------------------

  SUBROUTINE random_permut_i4(vec, save_state)
    !
    ! Generate a random permuation of the inout vector vec
    !
    use mo_xor4096, only: xor4096, n_save_state

    implicit none

    INTEGER(i4), DIMENSION(:),                      INTENT(INOUT) :: vec
    INTEGER(i4), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    !     Local variables
    INTEGER(i4) :: l, m, nvec
    INTEGER(i4) :: iSeed
    INTEGER(i4) :: L1
    REAL(sp)    :: rn

    iSeed = 0_i4

    nvec = size(vec)
    do m = 1, nvec
       if (present(save_state)) then
          call xor4096(iSeed, rn, save_state=save_state)
       else
          call xor4096(iSeed, rn)
       end if
       l      = m + int(rn * real(nvec + 1 - m, sp), i4)
       L1     = vec(l)
       vec(l) = vec(m)
       vec(m) = L1
    end do

  END SUBROUTINE random_permut_i4

  SUBROUTINE random_permut_i8(vec, save_state)
    !
    ! Generate a random permuation of the inout vector vec
    !
    use mo_xor4096, only: xor4096, n_save_state

    implicit none

    INTEGER(i8), DIMENSION(:),                      INTENT(INOUT) :: vec
    INTEGER(i8), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    !     Local variables
    INTEGER(i4) :: l, m, nvec
    INTEGER(i8) :: iSeed
    INTEGER(i8) :: L1
    REAL(dp)    :: rn

    iSeed = 0_i8

    nvec = size(vec)
    do m = 1, nvec
       if (present(save_state)) then
          call xor4096(iSeed, rn, save_state=save_state)
       else
          call xor4096(iSeed, rn)
       end if
       l      = m + int(rn * real(nvec + 1 - m, dp), i4)
       L1     = vec(l)
       vec(l) = vec(m)
       vec(m) = L1
    end do

  END SUBROUTINE random_permut_i8

  SUBROUTINE random_permut_sp(vec, save_state)
    !
    ! Generate a random permuation of the inout vector vec
    !
    use mo_xor4096, only: xor4096, n_save_state

    implicit none

    REAL(sp),    DIMENSION(:),                      INTENT(INOUT) :: vec
    INTEGER(i4), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    !     Local variables
    INTEGER(i4) :: l, m, nvec
    INTEGER(i4) :: iSeed
    REAL(sp)    :: L1
    REAL(sp)    :: rn

    iSeed = 0_i4

    nvec = size(vec)
    do m = 1, nvec
       if (present(save_state)) then
          call xor4096(iSeed, rn, save_state=save_state)
       else
          call xor4096(iSeed, rn)
       end if
       l      = m + int(rn * real(nvec + 1 - m, sp), i4)
       L1     = vec(l)
       vec(l) = vec(m)
       vec(m) = L1
    end do

  END SUBROUTINE random_permut_sp

  SUBROUTINE random_permut_dp(vec, save_state)
    !
    ! Generate a random permuation of the inout vector vec
    !
    use mo_xor4096, only: xor4096, n_save_state

    implicit none

    REAL(dp),    DIMENSION(:),                      INTENT(INOUT) :: vec
    INTEGER(i8), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    !     Local variables
    INTEGER(i4) :: l, m, nvec
    INTEGER(i8) :: iSeed
    REAL(dp)    :: L1
    REAL(dp)    :: rn

    iSeed = 0_i8

    nvec = size(vec)
    do m = 1, nvec
       if (present(save_state)) then
          call xor4096(iSeed, rn, save_state=save_state)
       else
          call xor4096(iSeed, rn)
       end if
       l      = m + int(rn * real(nvec + 1 - m, dp), i4)
       L1     = vec(l)
       vec(l) = vec(m)
       vec(m) = L1
    end do

  END SUBROUTINE random_permut_dp

  ! ------------------------------------------------------------------

  FUNCTION random_index_permut_i4(n, save_state)
    !
    ! Generate a random ordering of the integers 1 ... n.
    !
    use mo_xor4096,      only: n_save_state

    implicit none

    INTEGER(i4),                                    INTENT(IN)    :: n
    INTEGER(i4), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state
    INTEGER(i4), DIMENSION(n)                                     :: random_index_permut_i4

    INTEGER(i4) :: i

    forall(i=1:n) random_index_permut_i4(i) = i
    if (present(save_state)) then
       call random_permut(random_index_permut_i4, save_state=save_state)
    else
       call random_permut(random_index_permut_i4)
    endif

  END FUNCTION random_index_permut_i4

  FUNCTION random_index_permut_i8(n, save_state)
    !
    ! Generate a random ordering of the integers 1 ... n.
    !
    use mo_xor4096,      only: n_save_state

    implicit none

    INTEGER(i8),                                    INTENT(IN)    :: n
    INTEGER(i8), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state
    INTEGER(i8), DIMENSION(n)                                     :: random_index_permut_i8

    INTEGER(i8) :: i

    forall(i=1_i8:n) random_index_permut_i8(i) = i
    if (present(save_state)) then
       call random_permut(random_index_permut_i8, save_state=save_state)
    else
       call random_permut(random_index_permut_i8)
    endif

  END FUNCTION random_index_permut_i8

  ! ------------------------------------------------------------------

  function nextpart_i4(part)

    implicit none

    integer(i4), dimension(:), intent(in) :: part   
    integer(i4), dimension(size(part))    :: nextpart_i4
    
    integer(i4)                           :: n, firstindex, summ, temp


    !Check if input is valid
    if (.not. isgoodpart(part))     stop 'Not a valid input partition.'
    if (part(1) .eq. size(part))    stop 'Input was the last possible partition.'

    nextpart_i4             =   part    

    !set the first entry to zero  and remember it
    temp                    =   nextpart_i4(1)
    nextpart_i4(1)          =   0_i4    

    !search for the first index that is not equal to zero
    firstindex  =   1_i4
    do n = 2_i4, size(nextpart_i4)
        if (nextpart_i4(n) .ne. 0_i4) then
            firstindex      =   n
            exit
        endif
    end do

    !decrease the first entry that is not equal to zero by one
    nextpart_i4(firstindex) =   nextpart_i4(firstindex)-1_i4

    !calculate the actual sum of the partition
    summ                    =   size(nextpart_i4)-firstindex-temp

    !fill up the partition below the index from above, until the sum equals n again       
    do n=firstindex-1_i4,1_i4,-1_i4
        temp                =   nextpart_i4(n)
        !set the next entry as big as possible
        nextpart_i4(n)      =   FLOOR(real(size(nextpart_i4)-summ,sp)/real(n,sp))
        !calculate the actual sum of the partition
        summ                =   summ + n*(nextpart_i4(n) - temp) 
        !if the sum is n we are done
        if (summ == size(nextpart_i4)) exit
    end do

  end function nextpart_i4

  function nextpart_i8(part)

    implicit none

    integer(i8), dimension(:), intent(in) :: part   
    integer(i8), dimension(size(part))    :: nextpart_i8
    
    integer(i8)                           :: n, firstindex, summ, temp


    !Check if input is valid
    if (.not. isgoodpart(part))     stop 'Not a valid input partition.'
    if (part(1) .eq. size(part))    stop 'Input was the last possible partition.'

    nextpart_i8             =   part    

    !set the first entry to zero  and remember it
    temp                    =   nextpart_i8(1)
    nextpart_i8(1)          =   0_i8    

    !search for the first index that is not equal to zero
    firstindex  =   1_i8
    do n = 2_i8, size(nextpart_i8)
        if (nextpart_i8(n) .ne. 0_i8) then
            firstindex      =   n
            exit
        endif
    end do

    !decrease the first entry that is not equal to zero by one
    nextpart_i8(firstindex) =   nextpart_i8(firstindex)-1_i8

    !calculate the actual sum of the partition
    summ                    =   size(nextpart_i8)-firstindex-temp

    !fill up the partition below the index from above, until the sum equals n again       
    do n=firstindex-1_i8,1_i8,-1_i8
        temp                =   nextpart_i8(n)
        !set the next entry as big as possible
        nextpart_i8(n)      =   FLOOR(real(size(nextpart_i8)-summ,dp)/real(n,dp))
        !calculate the actual sum of the partition
        summ                =   summ + n*(nextpart_i8(n) - temp) 
        !if the sum is n we are done
        if (summ == size(nextpart_i8)) exit
    end do

  end function nextpart_i8

  ! ------------------------------------------------------------------

  function prevpart_i4(part)

    implicit none

    integer(i4), dimension(:), intent(in) :: part   
    integer(i4), dimension(size(part))    :: prevpart_i4
    
    integer(i4)                           :: n, summ


    !Check if input is valid
    if (.not. isgoodpart(part))     stop 'Not a valid input partition.'
    if (part(size(part)) .eq. 1_i4) stop 'Input was the first possible partition.'

    prevpart_i4             =   part    

    summ                    =   size(prevpart_i4)

    do n = 1_i4, size(prevpart_i4)
        if (summ + n .gt. size(prevpart_i4)) then
            if (prevpart_i4(n) .ne. 0_i4) then
                summ            =   summ - prevpart_i4(n)*n
                prevpart_i4(n)  =   0_i4
            endif
        else
            !increase the first entry that is increasable
            prevpart_i4(n)  =   prevpart_i4(n)+1_i4
            !correct the remainder by filling up 1's
            prevpart_i4(1)  =   size(prevpart_i4) - summ - n
            return
        endif
    end do

  end function prevpart_i4

  function prevpart_i8(part)

    implicit none

    integer(i8), dimension(:), intent(in) :: part   
    integer(i8), dimension(size(part))    :: prevpart_i8
    
    integer(i8)                           :: n, summ


    !Check if input is valid
    if (.not. isgoodpart(part))     stop 'Not a valid input partition.'
    if (part(size(part)) .eq. 1_i8) stop 'Input was the first possible partition.'

    prevpart_i8             =   part    

    summ                    =   size(prevpart_i8)

    do n = 1_i8, size(prevpart_i8)
        if (summ + n .gt. size(prevpart_i8)) then
            if (prevpart_i8(n) .ne. 0_i8) then
                summ            =   summ - prevpart_i8(n)*n
                prevpart_i8(n)  =   0_i8
            endif
        else
            !increase the first entry that is increasable
            prevpart_i8(n)  =   prevpart_i8(n)+1_i8
            !correct the remainder by filling up 1's
            prevpart_i8(1)  =   size(prevpart_i8) - summ - n
            return
        endif
    end do

  end function prevpart_i8

  ! ------------------------------------------------------------------

  function isgoodpart_i4(part)

    implicit none

    integer(i4), dimension(:), intent(in) :: part    
    logical                               :: isgoodpart_i4
    
    integer(i4)                           :: n, testn
    
    isgoodpart_i4   =   .false.
    testn           =   0_i4

    if (size(part) .le. 1_i4)    stop 'n must be at least 2'

    !check if p_i >= 0
    if (any(part<0_i4)) return

    !Check if it is a partition: sum(i*p_i)=n
    do n=size(part),1_i4,-1_i4
        testn       =   testn+n*part(n)
        if (testn > size(part)) return
    end do
    
    if (testn==size(part))  isgoodpart_i4  =   .true.
   
  end function isgoodpart_i4
  
  function isgoodpart_i8(part)

    implicit none

    integer(i8), dimension(:), intent(in) :: part    
    logical                               :: isgoodpart_i8
    
    integer(i8)                           :: n, testn
    
    isgoodpart_i8   =   .false.
    testn           =   0_i8

    if (size(part) .le. 1_i8)    stop 'n must be at least 2'

    !check if p_i >= 0
    if (any(part<0_i8)) return

    !Check if it is a partition: sum(i*p_i)=n
    do n=size(part),1_i8,-1_i8
        testn       =   testn+n*part(n)
        if (testn > size(part)) return
    end do
    
    if (testn==size(part))  isgoodpart_i8  =   .true.
   
  end function isgoodpart_i8
  
END MODULE mo_combinatorics

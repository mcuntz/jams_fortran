!> \file mo_sort.f90

!> \brief Quicksort array or indices

!> \details This module contains sorting routines from numerical recipes
!> 1. quicksort one array
!> 2. return indeces that would sort the array

!> \authors Matthias Cuntz - modified Numerical Recipes routines sort and indexx
!> \date March 2011

MODULE mo_sort

  ! This module contains sorting routines from numerical recipes:
  ! 1. quicksort one array
  ! 2. return indeces that would sort the array

  ! Usage:
  !   USE mo_sort, ONLY: sort
  !   call sort(arr)
  ! Or if you want to sort more than one array
  !   USE mo_sort, ONLY: sort_index
  !   integer(i4), dimension(size(arr)) :: ii
  !   ii = sort_index(arr)
  !   arr = arr(ii)

  ! Written March 2011, Matthias Cuntz
  !   - modified numerical recipes: sort and indexx
  ! Modified May 2016, Matthias Cuntz - use swap from mo_utils

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2011 Matthias Cuntz - mc (at) macu (dot) de
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

  ! Note on Numerical Recipes License
  ! ---------------------------------
  ! Be aware that some code is under the Numerical Recipes License 3rd
  ! edition <http://www.nr.com/aboutNR3license.html>

  ! The Numerical Recipes Personal Single-User License lets you personally
  ! use Numerical Recipes code ("the code") on any number of computers,
  ! but only one computer at a time. You are not permitted to allow anyone
  ! else to access or use the code. You may, under this license, transfer
  ! precompiled, executable applications incorporating the code to other,
  ! unlicensed, persons, providing that (i) the application is
  ! noncommercial (i.e., does not involve the selling or licensing of the
  ! application for a fee), and (ii) the application was first developed,
  ! compiled, and successfully run by you, and (iii) the code is bound
  ! into the application in such a manner that it cannot be accessed as
  ! individual routines and cannot practicably be unbound and used in
  ! other programs. That is, under this license, your application user
  ! must not be able to use Numerical Recipes code as part of a program
  ! library or "mix and match" workbench.

  ! Businesses and organizations that purchase the disk or code download,
  ! and that thus acquire one or more Numerical Recipes Personal
  ! Single-User Licenses, may permanently assign those licenses, in the
  ! number acquired, to individual employees. Such an assignment must be
  ! made before the code is first used and, once made, it is irrevocable
  ! and can not be transferred.

  ! If you do not hold a Numerical Recipes License, this code is only for
  ! informational and educational purposes but cannot be used.

  USE mo_kind,   ONLY: i4, sp, dp
  USE mo_utils,  ONLY: swap

  IMPLICIT NONE

  PRIVATE

#ifndef __PYTHON__
  PUBLIC :: sort       ! Quicksorts one array
  PUBLIC :: sort_index ! Returns indeces that sort the input array
#else
  PUBLIC :: nrsort       ! Quicksorts one array
  PUBLIC :: nrsort_index ! Returns indeces that sort the input array
#endif

  ! ------------------------------------------------------------------

  !     NAME
  !         sort

  !     PURPOSE
  !>        \brief Sorts the input array in ascending order.

  !     CALLING SEQUENCE
  !         call sort(arr)
  
  !     INDENT(IN)
  !         None

  !     INDENT(INOUT)
  !>        \param[in,out] "real(sp/dp) :: arr(:)"     On input:  unsorted 1D-array\n
  !>                                                   On output: ascendingly sorted input array

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         No mask or undefined value provided.

  !     EXAMPLE
  !         vec = (/ 3., 2, 1., 4., 6., 5. /)
  !         call sort(vec)
  !         -> see also example in test directory

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Matthias Cuntz - adapted routine sort from numerical recipes
  !>        \date Nov 2011
#ifndef __PYTHON__
  INTERFACE sort
     MODULE PROCEDURE isort_i4, isort_sp, isort_dp
  END INTERFACE sort
#else
  INTERFACE nrsort
     MODULE PROCEDURE isort_i4, isort_sp, isort_dp
  END INTERFACE nrsort
#endif

  ! ------------------------------------------------------------------

  !     NAME
  !         sort_index

  !     PURPOSE
  !>        \brief Gives the indeces that would sort an array in ascending order.

  !     CALLING SEQUENCE
  !         ii = sort_index(arr)
  
  !     INDENT(IN)
  !>        \param[in] "real(sp/dp) :: arr(:)"     Unsorted 1D-array

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RETURNS
  !>       \return integer(i4) :: indices(:) &mdash; Indices that would sort arr in ascending order

  !     RESTRICTIONS
  !         No mask or undefined value provided.

  !     EXAMPLE
  !         vec = (/ 3., 2, 1., 4., 6., 5. /)
  !         ii = sort_index(arr)
  !         arr = arr(ii)
  !         -> see also example in test directory

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Matthias Cuntz - adapted routine indexx from numerical recipes
  !>        \date Nov 2011
#ifndef __PYTHON__
  INTERFACE sort_index
     MODULE PROCEDURE isort_index_i4, isort_index_sp, isort_index_dp
  END INTERFACE sort_index
#else
  INTERFACE nrsort_index
     MODULE PROCEDURE isort_index_i4, isort_index_sp, isort_index_dp
  END INTERFACE nrsort_index
#endif

  ! Local routine for indeces exchange, from numerical recipes
  INTERFACE arth
     MODULE PROCEDURE arth_r, arth_d, arth_i
  END INTERFACE arth
  INTERFACE icomp_xchg
     MODULE PROCEDURE icomp_xchg_i4, icomp_xchg_sp, icomp_xchg_dp
  END INTERFACE icomp_xchg

  INTEGER(I4), PARAMETER :: NPAR_ARTH=16, NPAR2_ARTH=8

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  ! Array function returning an arithmetic progression.
  FUNCTION arth_r(first,increment,n)
    REAL(SP), INTENT(IN) :: first,increment
    INTEGER(I4), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: arth_r
    INTEGER(I4) :: k,k2
    REAL(SP) :: temp
    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_r(k)=arth_r(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_r

  FUNCTION arth_d(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER(I4) :: k,k2
    REAL(DP) :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_d(k)=arth_d(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_d(k)=arth_d(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_d

  FUNCTION arth_i(first,increment,n)
    INTEGER(I4), INTENT(IN) :: first,increment,n
    INTEGER(I4), DIMENSION(n) :: arth_i
    INTEGER(I4) :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_i

  ! ------------------------------------------------------------------

  ! From numerical recipes documentation
  ! Swap or don''t swap integer arguments, depending on the ordering of their corresponding
  ! elements in an array arr.
  SUBROUTINE icomp_xchg_i4(arr, i, j)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:), INTENT(IN)    :: arr
    INTEGER(i4),               INTENT(INOUT) :: i, j

    INTEGER(i4) :: swp

    if (arr(j) < arr(i)) then
       swp = i
       i   = j
       j   = swp
    end if

  END SUBROUTINE icomp_xchg_i4


  SUBROUTINE icomp_xchg_sp(arr, i, j)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:), INTENT(IN)    :: arr
    INTEGER(i4),               INTENT(INOUT) :: i, j

    INTEGER(i4) :: swp

    if (arr(j) < arr(i)) then
       swp = i
       i   = j
       j   = swp
    end if

  END SUBROUTINE icomp_xchg_sp


  SUBROUTINE icomp_xchg_dp(arr, i, j)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:), INTENT(IN)    :: arr
    INTEGER(i4),               INTENT(INOUT) :: i, j

    INTEGER(i4) :: swp

    if (arr(j) < arr(i)) then
       swp = i
       i   = j
       j   = swp
    end if

  END SUBROUTINE icomp_xchg_dp

! ------------------------------------------------------------------

  SUBROUTINE isort_i4(arr)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:), INTENT(INOUT) :: arr

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    INTEGER(i4) :: a
    INTEGER(i4) :: n, k, i, j, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n      = size(arr)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             a = arr(j)
             do i=j-1, l, -1
                if (arr(i) <= a) exit
                arr(i+1) = arr(i)
             end do
             arr(i+1) = a
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(arr(k),   arr(l+1))
          call swap(arr(l),   arr(r),   arr(l)>arr(r))
          call swap(arr(l+1), arr(r),   arr(l+1)>arr(r))
          call swap(arr(l),   arr(l+1), arr(l)>arr(l+1))
          i = l+1
          j = r
          a = arr(l+1)
          do
             do
                i = i+1
                if (arr(i) >= a) exit
             end do
             do
                j = j-1
                if (arr(j) <= a) exit
             end do
             if (j < i) exit
             call swap(arr(i), arr(j))
          end do
          arr(l+1) = arr(j)
          arr(j) = a
          jstack = jstack+2
          if (jstack > NSTACK) stop 'isort_sp: NSTACK too small'
          if (r-i+1 >= j-l) then
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
          else
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l = i
          end if
       end if
    end do

  END SUBROUTINE isort_i4

  ! ------------------------------------------------------------------

  ! quicksort
  SUBROUTINE isort_sp(arr)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(INOUT) :: arr

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    REAL(sp)    :: a
    INTEGER(i4) :: n, k, i, j, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n      = size(arr)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             a = arr(j)
             do i=j-1, l, -1
                if (arr(i) <= a) exit
                arr(i+1) = arr(i)
             end do
             arr(i+1) = a
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(arr(k),   arr(l+1))
          call swap(arr(l),   arr(r),   arr(l)>arr(r))
          call swap(arr(l+1), arr(r),   arr(l+1)>arr(r))
          call swap(arr(l),   arr(l+1), arr(l)>arr(l+1))
          i = l+1
          j = r
          a = arr(l+1)
          do
             do
                i = i+1
                if (arr(i) >= a) exit
             end do
             do
                j = j-1
                if (arr(j) <= a) exit
             end do
             if (j < i) exit
             call swap(arr(i), arr(j))
          end do
          arr(l+1) = arr(j)
          arr(j) = a
          jstack = jstack+2
          if (jstack > NSTACK) stop 'isort_sp: NSTACK too small'
          if (r-i+1 >= j-l) then
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
          else
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l = i
          end if
       end if
    end do

  END SUBROUTINE isort_sp


  SUBROUTINE isort_dp(arr)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(INOUT) :: arr

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    REAL(dp)    :: a
    INTEGER(i4) :: n, k, i, j, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n      = size(arr)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             a = arr(j)
             do i=j-1, l, -1
                if (arr(i) <= a) exit
                arr(i+1) = arr(i)
             end do
             arr(i+1) = a
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(arr(k),   arr(l+1))
          call swap(arr(l),   arr(r),   arr(l)>arr(r))
          call swap(arr(l+1), arr(r),   arr(l+1)>arr(r))
          call swap(arr(l),   arr(l+1), arr(l)>arr(l+1))
          i = l+1
          j = r
          a = arr(l+1)
          do
             do
                i = i+1
                if (arr(i) >= a) exit
             end do
             do
                j = j-1
                if (arr(j) <= a) exit
             end do
             if (j < i) exit
             call swap(arr(i), arr(j))
          end do
          arr(l+1) = arr(j)
          arr(j) = a
          jstack = jstack+2
          if (jstack > NSTACK) stop 'isort_dp: NSTACK too small'
          if (r-i+1 >= j-l) then
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
          else
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l = i
          end if
       end if
    end do

  END SUBROUTINE isort_dp


  ! ------------------------------------------------------------------

  ! quicksort but return indices
  FUNCTION isort_index_sp(arr)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:), INTENT(IN) :: arr
    INTEGER(i4), DIMENSION(size(arr))     :: isort_index_sp

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    REAL(sp)    :: a
    INTEGER(i4) :: n, k, i, j, indext, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n = size(arr)
    isort_index_sp = arth(1,1,n)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             indext = isort_index_sp(j)
             a = arr(indext)
             do i=j-1, l, -1
                if (arr(isort_index_sp(i)) <= a) exit
                isort_index_sp(i+1) = isort_index_sp(i)
             end do
             isort_index_sp(i+1) = indext
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(isort_index_sp(k), isort_index_sp(l+1))
          call icomp_xchg(arr, isort_index_sp(l),   isort_index_sp(r))
          call icomp_xchg(arr, isort_index_sp(l+1), isort_index_sp(r))
          call icomp_xchg(arr, isort_index_sp(l),   isort_index_sp(l+1))
          i = l+1
          j = r
          indext = isort_index_sp(l+1)
          a = arr(indext)
          do
             do
                i = i+1
                if (arr(isort_index_sp(i)) >= a) exit
             end do
             do
                j = j-1
                if (arr(isort_index_sp(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(isort_index_sp(i), isort_index_sp(j))
          end do
          isort_index_sp(l+1) = isort_index_sp(j)
          isort_index_sp(j)   = indext
          jstack = jstack+2
          if (jstack > NSTACK) stop 'isort_index_sp: NSTACK too small'
          if (r-i+1 >= j-l) then
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
          else
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l=i
          end if
       end if
    end do

  END FUNCTION isort_index_sp


  FUNCTION isort_index_dp(arr)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:), INTENT(IN) :: arr
    INTEGER(i4), DIMENSION(size(arr))     :: isort_index_dp

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    REAL(dp)    :: a
    INTEGER(i4) :: n, k, i, j, indext, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n = size(arr)
    isort_index_dp = arth(1,1,n)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             indext = isort_index_dp(j)
             a = arr(indext)
             do i=j-1, l, -1
                if (arr(isort_index_dp(i)) <= a) exit
                isort_index_dp(i+1) = isort_index_dp(i)
             end do
             isort_index_dp(i+1) = indext
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(isort_index_dp(k), isort_index_dp(l+1))
          call icomp_xchg(arr, isort_index_dp(l),   isort_index_dp(r))
          call icomp_xchg(arr, isort_index_dp(l+1), isort_index_dp(r))
          call icomp_xchg(arr, isort_index_dp(l),   isort_index_dp(l+1))
          i = l+1
          j = r
          indext = isort_index_dp(l+1)
          a = arr(indext)
          do
             do
                i = i+1
                if (arr(isort_index_dp(i)) >= a) exit
             end do
             do
                j = j-1
                if (arr(isort_index_dp(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(isort_index_dp(i), isort_index_dp(j))
          end do
          isort_index_dp(l+1) = isort_index_dp(j)
          isort_index_dp(j)   = indext
          jstack = jstack+2
          if (jstack > NSTACK) stop 'isort_index_dp: NSTACK too small'
          if (r-i+1 >= j-l) then
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
          else
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l=i
          end if
       end if
    end do

  END FUNCTION isort_index_dp


  FUNCTION isort_index_i4(iarr)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:), INTENT(IN) :: iarr
    INTEGER(i4), DIMENSION(size(iarr))    :: isort_index_i4

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    INTEGER(i4) :: a
    INTEGER(i4) :: n, k, i, j, indext, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n = size(iarr)
    isort_index_i4 = arth(1,1,n)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             indext = isort_index_i4(j)
             a = iarr(indext)
             ! Corrected from Numerical Recipes Forum
             ! do i=j-1, 1, -1
             do i=j-1, l, -1
                if (iarr(isort_index_i4(i)) <= a) exit
                isort_index_i4(i+1) = isort_index_i4(i)
             end do
             isort_index_i4(i+1) = indext
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(isort_index_i4(k), isort_index_i4(l+1))
          call icomp_xchg(iarr, isort_index_i4(l),   isort_index_i4(r))
          call icomp_xchg(iarr, isort_index_i4(l+1), isort_index_i4(r))
          call icomp_xchg(iarr, isort_index_i4(l),   isort_index_i4(l+1))
          i = l+1
          j = r
          indext = isort_index_i4(l+1)
          a = iarr(indext)
          do
             do
                i = i+1
                if (iarr(isort_index_i4(i)) >= a) exit
             end do
             do
                j = j-1
                if (iarr(isort_index_i4(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(isort_index_i4(i), isort_index_i4(j))
          end do
          isort_index_i4(l+1) = isort_index_i4(j)
          isort_index_i4(j)   = indext
          jstack = jstack+2
          if (jstack > NSTACK) stop 'isort_index_i4: NSTACK too small'
          if (r-i+1 >= j-l) then
             istack(jstack)   = r
             istack(jstack-1) = i
             r = j-1
          else
             istack(jstack)   = j-1
             istack(jstack-1) = l
             l = i
          end if
       end if
    end do

  END FUNCTION isort_index_i4

  ! ------------------------------------------------------------------

END MODULE mo_sort

MODULE mo_sort

  ! This module contains a sorting routines from numerical recipes:
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

  ! Copyright 2011 Matthias Cuntz


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
  USE mo_nrutil, ONLY: arth, swap

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sort       ! Quicksorts one array
  PUBLIC :: sort_index ! Returns indeces that sort the input array


  INTERFACE sort
     MODULE PROCEDURE sort_sp, sort_dp
  END INTERFACE sort
  INTERFACE sort_index
     MODULE PROCEDURE sort_index_i4, sort_index_sp, sort_index_dp
  END INTERFACE sort_index

  ! Local routine for indeces exchange, from numerical recipes
  INTERFACE icomp_xchg
     MODULE PROCEDURE icomp_xchg_i4, icomp_xchg_sp, icomp_xchg_dp
  END INTERFACE icomp_xchg

  ! ------------------------------------------------------------------

CONTAINS

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

  !     NAME
  !         sort

  !     PURPOSE
  !         Sorts the input array in ascending order.

  !     CALLING SEQUENCE
  !         call sort(arr)
  
  !     INDENT(IN)
  !         None

  !     INDENT(INOUT)
  !         real(sp/dp) :: arr(:)     On input:  unsorted 1D-array
  !                                   On output: ascendingly sorted input array

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
  !         Written,  Matthias Cuntz, Nov 2011 - adapted routine sort from numerical recipes

  SUBROUTINE sort_sp(arr)

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
          if (jstack > NSTACK) stop 'sort_sp: NSTACK too small'
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

  END SUBROUTINE sort_sp


  SUBROUTINE sort_dp(arr)

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
          if (jstack > NSTACK) stop 'sort_dp: NSTACK too small'
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

  END SUBROUTINE sort_dp


  ! ------------------------------------------------------------------

  !     NAME
  !         sort_index

  !     PURPOSE
  !         Gives the indeces that would sort an array in ascending order.

  !     CALLING SEQUENCE
  !         ii = sort_index(arr)
  
  !     INDENT(IN)
  !         real(sp/dp) :: arr(:)     Unsorted 1D-array

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         integer(i4)               Indeces that would sort arr in ascending order

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
  !         ii = sort_index(arr)
  !         arr = arr(ii)
  !         -> see also example in test directory

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written,  Matthias Cuntz, Nov 2011 - adapted routine indexx from numerical recipes

  FUNCTION sort_index_sp(arr)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:), INTENT(IN) :: arr
    INTEGER(i4), DIMENSION(size(arr))     :: sort_index_sp

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    REAL(sp)    :: a
    INTEGER(i4) :: n, k, i, j, indext, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n = size(arr)
    sort_index_sp = arth(1,1,n)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             indext = sort_index_sp(j)
             a = arr(indext)
             do i=j-1, l, -1
                if (arr(sort_index_sp(i)) <= a) exit
                sort_index_sp(i+1) = sort_index_sp(i)
             end do
             sort_index_sp(i+1) = indext
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(sort_index_sp(k), sort_index_sp(l+1))
          call icomp_xchg(arr, sort_index_sp(l),   sort_index_sp(r))
          call icomp_xchg(arr, sort_index_sp(l+1), sort_index_sp(r))
          call icomp_xchg(arr, sort_index_sp(l),   sort_index_sp(l+1))
          i = l+1
          j = r
          indext = sort_index_sp(l+1)
          a = arr(indext)
          do
             do
                i = i+1
                if (arr(sort_index_sp(i)) >= a) exit
             end do
             do
                j = j-1
                if (arr(sort_index_sp(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(sort_index_sp(i), sort_index_sp(j))
          end do
          sort_index_sp(l+1) = sort_index_sp(j)
          sort_index_sp(j)   = indext
          jstack = jstack+2
          if (jstack > NSTACK) stop 'sort_index_sp: NSTACK too small'
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

  END FUNCTION sort_index_sp


  FUNCTION sort_index_dp(arr)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:), INTENT(IN) :: arr
    INTEGER(i4), DIMENSION(size(arr))     :: sort_index_dp

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    REAL(dp)    :: a
    INTEGER(i4) :: n, k, i, j, indext, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n = size(arr)
    sort_index_dp = arth(1,1,n)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             indext = sort_index_dp(j)
             a = arr(indext)
             do i=j-1, l, -1
                if (arr(sort_index_dp(i)) <= a) exit
                sort_index_dp(i+1) = sort_index_dp(i)
             end do
             sort_index_dp(i+1) = indext
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(sort_index_dp(k), sort_index_dp(l+1))
          call icomp_xchg(arr, sort_index_dp(l),   sort_index_dp(r))
          call icomp_xchg(arr, sort_index_dp(l+1), sort_index_dp(r))
          call icomp_xchg(arr, sort_index_dp(l),   sort_index_dp(l+1))
          i = l+1
          j = r
          indext = sort_index_dp(l+1)
          a = arr(indext)
          do
             do
                i = i+1
                if (arr(sort_index_dp(i)) >= a) exit
             end do
             do
                j = j-1
                if (arr(sort_index_dp(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(sort_index_dp(i), sort_index_dp(j))
          end do
          sort_index_dp(l+1) = sort_index_dp(j)
          sort_index_dp(j)   = indext
          jstack = jstack+2
          if (jstack > NSTACK) stop 'sort_index_dp: NSTACK too small'
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

  END FUNCTION sort_index_dp


  FUNCTION sort_index_i4(iarr)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:), INTENT(IN) :: iarr
    INTEGER(i4), DIMENSION(size(iarr))    :: sort_index_i4

    INTEGER(i4), PARAMETER :: NN=15
    INTEGER(i4), PARAMETER :: NSTACK=50
    INTEGER(i4) :: a
    INTEGER(i4) :: n, k, i, j, indext, jstack, l, r
    INTEGER(i4), DIMENSION(NSTACK) :: istack

    n = size(iarr)
    sort_index_i4 = arth(1,1,n)
    jstack = 0
    l = 1
    r = n
    do
       if (r-l < NN) then
          do j=l+1, r
             indext = sort_index_i4(j)
             a = iarr(indext)
             ! Corrected from Numerical Recipes Forum
             ! do i=j-1, 1, -1
             do i=j-1, l, -1
                if (iarr(sort_index_i4(i)) <= a) exit
                sort_index_i4(i+1) = sort_index_i4(i)
             end do
             sort_index_i4(i+1) = indext
          end do
          if (jstack == 0) return
          r = istack(jstack)
          l = istack(jstack-1)
          jstack = jstack-2
       else
          k = (l+r)/2
          call swap(sort_index_i4(k), sort_index_i4(l+1))
          call icomp_xchg(iarr, sort_index_i4(l),   sort_index_i4(r))
          call icomp_xchg(iarr, sort_index_i4(l+1), sort_index_i4(r))
          call icomp_xchg(iarr, sort_index_i4(l),   sort_index_i4(l+1))
          i = l+1
          j = r
          indext = sort_index_i4(l+1)
          a = iarr(indext)
          do
             do
                i = i+1
                if (iarr(sort_index_i4(i)) >= a) exit
             end do
             do
                j = j-1
                if (iarr(sort_index_i4(j)) <= a) exit
             end do
             if (j < i) exit
             call swap(sort_index_i4(i), sort_index_i4(j))
          end do
          sort_index_i4(l+1) = sort_index_i4(j)
          sort_index_i4(j)   = indext
          jstack = jstack+2
          if (jstack > NSTACK) stop 'sort_index_i4: NSTACK too small'
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

  END FUNCTION sort_index_i4

  ! ------------------------------------------------------------------

END MODULE mo_sort

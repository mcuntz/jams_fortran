!> \file mo_quicksort.f90

!> \brief Different quicksort algorithms

!> \details This module contains different implementations of quicksort
!>          including an OpenMP version

!> \authors Matthias Cuntz
!> \date May 2014

MODULE mo_quicksort

  ! Module with different implementations of quicksort

  ! Written May 2014, Matthias Cuntz

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2014 Matthias Cuntz - mc (at) macu (dot) de
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

  USE mo_kind,  ONLY: i4, sp, dp
  use mo_utils, only: swap

  IMPLICIT NONE

  ! Aliases/wrapper
#ifndef __PYTHON__
  PUBLIC :: sort          ! qsort
  PUBLIC :: sort_index    ! qsort_index
  PUBLIC :: sortmp        ! qsortmp
  PUBLIC :: sortmp_index  ! wrapper for qsortmp
#endif

  ! Standard
  PUBLIC :: qsort         ! Non-recursive quicksort from Wikibooks
  PUBLIC :: qsort_index   ! Same but returning indexes
  PUBLIC :: qsortmp       ! Recursive algorithm with OpenMP support of David Bal
  PUBLIC :: qsortmp_index ! wrapper for qsortmp

  ! For comparison
  PUBLIC :: qsortc        ! Juli Rew's implementation of the recursive algorithm of Cormen et al. (1997)
  PUBLIC :: quick_sort    ! Alan Miller's implementation of the recursive algorithm of Brainerd et al. (1990)
  PUBLIC :: qsortb        ! Recursive algorithm of David Bal

  ! ------------------------------------------------------------------

  !     NAME
  !         qsort

  !     PURPOSE
  !         Non-recursive quicksort from Wikibooks
  !
  !>        \brief Quicksort.
  !
  !>        \details Non-recursive quicksort from Wikibooks.
  !>
  !>                 This implementation of quicksort is non-recursive, and choose as pivot element
  !>                 the median of three (the first, last and middle element of the list).
  !>                 It also uses insertion sort to sort lists with less than 10 elements.
  !>
  !>                 This implementation of Quicksort closely follows the one that can be found
  !>                 in the FORTRAN 90/95 GPL library AFNL (see references).
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !>        \param[inout] "real(sp/dp/i4) :: arr(:)"         The array to be sorted

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
  !>        \param[out] "integer(i4), optional :: indx(:)"   Index vector that would have sorted the array
  !
  !     RESTRICTIONS
  !         No mask or undefined value.
  !
  !     EXAMPLE
  !         vec  = (/ 3., 2., 1., 4., 6., 5. /)
  !         vec2 = (/ 1., 2., 3., 4., 6., 5. /)
  !         call qsort(vec, ii)
  !         vec2 = vec2(ii)
  !         -> see also example in test directory

  !     LITERATURE
  !         Wikibooks
  !             http://en.wikibooks.org/wiki/Algorithm_Implementation/Sorting/Quicksort#FORTRAN_90.2F95
  !         A Fortran Numerical Library (AFNL)
  !             http://sourceforge.net/projects/afnl

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date May 2014
  INTERFACE qsort
     MODULE PROCEDURE qsort_dp, qsort_sp, qsort_i4
  END INTERFACE qsort

  ! ------------------------------------------------------------------

  !     NAME
  !         qsort_index

  !     PURPOSE
  !         Non-recursive quicksort from Wikibooks returning the indexes to sort array.
  !
  !>        \brief Quicksort on indexes.
  !
  !>        \details Non-recursive quicksort from Wikibooks returning the indexes to sort the array.
  !>
  !>                 This implementation of quicksort is non-recursive, and choose as pivot element
  !>                 the median of three (the first, last and middle element of the list).
  !>                 It also uses insertion sort to sort lists with less than 10 elements.
  !>
  !>                 This implementation of Quicksort closely follows the one that can be found
  !>                 in the FORTRAN 90/95 GPL library AFNL (see references).
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp/i4) :: arr(:)"    The array which would be sorted given the returned indexes
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
  !     RETURNS
  !>       \return integer(i4) :: indixes(:) &mdash; Indexes that sort arr in ascending order
  !
  !     RESTRICTIONS
  !         No mask or undefined value.
  !
  !     EXAMPLE
  !         vec = (/ 3., 2, 1., 4., 6., 5. /)
  !         ii = qsort_index(vec)
  !         vec = vec(ii)
  !         -> see also example in test directory

  !     LITERATURE
  !         Wikibooks
  !             http://en.wikibooks.org/wiki/Algorithm_Implementation/Sorting/Quicksort#FORTRAN_90.2F95
  !         A Fortran Numerical Library (AFNL)
  !             http://sourceforge.net/projects/afnl

  !     HISTORY
  !>        \author  Matthias Cuntz - adapted qsort
  !>        \date May 2014
  !         modified Stephan Thober - added if clause when given array has size 1
  !                                 - modified documentation in subroutine
  INTERFACE qsort_index
     MODULE PROCEDURE qsort_index_dp, qsort_index_sp, qsort_index_i4
  END INTERFACE qsort_index

  ! ------------------------------------------------------------------

  !     NAME
  !         qsortmp

  !     PURPOSE
  !         Recursive quicksort of David Bal with OpenMP support
  !
  !>        \brief OpenMP Quicksort.
  !
  !>        \details Recursive quicksort of David Bal with OpenMP support.
  !>
  !>                 The input array is divided into groups and groups are sorted in different threads.
  !>                 The sorted groups are merged together using a bottom-up merge sort approach.
  !>                 Neighboring groups are merged as soon as they are finished. Also a multi-threaded merge
  !>                 is used to speed up the final stages of the merge.
  !>
  !>                 The individual quicksorts use the median of 9 points to calculate the pivot element.
  !>                 It also uses insertion sort below 50 elements.
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !>        \param[inout] "real(sp/dp/i4) :: arr(:)"         The array to be sorted

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4), optional :: nt"         Number of threads to use. Default: all available threads.
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "integer(i4), optional :: indx(:)"   Index vector that would have sorted the array
  !
  !     RESTRICTIONS
  !         No mask or undefined value.
  !
  !     EXAMPLE
  !         vec  = (/ 3., 2., 1., 4., 6., 5. /)
  !         vec2 = (/ 1., 2., 3., 4., 6., 5. /)
  !         call qsortmp(vec, ii)
  !         vec2 = vec2(ii)
  !         -> see also example in test directory

  !     LITERATURE
  !         David Bal
  !             Quicksort - http://balfortran.org/
  !             Multithreaded Quicksort - http://balfortran.org/multithreaded_sort.html

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date May 2014
  INTERFACE qsortmp
     MODULE PROCEDURE qsortmp_dp, qsortmp_sp, qsortmp_i4
  END INTERFACE qsortmp

  ! ------------------------------------------------------------------

  !     NAME
  !         sortmp_index

  !     PURPOSE
  !>        \brief Gives the indeces that would sort an array in ascending order.
  !
  !>        \details Wrapper for recursive quicksort of David Bal with OpenMP support.

  !     CALLING SEQUENCE
  !         ii = sortmp_index(arr)
  
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
  !         ii = sortmp_index(arr)
  !         arr = arr(ii)
  !         -> see also example in test directory

  !     LITERATURE
  !         David Bal
  !             Quicksort - http://balfortran.org/
  !             Multithreaded Quicksort - http://balfortran.org/multithreaded_sort.html

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date May 2014
  INTERFACE qsortmp_index
     MODULE PROCEDURE sortmp_index_i4, sortmp_index_sp, sortmp_index_dp
  END INTERFACE qsortmp_index
  
  ! aliases
#ifndef __PYTHON__
  INTERFACE sort
     MODULE PROCEDURE qsort_dp, qsort_sp, qsort_i4
  END INTERFACE sort
  INTERFACE sort_index
     MODULE PROCEDURE qsort_index_dp, qsort_index_sp, qsort_index_i4
  END INTERFACE sort_index
  INTERFACE sortmp
     MODULE PROCEDURE qsortmp_dp, qsortmp_sp, qsortmp_i4
  END INTERFACE sortmp
  INTERFACE sortmp_index
     MODULE PROCEDURE sortmp_index_i4, sortmp_index_sp, sortmp_index_dp
  END INTERFACE sortmp_index
#endif

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION sortmp_index_dp(arr)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:), INTENT(IN) :: arr
    INTEGER(i4), DIMENSION(size(arr))     :: sortmp_index_dp

    REAL(dp), DIMENSION(size(arr)) :: iarr

    iarr = arr
    call qsortmp(iarr, sortmp_index_dp)

  END FUNCTION sortmp_index_dp

  FUNCTION sortmp_index_sp(arr)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:), INTENT(IN) :: arr
    INTEGER(i4), DIMENSION(size(arr))     :: sortmp_index_sp

    REAL(sp), DIMENSION(size(arr)) :: iarr

    iarr = arr
    call qsortmp(iarr, sortmp_index_sp)

  END FUNCTION sortmp_index_sp

  FUNCTION sortmp_index_i4(arr)

    IMPLICIT NONE

    integer(i4), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(i4), DIMENSION(size(arr))     :: sortmp_index_i4

    integer(i4), DIMENSION(size(arr)) :: iarr

    iarr = arr
    call qsortmp(iarr, sortmp_index_i4)

  END FUNCTION sortmp_index_i4

  ! ------------------------------------------------------------------
  ! http://en.wikibooks.org/wiki/Algorithm_Implementation/Sorting/Quicksort#FORTRAN_90.2F95

  Subroutine qsort_dp(X, Ipt)

    ! ***********************************
    ! * Sort Array X(:) in ascendent order
    ! * If present Ipt, a pointer with the
    ! * changes is returned in Ipt.
    ! ***********************************

    IMPLICIT NONE

    Type Limits
       Integer(i4) :: Ileft, Iright
    End Type Limits

    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer(i4), Parameter :: Isw = 10

    Real(dp), Intent (inout) :: X(:)
    Integer(i4), Intent (out), Optional :: Ipt(:)

    Integer(i4) :: I, Ipvn, Ileft, Iright, ISpos, ISmax
    Integer(i4), Allocatable :: IIpt(:)
    Type(Limits), Allocatable :: Stack(:)

    Allocate(Stack(Size(X)))

    Stack(:)%Ileft = 0
    If (Present(Ipt)) Then
       Forall (I=1:Size(Ipt)) Ipt(I) = I

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)

       Do While (Stack(ISpos)%Ileft /= 0)

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL Insrtlc_Dp(X, Ipt, Ileft,Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = Choosepiv_Dp(X, Ileft, Iright)
             Ipvn = Partition_Dp(X, Ileft, Iright, Ipvn, Ipt)

             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1)%Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do

    Else

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)

       !Allocate(IIpt(10))
       Allocate(IIpt(size(X)))
       Do While (Stack(ISpos)%Ileft /= 0)

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL Insrtlc_Dp(X, IIpt, Ileft, Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = Choosepiv_Dp(X, Ileft, Iright)
             Ipvn = Partition_Dp(X, Ileft, Iright, Ipvn)

             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1)%Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do
       Deallocate(IIpt)

    End If

    Deallocate(Stack)

    Return

  CONTAINS

    ! ***********************************

    Integer(i4) Function Choosepiv_Dp(XX, IIleft, IIright) Result (IIpv)

      ! ***********************************
      ! * Choose a Pivot element from XX(Ileft:Iright)
      ! * for Qsort. This routine chooses the median
      ! * of the first, last and mid element of the
      ! * list.
      ! ***********************************

      IMPLICIT NONE

      Real(dp), Intent (in) :: XX(:)
      Integer(i4), Intent (in) :: IIleft, IIright

      Real(dp) :: XXcp(3)
      Integer(i4) :: IIpt(3), IImd

      IImd = Int((IIleft+IIright)/2)
      XXcp(1) = XX(IIleft)
      XXcp(2) = XX(IImd)
      XXcp(3) = XX(IIright)
      IIpt = (/1,2,3/)

      CALL Insrtlc_Dp(XXcp, IIpt, 1, 3)

      Select Case (IIpt(2))
      Case (1)
         IIpv = IIleft
      Case (2)
         IIpv = IImd
      Case (3)
         IIpv = IIright
      End Select

      Return

    End Function Choosepiv_Dp

    ! ***********************************

    Subroutine Insrtlc_Dp(XX, IIpt, IIl, IIr)

      ! ***********************************
      ! * Perform an insertion sort of the list
      ! * XX(:) between index values IIl and IIr.
      ! * IIpt(:) returns the permutations
      ! * made to sort.
      ! ***********************************

      IMPLICIT NONE

      Real(dp), Intent(inout) :: XX(:)
      Integer(i4), Intent(inout) :: IIpt(:)
      Integer(i4), Intent(in) :: IIl, IIr

      Real(dp) :: RRtmp
      Integer(i4) :: II, JJ


      Do II = IIl+1, IIr
         RRtmp = XX(II)
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(JJ)) Then
               XX(JJ+1) = XX(JJ)
               CALL swap(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
         XX(JJ+1) = RRtmp
      End Do

      Return

    End Subroutine Insrtlc_Dp

  End Subroutine qsort_dp


  Integer Function Partition_Dp(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)

    ! ***********************************
    ! * This routine arranges the array X
    ! * between the index values Ileft and Iright
    ! * positioning elements smallers than
    ! * X(Ipv) at the left and the others
    ! * at the right.
    ! * Internal routine used by Qsort.
    ! ***********************************

    IMPLICIT NONE

    Real(dp), Intent (inout) :: X(:)
    Integer(i4), Intent (in) :: Ileft, Iright, Ipv
    Integer(i4), Intent (inout), Optional :: Ipt(:)

    Real(dp) :: Rpv
    Integer(i4) :: I

    Rpv = X(Ipv)
    CALL swap(X, Ipv, Iright)
    If (Present(Ipt)) CALL swap(Ipt, Ipv, Iright)
    Ipvfn = Ileft

    If (Present(Ipt))  Then
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL swap(X, I, Ipvfn)
             CALL swap(Ipt, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    Else
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL swap(X, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    End If

    CALL swap(X, Ipvfn, Iright)
    If (Present(Ipt)) CALL swap(Ipt, Ipvfn, Iright)

    Return

  End Function Partition_Dp


  Subroutine qsort_sp(X, Ipt)

    ! ***********************************
    ! * Sort Array X(:) in ascendent order
    ! * If present Ipt, a pointer with the
    ! * changes is returned in Ipt.
    ! ***********************************

    IMPLICIT NONE

    Type Limits
       Integer(i4) :: Ileft, Iright
    End Type Limits

    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer(i4), Parameter :: Isw = 10

    Real(sp), Intent (inout) :: X(:)
    Integer(i4), Intent (out), Optional :: Ipt(:)

    Integer(i4) :: I, Ipvn, Ileft, Iright, ISpos, ISmax
    Integer(i4), Allocatable :: IIpt(:)
    Type(Limits), Allocatable :: Stack(:)

    Allocate(Stack(Size(X)))

    Stack(:)%Ileft = 0
    If (Present(Ipt)) Then
       Forall (I=1:Size(Ipt)) Ipt(I) = I

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)

       Do While (Stack(ISpos)%Ileft /= 0)

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL Insrtlc_Sp(X, Ipt, Ileft,Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = Choosepiv_Sp(X, Ileft, Iright)
             Ipvn = Partition_Sp(X, Ileft, Iright, Ipvn, Ipt)

             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1)%Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do

    Else

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)

       !Allocate(IIpt(10))
       Allocate(IIpt(size(X)))
       Do While (Stack(ISpos)%Ileft /= 0)

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL Insrtlc_Sp(X, IIpt, Ileft, Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = Choosepiv_Sp(X, Ileft, Iright)
             Ipvn = Partition_Sp(X, Ileft, Iright, Ipvn)

             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1)%Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do
       Deallocate(IIpt)

    End If

    Deallocate(Stack)

    Return

  CONTAINS

    ! ***********************************

    Integer(i4) Function Choosepiv_Sp(XX, IIleft, IIright) Result (IIpv)

      ! ***********************************
      ! * Choose a Pivot element from XX(Ileft:Iright)
      ! * for Qsort. This routine chooses the median
      ! * of the first, last and mid element of the
      ! * list.
      ! ***********************************

      IMPLICIT NONE

      Real(sp), Intent (in) :: XX(:)
      Integer(i4), Intent (in) :: IIleft, IIright

      Real(sp) :: XXcp(3)
      Integer(i4) :: IIpt(3), IImd

      IImd = Int((IIleft+IIright)/2)
      XXcp(1) = XX(IIleft)
      XXcp(2) = XX(IImd)
      XXcp(3) = XX(IIright)
      IIpt = (/1,2,3/)

      CALL Insrtlc_Sp(XXcp, IIpt, 1, 3)

      Select Case (IIpt(2))
      Case (1)
         IIpv = IIleft
      Case (2)
         IIpv = IImd
      Case (3)
         IIpv = IIright
      End Select

      Return

    End Function Choosepiv_Sp

    ! ***********************************

    Subroutine Insrtlc_Sp(XX, IIpt, IIl, IIr)

      ! ***********************************
      ! * Perform an insertion sort of the list
      ! * XX(:) between index values IIl and IIr.
      ! * IIpt(:) returns the permutations
      ! * made to sort.
      ! ***********************************

      IMPLICIT NONE

      Real(sp), Intent(inout) :: XX(:)
      Integer(i4), Intent(inout) :: IIpt(:)
      Integer(i4), Intent(in) :: IIl, IIr

      Real(sp) :: RRtmp
      Integer(i4) :: II, JJ


      Do II = IIl+1, IIr
         RRtmp = XX(II)
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(JJ)) Then
               XX(JJ+1) = XX(JJ)
               CALL swap(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
         XX(JJ+1) = RRtmp
      End Do

      Return

    End Subroutine Insrtlc_Sp

  End Subroutine qsort_sp


  Integer Function Partition_Sp(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)

    ! ***********************************
    ! * This routine arranges the array X
    ! * between the index values Ileft and Iright
    ! * positioning elements smallers than
    ! * X(Ipv) at the left and the others
    ! * at the right.
    ! * Internal routine used by Qsort.
    ! ***********************************

    IMPLICIT NONE

    Real(sp), Intent (inout) :: X(:)
    Integer(i4), Intent (in) :: Ileft, Iright, Ipv
    Integer(i4), Intent (inout), Optional :: Ipt(:)

    Real(sp) :: Rpv
    Integer(i4) :: I

    Rpv = X(Ipv)
    CALL swap(X, Ipv, Iright)
    If (Present(Ipt)) CALL swap(Ipt, Ipv, Iright)
    Ipvfn = Ileft

    If (Present(Ipt))  Then
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL swap(X, I, Ipvfn)
             CALL swap(Ipt, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    Else
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL swap(X, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    End If

    CALL swap(X, Ipvfn, Iright)
    If (Present(Ipt)) CALL swap(Ipt, Ipvfn, Iright)

    Return

  End Function Partition_Sp


  Subroutine qsort_i4(X, Ipt)

    ! ***********************************
    ! * Sort Array X(:) in ascendent order
    ! * If present Ipt, a pointer with the
    ! * changes is returned in Ipt.
    ! ***********************************

    IMPLICIT NONE

    Type Limits
       Integer(i4) :: Ileft, Iright
    End Type Limits

    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer(i4), Parameter :: Isw = 10

    integer(i4), Intent (inout) :: X(:)
    Integer(i4), Intent (out), Optional :: Ipt(:)

    Integer(i4) :: I, Ipvn, Ileft, Iright, ISpos, ISmax
    Integer(i4), Allocatable :: IIpt(:)
    Type(Limits), Allocatable :: Stack(:)

    Allocate(Stack(Size(X)))

    Stack(:)%Ileft = 0
    If (Present(Ipt)) Then
       Forall (I=1:Size(Ipt)) Ipt(I) = I

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)

       Do While (Stack(ISpos)%Ileft /= 0)

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL Insrtlc_I4(X, Ipt, Ileft,Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = Choosepiv_I4(X, Ileft, Iright)
             Ipvn = Partition_I4(X, Ileft, Iright, Ipvn, Ipt)

             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1)%Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do

    Else

       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)

       !Allocate(IIpt(10))
       Allocate(IIpt(size(X)))
       Do While (Stack(ISpos)%Ileft /= 0)

          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL Insrtlc_I4(X, IIpt, Ileft, Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = Choosepiv_I4(X, Ileft, Iright)
             Ipvn = Partition_I4(X, Ileft, Iright, Ipvn)

             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1)%Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do
       Deallocate(IIpt)

    End If

    Deallocate(Stack)

    Return

  CONTAINS

    ! ***********************************

    Integer(i4) Function Choosepiv_I4(XX, IIleft, IIright) Result (IIpv)

      ! ***********************************
      ! * Choose a Pivot element from XX(Ileft:Iright)
      ! * for Qsort. This routine chooses the median
      ! * of the first, last and mid element of the
      ! * list.
      ! ***********************************

      IMPLICIT NONE

      integer(i4), Intent (in) :: XX(:)
      Integer(i4), Intent (in) :: IIleft, IIright

      integer(i4) :: XXcp(3)
      Integer(i4) :: IIpt(3), IImd

      IImd = Int((IIleft+IIright)/2)
      XXcp(1) = XX(IIleft)
      XXcp(2) = XX(IImd)
      XXcp(3) = XX(IIright)
      IIpt = (/1,2,3/)

      CALL Insrtlc_I4(XXcp, IIpt, 1, 3)

      Select Case (IIpt(2))
      Case (1)
         IIpv = IIleft
      Case (2)
         IIpv = IImd
      Case (3)
         IIpv = IIright
      End Select

      Return

    End Function Choosepiv_I4

    ! ***********************************

    Subroutine Insrtlc_I4(XX, IIpt, IIl, IIr)

      ! ***********************************
      ! * Perform an insertion sort of the list
      ! * XX(:) between index values IIl and IIr.
      ! * IIpt(:) returns the permutations
      ! * made to sort.
      ! ***********************************

      IMPLICIT NONE

      integer(i4), Intent(inout) :: XX(:)
      Integer(i4), Intent(inout) :: IIpt(:)
      Integer(i4), Intent(in) :: IIl, IIr

      integer(i4) :: RRtmp
      Integer(i4) :: II, JJ


      Do II = IIl+1, IIr
         RRtmp = XX(II)
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(JJ)) Then
               XX(JJ+1) = XX(JJ)
               CALL swap(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
         XX(JJ+1) = RRtmp
      End Do

      Return

    End Subroutine Insrtlc_I4

  End Subroutine qsort_i4


  Integer Function Partition_I4(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)

    ! ***********************************
    ! * This routine arranges the array X
    ! * between the index values Ileft and Iright
    ! * positioning elements smallers than
    ! * X(Ipv) at the left and the others
    ! * at the right.
    ! * Internal routine used by Qsort.
    ! ***********************************

    IMPLICIT NONE

    integer(i4), Intent (inout) :: X(:)
    Integer(i4), Intent (in) :: Ileft, Iright, Ipv
    Integer(i4), Intent (inout), Optional :: Ipt(:)

    integer(i4) :: Rpv
    Integer(i4) :: I

    Rpv = X(Ipv)
    CALL swap(X, Ipv, Iright)
    If (Present(Ipt)) CALL swap(Ipt, Ipv, Iright)
    Ipvfn = Ileft

    If (Present(Ipt))  Then
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL swap(X, I, Ipvfn)
             CALL swap(Ipt, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    Else
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL swap(X, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    End If

    CALL swap(X, Ipvfn, Iright)
    If (Present(Ipt)) CALL swap(Ipt, Ipvfn, Iright)

    Return

  End Function Partition_I4

  ! ------------------------------------------------------------------

  function Qsort_index_dp(X) result(Ipt)

    ! ***********************************
    ! * Sort Array X(:) in ascendent order
    ! * and return permutation in Ipt
    ! ***********************************

    IMPLICIT NONE

    Type Limits
       Integer(i4) :: Ileft, Iright
    End Type Limits

    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer(i4), Parameter :: Isw = 10

    Real(dp),    dimension(:),      Intent(in) :: X
    Integer(i4), dimension(size(X))            :: Ipt

    Integer(i4) :: I, Ipvn, Ileft, Iright, ISpos, ISmax, Imd
    Type(Limits), Allocatable :: Stack(:)
    Real(dp)    :: XXcp(3)
    Integer(i4) :: IIpt(3)

    ! dont sort if size equal to one
    if ( size(X) .eq. 1 ) Then
       Ipt = 1
       return
    end if

    Allocate(Stack(Size(X)))

    Stack(:)%Ileft = 0
    Forall (I=1:Size(Ipt)) Ipt(I) = I

    ! Iniitialize the stack
    Ispos = 1
    Ismax = 1
    Stack(ISpos)%Ileft  = 1
    Stack(ISpos)%Iright = Size(X)

    Do While (Stack(ISpos)%Ileft /= 0)

       Ileft = Stack(ISPos)%Ileft
       Iright = Stack(ISPos)%Iright
       If (Iright-Ileft <= Isw) Then
          CALL InsrtLC_index_dp(X, Ipt, Ileft, Iright)
          ISpos = ISPos + 1
       Else
          ! Ipvn = ChoosePiv_index_dp(X(Ipt), Ileft, Iright)
          ! Choose Pivot directly because copy of X takes long
          Imd = Ileft + (Iright-Ileft)/2
          XXcp(1) = X(Ipt(Ileft))
          XXcp(2) = X(Ipt(Imd))
          XXcp(3) = X(Ipt(Iright))
          IIpt = (/1,2,3/)
          CALL InsrtLC_index_dp(XXcp, IIpt, 1, 3)
          Select Case (IIpt(2))
          Case (1)
             Ipvn = Ileft
          Case (2)
             Ipvn = Imd
          Case (3)
             Ipvn = Iright
          End Select
          Ipvn = Partition_index_dp(X, Ileft, Iright, Ipvn, Ipt)

          Stack(ISmax+1)%Ileft = Ileft
          Stack(ISmax+1)%Iright = Ipvn-1
          Stack(ISmax+2)%Ileft = Ipvn + 1
          Stack(ISmax+2)%Iright = Iright
          ISpos = ISpos + 1
          ISmax = ISmax + 2
       End If
    End Do

    Deallocate(Stack)

    Return

  CONTAINS

    ! ***********************************

    Subroutine InsrtLC_index_dp(XX, IIpt, IIl, IIr)

      ! ***********************************
      ! * Perform an insertion sort of the list
      ! * XX(:) between index values IIl and IIr
      ! * but returns only IIpt(:) the permutations
      ! * made to sort.
      ! ***********************************

      IMPLICIT NONE

      Real(dp), Intent(in) :: XX(:)
      Integer(i4), Intent(inout) :: IIpt(:)
      Integer(i4), Intent(in) :: IIl, IIr

      Real(dp) :: RRtmp
      Integer(i4) :: II, JJ


      Do II = IIl+1, IIr
         RRtmp = XX(IIpt(II))
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(IIpt(JJ))) Then
               CALL swap(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
      End Do

      Return

    End Subroutine InsrtLC_index_dp

  end function Qsort_index_dp

  ! ***********************************

  Integer(i4) Function Partition_index_dp(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)

    ! ***********************************
    ! * This routine arranges the indexes Ipt of the array X
    ! * between the index values Ileft and Iright
    ! * positioning elements smallers than
    ! * X(Ipv) at the left and the others
    ! * at the right.
    ! * Internal routine used by Qsort.
    ! ***********************************

    IMPLICIT NONE

    Real(dp), Intent (in) :: X(:)
    Integer(i4), Intent (in) :: Ileft, Iright, Ipv
    Integer(i4), Intent (inout) :: Ipt(:)

    Real(dp) :: Rpv
    Integer(i4) :: I

    Rpv = X(Ipt(Ipv))
    CALL swap(Ipt, Ipv, Iright)
    Ipvfn = Ileft

    Do I = Ileft, Iright-1
       If (X(Ipt(I)) <= Rpv) Then
          CALL swap(Ipt, I, Ipvfn)
          Ipvfn = Ipvfn + 1
       End If
    End Do
    CALL swap(Ipt, Ipvfn, Iright)

    Return

  End Function Partition_index_dp



  function Qsort_index_sp(X) result(Ipt)

    ! ***********************************
    ! * Sort Array X(:) in ascendent order
    ! * and return permutation in Ipt
    ! ***********************************

    IMPLICIT NONE

    Type Limits
       Integer(i4) :: Ileft, Iright
    End Type Limits

    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer(i4), Parameter :: Isw = 10

    Real(sp),    dimension(:),      Intent(in) :: X
    Integer(i4), dimension(size(X))            :: Ipt

    Integer(i4) :: I, Ipvn, Ileft, Iright, ISpos, ISmax, Imd
    Type(Limits), Allocatable :: Stack(:)
    Real(sp)    :: XXcp(3)
    Integer(i4) :: IIpt(3)

    ! dont sort if size equal to one
    if ( size(X) .eq. 1 ) Then
       Ipt = 1
       Return
    end if

    Allocate(Stack(Size(X)))

    Stack(:)%Ileft = 0
    Forall (I=1:Size(Ipt)) Ipt(I) = I

    ! Iniitialize the stack
    Ispos = 1
    Ismax = 1
    Stack(ISpos)%Ileft  = 1
    Stack(ISpos)%Iright = Size(X)

    Do While (Stack(ISpos)%Ileft /= 0)

       Ileft = Stack(ISPos)%Ileft
       Iright = Stack(ISPos)%Iright
       If (Iright-Ileft <= Isw) Then
          CALL InsrtLC_index_sp(X, Ipt, Ileft, Iright)
          ISpos = ISPos + 1
       Else
          ! Ipvn = ChoosePiv_index_sp(X(Ipt), Ileft, Iright)
          ! Choose Pivot directly because copy of X takes long
          Imd = Ileft + (Iright-Ileft)/2
          XXcp(1) = X(Ipt(Ileft))
          XXcp(2) = X(Ipt(Imd))
          XXcp(3) = X(Ipt(Iright))
          IIpt = (/1,2,3/)
          CALL InsrtLC_index_sp(XXcp, IIpt, 1, 3)
          Select Case (IIpt(2))
          Case (1)
             Ipvn = Ileft
          Case (2)
             Ipvn = Imd
          Case (3)
             Ipvn = Iright
          End Select
          Ipvn = Partition_index_sp(X, Ileft, Iright, Ipvn, Ipt)

          Stack(ISmax+1)%Ileft = Ileft
          Stack(ISmax+1)%Iright = Ipvn-1
          Stack(ISmax+2)%Ileft = Ipvn + 1
          Stack(ISmax+2)%Iright = Iright
          ISpos = ISpos + 1
          ISmax = ISmax + 2
       End If
    End Do

    Deallocate(Stack)

    Return

  CONTAINS

    ! ***********************************

    Subroutine InsrtLC_index_sp(XX, IIpt, IIl, IIr)

      ! ***********************************
      ! * Perform an insertion sort of the list
      ! * XX(:) between index values IIl and IIr
      ! * but returns only IIpt(:) the permutations
      ! * made to sort.
      ! ***********************************

      IMPLICIT NONE

      Real(sp), Intent(in) :: XX(:)
      Integer(i4), Intent(inout) :: IIpt(:)
      Integer(i4), Intent(in) :: IIl, IIr

      Real(sp) :: RRtmp
      Integer(i4) :: II, JJ


      Do II = IIl+1, IIr
         RRtmp = XX(IIpt(II))
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(IIpt(JJ))) Then
               CALL swap(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
      End Do

      Return

    End Subroutine InsrtLC_index_sp

  end function Qsort_index_sp

  ! ***********************************

  Integer(i4) Function Partition_index_sp(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)

    ! ***********************************
    ! * This routine arranges the indexes Ipt of the array X
    ! * between the index values Ileft and Iright
    ! * positioning elements smallers than
    ! * X(Ipv) at the left and the others
    ! * at the right.
    ! * Internal routine used by Qsort.
    ! ***********************************

    IMPLICIT NONE

    Real(sp), Intent (in) :: X(:)
    Integer(i4), Intent (in) :: Ileft, Iright, Ipv
    Integer(i4), Intent (inout) :: Ipt(:)

    Real(sp) :: Rpv
    Integer(i4) :: I

    Rpv = X(Ipt(Ipv))
    CALL swap(Ipt, Ipv, Iright)
    Ipvfn = Ileft

    Do I = Ileft, Iright-1
       If (X(Ipt(I)) <= Rpv) Then
          CALL swap(Ipt, I, Ipvfn)
          Ipvfn = Ipvfn + 1
       End If
    End Do
    CALL swap(Ipt, Ipvfn, Iright)

    Return

  End Function Partition_index_sp


  function Qsort_index_i4(X) result(Ipt)

    ! ***********************************
    ! * Sort Array X(:) in ascendent order
    ! * and return permutation in Ipt
    ! ***********************************

    IMPLICIT NONE

    Type Limits
       Integer(i4) :: Ileft, Iright
    End Type Limits

    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer(i4), Parameter :: Isw = 10

    integer(i4),    dimension(:),      Intent(in) :: X
    Integer(i4), dimension(size(X))            :: Ipt

    Integer(i4) :: I, Ipvn, Ileft, Iright, ISpos, ISmax, Imd
    Type(Limits), Allocatable :: Stack(:)
    integer(i4)    :: XXcp(3)
    Integer(i4) :: IIpt(3)

    ! dont sort if size equal to one
    if ( size(X) .eq. 1 ) Then
       Ipt = 1
       Return
    end if

    Allocate(Stack(Size(X)))

    Stack(:)%Ileft = 0
    Forall (I=1:Size(Ipt)) Ipt(I) = I

    ! Iniitialize the stack
    Ispos = 1
    Ismax = 1
    Stack(ISpos)%Ileft  = 1
    Stack(ISpos)%Iright = Size(X)

    Do While (Stack(ISpos)%Ileft /= 0)

       Ileft = Stack(ISPos)%Ileft
       Iright = Stack(ISPos)%Iright
       If (Iright-Ileft <= Isw) Then
          CALL InsrtLC_index_i4(X, Ipt, Ileft, Iright)
          ISpos = ISPos + 1
       Else
          ! Ipvn = ChoosePiv_index_i4(X(Ipt), Ileft, Iright)
          ! Choose Pivot directly because copy of X takes long
          Imd = Ileft + (Iright-Ileft)/2
          XXcp(1) = X(Ipt(Ileft))
          XXcp(2) = X(Ipt(Imd))
          XXcp(3) = X(Ipt(Iright))
          IIpt = (/1,2,3/)
          CALL InsrtLC_index_i4(XXcp, IIpt, 1, 3)
          Select Case (IIpt(2))
          Case (1)
             Ipvn = Ileft
          Case (2)
             Ipvn = Imd
          Case (3)
             Ipvn = Iright
          End Select
          Ipvn = Partition_index_i4(X, Ileft, Iright, Ipvn, Ipt)

          Stack(ISmax+1)%Ileft = Ileft
          Stack(ISmax+1)%Iright = Ipvn-1
          Stack(ISmax+2)%Ileft = Ipvn + 1
          Stack(ISmax+2)%Iright = Iright
          ISpos = ISpos + 1
          ISmax = ISmax + 2
       End If
    End Do

    Deallocate(Stack)

    Return

  CONTAINS

    ! ***********************************

    Subroutine InsrtLC_index_i4(XX, IIpt, IIl, IIr)

      ! ***********************************
      ! * Perform an insertion sort of the list
      ! * XX(:) between index values IIl and IIr
      ! * but returns only IIpt(:) the permutations
      ! * made to sort.
      ! ***********************************

      IMPLICIT NONE

      integer(i4), Intent(in) :: XX(:)
      Integer(i4), Intent(inout) :: IIpt(:)
      Integer(i4), Intent(in) :: IIl, IIr

      integer(i4) :: RRtmp
      Integer(i4) :: II, JJ


      Do II = IIl+1, IIr
         RRtmp = XX(IIpt(II))
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(IIpt(JJ))) Then
               CALL swap(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
      End Do

      Return

    End Subroutine InsrtLC_index_i4

  end function Qsort_index_i4

  ! ***********************************

  Integer(i4) Function Partition_index_i4(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)

    ! ***********************************
    ! * This routine arranges the indexes Ipt of the array X
    ! * between the index values Ileft and Iright
    ! * positioning elements smallers than
    ! * X(Ipv) at the left and the others
    ! * at the right.
    ! * Internal routine used by Qsort.
    ! ***********************************

    IMPLICIT NONE

    integer(i4), Intent (in) :: X(:)
    Integer(i4), Intent (in) :: Ileft, Iright, Ipv
    Integer(i4), Intent (inout) :: Ipt(:)

    integer(i4) :: Rpv
    Integer(i4) :: I

    Rpv = X(Ipt(Ipv))
    CALL swap(Ipt, Ipv, Iright)
    Ipvfn = Ileft

    Do I = Ileft, Iright-1
       If (X(Ipt(I)) <= Rpv) Then
          CALL swap(Ipt, I, Ipvfn)
          Ipvfn = Ipvfn + 1
       End If
    End Do
    CALL swap(Ipt, Ipvfn, Iright)

    Return

  End Function Partition_index_i4



  ! ------------------------------------------------------------------
  ! http://www.fortran.com/qsort_c.f95

  ! Recursive Fortran 95 quicksort routine
  ! sorts real numbers into ascending numerical order
  ! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
  ! Based on algorithm from Cormen et al., Introduction to Algorithms, 1997 printing
  ! Made F conformant by Walt Brainerd

  recursive subroutine QsortC(A)

    IMPLICIT NONE

    real(dp), intent(in out), dimension(:) :: A
    integer(i4) :: iq

    if(size(A) > 1) then
       call Partitionc(A, iq)
       call QsortC(A(:iq-1))
       call QsortC(A(iq:))
    endif

  end subroutine QsortC

  subroutine Partitionc(A, marker)

    IMPLICIT NONE

    real(dp), intent(in out), dimension(:) :: A
    integer(i4), intent(out) :: marker
    integer(i4) :: i, j
    real(dp) :: temp
    real(dp) :: x      ! pivot point
    x = A(1)
    i= 0
    j= size(A) + 1

    do
       j = j-1
       do
          if (A(j) <= x) exit
          j = j-1
       end do
       i = i+1
       do
          if (A(i) >= x) exit
          i = i+1
       end do
       if (i < j) then
          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp
       elseif (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do

  end subroutine Partitionc




  ! ------------------------------------------------------------------
  ! http://jblevins.org/mirror/amiller/qsort.f90

  RECURSIVE SUBROUTINE quick_sort(list, order)

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.

    IMPLICIT NONE

    REAL(dp), DIMENSION (:), INTENT(INOUT)  :: list
    INTEGER(i4), DIMENSION (:), INTENT(OUT)  :: order

    ! Local variable
    INTEGER(i4) :: i

    DO i = 1, SIZE(list)
       order(i) = i
    END DO

    CALL quick_sort_1(1, SIZE(list))

  CONTAINS

    RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

      INTEGER(i4), INTENT(IN) :: left_end, right_end

      !     Local variables
      INTEGER(i4)             :: i, j, itemp
      REAL(dp)                :: reference, temp
      INTEGER(i4), PARAMETER  :: max_simple_sort_size = 6

      IF (right_end < left_end + max_simple_sort_size) THEN
         ! Use interchange sort for small lists
         CALL interchange_sort(left_end, right_end)

      ELSE
         ! Use partition ("quick") sort
         reference = list((left_end + right_end)/2)
         i = left_end - 1
         j = right_end + 1

         DO
            ! Scan list from left end until element >= reference is found
            DO
               i = i + 1
               IF (list(i) >= reference) EXIT
            END DO
            ! Scan list from right end until element <= reference is found
            DO
               j = j - 1
               IF (list(j) <= reference) EXIT
            END DO


            IF (i < j) THEN
               ! Swap two out-of-order elements
               temp = list(i)
               list(i) = list(j)
               list(j) = temp
               itemp = order(i)
               order(i) = order(j)
               order(j) = itemp
            ELSE IF (i == j) THEN
               i = i + 1
               EXIT
            ELSE
               EXIT
            END IF
         END DO

         IF (left_end < j) CALL quick_sort_1(left_end, j)
         IF (i < right_end) CALL quick_sort_1(i, right_end)
      END IF

    END SUBROUTINE quick_sort_1


    SUBROUTINE interchange_sort(left_end, right_end)

      INTEGER(i4), INTENT(IN) :: left_end, right_end

      !     Local variables
      INTEGER(i4)             :: i, j, itemp
      REAL(dp)                :: temp

      DO i = left_end, right_end - 1
         DO j = i+1, right_end
            IF (list(i) > list(j)) THEN
               temp = list(i)
               list(i) = list(j)
               list(j) = temp
               itemp = order(i)
               order(i) = order(j)
               order(j) = itemp
            END IF
         END DO
      END DO

    END SUBROUTINE interchange_sort

  END SUBROUTINE quick_sort




  ! ------------------------------------------------------------------
  ! http://balfortran.org

  recursive subroutine qsortb(a, limit)

    ! DUMMY ARGUMENTS
    real(dp),    dimension(:), intent(inout) :: A
    integer(i4), optional,     intent(in)    :: limit

    ! LOCAL VARIABLES
    integer(i4) :: nA
    integer(i4) :: ilimit
    integer(i4) :: left, right
    integer(i4) :: eighth
    real(dp)    :: pivot, temp
    integer(i4) :: marker
    real(dp), dimension(9) :: sample

    nA = size(A)
    if (present(limit)) then
       ilimit = limit
    else
       ilimit = 50_i4 ! http://balfortran.org/optimized_quicksort.html
    endif

    if (nA > 1) then
       if (nA > ilimit) then ! Do quicksort for large groups

          ! 9-SAMPLE PIVOT METHOD
          ! eighth = nA/8_i4
          ! sample = a(1:nA:eighth)
          if (mod(nA,8_i4) == 0) then
             eighth = nA/8_i4 - 1
          else
             eighth = nA/8_i4
          endif
          sample = a(1:8_i4*eighth+1:eighth)

          ! Sort Network for N=9, using Batcher's Merge-Exchange.
          ! Skip some steps because we only care about the median (5)
          if (sample(1) > sample(9)) then
             temp = sample(1)
             sample(1) = sample(9)
             sample(9) = temp
          end if
          if (sample(1) > sample(5)) then
             temp = sample(1)
             sample(1) = sample(5)
             sample(5) = temp
          end if
          if (sample(2) > sample(6)) then
             temp = sample(2)
             sample(2) = sample(6)
             sample(6) = temp
          end if
          if (sample(3) > sample(7)) then
             temp = sample(3)
             sample(3) = sample(7)
             sample(7) = temp
          end if
          if (sample(4) > sample(8)) then
             temp = sample(4)
             sample(4) = sample(8)
             sample(8) = temp
          end if
          if (sample(5) > sample(9)) then
             temp = sample(5)
             sample(5) = sample(9)
             sample(9) = temp
          end if
          if (sample(1) > sample(3)) then
             temp = sample(1)
             sample(1) = sample(3)
             sample(3) = temp
          end if
          if (sample(2) > sample(4)) then
             temp = sample(2)
             sample(2) = sample(4)
             sample(4) = temp
          end if
          if (sample(5) > sample(7)) then
             temp = sample(5)
             sample(5) = sample(7)
             sample(7) = temp
          end if
          if (sample(6) > sample(8)) then
             temp = sample(6)
             sample(6) = sample(8)
             sample(8) = temp
          end if
          if (sample(3) > sample(9)) then
             temp = sample(3)
             sample(3) = sample(9)
             sample(9) = temp
          end if
          if (sample(3) > sample(5)) then
             temp = sample(3)
             sample(3) = sample(5)
             sample(5) = temp
          end if
          if (sample(4) > sample(6)) then
             temp = sample(4)
             sample(4) = sample(6)
             sample(6) = temp
          end if
          if (sample(7) > sample(9)) then
             temp = sample(7)
             sample(7) = sample(9)
             sample(9) = temp
          end if
          if (sample(1) > sample(2)) then
             temp = sample(1)
             sample(1) = sample(2)
             sample(2) = temp
          end if
          if (sample(3) > sample(4)) then
             temp = sample(3)
             sample(3) = sample(4)
             sample(4) = temp
          end if
          if (sample(5) > sample(6)) then
             temp = sample(5)
             sample(5) = sample(6)
             sample(6) = temp
          end if
          if (sample(7) > sample(8)) then
             temp = sample(7)
             sample(7) = sample(8)
             sample(8) = temp
          end if
          if (sample(2) > sample(9)) then
             temp = sample(2)
             sample(2) = sample(9)
             sample(9) = temp
          end if
          if (sample(2) > sample(5)) then
             temp = sample(2)
             sample(2) = sample(5)
             sample(5) = temp
          end if
          if (sample(4) > sample(7)) then
             temp = sample(4)
             sample(4) = sample(7)
             sample(7) = temp
          end if
          !if (sample(6) > sample(9)) then
          ! temp = sample(6)
          ! sample(6) = sample(9)
          ! sample(9) = temp
          ! end if ! skipped
          !if (sample(2) > sample(3)) then
          ! temp = sample(2)
          ! sample(2) = sample(3)
          ! sample(3) = temp
          ! end if ! skipped
          if (sample(4) > sample(5)) then
             temp = sample(4)
             sample(4) = sample(5)
             sample(5) = temp
          end if
          !if (sample(6) > sample(7)) then
          ! temp = sample(6)
          ! sample(6) = sample(7)
          ! sample(7) = temp
          ! end if ! skipped
          !if (sample(8) > sample(9)) then
          ! temp = sample(8)
          ! sample(8) = sample(9)
          ! sample(9) = temp
          ! end if ! skipped
          pivot = sample(5)
          ! ************************
          left = 0
          right = nA + 1
          do while (left < right)
             right = right - 1
             do while (A(right) > pivot)
                right = right - 1
             end do
             left = left + 1
             do while (A(left) < pivot)
                left = left + 1
             end do
             if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
             end if
          end do

          if (left == right) then
             marker = left + 1
          else
             marker = left
          end if

          call qsortb(A(:marker-1), ilimit)
          call qsortb(A(marker:),   ilimit)

       else

          call InsertionSort(nA,A)    ! Insertion sort for small groups is faster than Quicksort

       end if

    end if

  end subroutine qsortb

  subroutine InsertionSort(nA,A)

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA
    real(dp), dimension(nA), intent(inout) :: A

    ! LOCAL VARIABLES
    real(dp) :: temp
    integer(i4) :: i, j

    outer: do i = 2, nA
       j = i - 1
       temp = A(i)
       inner: do
          if (j == 0) exit inner
          if (a(j) <= temp) exit inner
          A(j+1) = A(j)
          j = j - 1
       end do inner
       a(j+1) = temp
    end do outer

  end subroutine InsertionSort


  ! ------------------------------------------------------------------
  ! Optimised OpenMP quicksort http://balfortran.org

  ! Better OpenMP multi-threaded sort.
  ! Sorts real arrays
  ! Sorts equal parts of the input array on N number of threads.
  ! Uses both quicksort & insertion sort.  In quicksort, when the array size is less than a defined limit insertion sort is used
  ! instead, as it is faster for small arrays.
  ! After the initial sort, it does a bottom-up merge sort to combine the sorted chunks into a single sorted list.  The code will
  ! merge multiple chunks are once, and uses multiple threads for each merge.  For example, when run with 4 threads all 4 will be
  ! used to sort 4 chunks.  Then each set of 2 chunks will use 2 threads to merge together to form  2 larger chunks.  Finally,
  ! 4 threads will merge the 2 chunks together.
  ! Optimally uses power of 2 number of threads (2, 4, 8, etc).  The initial sort can use non-power of 2 threads, but not the merge.
  ! Uses double the memory of input array A as temp storage.

  ! Note: Subroutines use explicit size arrays instead of allocatable arrays.  The code runs faster this way.  Segmentation fault
  ! may occur for large arrays (i.e. 1 billion+ elements).  If a problem occurs either change to allocatable arrays or set memory
  ! stack size to unlimited.  In Linux the BASH command is "ulimit -s unlimited"

  ! Main sorting subroutine
  subroutine QSortMP_dp(A, ii, nt)

    ! USED MODULES
    !use, intrinsic :: iso_fortran_env, only: int64, real32, real64
    USE mo_kind, ONLY: i4, dp
    !$ use omp_lib

    ! DUMMY ARGUMENTS
    ! a: input data to be sorted
    ! ii: optional output of ranks (indexes to sort array)
    ! nt: number of threads, if omitted then uses all available
    real(dp),    dimension(:),           intent(inout) :: A
    integer(i4), dimension(:), optional, intent(out)   :: ii
    integer(i4),               optional, intent(in)    :: nt

    ! LOCAL VARIABLES
    ! A2: temp storage for partially sorted copies of A.  Program toggles between using arrays A and A2 as input and output to
    !   progressively sort the data.  The final sorted data is in array A.
    ! nt1: number of threads available.  Either a copy of input variable nt or omp_get_max_threads() if nt is not present.
    ! nt2: largest power of 2 number of threads (i.e. 2,4,8,16...)
    ! nt3: number of threads for each merge job.
    ! chunks: number of sort chunks for multithreaded quicksort.  Uses chunks_per_thread for each thread.  Using more smaller chunks
    !      helps load balancing.
    ! swap_flag: flag used to toggle between array a and array a2
    ! t: thread number
    ! s, f: Chunk indexes (start and finish) for each thread used by the initial quicksort/insertion sort combo.
    ! i, j: loop counters to merge sorted chunks together in order For example, with 4 threads first merge chunks 1 & 2 and 3 & 4 at the
    !       same time.  Next merge the previously made 1+2 chunk and 3+4 chunk to generate the completely sorted result.
    ! levels: Log(base 2) of nt2 (i.e 8 threads has 3 levels 2**3).  Extent of outer merge loop (level_loop).
    ! step: number of threads available for each merge.
    ! gap: number of pieces between left and right chunks to be merged.
    ! span: number of pieces in each left and right chunk.  As level_loop progresses each left and right chunk spans multiple chunks.
    ! l1, l2, r1, r2: index locations of left and right chunks to be merged.  1 is start, 2 is finish.
    ! left_part: temp array half the size of chunks being merged
    ! insert_limit: array size limit when quicksort changes to insertion sort.  50 is a good value for small data types.  As derived
    !       derived type increases in bytes the insert_limit should become lower.
    ! chunks_per_thread: number of work chunks for each thread.  More chunks improve load balancing.  However, too many can be slow.
    !       Quicksort is more efficient than then merge.  For quicksort 1 is optimal. For merge sort 64 is optimal.  Must be a power of
    !       2 (1, 2, 4, 8...)
    ! s2, f2: start and finish of sections of array a2 to be copied into array a
    ! verbose: T means to output messages, F means to output no messages
    integer(i4) :: nA
    real(dp), dimension(:), allocatable :: A2
    integer(i4) :: nt1
    integer(i4) :: nt2
    integer(i4) :: nt3
    integer(i4) :: chunks
    logical :: swap_flag
    integer(i4), dimension(:), allocatable :: s, f
    integer(i4) :: i, j
    integer(i4) :: levels
    integer(i4) :: step
    integer(i4) :: gap
    integer(i4) :: span
    integer(i4) :: l1, l2, r1, r2
    integer(i4), parameter :: insert_limit = 50
    integer(i4) :: chunks_per_thread = 1
    logical, parameter :: verbose = .false.
    integer(i4) :: indx(size(A))
    integer(i4), dimension(:), allocatable :: indx2

    nA = size(A)
    ! FIND AVAILABLE THREADS
    nt1 = 1  ! default to 1 thread in not using openmp
    !$ nt1 = omp_get_max_threads()  ! get max threads available if using OpenMP
    ! if (nt1 == 1) then
    !    write (*,*) "WARNING: Multi-threaded sort requested, but either system has only one CPU core or OpenMP is not enabled."
    ! end if
    if (present(nt)) nt1 = nt

    forall(i=1:nA) indx(i) = i

    multithread: if (nA < 10 .or. nt1 == 1) then

       ! Single-threaded
       if (verbose) write (*,*) "Single threaded"
       call QuickSortA_dp(A, nA, indx, insert_limit)

    else multithread ! PARALLEL SORT

       ! THREAD AND CHUNK CALCULATIONS
       nt2 = 2 ** int(log(real(nt1))/log(2.0)) ! get largest power of 2 number of threads (i.e. 2,4,8,16...)
       chunks = nt2 * chunks_per_thread           ! set number of smaller work chunks to increase load balancing
       levels = nint(log(real(chunks))/log(2.0))  ! set number of levels of merging
       swap_flag = .true.
       ! Make adjustments to ensure results end in array A for 2, 8, 32... threads and avoid an extra data copy.
       if (mod(levels,2) /= 0) then
          if (chunks_per_thread == 1) then    ! Ensure results end in array A by doubling chunks and possibly doubling threads
             ! If no thread count specified then double threads if more threads available.  This may result in more than one
             ! thread per core, but execution is faster.
             if (.not. present(nt) .and. nt1 > nt2) then
                nt2 = nt2 * 2
                nt1 = nt2
             end if
             chunks_per_thread = chunks_per_thread * 2
             chunks = chunks * 2
             levels = levels + 1
          else
             chunks_per_thread = chunks_per_thread / 2
             chunks = chunks / 2
             levels = levels - 1
          end if
       end if
       if (verbose) then
          write (*,"(A,I3)") "Threads used =", nt1
          write (*,"(A,I5)") "Chunks =", chunks
          if (nt2 /= nt1) write (*,"(A,I3,a)") "Only efficiently using",nt2," threads."
       end if

       allocate (s(chunks),f(chunks))        ! allocate start and finish locations of chunks
       allocate (A2(nA))                     ! allocate copy of A for temporary storage in merging
       allocate (indx2(nA))                  ! allocate copy of indx for temporary storage in merging

       ! SORT PIECES
       !$omp parallel do &
       !$omp num_threads(nt1) &
       ! !$omp private(t) &
       !$omp shared(A, A2, nA, indx, indx2, s, f, nt2, chunks, swap_flag)
       do i = 1, chunks
          s(i) = nA * (i-1) / chunks + 1     ! start
          f(i) = nA * i / chunks             ! finish
          call QuickSortA_dp(a(s(i):f(i)), &    ! section to be sorted
               &               f(i)-s(i)+1, &    ! size of section
               &               indx(s(i):f(i)), & ! indexes of section to be sorted
               &               insert_limit)     ! Insertion sort limit (50 is a good)
          if (.not. swap_flag) then
             a2(s(i):f(i)) = a(s(i):f(i))
             indx2(s(i):f(i)) = indx(s(i):f(i))
          endif
       end do
       !$omp end parallel do

       ! MERGE SORTED CHUNKS
       level_loop: do i = levels, 1, -1
          step = 2 ** (levels - i + 1)
          gap = 2 ** (levels - i)
          span = 2 ** (levels - i) - 1
          nt3 = max(2**(levels-i+1)/chunks_per_thread,1)                ! threads available for each merge
          !$omp parallel do &
          !$omp num_threads(nt2/nt3) &
          !            !$omp num_threads(nt2/nt3) &
          !$omp private (l1, l2, r1, r2) &
          !$omp shared (A, A2, indx, indx2, chunks, nt3, s, f, gap, span, step, swap_flag, nA, i)
          merge_loop: do j = 1, chunks, step
             l1 = s(j)
             l2 = f(j+span)
             r1 = s(j+gap)
             r2 = f(j+gap+span)

             ! MERGE 2 CHUNKS
             ! Merges sorted chunks.  Toggles input and output between arrays A and A2.  The duplicate memory is needed for the
             ! multithreaded merge.  Provisions are made above to ensure final output is in array A.
             if (swap_flag) then
                call Merge8A_mt_dp(A(l1:l2), &                  ! left part
                     &          l2-l1+1, &                  ! size of left part
                     &          A(r1:r2), &                 ! right part
                     &          r2-r1+1, &                  ! size of right part
                     &          A2(l1:r2), &                ! output array
                     &          r2-l1+1, &                  ! size of output array
                     &          indx(l1:l2), &              ! left index part
                     &          indx(r1:r2), &              ! right index part
                     &          indx2(l1:r2), &             ! index output array
                     &          nt3)                        ! number of threads
             else
                call Merge8A_mt_dp(A2(l1:l2), &                 ! left part
                     &          l2-l1+1, &                  ! size of left part
                     &          A2(r1:r2), &                ! right part
                     &          r2-r1+1, &                  ! size of right part
                     &          A(l1:r2), &                 ! output array
                     &          r2-l1+1, &                  ! size of output array
                     &          indx2(l1:l2), &             ! left index part
                     &          indx2(r1:r2), &             ! right index part
                     &          indx(l1:r2), &              ! index output array
                     &          nt3)                        ! number of threads
             end if

          enddo merge_loop
          !$omp end parallel do
          swap_flag = .not. swap_flag   ! toggle swap_flag to swap between array A and A2 as input and output arrays.
       enddo level_loop

       deallocate(s)
       deallocate(A2)
       deallocate(indx2)

    end if multithread

    if (present(ii)) ii = indx

  end subroutine QSortMP_dp

  ! Ascending Quicksort/Insertion Sort Combo for derived type
  recursive subroutine QuickSortA_dp(A,nA,indx,limit)
    ! USED MODULES
    use mo_kind, only: i4, dp

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA
    real(dp), dimension(nA), intent(inout) :: A
    integer(i4), dimension(nA), intent(inout) :: indx
    integer(i4), intent(in) :: limit

    ! LOCAL VARIABLES
    integer(i4) :: left, right
    real(dp) :: pivot
    real(dp) :: temp
    real(dp) :: temp2
    integer(i4) :: marker
    integer(i4) :: eighth
    real(dp), dimension(9) :: sample
    integer(i4) :: itemp2

    if (nA > 1) then
       if (nA > limit) then ! Do quicksort for large groups
          ! ************************
          ! 9-SAMPLE PIVOT METHOD
          eighth = nA/8
          !            sample = a(1:nA:eighth)

          sample(1) = a(1)
          sample(2) = a(1 + eighth)
          sample(3) = a(1 + eighth*2)
          sample(4) = a(1 + eighth*3)
          sample(5) = a(1 + eighth*4)
          sample(6) = a(1 + eighth*5)
          sample(7) = a(1 + eighth*6)
          sample(8) = a(1 + eighth*7)
          sample(9) = a(nA)
          ! Sort Network for N=9, using Batcher's Merge-Exchange. Skip some steps because I only care about the median (5)
          if (sample(1) > sample(9)) then
             temp = sample(1)
             sample(1) = sample(9)
             sample(9) = temp
          end if
          if (sample(1) > sample(5)) then
             temp = sample(1)
             sample(1) = sample(5)
             sample(5) = temp
          end if
          if (sample(2) > sample(6)) then
             temp = sample(2)
             sample(2) = sample(6)
             sample(6) = temp
          end if
          if (sample(3) > sample(7)) then
             temp = sample(3)
             sample(3) = sample(7)
             sample(7) = temp
          end if
          if (sample(4) > sample(8)) then
             temp = sample(4)
             sample(4) = sample(8)
             sample(8) = temp
          end if
          if (sample(5) > sample(9)) then
             temp = sample(5)
             sample(5) = sample(9)
             sample(9) = temp
          end if
          if (sample(1) > sample(3)) then
             temp = sample(1)
             sample(1) = sample(3)
             sample(3) = temp
          end if
          if (sample(2) > sample(4)) then
             temp = sample(2)
             sample(2) = sample(4)
             sample(4) = temp
          end if
          if (sample(5) > sample(7)) then
             temp = sample(5)
             sample(5) = sample(7)
             sample(7) = temp
          end if
          if (sample(6) > sample(8)) then
             temp = sample(6)
             sample(6) = sample(8)
             sample(8) = temp
          end if
          if (sample(3) > sample(9)) then
             temp = sample(3)
             sample(3) = sample(9)
             sample(9) = temp
          end if
          if (sample(3) > sample(5)) then
             temp = sample(3)
             sample(3) = sample(5)
             sample(5) = temp
          end if
          if (sample(4) > sample(6)) then
             temp = sample(4)
             sample(4) = sample(6)
             sample(6) = temp
          end if
          if (sample(7) > sample(9)) then
             temp = sample(7)
             sample(7) = sample(9)
             sample(9) = temp
          end if
          if (sample(1) > sample(2)) then
             temp = sample(1)
             sample(1) = sample(2)
             sample(2) = temp
          end if
          if (sample(3) > sample(4)) then
             temp = sample(3)
             sample(3) = sample(4)
             sample(4) = temp
          end if
          if (sample(5) > sample(6)) then
             temp = sample(5)
             sample(5) = sample(6)
             sample(6) = temp
          end if
          if (sample(7) > sample(8)) then
             temp = sample(7)
             sample(7) = sample(8)
             sample(8) = temp
          end if
          if (sample(2) > sample(9)) then
             temp = sample(2)
             sample(2) = sample(9)
             sample(9) = temp
          end if
          if (sample(2) > sample(5)) then
             temp = sample(2)
             sample(2) = sample(5)
             sample(5) = temp
          end if
          if (sample(4) > sample(7)) then
             temp = sample(4)
             sample(4) = sample(7)
             sample(7) = temp
          end if
          !if (sample(6) > sample(9)) then
          ! temp = sample(6)
          ! sample(6) = sample(9)
          ! sample(9) = temp
          ! end if ! skipped
          !if (sample(2) > sample(3)) then
          ! temp = sample(2)
          ! sample(2) = sample(3)
          ! sample(3) = temp
          ! end if ! skipped
          if (sample(4) > sample(5)) then
             temp = sample(4)
             sample(4) = sample(5)
             sample(5) = temp
          end if
          !if (sample(6) > sample(7)) then
          ! temp = sample(6)
          ! sample(6) = sample(7)
          ! sample(7) = temp
          ! end if ! skipped
          !if (sample(8) > sample(9)) then
          ! temp = sample(8)
          ! sample(8) = sample(9)
          ! sample(9) = temp
          ! end if ! skipped
          pivot = sample(5)
          ! ************************
          left = 0
          right = nA + 1
          do while (left < right)
             right = right - 1
             do while (A(right) > pivot)
                right = right - 1
             end do
             left = left + 1
             do while (A(left) < pivot)
                left = left + 1
             end do
             if (left < right) then
                temp2 = A(left)
                A(left) = A(right)
                A(right) = temp2
                itemp2 = indx(left)
                indx(left) = indx(right)
                indx(right) = itemp2
             end if
          end do

          if (left == right) then
             marker = left + 1
          else
             marker = left
          end if

          call QuickSortA_dp(A(:marker-1),marker-1,indx(:marker-1),limit)
          call QuickSortA_dp(A(marker:),nA-marker+1,indx(marker:),limit)

       else
          call InsertionSortA_dp(A,nA,indx)    ! Insertion sort for small groups is faster than Quicksort
       end if
    end if

  end subroutine QuickSortA_dp

  subroutine InsertionSortA_dp(A,nA,indx)
    ! USED MODULES
    use mo_kind, only: i4, dp

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA
    real(dp), dimension(nA), intent(inout) :: A
    integer(i4), dimension(nA), intent(inout) :: indx

    ! LOCAL VARIABLES
    real(dp) :: temp
    integer(i4) :: i, j
    integer(i4) :: itemp

    outter: do i = 2, nA
       j = i - 1
       temp = A(i)
       itemp = indx(i)
       inner: do
          if (j == 0) exit inner
          if (a(j) <= temp) exit inner
          A(j+1) = A(j)
          indx(j+1) = indx(j)
          j = j - 1
       end do inner
       a(j+1) = temp
       indx(j+1) = itemp
    end do outter

  end subroutine InsertionSortA_dp

  ! MULTI-THREAD MERGE
  ! Divides array a into nt number of parts.  Lookup matching values in array b.  Merges the parts together into new memory array c.
  ! Ascending sort
  subroutine Merge8a_mt_dp(A,nA,B,nB,C,nC,iA,iB,iC,nt)

    ! USED MODULES
    use mo_kind,  only: i4, dp
    use mo_utils, only: eq
    !$ use omp_lib, only: omp_get_thread_num

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA, nB, nC                      ! Size of arrays.  Normal usage: nA+nB = nC
    real(dp), dimension(nA), intent(in) :: A     ! left part
    real(dp), dimension(nB), intent(in) :: B     ! right part
    real(dp), dimension(nC), intent(out) :: C    ! sorted output must be new memory (not overlapping a or b)
    integer(i4), dimension(nA), intent(in) :: iA     ! left index part
    integer(i4), dimension(nB), intent(in) :: iB     ! right index part
    integer(i4), dimension(nC), intent(out) :: iC    ! sorted output must be new memory (not overlapping a or b)
    integer(i4), intent(in) :: nt                                           ! number of threads

    ! LOCAL VARIABLES
    integer(i4) :: t                                             ! thread number
    integer(i4), dimension(0:nt) :: dividerA        ! section divider of array a
    integer(i4), dimension(0:nt) :: dividerB        ! section divider of array b
    integer(i4), dimension(0:nt) :: dividerC        ! section divider of array c
    real(dp) :: divider_value                      ! value of section dividers
    integer(i4) :: nA_part, nB_part                 ! number of left and right elements of group
    integer(i4) :: nC_part                          ! number of sorted elements in group

    t=1
    !$ call omp_set_nested(.true.)

    ! initialize first and last divider
    dividerA(0) = 0
    dividerB(0) = 0
    dividerB(nt) = nB
    dividerC(0) = 0
    dividerC(nt) = nC

    !$omp parallel &
    !$omp num_threads(nt) &
    !$omp default(none) &
    !$omp private (t, divider_value, nA_part, nB_part, nC_part) &
    !$omp shared (a, nA, b, nB, c, nC, iA, iB, iC, nt, dividerA, dividerB, dividerC)

    !$ t = omp_get_thread_num() + 1  ! add 1 to make first thread 1 instead of 0
    dividerA(t) = nA * t / nt   ! starting value for dividerA.  Will change if not last in group.
    if (t < nt) then
       if (eq(a(dividerA(t)),a(dividerA(t)-1))) then
          ! Get new DividerA because starting value is in group, and not the last in the group.
          dividerA(t) = LookupAscending_dp(a,nA,a(dividerA(t))) - 1
       end if
    end if
    ! POSSIBILITIES
    ! dividerA is zero: Previous groups took all A values (same number extends from one dividerA to the end of A)
    ! dividerA if equal to nA: All values after dividerA starting value are same for the rest of A.
    ! dividerA is between 0 and nA: There are values that exceed the dividerA starting value (normal operation)

    divider_value = a(dividerA(t))
    if (t < nt) then
       ! find closest match that is just less than divider_value in array b
       dividerB(t) = LookupAscending_dp(b,nB,divider_value) - 1
       dividerC(t) = dividerA(t) + dividerB(t)
    end if
    !$omp barrier   ! make sure all divider values are calculated because each thread needs previous thread's divider values

    nA_part = dividerA(t) - (dividerA(t-1)+1) + 1
    nB_part = dividerB(t) - (dividerB(t-1)+1) + 1
    nC_part = dividerC(t) - (dividerC(t-1)+1) + 1

    ! POSSIBILITIES
    ! 1. nA_part == 0: Caused by dividerA being equal to nA.  This happens when all values in A beyond starting value of dividerA are
    !    equal.  After first thread, copy none of A-part and all of B-part into C.  First thead does normal (3) operation with all of A
    !    and a part of B (nA_part is not equal to 0 for the first thread, only for all threads greater than one)
    ! 2. nB_part == 0: Happens if Lookup finds no value in B that is less than divider_value in A.  Copy A-part into C.
    ! 3. Neither nA_part or nB_part are equal to 0.  Merge A-part and B-part (normal operation).

    if (nA_part == 0) then      ! possibility 1
       if (nB_part > 0) then
          C(dividerC(t-1)+1:dividerC(t-1)+nB_part) = B(dividerB(t-1)+1:dividerB(t)) ! copy only B-part
          iC(dividerC(t-1)+1:dividerC(t-1)+nB_part) = iB(dividerB(t-1)+1:dividerB(t)) ! copy only B-part
       endif
    else if (nB_part == 0) then ! possibility 2
       C(dividerC(t-1)+1:dividerC(t)) = A(dividerA(t-1)+1:dividerA(t)) ! copy only A-part
       iC(dividerC(t-1)+1:dividerC(t)) = iA(dividerA(t-1)+1:dividerA(t)) ! copy only A-part
    else                        ! possibility 3
       call Merge8a_dp(   A(dividerA(t-1)+1:dividerA(t)), nA_part, &   ! A-part
            &               B(dividerB(t-1)+1:dividerB(t)), nB_part, &   ! B-part
            &               C(dividerC(t-1)+1:dividerC(t)), nC_part, &  ! sorted part
            &               iA(dividerA(t-1)+1:dividerA(t)), &   ! A-part
            &               iB(dividerB(t-1)+1:dividerB(t)), &   ! B-part
            &               iC(dividerC(t-1)+1:dividerC(t)))     ! sorted part
    end if

    !$omp end parallel

  end subroutine Merge8a_mt_dp

  ! find first location where value is exceeded in array a.  Array a must be sorted in ascending order.  Uses binary search.
  function LookupAscending_dp(a,nA,value)

    ! USED MODULES
    use mo_kind, only: i4, dp

    ! DUMMY ARGUMENTS
    integer(i4) :: LookupAscending_dp
    integer(i4), intent(in) :: nA         ! Size of a
    real(dp), dimension(nA), intent(in) :: a
    real(dp), intent(in) :: value

    ! LOCAL VARIABLES
    integer(i4) :: half ! half the distance between binary searches

    if (a(1) > value) then
       ! Note: If using this function for other purposes, you may want to set LookupAscending to zero when the first value exceeds
       ! the value to make the result distinct from the first value being the match, but returning that value causes the
       ! multi-threaded merge to fail.
       LookupAscending_dp = 1    ! first value exceeds
    else
       if (a(nA) <= value) then
          LookupAscending_dp = nA + 1 ! last value is too small
       else    ! value is in between.  Find it.
          LookupAscending_dp = (nA+1)/2    ! starting location in middle
          half = LookupAscending_dp
          do
             half = (half+1) / 2
             if (a(LookupAscending_dp) > value) then
                if (a(LookupAscending_dp-1) <= value) then
                   exit    ! found
                else
                   LookupAscending_dp = max(LookupAscending_dp - half, 1) ! move down half
                end if
             else
                LookupAscending_dp = min(LookupAscending_dp + half, nA) ! move up half
             end if
          end do
       end if
    end if

  end function LookupAscending_dp

  ! Ascending merge (merges 2 ascending sorted lists into 1 ascending sorted list)
  subroutine Merge8a_dp(a,nA,b,nB,c,nC,ia,ib,ic)

    ! USED MODULES
    use mo_kind, only: i4, dp

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA, nB, nC         ! Size of arrays.  Normal usage: nC = nA+nB
    !   real(dp), dimension(:), intent(in) :: A ! Make sure A is a copy of any subset of C or sort overwrites as it goes
    !   real(dp), dimension(:), intent(in) :: B        ! B overlays C(nA+1:nC)
    !   real(dp), dimension(:), intent(out) :: C
    !    ! ***NOTE: USING EXPLICIT ARRAY SIZES AS SHOW BELOW SOMETIMES CAUSES SEGMENTATION FAULTS WITH LARGE ARRAY SIZES
    !   automatic arrays are faster than assumed size arrays.  Switch to above variable declarations if problems happen.
    real(dp), dimension(nA), intent(in) :: A  ! left part
    real(dp), dimension(nB), intent(in) :: B ! right part
    real(dp), dimension(nC), intent(out) :: C ! output array
    integer(i4), dimension(nA), intent(in) :: iA  ! left part
    integer(i4), dimension(nB), intent(in) :: iB ! right part
    integer(i4), dimension(nC), intent(out) :: iC ! output array
    ! Note: for single-threaded merge, such as the call from MergeSort, array a is a copy of the left part of array c.  Also, B overlays
    ! array c(nA+1:nC).  As multi-threaded merge array a and b are parts of the same input array.  Multi-threaded usage also requires
    ! array c to be different memory not overlapping arrays a or b.

    ! LOCAL VARIABLES
    integer(i4) :: i, j, k

    i = 1
    j = 1
    k = 1
    do
       if (i > nA .or. j > nB) exit
       if (A(i) <= B(j)) then
          C(k) = A(i)
          iC(k) = iA(i)
          i = i + 1
       else
          C(k) = B(j)
          iC(k) = iB(j)
          j = j + 1
       end if
       k = k + 1
    end do
    if (i <= nA) then
       C(k:) = A(i:)
       iC(k:) = iA(i:)
       return
    end if
    if (j <= nB) then
       C(k:) = B(j:)    ! This statement is only necessary for multi-threaded merge
       iC(k:) = iB(j:)    ! This statement is only necessary for multi-threaded merge
    endif

  end subroutine Merge8a_dp


  subroutine QSortMP_sp(A, ii, nt)

    ! USED MODULES
    !use, intrinsic :: iso_fortran_env, only: int64, real32, real64
    USE mo_kind, ONLY: i4, sp
    !$ use omp_lib

    ! DUMMY ARGUMENTS
    ! a: input data to be sorted
    ! ii: optional output of ranks (indexes to sort array)
    ! nt: number of threads, if omitted then uses all available
    real(sp),    dimension(:),           intent(inout) :: A
    integer(i4), dimension(:), optional, intent(out)   :: ii
    integer(i4),               optional, intent(in)    :: nt

    ! LOCAL VARIABLES
    ! A2: temp storage for partially sorted copies of A.  Program toggles between using arrays A and A2 as input and output to
    !   progressively sort the data.  The final sorted data is in array A.
    ! nt1: number of threads available.  Either a copy of input variable nt or omp_get_max_threads() if nt is not present.
    ! nt2: largest power of 2 number of threads (i.e. 2,4,8,16...)
    ! nt3: number of threads for each merge job.
    ! chunks: number of sort chunks for multithreaded quicksort.  Uses chunks_per_thread for each thread.  Using more smaller chunks
    !      helps load balancing.
    ! swap_flag: flag used to toggle between array a and array a2
    ! t: thread number
    ! s, f: Chunk indexes (start and finish) for each thread used by the initial quicksort/insertion sort combo.
    ! i, j: loop counters to merge sorted chunks together in order For example, with 4 threads first merge chunks 1 & 2 and 3 & 4 at the
    !       same time.  Next merge the previously made 1+2 chunk and 3+4 chunk to generate the completely sorted result.
    ! levels: Log(base 2) of nt2 (i.e 8 threads has 3 levels 2**3).  Extent of outer merge loop (level_loop).
    ! step: number of threads available for each merge.
    ! gap: number of pieces between left and right chunks to be merged.
    ! span: number of pieces in each left and right chunk.  As level_loop progresses each left and right chunk spans multiple chunks.
    ! l1, l2, r1, r2: index locations of left and right chunks to be merged.  1 is start, 2 is finish.
    ! left_part: temp array half the size of chunks being merged
    ! insert_limit: array size limit when quicksort changes to insertion sort.  50 is a good value for small data types.  As derived
    !       derived type increases in bytes the insert_limit should become lower.
    ! chunks_per_thread: number of work chunks for each thread.  More chunks improve load balancing.  However, too many can be slow.
    !       Quicksort is more efficient than then merge.  For quicksort 1 is optimal. For merge sort 64 is optimal.  Must be a power of
    !       2 (1, 2, 4, 8...)
    ! s2, f2: start and finish of sections of array a2 to be copied into array a
    ! verbose: T means to output messages, F means to output no messages
    integer(i4) :: nA
    real(sp), dimension(:), allocatable :: A2
    integer(i4) :: nt1
    integer(i4) :: nt2
    integer(i4) :: nt3
    integer(i4) :: chunks
    logical :: swap_flag
    integer(i4), dimension(:), allocatable :: s, f
    integer(i4) :: i, j
    integer(i4) :: levels
    integer(i4) :: step
    integer(i4) :: gap
    integer(i4) :: span
    integer(i4) :: l1, l2, r1, r2
    integer(i4), parameter :: insert_limit = 50
    integer(i4) :: chunks_per_thread = 1
    logical, parameter :: verbose = .false.
    integer(i4) :: indx(size(A))
    integer(i4), dimension(:), allocatable :: indx2

    nA = size(A)
    ! FIND AVAILABLE THREADS
    nt1 = 1  ! default to 1 thread in not using openmp
    !$ nt1 = omp_get_max_threads()  ! get max threads available if using OpenMP
    ! if (nt1 == 1) then
    !    write (*,*) "WARNING: Multi-threaded sort requested, but either system has only one CPU core or OpenMP is not enabled."
    ! end if
    if (present(nt)) nt1 = nt

    forall(i=1:nA) indx(i) = i

    multithread: if (nA < 10 .or. nt1 == 1) then

       ! Single-threaded
       if (verbose) write (*,*) "Single threaded"
       call QuickSortA_sp(A, nA, indx, insert_limit)

    else multithread ! PARALLEL SORT

       ! THREAD AND CHUNK CALCULATIONS
       nt2 = 2 ** int(log(real(nt1))/log(2.0)) ! get largest power of 2 number of threads (i.e. 2,4,8,16...)
       chunks = nt2 * chunks_per_thread           ! set number of smaller work chunks to increase load balancing
       levels = nint(log(real(chunks))/log(2.0))  ! set number of levels of merging
       swap_flag = .true.
       ! Make adjustments to ensure results end in array A for 2, 8, 32... threads and avoid an extra data copy.
       if (mod(levels,2) /= 0) then
          if (chunks_per_thread == 1) then    ! Ensure results end in array A by doubling chunks and possibly doubling threads
             ! If no thread count specified then double threads if more threads available.  This may result in more than one
             ! thread per core, but execution is faster.
             if (.not. present(nt) .and. nt1 > nt2) then
                nt2 = nt2 * 2
                nt1 = nt2
             end if
             chunks_per_thread = chunks_per_thread * 2
             chunks = chunks * 2
             levels = levels + 1
          else
             chunks_per_thread = chunks_per_thread / 2
             chunks = chunks / 2
             levels = levels - 1
          end if
       end if
       if (verbose) then
          write (*,"(A,I3)") "Threads used =", nt1
          write (*,"(A,I5)") "Chunks =", chunks
          if (nt2 /= nt1) write (*,"(A,I3,a)") "Only efficiently using",nt2," threads."
       end if

       allocate (s(chunks),f(chunks))        ! allocate start and finish locations of chunks
       allocate (A2(nA))                     ! allocate copy of A for temporary storage in merging
       allocate (indx2(nA))                  ! allocate copy of indx for temporary storage in merging

       ! SORT PIECES
       !$omp parallel do &
       !$omp num_threads(nt1) &
       ! !$omp private(t) &
       !$omp shared(A, A2, nA, indx, indx2, s, f, nt2, chunks, swap_flag)
       do i = 1, chunks
          s(i) = nA * (i-1) / chunks + 1     ! start
          f(i) = nA * i / chunks             ! finish
          call QuickSortA_sp(a(s(i):f(i)), &    ! section to be sorted
               &               f(i)-s(i)+1, &    ! size of section
               &               indx(s(i):f(i)), & ! indexes of section to be sorted
               &               insert_limit)     ! Insertion sort limit (50 is a good)
          if (.not. swap_flag) then
             a2(s(i):f(i)) = a(s(i):f(i))
             indx2(s(i):f(i)) = indx(s(i):f(i))
          endif
       end do
       !$omp end parallel do

       ! MERGE SORTED CHUNKS
       level_loop: do i = levels, 1, -1
          step = 2 ** (levels - i + 1)
          gap = 2 ** (levels - i)
          span = 2 ** (levels - i) - 1
          nt3 = max(2**(levels-i+1)/chunks_per_thread,1)                ! threads available for each merge
          !$omp parallel do &
          !$omp num_threads(nt2/nt3) &
          !            !$omp num_threads(nt2/nt3) &
          !$omp private (l1, l2, r1, r2) &
          !$omp shared (A, A2, indx, indx2, chunks, nt3, s, f, gap, span, step, swap_flag, nA, i)
          merge_loop: do j = 1, chunks, step
             l1 = s(j)
             l2 = f(j+span)
             r1 = s(j+gap)
             r2 = f(j+gap+span)

             ! MERGE 2 CHUNKS
             ! Merges sorted chunks.  Toggles input and output between arrays A and A2.  The duplicate memory is needed for the
             ! multithreaded merge.  Provisions are made above to ensure final output is in array A.
             if (swap_flag) then
                call Merge8A_mt_sp(A(l1:l2), &                  ! left part
                     &          l2-l1+1, &                  ! size of left part
                     &          A(r1:r2), &                 ! right part
                     &          r2-r1+1, &                  ! size of right part
                     &          A2(l1:r2), &                ! output array
                     &          r2-l1+1, &                  ! size of output array
                     &          indx(l1:l2), &              ! left index part
                     &          indx(r1:r2), &              ! right index part
                     &          indx2(l1:r2), &             ! index output array
                     &          nt3)                        ! number of threads
             else
                call Merge8A_mt_sp(A2(l1:l2), &                 ! left part
                     &          l2-l1+1, &                  ! size of left part
                     &          A2(r1:r2), &                ! right part
                     &          r2-r1+1, &                  ! size of right part
                     &          A(l1:r2), &                 ! output array
                     &          r2-l1+1, &                  ! size of output array
                     &          indx2(l1:l2), &             ! left index part
                     &          indx2(r1:r2), &             ! right index part
                     &          indx(l1:r2), &              ! index output array
                     &          nt3)                        ! number of threads
             end if

          enddo merge_loop
          !$omp end parallel do
          swap_flag = .not. swap_flag   ! toggle swap_flag to swap between array A and A2 as input and output arrays.
       enddo level_loop

       deallocate(s)
       deallocate(A2)
       deallocate(indx2)

    end if multithread

    if (present(ii)) ii = indx

  end subroutine QSortMP_sp

  ! Ascending Quicksort/Insertion Sort Combo for derived type
  recursive subroutine QuickSortA_sp(A,nA,indx,limit)
    ! USED MODULES
    use mo_kind, only: i4, sp

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA
    real(sp), dimension(nA), intent(inout) :: A
    integer(i4), dimension(nA), intent(inout) :: indx
    integer(i4), intent(in) :: limit

    ! LOCAL VARIABLES
    integer(i4) :: left, right
    real(sp) :: pivot
    real(sp) :: temp
    real(sp) :: temp2
    integer(i4) :: marker
    integer(i4) :: eighth
    real(sp), dimension(9) :: sample
    integer(i4) :: itemp2

    if (nA > 1) then
       if (nA > limit) then ! Do quicksort for large groups
          ! ************************
          ! 9-SAMPLE PIVOT METHOD
          eighth = nA/8
          !            sample = a(1:nA:eighth)

          sample(1) = a(1)
          sample(2) = a(1 + eighth)
          sample(3) = a(1 + eighth*2)
          sample(4) = a(1 + eighth*3)
          sample(5) = a(1 + eighth*4)
          sample(6) = a(1 + eighth*5)
          sample(7) = a(1 + eighth*6)
          sample(8) = a(1 + eighth*7)
          sample(9) = a(nA)
          ! Sort Network for N=9, using Batcher's Merge-Exchange. Skip some steps because I only care about the median (5)
          if (sample(1) > sample(9)) then
             temp = sample(1)
             sample(1) = sample(9)
             sample(9) = temp
          end if
          if (sample(1) > sample(5)) then
             temp = sample(1)
             sample(1) = sample(5)
             sample(5) = temp
          end if
          if (sample(2) > sample(6)) then
             temp = sample(2)
             sample(2) = sample(6)
             sample(6) = temp
          end if
          if (sample(3) > sample(7)) then
             temp = sample(3)
             sample(3) = sample(7)
             sample(7) = temp
          end if
          if (sample(4) > sample(8)) then
             temp = sample(4)
             sample(4) = sample(8)
             sample(8) = temp
          end if
          if (sample(5) > sample(9)) then
             temp = sample(5)
             sample(5) = sample(9)
             sample(9) = temp
          end if
          if (sample(1) > sample(3)) then
             temp = sample(1)
             sample(1) = sample(3)
             sample(3) = temp
          end if
          if (sample(2) > sample(4)) then
             temp = sample(2)
             sample(2) = sample(4)
             sample(4) = temp
          end if
          if (sample(5) > sample(7)) then
             temp = sample(5)
             sample(5) = sample(7)
             sample(7) = temp
          end if
          if (sample(6) > sample(8)) then
             temp = sample(6)
             sample(6) = sample(8)
             sample(8) = temp
          end if
          if (sample(3) > sample(9)) then
             temp = sample(3)
             sample(3) = sample(9)
             sample(9) = temp
          end if
          if (sample(3) > sample(5)) then
             temp = sample(3)
             sample(3) = sample(5)
             sample(5) = temp
          end if
          if (sample(4) > sample(6)) then
             temp = sample(4)
             sample(4) = sample(6)
             sample(6) = temp
          end if
          if (sample(7) > sample(9)) then
             temp = sample(7)
             sample(7) = sample(9)
             sample(9) = temp
          end if
          if (sample(1) > sample(2)) then
             temp = sample(1)
             sample(1) = sample(2)
             sample(2) = temp
          end if
          if (sample(3) > sample(4)) then
             temp = sample(3)
             sample(3) = sample(4)
             sample(4) = temp
          end if
          if (sample(5) > sample(6)) then
             temp = sample(5)
             sample(5) = sample(6)
             sample(6) = temp
          end if
          if (sample(7) > sample(8)) then
             temp = sample(7)
             sample(7) = sample(8)
             sample(8) = temp
          end if
          if (sample(2) > sample(9)) then
             temp = sample(2)
             sample(2) = sample(9)
             sample(9) = temp
          end if
          if (sample(2) > sample(5)) then
             temp = sample(2)
             sample(2) = sample(5)
             sample(5) = temp
          end if
          if (sample(4) > sample(7)) then
             temp = sample(4)
             sample(4) = sample(7)
             sample(7) = temp
          end if
          !if (sample(6) > sample(9)) then
          ! temp = sample(6)
          ! sample(6) = sample(9)
          ! sample(9) = temp
          ! end if ! skipped
          !if (sample(2) > sample(3)) then
          ! temp = sample(2)
          ! sample(2) = sample(3)
          ! sample(3) = temp
          ! end if ! skipped
          if (sample(4) > sample(5)) then
             temp = sample(4)
             sample(4) = sample(5)
             sample(5) = temp
          end if
          !if (sample(6) > sample(7)) then
          ! temp = sample(6)
          ! sample(6) = sample(7)
          ! sample(7) = temp
          ! end if ! skipped
          !if (sample(8) > sample(9)) then
          ! temp = sample(8)
          ! sample(8) = sample(9)
          ! sample(9) = temp
          ! end if ! skipped
          pivot = sample(5)
          ! ************************
          left = 0
          right = nA + 1
          do while (left < right)
             right = right - 1
             do while (A(right) > pivot)
                right = right - 1
             end do
             left = left + 1
             do while (A(left) < pivot)
                left = left + 1
             end do
             if (left < right) then
                temp2 = A(left)
                A(left) = A(right)
                A(right) = temp2
                itemp2 = indx(left)
                indx(left) = indx(right)
                indx(right) = itemp2
             end if
          end do

          if (left == right) then
             marker = left + 1
          else
             marker = left
          end if

          call QuickSortA_sp(A(:marker-1),marker-1,indx(:marker-1),limit)
          call QuickSortA_sp(A(marker:),nA-marker+1,indx(marker:),limit)

       else
          call InsertionSortA_sp(A,nA,indx)    ! Insertion sort for small groups is faster than Quicksort
       end if
    end if

  end subroutine QuickSortA_sp

  subroutine InsertionSortA_sp(A,nA,indx)
    ! USED MODULES
    use mo_kind, only: i4, sp

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA
    real(sp), dimension(nA), intent(inout) :: A
    integer(i4), dimension(nA), intent(inout) :: indx

    ! LOCAL VARIABLES
    real(sp) :: temp
    integer(i4) :: i, j
    integer(i4) :: itemp

    outter: do i = 2, nA
       j = i - 1
       temp = A(i)
       itemp = indx(i)
       inner: do
          if (j == 0) exit inner
          if (a(j) <= temp) exit inner
          A(j+1) = A(j)
          indx(j+1) = indx(j)
          j = j - 1
       end do inner
       a(j+1) = temp
       indx(j+1) = itemp
    end do outter

  end subroutine InsertionSortA_sp

  ! MULTI-THREAD MERGE
  ! Divides array a into nt number of parts.  Lookup matching values in array b.  Merges the parts together into new memory array c.
  ! Ascending sort
  subroutine Merge8a_mt_sp(A,nA,B,nB,C,nC,iA,iB,iC,nt)

    ! USED MODULES
    use mo_kind,  only: i4, sp
    use mo_utils, only: eq
    !$ use omp_lib, only: omp_get_thread_num

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA, nB, nC                      ! Size of arrays.  Normal usage: nA+nB = nC
    real(sp), dimension(nA), intent(in) :: A     ! left part
    real(sp), dimension(nB), intent(in) :: B     ! right part
    real(sp), dimension(nC), intent(out) :: C    ! sorted output must be new memory (not overlapping a or b)
    integer(i4), dimension(nA), intent(in) :: iA     ! left index part
    integer(i4), dimension(nB), intent(in) :: iB     ! right index part
    integer(i4), dimension(nC), intent(out) :: iC    ! sorted output must be new memory (not overlapping a or b)
    integer(i4), intent(in) :: nt                                           ! number of threads

    ! LOCAL VARIABLES
    integer(i4) :: t                                             ! thread number
    integer(i4), dimension(0:nt) :: dividerA        ! section divider of array a
    integer(i4), dimension(0:nt) :: dividerB        ! section divider of array b
    integer(i4), dimension(0:nt) :: dividerC        ! section divider of array c
    real(sp) :: divider_value                      ! value of section dividers
    integer(i4) :: nA_part, nB_part                 ! number of left and right elements of group
    integer(i4) :: nC_part                          ! number of sorted elements in group

    t=1
    !$ call omp_set_nested(.true.)

    ! initialize first and last divider
    dividerA(0) = 0
    dividerB(0) = 0
    dividerB(nt) = nB
    dividerC(0) = 0
    dividerC(nt) = nC

    !$omp parallel &
    !$omp num_threads(nt) &
    !$omp default(none) &
    !$omp private (t, divider_value, nA_part, nB_part, nC_part) &
    !$omp shared (a, nA, b, nB, c, nC, iA, iB, iC, nt, dividerA, dividerB, dividerC)

    !$ t = omp_get_thread_num() + 1  ! add 1 to make first thread 1 instead of 0
    dividerA(t) = nA * t / nt   ! starting value for dividerA.  Will change if not last in group.
    if (t < nt) then
       if (eq(a(dividerA(t)),a(dividerA(t)-1))) then
          ! Get new DividerA because starting value is in group, and not the last in the group.
          dividerA(t) = LookupAscending_sp(a,nA,a(dividerA(t))) - 1
       end if
    end if
    ! POSSIBILITIES
    ! dividerA is zero: Previous groups took all A values (same number extends from one dividerA to the end of A)
    ! dividerA if equal to nA: All values after dividerA starting value are same for the rest of A.
    ! dividerA is between 0 and nA: There are values that exceed the dividerA starting value (normal operation)

    divider_value = a(dividerA(t))
    if (t < nt) then
       ! find closest match that is just less than divider_value in array b
       dividerB(t) = LookupAscending_sp(b,nB,divider_value) - 1
       dividerC(t) = dividerA(t) + dividerB(t)
    end if
    !$omp barrier   ! make sure all divider values are calculated because each thread needs previous thread's divider values

    nA_part = dividerA(t) - (dividerA(t-1)+1) + 1
    nB_part = dividerB(t) - (dividerB(t-1)+1) + 1
    nC_part = dividerC(t) - (dividerC(t-1)+1) + 1

    ! POSSIBILITIES
    ! 1. nA_part == 0: Caused by dividerA being equal to nA.  This happens when all values in A beyond starting value of dividerA are
    !    equal.  After first thread, copy none of A-part and all of B-part into C.  First thead does normal (3) operation with all of A
    !    and a part of B (nA_part is not equal to 0 for the first thread, only for all threads greater than one)
    ! 2. nB_part == 0: Happens if Lookup finds no value in B that is less than divider_value in A.  Copy A-part into C.
    ! 3. Neither nA_part or nB_part are equal to 0.  Merge A-part and B-part (normal operation).

    if (nA_part == 0) then      ! possibility 1
       if (nB_part > 0) then
          C(dividerC(t-1)+1:dividerC(t-1)+nB_part) = B(dividerB(t-1)+1:dividerB(t)) ! copy only B-part
          iC(dividerC(t-1)+1:dividerC(t-1)+nB_part) = iB(dividerB(t-1)+1:dividerB(t)) ! copy only B-part
       endif
    else if (nB_part == 0) then ! possibility 2
       C(dividerC(t-1)+1:dividerC(t)) = A(dividerA(t-1)+1:dividerA(t)) ! copy only A-part
       iC(dividerC(t-1)+1:dividerC(t)) = iA(dividerA(t-1)+1:dividerA(t)) ! copy only A-part
    else                        ! possibility 3
       call Merge8a_sp(   A(dividerA(t-1)+1:dividerA(t)), nA_part, &   ! A-part
            &               B(dividerB(t-1)+1:dividerB(t)), nB_part, &   ! B-part
            &               C(dividerC(t-1)+1:dividerC(t)), nC_part, &  ! sorted part
            &               iA(dividerA(t-1)+1:dividerA(t)), &   ! A-part
            &               iB(dividerB(t-1)+1:dividerB(t)), &   ! B-part
            &               iC(dividerC(t-1)+1:dividerC(t)))     ! sorted part
    end if

    !$omp end parallel

  end subroutine Merge8a_mt_sp

  ! find first location where value is exceeded in array a.  Array a must be sorted in ascending order.  Uses binary search.
  function LookupAscending_sp(a,nA,value)

    ! USED MODULES
    use mo_kind, only: i4, sp

    ! DUMMY ARGUMENTS
    integer(i4) :: LookupAscending_sp
    integer(i4), intent(in) :: nA         ! Size of a
    real(sp), dimension(nA), intent(in) :: a
    real(sp), intent(in) :: value

    ! LOCAL VARIABLES
    integer(i4) :: half ! half the distance between binary searches

    if (a(1) > value) then
       ! Note: If using this function for other purposes, you may want to set LookupAscending to zero when the first value exceeds
       ! the value to make the result distinct from the first value being the match, but returning that value causes the
       ! multi-threaded merge to fail.
       LookupAscending_sp = 1    ! first value exceeds
    else
       if (a(nA) <= value) then
          LookupAscending_sp = nA + 1 ! last value is too small
       else    ! value is in between.  Find it.
          LookupAscending_sp = (nA+1)/2    ! starting location in middle
          half = LookupAscending_sp
          do
             half = (half+1) / 2
             if (a(LookupAscending_sp) > value) then
                if (a(LookupAscending_sp-1) <= value) then
                   exit    ! found
                else
                   LookupAscending_sp = max(LookupAscending_sp - half, 1) ! move down half
                end if
             else
                LookupAscending_sp = min(LookupAscending_sp + half, nA) ! move up half
             end if
          end do
       end if
    end if

  end function LookupAscending_sp

  ! Ascending merge (merges 2 ascending sorted lists into 1 ascending sorted list)
  subroutine Merge8a_sp(a,nA,b,nB,c,nC,ia,ib,ic)

    ! USED MODULES
    use mo_kind, only: i4, sp

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA, nB, nC         ! Size of arrays.  Normal usage: nC = nA+nB
    !   real(sp), dimension(:), intent(in) :: A ! Make sure A is a copy of any subset of C or sort overwrites as it goes
    !   real(sp), dimension(:), intent(in) :: B        ! B overlays C(nA+1:nC)
    !   real(sp), dimension(:), intent(out) :: C
    !    ! ***NOTE: USING EXPLICIT ARRAY SIZES AS SHOW BELOW SOMETIMES CAUSES SEGMENTATION FAULTS WITH LARGE ARRAY SIZES
    !   automatic arrays are faster than assumed size arrays.  Switch to above variable declarations if problems happen.
    real(sp), dimension(nA), intent(in) :: A  ! left part
    real(sp), dimension(nB), intent(in) :: B ! right part
    real(sp), dimension(nC), intent(out) :: C ! output array
    integer(i4), dimension(nA), intent(in) :: iA  ! left part
    integer(i4), dimension(nB), intent(in) :: iB ! right part
    integer(i4), dimension(nC), intent(out) :: iC ! output array
    ! Note: for single-threaded merge, such as the call from MergeSort, array a is a copy of the left part of array c.  Also, B overlays
    ! array c(nA+1:nC).  As multi-threaded merge array a and b are parts of the same input array.  Multi-threaded usage also requires
    ! array c to be different memory not overlapping arrays a or b.

    ! LOCAL VARIABLES
    integer(i4) :: i, j, k

    i = 1
    j = 1
    k = 1
    do
       if (i > nA .or. j > nB) exit
       if (A(i) <= B(j)) then
          C(k) = A(i)
          iC(k) = iA(i)
          i = i + 1
       else
          C(k) = B(j)
          iC(k) = iB(j)
          j = j + 1
       end if
       k = k + 1
    end do
    if (i <= nA) then
       C(k:) = A(i:)
       iC(k:) = iA(i:)
       return
    end if
    if (j <= nB) then
       C(k:) = B(j:)    ! This statement is only necessary for multi-threaded merge
       iC(k:) = iB(j:)    ! This statement is only necessary for multi-threaded merge
    endif

  end subroutine Merge8a_sp


  subroutine QSortMP_i4(A, ii, nt)

    ! USED MODULES
    !use, intrinsic :: iso_fortran_env, only: int64, real32, real64
    USE mo_kind, ONLY: i4, i4
    !$ use omp_lib

    ! DUMMY ARGUMENTS
    ! a: input data to be sorted
    ! ii: optional output of ranks (indexes to sort array)
    ! nt: number of threads, if omitted then uses all available
    integer(i4),    dimension(:),           intent(inout) :: A
    integer(i4), dimension(:), optional, intent(out)   :: ii
    integer(i4),               optional, intent(in)    :: nt

    ! LOCAL VARIABLES
    ! A2: temp storage for partially sorted copies of A.  Program toggles between using arrays A and A2 as input and output to
    !   progressively sort the data.  The final sorted data is in array A.
    ! nt1: number of threads available.  Either a copy of input variable nt or omp_get_max_threads() if nt is not present.
    ! nt2: largest power of 2 number of threads (i.e. 2,4,8,16...)
    ! nt3: number of threads for each merge job.
    ! chunks: number of sort chunks for multithreaded quicksort.  Uses chunks_per_thread for each thread.  Using more smaller chunks
    !      helps load balancing.
    ! swap_flag: flag used to toggle between array a and array a2
    ! t: thread number
    ! s, f: Chunk indexes (start and finish) for each thread used by the initial quicksort/insertion sort combo.
    ! i, j: loop counters to merge sorted chunks together in order For example, with 4 threads first merge chunks 1 & 2 and 3 & 4 at the
    !       same time.  Next merge the previously made 1+2 chunk and 3+4 chunk to generate the completely sorted result.
    ! levels: Log(base 2) of nt2 (i.e 8 threads has 3 levels 2**3).  Extent of outer merge loop (level_loop).
    ! step: number of threads available for each merge.
    ! gap: number of pieces between left and right chunks to be merged.
    ! span: number of pieces in each left and right chunk.  As level_loop progresses each left and right chunk spans multiple chunks.
    ! l1, l2, r1, r2: index locations of left and right chunks to be merged.  1 is start, 2 is finish.
    ! left_part: temp array half the size of chunks being merged
    ! insert_limit: array size limit when quicksort changes to insertion sort.  50 is a good value for small data types.  As derived
    !       derived type increases in bytes the insert_limit should become lower.
    ! chunks_per_thread: number of work chunks for each thread.  More chunks improve load balancing.  However, too many can be slow.
    !       Quicksort is more efficient than then merge.  For quicksort 1 is optimal. For merge sort 64 is optimal.  Must be a power of
    !       2 (1, 2, 4, 8...)
    ! s2, f2: start and finish of sections of array a2 to be copied into array a
    ! verbose: T means to output messages, F means to output no messages
    integer(i4) :: nA
    integer(i4), dimension(:), allocatable :: A2
    integer(i4) :: nt1
    integer(i4) :: nt2
    integer(i4) :: nt3
    integer(i4) :: chunks
    logical :: swap_flag
    integer(i4), dimension(:), allocatable :: s, f
    integer(i4) :: i, j
    integer(i4) :: levels
    integer(i4) :: step
    integer(i4) :: gap
    integer(i4) :: span
    integer(i4) :: l1, l2, r1, r2
    integer(i4), parameter :: insert_limit = 50
    integer(i4) :: chunks_per_thread = 1
    logical, parameter :: verbose = .false.
    integer(i4) :: indx(size(A))
    integer(i4), dimension(:), allocatable :: indx2

    nA = size(A)
    ! FIND AVAILABLE THREADS
    nt1 = 1  ! default to 1 thread in not using openmp
    !$ nt1 = omp_get_max_threads()  ! get max threads available if using OpenMP
    ! if (nt1 == 1) then
    !    write (*,*) "WARNING: Multi-threaded sort requested, but either system has only one CPU core or OpenMP is not enabled."
    ! end if
    if (present(nt)) nt1 = nt

    forall(i=1:nA) indx(i) = i

    multithread: if (nA < 10 .or. nt1 == 1) then

       ! Single-threaded
       if (verbose) write (*,*) "Single threaded"
       call QuickSortA_i4(A, nA, indx, insert_limit)

    else multithread ! PARALLEL SORT

       ! THREAD AND CHUNK CALCULATIONS
       nt2 = 2 ** int(log(real(nt1))/log(2.0)) ! get largest power of 2 number of threads (i.e. 2,4,8,16...)
       chunks = nt2 * chunks_per_thread           ! set number of smaller work chunks to increase load balancing
       levels = nint(log(real(chunks))/log(2.0))  ! set number of levels of merging
       swap_flag = .true.
       ! Make adjustments to ensure results end in array A for 2, 8, 32... threads and avoid an extra data copy.
       if (mod(levels,2) /= 0) then
          if (chunks_per_thread == 1) then    ! Ensure results end in array A by doubling chunks and possibly doubling threads
             ! If no thread count specified then double threads if more threads available.  This may result in more than one
             ! thread per core, but execution is faster.
             if (.not. present(nt) .and. nt1 > nt2) then
                nt2 = nt2 * 2
                nt1 = nt2
             end if
             chunks_per_thread = chunks_per_thread * 2
             chunks = chunks * 2
             levels = levels + 1
          else
             chunks_per_thread = chunks_per_thread / 2
             chunks = chunks / 2
             levels = levels - 1
          end if
       end if
       if (verbose) then
          write (*,"(A,I3)") "Threads used =", nt1
          write (*,"(A,I5)") "Chunks =", chunks
          if (nt2 /= nt1) write (*,"(A,I3,a)") "Only efficiently using",nt2," threads."
       end if

       allocate (s(chunks),f(chunks))        ! allocate start and finish locations of chunks
       allocate (A2(nA))                     ! allocate copy of A for temporary storage in merging
       allocate (indx2(nA))                  ! allocate copy of indx for temporary storage in merging

       ! SORT PIECES
       !$omp parallel do &
       !$omp num_threads(nt1) &
       ! !$omp private(t) &
       !$omp shared(A, A2, nA, indx, indx2, s, f, nt2, chunks, swap_flag)
       do i = 1, chunks
          s(i) = nA * (i-1) / chunks + 1     ! start
          f(i) = nA * i / chunks             ! finish
          call QuickSortA_i4(a(s(i):f(i)), &    ! section to be sorted
               &               f(i)-s(i)+1, &    ! size of section
               &               indx(s(i):f(i)), & ! indexes of section to be sorted
               &               insert_limit)     ! Insertion sort limit (50 is a good)
          if (.not. swap_flag) then
             a2(s(i):f(i)) = a(s(i):f(i))
             indx2(s(i):f(i)) = indx(s(i):f(i))
          endif
       end do
       !$omp end parallel do

       ! MERGE SORTED CHUNKS
       level_loop: do i = levels, 1, -1
          step = 2 ** (levels - i + 1)
          gap = 2 ** (levels - i)
          span = 2 ** (levels - i) - 1
          nt3 = max(2**(levels-i+1)/chunks_per_thread,1)                ! threads available for each merge
          !$omp parallel do &
          !$omp num_threads(nt2/nt3) &
          !            !$omp num_threads(nt2/nt3) &
          !$omp private (l1, l2, r1, r2) &
          !$omp shared (A, A2, indx, indx2, chunks, nt3, s, f, gap, span, step, swap_flag, nA, i)
          merge_loop: do j = 1, chunks, step
             l1 = s(j)
             l2 = f(j+span)
             r1 = s(j+gap)
             r2 = f(j+gap+span)

             ! MERGE 2 CHUNKS
             ! Merges sorted chunks.  Toggles input and output between arrays A and A2.  The duplicate memory is needed for the
             ! multithreaded merge.  Provisions are made above to ensure final output is in array A.
             if (swap_flag) then
                call Merge8A_mt_i4(A(l1:l2), &                  ! left part
                     &          l2-l1+1, &                  ! size of left part
                     &          A(r1:r2), &                 ! right part
                     &          r2-r1+1, &                  ! size of right part
                     &          A2(l1:r2), &                ! output array
                     &          r2-l1+1, &                  ! size of output array
                     &          indx(l1:l2), &              ! left index part
                     &          indx(r1:r2), &              ! right index part
                     &          indx2(l1:r2), &             ! index output array
                     &          nt3)                        ! number of threads
             else
                call Merge8A_mt_i4(A2(l1:l2), &                 ! left part
                     &          l2-l1+1, &                  ! size of left part
                     &          A2(r1:r2), &                ! right part
                     &          r2-r1+1, &                  ! size of right part
                     &          A(l1:r2), &                 ! output array
                     &          r2-l1+1, &                  ! size of output array
                     &          indx2(l1:l2), &             ! left index part
                     &          indx2(r1:r2), &             ! right index part
                     &          indx(l1:r2), &              ! index output array
                     &          nt3)                        ! number of threads
             end if

          enddo merge_loop
          !$omp end parallel do
          swap_flag = .not. swap_flag   ! toggle swap_flag to swap between array A and A2 as input and output arrays.
       enddo level_loop

       deallocate(s)
       deallocate(A2)
       deallocate(indx2)

    end if multithread

    if (present(ii)) ii = indx

  end subroutine QSortMP_i4

  ! Ascending Quicksort/Insertion Sort Combo for derived type
  recursive subroutine QuickSortA_i4(A,nA,indx,limit)
    ! USED MODULES
    use mo_kind, only: i4, i4

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA
    integer(i4), dimension(nA), intent(inout) :: A
    integer(i4), dimension(nA), intent(inout) :: indx
    integer(i4), intent(in) :: limit

    ! LOCAL VARIABLES
    integer(i4) :: left, right
    integer(i4) :: pivot
    integer(i4) :: temp
    integer(i4) :: temp2
    integer(i4) :: marker
    integer(i4) :: eighth
    integer(i4), dimension(9) :: sample
    integer(i4) :: itemp2

    if (nA > 1) then
       if (nA > limit) then ! Do quicksort for large groups
          ! ************************
          ! 9-SAMPLE PIVOT METHOD
          eighth = nA/8
          !            sample = a(1:nA:eighth)

          sample(1) = a(1)
          sample(2) = a(1 + eighth)
          sample(3) = a(1 + eighth*2)
          sample(4) = a(1 + eighth*3)
          sample(5) = a(1 + eighth*4)
          sample(6) = a(1 + eighth*5)
          sample(7) = a(1 + eighth*6)
          sample(8) = a(1 + eighth*7)
          sample(9) = a(nA)
          ! Sort Network for N=9, using Batcher's Merge-Exchange. Skip some steps because I only care about the median (5)
          if (sample(1) > sample(9)) then
             temp = sample(1)
             sample(1) = sample(9)
             sample(9) = temp
          end if
          if (sample(1) > sample(5)) then
             temp = sample(1)
             sample(1) = sample(5)
             sample(5) = temp
          end if
          if (sample(2) > sample(6)) then
             temp = sample(2)
             sample(2) = sample(6)
             sample(6) = temp
          end if
          if (sample(3) > sample(7)) then
             temp = sample(3)
             sample(3) = sample(7)
             sample(7) = temp
          end if
          if (sample(4) > sample(8)) then
             temp = sample(4)
             sample(4) = sample(8)
             sample(8) = temp
          end if
          if (sample(5) > sample(9)) then
             temp = sample(5)
             sample(5) = sample(9)
             sample(9) = temp
          end if
          if (sample(1) > sample(3)) then
             temp = sample(1)
             sample(1) = sample(3)
             sample(3) = temp
          end if
          if (sample(2) > sample(4)) then
             temp = sample(2)
             sample(2) = sample(4)
             sample(4) = temp
          end if
          if (sample(5) > sample(7)) then
             temp = sample(5)
             sample(5) = sample(7)
             sample(7) = temp
          end if
          if (sample(6) > sample(8)) then
             temp = sample(6)
             sample(6) = sample(8)
             sample(8) = temp
          end if
          if (sample(3) > sample(9)) then
             temp = sample(3)
             sample(3) = sample(9)
             sample(9) = temp
          end if
          if (sample(3) > sample(5)) then
             temp = sample(3)
             sample(3) = sample(5)
             sample(5) = temp
          end if
          if (sample(4) > sample(6)) then
             temp = sample(4)
             sample(4) = sample(6)
             sample(6) = temp
          end if
          if (sample(7) > sample(9)) then
             temp = sample(7)
             sample(7) = sample(9)
             sample(9) = temp
          end if
          if (sample(1) > sample(2)) then
             temp = sample(1)
             sample(1) = sample(2)
             sample(2) = temp
          end if
          if (sample(3) > sample(4)) then
             temp = sample(3)
             sample(3) = sample(4)
             sample(4) = temp
          end if
          if (sample(5) > sample(6)) then
             temp = sample(5)
             sample(5) = sample(6)
             sample(6) = temp
          end if
          if (sample(7) > sample(8)) then
             temp = sample(7)
             sample(7) = sample(8)
             sample(8) = temp
          end if
          if (sample(2) > sample(9)) then
             temp = sample(2)
             sample(2) = sample(9)
             sample(9) = temp
          end if
          if (sample(2) > sample(5)) then
             temp = sample(2)
             sample(2) = sample(5)
             sample(5) = temp
          end if
          if (sample(4) > sample(7)) then
             temp = sample(4)
             sample(4) = sample(7)
             sample(7) = temp
          end if
          !if (sample(6) > sample(9)) then
          ! temp = sample(6)
          ! sample(6) = sample(9)
          ! sample(9) = temp
          ! end if ! skipped
          !if (sample(2) > sample(3)) then
          ! temp = sample(2)
          ! sample(2) = sample(3)
          ! sample(3) = temp
          ! end if ! skipped
          if (sample(4) > sample(5)) then
             temp = sample(4)
             sample(4) = sample(5)
             sample(5) = temp
          end if
          !if (sample(6) > sample(7)) then
          ! temp = sample(6)
          ! sample(6) = sample(7)
          ! sample(7) = temp
          ! end if ! skipped
          !if (sample(8) > sample(9)) then
          ! temp = sample(8)
          ! sample(8) = sample(9)
          ! sample(9) = temp
          ! end if ! skipped
          pivot = sample(5)
          ! ************************
          left = 0
          right = nA + 1
          do while (left < right)
             right = right - 1
             do while (A(right) > pivot)
                right = right - 1
             end do
             left = left + 1
             do while (A(left) < pivot)
                left = left + 1
             end do
             if (left < right) then
                temp2 = A(left)
                A(left) = A(right)
                A(right) = temp2
                itemp2 = indx(left)
                indx(left) = indx(right)
                indx(right) = itemp2
             end if
          end do

          if (left == right) then
             marker = left + 1
          else
             marker = left
          end if

          call QuickSortA_i4(A(:marker-1),marker-1,indx(:marker-1),limit)
          call QuickSortA_i4(A(marker:),nA-marker+1,indx(marker:),limit)

       else
          call InsertionSortA_i4(A,nA,indx)    ! Insertion sort for small groups is faster than Quicksort
       end if
    end if

  end subroutine QuickSortA_i4

  subroutine InsertionSortA_i4(A,nA,indx)
    ! USED MODULES
    use mo_kind, only: i4, i4

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA
    integer(i4), dimension(nA), intent(inout) :: A
    integer(i4), dimension(nA), intent(inout) :: indx

    ! LOCAL VARIABLES
    integer(i4) :: temp
    integer(i4) :: i, j
    integer(i4) :: itemp

    outter: do i = 2, nA
       j = i - 1
       temp = A(i)
       itemp = indx(i)
       inner: do
          if (j == 0) exit inner
          if (a(j) <= temp) exit inner
          A(j+1) = A(j)
          indx(j+1) = indx(j)
          j = j - 1
       end do inner
       a(j+1) = temp
       indx(j+1) = itemp
    end do outter

  end subroutine InsertionSortA_i4

  ! MULTI-THREAD MERGE
  ! Divides array a into nt number of parts.  Lookup matching values in array b.  Merges the parts together into new memory array c.
  ! Ascending sort
  subroutine Merge8a_mt_i4(A,nA,B,nB,C,nC,iA,iB,iC,nt)

    ! USED MODULES
    use mo_kind, only: i4, i4
    !$ use omp_lib, only: omp_get_thread_num

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA, nB, nC                      ! Size of arrays.  Normal usage: nA+nB = nC
    integer(i4), dimension(nA), intent(in) :: A     ! left part
    integer(i4), dimension(nB), intent(in) :: B     ! right part
    integer(i4), dimension(nC), intent(out) :: C    ! sorted output must be new memory (not overlapping a or b)
    integer(i4), dimension(nA), intent(in) :: iA     ! left index part
    integer(i4), dimension(nB), intent(in) :: iB     ! right index part
    integer(i4), dimension(nC), intent(out) :: iC    ! sorted output must be new memory (not overlapping a or b)
    integer(i4), intent(in) :: nt                                           ! number of threads

    ! LOCAL VARIABLES
    integer(i4) :: t                                             ! thread number
    integer(i4), dimension(0:nt) :: dividerA        ! section divider of array a
    integer(i4), dimension(0:nt) :: dividerB        ! section divider of array b
    integer(i4), dimension(0:nt) :: dividerC        ! section divider of array c
    integer(i4) :: divider_value                      ! value of section dividers
    integer(i4) :: nA_part, nB_part                 ! number of left and right elements of group
    integer(i4) :: nC_part                          ! number of sorted elements in group

    t=1
    !$ call omp_set_nested(.true.)

    ! initialize first and last divider
    dividerA(0) = 0
    dividerB(0) = 0
    dividerB(nt) = nB
    dividerC(0) = 0
    dividerC(nt) = nC

    !$omp parallel &
    !$omp num_threads(nt) &
    !$omp default(none) &
    !$omp private (t, divider_value, nA_part, nB_part, nC_part) &
    !$omp shared (a, nA, b, nB, c, nC, iA, iB, iC, nt, dividerA, dividerB, dividerC)

    !$ t = omp_get_thread_num() + 1  ! add 1 to make first thread 1 instead of 0
    dividerA(t) = nA * t / nt   ! starting value for dividerA.  Will change if not last in group.
    if (t < nt) then
       if (a(dividerA(t)) == a(dividerA(t)-1)) then
          ! Get new DividerA because starting value is in group, and not the last in the group.
          dividerA(t) = LookupAscending_i4(a,nA,a(dividerA(t))) - 1
       end if
    end if
    ! POSSIBILITIES
    ! dividerA is zero: Previous groups took all A values (same number extends from one dividerA to the end of A)
    ! dividerA if equal to nA: All values after dividerA starting value are same for the rest of A.
    ! dividerA is between 0 and nA: There are values that exceed the dividerA starting value (normal operation)

    divider_value = a(dividerA(t))
    if (t < nt) then
       ! find closest match that is just less than divider_value in array b
       dividerB(t) = LookupAscending_i4(b,nB,divider_value) - 1
       dividerC(t) = dividerA(t) + dividerB(t)
    end if
    !$omp barrier   ! make sure all divider values are calculated because each thread needs previous thread's divider values

    nA_part = dividerA(t) - (dividerA(t-1)+1) + 1
    nB_part = dividerB(t) - (dividerB(t-1)+1) + 1
    nC_part = dividerC(t) - (dividerC(t-1)+1) + 1

    ! POSSIBILITIES
    ! 1. nA_part == 0: Caused by dividerA being equal to nA.  This happens when all values in A beyond starting value of dividerA are
    !    equal.  After first thread, copy none of A-part and all of B-part into C.  First thead does normal (3) operation with all of A
    !    and a part of B (nA_part is not equal to 0 for the first thread, only for all threads greater than one)
    ! 2. nB_part == 0: Happens if Lookup finds no value in B that is less than divider_value in A.  Copy A-part into C.
    ! 3. Neither nA_part or nB_part are equal to 0.  Merge A-part and B-part (normal operation).

    if (nA_part == 0) then      ! possibility 1
       if (nB_part > 0) then
          C(dividerC(t-1)+1:dividerC(t-1)+nB_part) = B(dividerB(t-1)+1:dividerB(t)) ! copy only B-part
          iC(dividerC(t-1)+1:dividerC(t-1)+nB_part) = iB(dividerB(t-1)+1:dividerB(t)) ! copy only B-part
       endif
    else if (nB_part == 0) then ! possibility 2
       C(dividerC(t-1)+1:dividerC(t)) = A(dividerA(t-1)+1:dividerA(t)) ! copy only A-part
       iC(dividerC(t-1)+1:dividerC(t)) = iA(dividerA(t-1)+1:dividerA(t)) ! copy only A-part
    else                        ! possibility 3
       call Merge8a_i4(   A(dividerA(t-1)+1:dividerA(t)), nA_part, &   ! A-part
            &               B(dividerB(t-1)+1:dividerB(t)), nB_part, &   ! B-part
            &               C(dividerC(t-1)+1:dividerC(t)), nC_part, &  ! sorted part
            &               iA(dividerA(t-1)+1:dividerA(t)), &   ! A-part
            &               iB(dividerB(t-1)+1:dividerB(t)), &   ! B-part
            &               iC(dividerC(t-1)+1:dividerC(t)))     ! sorted part
    end if

    !$omp end parallel

  end subroutine Merge8a_mt_i4

  ! find first location where value is exceeded in array a.  Array a must be sorted in ascending order.  Uses binary search.
  function LookupAscending_i4(a,nA,value)

    ! USED MODULES
    use mo_kind, only: i4, i4

    ! DUMMY ARGUMENTS
    integer(i4) :: LookupAscending_i4
    integer(i4), intent(in) :: nA         ! Size of a
    integer(i4), dimension(nA), intent(in) :: a
    integer(i4), intent(in) :: value

    ! LOCAL VARIABLES
    integer(i4) :: half ! half the distance between binary searches

    if (a(1) > value) then
       ! Note: If using this function for other purposes, you may want to set LookupAscending to zero when the first value exceeds
       ! the value to make the result distinct from the first value being the match, but returning that value causes the
       ! multi-threaded merge to fail.
       LookupAscending_i4 = 1    ! first value exceeds
    else
       if (a(nA) <= value) then
          LookupAscending_i4 = nA + 1 ! last value is too small
       else    ! value is in between.  Find it.
          LookupAscending_i4 = (nA+1)/2    ! starting location in middle
          half = LookupAscending_i4
          do
             half = (half+1) / 2
             if (a(LookupAscending_i4) > value) then
                if (a(LookupAscending_i4-1) <= value) then
                   exit    ! found
                else
                   LookupAscending_i4 = max(LookupAscending_i4 - half, 1) ! move down half
                end if
             else
                LookupAscending_i4 = min(LookupAscending_i4 + half, nA) ! move up half
             end if
          end do
       end if
    end if

  end function LookupAscending_i4

  ! Ascending merge (merges 2 ascending sorted lists into 1 ascending sorted list)
  subroutine Merge8a_i4(a,nA,b,nB,c,nC,ia,ib,ic)

    ! USED MODULES
    use mo_kind, only: i4, i4

    ! DUMMY ARGUMENTS
    integer(i4), intent(in) :: nA, nB, nC         ! Size of arrays.  Normal usage: nC = nA+nB
    !   integer(i4), dimension(:), intent(in) :: A ! Make sure A is a copy of any subset of C or sort overwrites as it goes
    !   integer(i4), dimension(:), intent(in) :: B        ! B overlays C(nA+1:nC)
    !   integer(i4), dimension(:), intent(out) :: C
    !    ! ***NOTE: USING EXPLICIT ARRAY SIZES AS SHOW BELOW SOMETIMES CAUSES SEGMENTATION FAULTS WITH LARGE ARRAY SIZES
    !   automatic arrays are faster than assumed size arrays.  Switch to above variable declarations if problems happen.
    integer(i4), dimension(nA), intent(in) :: A  ! left part
    integer(i4), dimension(nB), intent(in) :: B ! right part
    integer(i4), dimension(nC), intent(out) :: C ! output array
    integer(i4), dimension(nA), intent(in) :: iA  ! left part
    integer(i4), dimension(nB), intent(in) :: iB ! right part
    integer(i4), dimension(nC), intent(out) :: iC ! output array
    ! Note: for single-threaded merge, such as the call from MergeSort, array a is a copy of the left part of array c.  Also, B overlays
    ! array c(nA+1:nC).  As multi-threaded merge array a and b are parts of the same input array.  Multi-threaded usage also requires
    ! array c to be different memory not overlapping arrays a or b.

    ! LOCAL VARIABLES
    integer(i4) :: i, j, k

    i = 1
    j = 1
    k = 1
    do
       if (i > nA .or. j > nB) exit
       if (A(i) <= B(j)) then
          C(k) = A(i)
          iC(k) = iA(i)
          i = i + 1
       else
          C(k) = B(j)
          iC(k) = iB(j)
          j = j + 1
       end if
       k = k + 1
    end do
    if (i <= nA) then
       C(k:) = A(i:)
       iC(k:) = iA(i:)
       return
    end if
    if (j <= nB) then
       C(k:) = B(j:)    ! This statement is only necessary for multi-threaded merge
       iC(k:) = iB(j:)    ! This statement is only necessary for multi-threaded merge
    endif

  end subroutine Merge8a_i4

END MODULE mo_quicksort

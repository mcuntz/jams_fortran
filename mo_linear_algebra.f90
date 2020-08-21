!> \file mo_linear_algebra.f90

!> \brief Wrapper for LAPACK's linear algebra routines.

!> \details This modules provides mostly wrappers for LAPACK's F77 linear algebra routines.
!>          It adds a few convenience functions such as diag.

!> \authors Matthias Cuntz, Sebastian Mueller
!> \date May 2014
MODULE mo_linear_algebra

  ! Wrapper for LAPACK's F77 linear algebra routines.

  ! Written  Matthias Cuntz,    May 2014
  ! Modified Matthias Cuntz,    May 2016 - calc single precision via double precision
  ! Modified Sebastian Mueller, Oct 2016 - solver for banded coefficent-matrices and some involved algorithms
  ! Modified Matthias Cuntz,    Mar 2020 - allocate out of sp routines calling dp routines only if not allocated
  !                                      - allocate out for Python.

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2014-2020 Matthias Cuntz, Sebastian Mueller - mc (at) macu (dot) de
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

  USE mo_kind, ONLY: i4, i8, sp, dp

  IMPLICIT NONE

  PUBLIC :: banded                      ! banded form of a given matrix
  PUBLIC :: diag                        ! diagonal of matrix
  PUBLIC :: inverse                     ! inverse of matrix
  PUBLIC :: min_diag                    ! minor diagonal of matrix
  PUBLIC :: solve_linear_equations      ! solve linear system of equations with LU decomposition
  PUBLIC :: solve_linear_equations_svd  ! solve linear system of equations with SVD
  PUBLIC :: solve_linear_equations_band ! solve linear system of equations for banded matrix


  ! ------------------------------------------------------------------
  !
  !     NAME
  !         banded
  !
  !     PURPOSE
  !         Banded form of squared matrix for the banded linear equations solver
  !
  !>        \brief Banded form of squared matrix.
  !
  !>        \details Converts a given squared matrix into the banded form needed for the Lapack-solver for banded matrices.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: matrix(:,:)"    Squared 2D-array
  !>        \param[in] "integer(i4) :: l          "    number of lower minor-diagonals
  !>        \param[in] "integer(i4) :: u          "    number of upper minor-diagonals
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !>       \return     real(sp/dp) :: banded(:) &mdash; Banded form of input matrix
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         band = banded(matrix,l,u)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Sebastian Mueller
  !>        \date October 2016
  INTERFACE banded
     MODULE PROCEDURE banded_dp, banded_sp
  END INTERFACE banded

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         diag
  !
  !     PURPOSE
  !         Returns the diagonal of a square 2D-array.
  !
  !>        \brief Diagonal elements of a squared matrix.
  !
  !>        \details Returns the diagonal of a square 2D-array.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp)/integer(i4/i8)/logical :: matrix(:,:)"    Squared 2D-array
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !>       \return     real(sp/dp)/integer(i4/i8)/logical :: diag(:) &mdash; Diagonal elements of matrix
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         m = diag(matrix)
  !         -> see also example in test directory
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date May 2014
  INTERFACE diag
     MODULE PROCEDURE diag_sp, diag_dp, diag_i4, diag_i8, diag_lgt
  END INTERFACE diag

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         inverse
  !
  !     PURPOSE
  !         Inverse of squared matrix
  !
  !>        \brief Inverse of squared matrix.
  !
  !>        \details Inverts squared matrix using LU decomposition.
  !>
  !>                 Uses standard Lapack routines using LU decomposition: dgetrf and dgetri.
  !>                 Conditions columns before decomposition.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: matrix(:,:)"    Squared 2D-array
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: condition" If true, condition matrix before inversion (default: true)
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURNS
  !>       \return     real(sp/dp) :: inverse(:) &mdash; Inverse of input matrix
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         inv = inverse(matrix)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date May 2014
  INTERFACE inverse
     MODULE PROCEDURE inverse_dp, inverse_sp
  END INTERFACE inverse

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         min_diag
  !
  !     PURPOSE
  !         Returns the n-th minor-diagonal of a square 2D-array, where positive n indicate super-diagonals and negativ n indicate
  !         sub-diagonals.
  !
  !>        \brief Minor-diagonal elements of a squared matrix.
  !
  !>        \details Returns the n-th minor-diagonal of a square 2D-array.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp)/integer(i4/i8)/logical :: matrix(:,:)"    Squared 2D-array
  !>        \param[in] "integer(i4)                        :: n          "    number of minor-diagonal
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !>       \return     real(sp/dp)/integer(i4/i8)/logical :: min_diag(:) &mdash; n-th minor-diagonal elements of matrix
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         m = min_diag(matrix,n)
  !         -> see also example in test directory
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Sebastian Mueller
  !>        \date October 2016
  INTERFACE min_diag
     MODULE PROCEDURE min_diag_sp, min_diag_dp, min_diag_i4, min_diag_i8, min_diag_lgt
  END INTERFACE min_diag

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         solve_linear_equations
  !
  !     PURPOSE
  !         Solve linear system of equations
  !
  !>        \brief Solve linear system of equations.
  !
  !>        \details Solve linear system of equations using LU decomposition
  !>                 \f[ A x = b \f]
  !>                 Returns \f$ x \f$.
  !>
  !>                 Uses standard Lapack routine using LU decomposition: dgesv.
  !>                 Conditions rows and columns before decomposition.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: lhs(:,:)"    Coefficients of left hand side \f$ A \f$
  !>        \param[in] "real(sp/dp) :: rhs(:)"      Right hand side \f$ b \f$
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: condition" If true, condition lhs and rhs before decomposition (default: true)
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURNS
  !>       \return     real(sp/dp) :: x(:) &mdash; Solution x to \f$ A x = b \f$
  !
  !     RESTRICTIONS
  !         Only one right-hand side.
  !
  !     EXAMPLE
  !         sol = solve_linear_equations(lhs, rhs)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date May 2014
  INTERFACE solve_linear_equations
     MODULE PROCEDURE solve_linear_equations_1_dp, solve_linear_equations_1_sp
  END INTERFACE solve_linear_equations

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         solve_linear_equations_svd
  !
  !     PURPOSE
  !         Solve linear system of equations
  !
  !>        \brief Solve linear system of equations.
  !
  !>        \details Solve linear system of equations using singular value decomposition
  !>                 \f[ A x = b \f]
  !>                 Returns \f$ x \f$.
  !>
  !>                 Uses standard Lapack routine to do SVD: dgesvd.
  !>                 then solving routine similar to Numerical Recipes: svdksb.
  !>                 It conditions rows and columns before decomposition.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: lhs(:,:)"    Coefficients of left hand side \f$ A \f$
  !>        \param[in] "real(sp/dp) :: rhs(:)"      Right hand side \f$ b \f$
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: condition" If true, condition lhs and rhs before decomposition (default: true)
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURNS
  !>       \return     real(sp/dp) :: x(:) &mdash; Solution x to \f$ A x = b \f$
  !
  !     RESTRICTIONS
  !         Only one right-hand side.
  !
  !     EXAMPLE
  !         sol = solve_linear_equations_svd(lhs, rhs)
  !         -> see also example in test directory

  !     LITERATURE
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !>        \author Matthias Cuntz
  !>        \date May 2014
  INTERFACE solve_linear_equations_svd
     MODULE PROCEDURE solve_linear_equations_svd_1_dp, solve_linear_equations_svd_1_sp
  END INTERFACE solve_linear_equations_svd

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         solve_linear_equations_band
  !
  !     PURPOSE
  !         Solve linear system of equations for a banded matrix A=(aij).
  !         The band storage scheme is illustrated by the following example, when
  !         N = 6, l = 2, u = 1:
  !
  !             *   a12  a23  a34  a45  a56
  !            a11  a22  a33  a44  a55  a66
  !            a21  a32  a43  a54  a65   *
  !            a31  a42  a53  a64   *    *
  !
  !         Array elements marked * are not used by the routine. You can set them zero.
  !
  !>        \brief Solve linear system of equations.
  !
  !>        \details Solve linear system of equations for banded coefficent-matrix
  !>                 \f[ A x = b \f]
  !>                 Returns \f$ x \f$.
  !>
  !>                 Uses standard Lapack routine for banded matrix: dgbsv.
  !>                 Conditions columns before decomposition.
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: lhsb(:,:)"    Coefficients of left hand side in banded form of \f$ A \f$
  !>        \param[in] "real(sp/dp) :: rhs(:)   "    Right hand side \f$ b \f$
  !>        \param[in] "integer(i4) :: l        "    number of lower minor-diagonals
  !>        \param[in] "integer(i4) :: u        "    number of upper minor-diagonals
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "logical, optional :: condition" If true, condition lhsb before decomposition (default: true)
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURNS
  !>       \return     real(sp/dp) :: x(:) &mdash; Solution x to \f$ A x = b \f$
  !
  !     RESTRICTIONS
  !         Only one right-hand side.
  !         Left-hand side needs to be in banded form. You can use the banded-function to convert a given squared matrix.
  !
  !     EXAMPLE
  !         sol = solve_linear_equations_band(lhsb, rhs, l, u)
  !         sol = solve_linear_equations_band(banded(lhs), rhs, l, u)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Sebastian Mueller
  !>        \date October 2016
  !         Modified, Matthias Cuntz, May 2018 - allocate output for solve_linear_equations_band_1_sp
  INTERFACE solve_linear_equations_band
     MODULE PROCEDURE solve_linear_equations_band_1_dp, solve_linear_equations_band_1_sp
  END INTERFACE solve_linear_equations_band

  PRIVATE

  INTERFACE svdksb
     MODULE PROCEDURE svdksb_dp, svdksb_sp
  END INTERFACE svdksb

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION banded_dp(matrix, l, u)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN)   :: matrix
    INTEGER(i4),              INTENT(IN)   :: l, u
#ifndef __PYTHON__
    REAL(dp), DIMENSION(:,:), allocatable  :: banded_dp
#else
    REAL(dp), DIMENSION(l+u+1,size(matrix,1)) :: banded_dp
#endif

    INTEGER(i4) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'banded_dp: array must be squared matrix.'
    if (l              >= size(matrix,1)) stop 'banded_dp: l is to big. You need to choose minor diagonals within the matrix.'
    if (u              >= size(matrix,1)) stop 'banded_dp: u is to big. You need to choose minor diagonals within the matrix.'
    if (l              <  0_i4          ) stop 'banded_dp: l needs to be non-negativ.'
    if (u              <  0_i4          ) stop 'banded_dp: u needs to be non-negativ.'
#ifndef __PYTHON__
    if (.not. allocated(banded_dp)) allocate(banded_dp(l+u+1,size(matrix,1)))
#endif

    banded_dp = 0.0_dp

    do i=0_i4,u
        banded_dp(u+1-i,i+1:size(matrix,1)) = min_diag_dp(matrix, i)
    end do

    do i=1_i4,l
        banded_dp(u+1+i,1:size(matrix,1)-i) = min_diag_dp(matrix,-i)
    end do

  END FUNCTION banded_dp

  FUNCTION banded_sp(matrix, l, u)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN)   :: matrix
    INTEGER(i4),              INTENT(IN)   :: l, u
#ifndef __PYTHON__
    REAL(sp), DIMENSION(:,:), allocatable  :: banded_sp
#else
    REAL(sp), DIMENSION(l+u+1,size(matrix,1)) :: banded_sp
#endif

    INTEGER(i4) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'banded_sp: array must be squared matrix.'
    if (l              >= size(matrix,1)) stop 'banded_sp: l is to big. You need to choose minor diagonals within the matrix.'
    if (u              >= size(matrix,1)) stop 'banded_sp: u is to big. You need to choose minor diagonals within the matrix.'
    if (l              <  0_i4          ) stop 'banded_sp: l needs to be non-negativ.'
    if (u              <  0_i4          ) stop 'banded_sp: u needs to be non-negativ.'
#ifndef __PYTHON__
    if (.not. allocated(banded_sp)) allocate(banded_sp(l+u+1,size(matrix,1)))
#endif

    banded_sp = 0.0_sp

    do i=0_i4,u
        banded_sp(u+1-i,i+1:size(matrix,1)) = min_diag_sp(matrix, i)
    end do

    do i=1_i4,l
        banded_sp(u+1+i,1:size(matrix,1)-i) = min_diag_sp(matrix,-i)
    end do

  END FUNCTION banded_sp

  ! ------------------------------------------------------------------

  FUNCTION diag_dp(matrix)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: matrix
#ifndef __PYTHON__
    REAL(dp), DIMENSION(:), allocatable  :: diag_dp
#else
    REAL(dp), DIMENSION(size(matrix,1)) :: diag_dp
#endif

    INTEGER(i4) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_dp: array must be squared matrix.'
#ifndef __PYTHON__
    if (.not. allocated(diag_dp)) allocate(diag_dp(size(matrix,1)))
#endif

    forall(i=1:size(matrix,1)) diag_dp(i) = matrix(i,i)

  END FUNCTION diag_dp

  FUNCTION diag_sp(matrix)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN) :: matrix
#ifndef __PYTHON__
    REAL(sp), DIMENSION(:), allocatable  :: diag_sp
#else
    REAL(sp), DIMENSION(size(matrix,1)) :: diag_sp
#endif

    INTEGER(i4) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_sp: array must be squared matrix.'
#ifndef __PYTHON__
    if (.not. allocated(diag_sp)) allocate(diag_sp(size(matrix,1)))
#endif

    forall(i=1:size(matrix,1)) diag_sp(i) = matrix(i,i)

  END FUNCTION diag_sp

  FUNCTION diag_i4(matrix)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:,:), INTENT(IN) :: matrix
#ifndef __PYTHON__
    INTEGER(i4), DIMENSION(:), allocatable  :: diag_i4
#else
    INTEGER(i4), DIMENSION(size(matrix,1)) :: diag_i4
#endif

    INTEGER(i4) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_i4: array must be squared matrix.'
#ifndef __PYTHON__
    if (.not. allocated(diag_i4)) allocate(diag_i4(size(matrix,1)))
#endif

    forall(i=1:size(matrix,1)) diag_i4(i) = matrix(i,i)

  END FUNCTION diag_i4

  FUNCTION diag_i8(matrix)

    IMPLICIT NONE

    INTEGER(i8), DIMENSION(:,:), INTENT(IN) :: matrix
#ifndef __PYTHON__
    INTEGER(i8), DIMENSION(:), allocatable  :: diag_i8
#else
    INTEGER(i8), DIMENSION(size(matrix,1)) :: diag_i8
#endif

    INTEGER(i8) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_i8: array must be squared matrix.'
#ifndef __PYTHON__
    if (.not. allocated(diag_i8)) allocate(diag_i8(size(matrix,1)))
#endif

    forall(i=1:size(matrix,1)) diag_i8(i) = matrix(i,i)

  END FUNCTION diag_i8

  FUNCTION diag_lgt(matrix)

    IMPLICIT NONE

    LOGICAL, DIMENSION(:,:), INTENT(IN) :: matrix
#ifndef __PYTHON__
    LOGICAL, DIMENSION(:), allocatable  :: diag_lgt
#else
    LOGICAL, DIMENSION(size(matrix,1)) :: diag_lgt
#endif

    INTEGER(i4) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_lgt: array must be squared matrix.'
#ifndef __PYTHON__
    if (.not. allocated(diag_lgt)) allocate(diag_lgt(size(matrix,1)))
#endif

    forall(i=1:size(matrix,1)) diag_lgt(i) = matrix(i,i)

  END FUNCTION diag_lgt

  ! ------------------------------------------------------------------

  FUNCTION inverse_dp(matrix, condition)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN)  :: matrix
    LOGICAL,  OPTIONAL,       INTENT(IN)  :: condition
#ifndef __PYTHON__
    REAL(dp), DIMENSION(:,:), allocatable :: inverse_dp
#else
    REAL(dp), DIMENSION(size(matrix,2),size(matrix,1)) :: inverse_dp
#endif

    INTEGER(i4) :: i, nn
    real(dp),    dimension(:),   allocatable :: scale_cols    ! scale matrix for better conditioning
    integer(i4), dimension(:),   allocatable :: ipiv          ! needed for dgetrf lapack routine
    integer(i4)                              :: info          !              "
    real(dp),    dimension(:),   allocatable :: work          ! needed for dgetri lapack routine
    integer(i4)                              :: lwork         !              "
    logical :: icondition

    external :: dgetrf ! Lapack routine to compute LU factorization of a general matrix
    external :: dgetri ! Lapack routine to compute inverse of matrix using the LU factorization computed by DGETRF

    nn = size(matrix,2)
    if (size(matrix,1) /= nn) stop 'inverse_dp: matrix must be square.'
#ifndef __PYTHON__
    if (.not. allocated(inverse_dp)) allocate(inverse_dp(nn,nn))
#endif
    inverse_dp = matrix

    if (present(condition)) then
       icondition = condition
    else
       icondition = .true.
    endif

    if (icondition) then
       ! Condition columns
       allocate(scale_cols(nn))
       scale_cols(:) = maxval(abs(inverse_dp(:,:)),1)
       where (scale_cols(:) < tiny(1.0_dp)) scale_cols(:) = 1.0_dp
       scale_cols(:) = 1.0_dp / scale_cols(:)
       forall(i=1:nn) inverse_dp(:,i) = inverse_dp(:,i) * scale_cols(i)
    endif

    ! LU factorisation of imatrix
    allocate(ipiv(nn))
    call dgetrf(nn, nn, inverse_dp, nn, ipiv, info)
    if (info /= 0) stop 'inverse_dp: LU factorisation did not work.'
    ! Inverse of LU factorisation of imatrix
    allocate(work(1))
    call dgetri(nn, inverse_dp, nn, ipiv, work, -1, info)
    lwork = int(work(1),i4)
    deallocate(work)
    allocate(work(lwork))
    call dgetri(nn, inverse_dp, nn, ipiv, work, lwork, info)
    if (info /= 0) stop 'hdmr_hessian: Inversion did not work.'

    if (icondition) then
       ! rescale result
       forall(i=1:nn) inverse_dp(i,:) = inverse_dp(i,:) * scale_cols(i)

       deallocate(scale_cols)
    endif

    deallocate(ipiv)

  END FUNCTION inverse_dp


  FUNCTION inverse_sp(matrix, condition)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN)  :: matrix
    LOGICAL,  OPTIONAL,       INTENT(IN)  :: condition
#ifndef __PYTHON__
    REAL(sp), DIMENSION(:,:), allocatable :: inverse_sp
#else
    REAL(sp), DIMENSION(size(matrix,2),size(matrix,1)) :: inverse_sp
#endif

#ifndef __PYTHON__
    if (.not. allocated(inverse_sp)) allocate(inverse_sp(size(matrix,2),size(matrix,1)))
#endif
    inverse_sp = real(inverse_dp(real(matrix,dp), condition), sp)

  END FUNCTION inverse_sp

  ! ------------------------------------------------------------------

  FUNCTION min_diag_dp(matrix, n)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: matrix
    INTEGER(i4),              INTENT(IN) :: n
#ifndef __PYTHON__
    REAL(dp), DIMENSION(:), allocatable  :: min_diag_dp
#else
    REAL(dp), DIMENSION(size(matrix,1)-abs(n)) :: min_diag_dp
#endif

    if (size(matrix,1) /= size(matrix,2)) stop 'min_diag_dp: array must be squared matrix.'
    if (abs(n)         >= size(matrix,1)) stop 'min_diag_dp: n is to big. You need to choose a minor diagonal within the matrix.'
#ifndef __PYTHON__
    if (.not. allocated(min_diag_dp)) allocate(min_diag_dp(size(matrix,1) - abs(n)))
#endif

    if (n >= 0_i4) then
        min_diag_dp = diag_dp(matrix(1_i4:size(matrix,1)-n,1_i4+n:size(matrix,1)))
    else
        min_diag_dp = diag_dp(matrix(1_i4-n:size(matrix,1),1_i4:size(matrix,1)+n))
    end if

  END FUNCTION min_diag_dp

  FUNCTION min_diag_sp(matrix, n)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN) :: matrix
    INTEGER(i4),              INTENT(IN) :: n
#ifndef __PYTHON__
    REAL(sp), DIMENSION(:), allocatable  :: min_diag_sp
#else
    REAL(sp), DIMENSION(size(matrix,1)-abs(n)) :: min_diag_sp
#endif

    if (size(matrix,1) /= size(matrix,2)) stop 'min_diag_sp: array must be squared matrix.'
    if (abs(n)         >= size(matrix,1)) stop 'min_diag_sp: n is to big. You need to choose a minor diagonal within the matrix.'
#ifndef __PYTHON__
    if (.not. allocated(min_diag_sp)) allocate(min_diag_sp(size(matrix,1) - abs(n)))
#endif

    if (n >= 0_i4) then
        min_diag_sp = diag_sp(matrix(1_i4:size(matrix,1)-n,1_i4+n:size(matrix,1)))
    else
        min_diag_sp = diag_sp(matrix(1_i4-n:size(matrix,1),1_i4:size(matrix,1)+n))
    end if

  END FUNCTION min_diag_sp

  FUNCTION min_diag_i4(matrix, n)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:,:), INTENT(IN) :: matrix
    INTEGER(i4),                 INTENT(IN) :: n
#ifndef __PYTHON__
    INTEGER(i4), DIMENSION(:), allocatable  :: min_diag_i4
#else
    INTEGER(i4), DIMENSION(size(matrix,1)-abs(n)) :: min_diag_i4
#endif

    if (size(matrix,1) /= size(matrix,2)) stop 'min_diag_i4: array must be squared matrix.'
    if (abs(n)         >= size(matrix,1)) stop 'min_diag_i4: n is to big. You need to choose a minor diagonal within the matrix.'
#ifndef __PYTHON__
    if (.not. allocated(min_diag_i4)) allocate(min_diag_i4(size(matrix,1) - abs(n)))
#endif

    if (n >= 0_i4) then
        min_diag_i4 = diag_i4(matrix(1_i4:size(matrix,1)-n,1_i4+n:size(matrix,1)))
    else
        min_diag_i4 = diag_i4(matrix(1_i4-n:size(matrix,1),1_i4:size(matrix,1)+n))
    end if

  END FUNCTION min_diag_i4

  FUNCTION min_diag_i8(matrix, n)

    IMPLICIT NONE

    INTEGER(i8), DIMENSION(:,:), INTENT(IN) :: matrix
    INTEGER(i4),                 INTENT(IN) :: n
#ifndef __PYTHON__
    INTEGER(i8), DIMENSION(:), allocatable  :: min_diag_i8
#else
    INTEGER(i8), DIMENSION(size(matrix,1)-abs(n)) :: min_diag_i8
#endif

    if (size(matrix,1) /= size(matrix,2)) stop 'min_diag_i8: array must be squared matrix.'
    if (abs(n)         >= size(matrix,1)) stop 'min_diag_i8: n is to big. You need to choose a minor diagonal within the matrix.'
#ifndef __PYTHON__
    if (.not. allocated(min_diag_i8)) allocate(min_diag_i8(size(matrix,1) - abs(n)))
#endif

    if (n >= 0_i4) then
        min_diag_i8 = diag_i8(matrix(1_i4:size(matrix,1)-n,1_i4+n:size(matrix,1)))
    else
        min_diag_i8 = diag_i8(matrix(1_i4-n:size(matrix,1),1_i4:size(matrix,1)+n))
    end if

  END FUNCTION min_diag_i8

  FUNCTION min_diag_lgt(matrix, n)

    IMPLICIT NONE

    LOGICAL, DIMENSION(:,:), INTENT(IN) :: matrix
    INTEGER(i4),             INTENT(IN) :: n
#ifndef __PYTHON__
    LOGICAL, DIMENSION(:), allocatable  :: min_diag_lgt
#else
    LOGICAL, DIMENSION(size(matrix,1)-abs(n)) :: min_diag_lgt
#endif

    if (size(matrix,1) /= size(matrix,2)) stop 'min_diag_lgt: array must be squared matrix.'
    if (abs(n)         >= size(matrix,1)) stop 'min_diag_lgt: n is to big. You need to choose a minor diagonal within the matrix.'
#ifndef __PYTHON__
    if (.not. allocated(min_diag_lgt)) allocate(min_diag_lgt(size(matrix,1) - abs(n)))
#endif

    if (n >= 0_i4) then
        min_diag_lgt = diag_lgt(matrix(1_i4:size(matrix,1)-n,1_i4+n:size(matrix,1)))
    else
        min_diag_lgt = diag_lgt(matrix(1_i4-n:size(matrix,1),1_i4:size(matrix,1)+n))
    end if

  END FUNCTION min_diag_lgt

  ! ------------------------------------------------------------------

  FUNCTION solve_linear_equations_1_dp(lhs, rhs, condition)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: lhs
    REAL(dp), DIMENSION(:),   INTENT(IN) :: rhs
    LOGICAL,  OPTIONAL,       INTENT(IN) :: condition
#ifndef __PYTHON__
    REAL(dp), DIMENSION(:), allocatable  :: solve_linear_equations_1_dp
#else
    REAL(dp), DIMENSION(size(lhs,1)) :: solve_linear_equations_1_dp
#endif

    INTEGER(i4) :: ii, nZeilen
    real(dp),    dimension(:,:), allocatable :: ilhs          ! internal lhs
    real(dp),    dimension(:),   allocatable :: irhs          ! internal rhs
    real(dp),    dimension(:),   allocatable :: scale_cols    ! scale matrix for better conditioning
    real(dp),    dimension(:),   allocatable :: scale_rows    !              "
    integer(i4), dimension(:),   allocatable :: ipiv          ! needed for dgesv lapack routine
    integer(i4)                              :: info          !              "
    logical :: icondition

    external :: dgesv  ! Lapack routine to compute solution of real system of linear equations

    nZeilen = size(lhs,1)
    if (size(lhs,2) /= nZeilen) stop 'solve_linear_equations_1_dp: left hand side must be squared matrix.'
    if (size(rhs,1) /= nZeilen) stop 'solve_linear_equations_1_dp: right hand side must have same size as left hand side.'
    ! internal arrays
    allocate(ilhs(nZeilen,nZeilen), irhs(nZeilen))
    ilhs = lhs
    irhs = rhs

    if (present(condition)) then
       icondition = condition
    else
       icondition = .true.
    endif

    if (icondition) then
       ! Condition of matrix
       allocate(scale_cols(nZeilen), scale_rows(nZeilen))
       ! Condition columns
       scale_cols(:) = maxval(abs(ilhs(:,:)),1)
       where (scale_cols(:) < tiny(1.0_dp)) scale_cols(:) = 1.0_dp
       scale_cols(:) = 1.0_dp / scale_cols(:)
       forall(ii=1:nZeilen) ilhs(:,ii) = ilhs(:,ii) * scale_cols(ii)
       ! Condition rows
       scale_rows(:) = maxval(abs(ilhs(:,:)),2)
       where (scale_rows(:) < tiny(1.0_dp)) scale_rows(:) = 1.0_dp
       scale_rows(:) = 1.0_dp / scale_rows(:)
       forall(ii=1:nZeilen) ilhs(ii,:) = ilhs(ii,:) * scale_rows(ii)
       irhs(:) = irhs(:) * scale_rows(:)
    endif

    ! solve linear system of equations
    allocate(ipiv(nZeilen))
    call dgesv(nZeilen, 1, ilhs, nZeilen, ipiv, irhs, nZeilen, info)
    if (info /= 0) stop 'solve_linear_equations_1_dp: Solving of linear system did not work.'

#ifndef __PYTHON__
    if (.not. allocated(solve_linear_equations_1_dp)) allocate(solve_linear_equations_1_dp(nZeilen))
#endif
    solve_linear_equations_1_dp = irhs

    if (icondition) then
       ! rescale result
       solve_linear_equations_1_dp(:) = solve_linear_equations_1_dp(:) * scale_cols(:)

       deallocate(scale_cols, scale_rows)
    endif

    deallocate(ilhs, irhs)
    deallocate(ipiv)

  END FUNCTION solve_linear_equations_1_dp


  FUNCTION solve_linear_equations_1_sp(lhs, rhs, condition)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN) :: lhs
    REAL(sp), DIMENSION(:),   INTENT(IN) :: rhs
    LOGICAL,  OPTIONAL,       INTENT(IN) :: condition
#ifndef __PYTHON__
    REAL(sp), DIMENSION(:), allocatable  :: solve_linear_equations_1_sp
#else
    REAL(sp), DIMENSION(size(rhs,1)) :: solve_linear_equations_1_sp
#endif

#ifndef __PYTHON__
    if (.not. allocated(solve_linear_equations_1_sp)) allocate(solve_linear_equations_1_sp(size(rhs,1)))
#endif
    solve_linear_equations_1_sp = real(solve_linear_equations_1_dp(real(lhs,dp), real(rhs,dp), condition), sp)

  END FUNCTION solve_linear_equations_1_sp

  ! ------------------------------------------------------------------

  FUNCTION solve_linear_equations_svd_1_dp(lhs, rhs, condition)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: lhs
    REAL(dp), DIMENSION(:),   INTENT(IN) :: rhs
    LOGICAL,  OPTIONAL,       INTENT(IN) :: condition
#ifndef __PYTHON__
    REAL(dp), DIMENSION(:), allocatable  :: solve_linear_equations_svd_1_dp
#else
    REAL(dp), DIMENSION(size(lhs,1)) :: solve_linear_equations_svd_1_dp
#endif

    INTEGER(i4) :: ii, nZeilen
    real(dp),    dimension(:,:), allocatable :: ilhs          ! internal lhs
    real(dp),    dimension(:),   allocatable :: irhs          ! internal rhs
    real(dp),    dimension(:),   allocatable :: scale_cols    ! scale matrix for better conditioning
    real(dp),    dimension(:),   allocatable :: scale_rows    !              "
    real(dp),    dimension(:),   allocatable :: work          ! needed for dgesvd lapack routine
    real(dp),    dimension(:),   allocatable :: svdw          !              "
    real(dp),    dimension(:,:), allocatable :: svdu          !              "
    real(dp),    dimension(:,:), allocatable :: svdv          !              "
    integer(i4)                              :: lwork         !              "
    integer(i4)                              :: info          !              "
    real(dp),    parameter                   :: svdtol = 1.0e-5_dp ! if <svdtol*maxval(svdw), then set svdw=0
    logical :: icondition

    external :: dgesvd ! Lapack routine to compute singular value decomposition of matrix

    nZeilen = size(lhs,1)
    if (size(lhs,2) /= nZeilen) stop 'solve_linear_equations_svd_1_dp: left hand side must be squared matrix.'
    if (size(rhs,1) /= nZeilen) stop 'solve_linear_equations_svd_1_dp: right hand side must have same size as left hand side.'
    ! internal arrays
    allocate(ilhs(nZeilen,nZeilen), irhs(nZeilen))
    ilhs = lhs
    irhs = rhs

    if (present(condition)) then
       icondition = condition
    else
       icondition = .true.
    endif

    if (icondition) then
       ! Condition of matrix
       allocate(scale_cols(nZeilen), scale_rows(nZeilen))
       ! Condition columns
       scale_cols(:) = maxval(abs(ilhs(:,:)),1)
       where (scale_cols(:) < tiny(1.0_dp)) scale_cols(:) = 1.0_dp
       scale_cols(:) = 1.0_dp / scale_cols(:)
       forall(ii=1:nZeilen) ilhs(:,ii) = ilhs(:,ii) * scale_cols(ii)
       ! Condition rows
       scale_rows(:) = maxval(abs(ilhs(:,:)),2)
       where (scale_rows(:) < tiny(1.0_dp)) scale_rows(:) = 1.0_dp
       scale_rows(:) = 1.0_dp / scale_rows(:)
       forall(ii=1:nZeilen) ilhs(ii,:) = ilhs(ii,:) * scale_rows(ii)
       irhs(:) = irhs(:) * scale_rows(:)
    endif

    ! solve linear system of equations
    ! Use Lapack singular value decomposition and Numerical Recipes solution
    allocate(work(1))
    allocate(svdu(nZeilen,nZeilen), svdv(nZeilen,nZeilen), svdw(nZeilen))
    call dgesvd('A', 'A', nZeilen, nZeilen, ilhs, nZeilen, svdw, svdu, nZeilen, svdv, nZeilen, work, -1, info)
    lwork = int(work(1),i4)
    deallocate(work)
    allocate(work(lwork))
    call dgesvd('A', 'A', nZeilen, nZeilen, ilhs, nZeilen, svdw, svdu, nZeilen, svdv, nZeilen, work, lwork, info)
    if (info /= 0) stop 'solve_linear_equations_svd_1_dp: Solving of linear system did not work.'
    ! write(*,*) 'Matrix condition ', maxval(abs(svdw))/minval(abs(svdw))
    where (svdw < svdtol*maxval(svdw)) svdw = 0.0_dp

#ifndef __PYTHON__
    if (.not. allocated(solve_linear_equations_svd_1_dp)) allocate(solve_linear_equations_svd_1_dp(nZeilen))
#endif
    ! svdv is V**T from dgesvd, svdksb wants V
    solve_linear_equations_svd_1_dp = svdksb(svdu,svdw,transpose(svdv),irhs)

    if (icondition) then
       ! rescale result
       solve_linear_equations_svd_1_dp(:) = solve_linear_equations_svd_1_dp(:) * scale_cols(:)

       deallocate(scale_cols, scale_rows)
    endif

    deallocate(ilhs, irhs)
    deallocate(work, svdu, svdv, svdw)

  END FUNCTION solve_linear_equations_svd_1_dp


  FUNCTION solve_linear_equations_svd_1_sp(lhs, rhs, condition)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN) :: lhs
    REAL(sp), DIMENSION(:),   INTENT(IN) :: rhs
    LOGICAL,  OPTIONAL,       INTENT(IN) :: condition
#ifndef __PYTHON__
    REAL(sp), DIMENSION(:), allocatable  :: solve_linear_equations_svd_1_sp
#else
    REAL(sp), DIMENSION(size(rhs,1)) :: solve_linear_equations_svd_1_sp
#endif

#ifndef __PYTHON__
    if (.not. allocated(solve_linear_equations_svd_1_sp)) allocate(solve_linear_equations_svd_1_sp(size(rhs,1)))
#endif
    solve_linear_equations_svd_1_sp = real(solve_linear_equations_svd_1_dp(real(lhs,dp), real(rhs,dp), condition), sp)

  END FUNCTION solve_linear_equations_svd_1_sp

  ! ------------------------------------------------------------------

  FUNCTION solve_linear_equations_band_1_dp(lhsb, rhs, l, u, condition)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: lhsb !banded form of the lhs
    REAL(dp), DIMENSION(:),   INTENT(IN) :: rhs
    INTEGER(i4),              INTENT(IN) :: l, u
    LOGICAL,  OPTIONAL,       INTENT(IN) :: condition
#ifndef __PYTHON__
    REAL(dp), DIMENSION(:), allocatable  :: solve_linear_equations_band_1_dp
#else
    REAL(dp), DIMENSION(size(rhs,1)) :: solve_linear_equations_band_1_dp
#endif

    INTEGER(i4) :: ii, nSpalten, nZeilen
    real(dp),    dimension(:,:), allocatable :: ilhsb         ! internal lhsb
    real(dp),    dimension(:),   allocatable :: irhs          ! internal rhs
    real(dp),    dimension(:),   allocatable :: scale_cols    ! scale matrix for better conditioning
!    real(dp),    dimension(:),   allocatable :: scale_rows    !              "
    integer(i4), dimension(:),   allocatable :: ipiv          ! needed for dgesv lapack routine
    integer(i4)                              :: info          !              "
    logical :: icondition

    external :: dgbsv  ! Lapack routine to compute solution of real system of linear equations with banded matrix

    nZeilen  = size(lhsb,1)
    nSpalten = size(lhsb,2)
    if (l+u+1_i4    /= nZeilen ) stop 'solve_linear_equations_1_dp: Given banded matrix must have l+u+1 rows.'
    if (size(rhs,1) /= nSpalten) stop 'solve_linear_equations_1_dp: right hand side must have same size as left hand side.'
    ! internal arrays
    allocate(ilhsb(nZeilen + l,nSpalten), irhs(nSpalten))
    ilhsb = 0.0_dp

    ilhsb(1+l:nZeilen + l,:) = lhsb
    irhs                     = rhs

    if (present(condition)) then
       icondition = condition
    else
       icondition = .true.
    endif

    if (icondition) then
       ! Condition of matrix
       allocate(scale_cols(nSpalten))
       ! Condition columns
       scale_cols(:) = maxval(abs(ilhsb(:,:)),1)
       where (scale_cols(:) < tiny(1.0_dp)) scale_cols(:) = 1.0_dp
       scale_cols(:) = 1.0_dp / scale_cols(:)
       forall(ii=1:nSpalten) ilhsb(:,ii) = ilhsb(:,ii) * scale_cols(ii)
    endif

    ! solve linear system of equations
    allocate(ipiv(nSpalten))
    call dgbsv(nSpalten, l, u, 1, ilhsb, nZeilen+l, ipiv, irhs, nSpalten, info)
    if (info /= 0) stop 'solve_linear_equations_band_1_dp: Solving of linear system did not work.'

#ifndef __PYTHON__
    if (.not. allocated(solve_linear_equations_band_1_dp)) allocate(solve_linear_equations_band_1_dp(nSpalten))
#endif
    solve_linear_equations_band_1_dp = irhs

    if (icondition) then
       ! rescale result
       solve_linear_equations_band_1_dp(:) = solve_linear_equations_band_1_dp(:) * scale_cols(:)

       deallocate(scale_cols)
    endif

    deallocate(ilhsb, irhs)
    deallocate(ipiv)

  END FUNCTION solve_linear_equations_band_1_dp

  FUNCTION solve_linear_equations_band_1_sp(lhsb, rhs, l, u, condition)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN) :: lhsb !banded form of the lhs
    REAL(sp), DIMENSION(:),   INTENT(IN) :: rhs
    INTEGER(i4),              INTENT(IN) :: l, u
    LOGICAL,  OPTIONAL,       INTENT(IN) :: condition
#ifndef __PYTHON__
    REAL(sp), DIMENSION(:), allocatable  :: solve_linear_equations_band_1_sp
#else
    REAL(sp), DIMENSION(size(rhs,1)) :: solve_linear_equations_band_1_sp
#endif

#ifndef __PYTHON__
    if (.not. allocated(solve_linear_equations_band_1_sp)) allocate(solve_linear_equations_band_1_sp(size(rhs,1)))
#endif
    solve_linear_equations_band_1_sp = real(solve_linear_equations_band_1_dp(real(lhsb,dp), real(rhs,dp), l, u, condition), sp)

  END FUNCTION solve_linear_equations_band_1_sp

  ! ------------------------------------------------------------------

  ! From Numerical Recipes in F77 book
  ! Solves A*X=B for a vector X, where A is specified by the arrays u, w, v as returned by
  ! svdcmp. m and n are the logical dimensions of a, and will be equal for square matrices. mp
  ! and np are the physical dimensions of a. b(1:m) is the input right-handside. x(1:n) is
  ! the output solution vector. No input quantities are destroyed, so the routine may be called
  ! sequentially with different bs.
  FUNCTION svdksb_dp(u,w,v,b)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: u,v
    REAL(dp), DIMENSION(:),   INTENT(IN) :: w,b
    REAL(dp), DIMENSION(size(u,2))       :: svdksb_dp

    REAL(dp), DIMENSION(size(u,2)) :: tmp
    REAL(dp), DIMENSION(size(w))   :: wtmp

    wtmp = w
    where (abs(w) < tiny(1.0_dp)) wtmp = 1.0_dp
    tmp = matmul(b,u) / wtmp
    where (abs(w) < tiny(1.0_dp)) tmp = 0.0_dp
    svdksb_dp = matmul(v,tmp)

  END FUNCTION svdksb_dp

  FUNCTION svdksb_sp(u,w,v,b)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN) :: u,v
    REAL(sp), DIMENSION(:),   INTENT(IN) :: w,b
    REAL(sp), DIMENSION(size(u,2))       :: svdksb_sp

    REAL(sp), DIMENSION(size(u,2)) :: tmp
    REAL(sp), DIMENSION(size(w))   :: wtmp

    wtmp = w
    where (abs(w) < tiny(1.0_sp)) wtmp = 1.0_sp
    tmp = matmul(b,u) / wtmp
    where (abs(w) < tiny(1.0_sp)) tmp = 0.0_sp
    svdksb_sp = matmul(v,tmp)

  END FUNCTION svdksb_sp

  ! ------------------------------------------------------------------

END MODULE mo_linear_algebra

!> \file mo_linear_algebra.f90

!> \brief Wrapper for LAPACK's linear algebra routines.

!> \details This modules provides mostly wrappers for LAPACK's F77 linear algebra routines.
!>          It adds a few convenience functions such as diag.

!> \authors Matthias Cuntz
!> \date May 2014
MODULE mo_linear_algebra

  ! Wrapper for LAPACK's F77 linear algebra routines.

  ! Written  Matthias Cuntz, May 2014

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
  ! along with the UFZ Fortran library (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2014 Matthias Cuntz

  USE mo_kind, ONLY: i4, i8, sp, dp

  IMPLICIT NONE

  PUBLIC :: diag                       ! diagonal of matrix
  PUBLIC :: inverse                    ! inverse of matrix
  PUBLIC :: solve_linear_equations     ! solve linear system of equations with LU decomposition
  PUBLIC :: solve_linear_equations_svd ! solve linear system of equations with SVD

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

  PRIVATE

  INTERFACE svdksb
     MODULE PROCEDURE svdksb_dp, svdksb_sp
  END INTERFACE svdksb

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION diag_dp(matrix)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: matrix
    REAL(dp), DIMENSION(:), allocatable  :: diag_dp

    INTEGER(i4) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_dp: array must be squared matrix.'
    if (.not. allocated(diag_dp)) allocate(diag_dp(size(matrix,1)))

    forall(i=1:size(matrix,1)) diag_dp(i) = matrix(i,i)

  END FUNCTION diag_dp

  FUNCTION diag_sp(matrix)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN) :: matrix
    REAL(sp), DIMENSION(:), allocatable  :: diag_sp

    INTEGER(i4) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_sp: array must be squared matrix.'
    if (.not. allocated(diag_sp)) allocate(diag_sp(size(matrix,1)))

    forall(i=1:size(matrix,1)) diag_sp(i) = matrix(i,i)

  END FUNCTION diag_sp

  FUNCTION diag_i4(matrix)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:,:), INTENT(IN) :: matrix
    INTEGER(i4), DIMENSION(:), allocatable  :: diag_i4

    INTEGER(i4) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_i4: array must be squared matrix.'
    if (.not. allocated(diag_i4)) allocate(diag_i4(size(matrix,1)))

    forall(i=1:size(matrix,1)) diag_i4(i) = matrix(i,i)

  END FUNCTION diag_i4

  FUNCTION diag_i8(matrix)

    IMPLICIT NONE

    INTEGER(i8), DIMENSION(:,:), INTENT(IN) :: matrix
    INTEGER(i8), DIMENSION(:), allocatable  :: diag_i8

    INTEGER(i8) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_i8: array must be squared matrix.'
    if (.not. allocated(diag_i8)) allocate(diag_i8(size(matrix,1)))

    forall(i=1:size(matrix,1)) diag_i8(i) = matrix(i,i)

  END FUNCTION diag_i8

  FUNCTION diag_lgt(matrix)

    IMPLICIT NONE

    LOGICAL, DIMENSION(:,:), INTENT(IN) :: matrix
    LOGICAL, DIMENSION(:), allocatable  :: diag_lgt

    INTEGER(i4) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_lgt: array must be squared matrix.'
    if (.not. allocated(diag_lgt)) allocate(diag_lgt(size(matrix,1)))

    forall(i=1:size(matrix,1)) diag_lgt(i) = matrix(i,i)

  END FUNCTION diag_lgt

  ! ------------------------------------------------------------------

  FUNCTION inverse_dp(matrix, condition)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN)  :: matrix
    LOGICAL,  OPTIONAL,       INTENT(IN)  :: condition
    REAL(dp), DIMENSION(:,:), allocatable :: inverse_dp

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
    if (.not. allocated(inverse_dp)) allocate(inverse_dp(nn,nn))
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


  ! FUNCTION inverse_sp(matrix, condition)

  !   IMPLICIT NONE

  !   REAL(sp), DIMENSION(:,:), INTENT(IN)  :: matrix
  !   LOGICAL,  OPTIONAL,       INTENT(IN)  :: condition
  !   REAL(sp), DIMENSION(:,:), allocatable :: inverse_sp

  !   INTEGER(i4) :: i, nn
  !   real(sp),    dimension(:), allocatable :: scale_cols    ! scale matrix for better conditioning
  !   integer(i4), dimension(:), allocatable :: ipiv          ! needed for dgetrf lapack routine
  !   integer(i4)                            :: info          !              "
  !   real(sp),    dimension(:), allocatable :: work          ! needed for dgetri lapack routine
  !   integer(i4)                            :: lwork         !              "
  !   logical :: icondition

  !   external :: sgetrf ! Lapack routine to compute LU factorization of a general matrix
  !   external :: sgetri ! Lapack routine to compute inverse of matrix using the LU factorization computed by DGETRF

  !   nn = size(matrix,1)
  !   if (size(matrix,2) /= nn) stop 'inverse_sp: matrix must be square.'
  !   if (.not. allocated(inverse_sp)) allocate(inverse_sp(nn,nn))
  !   inverse_sp = matrix

  !   if (present(condition)) then
  !      icondition = condition
  !   else
  !      icondition = .true.
  !   endif

  !   if (icondition) then
  !      ! Condition columns
  !      allocate(scale_cols(nn))
  !      scale_cols(:) = maxval(abs(inverse_sp(:,:)),1)
  !      where (scale_cols(:) < tiny(1.0_sp)) scale_cols(:) = 1.0_sp
  !      scale_cols(:) = 1.0_sp / scale_cols(:)
  !      forall(i=1:nn) inverse_sp(:,i) = inverse_sp(:,i) * scale_cols(i)
  !   endif

  !   ! LU factorisation of imatrix
  !   allocate(ipiv(nn))
  !   call sgetrf(nn, nn, inverse_sp, nn, ipiv, info)
  !   if (info /= 0) stop 'inverse_sp: LU factorisation did not work.'
  !   ! Inverse of LU factorisation of imatrix
  !   allocate(work(1))
  !   call sgetri(nn, inverse_sp, nn, ipiv, work, -1, info)
  !   if (info /= 0) stop 'inverse_sp: determination of working matrix dimensions did not work.'
  !   lwork = int(work(1),i4)
  !   deallocate(work)
  !   allocate(work(lwork))
  !   call sgetri(nn, inverse_sp, nn, ipiv, work, lwork, info)
  !   if (info /= 0) stop 'inverse_sp: Inversion did not work.'
  !   deallocate(work)

  !   if (icondition) then
  !      ! rescale result
  !      forall(i=1:nn) inverse_sp(i,:) = inverse_sp(i,:) * scale_cols(i)
  !      deallocate(scale_cols)
  !   endif

  !   deallocate(ipiv)

  ! END FUNCTION inverse_sp

  FUNCTION inverse_sp(matrix, condition)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN)  :: matrix
    LOGICAL,  OPTIONAL,       INTENT(IN)  :: condition
    REAL(sp), DIMENSION(:,:), allocatable :: inverse_sp

    inverse_sp = real(inverse_dp(real(matrix,dp), condition), sp)
    
  END FUNCTION inverse_sp

  ! ------------------------------------------------------------------

  FUNCTION solve_linear_equations_1_dp(lhs, rhs, condition)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: lhs
    REAL(dp), DIMENSION(:),   INTENT(IN) :: rhs
    LOGICAL,  OPTIONAL,       INTENT(IN) :: condition
    REAL(dp), DIMENSION(:), allocatable  :: solve_linear_equations_1_dp

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

    if (.not. allocated(solve_linear_equations_1_dp)) allocate(solve_linear_equations_1_dp(nZeilen))
    solve_linear_equations_1_dp = irhs

    if (icondition) then
       ! rescale result
       solve_linear_equations_1_dp(:) = solve_linear_equations_1_dp(:) * scale_cols(:)

       deallocate(scale_cols, scale_rows)
    endif

    deallocate(ilhs, irhs)
    deallocate(ipiv)

  END FUNCTION solve_linear_equations_1_dp


  ! FUNCTION solve_linear_equations_1_sp(lhs, rhs, condition)

  !   IMPLICIT NONE

  !   REAL(sp), DIMENSION(:,:), INTENT(IN) :: lhs
  !   REAL(sp), DIMENSION(:),   INTENT(IN) :: rhs
  !   LOGICAL,  OPTIONAL,       INTENT(IN) :: condition
  !   REAL(sp), DIMENSION(:), allocatable  :: solve_linear_equations_1_sp

  !   INTEGER(i4) :: ii, nZeilen
  !   real(sp),    dimension(:,:), allocatable :: ilhs          ! internal lhs
  !   real(sp),    dimension(:),   allocatable :: irhs          ! internal rhs
  !   real(sp),    dimension(:),   allocatable :: scale_cols    ! scale matrix for better conditioning
  !   real(sp),    dimension(:),   allocatable :: scale_rows    !              "
  !   integer(i4), dimension(:),   allocatable :: ipiv          ! needed for dgesv lapack routine
  !   integer(i4)                              :: info          !              "
  !   logical :: icondition

  !   external :: sgesv  ! Lapack routine to compute solution of real system of linear equations

  !   nZeilen = size(lhs,1)
  !   if (size(lhs,2) /= nZeilen) stop 'solve_linear_equations_1_sp: left hand side must be squared matrix.'
  !   if (size(rhs,1) /= nZeilen) stop 'solve_linear_equations_1_sp: right hand side must have same size as left hand side.'
  !   ! internal arrays
  !   allocate(ilhs(nZeilen,nZeilen), irhs(nZeilen))
  !   ilhs = lhs
  !   irhs = rhs

  !   if (present(condition)) then
  !      icondition = condition
  !   else
  !      icondition = .true.
  !   endif

  !   if (icondition) then
  !      ! Condition of matrix
  !      allocate(scale_cols(nZeilen), scale_rows(nZeilen))
  !      ! Condition columns
  !      scale_cols(:) = maxval(abs(ilhs(:,:)),1)
  !      where (scale_cols(:) < tiny(1.0_sp)) scale_cols(:) = 1.0_sp
  !      scale_cols(:) = 1.0_sp / scale_cols(:)
  !      forall(ii=1:nZeilen) ilhs(:,ii) = ilhs(:,ii) * scale_cols(ii)
  !      ! Condition rows
  !      scale_rows(:) = maxval(abs(ilhs(:,:)),2)
  !      where (scale_rows(:) < tiny(1.0_sp)) scale_rows(:) = 1.0_sp
  !      scale_rows(:) = 1.0_sp / scale_rows(:)
  !      forall(ii=1:nZeilen) ilhs(ii,:) = ilhs(ii,:) * scale_rows(ii)
  !      irhs(:) = irhs(:) * scale_rows(:)
  !   endif

  !   ! solve linear system of equations
  !   allocate(ipiv(nZeilen))
  !   call sgesv(nZeilen, 1, ilhs, nZeilen, ipiv, irhs, nZeilen, info)
  !   if (info /= 0) stop 'solve_linear_equations_1_sp: Solving of linear system did not work.'

  !   if (.not. allocated(solve_linear_equations_1_sp)) allocate(solve_linear_equations_1_sp(nZeilen))
  !   solve_linear_equations_1_sp = irhs

  !   if (icondition) then
  !      ! rescale result
  !      solve_linear_equations_1_sp(:) = solve_linear_equations_1_sp(:) * scale_cols(:)

  !      deallocate(scale_cols, scale_rows)
  !   endif

  !   deallocate(ilhs, irhs)
  !   deallocate(ipiv)

  ! END FUNCTION solve_linear_equations_1_sp

  FUNCTION solve_linear_equations_1_sp(lhs, rhs, condition)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN) :: lhs
    REAL(sp), DIMENSION(:),   INTENT(IN) :: rhs
    LOGICAL,  OPTIONAL,       INTENT(IN) :: condition
    REAL(sp), DIMENSION(:), allocatable  :: solve_linear_equations_1_sp

    solve_linear_equations_1_sp = real(solve_linear_equations_1_dp(real(lhs,dp), real(rhs,dp), condition), sp)

  END FUNCTION solve_linear_equations_1_sp

  ! ------------------------------------------------------------------

  FUNCTION solve_linear_equations_svd_1_dp(lhs, rhs, condition)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: lhs
    REAL(dp), DIMENSION(:),   INTENT(IN) :: rhs
    LOGICAL,  OPTIONAL,       INTENT(IN) :: condition
    REAL(dp), DIMENSION(:), allocatable  :: solve_linear_equations_svd_1_dp

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

    if (.not. allocated(solve_linear_equations_svd_1_dp)) allocate(solve_linear_equations_svd_1_dp(nZeilen))
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


  ! FUNCTION solve_linear_equations_svd_1_sp(lhs, rhs, condition)

  !   IMPLICIT NONE

  !   REAL(sp), DIMENSION(:,:), INTENT(IN) :: lhs
  !   REAL(sp), DIMENSION(:),   INTENT(IN) :: rhs
  !   LOGICAL,  OPTIONAL,       INTENT(IN) :: condition
  !   REAL(sp), DIMENSION(:), allocatable  :: solve_linear_equations_svd_1_sp

  !   INTEGER(i4) :: ii, nZeilen
  !   real(sp),    dimension(:,:), allocatable :: ilhs          ! internal lhs
  !   real(sp),    dimension(:),   allocatable :: irhs          ! internal rhs
  !   real(sp),    dimension(:),   allocatable :: scale_cols    ! scale matrix for better conditioning
  !   real(sp),    dimension(:),   allocatable :: scale_rows    !              "
  !   real(sp),    dimension(:),   allocatable :: work          ! needed for dgesvd lapack routine
  !   real(sp),    dimension(:),   allocatable :: svdw          !              "
  !   real(sp),    dimension(:,:), allocatable :: svdu          !              "
  !   real(sp),    dimension(:,:), allocatable :: svdv          !              "
  !   integer(i4)                              :: lwork         !              "
  !   integer(i4)                              :: info          !              "
  !   real(sp),    parameter                   :: svdtol = 1.0e-5_sp ! if <svdtol*maxval(svdw), then set svdw=0
  !   logical :: icondition

  !   external :: sgesvd ! Lapack routine to compute singular value decomposition of matrix

  !   nZeilen = size(lhs,1)
  !   if (size(lhs,2) /= nZeilen) stop 'solve_linear_equations_svd_1_sp: left hand side must be squared matrix.'
  !   if (size(rhs,1) /= nZeilen) stop 'solve_linear_equations_svd_1_sp: right hand side must have same size as left hand side.'
  !   ! internal arrays
  !   allocate(ilhs(nZeilen,nZeilen), irhs(nZeilen))
  !   ilhs = lhs
  !   irhs = rhs

  !   if (present(condition)) then
  !      icondition = condition
  !   else
  !      icondition = .true.
  !   endif

  !   if (icondition) then
  !      ! Condition of matrix
  !      allocate(scale_cols(nZeilen), scale_rows(nZeilen))
  !      ! Condition columns
  !      scale_cols(:) = maxval(abs(ilhs(:,:)),1)
  !      where (scale_cols(:) < tiny(1.0_sp)) scale_cols(:) = 1.0_sp
  !      scale_cols(:) = 1.0_sp / scale_cols(:)
  !      forall(ii=1:nZeilen) ilhs(:,ii) = ilhs(:,ii) * scale_cols(ii)
  !      ! Condition rows
  !      scale_rows(:) = maxval(abs(ilhs(:,:)),2)
  !      where (scale_rows(:) < tiny(1.0_sp)) scale_rows(:) = 1.0_sp
  !      scale_rows(:) = 1.0_sp / scale_rows(:)
  !      forall(ii=1:nZeilen) ilhs(ii,:) = ilhs(ii,:) * scale_rows(ii)
  !      irhs(:) = irhs(:) * scale_rows(:)
  !   endif

  !   ! solve linear system of equations
  !   ! Use Lapack singular value decomposition and Numerical Recipes solution
  !   allocate(work(1))
  !   allocate(svdu(nZeilen,nZeilen), svdv(nZeilen,nZeilen), svdw(nZeilen))
  !   call sgesvd('A', 'A', nZeilen, nZeilen, ilhs, nZeilen, svdw, svdu, nZeilen, svdv, nZeilen, work, -1, info)
  !   lwork = int(work(1),i4)
  !   deallocate(work)
  !   allocate(work(lwork))
  !   call sgesvd('A', 'A', nZeilen, nZeilen, ilhs, nZeilen, svdw, svdu, nZeilen, svdv, nZeilen, work, lwork, info)
  !   if (info /= 0) stop 'solve_linear_equations_svd_1_sp: Solving of linear system did not work.'
  !   ! write(*,*) 'Matrix condition ', maxval(abs(svdw))/minval(abs(svdw))
  !   where (svdw < svdtol*maxval(svdw)) svdw = 0.0_sp

  !   if (.not. allocated(solve_linear_equations_svd_1_sp)) allocate(solve_linear_equations_svd_1_sp(nZeilen))
  !   ! svdv is V**T from dgesvd, svdksb wants V
  !   solve_linear_equations_svd_1_sp = svdksb(svdu,svdw,transpose(svdv),irhs)

  !   if (icondition) then
  !      ! rescale result
  !      solve_linear_equations_svd_1_sp(:) = solve_linear_equations_svd_1_sp(:) * scale_cols(:)

  !      deallocate(scale_cols, scale_rows)
  !   endif

  !   deallocate(ilhs, irhs)
  !   deallocate(work, svdu, svdv, svdw)

  ! END FUNCTION solve_linear_equations_svd_1_sp

  FUNCTION solve_linear_equations_svd_1_sp(lhs, rhs, condition)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:), INTENT(IN) :: lhs
    REAL(sp), DIMENSION(:),   INTENT(IN) :: rhs
    LOGICAL,  OPTIONAL,       INTENT(IN) :: condition
    REAL(sp), DIMENSION(:), allocatable  :: solve_linear_equations_svd_1_sp

    solve_linear_equations_svd_1_sp = real(solve_linear_equations_svd_1_dp(real(lhs,dp), real(rhs,dp), condition), sp)

  END FUNCTION solve_linear_equations_svd_1_sp

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

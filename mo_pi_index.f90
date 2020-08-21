!> \file mo_pi_index.f90

!> \brief Parameter sensitivity estimation using Parameter importance (PI) index.

!> \details This module contains routines to determine parameter sensitivities using
!>          the Parameter importance index (PI) based on the eigendecomposition of the sensitivity matrix.\n
!>          Based on pi.py which is part of the JAMS Python library.

!> \authors Juliane Mai
!> \date Jul 2013
MODULE mo_pi_index

  ! Written  Juliane Mai, Jul 2013

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2013 Juliane Mai
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

  ! If you use this routine in your work, you should cite the following reference
  ! Goehler M, J Mai, and M Cuntz (2013)
  !     Use of eigendecomposition in a parameter sensitivity analysis of the Community Land Model,
  !     J Geophys Res 188, 904-921, doi:10.1002/jgrg.20072

  USE mo_kind,   ONLY: i4, sp, dp
  USE mo_mad,    ONLY: mad

  IMPLICIT NONE

  PUBLIC :: pi_index         ! Parameter importance index

  ! ------------------------------------------------------------------

  !     NAME
  !         pi_index

  !     PURPOSE
  !>        \brief The parameter importance index (PI).
  !
  !>        \details Calculates the parameter importance index which is based on the eigendecomposition of \f$ S^T S \f$
  !>                 where S is the sensitivity matrix. The sensitivity matrix S contains nomalized approximates of
  !>                 the first partial derivatives, i.e. the change of the model due to a small change of the parameter.\n
  !>                 Either the sensitivity matrix S or the matrix \f$ M = S^T S \f$ can be given. \n
  !>                 If S is given, an apriori outlier test (MAD: domad) can be applied to the columns of the S matrix and
  !>                 the columnwise number of valid entries of matrix S (counter) can be returned.\n
  !>                 Optionally, the eigenvalues and eigenvectors can be reurned. \n
  !>                 There are different options to normalize the estimated sensitivity index. \n
  !>                 Optionally, the B index (b_index) which is a sensitivity index neglecting covariations between
  !>                 parameters can be returned.
  !>                 The same normalization (norm) will be applied to both the B index and the PI index.
  !>
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         \param[out] "real(sp/dp)           :: pi(:)"       Parameter importance index per parameter
  !
  !     INTENT(IN), OPTIONAL
  !>       \param[in] "real(sp/dp), optional   :: s(:,:)"      Sensitivity matrix \n
  !>                                                              dim_1 = number of parameter sets \n
  !>                                                              dim_2 = number of parameters
  !>       \param[in] "real(sp/dp), optional   :: m(:,:)"      \f$ M = S^T S \f$ \n
  !>                                                              dim_1 = dim_2 = number of parameters \n
  !>       \param[in] "integer(i4), optional   :: norm"        Normalization of index \n
  !>                                                              0 - no normalization (default) \n
  !>                                                              1 - normalized such that all PI sum up to 1 \n
  !>                                                              2 - normalized by sum of eigenvalues \n
  !>                                                              3 - normalized by sum of eigenvalues and
  !>                                                                   subsequently by sum of PI indexes
  !>       \param[in] "logical,     optional   :: domad"       If MAD outlier test should be performed on each
  !>                                                              column of matrix S. \n
  !>                                                              Only applicable if s is present. \n
  !>                                                              Default: .false.
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !>       \param[out] "integer(i4), allocatable, optional :: counter(:)"    Returns number of valid entries in
  !>                                                                            matrix S per column. \n
  !>                                                                            Only applicable if s is present. \n
  !>       \param[out] "real(sp/dp), optional              :: evalues(:)"    Eigenvalues of \f$ M = S^T S \f$ \n
  !>       \param[out] "real(sp/dp), optional              :: evectors(:,:)" Eigenvectors of \f$ M = S^T S \f$ \n
  !>       \param[out] "real(sp/dp), optional              :: b_index(:)"    B index, i.e. diagonal elements of M \n
  !>                                                                         normalized using norm \n

  !
  !     RETURN
  !        None
  !
  !     RESTRICTIONS
  !>       \note Needs Lapack. \n
  !>             Some options are only applicable if matrix S is present. \n
  !
  !     EXAMPLE
  !         s(:,1) = (/ 1.0, 2.0, 3.0, -999., 5.0, 6.0 /)
  !         s(:,2) = (/ 1.2, 2.5, 3.4, -999., 5.2, 6.8 /)
  !         call pi_index(pi, s=s)
  !         -> see also example in test directory

  !     LITERATURE
  !         Vajda, S., Valko, P., & Turanyi, T. (2004).
  !              Principal component analysis of kinetic models.
  !              International Journal of Chemical Kinetics, 17(1), 55-81.
  !         Goehler, M., Mai, J., & Cuntz, M. (2013).
  !              Use of eigendecomposition in a parameter sensitivity analysis of the Community Land Model.
  !              Journal of Geophysical Research: Biogeosciences, 118(2), 904-921. doi:10.1002/jgrg.20072

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Jul 2013
  !         Modified,

  !     The code is based on pi.py by Matthias Cuntz. pi.py is part of the JAMS Python library.

  INTERFACE pi_index
     MODULE PROCEDURE pi_index_dp, pi_index_sp
  END INTERFACE pi_index

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

#ifndef __PYTHON__
  subroutine pi_index_dp(pi_index_out, s, m, norm, domad, counter, evalues, evectors, b_index)
#else
  subroutine pi_index_dp(pi_index_out, s, m, norm, domad)
#endif

    real(dp), dimension(:),                   intent(out)           :: pi_index_out     ! Parameter importance index
    !                                         ! either S or M has to be given
    real(dp), dimension(:,:),                 intent(in),  optional :: s            ! sensitivity matrix S
    real(dp), dimension(:,:),                 intent(in),  optional :: m            ! M = S^T.S
    integer(i4),                              intent(in),  optional :: norm         ! normalization of PI index
    !                                                                               ! 0 - no normalization (DEFAULT)
    !                                                                               ! 1 - normalized such that all PI sum up to 1
    !                                                                               ! 2 - normalized by sum of eigenvalues
    !                                                                               ! 3 - normalized by sum of eigenvalues and
    !                                                                               !     subsequently by sum of PI indexes
    logical,                                  intent(in),  optional :: domad        ! .true. : prior mad oulier test on S matrix
    !                                                                               ! restriction: S matrix must to be given
    !                                                                               ! DEFAULT: .false.
#ifndef __PYTHON__
    integer(i4), dimension(:),   allocatable, intent(out), optional :: counter      ! number of valid entries in column of S matrix
    real(dp),    dimension(:),   allocatable, intent(out), optional :: evalues      ! eigenvalues of M
    real(dp),    dimension(:,:), allocatable, intent(out), optional :: evectors     ! eigenvectors of M
    real(dp),    dimension(:),   allocatable, intent(out), optional :: b_index      ! B index = sensitivity index neglecting
    !                                                                               ! covariations, i.e. diagonal elements of
    !                                                                               ! matrix M normalized using norm
#endif

    ! local variables
    integer(i4)                              :: ii, jj, ii_count, jj_count
    integer(i4)                              :: ipar
    integer(i4)                              :: npar        ! number of parameters
    logical                                  :: my_mad      ! if mad-test should performed or not
    integer(i4)                              :: my_norm
    logical,     dimension(:,:), allocatable :: mask
    logical,     dimension(:),   allocatable :: maskpara
    real(dp),    dimension(:,:), allocatable :: my_s
    real(dp),    dimension(:,:), allocatable :: my_m
    integer(i4), dimension(:,:), allocatable :: imask
    integer(i4), dimension(:,:), allocatable :: ntrue
    integer(i4), dimension(:),   allocatable :: my_counter  ! number of valid entries in ith column of matrix S
    real(dp),    dimension(:),   allocatable :: my_evalues  ! eigenvalues of M
    real(dp),    dimension(:,:), allocatable :: my_evectors ! eigenvectors of M
    integer(i4)                              :: lwork, info
    real(dp),    dimension(:),   allocatable :: work
    logical,     dimension(:),   allocatable :: tmpmad

    external :: dsyev

    if ( (.not. present(s) .and. .not. present(m)) .or.  (present(s) .and. present(m)) ) then
       stop 'pi_index_dp: either s or m has to be given'
    end if

    if (present(domad)) then
       my_mad = domad
    else
       my_mad = .false.
    end if

    npar = size(pi_index_out)
    allocate( maskpara(npar) )
    maskpara = .true.

    if (present(s)) then
       if ( npar .ne. size(s,2) ) stop 'pi_index_dp: size of s is not matching'
       allocate(my_counter(npar))
       my_counter = size(s,1)

       if (my_mad) then
          allocate( my_s(size(s,1),npar) )
          allocate( my_m(npar,npar) )
          allocate( imask(size(s,1),npar) )
          allocate( mask(size(s,1),npar) )
          allocate( ntrue(npar,npar) )
          allocate( tmpmad(size(s,1)) )

          mask     = .true.
          ! mask = .false. if entry in S matrix is zero (i.e. no response of model to change)
          where (abs(s) .lt. tiny(1.0_dp))
             mask = .false.
          end where
          ! maskpara = true,  if there are not only zeros in column of matrix S
          ! maskpara = false, if there are only zeros in column of matrix S
          maskpara = (count(mask,dim=1) .ne. 0_i4)
          do ipar=1, npar
             if (maskpara(ipar)) then
                tmpmad = mad(s(:,ipar), mask=mask(:,ipar), z=15.0_dp)
                mask(:,ipar)   = mask(:,ipar) .and. tmpmad
                my_counter(ipar) = count(mask(:,ipar))
                maskpara(ipar) = ( my_counter(ipar) .ne. 0_i4)
             else
                my_counter(ipar) = 0_i4
             end if
          end do

          imask = 1_i4
          my_s  = s
          ! set outliers in matrix S to zero
          where (.not. mask)
             imask = 0_i4
             my_s  = 0.0_dp
          end where
          ntrue = matmul(transpose(imask),imask)

          ! matrix M
          my_m = matmul(transpose(my_s),my_s)

          ! delete colums and rows in matrix M where maskpara=false
          allocate( my_evalues(count(maskpara)) )
          allocate (my_evectors(count(maskpara),count(maskpara)))
          ii_count = 0_i4
          do ii=1,npar
             if (maskpara(ii)) ii_count = ii_count + 1
             jj_count = 0_i4
             do jj=1,npar
                if (maskpara(jj)) jj_count = jj_count + 1
                if (maskpara(ii) .and. maskpara(jj)) then
                   my_evectors(ii_count, jj_count) = my_m(ii,jj)
                   ! and scale entries
                   my_evectors(ii_count, jj_count) = my_evectors(ii_count, jj_count) * real(size(my_s,1)) / real(ntrue(ii,jj),dp)
                end if
             end do
          end do

          deallocate(my_s)
          deallocate(my_m)
          deallocate(imask)
          deallocate(mask)
          deallocate(ntrue)
       else
          allocate( my_evalues(count(maskpara)) )
          allocate( my_evectors(count(maskpara),count(maskpara)) )
          my_evectors = matmul(transpose(s),s)
       end if
    else
       if ( (npar .ne. size(m,1)) .or. (npar .ne. size(m,2)) ) stop 'pi_index_dp: size of m is not matching'
#ifndef __PYTHON__
       if ( present(counter) ) stop 'pi_index_dp: argument counter only applicable if matrix S is given'
#endif
       allocate( my_evalues(count(maskpara)) )
       allocate( my_evectors(count(maskpara),count(maskpara)) )
       my_evectors = m
    end if

    if (present(norm)) then
       my_norm = norm
    else
       my_norm =  0
    end if

#ifndef __PYTHON__
    if (present(b_index)) then
       allocate(b_index(npar))
       forall(ii=1:npar) b_index(ii) = my_evectors(ii,ii)
       ! Normalize B index
       select case(my_norm)
       case(0_i4)  ! No normalization
          b_index = b_index
       case(1_i4)  ! Normalized such that B sum up to one
          b_index = b_index / sum(b_index)
       case(2_i4)  ! Normalized by sum of eigenvalues
          b_index = b_index / sum(my_evalues)
       case(3_i4)  ! Normalized by sum of eigenvalues and subsequent by sum of B
          b_index = b_index / sum(my_evalues)
          b_index = b_index / sum(b_index)
       case default
          stop 'pi_index_dp: This normalization method is not implemented.'
       end select
    end if
#endif

    ! Eigenvalues of Covariance Matrix
    ! (1) Query for optimal workspace: lwork=-1
    allocate(work(1))
    lwork = -1_i4
    call dsyev( 'V', 'L', count(maskpara), my_evectors, count(maskpara), my_evalues, work, lwork, info )
    lwork =int( work(1) )
    deallocate(work)

    ! (2) allocate "work" with proper length
    allocate(work(lwork))

    ! (3) calculate eigenvalues
    call dsyev( 'V', 'L', count(maskpara), my_evectors, count(maskpara), my_evalues, work, lwork, info )

    ! (4) Info statement
    if ( info .gt. 0_i4 ) then
       stop 'pi_index_dp: The algorithm failed to compute eigenvalues.'
    end if
    deallocate(work)

    ! Make Eigenvalues all non-negative
    my_evalues(:) = Max(0.0_dp, my_evalues(:))

    ! Calculate PI index
    ii_count = 0_i4
    do ipar=1, npar
       if (maskpara(ipar)) then
          ii_count = ii_count + 1
          pi_index_out(ipar) = dot_product(abs(my_evectors(ii_count,:)), my_evalues)
       else
          pi_index_out(ipar) = 0.0_dp
       end if
    end do

    ! Normalize PI index
    select case(my_norm)
    case(0_i4)  ! No normalization
       pi_index_out = pi_index_out
    case(1_i4)  ! Normalized such that PI sum up to one
       pi_index_out = pi_index_out / sum(pi_index_out)
    case(2_i4)  ! Normalized by sum of eigenvalues
       pi_index_out = pi_index_out / sum(my_evalues)
    case(3_i4)  ! Normalized by sum of eigenvalues and subsequent by sum of PI
       pi_index_out = pi_index_out / sum(my_evalues)
       pi_index_out = pi_index_out / sum(pi_index_out)
    case default
       stop 'pi_index_dp: This normalization method is not implemented.'
    end select

#ifndef __PYTHON__
    if (present(counter)) then
       allocate(counter(npar))
       counter = my_counter
    end if

    if (present(evalues)) then
       allocate(evalues(npar))
       evalues = my_evalues
    end if

    if (present(evectors)) then
       allocate(evectors(npar,npar))
       evectors = my_evectors
    end if
#endif

    deallocate(my_counter)
    deallocate(my_evalues)
    deallocate(my_evectors)

  end subroutine pi_index_dp

#ifndef __PYTHON__
  subroutine pi_index_sp(pi_index_out, s, m, norm, domad, counter, evalues, evectors)
#else
  subroutine pi_index_sp(pi_index_out, s, m, norm, domad)
#endif

    real(sp), dimension(:),                   intent(out)           :: pi_index_out     ! Parameter importance index
    ! either S or M has to be given
    real(sp), dimension(:,:),                 intent(in),  optional :: s            ! sensitivity matrix S
    real(sp), dimension(:,:),                 intent(in),  optional :: m            ! M = S^T.S
    integer(i4),                              intent(in),  optional :: norm         ! normalized PI index
    !                                                                               ! 0 - no normalization (DEFAULT)
    !                                                                               ! 1 - normalized such that all PI sum up to 1
    !                                                                               ! 2 - normalized by sum of eigenvalues
    !                                                                               ! 3 - noramlized by sum of eigenvalues and
    !                                                                               !     subsequently by sum of PI indexes
    logical,                                  intent(in),  optional :: domad        ! .true. : prior mad oulier test on S matrix
    !                                                                               ! restriction: S matrix must to be given
    !                                                                               ! DEFAULT: .false.
#ifndef __PYTHON__
    integer(i4), dimension(:),   allocatable, intent(out), optional :: counter      ! number of valid entries in column of S matrix
    real(sp),    dimension(:),   allocatable, intent(out), optional :: evalues      ! eigenvalues of M
    real(sp),    dimension(:,:), allocatable, intent(out), optional :: evectors     ! eigenvectors of M
#endif

    ! local variables
    integer(i4)                              :: ii, jj, ii_count, jj_count
    integer(i4)                              :: ipar
    integer(i4)                              :: npar        ! number of parameters
    logical                                  :: my_mad      ! if mad-test should performed or not
    integer(i4)                              :: my_norm
    logical,     dimension(:,:), allocatable :: mask
    logical,     dimension(:),   allocatable :: maskpara
    real(sp),    dimension(:,:), allocatable :: my_s
    real(sp),    dimension(:,:), allocatable :: my_m
    integer(i4), dimension(:,:), allocatable :: imask
    integer(i4), dimension(:,:), allocatable :: ntrue
    integer(i4), dimension(:),   allocatable :: my_counter  ! number of valid entries in ith column of matrix S
    real(sp),    dimension(:),   allocatable :: my_evalues  ! eigenvalues of M
    real(sp),    dimension(:,:), allocatable :: my_evectors ! eigenvectors of M
    integer(i4)                              :: lwork, info
    real(sp),    dimension(:),   allocatable :: work
    logical,     dimension(:),   allocatable :: tmpmad

    external :: ssyev

    if ( (.not. present(s) .and. .not. present(m)) .or.  (present(s) .and. present(m)) ) then
       stop 'pi_index_sp: either s or m has to be given'
    end if

    if (present(domad)) then
       my_mad = domad
    else
       my_mad = .false.
    end if

    npar = size(pi_index_out)
    allocate( maskpara(npar) )
    maskpara = .true.

    if (present(s)) then
       if ( npar .ne. size(s,2) ) stop 'pi_index_sp: size of s is not matching'
       allocate(my_counter(npar))
       my_counter = size(s,1)

       if (my_mad) then
          allocate( my_s(size(s,1),npar) )
          allocate( my_m(npar,npar) )
          allocate( imask(size(s,1),npar) )
          allocate( mask(size(s,1),npar) )
          allocate( ntrue(npar,npar) )
          allocate( tmpmad(size(s,1)) )

          mask     = .true.
          ! mask = .false. if entry in S matrix is zero (i.e. no response of model to change)
          where (abs(s) .lt. tiny(1.0_sp))
             mask = .false.
          end where
          ! maskpara = true,  if there are not only zeros in column of matrix S
          ! maskpara = false, if there are only zeros in column of matrix S
          maskpara = (count(mask,dim=1) .ne. 0_i4)
          do ipar=1, npar
             if (maskpara(ipar)) then
                tmpmad = mad(s(:,ipar), mask=mask(:,ipar), z=15.0_sp)
                mask(:,ipar)   = mask(:,ipar) .and. tmpmad
                my_counter(ipar) = count(mask(:,ipar))
                maskpara(ipar) = ( my_counter(ipar) .ne. 0_i4)
             else
                my_counter(ipar) = 0_i4
             end if
          end do

          imask = 1_i4
          my_s  = s
          ! set outliers in matrix S to zero
          where (.not. mask)
             imask = 0_i4
             my_s  = 0.0_sp
          end where
          ntrue = matmul(transpose(imask),imask)

          ! matrix M
          my_m = matmul(transpose(my_s),my_s)

          ! delete colums and rows in matrix M where maskpara=false
          allocate( my_evalues(count(maskpara)) )
          allocate (my_evectors(count(maskpara),count(maskpara)))
          ii_count = 0_i4
          do ii=1,npar
             if (maskpara(ii)) ii_count = ii_count + 1
             jj_count = 0_i4
             do jj=1,npar
                if (maskpara(jj)) jj_count = jj_count + 1
                if (maskpara(ii) .and. maskpara(jj)) then
                   my_evectors(ii_count, jj_count) = my_m(ii,jj)
                   ! and scale entries
                   my_evectors(ii_count, jj_count) = my_evectors(ii_count, jj_count) * real(size(my_s,1)) / real(ntrue(ii,jj),sp)
                end if
             end do
          end do

          deallocate(my_s)
          deallocate(my_m)
          deallocate(imask)
          deallocate(mask)
          deallocate(ntrue)
       else
          allocate( my_evalues(count(maskpara)) )
          allocate( my_evectors(count(maskpara),count(maskpara)) )
          my_evectors = matmul(transpose(s),s)
       end if
    else
       if ( (npar .ne. size(m,1)) .or. (npar .ne. size(m,2)) ) stop 'pi_index_sp: size of m is not matching'
#ifndef __PYTHON__
       if ( present(counter) ) stop 'pi_index_sp: argument counter only applicable if matrix S is given'
#endif
       allocate( my_evalues(count(maskpara)) )
       allocate( my_evectors(count(maskpara),count(maskpara)) )
       my_evectors = m
    end if

    if (present(norm)) then
       my_norm = norm
    else
       my_norm =  0
    end if

    ! Eigenvalues of Covariance Matrix
    ! (1) Query for optimal workspace: lwork=-1
    allocate(work(1))
    lwork = -1_i4
    call ssyev( 'Vectors', 'Lower', count(maskpara), my_evectors, count(maskpara), my_evalues, work, lwork, info )
    lwork =int( work(1) )
    deallocate(work)

    ! (2) allocate "work" with proper length
    allocate(work(lwork))

    ! (3) calculate eigenvalues
    call ssyev( 'Vectors', 'Lower', count(maskpara), my_evectors, count(maskpara), my_evalues, work, lwork, info )

    ! (4) Info statement
    if ( info .gt. 0_i4 ) then
       stop 'pi_index_sp: The algorithm failed to compute eigenvalues.'
    end if
    deallocate(work)

    ! Make Eigenvalues all non-negative
    my_evalues(:) = Max(0.0_sp, my_evalues(:))

    ! Calculate PI index
    ii_count = 0_i4
    do ipar=1, npar
       if (maskpara(ipar)) then
          ii_count = ii_count + 1
          pi_index_out(ipar) = dot_product(abs(my_evectors(ii_count,:)), my_evalues)
       else
          pi_index_out(ipar) = 0.0_sp
       end if
    end do

    ! Normalize PI index
    select case(my_norm)
    case(0_i4)  ! No normalization
       pi_index_out = pi_index_out
    case(1_i4)  ! Normalized such that PI sum up to one
       pi_index_out = pi_index_out / sum(pi_index_out)
    case(2_i4)  ! Normalized by sum of eigenvalues
       pi_index_out = pi_index_out / sum(my_evalues)
    case(3_i4)  ! Normalized by sum of eigenvalues and subsequent by sum of PI
       pi_index_out = pi_index_out / sum(my_evalues)
       pi_index_out = pi_index_out / sum(pi_index_out)
    case default
       stop 'pi_index_sp: This normalization method is not implemented.'
    end select

#ifndef __PYTHON__
    if (present(counter)) then
       allocate(counter(npar))
       counter = my_counter
    end if

    if (present(evalues)) then
       allocate(evalues(npar))
       evalues = my_evalues
    end if

    if (present(evectors)) then
       allocate(evectors(npar,npar))
       evectors = my_evectors
    end if
#endif

    deallocate(my_counter)
    deallocate(my_evalues)
    deallocate(my_evectors)

  end subroutine pi_index_sp

END MODULE mo_pi_index

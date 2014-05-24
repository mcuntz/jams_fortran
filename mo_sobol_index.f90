!> \file mo_sobol_index.f90

!> \brief Parameter sensitivity estimation using Sobol index (main and total effect).

!> \details This module contains routines to determine parameter sensitivities using 
!>          the Sobol index (SI) based on the variance of model output changing only one parameter at a time and 
!>          all parameter except one at a time. \n
!>          Based on sobol_index.py which is part of the Python CHS library.

!> \authors Juliane Mai
!> \date Jul 2013

MODULE mo_sobol_index

  ! Written  Juliane Mai,    Jul 2013
  ! Modified 

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! It is NOT released under the GNU Lesser General Public License, yet.

  ! If you use this routine, please contact Juliane Mai.

  ! Copyright 2013 Juliane Mai

  USE mo_kind,   ONLY: i4, sp, dp
  USE mo_moment, ONLY: variance

  IMPLICIT NONE

  PUBLIC :: sobol_index      ! Sobol index

  ! ------------------------------------------------------------------

  !     NAME
  !         sobol_index

  !     PURPOSE
  !>        \brief The Sobol index (SI).
  !
  !>        \details Calculates the Sobol index which is based on the change of variance of 
  !>                 model response due to parameter changes.\n
  !>                 Here, the main effect (Si) and the total effect (STi) are implemented. \n
  !>                 The model output has to be provided for:\n
  !>                   - a parameter set A = f(A)\n
  !>                   - a parameter set B which is independent to A  = f(B)\n
  !>                 For method #3, #4, #5, #6, #7, #8, #9=Default: \n
  !>                   - n parameter sets Ci where Ci equals parameter set A except parameter i 
  !>                     which is the value of parameter set B = f(A_B)\n
  !>                 For method #1, #2, #6: \n
  !>                   - n parameter sets Ci where Ci equals parameter set B except parameter i 
  !>                     which is the value of parameter set A = f(B_A)\n
  !>                 This has to be repeated several times. \n
  !>                 The model response can be either a scalar or a vector (e.g. a timeseries). 
  !>                 If the model response is a vector, the Si and STi are calculated for each time point 
  !>                 independently such that the returned Si and STi are also vectors. \n
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp)      :: ya(:)/ya(:,:)"       Model output for parameter sets A \n
  !>                                                               dim_1 = number of parameter sets \n
  !>                                                               dim_2 = number of model outputs, e.g. time steps
  !>        \param[in] "real(sp/dp)      :: yb(:)/yb(:,:)"       Model output for parameter sets B \n
  !>                                                               dim_1 = number of parameter sets \n
  !>                                                               dim_2 = number of model outputs, e.g. time steps
  !>        \param[in] "real(sp/dp)      :: yc(:,:)/yc(:,:,:)"   Model output for parameter sets C(i), i=1,npara \n
  !>                                                               dim_1 = number of parameter sets \n
  !>                                                               dim_2 = number of parameters \n
  !>                                                               dim_3 = number of model outputs, e.g. time steps
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(sp/dp)      :: si(:)/si(:,:)"   Sobol index - Main effect \n
  !>                                                             dim_1 = number of parameters \n
  !>                                                             dim_2 = number of model outputs, e.g. time steps
  !>        \param[out] "real(sp/dp)      :: sti(:)/sti(:,:)" Sobol index - Total effect \n
  !>                                                             dim_1 = number of parameters \n
  !>                                                             dim_2 = number of model outputs, e.g. time steps
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(sp/dp)      :: smean(:,:)"      Average index per parameter \n
  !>                                                             dim_1 = number of parameters \n
  !>                                                             dim_2 = 2, i.e. SI and STI
  !>        \param[out] "real(sp/dp)      :: wmean(:,:)"      Variance weighted average index per parameter \n
  !>                                                             dim_1 = number of parameters \n
  !>                                                             dim_2 = 2, i.e. SI and STI
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !>       \note Input values must be floating points.\n
  !>             Average indexes are only calculated when ya has more than one dimension, i.e. modeloutput is not scalar. 
  !>             Else average indexes are equal si and sti.
  !
  !     EXAMPLE
  !         ya      = (/ 1., 2, 3., -999., 5., 6. /)
  !         yb      = (/ 2., 3, 1., 7., 8., 9. /)
  !         yc(:,1) = (/ 2., 3, 1., 7., 8., 9. /)
  !         yc(:,2) = (/ 3., 3, 3., 3., 3., 3. /)
  !         call sobol_index(ya, yb, yc, si, sti)
  !         -> see also example in test directory

  !     LITERATURE
  !         Homma, T., & Saltelli, A. (1996). 
  !                Importance measures in global sensitivity analysis of nonlinear models. 
  !                Reliability Engineering & System Safety, 52(1), 1–17.
  !         Saltelli, A. (2002). 
  !                Making best use of model evaluations to compute sensitivity indices. 
  !                Computer Physics Communications, 145(2), 280–297. doi:10.1016/S0010-4655(02)00280-1
  !         Saltelli, A., Ratto, M., Andres, T., Campolongo, F., Cariboni, J., Gatelli, D., et al. (2008).
  !                Global Sensitivity Analysis. Wiley-Interscience.

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Jul 2013
  !         Modified, 

  !     The code is based on sobol_index.py by Matthias Cuntz. Sobol_index.py is part of the Python CHS library.

  INTERFACE sobol_index
     MODULE PROCEDURE sobol_index_0d_dp, sobol_index_0d_sp, &
          sobol_index_1d_dp, sobol_index_1d_sp
  END INTERFACE sobol_index

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  subroutine sobol_index_0d_dp(ya, yb, yc, si, sti, method, smean, wmean)

    real(dp), dimension(:),                       intent(in)  :: ya     ! Output running model with parameter sets A
    real(dp), dimension(:),                       intent(in)  :: yb     ! Output running model with parameter sets B
    real(dp), dimension(:,:),                     intent(in)  :: yc     ! Output running model with parameter sets C
    !                                                                   ! dim1 = number of parameter sets 
    !                                                                   !      = size(ya) = size(yb)
    !                                                                   ! dim2 = number of parameters
    real(dp), dimension(size(yc,2)),              intent(out) :: si     ! Sobol index (main effect)
    real(dp), dimension(size(yc,2)),              intent(out) :: sti    ! Sobol index (total effect)
    integer(i4),                        optional, intent(in)  :: method ! (1) 'Saltelli2008'  needs Ci = B_A
                                                                        ! (2) 'Homma1996'     needs Ci = B_A
                                                                        ! (3) 'Sobol2007'     needs Ci = A_B and B_A --> NA
                                                                        ! (4) 'Saltelli2010'  needs Ci = A_B
                                                                        ! (5) 'Jansen1999'    needs Ci = A_B
                                                                        ! (6) 'Mai2012'       needs Ci = B_A
                                                                        ! (7) 'Mai2013'       needs Ci = A_B
                                                                        ! (8) 'Mai2014'       needs Ci = A_B
                                                                        ! (9) 'Mai1999'       needs Ci = A_B  ----> Default
    real(dp), dimension(size(yc,2), 2), optional, intent(out) :: smean  ! Mean Sobol index (main and total effect)
    !                                                                   ! dim1 = number of parameters
    !                                                                   ! dim2 = 2 (i.e. SI and STI)
    !                                                                   ! == (/ si, sti /) in 0d version !!
    real(dp), dimension(size(yc,2), 2), optional, intent(out) :: wmean  ! Variance weighted mean Sobol index 
    !                                                                   ! (main and total effect)
    !                                                                   ! dim1 = number of parameters
    !                                                                   ! dim2 = 2 (i.e. SI and STI)
    !                                                                   ! == (/ si, sti /) in 0d version !!

    ! local variables
    integer(i4)                          :: ii
    integer(i4)                          :: nsets, npara
    integer(i4)                          :: meth
    real(dp), dimension(:), allocatable  :: yab
    real(dp)                             :: var_a, var_b, var_ab
    real(dp)                             :: mean_a, mean_b, mean_ab

    if (present(method)) then
       meth = method
    else
       meth = 9
    end if

    nsets = size(yc,1)
    npara = size(yc,2)

    select case(meth)
       case(1)
          ! 'Saltelli2008' - The formulation presented in 'The Primer'. (yc=f(B_A))
          !                  Si  = (1/n*sum_j(f(A)_j*f(B_A^i)_j) - mean(f(A))^2)/var(f(A))
          !                  STi = (var(f(A))-(1/n*sum_j(f(B)_j*f(B_A^i)_j) - mean(f(A))^2))/var(f(A))
          var_a  = variance(ya)
          mean_a = sum( ya(:) ) / real(nsets,dp)
          if ( var_a .gt. tiny(1.0_dp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  =          ( dot_product( ya(:) , (yc(:,ii)) ) / real(nsets,dp) - mean_a**2 ) / var_a
                sti(ii) = 1.0_dp - ( dot_product( yb(:) , (yc(:,ii)) ) / real(nsets,dp) - mean_a**2 ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_dp
             sti(:) = 0.0_dp
          end if
       case(2)
          ! 'Homma1996'    - Si takes mixed A, B term for squared mean in nominator. (yc=f(B_A))
          !                  Si  = 1/n*sum_j(f(A)_j*(f(B_A^i)_j - f(B)_j))/var(f(A))
          !                  STi as Saltelli2008
          var_a  = variance(ya)
          mean_a = sum( ya(:) ) / real(nsets,dp)
          if ( var_a .gt. tiny(1.0_dp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  =          ( dot_product( ya(:) , (yc(:,ii)-yb(:)) ) / real(nsets,dp) ) / var_a
                sti(ii) = 1.0_dp - ( dot_product( yb(:) , (yc(:,ii)) ) / real(nsets,dp) - mean_a**2 ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_dp
             sti(:) = 0.0_dp
          end if
       case(3)
          stop 'mo_sobol_index: Sobol2007 would need f(A_B) and f(B_A). It is thus not implemented here.'
       case(4)
          ! 'Saltelli2010' - Si takes B, A_B instead of A, B_A (yc=f(A_B))
          !                  Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var(f(A))
          !                  STi = 1/n*sum_j(f(A)_j*(f(A)_j-f(A_B^i)_j))/var(f(A))
          var_a  = variance(ya)
          if ( var_a .gt. tiny(1.0_dp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  = ( dot_product( yb(:) , (yc(:,ii)-ya(:)) ) / real(nsets,dp) ) / var_a
                sti(ii) = ( dot_product( ya(:) , (ya(:)-yc(:,ii)) ) / real(nsets,dp) ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_dp
             sti(:) = 0.0_dp
          end if
       case(5)
          ! 'Jansen1999'   - Calculate Si and STi by expectation(variance) instead of variance(expectation) (yc=f(A_B))
          !                  Si  = (var(f(A)) - 1/2n*sum_j(f(B)_j - f(A_B^i)_j)^2)/var(f(A))
          !                  STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
          var_a  = variance(ya)
          if ( var_a .gt. tiny(1.0_dp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  = 1.0_dp - ( sum( (yb(:)-yc(:,ii))**2 ) / real(2*nsets,dp) ) / var_a
                sti(ii) =          ( sum( (ya(:)-yc(:,ii))**2 ) / real(2*nsets,dp) ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_dp
             sti(:) = 0.0_dp
          end if
       case(6)
          ! 'Mai2012'      - Si of Homma 1996 but denominator is variance of A and B (yc=f(B_A)) and
          !                  STi of Homma 1996 but mean and variance of f(A) replaced by f(B)
          !                  Si  = 1/n*sum_j(f(A)_j*(f(B_A^i)_j - f(B)_j))/var([f(A),f(B)])
          !                  STi = (var(f(A))-(1/n*sum_j(f(B)_j*f(B_A^i)_j) - mean(f(B)^2)))/var(f(B))
          allocate(yab(2*nsets))
          yab(      1:  nsets) = ya
          yab(nsets+1:2*nsets) = yb
          var_ab = variance(yab)
          var_b  = variance(yb)

          mean_ab = dot_product(ya,yb) / real(nsets,dp)
          mean_b  = sum(yb) / real(nsets,dp)

          if ( var_ab .gt. tiny(1.0_dp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  =          ( dot_product(ya(:),yc(:,ii)) / real(nsets,dp) - mean_ab)   / var_ab
                sti(ii) = 1.0_dp - ( dot_product(yb(:),yc(:,ii)) / real(nsets,dp) - mean_b**2) / var_b
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_dp
             sti(:) = 0.0_dp
          end if
          deallocate(yab)
       case(7)
          ! 'Mai2013'      - Si of Saltelli2010 but denominator is variance of A and B (yc=f(A_B))
          !                  Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var([f(A),f(B)])
          !                  STi = 1/n*sum_j(f(A)_j*(f(A)_j)-f(A_B^i)_j)/var(f(A))
          allocate(yab(2*nsets))
          yab(      1:  nsets) = ya
          yab(nsets+1:2*nsets) = yb
          var_ab = variance(yab)
          var_a  = variance(ya)

          if ( var_ab .gt. tiny(1.0_dp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  = ( dot_product( yb(:) , (yc(:,ii)-ya(:)) ) / real(nsets,dp) ) / var_ab
                sti(ii) = ( dot_product( ya(:) , (ya(:)-yc(:,ii)) ) / real(nsets,dp) ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_dp
             sti(:) = 0.0_dp
          end if
          deallocate(yab)
       case(8)
          ! 'Mai2014'      - SI of Jansen1999  but denominator is variance of A and B (yc=f(A_B))
          !                  Si  = (var([f(A),f(B)]) - 1/2n*sum_j(f(B)_j - f(A_B^i)_j)^2)/var([f(A),f(B)])
          !                  STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
          allocate(yab(2*nsets))
          yab(      1:  nsets) = ya
          yab(nsets+1:2*nsets) = yb
          var_ab = variance(yab)
          var_a  = variance(ya)

          if ( var_ab .gt. tiny(1.0_dp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  = 1.0_dp - ( sum( (yb(:)-yc(:,ii))**2 ) / real(2*nsets,dp) ) / var_ab
                sti(ii) =          ( sum( (ya(:)-yc(:,ii))**2 ) / real(2*nsets,dp) ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_dp
             sti(:) = 0.0_dp
          end if
          deallocate(yab)
       case(9)
          ! 'Mai1999'      - SI of Mai2013 (Saltelli2010 with var([f(A),f(B)])) and STi of Jansen1999 (yc=f(A_B))
          !                  Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var([f(A),f(B)])
          !                  STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
          allocate(yab(2*nsets))
          yab(      1:  nsets) = ya
          yab(nsets+1:2*nsets) = yb
          var_ab = variance(yab)
          var_a  = variance(ya)

          if ( var_ab .gt. tiny(1.0_dp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  = ( dot_product( yb(:) , (yc(:,ii)-ya(:)) ) / real(nsets,dp) ) / var_ab
                sti(ii) = ( sum( (ya(:)-yc(:,ii))**2 ) / real(2*nsets,dp) ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_dp
             sti(:) = 0.0_dp
          end if
          deallocate(yab)
       case default
          stop 'mo_sobol_index: This method is not implemented!'
       end select

   if ( present(smean) ) then
       smean(:,1) = si(:)
       smean(:,2) = sti(:)
    end if

    if ( present(wmean) ) then
       wmean(:,1) = si(:)
       wmean(:,2) = sti(:)
    end if

  end subroutine sobol_index_0d_dp

  subroutine sobol_index_0d_sp(ya, yb, yc, si, sti, method, smean, wmean)

    real(sp), dimension(:),                       intent(in)  :: ya     ! Output running model with parameter sets A
    real(sp), dimension(:),                       intent(in)  :: yb     ! Output running model with parameter sets B
    real(sp), dimension(:,:),                     intent(in)  :: yc     ! Output running model with parameter sets C
    !                                                                   ! dim1 = number of parameter sets 
    !                                                                   !      = size(ya) = size(yb)
    !                                                                   ! dim2 = number of parameters
    real(sp), dimension(size(yc,2)),              intent(out) :: si     ! Sobol index (main effect)
    real(sp), dimension(size(yc,2)),              intent(out) :: sti    ! Sobol index (total effect)
    integer(i4),                        optional, intent(in)  :: method ! (1) 'Saltelli2008'  needs Ci = B_A
                                                                        ! (2) 'Homma1996'     needs Ci = B_A
                                                                        ! (3) 'Sobol2007'     needs Ci = A_B and B_A --> NA
                                                                        ! (4) 'Saltelli2010'  needs Ci = A_B
                                                                        ! (5) 'Jansen1999'    needs Ci = A_B
                                                                        ! (6) 'Mai2012'       needs Ci = B_A
                                                                        ! (7) 'Mai2013'       needs Ci = A_B
                                                                        ! (8) 'Mai2014'       needs Ci = A_B
                                                                        ! (9) 'Mai1999'       needs Ci = A_B  ----> Default
    real(sp), dimension(size(yc,2), 2), optional, intent(out) :: smean  ! Mean Sobol index (main and total effect)
    !                                                                   ! dim1 = number of parameters
    !                                                                   ! dim2 = 2 (i.e. SI and STI)
    !                                                                   ! == (/ si, sti /) in 0d version !!
    real(sp), dimension(size(yc,2), 2), optional, intent(out) :: wmean  ! Variance weighted mean Sobol index 
    !                                                                   ! (main and total effect)
    !                                                                   ! dim1 = number of parameters
    !                                                                   ! dim2 = 2 (i.e. SI and STI)
    !                                                                   ! == (/ si, sti /) in 0d version !!

    ! local variables
    integer(i4)                          :: ii
    integer(i4)                          :: nsets, npara
    integer(i4)                          :: meth
    real(sp), dimension(:), allocatable  :: yab
    real(sp)                             :: var_a, var_b, var_ab
    real(sp)                             :: mean_a, mean_b, mean_ab

    if (present(method)) then
       meth = method
    else
       meth = 9
    end if

    nsets = size(yc,1)
    npara = size(yc,2)

    select case(meth)
       case(1)
          ! 'Saltelli2008' - The formulation presented in 'The Primer'. (yc=f(B_A))
          !                  Si  = (1/n*sum_j(f(A)_j*f(B_A^i)_j) - mean(f(A))^2)/var(f(A))
          !                  STi = (var(f(A))-(1/n*sum_j(f(B)_j*f(B_A^i)_j) - mean(f(A))^2))/var(f(A))
          var_a  = variance(ya)
          mean_a = sum( ya(:) ) / real(nsets,sp)
          if ( var_a .gt. tiny(1.0_sp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  =          ( dot_product( ya(:) , (yc(:,ii)) ) / real(nsets,sp) - mean_a**2 ) / var_a
                sti(ii) = 1.0_sp - ( dot_product( yb(:) , (yc(:,ii)) ) / real(nsets,sp) - mean_a**2 ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_sp
             sti(:) = 0.0_sp
          end if
       case(2)
          ! 'Homma1996'    - Si takes mixed A, B term for squared mean in nominator. (yc=f(B_A))
          !                  Si  = 1/n*sum_j(f(A)_j*(f(B_A^i)_j - f(B)_j))/var(f(A))
          !                  STi as Saltelli2008
          var_a  = variance(ya)
          mean_a = sum( ya(:) ) / real(nsets,sp)
          if ( var_a .gt. tiny(1.0_sp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  =          ( dot_product( ya(:) , (yc(:,ii)-yb(:)) ) / real(nsets,sp) ) / var_a
                sti(ii) = 1.0_sp - ( dot_product( yb(:) , (yc(:,ii)) ) / real(nsets,sp) - mean_a**2 ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_sp
             sti(:) = 0.0_sp
          end if
       case(3)
          stop 'mo_sobol_index: Sobol2007 would need f(A_B) and f(B_A). It is thus not implemented here.'
       case(4)
          ! 'Saltelli2010' - Si takes B, A_B instead of A, B_A (yc=f(A_B))
          !                  Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var(f(A))
          !                  STi = 1/n*sum_j(f(A)_j*(f(A)_j-f(A_B^i)_j))/var(f(A))
          var_a  = variance(ya)
          if ( var_a .gt. tiny(1.0_sp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  = ( dot_product( yb(:) , (yc(:,ii)-ya(:)) ) / real(nsets,sp) ) / var_a
                sti(ii) = ( dot_product( ya(:) , (ya(:)-yc(:,ii)) ) / real(nsets,sp) ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_sp
             sti(:) = 0.0_sp
          end if
       case(5)
          ! 'Jansen1999'   - Calculate Si and STi by expectation(variance) instead of variance(expectation) (yc=f(A_B))
          !                  Si  = (var(f(A)) - 1/2n*sum_j(f(B)_j - f(A_B^i)_j)^2)/var(f(A))
          !                  STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
          var_a  = variance(ya)
          if ( var_a .gt. tiny(1.0_sp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  = 1.0_sp - ( sum( (yb(:)-yc(:,ii))**2 ) / real(2*nsets,sp) ) / var_a
                sti(ii) =          ( sum( (ya(:)-yc(:,ii))**2 ) / real(2*nsets,sp) ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_sp
             sti(:) = 0.0_sp
          end if
       case(6)
          ! 'Mai2012'      - Si of Homma 1996 but denominator is variance of A and B (yc=f(B_A)) and
          !                  STi of Homma 1996 but mean and variance of f(A) replaced by f(B)
          !                  Si  = 1/n*sum_j(f(A)_j*(f(B_A^i)_j - f(B)_j))/var([f(A),f(B)])
          !                  STi = (var(f(A))-(1/n*sum_j(f(B)_j*f(B_A^i)_j) - mean(f(B)^2)))/var(f(B))
          allocate(yab(2*nsets))
          yab(      1:  nsets) = ya
          yab(nsets+1:2*nsets) = yb
          var_ab = variance(yab)
          var_b  = variance(yb)

          mean_ab = dot_product(ya,yb) / real(nsets,sp)
          mean_b  = sum(yb) / real(nsets,sp)

          if ( var_ab .gt. tiny(1.0_sp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  =          ( dot_product(ya(:),yc(:,ii)) / real(nsets,sp) - mean_ab)   / var_ab
                sti(ii) = 1.0_sp - ( dot_product(yb(:),yc(:,ii)) / real(nsets,sp) - mean_b**2) / var_b
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_sp
             sti(:) = 0.0_sp
          end if
          deallocate(yab)
       case(7)
          ! 'Mai2013'      - Si of Saltelli2010 but denominator is variance of A and B (yc=f(A_B))
          !                  Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var([f(A),f(B)])
          !                  STi = 1/n*sum_j(f(A)_j*(f(A)_j)-f(A_B^i)_j)/var(f(A))
          allocate(yab(2*nsets))
          yab(      1:  nsets) = ya
          yab(nsets+1:2*nsets) = yb
          var_ab = variance(yab)
          var_a  = variance(ya)

          if ( var_ab .gt. tiny(1.0_sp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  = ( dot_product( yb(:) , (yc(:,ii)-ya(:)) ) / real(nsets,sp) ) / var_ab
                sti(ii) = ( dot_product( ya(:) , (ya(:)-yc(:,ii)) ) / real(nsets,sp) ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_sp
             sti(:) = 0.0_sp
          end if
          deallocate(yab)
       case(8)
          ! 'Mai2014'      - SI of Jansen1999  but denominator is variance of A and B (yc=f(A_B))
          !                  Si  = (var([f(A),f(B)]) - 1/2n*sum_j(f(B)_j - f(A_B^i)_j)^2)/var([f(A),f(B)])
          !                  STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
          allocate(yab(2*nsets))
          yab(      1:  nsets) = ya
          yab(nsets+1:2*nsets) = yb
          var_ab = variance(yab)
          var_a  = variance(ya)

          if ( var_ab .gt. tiny(1.0_sp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  = 1.0_sp - ( sum( (yb(:)-yc(:,ii))**2 ) / real(2*nsets,sp) ) / var_ab
                sti(ii) =          ( sum( (ya(:)-yc(:,ii))**2 ) / real(2*nsets,sp) ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_sp
             sti(:) = 0.0_sp
          end if
          deallocate(yab)
       case(9)
          ! 'Mai1999'      - SI of Mai2013 (Saltelli2010 with var([f(A),f(B)])) and STi of Jansen1999 (yc=f(A_B))
          !                  Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var([f(A),f(B)])
          !                  STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
          allocate(yab(2*nsets))
          yab(      1:  nsets) = ya
          yab(nsets+1:2*nsets) = yb
          var_ab = variance(yab)
          var_a  = variance(ya)

          if ( var_ab .gt. tiny(1.0_sp) ) then
             ! model outputs are different (usual case)
             do ii=1, npara
                si(ii)  = ( dot_product( yb(:) , (yc(:,ii)-ya(:)) ) / real(nsets,sp) ) / var_ab
                sti(ii) = ( sum( (ya(:)-yc(:,ii))**2 ) / real(2*nsets,sp) ) / var_a
             end do
          else
             ! model outputs are equal (should never happen)
             si(:)  = 0.0_sp
             sti(:) = 0.0_sp
          end if
          deallocate(yab)
       case default
          stop 'mo_sobol_index: This method is not implemented!'
       end select

   if ( present(smean) ) then
       smean(:,1) = si(:)
       smean(:,2) = sti(:)
    end if

    if ( present(wmean) ) then
       wmean(:,1) = si(:)
       wmean(:,2) = sti(:)
    end if

  end subroutine sobol_index_0d_sp

  subroutine sobol_index_1d_dp(ya, yb, yc, si, sti, method, smean, wmean)

    real(dp), dimension(:,:),                     intent(in)  :: ya     ! Output running model with parameter sets A
    !                                                                   ! dim1 = number of parameter sets 
    !                                                                   ! dim2 = number of modeloutputs 
    !                                                                   !        (e.g. number of timesteps)
    real(dp), dimension(:,:),                     intent(in)  :: yb     ! Output running model with parameter sets B
    !                                                                   ! dim1 = number of parameter sets 
    !                                                                   ! dim2 = number of modeloutputs 
    !                                                                   !        (e.g. number of timesteps)
    real(dp), dimension(:,:,:),                   intent(in)  :: yc     ! Output running model with parameter sets C
    !                                                                   ! dim1 = number of parameter sets 
    !                                                                   !      = size(ya) = size(yb)
    !                                                                   ! dim2 = number of parameters
    !                                                                   ! dim3 = number of modeloutputs 
    !                                                                   !        (e.g. number of timesteps)
    real(dp), dimension(size(yc,2),size(yc,3)),   intent(out) :: si     ! Sobol index (main effect)
    !                                                                   ! dim1 = number of parameters
    !                                                                   ! dim2 = number of modeloutputs 
    !                                                                   !        (e.g. number of timesteps)
    real(dp), dimension(size(yc,2),size(yc,3)),   intent(out) :: sti    ! Sobol index (total effect)
    !                                                                   ! dim1 = number of parameters
    !                                                                   ! dim2 = number of modeloutputs 
    !                                                                   !        (e.g. number of timesteps)
    integer(i4),                        optional, intent(in)  :: method ! (1) 'Saltelli2008'  needs Ci = B_A
                                                                        ! (2) 'Homma1996'     needs Ci = B_A
                                                                        ! (3) 'Sobol2007'     needs Ci = A_B and B_A --> NA
                                                                        ! (4) 'Saltelli2010'  needs Ci = A_B
                                                                        ! (5) 'Jansen1999'    needs Ci = A_B
                                                                        ! (6) 'Mai2012'       needs Ci = B_A
                                                                        ! (7) 'Mai2013'       needs Ci = A_B
                                                                        ! (8) 'Mai2014'       needs Ci = A_B
                                                                        ! (9) 'Mai1999'       needs Ci = A_B  ----> Default
    real(dp), dimension(size(yc,2), 2), optional, intent(out) :: smean  ! Mean Sobol index (main and total effect)
    !                                                                   ! dim1 = number of parameters
    !                                                                   ! dim2 = 2 (i.e. SI and STI)
    real(dp), dimension(size(yc,2), 2), optional, intent(out) :: wmean  ! Variance weighted mean Sobol index 
    !                                                                   ! (main and total effect)
    !                                                                   ! dim1 = number of parameters
    !                                                                   ! dim2 = 2 (i.e. SI and STI)

    ! local variables
    integer(i4)                          :: ii, iout
    integer(i4)                          :: nsets, npara, nout
    integer(i4)                          :: meth
    real(dp), dimension(:), allocatable  :: yab
    real(dp)                             :: var_a, var_b, var_ab
    real(dp)                             :: mean_a, mean_b, mean_ab
    ! for averaging
    real(dp), dimension(:), allocatable  :: si_num, sti_num
    real(dp)                             :: si_denom, sti_denom
    real(dp)                             :: si_num_ii, sti_num_ii

    if (present(method)) then
       meth = method
    else
       meth = 9
    end if

    nsets = size(yc,1)
    npara = size(yc,2)
    nout  = size(yc,3)

    si(:,:)  = 0.0_dp
    sti(:,:) = 0.0_dp

    allocate(si_num(npara))
    allocate(sti_num(npara))

    si_num(:)  = 0.0_dp
    sti_num(:) = 0.0_dp
    si_denom   = 0.0_dp
    sti_denom  = 0.0_dp

    select case(meth)
    case(1)
          ! 'Saltelli2008' - The formulation presented in 'The Primer'. (yc=f(B_A))
          !                  Si  = (1/n*sum_j(f(A)_j*f(B_A^i)_j) - mean(f(A))^2)/var(f(A))
          !                  STi = (var(f(A))-(1/n*sum_j(f(B)_j*f(B_A^i)_j) - mean(f(A))^2))/var(f(A))
          do iout=1, nout
             var_a  = variance(ya(:,iout))
             mean_a  = sum(ya(:,iout)) / real(nsets,dp)

             si_denom   = si_denom  + var_a
             sti_denom  = sti_denom + var_a
             
             if ( var_a .gt. tiny(1.0_dp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( dot_product( ya(:,iout) , (yc(:,ii,iout)) ) / real(nsets,dp) - mean_a**2 )
                   sti_num_ii = ( dot_product( yb(:,iout) , (yc(:,ii,iout)) ) / real(nsets,dp) - mean_a**2 )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  =          si_num_ii  / var_a
                   sti(ii,iout) = 1.0_dp - sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_dp
                sti(:,iout) = 0.0_dp
             end if
          end do
       case(2)
          ! 'Homma1996'    - Si takes mixed A, B term for squared mean in nominator. (yc=f(B_A))
          !                  Si  = 1/n*sum_j(f(A)_j*(f(B_A^i)_j - f(B)_j))/var(f(A))
          !                  STi as Saltelli2008
          do iout=1, nout
             var_a  = variance(ya(:,iout))
             mean_a  = sum(ya(:,iout)) / real(nsets,dp)

             si_denom   = si_denom  + var_a
             sti_denom  = sti_denom + var_a
             
             if ( var_a .gt. tiny(1.0_dp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( dot_product( ya(:,iout) , (yc(:,ii,iout)-yb(:,iout)) ) / real(nsets,dp) )
                   sti_num_ii = ( dot_product( yb(:,iout) , (yc(:,ii,iout)) ) / real(nsets,dp) - mean_a**2 )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  =          si_num_ii  / var_a
                   sti(ii,iout) = 1.0_dp - sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_dp
                sti(:,iout) = 0.0_dp
             end if
          end do
       case(3)
          stop 'mo_sobol_index: Sobol2007 would need f(A_B) and f(B_A). It is thus not implemented here.'
       case(4)
          ! 'Saltelli2010' - Si takes B, A_B instead of A, B_A (yc=f(A_B))
          !                  Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var(f(A))
          !                  STi = 1/n*sum_j(f(A)_j*(f(A)_j-f(A_B^i)_j))/var(f(A))
          do iout=1, nout
             var_a  = variance(ya(:,iout))

             si_denom   = si_denom  + var_a
             sti_denom  = sti_denom + var_a
             
             if ( var_a .gt. tiny(1.0_dp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( dot_product( yb(:,iout) , (yc(:,ii,iout)-ya(:,iout)) ) / real(nsets,dp) )
                   sti_num_ii = ( dot_product( ya(:,iout) , (ya(:,iout)-yc(:,ii,iout)) ) / real(nsets,dp) )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  =          si_num_ii  / var_a
                   sti(ii,iout) =          sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_dp
                sti(:,iout) = 0.0_dp
             end if
          end do
       case(5)
          ! 'Jansen1999'   - Calculate Si and STi by expectation(variance) instead of variance(expectation) (yc=f(A_B))
          !                  Si  = (var(f(A)) - 1/2n*sum_j(f(B)_j - f(A_B^i)_j)^2)/var(f(A))
          !                  STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
          do iout=1, nout
             var_a  = variance(ya(:,iout))

             si_denom   = si_denom  + var_a
             sti_denom  = sti_denom + var_a
             
             if ( var_a .gt. tiny(1.0_dp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( sum( (yb(:,iout)-yc(:,ii,iout))**2 ) / real(2*nsets,dp) )
                   sti_num_ii = ( sum( (ya(:,iout)-yc(:,ii,iout))**2 ) / real(2*nsets,dp) )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  = 1.0_dp - si_num_ii  / var_a
                   sti(ii,iout) =          sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_dp
                sti(:,iout) = 0.0_dp
             end if
          end do
       case(6)
          ! 'Mai2012'      - Si of Homma 1996 but denominator is variance of A and B (yc=f(B_A)) and
          !                  STi of Homma 1996 but mean and variance of f(A) replaced by f(B)
          !                  Si  = 1/n*sum_j(f(A)_j*(f(B_A^i)_j - f(B)_j))/var([f(A),f(B)])
          !                  STi = (var(f(A))-(1/n*sum_j(f(B)_j*f(B_A^i)_j) - mean(f(B)^2)))/var(f(B))
          allocate(yab(2*nsets))
          do iout=1, nout
             yab(      1:  nsets) = ya(:,iout)
             yab(nsets+1:2*nsets) = yb(:,iout)
             var_ab = variance(yab)
             var_b  = variance(yb(:,iout))
             
             mean_ab = dot_product(ya(:,iout),yb(:,iout)) / real(nsets,dp)
             mean_b  = sum(yb(:,iout)) / real(nsets,dp)

             si_denom   = si_denom  + var_ab
             sti_denom  = sti_denom + var_b
             
             if ( var_ab .gt. tiny(1.0_dp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( dot_product(ya(:,iout),yc(:,ii,iout)) / real(nsets,dp) - mean_ab ) 
                   sti_num_ii = ( dot_product(yb(:,iout),yc(:,ii,iout)) / real(nsets,dp) - mean_b**2 )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  =          si_num_ii  / var_ab
                   sti(ii,iout) = 1.0_dp - sti_num_ii / var_b
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_dp
                sti(:,iout) = 0.0_dp
             end if
          end do
          deallocate(yab)
       case(7)
          ! 'Mai2013'      - Si of Saltelli2010 but denominator is variance of A and B (yc=f(A_B))
          !                  Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var([f(A),f(B)])
          !                  STi = 1/n*sum_j(f(A)_j*(f(A)_j)-f(A_B^i)_j)/var(f(A))
          allocate(yab(2*nsets))
          do iout=1, nout
             yab(      1:  nsets) = ya(:,iout)
             yab(nsets+1:2*nsets) = yb(:,iout)
             var_ab = variance(yab)
             var_a  = variance(ya(:,iout))

             si_denom   = si_denom  + var_ab
             sti_denom  = sti_denom + var_a
             
             if ( var_ab .gt. tiny(1.0_dp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( dot_product( yb(:,iout) , (yc(:,ii,iout)-ya(:,iout)) ) / real(nsets,dp) )
                   sti_num_ii = ( dot_product( ya(:,iout) , (ya(:,iout)-yc(:,ii,iout)) ) / real(nsets,dp) )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  =          si_num_ii  / var_ab
                   sti(ii,iout) =          sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_dp
                sti(:,iout) = 0.0_dp
             end if
          end do
          deallocate(yab)
       case(8)
          ! 'Mai2014'      - SI of Jansen1999  but denominator is variance of A and B (yc=f(A_B))
          !                  Si  = (var([f(A),f(B)]) - 1/2n*sum_j(f(B)_j - f(A_B^i)_j)^2)/var([f(A),f(B)])
          !                  STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
          allocate(yab(2*nsets))
          do iout=1, nout
             yab(      1:  nsets) = ya(:,iout)
             yab(nsets+1:2*nsets) = yb(:,iout)
             var_ab = variance(yab)
             var_a  = variance(ya(:,iout))

             si_denom   = si_denom  + var_ab
             sti_denom  = sti_denom + var_a
             
             if ( var_ab .gt. tiny(1.0_dp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( sum( (yb(:,iout)-yc(:,ii,iout))**2 ) / real(2*nsets,dp) )
                   sti_num_ii = ( sum( (ya(:,iout)-yc(:,ii,iout))**2 ) / real(2*nsets,dp) )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  = 1.0_dp - si_num_ii  / var_ab
                   sti(ii,iout) =          sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_dp
                sti(:,iout) = 0.0_dp
             end if
          end do
          deallocate(yab)
       case(9)
          ! 'Mai1999'      - SI of Mai2013 (Saltelli2010 with var([f(A),f(B)])) and STi of Jansen1999 (yc=f(A_B))
          !                  Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var([f(A),f(B)])
          !                  STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
          allocate(yab(2*nsets))
          do iout=1, nout
             yab(      1:  nsets) = ya(:,iout)
             yab(nsets+1:2*nsets) = yb(:,iout)
             var_ab = variance(yab)
             var_a  = variance(ya(:,iout))

             si_denom   = si_denom  + var_ab
             sti_denom  = sti_denom + var_a
             
             if ( var_ab .gt. tiny(1.0_dp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( dot_product( yb(:,iout) , (yc(:,ii,iout)-ya(:,iout)) ) / real(nsets,dp) )
                   sti_num_ii = ( sum( (ya(:,iout)-yc(:,ii,iout))**2 ) / real(2*nsets,dp) )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  =          si_num_ii  / var_ab
                   sti(ii,iout) =          sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_dp
                sti(:,iout) = 0.0_dp
             end if
          end do
          deallocate(yab)
       case default
          stop 'mo_sobol_index: This method is not implemented!'
    end select

    if (present(smean)) then
       smean = 0.0_dp
       smean(:,1) = sum(si, dim=2) / real(nout,dp)
       smean(:,2) = sum(sti,dim=2) / real(nout,dp)
    end if

    if (present(wmean)) then
       wmean = 0.0_dp
       select case(meth)
          case(1,2,4,6,7,9)
             if ( si_denom .gt. tiny(1.0_dp) ) then
                wmean(:,1) =          si_num(:)  / si_denom
             end if
          case(5,8)
             if ( si_denom .gt. tiny(1.0_dp) ) then
                wmean(:,1) = 1.0_dp - si_num(:)  / si_denom
             end if
          case default
             stop 'mo_sobol_index: This averaging method is not implemented!'
       end select
       select case(meth)
          case(1,2,6)
             if ( sti_denom .gt. tiny(1.0_dp) ) then
                wmean(:,2) = 1.0_dp - sti_num(:) / sti_denom
             end if
          case(4,5,7,8,9)
             if ( sti_denom .gt. tiny(1.0_dp) ) then
                wmean(:,2) =          sti_num(:) / sti_denom
             end if
          case default
             stop 'mo_sobol_index: This averaging method is not implemented!'
       end select
    end if
       
    deallocate(si_num, sti_num)

  end subroutine sobol_index_1d_dp

  subroutine sobol_index_1d_sp(ya, yb, yc, si, sti, method, smean, wmean)

    real(sp), dimension(:,:),                     intent(in)  :: ya     ! Output running model with parameter sets A
    !                                                                   ! dim1 = number of parameter sets 
    !                                                                   ! dim2 = number of modeloutputs 
    !                                                                   !        (e.g. number of timesteps)
    real(sp), dimension(:,:),                     intent(in)  :: yb     ! Output running model with parameter sets B
    !                                                                   ! dim1 = number of parameter sets 
    !                                                                   ! dim2 = number of modeloutputs 
    !                                                                   !        (e.g. number of timesteps)
    real(sp), dimension(:,:,:),                   intent(in)  :: yc     ! Output running model with parameter sets C
    !                                                                   ! dim1 = number of parameter sets 
    !                                                                   !      = size(ya) = size(yb)
    !                                                                   ! dim2 = number of parameters
    !                                                                   ! dim3 = number of modeloutputs 
    !                                                                   !        (e.g. number of timesteps)
    real(sp), dimension(size(yc,2),size(yc,3)),   intent(out) :: si     ! Sobol index (main effect)
    !                                                                   ! dim1 = number of parameters
    !                                                                   ! dim2 = number of modeloutputs 
    !                                                                   !        (e.g. number of timesteps)
    real(sp), dimension(size(yc,2),size(yc,3)),   intent(out) :: sti    ! Sobol index (total effect)
    !                                                                   ! dim1 = number of parameters
    !                                                                   ! dim2 = number of modeloutputs 
    !                                                                   !        (e.g. number of timesteps)
    integer(i4),                        optional, intent(in)  :: method ! (1) 'Saltelli2008'  needs Ci = B_A
                                                                        ! (2) 'Homma1996'     needs Ci = B_A
                                                                        ! (3) 'Sobol2007'     needs Ci = A_B and B_A --> NA
                                                                        ! (4) 'Saltelli2010'  needs Ci = A_B
                                                                        ! (5) 'Jansen1999'    needs Ci = A_B
                                                                        ! (6) 'Mai2012'       needs Ci = B_A
                                                                        ! (7) 'Mai2013'       needs Ci = A_B
                                                                        ! (8) 'Mai2014'       needs Ci = A_B
                                                                        ! (9) 'Mai1999'       needs Ci = A_B  ----> Default
    real(sp), dimension(size(yc,2), 2), optional, intent(out) :: smean  ! Mean Sobol index (main and total effect)
    !                                                                   ! dim1 = number of parameters
    !                                                                   ! dim2 = 2 (i.e. SI and STI)
    real(sp), dimension(size(yc,2), 2), optional, intent(out) :: wmean  ! Variance weighted mean Sobol index 
    !                                                                   ! (main and total effect)
    !                                                                   ! dim1 = number of parameters
    !                                                                   ! dim2 = 2 (i.e. SI and STI)

    ! local variables
    integer(i4)                          :: ii, iout
    integer(i4)                          :: nsets, npara, nout
    integer(i4)                          :: meth
    real(sp), dimension(:), allocatable  :: yab
    real(sp)                             :: var_a, var_b, var_ab
    real(sp)                             :: mean_a, mean_b, mean_ab
    ! for averaging
    real(sp), dimension(:), allocatable  :: si_num, sti_num
    real(sp)                             :: si_denom, sti_denom
    real(sp)                             :: si_num_ii, sti_num_ii

    if (present(method)) then
       meth = method
    else
       meth = 9
    end if

    nsets = size(yc,1)
    npara = size(yc,2)
    nout  = size(yc,3)

    si(:,:)  = 0.0_sp
    sti(:,:) = 0.0_sp

    allocate(si_num(npara))
    allocate(sti_num(npara))

    si_num(:)  = 0.0_sp
    sti_num(:) = 0.0_sp
    si_denom   = 0.0_sp
    sti_denom  = 0.0_sp

    select case(meth)
    case(1)
          ! 'Saltelli2008' - The formulation presented in 'The Primer'. (yc=f(B_A))
          !                  Si  = (1/n*sum_j(f(A)_j*f(B_A^i)_j) - mean(f(A))^2)/var(f(A))
          !                  STi = (var(f(A))-(1/n*sum_j(f(B)_j*f(B_A^i)_j) - mean(f(A))^2))/var(f(A))
          do iout=1, nout
             var_a  = variance(ya(:,iout))
             mean_a  = sum(ya(:,iout)) / real(nsets,sp)

             si_denom   = si_denom  + var_a
             sti_denom  = sti_denom + var_a
             
             if ( var_a .gt. tiny(1.0_sp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( dot_product( ya(:,iout) , (yc(:,ii,iout)) ) / real(nsets,sp) - mean_a**2 )
                   sti_num_ii = ( dot_product( yb(:,iout) , (yc(:,ii,iout)) ) / real(nsets,sp) - mean_a**2 )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  =          si_num_ii  / var_a
                   sti(ii,iout) = 1.0_sp - sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_sp
                sti(:,iout) = 0.0_sp
             end if
          end do
       case(2)
          ! 'Homma1996'    - Si takes mixed A, B term for squared mean in nominator. (yc=f(B_A))
          !                  Si  = 1/n*sum_j(f(A)_j*(f(B_A^i)_j - f(B)_j))/var(f(A))
          !                  STi as Saltelli2008
          do iout=1, nout
             var_a  = variance(ya(:,iout))
             mean_a  = sum(ya(:,iout)) / real(nsets,sp)

             si_denom   = si_denom  + var_a
             sti_denom  = sti_denom + var_a
             
             if ( var_a .gt. tiny(1.0_sp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( dot_product( ya(:,iout) , (yc(:,ii,iout)-yb(:,iout)) ) / real(nsets,sp) )
                   sti_num_ii = ( dot_product( yb(:,iout) , (yc(:,ii,iout)) ) / real(nsets,sp) - mean_a**2 )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  =          si_num_ii  / var_a
                   sti(ii,iout) = 1.0_sp - sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_sp
                sti(:,iout) = 0.0_sp
             end if
          end do
       case(3)
          stop 'mo_sobol_index: Sobol2007 would need f(A_B) and f(B_A). It is thus not implemented here.'
       case(4)
          ! 'Saltelli2010' - Si takes B, A_B instead of A, B_A (yc=f(A_B))
          !                  Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var(f(A))
          !                  STi = 1/n*sum_j(f(A)_j*(f(A)_j-f(A_B^i)_j))/var(f(A))
          do iout=1, nout
             var_a  = variance(ya(:,iout))

             si_denom   = si_denom  + var_a
             sti_denom  = sti_denom + var_a
             
             if ( var_a .gt. tiny(1.0_sp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( dot_product( yb(:,iout) , (yc(:,ii,iout)-ya(:,iout)) ) / real(nsets,sp) )
                   sti_num_ii = ( dot_product( ya(:,iout) , (ya(:,iout)-yc(:,ii,iout)) ) / real(nsets,sp) )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  =          si_num_ii  / var_a
                   sti(ii,iout) =          sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_sp
                sti(:,iout) = 0.0_sp
             end if
          end do
       case(5)
          ! 'Jansen1999'   - Calculate Si and STi by expectation(variance) instead of variance(expectation) (yc=f(A_B))
          !                  Si  = (var(f(A)) - 1/2n*sum_j(f(B)_j - f(A_B^i)_j)^2)/var(f(A))
          !                  STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
          do iout=1, nout
             var_a  = variance(ya(:,iout))

             si_denom   = si_denom  + var_a
             sti_denom  = sti_denom + var_a
             
             if ( var_a .gt. tiny(1.0_sp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( sum( (yb(:,iout)-yc(:,ii,iout))**2 ) / real(2*nsets,sp) )
                   sti_num_ii = ( sum( (ya(:,iout)-yc(:,ii,iout))**2 ) / real(2*nsets,sp) )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  = 1.0_sp - si_num_ii  / var_a
                   sti(ii,iout) =          sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_sp
                sti(:,iout) = 0.0_sp
             end if
          end do
       case(6)
          ! 'Mai2012'      - Si of Homma 1996 but denominator is variance of A and B (yc=f(B_A)) and
          !                  STi of Homma 1996 but mean and variance of f(A) replaced by f(B)
          !                  Si  = 1/n*sum_j(f(A)_j*(f(B_A^i)_j - f(B)_j))/var([f(A),f(B)])
          !                  STi = (var(f(A))-(1/n*sum_j(f(B)_j*f(B_A^i)_j) - mean(f(B)^2)))/var(f(B))
          allocate(yab(2*nsets))
          do iout=1, nout
             yab(      1:  nsets) = ya(:,iout)
             yab(nsets+1:2*nsets) = yb(:,iout)
             var_ab = variance(yab)
             var_b  = variance(yb(:,iout))
             
             mean_ab = dot_product(ya(:,iout),yb(:,iout)) / real(nsets,sp)
             mean_b  = sum(yb(:,iout)) / real(nsets,sp)

             si_denom   = si_denom  + var_ab
             sti_denom  = sti_denom + var_b
             
             if ( var_ab .gt. tiny(1.0_sp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( dot_product(ya(:,iout),yc(:,ii,iout)) / real(nsets,sp) - mean_ab ) 
                   sti_num_ii = ( dot_product(yb(:,iout),yc(:,ii,iout)) / real(nsets,sp) - mean_b**2 )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  =          si_num_ii  / var_ab
                   sti(ii,iout) = 1.0_sp - sti_num_ii / var_b
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_sp
                sti(:,iout) = 0.0_sp
             end if
          end do
          deallocate(yab)
       case(7)
          ! 'Mai2013'      - Si of Saltelli2010 but denominator is variance of A and B (yc=f(A_B))
          !                  Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var([f(A),f(B)])
          !                  STi = 1/n*sum_j(f(A)_j*(f(A)_j)-f(A_B^i)_j)/var(f(A))
          allocate(yab(2*nsets))
          do iout=1, nout
             yab(      1:  nsets) = ya(:,iout)
             yab(nsets+1:2*nsets) = yb(:,iout)
             var_ab = variance(yab)
             var_a  = variance(ya(:,iout))

             si_denom   = si_denom  + var_ab
             sti_denom  = sti_denom + var_a
             
             if ( var_ab .gt. tiny(1.0_sp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( dot_product( yb(:,iout) , (yc(:,ii,iout)-ya(:,iout)) ) / real(nsets,sp) )
                   sti_num_ii = ( dot_product( ya(:,iout) , (ya(:,iout)-yc(:,ii,iout)) ) / real(nsets,sp) )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  =          si_num_ii  / var_ab
                   sti(ii,iout) =          sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_sp
                sti(:,iout) = 0.0_sp
             end if
          end do
          deallocate(yab)
       case(8)
          ! 'Mai2014'      - SI of Jansen1999  but denominator is variance of A and B (yc=f(A_B))
          !                  Si  = (var([f(A),f(B)]) - 1/2n*sum_j(f(B)_j - f(A_B^i)_j)^2)/var([f(A),f(B)])
          !                  STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
          allocate(yab(2*nsets))
          do iout=1, nout
             yab(      1:  nsets) = ya(:,iout)
             yab(nsets+1:2*nsets) = yb(:,iout)
             var_ab = variance(yab)
             var_a  = variance(ya(:,iout))

             si_denom   = si_denom  + var_ab
             sti_denom  = sti_denom + var_a
             
             if ( var_ab .gt. tiny(1.0_sp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( sum( (yb(:,iout)-yc(:,ii,iout))**2 ) / real(2*nsets,sp) )
                   sti_num_ii = ( sum( (ya(:,iout)-yc(:,ii,iout))**2 ) / real(2*nsets,sp) )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  = 1.0_sp - si_num_ii  / var_ab
                   sti(ii,iout) =          sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_sp
                sti(:,iout) = 0.0_sp
             end if
          end do
          deallocate(yab)
       case(9)
          ! 'Mai1999'      - SI of Mai2013 (Saltelli2010 with var([f(A),f(B)])) and STi of Jansen1999 (yc=f(A_B))
          !                  Si  = 1/n*sum_j(f(B)_j*(f(A_B^i)_j - f(A)_j))/var([f(A),f(B)])
          !                  STi = 1/2n*sum_j(f(A)_j - f(A_B^i)_j)^2/var(f(A))
          allocate(yab(2*nsets))
          do iout=1, nout
             yab(      1:  nsets) = ya(:,iout)
             yab(nsets+1:2*nsets) = yb(:,iout)
             var_ab = variance(yab)
             var_a  = variance(ya(:,iout))

             si_denom   = si_denom  + var_ab
             sti_denom  = sti_denom + var_a
             
             if ( var_ab .gt. tiny(1.0_sp) ) then
                ! model outputs are different (usual case)
                do ii=1, npara
                   si_num_ii  = ( dot_product( yb(:,iout) , (yc(:,ii,iout)-ya(:,iout)) ) / real(nsets,sp) )
                   sti_num_ii = ( sum( (ya(:,iout)-yc(:,ii,iout))**2 ) / real(2*nsets,sp) )

                   si_num(ii)  = si_num(ii)  + si_num_ii
                   sti_num(ii) = sti_num(ii) + sti_num_ii
                   
                   si(ii,iout)  =          si_num_ii  / var_ab
                   sti(ii,iout) =          sti_num_ii / var_a
                end do
             else
                ! model outputs are equal (should never happen)
                si(:,iout)  = 0.0_sp
                sti(:,iout) = 0.0_sp
             end if
          end do
          deallocate(yab)
       case default
          stop 'mo_sobol_index: This method is not implemented!'
    end select

    if (present(smean)) then
       smean = 0.0_sp
       smean(:,1) = sum(si, dim=2) / real(nout,sp)
       smean(:,2) = sum(sti,dim=2) / real(nout,sp)
    end if

    if (present(wmean)) then
       wmean = 0.0_sp
       select case(meth)
          case(1,2,4,6,7,9)
             if ( si_denom .gt. tiny(1.0_sp) ) then
                wmean(:,1) =          si_num(:)  / si_denom
             end if
          case(5,8)
             if ( si_denom .gt. tiny(1.0_sp) ) then
                wmean(:,1) = 1.0_sp - si_num(:)  / si_denom
             end if
          case default
             stop 'mo_sobol_index: This averaging method is not implemented!'
       end select
       select case(meth)
          case(1,2,6)
             if ( sti_denom .gt. tiny(1.0_sp) ) then
                wmean(:,2) = 1.0_sp - sti_num(:) / sti_denom
             end if
          case(4,5,7,8,9)
             if ( sti_denom .gt. tiny(1.0_sp) ) then
                wmean(:,2) =          sti_num(:) / sti_denom
             end if
          case default
             stop 'mo_sobol_index: This averaging method is not implemented!'
       end select
    end if
       
    deallocate(si_num, sti_num)

  end subroutine sobol_index_1d_sp

END MODULE mo_sobol_index

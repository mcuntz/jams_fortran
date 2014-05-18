!> \file mo_delsa.f90

!> \brief Distributed Evaluation of Local Sensitivity Analysis

!> \details This module calculates first order parameter sensitivity
!>           using DELSA (Rakovec et al., 2014, WRR).

!> \authors Oldrich Rakovec
!> \date May 2014

MODULE mo_delsa

  ! This module calculates first order parameter sensitivity using DELSA (Rakovec et al., 2014, WRR)

  ! Written  Oldrich Rakovec, May 2014

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

  ! Copyright 2014 Oldrich Rakovec

  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE

  PUBLIC :: delsa   ! delsa = distributed evaluation of local sensitivity analysis

  ! ------------------------------------------------------------------
  !     NAME
  !         delsa

  !     PURPOSE
  !         Calculates first order parameter sensitivity using the DELSA method thoughout the parameter space
  !
  !>        \brief DELSA.
  !
  !>        \details Calculates first order parameter sensitivity using the DELSA method thoughout the parameter space:
  !>        S_{L1}^j=\frac{\left(\frac{\partial \Psi_l}{\partial \theta_j} \right)^2 s_j^2}{V_L(\Psi)}
  !>        for more details see derivation of equation (13) in Rakovec et al. (2014, WRR)
  !
  !
  !     INTENT(IN)
  !>        \param[in] "real(sp/dp) :: parbase(:,:)"    array of base model parameters, first dimension is the sample size,
  !>                                                    i.e. number of parameter sets, second dimension is the number of
  !>                                                    parameters, which are recommended to be sampled using the Sobol'
  !>                                                    quasi random sequence (see for example mo_sobol in FORTRAN_chs_lib)
  !>        \param[in] "real(sp/dp) :: parpert(:,:)"    array of perturbed model parameters, first dimension is the sample size,
  !>                                                    i.e. number of parameter sets, second dimension is the number of
  !>                                                    parameters, which are recommended to be sampled using the Sobol'
  !>                                                    quasi random sequence (see mo_sobol), perturbation e.g. 1% parameter change
  !>                                                    as often done in literature (e.g. Hill, M. C. and C. R. Tiedeman; 2007)
  !>        \param[in] "real(sp/dp) :: outbase(:)"      baserun model output obtained using parbase(:,:)
  !>        \param[in] "real(sp/dp) :: outpert(:,:)"    perturbed model output obtained using parpert(:,:)
  !>        \param[in] "real(sp/dp) :: varprior(:)"     prior variance obtained using equation for uniform distribution:
  !>                                                    $\frac{1}{12}(\theta_{j,max}-\theta_{j,min})^2$

  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(sp/dp) :: delsafirst(:,:)" delsa, first order Local first order sensitivity index
  !>                                                     (see eq. 13 in Rakovec et al 2014), dimensions correspond to parbase(:,:)

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
  !         None
  !
  !     RESTRICTIONS
  !>       \note Works only with scalar model output.\n
  !>             No averaging of DELSA results implemented yet. \n
  !>             No check for outliers of the model output. \n
  !>             No masking of parameters possible.
  !
  !     EXAMPLE
  !         -> see example in test directory

  !     LITERATURE
  !         Rakovec, O., M. C. Hill, M. P. Clark, A. H. Weerts, A. J. Teuling, R. Uijlenhoet (2014), Distributed
  !                Evaluation of Local Sensitivity Analysis (DELSA), with application to hydrologic models,
  !                Water Resour. Res., 50, 1-18, doi:10.1002/2013WR014063.
  !         Hill, M.~C., and C.~R. Tiedeman (2007), Effective groundwater model calibration: with analysis of data,
  !                sensitivities, prediction and uncertainty, 455 pp., Wiley.

  !     HISTORY
  !>        \author Oldrich Rakovec
  !>        \date May 2014
  !         Modified,

  INTERFACE delsa
     MODULE PROCEDURE delsa_sp, delsa_dp
  END INTERFACE delsa

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  SUBROUTINE delsa_sp(parbase,parpert,outbase,outpert,varprior,delsafirst)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),                             INTENT(IN)  :: parbase      ! array of base model parameters
    REAL(sp), DIMENSION(:,:),                             INTENT(IN)  :: parpert      ! array of perturbed model parameters
    REAL(sp), DIMENSION(:),                               INTENT(IN)  :: outbase      ! model output using the base runs
    REAL(sp), DIMENSION(:,:),                             INTENT(IN)  :: outpert      ! model output using the perturbed runs
    REAL(sp), DIMENSION(:),                               INTENT(IN)  :: varprior     ! prior parameter variance
    REAL(sp), DIMENSION(size(parbase,1),size(parbase,2)), INTENT(OUT) :: delsafirst   ! first-order DELSA results

    ! local variables
    INTEGER(I4)                                                       :: Nsamp        ! sample size
    INTEGER(I4)                                                       :: Kpar         ! number of model parameters
    INTEGER(I4)                                                       :: rsamp,jpar   ! counter
    REAL(sp), DIMENSION(size(parbase,1),size(parbase,2))              :: deriv        ! array of parameter derivatives (eq.10)
    !                                                                                 ! see nominator of eq. 13 in Rakovec et al 2014)
    REAL(sp), DIMENSION(size(parbase,1),size(parbase,2))              :: varfir       ! array of lcal first order variance,
    !                                                                                 ! see nominator of eq. 13 in Rakovec et al 2014)
    REAL(sp), DIMENSION(size(parbase,1))                              :: vartot       ! local equation total variance (see eq. 12)

    Nsamp = size(parbase,1)
    Kpar  = size(parbase,2)

    ! initialization
    vartot(:)   = 0.0_sp
    deriv(:,:)  = 0.0_sp
    varfir(:,:) = 0.0_sp

    loop_sets: do rsamp = 1,Nsamp ! looping over parameter sets

       loop_para: do jpar = 1,Kpar ! looping over model parameters

          ! calculate derivative (see eq. 10 in Rakovec et al 2014)
          deriv(rsamp,jpar) = (outpert(rsamp,jpar) -  outbase(rsamp)) / (parpert(rsamp,jpar) - parbase(rsamp,jpar))

          ! calculate local first order variance (see nominator of eq. 13 in Rakovec et al 2014)
          varfir(rsamp,jpar) =  (deriv(rsamp,jpar)**2.0_i4)*(varprior(jpar))

          ! calculate local equation total variance (see eq. 12 in Rakovec et al 2014)
          vartot(rsamp) = vartot(rsamp) + varfir(rsamp,jpar)

       end do loop_para ! End loop over model parameters

       ! finally when vartot is complete, we can calculate delsa,
       ! here as the first order Local first order sensitivity index (see eq. 13 in Rakovec et al 2014)
       delsafirst(rsamp,:) =  varfir(rsamp,:) /  vartot(rsamp)

    end do loop_sets ! End loop over parameter sets

  END SUBROUTINE delsa_sp


  SUBROUTINE delsa_dp(parbase,parpert,outbase,outpert,varprior,delsafirst)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),                             INTENT(IN)  :: parbase      ! array of base model parameters
    REAL(dp), DIMENSION(:,:),                             INTENT(IN)  :: parpert      ! array of perturbed model parameters
    REAL(dp), DIMENSION(:),                               INTENT(IN)  :: outbase      ! model output using the base runs
    REAL(dp), DIMENSION(:,:),                             INTENT(IN)  :: outpert      ! model output using the perturbed runs
    REAL(dp), DIMENSION(:),                               INTENT(IN)  :: varprior     ! prior parameter variance
    REAL(dp), DIMENSION(size(parbase,1),size(parbase,2)), INTENT(OUT) :: delsafirst   ! first-order DELSA results

    ! local variables
    INTEGER(I4)                                                       :: Nsamp        ! sample size
    INTEGER(I4)                                                       :: Kpar         ! number of model parameters
    INTEGER(I4)                                                       :: rsamp,jpar   ! counter
    REAL(dp), DIMENSION(size(parbase,1),size(parbase,2))              :: deriv        ! array of parameter derivatives (eq.10)
    !                                                                                 ! see nominator of eq. 13 in Rakovec et al 2014)
    REAL(dp), DIMENSION(size(parbase,1),size(parbase,2))              :: varfir       ! array of lcal first order variance,
    !                                                                                 ! see nominator of eq. 13 in Rakovec et al 2014)
    REAL(dp), DIMENSION(size(parbase,1))                              :: vartot       ! local equation total variance (see eq. 12)

    Nsamp = size(parbase,1)
    Kpar  = size(parbase,2)

    ! initialization
    vartot(:)   = 0.0_dp
    deriv(:,:)  = 0.0_dp
    varfir(:,:) = 0.0_dp

    loop_sets: do rsamp = 1,Nsamp ! looping over parameter sets

       loop_para: do jpar = 1,Kpar ! looping over model parameters

          ! calculate derivative (see eq. 10 in Rakovec et al 2014)
          deriv(rsamp,jpar) = (outpert(rsamp,jpar) -  outbase(rsamp)) / (parpert(rsamp,jpar) - parbase(rsamp,jpar))

          ! calculate local first order variance (see nominator of eq. 13 in Rakovec et al 2014)
          varfir(rsamp,jpar) =  (deriv(rsamp,jpar)**2.0_i4)*(varprior(jpar))

          ! calculate local equation total variance (see eq. 12 in Rakovec et al 2014)
          vartot(rsamp) = vartot(rsamp) + varfir(rsamp,jpar)

       end do loop_para ! End loop over model parameters

       ! finally when vartot is complete, we can calculate delsa,
       ! here as the first order Local first order sensitivity index (see eq. 13 in Rakovec et al 2014)
       delsafirst(rsamp,:) =  varfir(rsamp,:) /  vartot(rsamp)

    end do loop_sets ! End loop over parameter sets

  END SUBROUTINE delsa_dp

END MODULE mo_delsa

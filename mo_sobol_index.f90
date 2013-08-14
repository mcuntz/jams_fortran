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

  ! Copyright 2013 Juliane Mai

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
  !>                   - a parameter set A \n
  !>                   - a parameter set B which is independent to A \n
  !>                   - n parameter sets Ci where Ci equals parameter set B except parameter i 
  !>                     which is the value of parameter set A \n
  !>                 This has to be repeated several times.  
  !>                 The model response can be either a scalar or a vector (e.g. a timeseries). 
  !>                 If the model response is a vector, the Si and STi are calculated for each time point 
  !>                 independently such that the returned Si and STi are also vectors.
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

  subroutine sobol_index_0d_dp(ya, yb, yc, si, sti, smean, wmean)

    real(dp), dimension(:),                       intent(in)  :: ya     ! Output running model with parameter sets A
    real(dp), dimension(:),                       intent(in)  :: yb     ! Output running model with parameter sets B
    real(dp), dimension(:,:),                     intent(in)  :: yc     ! Output running model with parameter sets C
    !                                                                   ! dim1 = number of parameter sets 
    !                                                                   !      = size(ya) = size(yb)
    !                                                                   ! dim2 = number of parameters
    real(dp), dimension(size(yc,2)),              intent(out) :: si     ! Sobol index (main effect)
    real(dp), dimension(size(yc,2)),              intent(out) :: sti    ! Sobol index (total effect)
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
    real(dp), dimension(:), allocatable  :: yab
    real(dp)                             :: var_ab, var_b
    real(dp)                             :: f0_ab, f0_b2

    nsets = size(yc,1)
    npara = size(yc,2)

    allocate(yab(2*nsets))
    yab(      1:  nsets) = ya
    yab(nsets+1:2*nsets) = yb

    var_ab = variance(yab)
    var_b  = variance(yb)

    f0_ab  = dot_product(ya,yb) / real(nsets,dp)
    f0_b2  = (sum(yb) / real(nsets,dp) )**2

    if ( var_ab .gt. tiny(1.0_dp) ) then
       ! model outputs are different (usual case)
       do ii=1, npara
          si(ii)  =          ( dot_product(ya(:),yc(:,ii)) / real(nsets,dp) - f0_ab) / var_ab
          sti(ii) = 1.0_dp - ( dot_product(yb(:),yc(:,ii)) / real(nsets,dp) - f0_b2) / var_b
       end do
    else
       ! model outputs are equal (should never happen)
       si(:)  = 0.0_dp
       sti(:) = 0.0_dp
    end if

   if ( present(smean) ) then
       smean(:,1) = si(:)
       smean(:,2) = sti(:)
    end if

    if ( present(wmean) ) then
       wmean(:,1) = si(:)
       wmean(:,2) = sti(:)
    end if

    deallocate(yab)

  end subroutine sobol_index_0d_dp

  subroutine sobol_index_0d_sp(ya, yb, yc, si, sti, smean, wmean)

    real(sp), dimension(:),                       intent(in)  :: ya     ! Output running model with parameter sets A
    real(sp), dimension(:),                       intent(in)  :: yb     ! Output running model with parameter sets B
    real(sp), dimension(:,:),                     intent(in)  :: yc     ! Output running model with parameter sets C
    !                                                                   ! dim1 = number of parameter sets 
    !                                                                   !      = size(ya) = size(yb)
    !                                                                   ! dim2 = number of parameters
    real(sp), dimension(size(yc,2)),              intent(out) :: si     ! Sobol index (main effect)
    real(sp), dimension(size(yc,2)),              intent(out) :: sti    ! Sobol index (total effect)
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
    real(sp), dimension(:), allocatable  :: yab
    real(sp)                             :: var_ab, var_b
    real(sp)                             :: f0_ab, f0_b2

    nsets = size(yc,1)
    npara = size(yc,2)

    allocate(yab(2*nsets))
    yab(      1:  nsets) = ya
    yab(nsets+1:2*nsets) = yb

    var_ab = variance(yab)
    var_b  = variance(yb)

    f0_ab  = dot_product(ya,yb) / real(nsets,sp)
    f0_b2  = (sum(yb) / real(nsets,sp) )**2

    if ( var_ab .gt. tiny(1.0_sp) ) then
       ! model outputs are different (usual case)
       do ii=1, npara
          si(ii)  =          ( dot_product(ya(:),yc(:,ii)) / real(nsets,sp) - f0_ab) / var_ab
          sti(ii) = 1.0_sp - ( dot_product(yb(:),yc(:,ii)) / real(nsets,sp) - f0_b2) / var_b
       end do
    else
       ! model outputs are equal (should never happen)
       si(:)  = 0.0_sp
       sti(:) = 0.0_sp
    end if

    if ( present(smean) ) then
       smean(:,1) = si(:)
       smean(:,2) = sti(:)
    end if

    if ( present(wmean) ) then
       wmean(:,1) = si(:)
       wmean(:,2) = sti(:)
    end if

    deallocate(yab)

  end subroutine sobol_index_0d_sp

  subroutine sobol_index_1d_dp(ya, yb, yc, si, sti, smean, wmean)

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
    real(dp), dimension(:), allocatable  :: yab
    real(dp)                             :: var_ab, var_b
    real(dp)                             :: f0_ab, f0_b2
    ! for averaging
    real(dp), dimension(:), allocatable  :: si_num, sti_num
    real(dp)                             :: si_denom, sti_denom
    real(dp)                             :: si_num_ii, sti_num_ii

    nsets = size(yc,1)
    npara = size(yc,2)
    nout  = size(yc,3)

    si(:,:)  = 0.0_dp
    sti(:,:) = 0.0_dp

    allocate(yab(2*nsets))

    allocate(si_num(npara))
    allocate(sti_num(npara))

    si_num(:)  = 0.0_dp
    sti_num(:) = 0.0_dp
    si_denom   = 0.0_dp
    sti_denom  = 0.0_dp

    do iout=1, nout
       yab(      1:  nsets) = ya(:,iout)
       yab(nsets+1:2*nsets) = yb(:,iout)
       
       f0_ab  = dot_product(ya(:,iout),yb(:,iout)) / real(nsets,dp)
       f0_b2  = (sum(yb(:,iout)) / real(nsets,dp) )**2

       var_ab = variance(yab)
       var_b  = variance(yb(:,iout))

       si_denom   = si_denom  + var_ab
       sti_denom  = sti_denom + var_b
       
       if ( var_ab .gt. tiny(1.0_dp) ) then
          ! model outputs are different (usual case)
          do ii=1, npara
             si_num_ii   = ( dot_product(ya(:,iout),yc(:,ii,iout)) / real(nsets,dp) - f0_ab )
             sti_num_ii  = ( dot_product(yb(:,iout),yc(:,ii,iout)) / real(nsets,dp) - f0_b2 )

             si_num(ii)  = si_num(ii)  + si_num_ii
             sti_num(ii) = sti_num(ii) + sti_num_ii
             
             si(ii,iout)  =          si_num_ii  / var_ab
             sti(ii,iout) = 1.0_dp - sti_num_ii / var_b
          end do
       end if
    end do

    if (present(smean)) then
       smean = 0.0_dp
       smean(:,1) = sum(si, dim=2) / real(nout,dp)
       smean(:,2) = sum(sti,dim=2) / real(nout,dp)
    end if

    if (present(wmean)) then
       wmean = 0.0_dp
       if ( si_denom .gt. tiny(1.0_dp) ) then
          wmean(:,1) =          si_num(:)  / si_denom
       end if
       if ( sti_denom .gt. tiny(1.0_dp) ) then
          wmean(:,2) = 1.0_dp - sti_num(:) / sti_denom
       end if
    end if
       
    deallocate(yab, si_num, sti_num)

  end subroutine sobol_index_1d_dp

  subroutine sobol_index_1d_sp(ya, yb, yc, si, sti, smean, wmean)

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
    real(sp), dimension(:), allocatable  :: yab
    real(sp)                             :: var_ab, var_b
    real(sp)                             :: f0_ab, f0_b2
    ! for averaging
    real(sp), dimension(:), allocatable  :: si_num, sti_num
    real(sp)                             :: si_denom, sti_denom
    real(sp)                             :: si_num_ii, sti_num_ii

    nsets = size(yc,1)
    npara = size(yc,2)
    nout  = size(yc,3)

    si(:,:)  = 0.0_sp
    sti(:,:) = 0.0_sp

    allocate(yab(2*nsets))

    allocate(si_num(npara))
    allocate(sti_num(npara))

    si_num(:)  = 0.0_sp
    sti_num(:) = 0.0_sp
    si_denom   = 0.0_sp
    sti_denom  = 0.0_sp

    do iout=1, nout
       yab(      1:  nsets) = ya(:,iout)
       yab(nsets+1:2*nsets) = yb(:,iout)
       
       f0_ab  = dot_product(ya(:,iout),yb(:,iout)) / real(nsets,sp)
       f0_b2  = (sum(yb(:,iout)) / real(nsets,sp) )**2

       var_ab = variance(yab)
       var_b  = variance(yb(:,iout))

       si_denom   = si_denom  + var_ab
       sti_denom  = sti_denom + var_b
       
       if ( var_ab .gt. tiny(1.0_sp) ) then
          ! model outputs are different (usual case)
          do ii=1, npara
             si_num_ii   = ( dot_product(ya(:,iout),yc(:,ii,iout)) / real(nsets,sp) - f0_ab )
             sti_num_ii  = ( dot_product(yb(:,iout),yc(:,ii,iout)) / real(nsets,sp) - f0_b2 )

             si_num(ii)  = si_num(ii)  + si_num_ii
             sti_num(ii) = sti_num(ii) + sti_num_ii
             
             si(ii,iout)  =          si_num_ii  / var_ab
             sti(ii,iout) = 1.0_sp - sti_num_ii / var_b
          end do
       end if
    end do

    if (present(smean)) then
       smean = 0.0_sp
       smean(:,1) = sum(si, dim=2) / real(nout,sp)
       smean(:,2) = sum(sti,dim=2) / real(nout,sp)
    end if

    if (present(wmean)) then
       wmean = 0.0_sp
       if ( si_denom .gt. tiny(1.0_sp) ) then
          wmean(:,1) =          si_num(:)  / si_denom
       end if
       if ( sti_denom .gt. tiny(1.0_sp) ) then
          wmean(:,2) = 1.0_sp - sti_num(:) / sti_denom
       end if
    end if
       
    deallocate(yab, si_num, sti_num)

  end subroutine sobol_index_1d_sp

END MODULE mo_sobol_index

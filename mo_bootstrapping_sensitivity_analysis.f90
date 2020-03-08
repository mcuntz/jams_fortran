!> \file mo_bootstrapping_sensitivity_analysis.f90

!> \brief Bootstrapping for sensitivity analysis.

!> \details This module generates new data sets using the bootstrapping method
!>          and calculates the sobol index for each data set.

!> \authors Leonie Bruckert
!> \date Oct 2014

MODULE mo_bootstrapping_sensitivity_analysis

  ! Written  Leonie Bruckert, Oct 2014

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2014 Leonie Bruckert
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

  USE mo_kind,          ONLY: i4, i8, sp, dp
  USE mo_xor4096,       ONLY: get_timeseed, xor4096, n_save_state
  USE mo_xor4096_apps,  ONLY: xor4096_range
  USE mo_sobol_index,   ONLY: sobol_index

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: bootstrap_si  ! generates new data sets and calculates the Sobol index for each set

  ! ------------------------------------------------------------------

  !     NAME
  !>        bootstrap_si

  !     PURPOSE
  !>        Artificially generates new data sets by randomly sampling with replacement and
  !>        then calculates the Sobol index for each new data set to estimate parameter sensitivity.
  !

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp)      :: ya(:)/ya(:,:)"       Model output for parameter sets A \n
  !>                                                               dim_1 = number of parameter sets \n
  !>                                                               dim_2 = number of model outputs, e.g. time steps
  !>
  !>        \param[in] "real(sp/dp)      :: yb(:)/yb(:,:)"       Model output for parameter sets B \n
  !>                                                               dim_1 = number of parameter sets \n
  !>                                                               dim_2 = number of model outputs, e.g. time steps
  !>
  !>        \param[in] "real(sp/dp)      :: yc(:,:)/yc(:,:,:)"   Model output for parameter sets C(i), i=1,npara \n
  !>                                                               dim_1 = number of parameter sets \n
  !>                                                               dim_2 = number of parameters \n
  !>                                                               dim_3 = number of model outputs, e.g. time steps
  !>
  !>        \param[in] "integer(i4)      :: n"                   Number of new data sets

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(sp/dp)      :: si(:,:)/si(:,:,:)"    Sobol index - Main effect \n
  !>                                                                 dim_1 = number of new data sets
  !>                                                                 dim_2 = number of parameters \n
  !>                                                                 dim_3 = number of model outputs, e.g. time steps \n
  !>
  !>        \param[out] "real(sp/dp)      :: sti(:,:)/sti(:,:,:)"  Sobol index - Main effect \n
  !>                                                                 dim_1 = number of new data sets
  !>                                                                 dim_2 = number of parameters \n
  !>                                                                 dim_3 = number of model outputs, e.g. time steps \n


  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4),    optional :: method"        Method for calculating the Sobol index \n
  !>                                                                 Default = 9
  !>        \param[in] "integer(i4/i8), optional :: seed"          Seed != 0

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       None
  !
  !     RESTRICTIONS
  !>       \note Input values must be floating points.
  !
  !     EXAMPLE
  !>        Model with 3 parameters:
  !>        ya    = (/ 1., 2., 3., 4. /)
  !>        yb    = (/ 3., 1., 3., 6. /)
  !>        yc(:,1)= (/ 1., 2.5, 2., 3. /)
  !>        yc(:,2)= (/ 4., 1., 4. 7. /)
  !>        yc(:,3)= (/ 2., 2., 3., 5. /)
  !>        n     = 5
  !>        call bootstrap_si(ya, yb, yc, n, si, sti)
  !>        -> see also example in test directory

  !     LITERATURE
  !

  !     HISTORY
  !>        \author Leonie Bruckert
  !>        \date Oct 2014

  INTERFACE bootstrap_si
     MODULE PROCEDURE bootstrap_si_0d_dp, bootstrap_si_0d_sp, bootstrap_si_1d_dp, bootstrap_si_1d_sp
  END INTERFACE bootstrap_si
  
  ! ------------------------------------------------------------------

CONTAINS

  SUBROUTINE bootstrap_si_0d_dp(ya,yb,yc,n,si,sti,method,seed)

    IMPLICIT NONE

    real(dp),     dimension(:),            intent(in)             :: ya     ! model output for parameter set A
    real(dp),     dimension(:),            intent(in)             :: yb     ! model output for parameter set B
    real(dp),     dimension(:,:),          intent(in)             :: yc     ! model output for parameter set C
    integer(i4),                           intent(in)             :: n      ! number of new data sets

    real(dp),     dimension(n,size(yc,2)), intent(out)            :: si     ! Sobol index (main effect)
    real(dp),     dimension(n,size(yc,2)), intent(out)            :: sti    ! Sobol index (total effect)

    integer(i4),                           intent(in),  optional  :: method ! method for calculating the sobol index
    integer(i8),                           intent(in),  optional  :: seed   ! seed

    ! local variables
    real(dp),     dimension(n)                         :: r
    integer(i8),  dimension(size(ya),n)                :: rn

    integer(i8),  dimension(n)                         :: seed_tmp
    integer(i8),  dimension(:,:), allocatable          :: save_state

    real(dp),     dimension(size(ya))                  :: ya_new, yb_new
    real(dp),     dimension(size(ya),size(yc,2))       :: yc_new

    integer(i4)                                        :: meth
    integer(i4)                                        :: i, j, k


    allocate(save_state(n,n_save_state))

    if (present(seed)) then
       seed_tmp(1)=seed
       do i=2,n
          seed_tmp(i)=seed_tmp(i-1)+1000_i8
       end do
    else
       call get_timeseed(seed_tmp)
    end if

    call xor4096(seed_tmp, r, save_state=save_state)
    call xor4096_range((/1_i8,int(size(ya),i8)/),rn, save_state=save_state) ! the numbers of a column represent the indices
    ! of the chosen parameter sets
    if (present(method)) then
       meth=method
    else
       meth=9
    end if

    loop1: do k=1,n

       ! generate new data set
       do j=1,size(ya)
          ya_new(j)=ya(rn(j,k))
          yb_new(j)=yb(rn(j,k))
       end do

       do i=1,size(yc,2)
          do j=1,size(ya)
             yc_new(j,i)=yc(rn(j,k),i)
          end do
       end do

       call sobol_index(ya_new, yb_new, yc_new, si(k,:), sti(k,:), method=meth)

    end do loop1

    deallocate(save_state)

  END SUBROUTINE bootstrap_si_0d_dp

  ! ------------------------------------------------------------------

  SUBROUTINE bootstrap_si_0d_sp(ya,yb,yc,n,si,sti,method,seed)

    IMPLICIT NONE

    real(sp),     dimension(:),            intent(in)             :: ya     ! model output for parameter set A
    real(sp),     dimension(:),            intent(in)             :: yb     ! model output for parameter set B
    real(sp),     dimension(:,:),          intent(in)             :: yc     ! model output for parameter set C
    integer(i4),                           intent(in)             :: n      ! number of new data sets

    real(sp),     dimension(n,size(yc,2)), intent(out)            :: si     ! Sobol index (main effect)
    real(sp),     dimension(n,size(yc,2)), intent(out)            :: sti    ! Sobol index (total effect)

    integer(i4),                           intent(in),  optional  :: method ! method for calculating the sobol index
    integer(i4),                           intent(in),  optional  :: seed   ! seed

    ! local variables
    real(sp),     dimension(n)                         :: r
    integer(i4),  dimension(size(ya),n)                :: rn

    integer(i4),  dimension(n)                         :: seed_tmp
    integer(i4),  dimension(:,:), allocatable          :: save_state

    real(sp),     dimension(size(ya))                  :: ya_new, yb_new
    real(sp),     dimension(size(ya),size(yc,2))       :: yc_new

    integer(i4)                                        :: meth
    integer(i4)                                        :: i, j, k


    allocate(save_state(n,n_save_state))

    if (present(seed)) then
       seed_tmp(1)=seed
       do i=2,n
          seed_tmp(i)=seed_tmp(i-1)+1000_i4
       end do
    else
       call get_timeseed(seed_tmp)
    end if

    call xor4096(seed_tmp, r, save_state=save_state)
    call xor4096_range((/1_i4,int(size(ya),i4)/),rn, save_state=save_state) ! the numbers of a column represent the indices
    ! of the chosen parameter sets
    if (present(method)) then
       meth=method
    else
       meth=9
    end if

    loop1: do k=1,n

       ! generate new data set
       do j=1,size(ya)
          ya_new(j)=ya(rn(j,k))
          yb_new(j)=yb(rn(j,k))
       end do

       do i=1,size(yc,2)
          do j=1,size(ya)
             yc_new(j,i)=yc(rn(j,k),i)
          end do
       end do

       call sobol_index(ya_new, yb_new, yc_new, si(k,:), sti(k,:), method=meth)

    end do loop1

    deallocate(save_state)

  END SUBROUTINE bootstrap_si_0d_sp

  ! ------------------------------------------------------------------

  SUBROUTINE bootstrap_si_1d_dp(ya,yb,yc,n,si,sti,method,seed)

    IMPLICIT NONE

    real(dp),     dimension(:,:),                     intent(in)             :: ya     ! model output for parameter set A
    real(dp),     dimension(:,:),                     intent(in)             :: yb     ! model output for parameter set B
    real(dp),     dimension(:,:,:),                   intent(in)             :: yc     ! model output for parameter set C
    integer(i4),                                      intent(in)             :: n      ! number of new data sets

    real(dp),     dimension(n,size(yc,2),size(ya,2)), intent(out)            :: si     ! Sobol index (main effect)
    real(dp),     dimension(n,size(yc,2),size(ya,2)), intent(out)            :: sti    ! Sobol index (total effect)

    integer(i4),                                      intent(in),  optional  :: method ! method for calculating the sobol index
    integer(i8),                                      intent(in),  optional  :: seed   ! seed

    ! local variables
    real(dp),     dimension(n)                                 :: r
    integer(i8),  dimension(size(ya,1),n)                      :: rn

    integer(i8),  dimension(n)                                 :: seed_tmp
    integer(i8),  dimension(:,:), allocatable                  :: save_state

    real(dp),     dimension(size(ya,1),size(ya,2))             :: ya_new, yb_new
    real(dp),     dimension(size(ya,1),size(yc,2),size(ya,2))  :: yc_new

    integer(i4)                                                :: meth
    integer(i4)                                                :: i, j, k


    allocate(save_state(n,n_save_state))

    if (present(seed)) then
       seed_tmp(1)=seed
       do i=2,n
          seed_tmp(i)=seed_tmp(i-1)+1000_i8
       end do
    else
       call get_timeseed(seed_tmp)
    end if

    call xor4096(seed_tmp, r, save_state=save_state)
    call xor4096_range((/1_i8,int(size(ya,1),i8)/),rn, save_state=save_state) ! the numbers of a column represent the indices
    ! of the chosen parameter sets
    if (present(method)) then
       meth=method
    else
       meth=9
    end if

    loop1: do k=1,n

       ! generate new data set
       do j=1,size(ya,1)
          ya_new(j,:)=ya(rn(j,k),:)
          yb_new(j,:)=yb(rn(j,k),:)
       end do

       do i=1,size(yc,2)
          do j=1,size(ya,1)
             yc_new(j,i,:)=yc(rn(j,k),i,:)
          end do
       end do

       call sobol_index(ya_new, yb_new, yc_new, si(k,:,:), sti(k,:,:), method=meth)

    end do loop1

    deallocate(save_state)

  END SUBROUTINE bootstrap_si_1d_dp

  ! ------------------------------------------------------------------

  SUBROUTINE bootstrap_si_1d_sp(ya,yb,yc,n,si,sti,method,seed)

    IMPLICIT NONE

    real(sp),     dimension(:,:),                     intent(in)             :: ya     ! model output for parameter set A
    real(sp),     dimension(:,:),                     intent(in)             :: yb     ! model output for parameter set B
    real(sp),     dimension(:,:,:),                   intent(in)             :: yc     ! model output for parameter set C
    integer(i4),                                      intent(in)             :: n      ! number of new data sets

    real(sp),     dimension(n,size(yc,2),size(ya,2)), intent(out)            :: si     ! Sobol index (main effect)
    real(sp),     dimension(n,size(yc,2),size(ya,2)), intent(out)            :: sti    ! Sobol index (total effect)

    integer(i4),                                      intent(in),  optional  :: method ! method for calculating the sobol index
    integer(i4),                                      intent(in),  optional  :: seed   ! seed

    ! local variables
    real(sp),     dimension(n)                                 :: r
    integer(i4),  dimension(size(ya,1),n)                      :: rn

    integer(i4),  dimension(n)                                 :: seed_tmp
    integer(i4),  dimension(:,:), allocatable                  :: save_state

    real(sp),     dimension(size(ya,1),size(ya,2))             :: ya_new, yb_new
    real(sp),     dimension(size(ya,1),size(yc,2),size(ya,2))  :: yc_new

    integer(i4)                                                :: meth
    integer(i4)                                                :: i, j, k


    allocate(save_state(n,n_save_state))

    if (present(seed)) then
       seed_tmp(1)=seed
       do i=2,n
          seed_tmp(i)=seed_tmp(i-1)+1000_i4
       end do
    else
       call get_timeseed(seed_tmp)
    end if

    call xor4096(seed_tmp, r, save_state=save_state)
    call xor4096_range((/1_i4,int(size(ya,1),i4)/),rn, save_state=save_state) ! the numbers of a column represent the indices
    ! of the chosen parameter sets
    if (present(method)) then
       meth=method
    else
       meth=9
    end if

    loop1: do k=1,n

       ! generate new data set
       do j=1,size(ya,1)
          ya_new(j,:)=ya(rn(j,k),:)
          yb_new(j,:)=yb(rn(j,k),:)
       end do

       do i=1,size(yc,2)
          do j=1,size(ya,1)
             yc_new(j,i,:)=yc(rn(j,k),i,:)
          end do
       end do

       call sobol_index(ya_new, yb_new, yc_new, si(k,:,:), sti(k,:,:), method=meth)

    end do loop1

    deallocate(save_state)

  END SUBROUTINE bootstrap_si_1d_sp

  ! ------------------------------------------------------------------

END MODULE MO_BOOTSTRAPPING_SENSITIVITY_ANALYSIS

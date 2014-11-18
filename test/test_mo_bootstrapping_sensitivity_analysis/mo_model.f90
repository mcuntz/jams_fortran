MODULE mo_model

  USE mo_kind,   only: dp, i4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: par_range
  PUBLIC :: model

  INTERFACE model
     MODULE PROCEDURE model_0d, model_1d
  END INTERFACE model

  integer(i4) :: npar = 4_i4

CONTAINS
!----------------------------------------------------------

  FUNCTION par_range()

    IMPLICIT NONE

    real(dp),   dimension(npar,2)   :: par_range

    par_range(1,:) = (/2.0_dp, 6.0_dp/)
    par_range(2,:) = (/1.0_dp, 2.0_dp/)
    par_range(3,:) = (/0.0_dp, 1.0_dp/)
    par_range(4,:) = (/-2.0_dp, -1.0_dp/)

  END FUNCTION par_range
!----------------------------------------------------------

  FUNCTION model_0d(p,t)

    IMPLICIT NONE

    real(dp),    dimension(:,:),    intent(in)      :: p              ! Parameter sets
    real(dp),                       intent(in)      :: t              ! Time

    real(dp),    dimension(size(p,1))               :: model_0d       ! Model output

    integer(i4)                                     :: i               

    do i = 1, size(p,1)              
          model_0d(i) = p(i,1)*t**2 + p(i,2)*t + p(i,4)
    end do

  END FUNCTION model_0d
!----------------------------------------------------------

  FUNCTION model_1d(p,t)

    IMPLICIT NONE

    real(dp),    dimension(:,:),        intent(in)    :: p              ! Parameter sets
    real(dp),    dimension(:),          intent(in)    :: t              ! Time series

    real(dp),    dimension(size(p,1),size(t))         :: model_1d       ! Model output

    integer(i4)                                       :: i, j               

    do i = 1, size(p,1)
       do j = 1, size(t)
          model_1d(i,j) = p(i,1)*t(j)**2 + p(i,2)*t(j) + p(i,4)
       end do
    end do

  END FUNCTION model_1d
!----------------------------------------------------------

END MODULE mo_model

Module mo_cost

use mo_kind, only: sp, dp

IMPLICIT NONE

PRIVATE

PUBLIC :: cost_sp, cost_dp

CONTAINS

FUNCTION cost_sp (paraset)

implicit none

REAL(SP), DIMENSION(:), INTENT(IN) :: paraset
REAL(SP)                           :: cost_sp
REAL(SP), DIMENSION(6,2)           :: meas
REAL(SP), DIMENSION(6)             :: calc

! function: f(x) = ax^3 + bx^2 + cx + d
! measurements: (0.5,5.73), (1.0, 20.70), (1.5, 1.80), (2.0, 9.73), (2.5, 35.01), (3.0, 115.43)
! --> a=1.0, b=20.0, c=0.2, d=0.5

meas(:,1) = (/0.5_sp, 1.0_sp, 1.5_sp, 2.0_sp, 2.5_sp, 3.0_sp/)
meas(:,2) = (/5.7250_sp, 21.7000_sp, 49.1750_sp, 88.9000_sp, 141.6250_sp, 208.1000_sp/)

calc(:) = paraset(1)*meas(:,1)**3+paraset(2)*meas(:,1)**2+paraset(3)*meas(:,1)+paraset(4)

! MSE  Mean Square Error
! cost_sp = dot_product(meas(:,2)-calc(:),meas(:,2)-calc(:)))/size(meas,1)

! RMSE  Root Mean Square Error
! cost_sp = sqrt(dot_product(meas(:,2)-calc(:),meas(:,2)-calc(:))/size(meas,1))

! MAE  Mean Absolute Error
cost_sp = sum(abs( meas(:,2)-calc(:) ))/size(meas,1)

RETURN
END FUNCTION cost_sp

FUNCTION cost_dp (paraset)

implicit none

REAL(DP), DIMENSION(:), INTENT(IN) :: paraset
REAL(DP)                           :: cost_dp
REAL(DP), DIMENSION(6,2)           :: meas
REAL(DP), DIMENSION(6)             :: calc

! function: f(x) = ax^3 + bx^2 + cx + d
! measurements: (0.5,5.73), (1.0, 20.70), (1.5, 1.80), (2.0, 9.73), (2.5, 35.01), (3.0, 115.43)
! --> a=1.0, b=20.0, c=0.2, d=0.5

meas(:,1) = (/0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp/)
meas(:,2) = (/5.7250_dp, 21.7000_dp, 49.1750_dp, 88.9000_dp, 141.6250_dp, 208.1000_dp/)

calc(:) = paraset(1)*meas(:,1)**3+paraset(2)*meas(:,1)**2+paraset(3)*meas(:,1)+paraset(4)

! MSE  Mean Square Error
! cost_dp = dot_product(meas(:,2)-calc(:),meas(:,2)-calc(:)))/size(meas,1)

! RMSE  Root Mean Square Error
! cost_dp = sqrt(dot_product(meas(:,2)-calc(:),meas(:,2)-calc(:))/size(meas,1))

! MAE  Mean Absolute Error
cost_dp = sum(abs( meas(:,2)-calc(:) ))/size(meas,1)

RETURN
END FUNCTION cost_dp

END MODULE mo_cost

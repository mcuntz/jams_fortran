module mo_likelihood

  USE mo_kind, ONLY: i4, dp

  Implicit NONE

  PRIVATE

  PUBLIC :: data        
  PUBLIC :: model
  PUBLIC :: likelihood_dp, loglikelihood_dp
  PUBLIC :: setmeas
  PUBLIC :: stdev_data_dp

  INTERFACE data
     MODULE PROCEDURE data_dp
  END INTERFACE data

  INTERFACE model
     MODULE PROCEDURE model_dp
  END INTERFACE model

  INTERFACE setmeas
     MODULE PROCEDURE setmeas_dp
  END INTERFACE setmeas

  REAL(DP),DIMENSION(100,2)  :: meas

  ! ------------------------------------------------------------------

CONTAINS

  ! -------------------------------
  ! A Likelihood function
  ! -------------------------------
  function likelihood_dp(paraset,stddev)
    REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset            ! parameter set
    REAL(DP),               INTENT(IN)  :: stddev             ! standard deviation of data
    REAL(DP)                            :: likelihood_dp

    ! local
    REAL(DP), DIMENSION(size(meas,1))   :: errors

    errors  = model(paraset)-data()
    likelihood_dp = sum( errors(:) * errors(:) / stddev**2 )
    likelihood_dp = exp(-0.5_dp * likelihood_dp)

  end function likelihood_dp

  ! -------------------------------
  ! A Log-Likelihood function
  ! -------------------------------
  function loglikelihood_dp(paraset,stddev)
    REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset            ! parameter set
    REAL(DP),               INTENT(IN)  :: stddev             ! standard deviation of data
    REAL(DP)                            :: loglikelihood_dp

    ! local
    REAL(DP), DIMENSION(size(meas,1))   :: errors

    errors  = model(paraset)-data()
    loglikelihood_dp = sum( errors(:) * errors(:) / stddev**2 )
    loglikelihood_dp = -0.5_dp * loglikelihood_dp

  end function loglikelihood_dp

  ! -------------------------------
  ! Standard Deviation of the data points
  ! -------------------------------
  function stdev_data_dp(paraset)
    use mo_moment, only: mean
    REAL(DP), DIMENSION(:), INTENT(IN)  :: paraset
    REAL(DP)                            :: stdev_data_dp
    
    ! local
    REAL(DP), DIMENSION(size(meas,1))   :: errors
    
    errors = model_dp(paraset)-data()
    stdev_data_dp  = sqrt(sum( (errors(:) - mean(errors)) * (errors(:) - mean(errors))) / real(size(meas,1)-1,dp))
    
  end function stdev_data_dp

  ! -------------------------------
  ! A Model: p1*x^2 + p2*x + p3
  ! -------------------------------
  function model_dp(paraset)
    use mo_kind, only: dp
    REAL(DP), DIMENSION(:), INTENT(IN)     :: paraset
    REAL(DP), DIMENSION(size(meas,1))      :: model_dp

    model_dp(:) = paraset(1) * meas(:,1) * meas(:,1) + paraset(2) * meas(:,1) + paraset(3)

  end function model_dp

  ! -------------------------------
  ! DATA
  ! -------------------------------
  function data_dp()
    use mo_kind
    REAL(DP), DIMENSION(size(meas,1))      :: data_dp

    data_dp = meas(:,2)

  end function data_dp

  subroutine setmeas_dp()

    integer(i4) :: i

    do i=1,100
       meas(i,1) = real(i,dp)
    end do

    meas(:,2) = (/ 5.49537_dp, 10.7835_dp, 17.6394_dp, 26.8661_dp, 36.9247_dp, 50.9517_dp, 66.2058_dp, &
         82.9703_dp, 101.26_dp, 123.076_dp, 145.457_dp, 171.078_dp, 198.349_dp, 227.23_dp, 257.922_dp, &
         290.098_dp, 325.775_dp, 362.724_dp, 402.669_dp, 442.461_dp, 486.122_dp, 531.193_dp, &
         577.931_dp, 627.091_dp, 678.096_dp, 731.364_dp, 786.039_dp, 843.531_dp, 903.126_dp, &
         963.037_dp, 1025.85_dp, 1091.85_dp, 1158.47_dp, 1226.65_dp, 1298.25_dp, 1370.87_dp, &
         1444.96_dp, 1522.6_dp, 1602.6_dp, 1684.38_dp, 1765.15_dp, 1850.74_dp, 1937.85_dp, 2027.41_dp, &
         2118.44_dp, 2210.62_dp, 2306.9_dp, 2403.27_dp, 2501.83_dp, 2602.96_dp, 2705.29_dp, &
         2811.44_dp, 2917.82_dp, 3027.09_dp, 3137.64_dp, 3250.86_dp, 3366.67_dp, 3482.56_dp, &
         3602.37_dp, 3722.55_dp, 3845.8_dp, 3970.15_dp, 4098.15_dp, 4227.27_dp, 4357.77_dp, &
         4491.49_dp, 4626.23_dp, 4762.95_dp, 4901.15_dp, 5042.95_dp, 5184.86_dp, 5330.4_dp, &
         5478.11_dp, 5626.97_dp, 5776.94_dp, 5931.15_dp, 6085.9_dp, 6242.7_dp, 6402.17_dp, 6563.46_dp, &
         6726.37_dp, 6890.34_dp, 7058.4_dp, 7227.17_dp, 7397.09_dp, 7570.77_dp, 7745.95_dp, &
         7923.72_dp, 8102.21_dp, 8282.98_dp, 8465.42_dp, 8651.37_dp, 8838.43_dp, 9027.57_dp, &
         9219.21_dp, 9411.24_dp, 9606.31_dp, 9802.88_dp, 10002.3_dp, 10202.6_dp /)
  end subroutine setmeas_dp

end module mo_likelihood

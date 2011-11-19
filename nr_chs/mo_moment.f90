MODULE mo_moment

  ! This module contains the routine for calculating optionally
  ! the 1st to 4th moment of a distribution plus some extras

  ! Usage:
  !   USE mo_moment, ONLY: moment
  !   call moment(arr, average=ave, variance=var, skewness=skew, curtosis=curt, &
  !               mean=mean, stddev=stddev, absdev=absdev)

  ! Written March 2011, Matthias Cuntz
  !   - modified numerical recipes: moment

  !USE kinds, ONLY: sp, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: moment       ! Moments of an array

  INTERFACE moment
     MODULE PROCEDURE moment_sp, moment_dp
  END INTERFACE moment

  ! if no kind module
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)  
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  SUBROUTINE moment_sp(dat, average, variance, skewness, curtosis, mean, stddev, absdev)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN)  :: dat
    REAL(sp), OPTIONAL,     INTENT(OUT) :: average
    REAL(sp), OPTIONAL,     INTENT(OUT) :: variance
    REAL(sp), OPTIONAL,     INTENT(OUT) :: skewness
    REAL(sp), OPTIONAL,     INTENT(OUT) :: curtosis
    REAL(sp), OPTIONAL,     INTENT(OUT) :: mean
    REAL(sp), OPTIONAL,     INTENT(OUT) :: stddev
    REAL(sp), OPTIONAL,     INTENT(OUT) :: absdev

    REAL(sp) :: n

    REAL(sp) :: ep, ave, var
    REAL(sp), DIMENSION(size(dat)) :: p, s

    n = real(size(dat),sp)
    if (n <= 1.0_sp) stop 'moment_sp: n must be at least 2'
    ! Any optional argument
    if (.not. (present(average) .or. present(variance) .or. present(skewness) .or. &
         present(curtosis) .or. present(mean) .or. present(stddev) .or. present(absdev))) return
    ! Average
    ave  = sum(dat(:))/n
    if (present(average)) average = ave
    if (present(mean))    mean    = ave
    if (.not. (present(variance) .or. present(skewness) .or. &
         present(curtosis) .or. present(stddev) .or. present(absdev))) return
    ! Absolute deviation
    s(:) = dat(:)-ave
    if (present(absdev)) absdev = sum(abs(s(:)))/n
    ! Variance / Standard deviation
    if (.not. (present(variance) .or. present(skewness) .or. &
         present(curtosis) .or. present(stddev))) return
    ep   = sum(s(:))
    p(:) = s(:)*s(:)
    var = sum(p(:))
    var = (var-ep*ep/n)/(n-1.0_sp)
    if (present(variance)) variance = var
    if (present(stddev))    stddev   = sqrt(var)
    if (.not. (present(skewness) .or. present(curtosis))) return
    ! Skewness
    if (var == 0.0_sp) stop 'moment_sp: no skewness or kurtosis when zero variance'
    p(:) = p(:)*s(:)
    if (present(skewness)) then
       skewness = sum(p(:))
       skewness = skewness/(n*stddev*stddev*stddev)
    endif
    ! Curtosis
    if (present(curtosis)) then
       p(:) = p(:)*s(:)
       curtosis = sum(p(:))
       curtosis = curtosis/(n*variance*variance) - 3.0_sp
    end if
    
  END SUBROUTINE moment_sp


  SUBROUTINE moment_dp(dat, average, variance, skewness, curtosis, mean, stddev, absdev)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN)  :: dat
    REAL(dp), OPTIONAL,     INTENT(OUT) :: average
    REAL(dp), OPTIONAL,     INTENT(OUT) :: variance
    REAL(dp), OPTIONAL,     INTENT(OUT) :: skewness
    REAL(dp), OPTIONAL,     INTENT(OUT) :: curtosis
    REAL(dp), OPTIONAL,     INTENT(OUT) :: mean
    REAL(dp), OPTIONAL,     INTENT(OUT) :: stddev
    REAL(dp), OPTIONAL,     INTENT(OUT) :: absdev

    REAL(dp) :: n

    REAL(dp) :: ep, ave, var
    REAL(dp), DIMENSION(size(dat)) :: p, s

    n = real(size(dat),dp)
    if (n <= 1.0_dp) stop 'moment_dp: n must be at least 2'
    ! Any optional argument
    if (.not. (present(average) .or. present(variance) .or. present(skewness) .or. &
         present(curtosis) .or. present(mean) .or. present(stddev) .or. present(absdev))) return
    ! Average
    ave  = sum(dat(:))/n
    if (present(average)) average = ave
    if (present(mean))    mean    = ave
    if (.not. (present(variance) .or. present(skewness) .or. &
         present(curtosis) .or. present(stddev) .or. present(absdev))) return
    ! Absolute deviation
    s(:) = dat(:)-ave
    if (present(absdev)) absdev = sum(abs(s(:)))/n
    ! Variance / Standard deviation
    if (.not. (present(variance) .or. present(skewness) .or. &
         present(curtosis) .or. present(stddev))) return
    ep   = sum(s(:))
    p(:) = s(:)*s(:)
    var = sum(p(:))
    var = (var-ep*ep/n)/(n-1.0_dp)
    if (present(variance)) variance = var
    if (present(stddev))    stddev   = sqrt(var)
    if (.not. (present(skewness) .or. present(curtosis))) return
    ! Skewness
    if (var == 0.0_dp) stop 'moment_dp: no skewness or kurtosis when zero variance'
    p(:) = p(:)*s(:)
    if (present(skewness)) then
       skewness = sum(p(:))
       skewness = skewness/(n*stddev*stddev*stddev)
    endif
    ! Curtosis
    if (present(curtosis)) then
       p(:) = p(:)*s(:)
       curtosis = sum(p(:))
       curtosis = curtosis/(n*variance*variance) - 3.0_dp
    end if
    
  END SUBROUTINE moment_dp

END MODULE mo_moment

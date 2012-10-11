program test

  ! use mo_moment   
  USE mo_kind         , ONLY: i4, sp, dp
  USE mo_errorMeasures, ONLY: BIAS, NSE, SAE, MAE, SSE, MSE, RMSE
  !
  IMPLICIT NONE
  !
  INTEGER(i4)                 :: i, j, k
  REAL(dp), DIMENSION(54)     :: vec1, vec2
  REAL(dp), DIMENSION(9,6,8)  :: arr1, arr2
  LOGICAL,  DIMENSION(54)     :: maskvec
  LOGICAL,  DIMENSION(9,6,8)  :: mask
  LOGICAL                     :: isgood

  write(*,*) ''
  write(*,*) 'Test mo_errormeasures.f90'
  !
  ! read random mask - shape = (6,9,5)
  !open(unit=20, file='field_maskf90.dat', action='read', status='old')
  open(unit=20, file='../FORTRAN_chs_lib/test/test_mo_errormeasures/field_maskf90.dat', action='read', status='old')
  do k = 1, size(arr1, dim=3)
     do i = 1, size(arr1, dim=1)
        read(20,*) (mask(i,j,k), j=1, size(arr1, dim=2))
     end do
  end do
  close(20)
  ! read random number array - shape = (6,9,5)
  !open(unit=20, file='field.dat', action='read', status='old')
  open(unit=20, file='../FORTRAN_chs_lib/test/test_mo_errormeasures/field.dat', action='read', status='old')
  do k = 1, size(arr1, dim=3)
     do i = 1, size(arr1, dim=1)
        read(20,*) (arr1(i,j,k), j=1, size(arr1, dim=2))
     end do
     arr2(:,:,k) = arr1(:,:,k) + real(k, dp) / 5.0_dp
  end do
  close(20)
  !
  ! create example for 1D
  vec1     = reshape(arr1(:,:,1), (/size(arr1, dim=1)*size(arr1, dim=2)/))
  vec2     = reshape(arr1(:,:,2), (/size(arr1, dim=1)*size(arr1, dim=2)/))
  maskvec  = reshape(mask(:,:,1), (/size(mask, dim=1)*size(mask, dim=2)/))
  !
  isgood = .TRUE.
  !
  ! Test masked arrays (sp / dp, 1D / 2D / 3D)
  ! BIAS
  isgood = isgood .and. (anint(10000._sp*BIAS(real(vec1, sp), real(vec2, sp), mask=maskvec)) .EQ. 14353._sp)
  isgood = isgood .and. (anint(10000._dp*BIAS(vec1, vec2, mask=maskvec)) .EQ. 14353._dp)
  isgood = isgood .and. (anint(10000._sp*BIAS(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))) .EQ. &
           20813._sp)
  isgood = isgood .and. (anint(10000._dp*BIAS(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))) .EQ. 20813._dp)
  isgood = isgood .and. (anint(10000._sp*BIAS(real(arr1, sp), real(arr2, sp), mask=mask)) .EQ. 9097._sp) 
  isgood = isgood .and. (anint(10000._dp*BIAS(arr1, arr2, mask=mask)) .EQ. 9097._dp)
  ! MAE
  isgood = isgood .and. (anint(10000._sp*MAE(real(vec1, sp), real(vec2, sp), mask=maskvec)) .EQ. 16123._sp)
  isgood = isgood .and. (anint(10000._dp*MAE(vec1, vec2, mask=maskvec)) .EQ. 16123._dp)
  isgood = isgood .and. (anint(10000._sp*MAE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))) .EQ. &
           21134._sp)
  isgood = isgood .and. (anint(10000._dp*MAE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))) .EQ. 21134._dp)
  isgood = isgood .and. (anint(10000._sp*MAE(real(arr1, sp), real(arr2, sp), mask=mask)) .EQ. 9097._sp) 
  isgood = isgood .and. (anint(10000._dp*MAE(arr1, arr2, mask=mask)) .EQ. 9097._dp)
  ! MSE
  isgood = isgood .and. (anint(10000._sp*MSE(real(vec1, sp), real(vec2, sp), mask=maskvec)) .EQ. 42183._sp)
  isgood = isgood .and. (anint(10000._dp*MSE(vec1, vec2, mask=maskvec)) .EQ. 42183._dp)
  isgood = isgood .and. (anint(10000._sp*MSE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))) .EQ. &
           56920._sp)
  isgood = isgood .and. (anint(10000._dp*MSE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))) .EQ. 56920._dp)
  isgood = isgood .and. (anint(10000._sp*MSE(real(arr1, sp), real(arr2, sp), mask=mask)) .EQ. 10287._sp) 
  isgood = isgood .and. (anint(10000._dp*MSE(arr1, arr2, mask=mask)) .EQ. 10287._dp)
  ! NSE
  isgood = isgood .and. (anint(10000._sp*NSE(real(vec1, sp), real(vec2, sp), mask=maskvec)) .EQ. -23368._sp)
  isgood = isgood .and. (anint(10000._dp*NSE(vec1, vec2, mask=maskvec)) .EQ. -23368._dp)
  isgood = isgood .and. (anint(10000._sp*NSE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))) .EQ. &
           -35026._sp)
  isgood = isgood .and. (anint(10000._dp*NSE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))) .EQ. -35026._dp)
  isgood = isgood .and. (anint(10000._sp*NSE(real(arr1, sp), real(arr2, sp), mask=mask)) .EQ. 8164._sp) 
  isgood = isgood .and. (anint(10000._dp*NSE(arr1, arr2, mask=mask)) .EQ. 8164._dp)
  ! SAE
  isgood = isgood .and. (anint(10000._sp*SAE(real(vec1, sp), real(vec2, sp), mask=maskvec)) .EQ. 403063._sp)
  isgood = isgood .and. (anint(10000._dp*SAE(vec1, vec2, mask=maskvec)) .EQ. 403063._dp)
  isgood = isgood .and. (anint(10000._sp*SAE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))) .EQ. &
           528359._sp)
  isgood = isgood .and. (anint(10000._dp*SAE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))) .EQ. 528359._dp)
  isgood = isgood .and. (anint(10000._sp*SAE(real(arr1, sp), real(arr2, sp), mask=mask)) .EQ. 2055999._sp) 
  isgood = isgood .and. (anint(10000._dp*SAE(arr1, arr2, mask=mask)) .EQ. 2056000._dp)
  ! SSE
  isgood = isgood .and. (anint(10000._sp*SSE(real(vec1, sp), real(vec2, sp), mask=maskvec)) .EQ. 1054575._sp)
  isgood = isgood .and. (anint(10000._dp*SSE(vec1, vec2, mask=maskvec)) .EQ. 1054575._dp)
  isgood = isgood .and. (anint(10000._sp*SSE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))) .EQ. &
           1423003._sp)
  isgood = isgood .and. (anint(10000._dp*SSE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))) .EQ. 1423003._dp)
  isgood = isgood .and. (anint(10000._sp*SSE(real(arr1, sp), real(arr2, sp), mask=mask)) .EQ. 2324799._sp) 
  isgood = isgood .and. (anint(10000._dp*SSE(arr1, arr2, mask=mask)) .EQ. 2324800._dp)
  ! RMSE
  isgood = isgood .and. (anint(10000._sp*RMSE(real(vec1, sp), real(vec2, sp), mask=maskvec)) .EQ. 20538._sp)
  isgood = isgood .and. (anint(10000._dp*RMSE(vec1, vec2, mask=maskvec)) .EQ. 20538._dp)
  isgood = isgood .and. (anint(10000._sp*RMSE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))) .EQ. &
           23858._sp)
  isgood = isgood .and. (anint(10000._dp*RMSE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))) .EQ. 23858._dp)
  isgood = isgood .and. (anint(10000._sp*RMSE(real(arr1, sp), real(arr2, sp), mask=mask)) .EQ. 10142._sp) 
  isgood = isgood .and. (anint(10000._dp*RMSE(arr1, arr2, mask=mask)) .EQ. 10142._dp)
  !
  if (isgood) then
     write(*,*) 'mo_errormeasures o.k.'
  else
     write(*,*) '***ERROR: mo_errormeasures failed! ***'
  endif
  !
end program test

program test

  ! use mo_moment
  USE mo_kind         , ONLY: i4, sp, dp
  USE mo_errorMeasures, ONLY: BIAS, KGE, LNNSE, MAE, MSE, NSE, RMSE, SAE, SSE
  use mo_ansi_colors, only: color, c_red, c_green
  !
  IMPLICIT NONE
  !
  INTEGER(i4)                 :: i, j, k
  REAL(dp), DIMENSION(54)     :: vec1, vec2
  REAL(dp), DIMENSION(9,6,8)  :: arr1, arr2
  LOGICAL,  DIMENSION(54)     :: maskvec
  LOGICAL,  DIMENSION(9,6,8)  :: mask
  LOGICAL                     :: isgood
  real(dp)                    :: stmp
  real(dp)                    :: dtmp

  write(*,*) ''
  write(*,*) 'Test mo_errormeasures.f90'
  !
  ! read random mask - shape = (6,9,5)
  !open(unit=20, file='field_maskf90.dat', action='read', status='old')
  open(unit=20, file='field_maskf90.dat', action='read', status='old')
  do k = 1, size(arr1, dim=3)
     do i = 1, size(arr1, dim=1)
        read(20,*) (mask(i,j,k), j=1, size(arr1, dim=2))
     end do
  end do
  close(20)
  ! read random number array - shape = (6,9,5)
  !open(unit=20, file='field.dat', action='read', status='old')
  open(unit=20, file='field.dat', action='read', status='old')
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
  stmp = BIAS(real(vec1, sp), real(vec2, sp), mask=maskvec)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 14353)
  dtmp = BIAS(vec1, vec2, mask=maskvec)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 14353)
  stmp = BIAS(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 20813)
  dtmp = BIAS(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 20813)
  stmp = BIAS(real(arr1, sp), real(arr2, sp), mask=mask)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 9097)
  dtmp = BIAS(arr1, arr2, mask=mask)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 9097)
  ! KGE
  stmp = KGE(real(vec1, sp), real(vec2, sp), mask=maskvec)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 1783)
  dtmp = KGE(vec1, vec2, mask=maskvec)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 1783)
  stmp = KGE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. -4730)
  dtmp = KGE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. -4730)
  stmp = KGE(real(arr1, sp), real(arr2, sp), mask=mask)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 7631)
  dtmp = KGE(arr1, arr2, mask=mask)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 7631)
  ! MAE
  stmp = MAE(real(vec1, sp), real(vec2, sp), mask=maskvec)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 16123)
  dtmp = MAE(vec1, vec2, mask=maskvec)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 16123)
  stmp = MAE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. &
           21134)
  dtmp = MAE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 21134)
  stmp = MAE(real(arr1, sp), real(arr2, sp), mask=mask)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 9097)
  dtmp = MAE(arr1, arr2, mask=mask)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 9097)
  ! MSE
  stmp = MSE(real(vec1, sp), real(vec2, sp), mask=maskvec)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 42183)
  dtmp = MSE(vec1, vec2, mask=maskvec)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 42183)
  stmp = MSE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 56920)
  dtmp = MSE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 56920)
  stmp = MSE(real(arr1, sp), real(arr2, sp), mask=mask)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 10287)
  dtmp = MSE(arr1, arr2, mask=mask)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 10287)
  ! NSE
  stmp = NSE(real(vec1, sp), real(vec2, sp), mask=maskvec)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. -23368)
  dtmp = NSE(vec1, vec2, mask=maskvec)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. -23368)
  stmp = NSE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. -35026)
  dtmp = NSE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. -35026)
  stmp = NSE(real(arr1, sp), real(arr2, sp), mask=mask)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 8164)
  dtmp = NSE(arr1, arr2, mask=mask)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 8164)
  ! SAE
  stmp = SAE(real(vec1, sp), real(vec2, sp), mask=maskvec)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 403063)
  dtmp = SAE(vec1, vec2, mask=maskvec)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 403063)
  stmp = SAE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 528359)
  dtmp = SAE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 528359)
  stmp = SAE(real(arr1, sp), real(arr2, sp), mask=mask)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 2055999)
  dtmp = SAE(arr1, arr2, mask=mask)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 2056000)
  ! SSE
  stmp = SSE(real(vec1, sp), real(vec2, sp), mask=maskvec)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 1054575)
  dtmp = SSE(vec1, vec2, mask=maskvec)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 1054575)
  stmp = SSE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 1423003)
  dtmp = SSE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 1423003)
  stmp = SSE(real(arr1, sp), real(arr2, sp), mask=mask)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 2324799)
  dtmp = SSE(arr1, arr2, mask=mask)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 2324800)
  ! RMSE
  stmp = RMSE(real(vec1, sp), real(vec2, sp), mask=maskvec)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 20538)
  dtmp = RMSE(vec1, vec2, mask=maskvec)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 20538)
  stmp = RMSE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 23858)
  dtmp = RMSE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 23858)
  stmp = RMSE(real(arr1, sp), real(arr2, sp), mask=mask)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 10142)
  dtmp = RMSE(arr1, arr2, mask=mask)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 10142)
  !
  ! LNNSE
  ! Attention masks might have changed since inout
  !
  stmp = LNNSE(real(vec1, sp), real(vec2, sp), mask=maskvec)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. -8993)
  dtmp = LNNSE(vec1, vec2, mask=maskvec)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. -8993)
  stmp = LNNSE(real(arr1(:,:,1), sp), real(arr2(:,:,4), sp), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. -21656)
  dtmp = LNNSE(arr1(:,:,1), arr2(:,:,4), mask=mask(:,:,1))
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. -21656)
  stmp = LNNSE(real(arr1, sp), real(arr2, sp), mask=mask)
  isgood = isgood .and. (nint(10000._sp*stmp) .EQ. 8915)
  dtmp = LNNSE(arr1, arr2, mask=mask)
  isgood = isgood .and. (nint(10000._dp*dtmp) .EQ. 8915)
  !
  if (isgood) then
     write(*,*) 'mo_errormeasures ', color('o.k.', c_green)
  else
     write(*,*) 'mo_errormeasures ', color('failed', c_red)
  endif
  !
end program test

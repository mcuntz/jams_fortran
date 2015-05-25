! ------------------------------------------------------------------------------
!
! test program for mo_spatialsimilarity
!
! author: Matthias Zink
!
! created: May 2015
!
! ------------------------------------------------------------------------------

program test
  !
  use mo_kind,              only: sp, dp
  use mo_spatialsimilarity, only: PD, NNDV
  !
  implicit none
  !
  real(dp), dimension(3,3)  :: mat1, mat2 ! data arrays
  logical,  dimension(3,3)  :: masking    ! mask for excluding nodata values
  logical                   :: validity   ! number of valid cells
  logical                   :: isgood
  !
  isgood = .true.
  !
  mat1 = reshape((/12.0_dp,  4.0_dp, 15.0_dp, 17.0_dp, 10.0_dp,  2.0_dp,  1.0_dp, 11.0_dp, 20.0_dp/),(/3,3/))
  mat2 = reshape((/ 7.0_dp, 12.0_dp,  5.0_dp,  9.0_dp, 11.0_dp, 13.0_dp, 12.0_dp, 11.0_dp,  7.0_dp/),(/3,3/))


  masking = .TRUE.  
  ! sp without mask
  isgood  = isgood .and. (nint(100000.0_sp * NNDV(real(mat1(:,:),sp), real(mat2(:,:),sp)              )) .EQ. 29815)
  isgood  = isgood .and. (nint(100000.0_sp *   PD(real(mat1(:,:),sp), real(mat2(:,:),sp), mask=masking)) .EQ.  8148)

  ! dp without mask
  isgood  = isgood .and. (nint(100000.0_dp * NNDV(mat1(:,:), mat2(:,:),               valid=validity)) .EQ. 29815)
  isgood  = isgood .and. validity
  isgood  = isgood .and. (nint(100000.0_dp *   PD(mat1(:,:), mat2(:,:), mask=masking, valid=validity)) .EQ.  8148)
  isgood  = isgood .and. validity

  masking = reshape((/.TRUE., .TRUE., .FALSE., .FALSE., .TRUE., .TRUE., .TRUE., .TRUE., .FALSE./),(/3,3/))
  ! sp with mask
  isgood  = isgood .and. (nint(1000000.0_sp * NNDV(real(mat1(:,:),sp), real(mat2(:,:),sp), mask=masking)) .EQ. 300000)
  isgood  = isgood .and. (nint(1000000.0_sp *   PD(real(mat1(:,:),sp), real(mat2(:,:),sp), mask=masking)) .EQ.  55556)

  ! dp with mask
  masking = reshape((/.TRUE., .TRUE., .FALSE., .FALSE., .TRUE., .TRUE., .TRUE., .TRUE., .FALSE./),(/3,3/))
  isgood  = isgood .and. (nint(1000000.0_dp * NNDV(mat1(:,:), mat2(:,:), mask=masking, valid=validity)) .EQ. 300000)
  isgood  = isgood .and. validity
  isgood  = isgood .and. (nint(1000000.0_dp *   PD(mat1(:,:), mat2(:,:), mask=masking, valid=validity)) .EQ.  55556)
  isgood  = isgood .and. validity

  ! entire mask .false. - check valaidity argument
  masking = .FALSE.
  isgood  = isgood .and. (nint(1000000.0_sp * NNDV(mat1(:,:), mat2(:,:), valid=validity, mask=masking)) .EQ. 0)
  isgood  = isgood .and. ( .not. validity)
  isgood  = isgood .and. (nint(1000000.0_sp *   PD(mat1(:,:), mat2(:,:), valid=validity, mask=masking)) .EQ. 0)
  isgood  = isgood .and. ( .not. validity)
  
  if (isgood) then
     write(*,*) 'mo_spatialsimilarity o.k.'
  else
     write(*,*) 'mo_spatialsimilarity failed'
  endif
  !  
end program test

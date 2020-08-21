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
  use mo_ansi_colors,       only: color, c_red, c_green
  use mo_spatialsimilarity, only: PD, NNDV
  !
  implicit none
  !
  real(dp), dimension(3,3)  :: mat1, mat2 ! data arrays
  logical,  dimension(3,3)  :: masking    ! mask for excluding nodata values
  logical                   :: validity   ! number of valid cells
  logical                   :: isgood
  real(sp) :: stmp
  real(dp) :: dtmp
  !
  isgood = .true.
  !
  mat1 = reshape((/12.0_dp,  4.0_dp, 15.0_dp, 17.0_dp, 10.0_dp,  2.0_dp,  1.0_dp, 11.0_dp, 20.0_dp/),(/3,3/))
  mat2 = reshape((/ 7.0_dp, 12.0_dp,  5.0_dp,  9.0_dp, 11.0_dp, 13.0_dp, 12.0_dp, 11.0_dp,  7.0_dp/),(/3,3/))


  masking = .TRUE.
  ! sp without mask
  stmp = NNDV(real(mat1(:,:),sp), real(mat2(:,:),sp)              )
  isgood  = isgood .and. (nint(100000.0_sp * stmp) .EQ. 29815)
  stmp = PD(real(mat1(:,:),sp), real(mat2(:,:),sp), mask=masking)
  isgood  = isgood .and. (nint(100000.0_sp * stmp) .EQ.  8148)

  ! dp without mask
  dtmp = NNDV(mat1(:,:), mat2(:,:),               valid=validity)
  isgood  = isgood .and. (nint(100000.0_dp * dtmp) .EQ. 29815)
  isgood  = isgood .and. validity
  dtmp =   PD(mat1(:,:), mat2(:,:), mask=masking, valid=validity)
  isgood  = isgood .and. (nint(100000.0_dp * dtmp) .EQ.  8148)
  isgood  = isgood .and. validity

  masking = reshape((/.TRUE., .TRUE., .FALSE., .FALSE., .TRUE., .TRUE., .TRUE., .TRUE., .FALSE./),(/3,3/))
  ! sp with mask
  stmp = NNDV(real(mat1(:,:),sp), real(mat2(:,:),sp), mask=masking)
  isgood  = isgood .and. (nint(1000000.0_sp * stmp) .EQ. 300000)
  stmp =   PD(real(mat1(:,:),sp), real(mat2(:,:),sp), mask=masking)
  isgood  = isgood .and. (nint(1000000.0_sp * stmp) .EQ.  55556)

  ! dp with mask
  masking = reshape((/.TRUE., .TRUE., .FALSE., .FALSE., .TRUE., .TRUE., .TRUE., .TRUE., .FALSE./),(/3,3/))
  dtmp = NNDV(mat1(:,:), mat2(:,:), mask=masking, valid=validity)
  isgood  = isgood .and. (nint(1000000.0_dp * dtmp) .EQ. 300000)
  isgood  = isgood .and. validity
  dtmp =   PD(mat1(:,:), mat2(:,:), mask=masking, valid=validity)
  isgood  = isgood .and. (nint(1000000.0_dp * dtmp) .EQ.  55556)
  isgood  = isgood .and. validity

  ! entire mask .false. - check validity argument
  masking = .FALSE.
  dtmp = NNDV(mat1(:,:), mat2(:,:), valid=validity, mask=masking)
  isgood  = isgood .and. (nint(1000000.0_dp * dtmp) .EQ. 0)
  isgood  = isgood .and. ( .not. validity)
  dtmp =   PD(mat1(:,:), mat2(:,:), valid=validity, mask=masking)
  isgood  = isgood .and. (nint(1000000.0_dp * dtmp) .EQ. 0)
  isgood  = isgood .and. ( .not. validity)

  if (isgood) then
     write(*,*) 'mo_spatialsimilarity ', color('o.k.', c_green)
  else
     write(*,*) 'mo_spatialsimilarity ', color('failed', c_red)
  endif
  !
end program test

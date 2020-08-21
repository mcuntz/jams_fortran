! ------------------------------------------------------------------------------
!
! test program for mo_standard_score
!
! author: Matthias Zink
!
! created: May 2015
!
! ------------------------------------------------------------------------------

program test

  use mo_kind,           only: i4, sp, dp
  use mo_ansi_colors, only: color, c_red, c_green
  use mo_standard_score, only: standard_score, classified_standard_score

  implicit none

  real(dp),     dimension(9)    :: vec1       ! data array
  integer(i4),  dimension(9)    :: control    ! array with contral values
  integer(i4),  dimension(9)    :: classes    ! classes for classified_standard_score
  logical,      dimension(9)    :: masking    ! mask for excluding nodata values
  logical                       :: isgood
  real(sp), dimension(9) :: stmp
  real(dp), dimension(9) :: dtmp

  isgood = .true.

  vec1(:)    = (/12.0_dp,  4.0_dp, 15.0_dp, 17.0_dp, 10.0_dp,  2.0_dp,  1.0_dp, 11.0_dp, 20.0_dp/)

  masking = .TRUE.
  ! sp without mask
  ! standard_score
  control(:) = (/26518, -92813, 71267, 101100, -3315, -122645, -137562, 11602, 145849/)
  stmp = standard_score(real(vec1(:),sp)              )
  isgood  = isgood .and. all(nint(100000.0_sp * stmp) .EQ. control)
  if (.not. isgood) write(*,*) 'mo_standard_score 1'
  stmp = standard_score(real(vec1(:),sp), mask=masking)
  isgood  = isgood .and. all(nint(100000.0_sp * stmp) .EQ. control)
  if (.not. isgood) write(*,*) 'mo_standard_score 2'

  ! classified_standard_score
  classes =(/3,2,3,2,2,3,1,1,3/)
  control(:) = (/-3295, -97340, 36240, 102463, -5123, -135075, -70711, 70711, 102130/)
  stmp = classified_standard_score(real(vec1(:),sp), classes)
  isgood  = isgood .and. all(nint(100000.0_sp * stmp) .EQ. control)
  if (.not. isgood) write(*,*) 'mo_standard_score 3'

  masking(:) = (/.TRUE., .TRUE., .FALSE., .FALSE., .TRUE., .TRUE., .TRUE., .TRUE., .FALSE./)

  ! dp with mask
  ! standard_score
  control(:) = (/109170, -54585, 170578, 211517, 68231, -95524, -115993, 88701, 272925/)
  dtmp =            standard_score(     vec1(:),     mask=masking)
  isgood  = isgood .and. all(nint(100000.0_dp * dtmp) .EQ. control)
  if (.not. isgood) write(*,*) 'mo_standard_score 4'

  ! classified_standard_score
  control = (/70711, -70711, 0, 0, 70711, -70711, -70711, 70711, 0/)
  dtmp = classified_standard_score(vec1, classes, mask=masking)
  isgood  = isgood .and. all(nint(100000.0_dp * dtmp) .EQ. control)
  if (.not. isgood) write(*,*) 'mo_standard_score 5'

  if (isgood) then
     write(*,*) 'mo_standard_score ', color('o.k.', c_green)
  else
     write(*,*) 'mo_standard_score ', color('failed', c_red)
  endif
  !
end program test

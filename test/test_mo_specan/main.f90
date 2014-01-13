! testprogram using Wolfer’s sunspot data (Anderson 1971)
! calculating periodogram
! find best parameter for linear reservoir of logarithmic
! periodogram


PROGRAM sunspot_test

  USE mo_kind,   only: i4, dp
  USE mo_specan, only: periodogram, linres

  IMPLICIT NONE

  INTEGER(i4) :: i

  REAL(dp), DIMENSION(176)     :: sunspot
  REAL(dp), DIMENSION(176)     :: time  
  ! REAL(dp), DIMENSION(176)     :: sun_trans
  REAL(dp), DIMENSION(176/2,2) :: output
  REAL(dp), DIMENSION(2)       :: bestpara

  logical :: isgood

  isgood = .true.

  sunspot = (/809.0_dp, 834.0_dp, 477.0_dp, 478.0_dp, 307.0_dp, 122.0_dp, 96.0_dp, &
       102.0_dp, 324.0_dp, 476.0_dp, 540.0_dp, 629.0_dp, 859.0_dp, 612.0_dp,       &
       451.0_dp, 364.0_dp, 209.0_dp, 114.0_dp, 378.0_dp, 698.0_dp, 1061.0_dp,      &
       1008.0_dp, 816.0_dp, 665.0_dp, 348.0_dp, 306.0_dp, 70.0_dp, 198.0_dp,       &
       925.0_dp, 1544.0_dp, 1259.0_dp, 848.0_dp, 681.0_dp, 385.0_dp, 228.0_dp,     &
       102.0_dp, 241.0_dp, 829.0_dp, 1320.0_dp, 1309.0_dp, 1181.0_dp, 899.0_dp,    &
       666.0_dp, 600.0_dp, 469.0_dp, 410.0_dp, 213.0_dp, 160.0_dp, 64.0_dp,        &
       41.0_dp, 68.0_dp, 145.0_dp, 340.0_dp, 450.0_dp, 431.0_dp, 475.0_dp,         &
       422.0_dp, 281.0_dp, 101.0_dp, 81.0_dp, 25.0_dp, 0.0_dp, 14.0_dp,            &
       50.0_dp, 122.0_dp, 139.0_dp, 354.0_dp, 458.0_dp, 411.0_dp, 304.0_dp,        &
       239.0_dp, 157.0_dp, 66.0_dp, 40.0_dp, 18.0_dp, 85.0_dp, 166.0_dp,           &
       363.0_dp, 497.0_dp, 625.0_dp, 670.0_dp, 710.0_dp, 478.0_dp, 275.0_dp,       &
       85.0_dp, 132.0_dp, 569.0_dp, 1215.0_dp, 1383.0_dp, 1032.0_dp, 858.0_dp,     &
       632.0_dp, 368.0_dp, 242.0_dp, 107.0_dp, 150.0_dp, 401.0_dp, 615.0_dp,       &
       985.0_dp, 1243.0_dp, 959.0_dp, 665.0_dp, 645.0_dp, 542.0_dp, 390.0_dp,      &
       206.0_dp, 67.0_dp, 43.0_dp, 228.0_dp, 548.0_dp, 938.0_dp, 957.0_dp,         &
       772.0_dp, 591.0_dp, 440.0_dp, 470.0_dp, 305.0_dp, 163.0_dp, 73.0_dp,        &
       373.0_dp, 739.0_dp, 1391.0_dp, 1112.0_dp, 1017.0_dp, 663.0_dp, 447.0_dp,    &
       171.0_dp, 113.0_dp, 123.0_dp, 34.0_dp, 60.0_dp, 323.0_dp, 543.0_dp,         &
       597.0_dp, 637.0_dp, 635.0_dp, 522.0_dp, 254.0_dp, 131.0_dp, 68.0_dp,        &
       63.0_dp, 71.0_dp, 356.0_dp, 730.0_dp, 849.0_dp, 780.0_dp, 640.0_dp,         &
       418.0_dp, 262.0_dp, 267.0_dp, 121.0_dp, 95.0_dp, 27.0_dp, 50.0_dp,          &
       244.0_dp, 420.0_dp, 635.0_dp, 538.0_dp, 620.0_dp, 485.0_dp, 439.0_dp,       &
       186.0_dp, 57.0_dp, 36.0_dp, 14.0_dp, 96.0_dp, 474.0_dp, 571.0_dp,           &
       1039.0_dp, 806.0_dp, 636.0_dp, 376.0_dp, 261.0_dp, 142.0_dp, 58.0_dp,       &
       167.0_dp/)

  time(1)=1749.0_dp

  do i = 2, 176
     time(i) = time(i-1) + 1.0_dp
  end do

  call periodogram(sunspot,output)

  write(*,*) 'printing first five values of periodogram:'

  do i = 1, 5
     write(*,*) output(i,:)
  end do

  call linres( &
       output(:,2),   & ! data(:)        Time series
       0.00000005_dp, & ! a_min          minimum of first parameter
       0.1_dp,        & ! a_max          maximum of first parameter
       0._dp,         & ! tc_min         minimum of second parameter
       100._dp,       & ! tc_max         maximum of second parameter
       0.0000005_dp,  & ! initial_a      first initial parameter
       50._dp,        & ! initial_tc     second initial parameter
       bestpara )       ! best parameter set

  write(*,*) 'best parameter found:'
  write(*,*) '   alpha:', bestpara(1)
  write(*,*) '   tc:   ', bestpara(2)

  ! check if program was running properly
  isgood = isgood .and. ( 35268373_i4 .eq. int(output(1,2),i4) )
  isgood = isgood .and. ( 1589_i4     .eq. int(output(2,2),i4) )
  isgood = isgood .and. ( 1217716_i4  .eq. int(output(3,2),i4) )
  isgood = isgood .and. ( 538747_i4   .eq. int(output(4,2),i4) )
  isgood = isgood .and. ( 245925_i4   .eq. int(output(5,2),i4) )

  isgood = isgood .and. ( 6368_i4   .eq. int(bestpara(1)*10000000.0_dp,i4) )
  isgood = isgood .and. ( 376437_i4 .eq. int(bestpara(2)*10000.0_dp,i4) )

  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_specan: o.k.'
  else
     write(*,*) 'mo_specan: failed'
  end if

END PROGRAM sunspot_test

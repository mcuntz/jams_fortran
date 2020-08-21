program test

  !----------------------------------------------------------------------------
  ! MODULS AND VARIABLS
  !----------------------------------------------------------------------------

  USE mo_kind,         ONLY: sp, dp, i4
  use mo_ansi_colors,  only: color, c_red, c_green
  USE mo_pumpingtests, ONLY: thiem, theis, ext_thiem2d, ext_thiem3d, ext_theis2d, ext_theis3d
  use mo_utils,        ONLY: eq

  Implicit none

  REAL(sp) :: Reftest, hreftest, Qwtest, Ktest, Ltest, Stest, TGtest, &
       vartest, corrtest, anisotest!, proptest
  REAL(dp) :: Reftestd, hreftestd, Qwtestd, Ktestd, Ltestd, Stestd, &
       TGtestd, vartestd, corrtestd, anisotestd
  !, proptestd
  LOGICAL  :: isgood

  real(dp), dimension(3,3)   :: fpoints_dp, fpoints2_dp, fpoints3_dp
  real(sp), dimension(3,3)   :: fpoints_sp, fpoints2_sp, fpoints3_sp

  real(dp), dimension(9,2)   :: grid_dp
  real(sp), dimension(9,2)   :: grid_sp

  integer(i4) :: m,n

  real(sp) :: stmp
  real(dp) :: dtmp

  !----------------------------------------------------------------------------
  ! test-parameters
  !----------------------------------------------------------------------------

  Reftest         = 1.0_sp
  hreftest        = 0.0_sp
  Qwtest          = -0.001_sp
  Ktest           = 0.0001_sp
  Ltest           = 1.0_sp
  Stest           = 0.001_sp
  TGtest          = Ltest*Ktest
  vartest         = 1.0_sp
  !proptest        = 0.8_sp
  corrtest        = 10.0_sp
  anisotest       = 1.0_sp

  Reftestd        = 1.0_dp
  hreftestd       = 0.0_dp
  Qwtestd         = -0.001_dp
  Ktestd          = 0.0001_dp
  Ltestd          = 1.0_dp
  Stestd          = 0.001_dp
  TGtestd         = Ltestd*Ktestd
  vartestd        = 1.0_dp
  !proptestd       = 0.8_dp
  corrtestd       = 10.0_dp
  anisotestd      = 1.0_dp

  grid_sp         =   reshape((/  1._sp,      2._sp,      3._sp,&
       1._sp,      2._sp,      3._sp,&
       1._sp,      2._sp,      3._sp,&
       50.0_sp,    50.0_sp,    50.0_sp,&
       100.0_sp,   100.0_sp,   100.0_sp,&
       1000.0_sp,  1000.0_sp,  1000.0_sp /), (/9,2/))

  grid_dp         =   reshape((/  1._dp,      2._dp,      3._dp,&
       1._dp,      2._dp,      3._dp,&
       1._dp,      2._dp,      3._dp,&
       50.0_dp,    50.0_dp,    50.0_dp,&
       100.0_dp,   100.0_dp,   100.0_dp,&
       1000.0_dp,  1000.0_dp,  1000.0_dp /), (/9,2/))

  isgood          = .TRUE.

  !----------------------------------------------------------------------------
  ! tests
  !----------------------------------------------------------------------------

  write (*,*) ''
  write (*,*) 'Test mo_pumpingtests.f90'
  write (*,*) ''

  write (*,*) "THIEM's SOLUTION"
  write (*,*) '-------------------------------------------------------------------------'

  ! Thiemsp(0.2)         =  -2.56150007
  ! Thiemsp(0.3)         =  -1.91618240
  ! Thiemdp(0.2)         =  -2.5614999936338809
  ! Thiemdp(0.3)         =  -1.9161822315668402

  stmp = thiem(0.2_sp, (/Ktest /), (/Reftest , hreftest , Qwtest , Ltest /))
  isgood = isgood .and. (nint(stmp*1000)       .EQ. -2562)
  stmp = thiem(0.3_sp, (/Ktest /), (/Reftest , hreftest , Qwtest , Ltest /))
  isgood = isgood .and. (nint(stmp*1000)       .EQ. -1916)
  dtmp = thiem(0.2_dp, (/Ktestd/), (/Reftestd, hreftestd, Qwtestd, Ltestd/))
  isgood = isgood .and. (nint(dtmp*1000000)    .EQ. -2561500)
  dtmp = thiem(0.3_dp, (/Ktestd/), (/Reftestd, hreftestd, Qwtestd, Ltestd/))
  isgood = isgood .and. (nint(dtmp*1000000)    .EQ. -1916182)
  !sp
  write (*,*) 'Thiemsp(0.2) =', thiem(0.2_sp, (/Ktest /), (/Reftest , hreftest , Qwtest , Ltest /))
  write (*,*) 'Thiemsp(0.3) =', thiem(0.3_sp, (/Ktest /), (/Reftest , hreftest , Qwtest , Ltest /))
  !dp
  write (*,*) 'Thiemdp(0.2) =', thiem(0.2_dp, (/Ktestd/), (/Reftestd, hreftestd, Qwtestd, Ltestd/))
  write (*,*) 'Thiemdp(0.3) =', thiem(0.3_dp, (/Ktestd/), (/Reftestd, hreftestd, Qwtestd, Ltestd/))

  write (*,*) ' '

  write (*,*) "EXTENDED THIEM's SOLUTION IN 3D"
  write (*,*) '-------------------------------------------------------------------------'

  ! ext_Thiem3d_sp(0.2)  =  -4.19452763
  ! ext_Thiem3d_sp(0.3)  =  -3.13209820
  ! ext_Thiem3d_sp(0.2)  =  -1.55417931           optional
  ! ext_Thiem3d_sp(0.3)  =  -1.16274548           optional
  ! ext_Thiem3d_dp(0.2)  =  -4.1945268260359985
  ! ext_Thiem3d_dp(0.3)  =  -3.1320977426835279
  ! ext_Thiem3d_dp(0.2)  =  -1.5541790860603071   optional
  ! ext_Thiem3d_dp(0.3)  =  -1.1627453616249193   optional

  stmp = ext_thiem3d(0.2_sp, (/Ktest , vartest , corrtest , anisotest/), &
       (/Reftest , hreftest , Qwtest , Ltest  /))
  isgood = isgood .and. (nint(stmp*1000)          .EQ. -4195)
  stmp = ext_thiem3d(0.3_sp, (/Ktest , vartest , corrtest , anisotest/), &
       (/Reftest , hreftest , Qwtest , Ltest  /))
  isgood = isgood .and. (nint(stmp*1000)          .EQ. -3132)
  stmp = ext_thiem3d(0.2_sp, (/Ktest , vartest , corrtest , anisotest/), &
       (/Reftest , hreftest , Qwtest , Ltest  /), .false., 1.0_sp)
  isgood = isgood .and. (nint(stmp*1000)                            .EQ. -1554)
  stmp = ext_thiem3d(0.3_sp, (/Ktest , vartest , corrtest , anisotest/), &
       (/Reftest , hreftest , Qwtest , Ltest  /), .false., 1.0_sp)
  isgood = isgood .and. (nint(stmp*1000)                            .EQ. -1163)
  dtmp = ext_thiem3d(0.2_dp, (/Ktestd, vartestd, corrtestd, anisotestd/), &
       (/Reftestd, hreftestd, Qwtestd, Ltestd/))
  isgood = isgood .and. (nint(dtmp*1000000).EQ. -4194527)
  dtmp = ext_thiem3d(0.3_dp, (/Ktestd, vartestd, corrtestd, anisotestd/), &
       (/Reftestd, hreftestd, Qwtestd, Ltestd/))
  isgood = isgood .and. (nint(dtmp*1000000).EQ. -3132098)
  dtmp = ext_thiem3d(0.2_dp, (/Ktestd, vartestd, corrtestd, anisotestd/), &
       (/Reftestd, hreftestd, Qwtestd, Ltestd/), .false., 1.0_dp)
  isgood = isgood .and. (nint(dtmp*1000000)                         .EQ. -1554179)
  dtmp = ext_thiem3d(0.3_dp, (/Ktestd, vartestd, corrtestd, anisotestd/), &
       (/Reftestd, hreftestd, Qwtestd, Ltestd/), .false., 1.0_dp)
  isgood = isgood .and. (nint(dtmp*1000000)                         .EQ. -1162745)
  !sp
  write (*,*) 'ext_Thiem3d_sp(0.2)_std=', ext_thiem3d(0.2_sp, (/Ktest , vartest , corrtest , anisotest /), &
       (/Reftest , hreftest , Qwtest , Ltest /))
  write (*,*) 'ext_Thiem3d_sp(0.3)_std=', ext_thiem3d(0.3_sp, (/Ktest , vartest , corrtest , anisotest /), &
       (/Reftest , hreftest , Qwtest , Ltest /))
  !with optional inputs
  write (*,*) 'ext_Thiem3d_sp(0.2)_opt=', ext_thiem3d(0.2_sp, (/Ktest , vartest , corrtest , anisotest /), &
       (/Reftest , hreftest , Qwtest , Ltest /), &
       .false., 1.0_sp)
  write (*,*) 'ext_Thiem3d_sp(0.3)_opt=', ext_thiem3d(0.3_sp, (/Ktest , vartest , corrtest , anisotest /), &
       (/Reftest , hreftest , Qwtest , Ltest /), &
       .false., 1.0_sp)
  !dp
  write (*,*) 'ext_Thiem3d_dp(0.2)_std=', ext_thiem3d(0.2_dp, (/Ktestd, vartestd, corrtestd, anisotestd/), &
       (/Reftestd, hreftestd, Qwtestd, Ltestd/))
  write (*,*) 'ext_Thiem3d_dp(0.3)_std=', ext_thiem3d(0.3_dp, (/Ktestd, vartestd, corrtestd, anisotestd/), &
       (/Reftestd, hreftestd, Qwtestd, Ltestd/))
  !with optional inputs
  write (*,*) 'ext_Thiem3d_dp(0.2)_opt=', ext_thiem3d(0.2_dp, (/Ktestd, vartestd, corrtestd, anisotestd/), &
       (/Reftestd, hreftestd, Qwtestd, Ltestd/), &
       .false., 1.0_dp)
  write (*,*) 'ext_Thiem3d_dp(0.3)_opt=', ext_thiem3d(0.3_dp, (/Ktestd, vartestd, corrtestd, anisotestd/), &
       (/Reftestd, hreftestd, Qwtestd, Ltestd/), &
       .false., 1.0_dp)

  write (*,*) ' '

  write (*,*) "EXTENDED THIEM's SOLUTION IN 2D"
  write (*,*) '-------------------------------------------------------------------------'

  ! ext_Thiem2d_sp(0.2)  =  -2.57654548
  ! ext_Thiem2d_sp(0.3)  =  -1.93043315
  ! ext_Thiem2d_sp(0.2)  =  -2.56530547           optional
  ! ext_Thiem2d_sp(0.3)  =  -1.91978836           optional
  ! ext_Thiem2d_dp(0.2)  =  -2.5765449512179259
  ! ext_Thiem2d_dp(0.3)  =  -1.9304329811242074
  ! ext_Thiem2d_dp(0.2)  =  -2.5653048154494607   optional
  ! ext_Thiem2d_dp(0.3)  =  -1.9197882246039830   optional

  stmp = ext_thiem2d(0.2_sp, (/TGtest , vartest , corrtest /), (/Reftest , hreftest , Qwtest /))
  isgood = isgood .and. (nint(stmp*1000)                    .EQ. -2577)
  stmp = ext_thiem2d(0.3_sp, (/TGtest , vartest , corrtest /), (/Reftest , hreftest , Qwtest /))
  isgood = isgood .and. (nint(stmp*1000)                    .EQ. -1930)
  stmp = ext_thiem2d(0.2_sp, (/TGtest , vartest , corrtest /), (/Reftest , hreftest , Qwtest /),1.0_sp)
  isgood = isgood .and. (nint(stmp*1000)             .EQ. -2565)
  stmp = ext_thiem2d(0.3_sp, (/TGtest , vartest , corrtest /), (/Reftest , hreftest , Qwtest /),1.0_sp)
  isgood = isgood .and. (nint(stmp*1000)             .EQ. -1920)
  dtmp = ext_thiem2d(0.2_dp, (/TGtestd, vartestd, corrtestd/), (/Reftestd, hreftestd, Qwtestd/))
  isgood = isgood .and. (nint(dtmp*1000000)                 .EQ. -2576545)
  dtmp = ext_thiem2d(0.3_dp, (/TGtestd, vartestd, corrtestd/), (/Reftestd, hreftestd, Qwtestd/))
  isgood = isgood .and. (nint(dtmp*1000000)                 .EQ. -1930433)
  dtmp = ext_thiem2d(0.2_dp, (/TGtestd,vartestd,corrtestd/), (/Reftestd,hreftestd,Qwtestd/),1.0_dp)
  isgood = isgood .and. (nint(dtmp*1000000)            .EQ. -2565305)
  dtmp = ext_thiem2d(0.3_dp, (/TGtestd,vartestd,corrtestd/), (/Reftestd,hreftestd,Qwtestd/),1.0_dp)
  isgood = isgood .and. (nint(dtmp*1000000)            .EQ. -1919788)
  !sp
  write (*,*) 'ext_Thiem2d_sp(0.2)_std=', ext_thiem2d(0.2_sp, (/TGtest , vartest , corrtest /), &
       (/Reftest , hreftest , Qwtest /))
  write (*,*) 'ext_Thiem2d_sp(0.3)_std=', ext_thiem2d(0.3_sp, (/TGtest , vartest , corrtest /), &
       (/Reftest , hreftest , Qwtest /))
  !with optional inputs
  write (*,*) 'ext_Thiem2d_sp(0.2)_opt=', ext_thiem2d(0.2_sp,(/TGtest ,vartest ,corrtest /), &
       (/Reftest ,hreftest ,Qwtest /),1.0_sp)
  write (*,*) 'ext_Thiem2d_sp(0.3)_opt=', ext_thiem2d(0.3_sp,(/TGtest ,vartest ,corrtest /), &
       (/Reftest ,hreftest ,Qwtest /),1.0_sp)
  !dp
  write (*,*) 'ext_Thiem2d_dp(0.2)_std=', ext_thiem2d(0.2_dp, (/TGtestd, vartestd, corrtestd/), &
       (/Reftestd, hreftestd, Qwtestd/))
  write (*,*) 'ext_Thiem2d_dp(0.3)_std=', ext_thiem2d(0.3_dp, (/TGtestd, vartestd, corrtestd/), &
       (/Reftestd, hreftestd, Qwtestd/))
  !with optional inputs
  write (*,*) 'ext_Thiem2d_dp(0.2)_opt=', ext_thiem2d(0.2_dp,(/TGtestd,vartestd,corrtestd/), &
       (/Reftestd,hreftestd,Qwtestd/),1.0_dp)
  write (*,*) 'ext_Thiem2d_dp(0.3)_opt=', ext_thiem2d(0.3_dp,(/TGtestd,vartestd,corrtestd/), &
       (/Reftestd,hreftestd,Qwtestd/),1.0_dp)

  write (*,*) ' '

  write (*,*) "THEIS' SOLUTION"
  write (*,*) '-------------------------------------------------------------------------'

  ! Theissp(r=.2,t=.0)   =  0.00000000
  ! Theissp(r=.2,t=.1)   = -0.174579859
  ! Theissp(r=.3,t=.2)   = -0.142127588
  ! Theissp(r=.3,t=.3)   = -0.270834625
  ! Theisdp(r=.2,t=.0)   =  0.0000000000000000
  ! Theisdp(r=.2,t=.1)   = -0.17458018796997515
  ! Theisdp(r=.3,t=.2)   = -0.14212753460648780
  ! Theisdp(r=.3,t=.3)   = -0.27083461355368116

  stmp = theis(0.2_sp, 0.0_sp, (/Stest , Ktest /),   (/Qwtest , Ltest /) )
  isgood = isgood .and. (nint(stmp*1000)         .EQ. 0)
  stmp = theis(0.2_sp, 0.1_sp, (/Stest , Ktest /),   (/Qwtest , Ltest /) )
  isgood = isgood .and. (nint(stmp*1000)         .EQ. -175)
  stmp = theis(0.3_sp, 0.2_sp, (/Stest , Ktest /),   (/Qwtest , Ltest /) )
  isgood = isgood .and. (nint(stmp*1000)         .EQ. -142)
  stmp = theis(0.3_sp, 0.3_sp, (/Stest , Ktest /),   (/Qwtest , Ltest /) )
  isgood = isgood .and. (nint(stmp*1000)         .EQ. -271)
  dtmp = theis(0.2_dp, 0.0_dp, (/Stestd, Ktestd/),   (/Qwtestd, Ltestd/) )
  isgood = isgood .and. (nint(dtmp*1000000)      .EQ. 0)
  dtmp = theis(0.2_dp, 0.1_dp, (/Stestd, Ktestd/),   (/Qwtestd, Ltestd/) )
  isgood = isgood .and. (nint(dtmp*1000000)      .EQ. -174580)
  dtmp = theis(0.3_dp, 0.2_dp, (/Stestd, Ktestd/),   (/Qwtestd, Ltestd/) )
  isgood = isgood .and. (nint(dtmp*1000000)      .EQ. -142128)
  dtmp = theis(0.3_dp, 0.3_dp, (/Stestd, Ktestd/),   (/Qwtestd, Ltestd/) )
  isgood = isgood .and. (nint(dtmp*1000000)      .EQ. -270835)
  !sp
  write (*,*) 'Theissp(r=.2,t=.0)     =', theis(0.2_sp, 0.0_sp, (/Stest , Ktest /),   (/Qwtest , Ltest /) )
  write (*,*) 'Theissp(r=.2,t=.1)     =', theis(0.2_sp, 0.1_sp, (/Stest , Ktest /),   (/Qwtest , Ltest /) )
  write (*,*) 'Theissp(r=.3,t=.2)     =', theis(0.3_sp, 0.2_sp, (/Stest , Ktest /),   (/Qwtest , Ltest /) )
  write (*,*) 'Theissp(r=.3,t=.3)     =', theis(0.3_sp, 0.3_sp, (/Stest , Ktest /),   (/Qwtest , Ltest /) )
  !dp
  write (*,*) 'Theisdp(r=.2,t=.0)     =', theis(0.2_dp, 0.0_dp, (/Stestd, Ktestd/),   (/Qwtestd, Ltestd/) )
  write (*,*) 'Theisdp(r=.2,t=.1)     =', theis(0.2_dp, 0.1_dp, (/Stestd, Ktestd/),   (/Qwtestd, Ltestd/) )
  write (*,*) 'Theisdp(r=.3,t=.2)     =', theis(0.3_dp, 0.2_dp, (/Stestd, Ktestd/),   (/Qwtestd, Ltestd/) )
  write (*,*) 'Theisdp(r=.3,t=.3)     =', theis(0.3_dp, 0.3_dp, (/Stestd, Ktestd/),   (/Qwtestd, Ltestd/) )

  write (*,*) ' '

  write (*,*) "EXTENDED THEIS' SOLUTION IN 2D"
  write (*,*) '-------------------------------------------------------------------------'

  ! ext_theis2d_dp_vec(r=           1 ,t=          50 )  =   -2.6724678924208556
  ! ext_theis2d_dp_vec(r=           1 ,t=         100 )  =   -3.4409275281733089
  ! ext_theis2d_dp_vec(r=           1 ,t=        1000 )  =   -5.7384083577267218
  ! ext_theis2d_dp_vec(r=           2 ,t=          50 )  =   -1.2250879935634749
  ! ext_theis2d_dp_vec(r=           2 ,t=         100 )  =   -1.8423690881291142
  ! ext_theis2d_dp_vec(r=           2 ,t=        1000 )  =   -4.0012094551152613
  ! ext_theis2d_dp_vec(r=           3 ,t=          50 )  =   -0.60506500331545221
  ! ext_theis2d_dp_vec(r=           3 ,t=         100 )  =   -1.0950539046041157
  ! ext_theis2d_dp_vec(r=           3 ,t=        1000 )  =   -3.0545777884298797

  fpoints_dp    =             ext_theis2d(&
       rad     =   (/1._dp,2._dp,3._dp/),&
       time    =   (/50.0_dp,100.0_dp,1000.0_dp/),&
       params  =   (/TGtestd, vartestd, corrtestd, Stestd/),&
       inits   =   (/Qwtestd/))

  do m=1_i4,3_i4
     do n=1_i4,3_i4
        fpoints2_dp(m,n) =  ext_theis2d(&
             rad     =   real(m,dp),&
             time    =   real(425_i4*n**2_i4-1225_i4*n+850_i4,dp),&  !polynomial with p(1)=50 p(2)=100 p(3)=1000
             params  =   (/TGtestd, vartestd, corrtestd, Stestd/),&
             inits   =   (/Qwtestd/))
     end do
  end do

  fpoints3_dp =               reshape(&
       ext_theis2d(&
       grid    =   grid_dp,&
       params  =   (/TGtestd, vartestd, corrtestd, Stestd/),&
       inits   =   (/Qwtestd/)),(/3,3/))


  do m=1_i4,3_i4
     do n=1_i4,3_i4
        isgood = isgood .and. eq(fpoints_dp(m,n),fpoints2_dp(m,n))
        isgood = isgood .and. eq(fpoints_dp(m,n),fpoints3_dp(m,n))
     end do
  end do

  dtmp = fpoints_dp(1,1)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -2672468)
  dtmp = fpoints_dp(1,2)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -3440928)
  dtmp = fpoints_dp(1,3)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -5738408)
  dtmp = fpoints_dp(2,1)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -1225088)
  dtmp = fpoints_dp(2,2)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -1842369)
  dtmp = fpoints_dp(2,3)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -4001209)
  dtmp = fpoints_dp(3,1)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -605065)
  dtmp = fpoints_dp(3,2)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -1095054)
  dtmp = fpoints_dp(3,3)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -3054578)

  do m=1_i4,3_i4
     do n=1_i4,3_i4
        write(*,*) "ext_theis2d_dp_vec(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints_dp(m,n)
        write(*,*) "ext_theis2d_dp_sgl(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints2_dp(m,n)
        write(*,*) "ext_theis2d_dp_gri(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints3_dp(m,n)
     end do
  end do

  !### sp #########################################################################
  ! ext_theis2d_sp_vec(r=           1 ,t=          50 )  =   -2.67246795
  ! ext_theis2d_sp_vec(r=           1 ,t=         100 )  =   -3.44092774
  ! ext_theis2d_sp_vec(r=           1 ,t=        1000 )  =   -5.73840857
  ! ext_theis2d_sp_vec(r=           2 ,t=          50 )  =   -1.22508800
  ! ext_theis2d_sp_vec(r=           2 ,t=         100 )  =   -1.84236920
  ! ext_theis2d_sp_vec(r=           2 ,t=        1000 )  =   -4.00120974
  ! ext_theis2d_sp_vec(r=           3 ,t=          50 )  =   -0.605064988
  ! ext_theis2d_sp_vec(r=           3 ,t=         100 )  =   -1.09505391
  ! ext_theis2d_sp_vec(r=           3 ,t=        1000 )  =   -3.05457783

  fpoints_sp             =    ext_theis2d(&
       rad     =   (/1._sp,2._sp,3._sp/),&
       time    =   (/50.0_sp,100.0_sp,1000.0_sp/),&
       params  =   (/TGtest, vartest, corrtest, Stest/),&
       inits   =   (/Qwtest/))

  do m=1_i4,3_i4
     do n=1_i4,3_i4
        fpoints2_sp(m,n) =  ext_theis2d(&
             rad     =   real(m,sp),&
             time    =   real(425_i4*n**2_i4-1225_i4*n+850_i4,sp),&  !polynomial with p(1)=50 p(2)=100 p(3)=1000
             params  =   (/TGtest, vartest, corrtest, Stest/),&
             inits   =   (/Qwtest/))
     end do
  end do

  fpoints3_sp             =   reshape(&
       ext_theis2d(&
       grid    =   grid_sp,&
       params  =   (/TGtest, vartest, corrtest, Stest/),&
       inits   =   (/Qwtest/)),(/3,3/))

  do m=1_i4,3_i4
     do n=1_i4,3_i4
        isgood = isgood .and. eq(fpoints_sp(m,n),fpoints2_sp(m,n))
        isgood = isgood .and. eq(fpoints_sp(m,n),fpoints3_sp(m,n))
     end do
  end do

  stmp = fpoints_sp(1,1)
  isgood = isgood .and. (nint(stmp*1000) .eq. -2672)
  stmp = fpoints_sp(1,2)
  isgood = isgood .and. (nint(stmp*1000) .eq. -3441)
  stmp = fpoints_sp(1,3)
  isgood = isgood .and. (nint(stmp*1000) .eq. -5738)
  stmp = fpoints_sp(2,1)
  isgood = isgood .and. (nint(stmp*1000) .eq. -1225)
  stmp = fpoints_sp(2,2)
  isgood = isgood .and. (nint(stmp*1000) .eq. -1842)
  stmp = fpoints_sp(2,3)
  isgood = isgood .and. (nint(stmp*1000) .eq. -4001)
  stmp = fpoints_sp(3,1)
  isgood = isgood .and. (nint(stmp*1000) .eq. - 605)
  stmp = fpoints_sp(3,2)
  isgood = isgood .and. (nint(stmp*1000) .eq. -1095)
  stmp = fpoints_sp(3,3)
  isgood = isgood .and. (nint(stmp*1000) .eq. -3055)

  do m=1_i4,3_i4
     do n=1_i4,3_i4
        write(*,*) "ext_theis2d_sp_vec(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints_sp(m,n)
        write(*,*) "ext_theis2d_sp_sgl(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints2_sp(m,n)
        write(*,*) "ext_theis2d_sp_gri(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints3_sp(m,n)
     end do
  end do


  write (*,*) ' '

  write (*,*) "EXTENDED THEIS' SOLUTION IN 3D"
  write (*,*) '-------------------------------------------------------------------------'

  ! ext_theis3d_dp_vec(r=           1 ,t=          50 )  =   -2.6719739969853862
  ! ext_theis3d_dp_vec(r=           1 ,t=         100 )  =   -3.4114799964304905
  ! ext_theis3d_dp_vec(r=           1 ,t=        1000 )  =   -5.4641666480313242
  ! ext_theis3d_dp_vec(r=           2 ,t=          50 )  =   -1.2452670522580507
  ! ext_theis3d_dp_vec(r=           2 ,t=         100 )  =   -1.8293707976344475
  ! ext_theis3d_dp_vec(r=           2 ,t=        1000 )  =   -3.7458081837503667
  ! ext_theis3d_dp_vec(r=           3 ,t=          50 )  =  -0.64353916624054508
  ! ext_theis3d_dp_vec(r=           3 ,t=         100 )  =   -1.1069066253083155
  ! ext_theis3d_dp_vec(r=           3 ,t=        1000 )  =   -2.8268598328483159

  fpoints_dp    =             ext_theis3d(&
       rad     =   (/1._dp,2._dp,3._dp/),&
       time    =   (/50.0_dp,100.0_dp,1000.0_dp/),&
       params  =   (/Ktestd, vartestd, corrtestd, Stestd, anisotestd/),&
       inits   =   (/Qwtestd, Ltestd/))

  do m=1_i4,3_i4
     do n=1_i4,3_i4
        fpoints2_dp(m,n) =   ext_theis3d(&
             rad     =   real(m,dp),&
             time    =   real(425_i4*n**2_i4-1225_i4*n+850_i4,dp),&  !polynomial with p(1)=50 p(2)=100 p(3)=1000
             params  =   (/Ktestd, vartestd, corrtestd, Stestd, anisotestd/),&
             inits   =   (/Qwtestd, Ltestd/))
     end do
  end do

  fpoints3_dp             =   reshape(&
       ext_theis3d(&
       grid    =   grid_dp,&
       params  =   (/Ktestd, vartestd, corrtestd, Stestd, anisotestd/),&
       inits   =   (/Qwtestd, Ltestd/)),(/3,3/))

  do m=1_i4,3_i4
     do n=1_i4,3_i4
        isgood = isgood .and. eq(fpoints_dp(m,n), fpoints2_dp(m,n))
        isgood = isgood .and. eq(fpoints_dp(m,n), fpoints3_dp(m,n))
     end do
  end do

  dtmp = fpoints_dp(1,1)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -2671974)
  dtmp = fpoints_dp(1,2)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -3411480)
  dtmp = fpoints_dp(1,3)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -5464167)
  dtmp = fpoints_dp(2,1)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -1245267)
  dtmp = fpoints_dp(2,2)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -1829371)
  dtmp = fpoints_dp(2,3)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -3745808)
  dtmp = fpoints_dp(3,1)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -643539)
  dtmp = fpoints_dp(3,2)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -1106907)
  dtmp = fpoints_dp(3,3)
  isgood = isgood .and. (nint(dtmp*1000000) .eq. -2826860)

  do m=1_i4,3_i4
     do n=1_i4,3_i4
        write(*,*) "ext_theis3d_dp_vec(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints_dp(m,n)
        write(*,*) "ext_theis3d_dp_sgl(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints2_dp(m,n)
        write(*,*) "ext_theis3d_dp_gri(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints3_dp(m,n)
     end do
  end do

  ! ext_theis3d_sp_vec(r=           1 ,t=          50 )  =   -2.67197418
  ! ext_theis3d_sp_vec(r=           1 ,t=         100 )  =   -3.41148019
  ! ext_theis3d_sp_vec(r=           1 ,t=        1000 )  =   -5.46416712
  ! ext_theis3d_sp_vec(r=           2 ,t=          50 )  =   -1.24526703
  ! ext_theis3d_sp_vec(r=           2 ,t=         100 )  =   -1.82937086
  ! ext_theis3d_sp_vec(r=           2 ,t=        1000 )  =   -3.74580836
  ! ext_theis3d_sp_vec(r=           3 ,t=          50 )  =  -0.643539190
  ! ext_theis3d_sp_vec(r=           3 ,t=         100 )  =   -1.10690665
  ! ext_theis3d_sp_vec(r=           3 ,t=        1000 )  =   -2.82685995

  fpoints_sp              =   ext_theis3d(&
       rad     =   (/1._sp,2._sp,3._sp/),&
       time    =   (/50.0_sp,100.0_sp,1000.0_sp/),&
       params  =   (/Ktest, vartest, corrtest, Stest, anisotest/),&
       inits   =   (/Qwtest, Ltest/))

  do m=1_i4,3_i4
     do n=1_i4,3_i4
        fpoints2_sp(m,n) =   ext_theis3d(&
             rad     =   real(m,sp),&
             time    =   real(425_i4*n**2_i4-1225_i4*n+850_i4,sp),&  !polynomial with p(1)=50 p(2)=100 p(3)=1000
             params  =   (/Ktest, vartest, corrtest, Stest, anisotest/),&
             inits   =   (/Qwtest, Ltest/))
     end do
  end do

  fpoints3_sp             =   reshape(&
       ext_theis3d(&
       grid    =   grid_sp,&
       params  =   (/Ktest, vartest, corrtest, Stest, anisotest/),&
       inits   =   (/Qwtest, Ltest/)),(/3,3/))

  do m=1_i4,3_i4
     do n=1_i4,3_i4
        isgood = isgood .and. eq(fpoints_sp(m,n), fpoints2_sp(m,n))
        isgood = isgood .and. eq(fpoints_sp(m,n), fpoints3_sp(m,n))
     end do
  end do

  stmp = fpoints_sp(1,1)
  isgood = isgood .and. (nint(stmp*1000) .eq. -2672)
  stmp = fpoints_sp(1,2)
  isgood = isgood .and. (nint(stmp*1000) .eq. -3411)
  stmp = fpoints_sp(1,3)
  isgood = isgood .and. (nint(stmp*1000) .eq. -5464)
  stmp = fpoints_sp(2,1)
  isgood = isgood .and. (nint(stmp*1000) .eq. -1245)
  stmp = fpoints_sp(2,2)
  isgood = isgood .and. (nint(stmp*1000) .eq. -1829)
  stmp = fpoints_sp(2,3)
  isgood = isgood .and. (nint(stmp*1000) .eq. -3746)
  stmp = fpoints_sp(3,1)
  isgood = isgood .and. (nint(stmp*1000) .eq. -644)
  stmp = fpoints_sp(3,2)
  isgood = isgood .and. (nint(stmp*1000) .eq. -1107)
  stmp = fpoints_sp(3,3)
  isgood = isgood .and. (nint(stmp*1000) .eq. -2827)

  do m=1_i4,3_i4
     do n=1_i4,3_i4
        write(*,*) "ext_theis3d_sp_vec(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints_sp(m,n)
        write(*,*) "ext_theis3d_sp_sgl(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints2_sp(m,n)
        write(*,*) "ext_theis3d_sp_gri(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints3_sp(m,n)
     end do
  end do


  write (*,*) ' '

  write (*,*) '-------------------------------------------------------------------------'
  write (*,*) '-------------------------------------------------------------------------'
  if (isgood) then
     write(*,*) 'mo_pumpingtests ', color('o.k.', c_green)
  else
     write(*,*) 'mo_pumpingtests ', color('failed', c_red)
  endif

end program test

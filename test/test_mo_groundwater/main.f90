program test

!---------------------------------------------------------------------------------------------------
! MODULS AND VARIABLS
!---------------------------------------------------------------------------------------------------


USE mo_kind,            ONLY: sp, dp
USE mo_groundwater,     ONLY: thiem, theis, ext_thiem2d, ext_thiem3d

Implicit none

REAL(sp)            :: Reftest, hreftest, Qwtest, Ktest, Ltest, Stest, TGtest, vartest, proptest, corrtest, anisotest
REAL(dp)            :: Reftestd, hreftestd, Qwtestd, Ktestd, Ltestd, Stestd, TGtestd, vartestd, proptestd, corrtestd, anisotestd
LOGICAL             :: isgood

!---------------------------------------------------------------------------------------------------
! test-parameters
!---------------------------------------------------------------------------------------------------

Reftest         = 1.0_sp
hreftest        = 0.0_sp
Qwtest          = -0.001_sp
Ktest           = 0.0001_sp
Ltest           = 1.0_sp
Stest           = 0.001_sp
TGtest          = Ltest*Ktest
vartest         = 1.0_sp
proptest        = 0.8_sp
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
proptestd       = 0.8_dp
corrtestd       = 10.0_dp
anisotestd      = 1.0_dp

isgood          = .TRUE.

!---------------------------------------------------------------------------------------------------
! tests
!---------------------------------------------------------------------------------------------------

write (*,*) ''
write (*,*) 'Test mo_groundwater.f90'
write (*,*) ''

write (*,*) 'THIEM''s SOLUTION'
write (*,*) '-------------------------------------------------------------------------'

! Thiemsp(0.2)         =  -2.56150007    
! Thiemsp(0.3)         =  -1.91618240    
! Thiemdp(0.2)         =  -2.5614999936338809     
! Thiemdp(0.3)         =  -1.9161822315668402     
  
    isgood = isgood .and. (nint(thiem(0.2_sp, (/Ktest , Ltest /), (/Reftest , hreftest , Qwtest /))*1000)       .EQ. -2562)
    isgood = isgood .and. (nint(thiem(0.3_sp, (/Ktest , Ltest /), (/Reftest , hreftest , Qwtest /))*1000)       .EQ. -1916)
    isgood = isgood .and. (nint(thiem(0.2_dp, (/Ktestd, Ltestd/), (/Reftestd, hreftestd, Qwtestd/))*1000000)    .EQ. -2561500)
    isgood = isgood .and. (nint(thiem(0.3_dp, (/Ktestd, Ltestd/), (/Reftestd, hreftestd, Qwtestd/))*1000000)    .EQ. -1916182)
!sp
    write (*,*) 'Thiemsp(0.2)       =', thiem(0.3_sp, (/Ktest , Ltest /), (/Reftest , hreftest , Qwtest /))
    write (*,*) 'Thiemsp(0.3)       =', thiem(0.3_sp, (/Ktest , Ltest /), (/Reftest , hreftest , Qwtest /))
!dp
    write (*,*) 'Thiemdp(0.2)       =', thiem(0.2_dp, (/Ktestd, Ltestd/), (/Reftestd, hreftestd, Qwtestd/))
    write (*,*) 'Thiemdp(0.3)       =', thiem(0.3_dp, (/Ktestd, Ltestd/), (/Reftestd, hreftestd, Qwtestd/))

    write (*,*) ' '

write (*,*) 'EXTENDED THIEM''s SOLUTION IN 3D'
write (*,*) '-------------------------------------------------------------------------'

! ext_Thiem3d_sp(0.2)  =  -4.19452763    
! ext_Thiem3d_sp(0.3)  =  -3.13209820    
! ext_Thiem3d_sp(0.2)  =  -1.55417931           optional
! ext_Thiem3d_sp(0.3)  =  -1.16274548           optional
! ext_Thiem3d_dp(0.2)  =  -4.1945268260359985     
! ext_Thiem3d_dp(0.3)  =  -3.1320977426835279     
! ext_Thiem3d_dp(0.2)  =  -1.5541790860603071   optional  
! ext_Thiem3d_dp(0.3)  =  -1.1627453616249193   optional

    isgood = isgood .and. (nint(ext_thiem3d(0.2_sp, (/Ktest , vartest , corrtest , anisotest , Ltest /),&
                                                    & (/Reftest , hreftest , Qwtest /))*1000)                   .EQ. -4195)
    isgood = isgood .and. (nint(ext_thiem3d(0.3_sp, (/Ktest , vartest , corrtest , anisotest , Ltest /),&
                                                    & (/Reftest , hreftest , Qwtest /))*1000)                   .EQ. -3132)
    isgood = isgood .and. (nint(ext_thiem3d(0.2_sp, (/Ktest , vartest , corrtest , anisotest , Ltest /),&
                                                            & (/Reftest , hreftest , Qwtest /),&
                                                            & .false., 1.0_sp)*1000)                            .EQ. -1554)
    isgood = isgood .and. (nint(ext_thiem3d(0.3_sp, (/Ktest , vartest , corrtest , anisotest , Ltest /),&
                                                            & (/Reftest , hreftest , Qwtest /),&
                                                            & .false., 1.0_sp)*1000)                            .EQ. -1163)
    isgood = isgood .and. (nint(ext_thiem3d(0.2_dp, (/Ktestd, vartestd, corrtestd, anisotestd, Ltestd/),&
                                                            & (/Reftestd, hreftestd, Qwtestd/))*1000000)        .EQ. -4194527)
    isgood = isgood .and. (nint(ext_thiem3d(0.3_dp, (/Ktestd, vartestd, corrtestd, anisotestd, Ltestd/),&
                                                            & (/Reftestd, hreftestd, Qwtestd/))*1000000)        .EQ. -3132098)
    isgood = isgood .and. (nint(ext_thiem3d(0.2_dp, (/Ktestd, vartestd, corrtestd, anisotestd, Ltestd/),&
                                                            & (/Reftestd, hreftestd, Qwtestd/),&
                                                            & .false., 1.0_dp)*1000000)                         .EQ. -1554179)
    isgood = isgood .and. (nint(ext_thiem3d(0.3_dp, (/Ktestd, vartestd, corrtestd, anisotestd, Ltestd/),&
                                                            & (/Reftestd, hreftestd, Qwtestd/),&
                                                            & .false., 1.0_dp)*1000000)                         .EQ. -1162745)
!sp
    write (*,*) 'ext_Thiem3d_sp(0.2)=', ext_thiem3d(0.2_sp, (/Ktest , vartest , corrtest , anisotest , Ltest /),&
                                                            & (/Reftest , hreftest , Qwtest /))
    write (*,*) 'ext_Thiem3d_sp(0.3)=', ext_thiem3d(0.3_sp, (/Ktest , vartest , corrtest , anisotest , Ltest /),&
                                                            & (/Reftest , hreftest , Qwtest /))
!with optional inputs
    write (*,*) 'ext_Thiem3d_sp(0.2)=', ext_thiem3d(0.2_sp, (/Ktest , vartest , corrtest , anisotest , Ltest /),&
                                                            & (/Reftest , hreftest , Qwtest /),&
                                                            & .false., 1.0_sp)
    write (*,*) 'ext_Thiem3d_sp(0.3)=', ext_thiem3d(0.3_sp, (/Ktest , vartest , corrtest , anisotest , Ltest /),&
                                                            & (/Reftest , hreftest , Qwtest /),&
                                                            & .false., 1.0_sp)
!dp
    write (*,*) 'ext_Thiem3d_dp(0.2)=', ext_thiem3d(0.2_dp, (/Ktestd, vartestd, corrtestd, anisotestd, Ltestd/),&
                                                            & (/Reftestd, hreftestd, Qwtestd/))
    write (*,*) 'ext_Thiem3d_dp(0.3)=', ext_thiem3d(0.3_dp, (/Ktestd, vartestd, corrtestd, anisotestd, Ltestd/),&
                                                            & (/Reftestd, hreftestd, Qwtestd/))
!with optional inputs
    write (*,*) 'ext_Thiem3d_dp(0.2)=', ext_thiem3d(0.2_dp, (/Ktestd, vartestd, corrtestd, anisotestd, Ltestd/),&
                                                            & (/Reftestd, hreftestd, Qwtestd/),&
                                                            & .false., 1.0_dp)
    write (*,*) 'ext_Thiem3d_dp(0.3)=', ext_thiem3d(0.3_dp, (/Ktestd, vartestd, corrtestd, anisotestd, Ltestd/),&
                                                            & (/Reftestd, hreftestd, Qwtestd/),&
                                                            & .false., 1.0_dp)

    write (*,*) ' '

write (*,*) 'EXTENDED THIEM''s SOLUTION IN 2D'
write (*,*) '-------------------------------------------------------------------------'

! ext_Thiem2d_sp(0.2)  =  -2.57654548    
! ext_Thiem2d_sp(0.3)  =  -1.93043315    
! ext_Thiem2d_sp(0.2)  =  -2.56530547           optional    
! ext_Thiem2d_sp(0.3)  =  -1.91978836           optional
! ext_Thiem2d_dp(0.2)  =  -2.5765449512179259     
! ext_Thiem2d_dp(0.3)  =  -1.9304329811242074     
! ext_Thiem2d_dp(0.2)  =  -2.5653048154494607   optional  
! ext_Thiem2d_dp(0.3)  =  -1.9197882246039830   optional  

    isgood = isgood .and. (nint(ext_thiem2d(0.2_sp, (/TGtest , vartest , corrtest /),&
                                                    &(/Reftest , hreftest , Qwtest /))*1000)                    .EQ. -2577)
    isgood = isgood .and. (nint(ext_thiem2d(0.3_sp, (/TGtest , vartest , corrtest /),&
                                                    &(/Reftest , hreftest , Qwtest /))*1000)                    .EQ. -1930)
    isgood = isgood .and. (nint(ext_thiem2d(0.2_sp, (/TGtest , vartest , corrtest /),&
                                                    &(/Reftest , hreftest , Qwtest /),1.0_sp)*1000)             .EQ. -2565)
    isgood = isgood .and. (nint(ext_thiem2d(0.3_sp, (/TGtest , vartest , corrtest /),&
                                                    &(/Reftest , hreftest , Qwtest /),1.0_sp)*1000)             .EQ. -1920)
    isgood = isgood .and. (nint(ext_thiem2d(0.2_dp, (/TGtestd, vartestd, corrtestd/),&
                                                    &(/Reftestd, hreftestd, Qwtestd/))*1000000)                 .EQ. -2576545)
    isgood = isgood .and. (nint(ext_thiem2d(0.3_dp, (/TGtestd, vartestd, corrtestd/),&
                                                    &(/Reftestd, hreftestd, Qwtestd/))*1000000)                 .EQ. -1930433)
    isgood = isgood .and. (nint(ext_thiem2d(0.2_dp, (/TGtestd,vartestd,corrtestd/),&
                                                    &(/Reftestd,hreftestd,Qwtestd/),1.0_dp)*1000000)            .EQ. -2565305)
    isgood = isgood .and. (nint(ext_thiem2d(0.3_dp, (/TGtestd,vartestd,corrtestd/),&
                                                    &(/Reftestd,hreftestd,Qwtestd/),1.0_dp)*1000000)            .EQ. -1919788)
!sp
    write (*,*) 'ext_Thiem2d_sp(0.2)=', ext_thiem2d(0.2_sp, (/TGtest , vartest , corrtest /), (/Reftest , hreftest , Qwtest /))
    write (*,*) 'ext_Thiem2d_sp(0.3)=', ext_thiem2d(0.3_sp, (/TGtest , vartest , corrtest /), (/Reftest , hreftest , Qwtest /))
!with optional inputs
    write (*,*) 'ext_Thiem2d_sp(0.2)=', ext_thiem2d(0.2_sp,(/TGtest , vartest , corrtest /),(/Reftest , hreftest , Qwtest /),1.0_sp)
    write (*,*) 'ext_Thiem2d_sp(0.3)=', ext_thiem2d(0.3_sp,(/TGtest , vartest , corrtest /),(/Reftest , hreftest , Qwtest /),1.0_sp)
!dp
    write (*,*) 'ext_Thiem2d_dp(0.2)=', ext_thiem2d(0.2_dp, (/TGtestd, vartestd, corrtestd/), (/Reftestd, hreftestd, Qwtestd/))
    write (*,*) 'ext_Thiem2d_dp(0.3)=', ext_thiem2d(0.3_dp, (/TGtestd, vartestd, corrtestd/), (/Reftestd, hreftestd, Qwtestd/))
!with optional inputs
    write (*,*) 'ext_Thiem2d_dp(0.2)=', ext_thiem2d(0.2_dp,(/TGtestd, vartestd, corrtestd/),(/Reftestd, hreftestd, Qwtestd/),1.0_dp)
    write (*,*) 'ext_Thiem2d_dp(0.3)=', ext_thiem2d(0.3_dp,(/TGtestd, vartestd, corrtestd/),(/Reftestd, hreftestd, Qwtestd/),1.0_dp)

    write (*,*) ' '

write (*,*) 'THEIS'' SOLUTION'
write (*,*) '-------------------------------------------------------------------------'

! Theissp(r=.2,t=.0)   =  0.00000000    
! Theissp(r=.2,t=.1)   = -0.174579859    
! Theissp(r=.3,t=.2)   = -0.142127588    
! Theissp(r=.3,t=.3)   = -0.270834625    
! Theisdp(r=.2,t=.0)   =  0.0000000000000000     
! Theisdp(r=.2,t=.1)   = -0.17458018796997515     
! Theisdp(r=.3,t=.2)   = -0.14212753460648780     
! Theisdp(r=.3,t=.3)   = -0.27083461355368116 

    isgood = isgood .and. (nint(theis(0.2_sp, 0.0_sp, (/Stest , Ktest , Ltest /),   (/Qwtest /) )*1000)         .EQ. 0)
    isgood = isgood .and. (nint(theis(0.2_sp, 0.1_sp, (/Stest , Ktest , Ltest /),   (/Qwtest /) )*1000)         .EQ. -175)
    isgood = isgood .and. (nint(theis(0.3_sp, 0.2_sp, (/Stest , Ktest , Ltest /),   (/Qwtest /) )*1000)         .EQ. -142)
    isgood = isgood .and. (nint(theis(0.3_sp, 0.3_sp, (/Stest , Ktest , Ltest /),   (/Qwtest /) )*1000)         .EQ. -271)
    isgood = isgood .and. (nint(theis(0.2_dp, 0.0_dp, (/Stestd, Ktestd, Ltestd/),   (/Qwtestd/) )*1000000)      .EQ. 0)
    isgood = isgood .and. (nint(theis(0.2_dp, 0.1_dp, (/Stestd, Ktestd, Ltestd/),   (/Qwtestd/) )*1000000)      .EQ. -174580)
    isgood = isgood .and. (nint(theis(0.3_dp, 0.2_dp, (/Stestd, Ktestd, Ltestd/),   (/Qwtestd/) )*1000000)      .EQ. -142128)
    isgood = isgood .and. (nint(theis(0.3_dp, 0.3_dp, (/Stestd, Ktestd, Ltestd/),   (/Qwtestd/) )*1000000)      .EQ. -270835)
!sp
    write (*,*) 'Theissp(r=.2,t=.0) =', theis(0.2_sp, 0.0_sp, (/Stest , Ktest , Ltest /),   (/Qwtest /) )
    write (*,*) 'Theissp(r=.2,t=.1) =', theis(0.2_sp, 0.1_sp, (/Stest , Ktest , Ltest /),   (/Qwtest /) )
    write (*,*) 'Theissp(r=.3,t=.2) =', theis(0.3_sp, 0.2_sp, (/Stest , Ktest , Ltest /),   (/Qwtest /) )
    write (*,*) 'Theissp(r=.3,t=.3) =', theis(0.3_sp, 0.3_sp, (/Stest , Ktest , Ltest /),   (/Qwtest /) )
!dp
    write (*,*) 'Theisdp(r=.2,t=.0) =', theis(0.2_dp, 0.0_dp, (/Stestd, Ktestd, Ltestd/),   (/Qwtestd/) )
    write (*,*) 'Theisdp(r=.2,t=.1) =', theis(0.2_dp, 0.1_dp, (/Stestd, Ktestd, Ltestd/),   (/Qwtestd/) )
    write (*,*) 'Theisdp(r=.3,t=.2) =', theis(0.3_dp, 0.2_dp, (/Stestd, Ktestd, Ltestd/),   (/Qwtestd/) )
    write (*,*) 'Theisdp(r=.3,t=.3) =', theis(0.3_dp, 0.3_dp, (/Stestd, Ktestd, Ltestd/),   (/Qwtestd/) )
    
    write (*,*) ' '

write (*,*) '-------------------------------------------------------------------------'
write (*,*) '-------------------------------------------------------------------------'
if (isgood) then
    write(*,*) 'mo_groundwater o.k.'
else
    write(*,*) 'mo_groundwater failed'
endif

end program test

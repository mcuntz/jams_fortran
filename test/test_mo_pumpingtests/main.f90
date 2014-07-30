program test

!---------------------------------------------------------------------------------------------------
! MODULS AND VARIABLS
!---------------------------------------------------------------------------------------------------


USE mo_kind,            ONLY: sp, dp, i4
USE mo_pumpingtests,    ONLY: thiem, theis, ext_thiem2d, ext_thiem3d, ext_theis2d, ext_theis3d
use mo_utils,           ONLY: eq

Implicit none

REAL(sp)            :: Reftest, hreftest, Qwtest, Ktest, Ltest, Stest, TGtest, vartest, corrtest, anisotest!, proptest
REAL(dp)            :: Reftestd, hreftestd, Qwtestd, Ktestd, Ltestd, Stestd, TGtestd, vartestd, corrtestd, anisotestd!, proptestd
LOGICAL             :: isgood

real(dp), dimension(3,3)   :: fpoints_dp
real(sp), dimension(3,3)   :: fpoints_sp

integer(i4)         :: m,n

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

isgood          = .TRUE.

!---------------------------------------------------------------------------------------------------
! tests
!---------------------------------------------------------------------------------------------------

write (*,*) ''
write (*,*) 'Test mo_pumpingtests.f90'
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
    write (*,*) 'Thiemsp(0.2)           =', thiem(0.2_sp, (/Ktest , Ltest /), (/Reftest , hreftest , Qwtest /))
    write (*,*) 'Thiemsp(0.3)           =', thiem(0.3_sp, (/Ktest , Ltest /), (/Reftest , hreftest , Qwtest /))
!dp
    write (*,*) 'Thiemdp(0.2)           =', thiem(0.2_dp, (/Ktestd, Ltestd/), (/Reftestd, hreftestd, Qwtestd/))
    write (*,*) 'Thiemdp(0.3)           =', thiem(0.3_dp, (/Ktestd, Ltestd/), (/Reftestd, hreftestd, Qwtestd/))

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
    write (*,*) 'ext_Thiem3d_sp(0.2)_std=', ext_thiem3d(0.2_sp, (/Ktest , vartest , corrtest , anisotest , Ltest /),&
                                                            & (/Reftest , hreftest , Qwtest /))
    write (*,*) 'ext_Thiem3d_sp(0.3)_std=', ext_thiem3d(0.3_sp, (/Ktest , vartest , corrtest , anisotest , Ltest /),&
                                                            & (/Reftest , hreftest , Qwtest /))
!with optional inputs
    write (*,*) 'ext_Thiem3d_sp(0.2)_opt=', ext_thiem3d(0.2_sp, (/Ktest , vartest , corrtest , anisotest , Ltest /),&
                                                            & (/Reftest , hreftest , Qwtest /),&
                                                            & .false., 1.0_sp)
    write (*,*) 'ext_Thiem3d_sp(0.3)_opt=', ext_thiem3d(0.3_sp, (/Ktest , vartest , corrtest , anisotest , Ltest /),&
                                                            & (/Reftest , hreftest , Qwtest /),&
                                                            & .false., 1.0_sp)
!dp
    write (*,*) 'ext_Thiem3d_dp(0.2)_std=', ext_thiem3d(0.2_dp, (/Ktestd, vartestd, corrtestd, anisotestd, Ltestd/),&
                                                            & (/Reftestd, hreftestd, Qwtestd/))
    write (*,*) 'ext_Thiem3d_dp(0.3)_std=', ext_thiem3d(0.3_dp, (/Ktestd, vartestd, corrtestd, anisotestd, Ltestd/),&
                                                            & (/Reftestd, hreftestd, Qwtestd/))
!with optional inputs
    write (*,*) 'ext_Thiem3d_dp(0.2)_opt=', ext_thiem3d(0.2_dp, (/Ktestd, vartestd, corrtestd, anisotestd, Ltestd/),&
                                                            & (/Reftestd, hreftestd, Qwtestd/),&
                                                            & .false., 1.0_dp)
    write (*,*) 'ext_Thiem3d_dp(0.3)_opt=', ext_thiem3d(0.3_dp, (/Ktestd, vartestd, corrtestd, anisotestd, Ltestd/),&
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
    write (*,*) 'ext_Thiem2d_sp(0.2)_std=', ext_thiem2d(0.2_sp, (/TGtest , vartest , corrtest /), (/Reftest , hreftest , Qwtest /))
    write (*,*) 'ext_Thiem2d_sp(0.3)_std=', ext_thiem2d(0.3_sp, (/TGtest , vartest , corrtest /), (/Reftest , hreftest , Qwtest /))
!with optional inputs
    write (*,*) 'ext_Thiem2d_sp(0.2)_opt=', ext_thiem2d(0.2_sp,(/TGtest ,vartest ,corrtest /),(/Reftest ,hreftest ,Qwtest /),1.0_sp)
    write (*,*) 'ext_Thiem2d_sp(0.3)_opt=', ext_thiem2d(0.3_sp,(/TGtest ,vartest ,corrtest /),(/Reftest ,hreftest ,Qwtest /),1.0_sp)
!dp
    write (*,*) 'ext_Thiem2d_dp(0.2)_std=', ext_thiem2d(0.2_dp, (/TGtestd, vartestd, corrtestd/), (/Reftestd, hreftestd, Qwtestd/))
    write (*,*) 'ext_Thiem2d_dp(0.3)_std=', ext_thiem2d(0.3_dp, (/TGtestd, vartestd, corrtestd/), (/Reftestd, hreftestd, Qwtestd/))
!with optional inputs
    write (*,*) 'ext_Thiem2d_dp(0.2)_opt=', ext_thiem2d(0.2_dp,(/TGtestd,vartestd,corrtestd/),(/Reftestd,hreftestd,Qwtestd/),1.0_dp)
    write (*,*) 'ext_Thiem2d_dp(0.3)_opt=', ext_thiem2d(0.3_dp,(/TGtestd,vartestd,corrtestd/),(/Reftestd,hreftestd,Qwtestd/),1.0_dp)

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
    write (*,*) 'Theissp(r=.2,t=.0)     =', theis(0.2_sp, 0.0_sp, (/Stest , Ktest , Ltest /),   (/Qwtest /) )
    write (*,*) 'Theissp(r=.2,t=.1)     =', theis(0.2_sp, 0.1_sp, (/Stest , Ktest , Ltest /),   (/Qwtest /) )
    write (*,*) 'Theissp(r=.3,t=.2)     =', theis(0.3_sp, 0.2_sp, (/Stest , Ktest , Ltest /),   (/Qwtest /) )
    write (*,*) 'Theissp(r=.3,t=.3)     =', theis(0.3_sp, 0.3_sp, (/Stest , Ktest , Ltest /),   (/Qwtest /) )
!dp
    write (*,*) 'Theisdp(r=.2,t=.0)     =', theis(0.2_dp, 0.0_dp, (/Stestd, Ktestd, Ltestd/),   (/Qwtestd/) )
    write (*,*) 'Theisdp(r=.2,t=.1)     =', theis(0.2_dp, 0.1_dp, (/Stestd, Ktestd, Ltestd/),   (/Qwtestd/) )
    write (*,*) 'Theisdp(r=.3,t=.2)     =', theis(0.3_dp, 0.2_dp, (/Stestd, Ktestd, Ltestd/),   (/Qwtestd/) )
    write (*,*) 'Theisdp(r=.3,t=.3)     =', theis(0.3_dp, 0.3_dp, (/Stestd, Ktestd, Ltestd/),   (/Qwtestd/) )
    
    write (*,*) ' '

write (*,*) 'EXTENDED THEIS'' SOLUTION IN 2D'
write (*,*) '-------------------------------------------------------------------------'

! ext_theis2d_dp(r=1 ,t=  50 )  =   -2.6490725920191469     
! ext_theis2d_dp(r=1 ,t= 100 )  =   -3.5064712797000333     
! ext_theis2d_dp(r=1 ,t=1000 )  =   -6.4798765369790310     
! ext_theis2d_dp(r=2 ,t=  50 )  =   -1.1832411219975087     
! ext_theis2d_dp(r=2 ,t= 100 )  =   -1.9058011979543912     
! ext_theis2d_dp(r=2 ,t=1000 )  =   -4.7465171853361445     
! ext_theis2d_dp(r=3 ,t=  50 )  =  -0.60054292194236258     
! ext_theis2d_dp(r=3 ,t= 100 )  =   -1.1578489936001926     
! ext_theis2d_dp(r=3 ,t=1000 )  =   -3.8058729235417590     

fpoints_dp    =   ext_theis2d(&
                    rad     =   (/1._dp,2._dp,3._dp/),&
                    time    =   (/50.0_dp,100.0_dp,1000.0_dp/),&
                    params  =   (/0.0001_dp, 1.0_dp, 10.0_dp, 0.001_dp/),&
                    inits   =   (/0.01_dp,-0.001_dp/)&
                    )

do m=1_i4,3_i4
    do n=1_i4,3_i4
        isgood = isgood .and. &
                    eq(ext_theis2d(&
                        rad=real(m,dp),&
                        time=real(425_i4*n**2_i4-1225_i4*n+850_i4,dp),&  !polynomial with p(1)=50 p(2)=100 p(3)=1000
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/)),&
                    fpoints_dp(m,n))
    end do
end do        

    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=1.0_dp,&
                        time=50.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -2649072)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=1.0_dp,&
                        time=100.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -3506471)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=1.0_dp,&
                        time=1000.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -6479876)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=2.0_dp,&
                        time=50.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -1183241)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=2.0_dp,&
                        time=100.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -1905801)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=2.0_dp,&
                        time=1000.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -4746517)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=3.0_dp,&
                        time=50.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -600542)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=3.0_dp,&
                        time=100.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -1157848)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=3.0_dp,&
                        time=1000.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -3805872)

do m=1_i4,3_i4
do n=1_i4,3_i4
    write(*,*) "ext_theis2d_dp_vec(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints_dp(m,n)
    write(*,*) "ext_theis2d_dp_sgl(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", ext_theis2d(&
                        rad=real(m,dp),&
                        time=real(425_i4*n**2_i4-1225_i4*n+850_i4,dp),&  !polynomial with p(1)=50 p(2)=100 p(3)=1000
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))
end do
end do

! ext_theis2d_sp(r=1 ,t=  50 )  =   -2.64904737    
! ext_theis2d_sp(r=1 ,t= 100 )  =   -3.50652599    
! ext_theis2d_sp(r=1 ,t=1000 )  =   -6.48010302    
! ext_theis2d_sp(r=2 ,t=  50 )  =   -1.18317485    
! ext_theis2d_sp(r=2 ,t= 100 )  =   -1.90604305    
! ext_theis2d_sp(r=2 ,t=1000 )  =   -4.74695349    
! ext_theis2d_sp(r=3 ,t=  50 )  =  -0.600474238    
! ext_theis2d_sp(r=3 ,t= 100 )  =   -1.15807271    
! ext_theis2d_sp(r=3 ,t=1000 )  =   -3.80635810    

fpoints_sp    =   ext_theis2d(&
                    rad     =   (/1._sp,2._sp,3._sp/),&
                    time    =   (/50.0_sp,100.0_sp,1000.0_sp/),&
                    params  =   (/0.0001_sp, 1.0_sp, 10.0_sp, 0.001_sp/),&
                    inits   =   (/0.01_sp,-0.001_sp/)&
                    )

do m=1_i4,3_i4
    do n=1_i4,3_i4
        isgood = isgood .and. &
                    eq(ext_theis2d(&
                        rad=real(m,sp),&
                        time=real(425_i4*n**2_i4-1225_i4*n+850_i4,sp),&  !polynomial with p(1)=50 p(2)=100 p(3)=1000
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/)),&
                    fpoints_sp(m,n))
    end do
end do        

    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=1.0_sp,&
                        time=50.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -2649) !-2648
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=1.0_sp,&
                        time=100.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -3506)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=1.0_sp,&
                        time=1000.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -6480)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=2.0_sp,&
                        time=50.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -1183)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=2.0_sp,&
                        time=100.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -1906)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=2.0_sp,&
                        time=1000.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -4746)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=3.0_sp,&
                        time=50.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -600)
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=3.0_sp,&
                        time=100.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -1158) !-1157
    isgood = isgood .and. &
                (int(ext_theis2d(&
                        rad=3.0_sp,&
                        time=1000.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -3806)

do m=1_i4,3_i4
do n=1_i4,3_i4
    write(*,*) "ext_theis2d_sp_vec(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints_sp(m,n)
    write(*,*) "ext_theis2d_sp_sgl(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", ext_theis2d(&
                        rad=real(m,sp),&
                        time=real(425_i4*n**2_i4-1225_i4*n+850_i4,sp),&  !polynomial with p(1)=50 p(2)=100 p(3)=1000
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))

end do
end do


    write (*,*) ' '

write (*,*) 'EXTENDED THEIS'' SOLUTION IN 3D'
write (*,*) '-------------------------------------------------------------------------'

! ext_theis3d_dp(r=1 ,t=  50 )  =   -2.6560270406458253     
! ext_theis3d_dp(r=1 ,t= 100 )  =   -3.5135670372984724     
! ext_theis3d_dp(r=1 ,t=1000 )  =   -6.4871018871446848     
! ext_theis3d_dp(r=2 ,t=  50 )  =   -1.2070941817599892     
! ext_theis3d_dp(r=2 ,t= 100 )  =   -1.9315210194279706     
! ext_theis3d_dp(r=2 ,t=1000 )  =   -4.7741294815995907     
! ext_theis3d_dp(r=3 ,t=  50 )  =  -0.64295549214806225     
! ext_theis3d_dp(r=3 ,t= 100 )  =   -1.2072969523759849     
! ext_theis3d_dp(r=3 ,t=1000 )  =   -3.8633829776289548     

fpoints_dp    =   ext_theis3d(&
                    rad     =   (/1._dp,2._dp,3._dp/),&
                    time    =   (/50.0_dp,100.0_dp,1000.0_dp/),&
                    params  =   (/0.0001_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 0.001_dp/),&
                    inits   =   (/0.01_dp,-0.001_dp/)&
                    )

do m=1_i4,3_i4
    do n=1_i4,3_i4
        isgood = isgood .and. &
                    eq(ext_theis3d(&
                        rad=real(m,dp),&
                        time=real(425_i4*n**2_i4-1225_i4*n+850_i4,dp),&  !polynomial with p(1)=50 p(2)=100 p(3)=1000
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/)),&
                    fpoints_dp(m,n))
    end do
end do        

    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=1.0_dp,&
                        time=50.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -2656027)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=1.0_dp,&
                        time=100.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -3513567)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=1.0_dp,&
                        time=1000.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -6487101)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=2.0_dp,&
                        time=50.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -1207094)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=2.0_dp,&
                        time=100.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -1931521)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=2.0_dp,&
                        time=1000.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -4774129)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=3.0_dp,&
                        time=50.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -642955)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=3.0_dp,&
                        time=100.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -1207296)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=3.0_dp,&
                        time=1000.0_dp,&
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))*1000000)&
                .eq. -3863382)

do m=1_i4,3_i4
do n=1_i4,3_i4
    write(*,*) "ext_theis3d_dp_vec(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints_dp(m,n)
    write(*,*) "ext_theis3d_dp_sgl(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", ext_theis3d(&
                        rad=real(m,dp),&
                        time=real(425_i4*n**2_i4-1225_i4*n+850_i4,dp),&  !polynomial with p(1)=50 p(2)=100 p(3)=1000
                        params=(/0.0001_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 0.001_dp/),&
                        inits=(/0.01_dp,-0.001_dp/))
end do
end do

! ext_theis3d_sp(r=1 ,t=  50 )  =   -2.65659475    
! ext_theis3d_sp(r=1 ,t= 100 )  =   -3.51366711    
! ext_theis3d_sp(r=1 ,t=1000 )  =   -6.48712111    
! ext_theis3d_sp(r=2 ,t=  50 )  =   -1.20790708    
! ext_theis3d_sp(r=2 ,t= 100 )  =   -1.93158758    
! ext_theis3d_sp(r=2 ,t=1000 )  =   -4.77424622    
! ext_theis3d_sp(r=3 ,t=  50 )  =  -0.644515157    
! ext_theis3d_sp(r=3 ,t= 100 )  =   -1.20738590    
! ext_theis3d_sp(r=3 ,t=1000 )  =   -3.86354828    

fpoints_sp    =   ext_theis3d(&
                    rad     =   (/1._sp,2._sp,3._sp/),&
                    time    =   (/50.0_sp,100.0_sp,1000.0_sp/),&
                    params  =   (/0.0001_sp, 1.0_sp, 10.0_sp, 1.0_sp, 1.0_sp, 0.001_sp/),&
                    inits   =   (/0.01_sp,-0.001_sp/)&
                    )

do m=1_i4,3_i4
    do n=1_i4,3_i4
        isgood = isgood .and. &
                    eq(ext_theis3d(&
                        rad=real(m,sp),&
                        time=real(425_i4*n**2_i4-1225_i4*n+850_i4,sp),&  !polynomial with p(1)=50 p(2)=100 p(3)=1000
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 1.0_sp, 1.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/)),&
                    fpoints_sp(m,n))
    end do
end do        

    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=1.0_sp,&
                        time=50.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 1.0_sp, 1.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -2656)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=1.0_sp,&
                        time=100.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 1.0_sp, 1.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -3513)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=1.0_sp,&
                        time=1000.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 1.0_sp, 1.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -6487)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=2.0_sp,&
                        time=50.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 1.0_sp, 1.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -1207)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=2.0_sp,&
                        time=100.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 1.0_sp, 1.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -1931)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=2.0_sp,&
                        time=1000.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 1.0_sp, 1.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -4774)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=3.0_sp,&
                        time=50.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 1.0_sp, 1.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -644)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=3.0_sp,&
                        time=100.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 1.0_sp, 1.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -1207)
    isgood = isgood .and. &
                (int(ext_theis3d(&
                        rad=3.0_sp,&
                        time=1000.0_sp,&
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 1.0_sp, 1.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))*1000)&
                .eq. -3863)

do m=1_i4,3_i4
do n=1_i4,3_i4
    write(*,*) "ext_theis3d_sp_vec(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", fpoints_sp(m,n)
    write(*,*) "ext_theis3d_sp_sgl(r=",m,",t=",425_i4*n**2_i4-1225_i4*n+850_i4,")  = ", ext_theis3d(&
                        rad=real(m,sp),&
                        time=real(425_i4*n**2_i4-1225_i4*n+850_i4,sp),&  !polynomial with p(1)=50 p(2)=100 p(3)=1000
                        params=(/0.0001_sp, 1.0_sp, 10.0_sp, 1.0_sp, 1.0_sp, 0.001_sp/),&
                        inits=(/0.01_sp,-0.001_sp/))
end do
end do


    write (*,*) ' '

write (*,*) '-------------------------------------------------------------------------'
write (*,*) '-------------------------------------------------------------------------'
if (isgood) then
    write(*,*) 'mo_pumpingtests o.k.'
else
    write(*,*) 'mo_pumpingtests failed'
endif

end program test

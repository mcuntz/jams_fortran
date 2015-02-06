program test

use mo_kind,                only: i4, sp, dp
use mo_laplace_inversion,   only: NLInvSteh
use mo_utils,               only: eq

implicit none

integer(i4)             :: n,m
real(dp), dimension(10) :: test_s, comp_func1, comp_func2
real(dp), dimension(10) :: testpara, compvalue 
real(dp), dimension(10) :: meansq_err

real(sp), dimension(10) :: test_s_sp, comp_func1_sp, comp_func2_sp
real(sp), dimension(10) :: testpara_sp, compvalue_sp 
real(sp), dimension(6)  :: meansq_err_sp 

LOGICAL                 :: isgood

!-----------------------------------------------------------------------------------------------------------------------------------
!with testpara[_sp] you can choose one of the served functions:
!       function 1 [ 1/(s+1)    ]   : (2,0,0,0,0,...)
!       function 2 [ 1/(s^2)    ]   : (0,2,0,0,0,...)
!       function 3 [24/(s^5)    ]   : (0,0,2,0,0,...)
!       function 4 [ 1/(s^2 + 1)]   : (0,0,0,2,0,...)
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!Interface for the test-functions
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE
    function testfuncvec(s, para)
        use mo_kind, only: dp
        implicit none
        real(dp), dimension(:), intent(in) :: s
        real(dp), dimension(:), intent(in) :: para
        real(dp), dimension(size(s))       :: testfuncvec
    end function testfuncvec
END INTERFACE

INTERFACE
    function testfuncsgl(s, para)
        use mo_kind, only: dp
        implicit none
        real(dp),               intent(in) :: s
        real(dp), dimension(:), intent(in) :: para
        real(dp)                           :: testfuncsgl
    end function testfuncsgl
END INTERFACE

INTERFACE
    function testfuncvecsp(s, para)
        use mo_kind, only: sp
        implicit none
        real(sp), dimension(:), intent(in) :: s
        real(sp), dimension(:), intent(in) :: para
        real(sp), dimension(size(s))       :: testfuncvecsp
    end function testfuncvecsp
END INTERFACE

INTERFACE
    function testfuncsglsp(s, para)
        use mo_kind, only: sp
        implicit none
        real(sp),               intent(in) :: s
        real(sp), dimension(:), intent(in) :: para
        real(sp)                           :: testfuncsglsp
    end function testfuncsglsp
END INTERFACE


do n=1_i4, size(test_s)
    test_s(n)   =   real(n,dp)!*0.1_dp
end do
do n=1_i4, size(test_s_sp)
    test_s_sp(n)=   real(n,sp)!*0.1_sp
end do

isgood          =   .true.

!-----------------------------------------------------------------------------------------------------------------------------------
! 1. function
!-----------------------------------------------------------------------------------------------------------------------------------

!-------
! dp
!-------

compvalue   =   (/&
  0.36786938930120000_dp,&     
  0.13539293114857939_dp,&
   4.9877996513226830E-002_dp,&
   1.8286058932031109E-002_dp,&
   6.5874045566482612E-003_dp,&
   2.2951255530484788E-003_dp,&
   7.6967353312535096E-004_dp,&
   2.6597684820688181E-004_dp,&
   1.2496193720515686E-004_dp,&
   1.0067495791304415E-004_dp/)

!write(*,*) compvalue(10)

testpara=(/2.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/)

comp_func1=exp(-test_s(:))
comp_func2=NLInvSteh(testfuncvec, testpara, test_s, 12)

write (*,*) " "
write (*,*) " "
write (*,*) "dp: Function 1/(s+1) with inverse function e^(-t) [Stehfest-limit=12]"
write (*,*) '----------------------------------------------------------------------------------------------------------------------'
write (*,*) "time t                     analytical inverse        Stehfest-inverted        relativ error"
write (*,*) '-------------------------------------------------------------------------------------------------------'
do n=1_i4,size(test_s)
    write (*,*) test_s(n), comp_func1(n), comp_func2(n), abs((comp_func2(n)-comp_func1(n))/comp_func1(n))

    isgood = isgood .and. eq(comp_func2(n),NLInvSteh(testfuncsgl, testpara, test_s(n), 12))
    isgood = isgood .and. (int((comp_func2(n)-compvalue(n))*1000000.0_dp, i4) .eq. 0_i4)
end do


meansq_err=0.0_dp
do m=1_i4, size(meansq_err)
        !write (*,*) " "
        !write (*,*) "stehfest-limit: ", 2*m
    comp_func2=NLInvSteh(testfuncvec, testpara, test_s,2*m)
        !write (*,*) " "
        !write (*,*) "Function 1/(s+1) with inverse function e^(-t)"
        !write (*,*) '-------------------------------------------------------------------------------------------------------------'
        !write (*,*) "time t                     analytical inverse        Stehfest-inverted        relativ error"
        !write (*,*) '-------------------------------------------------------------------------------------------------------'
    do n=1_i4,size(test_s)
            !write (*,*) test_s(n), comp_func1(n), comp_func2(n), abs((comp_func2(n)-comp_func1(n))/comp_func1(n))
        meansq_err(m)=meansq_err(m)+(comp_func2(n)-comp_func1(n))**2
    end do
end do

write (*,*) " "
write (*,*) "MSE between the stehfest-inverted- and analytical-function in dependency of the stehfestlimit"
write (*,*) '-------------------------------------------------------------------------------------------------------'
do m=1_i4, size(meansq_err)
    write (*,*) "stehfest-limit with MSE: ", 2*m,":", meansq_err(m)
end do
write (*,*) "least absolut MSE with Stehfest-limit: ", minloc(meansq_err)*2

!-------
! sp
!-------

compvalue_sp   =   (/&
    0.262637794_sp,&
    0.138087913_sp,&
    1.62456371E-02_sp,&
    6.90439567E-02_sp,&
    6.49825484E-03_sp,&
    -3.61014158E-03_sp,&
    -2.93968674E-02_sp,&
    2.30146535E-02_sp,&
    4.81352210E-03_sp,&
    4.22386564E-02_sp/)


testpara_sp=(/2.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp/)

comp_func1_sp=exp(-test_s_sp(:))
comp_func2_sp=NLInvSteh(testfuncvecsp, testpara_sp, test_s_sp, 12_i4)

write (*,*) " "
write (*,*) "sp: Function 1/(s+1) with inverse function e^(-t) [Stehfest-limit=12]"
write (*,*) '----------------------------------------------------------------------------------------------------------------------'
write (*,*) "time t         analytical inverse  Stehfest-inverted   relativ error"
write (*,*) '-------------------------------------------------------------------------------------------------------'
do n=1_i4,size(test_s_sp)
    write (*,*) test_s_sp(n), comp_func1_sp(n), comp_func2_sp(n), abs((comp_func2_sp(n)-comp_func1_sp(n))/comp_func1_sp(n))

    isgood = isgood .and. eq(comp_func2_sp(n),NLInvSteh(testfuncsglsp, testpara_sp, test_s_sp(n), 12))
    isgood = isgood .and. (int((comp_func2_sp(n)-compvalue_sp(n))*1000.0_sp, i4) .eq. 0_i4)
end do


meansq_err_sp=0.0_sp
do m=1_i4, size(meansq_err_sp)
!        write (*,*) " "
!        write (*,*) "stehfest-limit: ", 2*m
    comp_func2_sp=NLInvSteh(testfuncvecsp, testpara_sp, test_s_sp,2*m)
!        write (*,*) " "
!        write (*,*) "Function 1/(s+1) with inverse function e^(-t)"
!        write (*,*) '-------------------------------------------------------------------------------------------------------------'
!        write (*,*) "time t        analytical inverse  Stehfest-inverted   relativ error"
!        write (*,*) '-------------------------------------------------------------------------------------------------------'
    do n=1_i4,size(test_s_sp)
!        write (*,*) test_s_sp(n), comp_func1_sp(n), comp_func2_sp(n), abs((comp_func2_sp(n)-comp_func1_sp(n))/comp_func1_sp(n))
        meansq_err_sp(m)=meansq_err_sp(m)+(comp_func2_sp(n)-comp_func1_sp(n))**2
    end do
end do

write (*,*) " "
write (*,*) "MSE between the stehfest-inverted- and analytical-function in dependency of the stehfestlimit"
write (*,*) '-------------------------------------------------------------------------------------------------------'
do m=1_i4, size(meansq_err_sp)
    write (*,*) "stehfest-limit with MSE: ", 2*m,":", meansq_err_sp(m)
end do
write (*,*) "least absolut MSE with Stehfest-limit: ", minloc(meansq_err_sp)*2


!-----------------------------------------------------------------------------------------------------------------------------------
! 2. function
!-----------------------------------------------------------------------------------------------------------------------------------

!-------
! dp
!-------

compvalue   =   (/&
   1.0000009622103812_dp,&
   2.0000019244207623_dp,&     
   3.0000028866118105_dp,&     
   4.0000038488415246_dp,&     
   5.0000048109464998_dp,&     
   6.0000057732236209_dp,&     
   7.0000067358527165_dp,&     
   8.0000076976830492_dp,&     
   9.0000086601285041_dp,&     
   10.000009621893000_dp/)     

testpara=(/0.0_dp,2.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/)

comp_func1=test_s
comp_func2=NLInvSteh(testfuncvec, testpara, test_s, 12)

write (*,*) " "
write (*,*) " "
write(*,*) "dp: function 1/(s^2) with inverse function t [Stehfest-limit=12]"
write (*,*) '----------------------------------------------------------------------------------------------------------------------'
write (*,*) "time t                     analytical inverse        Stehfest-inverted        relativ error"
write (*,*) '-------------------------------------------------------------------------------------------------------'
do n=1_i4,size(test_s)
    write (*,*) test_s(n), comp_func1(n), comp_func2(n), abs((comp_func2(n)-comp_func1(n))/comp_func1(n))

    isgood = isgood .and. eq(comp_func2(n),NLInvSteh(testfuncsgl, testpara, test_s(n), 12))
    isgood = isgood .and. (int((comp_func2(n)-compvalue(n))*1000000.0_dp, i4) .eq. 0_i4)
end do

!-------
! sp
!-------

compvalue_sp   =   (/&
    1.00012207_sp,&
    2.00024414_sp,&
    3.00002766_sp,&
    4.00048828_sp,&
    5.02748299_sp,&
    6.00005531_sp,&
    6.81852579_sp,&
    8.00097656_sp,&
    9.01331997_sp,&
    10.0549660_sp/)

testpara_sp=(/0.0_sp,2.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp/)

comp_func1_sp=test_s_sp
comp_func2_sp=NLInvSteh(testfuncvecsp, testpara_sp, test_s_sp, 12)

write (*,*) " "
write(*,*) "sp: function 1/(s^2) with inverse function t [Stehfest-limit=12]"
write (*,*) '----------------------------------------------------------------------------------------------------------------------'
write (*,*) "time t         analytical inverse  Stehfest-inverted   relativ error"
write (*,*) '-------------------------------------------------------------------------------------------------------'
do n=1_i4,size(test_s_sp)
    write (*,*) test_s_sp(n), comp_func1_sp(n), comp_func2_sp(n), abs((comp_func2_sp(n)-comp_func1_sp(n))/comp_func1_sp(n))

    isgood = isgood .and. eq(comp_func2_sp(n),NLInvSteh(testfuncsglsp, testpara_sp, test_s_sp(n), 12))
    isgood = isgood .and. (int((comp_func2_sp(n)-compvalue_sp(n))*1000.0_sp, i4) .eq. 0_i4)
end do

!!-----------------------------------------------------------------------------------------------------------------------------------
!! 3. function makes --> problems with different compiliers
!!-----------------------------------------------------------------------------------------------------------------------------------

!!-------
!! dp
!!-------

!compvalue   =   (/&
!   1.0072917132281780_dp,&     
!   16.116667411650848_dp,&     
!   81.590628772012082_dp,&     
!   257.86667858641357_dp,&     
!   629.55732075484286_dp,&     
!   1305.4500603521933_dp,&     
!   2418.5074034629538_dp,&     
!   4125.8668573826171_dp,&     
!   6608.8409304767702_dp,&     
!   10072.917132077486_dp/)     

!testpara=(/0.0_dp,0.0_dp,2.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/)

!comp_func1=test_s(:)**4_i4
!comp_func2=NLInvSteh(testfuncvec, testpara, test_s, 12)

!write (*,*) " "
!write (*,*) " "
!write (*,*) "dp: function 24/(s^5) with inverse function t^4 [Stehfest-limit=12]"
!write (*,*) '----------------------------------------------------------------------------------------------------------------------'
!write (*,*) "time t                     analytical inverse        Stehfest-inverted        relativ error"
!write (*,*) '-------------------------------------------------------------------------------------------------------'
!do n=1_i4,size(test_s)
!    write (*,*) test_s(n), comp_func1(n), comp_func2(n), abs((comp_func2(n)-comp_func1(n))/comp_func1(n))

!    isgood = isgood .and. eq(comp_func2(n),NLInvSteh(testfuncsgl, testpara, test_s(n), 12))
!    isgood = isgood .and. (int((comp_func2(n)-compvalue(n))*1000000.0_dp, i4) .eq. 0_i4)
!end do

!!-------
!! sp
!!-------

!compvalue_sp   =   (/&
!    1.01060343_sp,&
!    16.1696548_sp,&
!    81.7020187_sp,&
!    258.714478_sp,&
!    625.435364_sp,&
!    1307.23230_sp,&
!    2387.24829_sp,&
!    4139.43164_sp,&
!    6627.02637_sp,&
!    10006.9658_sp/)

!testpara_sp=(/0.0_sp,0.0_sp,2.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp/)

!comp_func1_sp=test_s_sp(:)**4_i4
!comp_func2_sp=NLInvSteh(testfuncvecsp, testpara_sp, test_s_sp, 12)

!write (*,*) " "
!write (*,*) "sp: function 24/(s^5) with inverse function t^4 [Stehfest-limit=12]"
!write (*,*) '----------------------------------------------------------------------------------------------------------------------'
!write (*,*) "time t         analytical inverse  Stehfest-inverted   relativ error"
!write (*,*) '-------------------------------------------------------------------------------------------------------'
!do n=1_i4,size(test_s_sp)
!    write (*,*) test_s_sp(n), comp_func1_sp(n), comp_func2_sp(n), abs((comp_func2_sp(n)-comp_func1_sp(n))/comp_func1_sp(n))

!    isgood = isgood .and. eq(comp_func2_sp(n),NLInvSteh(testfuncsglsp, testpara_sp, test_s_sp(n), 12))
!    isgood = isgood .and. (int((comp_func2_sp(n)-compvalue_sp(n))*1000.0_sp, i4) .eq. 0_i4)
!end do

!-----------------------------------------------------------------------------------------------------------------------------------
! 4. function
!-----------------------------------------------------------------------------------------------------------------------------------

!-------
! dp
!-------

compvalue   =   (/&
  0.83959871008261233_dp,&    
  0.94458688164622640_dp,&     
  3.2028647167191657E-002_dp,&
 -0.82004733974486821_dp,&     
 -0.55999021955903139_dp,&     
 -2.2709467812679900E-002_dp,&
  0.22006692791712928_dp,&     
  0.21684512943141837_dp,&     
  0.13265657517960888_dp,&     
  5.5232344654943893E-002_dp/)

testpara=(/0.0_dp,0.0_dp,0.0_dp,2.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/)

comp_func1=sin(test_s(:))
comp_func2=NLInvSteh(testfuncvec, testpara, test_s, 12)

write (*,*) " "
write (*,*) " "
write (*,*) "dp: function 1/(s^2 + 1) with inverse function sin(t) [Stehfest-limit=12]"
write (*,*) '----------------------------------------------------------------------------------------------------------------------'
write (*,*) "time t                     analytical inverse        Stehfest-inverted        relativ error"
write (*,*) '-------------------------------------------------------------------------------------------------------'
do n=1_i4,size(test_s)
    write (*,*) test_s(n), comp_func1(n), comp_func2(n), abs((comp_func2(n)-comp_func1(n))/comp_func1(n))

    isgood = isgood .and. eq(comp_func2(n),NLInvSteh(testfuncsgl, testpara, test_s(n), 12))
    isgood = isgood .and. (int((comp_func2(n)-compvalue(n))*1000000.0_dp, i4) .eq. 0_i4)
end do


meansq_err=0.0_dp
do m=1_i4, size(meansq_err)
!        write (*,*) " "
!        write (*,*) "stehfest-limit: ", 2*m
    comp_func2=NLInvSteh(testfuncvec, testpara, test_s, 2*m)
!        write (*,*) " "
!        write (*,*) "function 1/(s^2 + 1) with inverse function sin(t)"
!        write (*,*) '-------------------------------------------------------------------------------------------------------------'
!        write (*,*) "time t                     analytical inverse        Stehfest-inverted        relativ error"
!        write (*,*) '-------------------------------------------------------------------------------------------------------'
    do n=1_i4,size(test_s)
!        write (*,*) test_s(n), comp_func1(n), comp_func2(n), abs((comp_func2(n)-comp_func1(n))/comp_func1(n))
        meansq_err(m)=meansq_err(m)+(comp_func2(n)-comp_func1(n))**2
    end do
end do

write (*,*) " "
write (*,*) "MSE between the stehfest-inverted- and analytical-function in dependency of the stehfestlimit"
write (*,*) '-------------------------------------------------------------------------------------------------------'
do m=1_i4,size(meansq_err)
    write (*,*) "stehfest-limit with MSE: ", 2*m,":", meansq_err(m)
end do
write (*,*) "least absolut MSE with Stehfest-limit: ", minloc(meansq_err)*2
write (*,*) ' '

!-------
! sp
!-------

compvalue_sp   =   (/&
    0.840711713_sp,&
    0.959169507_sp,&
    6.67876154E-02_sp,&
    -0.852895975_sp,&
    -0.582676828_sp,&
    6.31774813E-02_sp,&
    0.236722142_sp,&
    0.222023711_sp,&
    0.134778619_sp,&
    6.28164634E-02_sp/)

testpara_sp=(/0.0_sp,0.0_sp,0.0_sp,2.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp/)

comp_func1_sp=sin(test_s_sp(:))
comp_func2_sp=NLInvSteh(testfuncvecsp, testpara_sp, test_s_sp, 12)

write (*,*) " "
write (*,*) " "
write (*,*) "sp: function 1/(s^2 + 1) with inverse function sin(t) [Stehfest-limit=12]"
write (*,*) '----------------------------------------------------------------------------------------------------------------------'
write (*,*) "time t         analytical inverse  Stehfest-inverted   relativ error"
write (*,*) '-------------------------------------------------------------------------------------------------------'
do n=1_i4,size(test_s_sp)
    write (*,*) test_s_sp(n), comp_func1_sp(n), comp_func2_sp(n), abs((comp_func2_sp(n)-comp_func1_sp(n))/comp_func1_sp(n))

    isgood = isgood .and. eq(comp_func2_sp(n),NLInvSteh(testfuncsglsp, testpara_sp, test_s_sp(n), 12))
    isgood = isgood .and. (int((comp_func2_sp(n)-compvalue_sp(n))*1000.0_sp, i4) .eq. 0_i4)
end do


meansq_err_sp=0.0_sp
do m=1_i4, size(meansq_err_sp)
!        write (*,*) " "
!        write (*,*) "stehfest-limit: ", 2*m
    comp_func2_sp=NLInvSteh(testfuncvecsp, testpara_sp, test_s_sp, 2*m)
!        write (*,*) " "
!        write (*,*) "function 1/(s^2 + 1) with inverse function sin(t)"
!        write (*,*) '-------------------------------------------------------------------------------------------------------------'
!        write (*,*) "time t                     analytical inverse        Stehfest-inverted        relativ error"
!        write (*,*) '-------------------------------------------------------------------------------------------------------'
    do n=1_i4,size(test_s_sp)
!        write (*,*) test_s_sp(n), comp_func1_sp(n), comp_func2_sp(n), abs((comp_func2_sp(n)-comp_func1_sp(n))/comp_func1_sp(n))
        meansq_err_sp(m)=meansq_err_sp(m)+(comp_func2_sp(n)-comp_func1_sp(n))**2
    end do
end do

write (*,*) " "
write (*,*) "MSE between the stehfest-inverted- and analytical-function in dependency of the stehfestlimit"
write (*,*) '-------------------------------------------------------------------------------------------------------'
do m=1_i4, size(meansq_err_sp)
    write (*,*) "stehfest-limit with MSE: ", 2*m,":", meansq_err_sp(m)
end do
write (*,*) "least absolut MSE with Stehfest-limit: ", minloc(meansq_err_sp)*2
write (*,*) ' '

!----------------------------------------------------------------------------------------------------------------

write (*,*) '-------------------------------------------------------------------------'
write (*,*) '-------------------------------------------------------------------------'
write (*,*) ' '
if (isgood) then
    write(*,*) 'mo_laplace_inversion o.k.'
else
    write(*,*) 'mo_laplace_inversion failed'
endif

end program test


function testfuncvec(s, para)
    use mo_kind, only: i4, dp
    implicit none
    real(dp), dimension(:), intent(in) :: s
    real(dp), dimension(:), intent(in) :: para
    real(dp), dimension(size(s))       :: testfuncvec

    testfuncvec =   0.0_dp
    
    if (para(1)>1) testfuncvec=1.0_dp/(s(:) + 1.0_dp)
    if (para(2)>1) testfuncvec=1.0_dp/(s(:)**2_i4)
    if (para(3)>1) testfuncvec=24.0_dp/(s(:)**5_i4)
    if (para(4)>1) testfuncvec=1.0_dp/(s(:)**2_i4 +1.0_dp)
    
end function testfuncvec

function testfuncsgl(s, para)
    use mo_kind, only: i4, dp
    implicit none
    real(dp),               intent(in) :: s
    real(dp), dimension(:), intent(in) :: para
    real(dp)                           :: testfuncsgl
    
    testfuncsgl =   0.0_dp
    
    if (para(1)>1) testfuncsgl=1.0_dp/(s + 1.0_dp)
    if (para(2)>1) testfuncsgl=1.0_dp/(s**2_i4)
    if (para(3)>1) testfuncsgl=24.0_dp/(s**5_i4)
    if (para(4)>1) testfuncsgl=1.0_dp/(s**2_i4 +1.0_dp)
    
end function testfuncsgl

function testfuncvecsp(s, para)
    use mo_kind, only: i4, sp
    implicit none
    real(sp), dimension(:), intent(in) :: s
    real(sp), dimension(:), intent(in) :: para
    real(sp), dimension(size(s))       :: testfuncvecsp

    testfuncvecsp =   0.0_sp
    
    if (para(1)>1) testfuncvecsp=1.0_sp/(s(:) + 1.0_sp)
    if (para(2)>1) testfuncvecsp=1.0_sp/(s(:)**2_i4)
    if (para(3)>1) testfuncvecsp=24.0_sp/(s(:)**5_i4)
    if (para(4)>1) testfuncvecsp=1.0_sp/(s(:)**2_i4 +1.0_sp)
    
end function testfuncvecsp

function testfuncsglsp(s, para)
    use mo_kind, only: i4, sp
    implicit none
    real(sp),               intent(in) :: s
    real(sp), dimension(:), intent(in) :: para
    real(sp)                           :: testfuncsglsp
    
    testfuncsglsp =   0.0_sp
    
    if (para(1)>1) testfuncsglsp=1.0_sp/(s + 1.0_sp)
    if (para(2)>1) testfuncsglsp=1.0_sp/(s**2_i4)
    if (para(3)>1) testfuncsglsp=24.0_sp/(s**5_i4)
    if (para(4)>1) testfuncsglsp=1.0_sp/(s**2_i4 +1.0_sp)
    
end function testfuncsglsp

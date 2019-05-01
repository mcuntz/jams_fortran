!> \file mo_c13o2_photosynthesis.f90

!> \brief Photosynthetic isotope discrimination

!> \details Photosynthetic isotope discrimination
!> after Farquhar et al. (1982), Evans et al. (1986), and Lloyd & Farquhar (1994).
!> It includes extensions similar to Wingate et al. (2007) so that discrimination
!> can be calculated over the whole diurnal cycle.

!> \author Matthias Cuntz
!> \date Apr 2019

MODULE mo_c13o2_photosynthesis

  use mo_kind, only: dp, i4

  implicit none

  private

  ! Public routines and variables
  
  public :: c13o2_discrimination        ! Photosynthetic 13C discrimination
  public :: c13o2_discrimination_simple ! Simple 13C discrimination = a+(b-a)*ci/ca
  public :: init_ratio_leaf_pools       ! Allocate and initialise isotope ratios of leaf carbohydrate pools

  ! Transitory starch concentration in leaf [mol(CO2)/m2]
  real(dp), dimension(:), allocatable, public :: Vstarch
  ! Isotopic composition if leaf sucrose
  real(dp), dimension(:), allocatable, public :: Rsucrose
  ! Isotopic composition if pool used for photorespiration
  real(dp), dimension(:), allocatable, public :: Rphoto
  ! Isotopic composition if transitory starch
  real(dp), dimension(:), allocatable, public :: Rstarch

  
  ! Private parameters
  
  ! Mesophyll conductance factor for C3; C4=2*C3: gm(co2) = meso_cond_fac*vcmax
  real(dp), parameter :: meso_cond_fac = 4000.0_dp ! from LSM

  ! fraction of leaf respiration in mesophyll in C4: Rdm = frdm*Rd
  real(dp), parameter :: frdm = 0.5_dp    ! von Caemmerer (2000), table 4.1

  ! mesophyll to bundle sheat conductance for CO2 in C4, von Caemmerer (2000) Section 4.3.2.2
  real(dp), parameter :: gbsc = 2.e-3_dp 

  ! Leakage to PEP carboxylation ratio in C4 plants during the day, von Caemmerer (2000) Section 4.3.3.2
  real(dp), parameter :: Phi  = 0.2_dp ! von caemmerer (2000), fig. 4.11

  ! Fractionation of leaf respiration
  ! Ghashghaie et al. (2003) reviews around -6
  ! Tcherkez et al. (2004) calculates -5.5
  ! Gessler et al. (2008) determines -2 for Ricinus
  real(dp), parameter :: eps_e = -6.0e-03_dp

  ! Fractionation of photo respiration
  ! Ghashghaie et al. (2003) reviews around +10
  ! Tcherkez et al. (2004) calculates -9.2
  real(dp), parameter :: eps_f = 10.0e-03_dp

  ! Parameters for photosynthesis discrimination: Lloyd and Farquhar (1994)
  real(dp), parameter :: &
       eps_a       =  4.4e-03_dp, & ! diff. frac.
       eps_al      =  0.7e-03_dp, & ! diff. in liquid
       eps_b3      = 30.0e-03_dp, & ! frac. during RUBISCO carbox.
       eps_b_ci    = 27.0e-03_dp, & ! effec. carbox. frac. in 'simple' model
       beta        =  0.05_dp,    & ! amount of PEP carbox. in C3
       Phi_lpj     =  0.40_dp       ! Phi used in LPJ

  ! Sucrose content of leaf
  ! 100-500 mol(C6)/gDW guess from figures in Grimmer & Komor (1999) and Komor(2000) for Ricinus
  ! Conversion to mol(CO2)/m2(leaf) based on Ricinus experiment with Claudia Keitel
  ! -> mol(C6)/gDW * C/C6 * gDW/leaf * leaf/cm2(leaf) * cm2(leaf)/m2(leaf)
  real(dp), parameter :: Vsucrose = 300.e-6_dp * 6._dp * 2.2_dp * 1.0_dp/500._dp * 1.e4_dp

  ! Photorespiration pool
  ! This is some mixed sugar/enzyme pool.
  ! Because 90% of the sugars in leaves are sucrose, this can be maximum Vsucrose/10.
  ! This is definitely too large. It is little, though, compared to the oxygenation flux and hence mixes rapidly.
  ! Mathematically, we only need that the isotope ratio is not equal the assimilation.
  real(dp), parameter :: Vphoto = Vsucrose/10._dp

  ! Fractionation of starch synthesis
  ! Tcherkez et al. (2004) models about Rchloroplast/Rinput=1.006,
  !   i.e. eps_starch=-6
  ! The measured equilibrium fractionation is -4.4 (according to Gerd Gleixner)
  real(dp), parameter :: eps_starch = -4.4e-03_dp

  ! Fractionation of sucrose
  ! Tcherkez et al. (2004) models about Rcytoplasm/Rinput=1.0019,
  !   i.e. eps_sucrose=-2
  ! Also in this simple model sucrose and starch are linked (ca.):
  ! A*RA = F_starch*(1-eps_starch)*RA + F_sucrose*(1-eps_sucrose)*RA
  ! so if F_starch = A/3 and F_sucrose=2A/3 then eps_sucrose=-eps_starch/2
  ! This is the 'wrong' direction, according to Tcherkez et al. (2004). 
  ! But the effect on sucrose is rather small so that we rather 'close the mass balance'.
  real(dp), parameter :: eps_sucrose = -eps_starch/2._dp

contains

  subroutine c13o2_discrimination(dt, isc3, &
       ! -- Input
       ! Photosynthesis variables
       Vcmax, GPP, Rd, Gammastar, &
       ! CO2 concentrations
       ca, ci, &
       ! Conductances
       ga, gb, gs, &
       ! leaf temperature
       Tl, &
       ! Ambient isotope ratio
       Rair, &
       ! Starch pool and isotope ratios of pools for respiration
       Vstarch, Rsucrose, Rphoto, Rstarch, &
       ! -- Output
       ! discrimination
       Disc, &
       ! 13CO2 flux
       Ass13)

    use mo_utils,     only: ne, le
    use mo_constants, only: twothird_dp, T0_dp
    use mo_isotope,   only: ratio_diff_air2vap, ratio_boundary_air2vap

    implicit none

    real(dp),               intent(in)    :: dt                   ! time step
    logical,  dimension(:), intent(in)    :: isc3                 ! C3 mask
    real(dp), dimension(:), intent(in)    :: Vcmax                ! Vcmax(Tl) [mol(air)/m2s]
    real(dp), dimension(:), intent(in)    :: GPP                  ! A+Rd [mol(co2)/m2s]
    real(dp), dimension(:), intent(in)    :: Rd                   ! Leaf respiration Rd [mol(co2)/m2s]
    real(dp), dimension(:), intent(in)    :: Gammastar            ! CO2 compensation point Gamma* [ppm]
    real(dp), dimension(:), intent(in)    :: ci                   ! Stomatal CO2 concentration [ppm]
    real(dp), dimension(:), intent(in)    :: ca                   ! Ambient CO2 concentration [ppm]
    real(dp), dimension(:), intent(in)    :: gs                   ! Stomatal conductance [mol(air)/m2s]
    real(dp), dimension(:), intent(in)    :: ga                   ! Aerodynamic conductance [mol(air)/m2s]
    real(dp), dimension(:), intent(in)    :: gb                   ! Boundary layer conductance [mol(air)/m2s]
    real(dp), dimension(:), intent(in)    :: Tl                   ! Leaf temperature [K]
    real(dp), dimension(:), intent(in)    :: Rair                 ! Isotopic composition of ambient CO2
    real(dp), dimension(:), intent(inout) :: Vstarch              ! Transitory starch pool [mol(CO2)/m2]
    real(dp), dimension(:), intent(inout) :: Rsucrose             ! Isotopic composition of sucrose pool
    real(dp), dimension(:), intent(inout) :: Rphoto               ! Isotopic composition of pool for photorespiration
    real(dp), dimension(:), intent(inout) :: Rstarch              ! Isotopic composition of transitory starch
    real(dp), dimension(:), intent(out)   :: Disc                 ! Discrimination
    real(dp), dimension(:), intent(out)   :: Ass13                ! 13CO2 flux [mol(13CO2)/m2s]

    ! Local variables
    integer(i4) :: nn  ! number of grid points
    integer(i4) :: jl  ! counter
    real(dp)    :: tmp ! temporary real
    real(dp), dimension(size(isc3)) :: Rdm, Rds, leakage, Vp, Phi1 ! for C4
    real(dp), dimension(size(isc3)) :: Ass                         ! net assimilation
    real(dp), dimension(size(isc3)) :: k                           ! carboxylation efficiency
    real(dp), dimension(size(isc3)) :: Vc                          ! carboxylation rate
    real(dp), dimension(size(isc3)) :: Photo                       ! Photorespiration
    ! resistances and conductances
    real(dp), dimension(size(isc3)) :: ra, rb, rs                            ! resistances for water
    real(dp), dimension(size(isc3)) :: rac, rbc, rsc, rmc, rwc, rmlc, rabsmc ! resistances for CO2
    real(dp), dimension(size(isc3)) :: gmc, gabsc, gabsmc                    ! conductances for CO2
    ! CO2 concentrations
    real(dp), dimension(size(isc3)) :: cc, cbs, ctmp ! , cs, cw
    ! fractionations
    real(dp), dimension(size(isc3)) :: eps_es, eps_b4, eps_b, eps_s, eps_a_eff
    real(dp)                        :: eps_ab
    ! for pools update
    real(dp) :: add_flux ! new flux to pool
    real(dp) :: add_r    ! isotope ratio of new flux

    nn = size(isc3,1)
    
    ! For C4
    Rdm = frdm*Rd  ! C4: respiration in mesophyll
    Rds = Rd - Rdm ! C4: respiration in bundle sheat
    ! L  - leakage of bundle sheat cells
    ! Vp - Rate of PEP-Carboxylations
    ! Phi = L/Vp seems to be rather constant in C4 plants = 0.2,
    ! probably because the amount of C3 to C4 photosynthesis is
    ! regulated (von Caemmerer 2000, section 4.3.3.2).
    ! However, Vp=0 (and L<>0) at night and Phi->inf.
    ! Phi should be 1 if GPP equals the amount of respiration in
    ! the bundle sheats (GPP=Rds).
    ! We calc PhiL=1/Phi because this goes to zero at night.
    ! We take a linear relationship that is 1 at Rds and capped at 5=1/0.2
    ! PhiL(0)=0; PhiL(Rds)=1 -> PhiL = (GPP/Rds < 5.) > 0.
    ! Rds is changing with time but this effect of a changing slope is rather small.
    ! This comes to vp=gpp and leakage=rds if gpp<5*rds otherwise phiL=1/phi.
    tmp = 1.0_dp / (1.0_dp-Phi)
    where (GPP > (Rds/Phi))
       leakage = (GPP-Rds)*(Phi*tmp)
       Vp      = (GPP-Rds)*tmp
    elsewhere
       leakage = Rds
       Vp      = GPP
    end where
    ! New continuous Phi
    where (Vp > 0.0_dp)
       Phi1    = leakage/Vp
    elsewhere
       Phi1    = 0.0_dp ! dummy
    end where

    ! net assimilation
    Ass = GPP - Rd


    !
    !-- Resistances & Conductances
    !
    
    ! Aerodynamic resistance
    ra  = 1.0_dp / ga ! [mol m-2 s-1]^-1
    rac = ra          ! [mol(CO2) m-2 s-1]^-1
    ! Leaf boundary layer conductance [mol m-2 s-1]
    rb  = 1.0_dp / gb ! [mol m-2 s-1]^-1
    rbc = rb * ratio_boundary_air2vap ! [mol(CO2) m-2 s-1]^-1
    ! Stomatal resistance
    rs  = 1.0_dp / gs ! [mol m-2 s-1]^-1
    rsc = rs * ratio_diff_air2vap ! [mol(CO2) m-2 s-1]^-1
    ! Mesophyll conductance [mol(CO2) m-2 s-1]
    ! ISOLSM takes 8000 but this is for C4
    ! John Evans suggests that it should be ca. 3000 for C3
    ! gm for C4 is about twice that of C3 (Evans & v.Caemmerer 1996).
    ! John Evans is critical about temp dependence of gm from Bernacchi et al. (2002)
    ! and suggests: if any temp depence, take the one of Vcmax.
    where (isc3)                         ! c3
       gmc =         meso_cond_fac*Vcmax
    elsewhere                            ! c4
       gmc =  2.0_dp*meso_cond_fac*Vcmax
    end where
    where (gmc > 0.0_dp)
       rmc = 1.0_dp / gmc ! [mol(co2) m-2 s-1]^-1
    elsewhere
       rmc = 0.0_dp       ! dummy
    end where
    ! resistance from stoma middle to the chloroplast surface [mol(CO2) m-2 s-1]^-1
    ! 0.25: Lloyd & Farquhar (1994)
    rwc  = 0.25_dp * rmc
    ! resistance from the chloroplast surface to inside the chloroplast [mol(CO2) m-2 s-1]^-1
    rmlc = (1.0_dp-0.25_dp) * rmc

    ! Total conductance from canopy to stomata for CO2
    gabsc = 1.0_dp / (rac+rbc+rsc)
    ! Total conductance from canopy to sites of carboxylation for CO2
    rabsmc   = rac + rbc + rsc + rwc + rmlc
    gabsmc   = 1.0_dp / rabsmc

    
    !
    !-- CO2 concentrations
    !

    ! CO2 at the site of carboxylation, i.e. in the chloroplasts
    where (gmc > 0.0_dp)
       cc = ci - ass*rmc
       cc = max(cc, 1.1*Gammastar) ! from isolsm
    elsewhere
       cc = ci
    end where
    ! C4: CO2 in the bundle sheat (von Caemmerer 2000: Eq. 4.12 and Table 4.1)
    cbs = cc + leakage/gbsc
    ! ! CO2 at the leaf surface
    ! cs   = ca - ass*(rac+rbc) ! same: cs = ci + ass*rsc
    ! ! CO2 at the chloroplast surface = 1/4 down from Ci->Cc
    ! ! cf. rwc, rmlc above (0.25 from Lloyd & Farquhar 1994)
    ! cw = ci - 0.25_dp*(ci-cc)

    
    !
    !-- Other photosynthesis variables
    !

    ! Total conductance from canopy to sites of carboxylation
    gabsmc = gabsc * (ca-ci)/(ca-cc)

    ! Carboxylation efficiency: initial slope of A vs Ci
    ! From Farquhar et al. (1982), eq. B11
    ! where (ci > Gammastar)
    !    k_ci = GPP/(ci-Gammastar)
    ! elsewhere
    !    k_ci = 0.0_dp
    ! end where
    ! Same for Cc in C3 and Cbs in C4
    where (isc3)
       ctmp = cc
    elsewhere
       ctmp = cbs
    end where
    where (ctmp > Gammastar)
       k = GPP/(ctmp-Gammastar)
    elsewhere
       k = 0.0_dp
    end where
    ! Carboxylation rate
    Vc = k*ctmp
    ! Photorespiration rate = 0.5*Vo
    photo = k*Gammastar

    
    !
    !-- Fractionations
    !
    
    ! diffusion fractionation through lamina
    eps_ab = (1.0_dp+eps_a)**twothird_dp - 1.0_dp
    ! frac. during CO2 dissolution
    eps_es = (1.18_dp - 0.0041_dp*(Tl-T0_dp))*1.e-3_dp ! Vogel et al. (1970), Szaran (1998)
    ! discrimination by PEP-c (<0)
    eps_b4 = (26.19_dp - 9483._dp/Tl)*1.e-3_dp         ! Henderson et al. (1992)         
    ! effective discrimination of carboxylation in C3 plants
    eps_b = eps_b3*(1.0_dp-beta) + eps_b4*beta         ! Brugnoli & Farquhar (2000)
    ! frac. during leakage of bundle sheets in C4
    eps_s = eps_es + eps_al
    ! effective fractionation from canopy canopy to sites of carboxylation (chloroplast interior)
    eps_a_eff = (rac*0.0_dp + rbc*eps_ab + rsc*eps_a + rwc*eps_a + rmlc*eps_s) / (rabsmc)

    
    !
    !-- 13CO2 fluxes
    !
    
    do jl=1, nn
       if (isc3(jl)) then
          ! C3
          ! Same for photorespiration as Wingate et al. (2007) for leaf respiration
          Ass13(jl) = (1.0_dp-eps_a_eff(jl))*gabsmc(jl) &
               / ((1.0_dp-eps_a_eff(jl))*gabsmc(jl) + (1.0_dp-eps_b(jl))*k(jl)) &
               * ((1.0_dp-eps_b(jl))*Rair(jl)*k(jl)*ca(jl) &
               - (1.0_dp-eps_f)*Rphoto(jl)*k(jl)*Gammastar(jl) &
               - (1.0_dp-eps_e)*Rsucrose(jl)*Rd(jl))
          if (ne(Ass(jl),0.0_dp)) then
             Disc(jl) = 1.0_dp-Ass13(jl)/(Ass(jl)*Rair(jl))
          else
             Disc(jl) = 0.0_dp
          end if
          ! Night
          if (le(k(jl),0.0_dp)) then
             Ass13(jl) = (1.0_dp-eps_e)*Rstarch(jl) ! no *ass
             Disc(jl)  = 1.0_dp-Ass13(jl)/Rair(jl)
             Ass13(jl) = Ass13(jl)*Ass(jl)          ! *ass
          end if
       else
          ! C4
          ! Same for photorespiration as Wingate et al. (2007) for leaf respiration
          if (k(jl) > 0.0_dp) then
             ! Day
             tmp = cc(jl)*(1.0_dp-(1.0_dp-eps_a_eff(jl))/(1.0_dp-eps_b4(jl))*Ass(jl)/Vp(jl) &
                  - (1.0_dp-eps_s(jl))/(1.0_dp-eps_b3)*(1.0_dp-eps_a_eff(jl))/(1.0_dp-eps_b4(jl)) &
                  *Phi1(jl)*Ass(jl)/vc(jl))
             if (abs(1.0_dp-tmp/ca(jl)) > 0.01_dp) then
                Ass13(jl) = (1.0_dp-eps_a_eff(jl))*(Rair(jl)*ca(jl) & ! no *ass
                     - (1.0_dp-eps_e)/(1.0_dp-eps_b4(jl))*Rdm(jl)/Vp(jl)*Rsucrose(jl)*cc(jl) &
                     - (1.0_dp-eps_s(jl))/(1.0_dp-eps_b3)/(1.0_dp-eps_b4(jl))*Phi1(jl)*cc(jl) &
                     * ((1.0_dp-eps_e)*Rsucrose(jl)*Rd(jl) &
                     +(1.0_dp-eps_f)*Photo(jl)*Rphoto(jl))/vc(jl)) &
                     / (ca(jl)-tmp)
             else
                Ass13(jl) = 0.0_dp
             end if
             Disc(jl)  = 1.0_dp - Ass13(jl)/Rair(jl)
             Ass13(jl) = Ass13(jl)*Ass(jl)                            ! *ass
          else
             ! Night
             Ass13(jl) = (1.0_dp-eps_e)*Rstarch(jl)  ! no *ass
             Disc(jl)  = 1.0_dp - Ass13(jl)/Rair(jl)
             Ass13(jl) = Ass13(jl)*Ass(jl)           ! *ass
          end if
       end if ! end C3/C4

       !
       !-- 13CO2 pools
       !
       
       ! Update substrate pools
       ! 2/3 of assimilation is sugar, 1/3 goes into starch pool
       if (Ass(jl) > 0.0_dp) then
          ! average sucrose pool (continuous)
          add_flux     = dt*twothird_dp               ! *Ass
          add_r        = (1.0_dp-eps_sucrose) * Ass13(jl) ! /Ass
          Rsucrose(jl) = (Vsucrose*Rsucrose(jl) + add_flux*add_r) / (Vsucrose + add_flux*Ass(jl))
          ! average photorespiration pool (continuous)
          add_flux    = dt*2.0_dp*Photo(jl) ! vo=2*photo
          add_r       = (1.0_dp-eps_sucrose) * Ass13(jl)/Ass(jl)
          Rphoto(jl)  = (Vphoto*Rphoto(jl) + add_flux*add_r) / (Vphoto + add_flux)
          ! integrated starch pool (new at every day)
          add_flux    = dt*(1.0_dp-twothird_dp)*Ass(jl)
          add_r       = (1.0_dp-eps_starch) * Ass13(jl)/Ass(jl)
          Rstarch(jl) = (Vstarch(jl)*Rstarch(jl) + add_flux*add_r) / (Vstarch(jl) + add_flux)
          Vstarch(jl) = Vstarch(jl) + add_flux
       else
          ! Rsucrose(jl) = Rsucrose(jl)
          ! Rphoto(jl)   = Rphoto(jl)
          ! Rstarch(jl)  = Rstarch(jl)
          Vstarch(jl) = 0.0_dp ! start new starch integration each day
       end if
       
    end do

    return

  end subroutine c13o2_discrimination


  subroutine c13o2_discrimination_simple(isc3, &
       ! -- Input
       ! GPP and Leaf respiration
       GPP, Rd, &
       ! Ambient and stomatal CO2 concentration 
       ca, ci, &
       ! leaf temperature
       Tl, &
       ! Ambient isotope ratio
       Ra, &
       ! -- Output
       ! discrimination
       Disc, &
       ! 13CO2 flux
       Ass13)

    use mo_constants, only: T0_dp

    implicit none

    logical,  dimension(:), intent(in)  :: isc3  ! C3 mask
    real(dp), dimension(:), intent(in)  :: GPP   ! A+Rd [mol(CO2)/m2s]
    real(dp), dimension(:), intent(in)  :: Rd    ! Leaf respiration Rd [mol(CO2)/m2s]
    real(dp), dimension(:), intent(in)  :: ci    ! Stomatal CO2 concentration [ppm]
    real(dp), dimension(:), intent(in)  :: ca    ! Ambient CO2 concentration [ppm]
    real(dp), dimension(:), intent(in)  :: Tl    ! Leaf temperature [K]
    real(dp), dimension(:), intent(in)  :: Ra    ! Isotopic composition of ambient CO2
    real(dp), dimension(:), intent(out) :: Disc  ! Discrimination
    real(dp), dimension(:), intent(out) :: Ass13 ! 13CO2 flux [mol(13CO2)/m2s]

    ! Local variables
    integer(i4) :: jl, nn
    real(dp)    :: tmp
    real(dp), dimension(size(isc3)) :: Ass    ! net assimilation
    ! fractionations
    real(dp), dimension(size(isc3)) :: eps_es, eps_b4, eps_s

    nn = size(isc3,1)
    
    ! net assimilation
    Ass = GPP - Rd

    ! frac. during CO2 dissolution
    eps_es = (1.18_dp - 0.0041_dp*(Tl-T0_dp))*1.e-3_dp ! Vogel et al. (1970), Szaran (1998)
    ! discrimination by PEP-c (<0)
    eps_b4 = (26.19_dp - 9483._dp/Tl)*1.e-3_dp         ! Henderson et al. (1992)         
    ! frac. during leakage of bundle sheets in C4
    eps_s = eps_es + eps_al

    !
    !-- 13CO2 fluxes
    !

    ! ToDo: Check nighttime a or b?
    do jl=1, nn
       if (isc3(jl)) then
          ! C3
          if (GPP(jl) > 0.0_dp) then
             Disc(jl)  = eps_a+(eps_b_ci-eps_a)*ci(jl)/ca(jl)
          else
             Disc(jl)  = eps_a
          end if
          Ass13(jl) = (1.0_dp-Disc(jl))*Ass(jl)*Ra(jl)           ! D = 1 - A'/A/Ra
       else
          ! C4
          tmp = ci(jl)*(1.0_dp-(1.0_dp-eps_a)/(1.0_dp-eps_b4(jl))*(1.0_dp-Phi) &
               - (1.0_dp-eps_s(jl))/(1.0_dp-eps_b3)*(1.0_dp-eps_a)/(1.0_dp-eps_b4(jl))*Phi)
          if (ca(jl) > ci(jl)) then
             Ass13(jl) = (1.0_dp-eps_a) * ca(jl)  / (ca(jl)-tmp) ! no *A*Ra
             Disc(jl)  = 1.0_dp - Ass13(jl)                      ! D = 1 - A'/A/Ra
             Ass13(jl) = Ass13(jl)*Ass(jl)*Ra(jl)                ! *A*Ra
          else
             Disc(jl)  = eps_a
             Ass13(jl) = (1.0_dp-eps_a)*Ra(jl)*Ass(jl)
          end if
       end if ! end C3/C4

    end do
    
    return

  end subroutine c13o2_discrimination_simple

  
  ! ------------------------------------------------------------------

  !     NAME
  !         init_ratio_leaf_pools

  !     PURPOSE
  !>        \brief Allocate and initialise transitory starch pool and isotope ratios of leaf carbohydrate pools.

  !>        \details Allocate and initialises the transitory starch pool and
  !>        the isotope ratios of the leaf carbohydrate pools for sucrose,
  !>        the pool for photorespiration, and transitory starch.

  !     CALLING SEQUENCE
  !         call init_ratio_leaf_pools(Vstarch, Rsucrose, Rphoto, Rstarch, Rinitc3, Rinitc4, isc3)

  !     PARAMETER
  !>        \param[inout] "real(dp) :: Vstarch(:)"     Transitort starch pool
  !>        \param[inout] "real(dp) :: Rsucrose(:)"    Isotope ratio of leaf sucrose
  !>        \param[inout] "real(dp) :: Rphoto(:)"      Isotope ratio of substrate for photorespiration
  !>        \param[inout] "real(dp) :: Rstarch(:)"     Isotope ratio of transitory starch
  !>        \param[in]    "real(dp) :: Rinitc3"        Initial isotope ratio for C3 plants
  !>        \param[in]    "real(dp) :: Rinitc4"        Initial isotope ratio for C4 plants
  !>        \param[in]    "logical  :: isc3(:)"        C3 mask with .true. for C3 and .false. for C4 plants

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  subroutine init_ratio_leaf_pools(Vstarch, Rsucrose, Rphoto, Rstarch, Rinitc3, Rinitc4, isc3)

    implicit none

    real(dp), dimension(:), allocatable, intent(inout) :: Vstarch
    real(dp), dimension(:), allocatable, intent(inout) :: Rsucrose
    real(dp), dimension(:), allocatable, intent(inout) :: Rphoto
    real(dp), dimension(:), allocatable, intent(inout) :: Rstarch
    real(dp),                            intent(in)    :: Rinitc3
    real(dp),                            intent(in)    :: Rinitc4
    logical,  dimension(:),              intent(in)    :: isc3

    integer(i4) :: nn

    nn = size(isc3,1)

    ! allocate
    allocate(Vstarch(nn))
    allocate(Rsucrose(nn))
    allocate(Rphoto(nn))
    allocate(Rstarch(nn))
    
    ! initialise
    Vstarch  = 0.0_dp
    Rsucrose = merge(Rinitc3, Rinitc4, isc3)
    Rphoto   = merge(Rinitc3, Rinitc4, isc3)
    Rstarch  = merge(Rinitc3, Rinitc4, isc3)

    return

  end subroutine init_ratio_leaf_pools

END MODULE mo_c13o2_photosynthesis

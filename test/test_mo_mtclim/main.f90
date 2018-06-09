! ------------------------------------------------------------------------------
!
! Main Program mtclim
!
! author: Johannes Brenner
!
! created: 03.08.2016
! last update: 30.05.2017
!
! ------------------------------------------------------------------------------
!
program main
!
use mo_kind,      only: i4, dp
use mo_mtclim
!
!----------------------------------------
! variable declarations
!
integer(i4)                         :: i
real(dp), dimension(:), allocatable :: parameterset
! path namelists
character(len=*), parameter :: fileini = "../fortran/test/test_mo_mtclim/ini"
character(len=*), parameter :: fileparameters = "../fortran/test/test_mo_mtclim/parameters"
! ini namefile
! input, output file names
character(265)        :: in_meteo, outprefix
namelist / io_ini / in_meteo, outprefix
!
! controls
integer(i4) :: nhead, ndays
integer(i4) :: indewpt, outhum, inyear, lwrad, netlwrad
!
namelist / control_ini / nhead, ndays, indewpt, outhum, inyear, lwrad, netlwrad
!
! site characteristics
real(dp)              :: base_elev,base_isoh,site_lat,site_elev,site_slp,site_asp
real(dp)              :: site_isoh,site_ehoriz,site_whoriz,tmax_lr,tmin_lr
!
namelist / parameters_ini / base_elev,base_isoh,site_lat,site_elev,site_slp,site_asp,&
site_isoh,site_ehoriz,site_whoriz,tmax_lr,tmin_lr
!
! parameters
!
real(dp)      :: TBASE,ABASE,C,B0,B1,B2,RAIN_SCALAR,DIF_ALB,SC_INT,SC_SLOPE, &
SNOW_TCRIT,SNOW_TRATE,TDAYCOEF,LWAVE_COR
namelist / parameters / TBASE,ABASE,C,B0,B1,B2,RAIN_SCALAR,DIF_ALB,SC_INT,SC_SLOPE, &
     SNOW_TCRIT,SNOW_TRATE,TDAYCOEF,LWAVE_COR
!
! read namelists in ini file
!
integer(i4),  dimension(:), allocatable :: yday         !/* array of yearday values */
real(dp),     dimension(:), allocatable :: tmax         !/* array of base maximum temperature values */
real(dp),     dimension(:), allocatable :: tmin         !/* array of base minimum temperature values */
real(dp),     dimension(:), allocatable :: prcp         !/* array of base daily precipitation values */
!real(dp),     dimension(:), allocatable :: swrad_obs    !/* array of observed incoming shortwave radiation values */
!real(dp),     dimension(:), allocatable :: lwrad_obs    !/* array of observed incoming longwave radiation values */
!real(dp),     dimension(:), allocatable :: netlwrad_obs !/* array of observed net outgoing longtwave radiation values */
!real(dp),     dimension(:), allocatable :: tdew        !/* array of base dewpoint temperature values */
real(dp),     dimension(:), allocatable :: s_tmax       !/* array of site tmax values */
real(dp),     dimension(:), allocatable :: s_tmin       !/* array of site tmin values */
real(dp),     dimension(:), allocatable :: s_tday       !/* array of site daylight temperature values */
real(dp),     dimension(:), allocatable :: s_prcp       !/* array of site prcp values */
real(dp),     dimension(:), allocatable :: s_hum        !/* array of site humidity values (VPD or VP, Pa) */
real(dp),     dimension(:), allocatable :: s_srad       !/* array of site shortwave radiation values */
real(dp),     dimension(:), allocatable :: s_srad_dif   !/* array of site shortwave diffuse radiation values */
real(dp),     dimension(:), allocatable :: s_lrad       !/* array of site longwave radiation values */
real(dp),     dimension(:), allocatable :: s_lradnet    !/* array of site net longwave radiation values */
real(dp),     dimension(:), allocatable :: s_dayl       !/* array of site daylength values */
real(dp),     dimension(:), allocatable :: s_swe        !/* array of site snow water equivalent values (cm) */
!
logical     :: isgood, allgood
allgood = .true.
!
Write(*,*) ''
Write(*,*) 'Test mo_mtclim.f90'
!
! read out namefiles
!
open(700, file=fileini)
read(700, nml=io_ini)
read(700, nml=control_ini)
read(700, nml=parameters_ini)
close(700)
!
open(700, file=fileparameters)
read(700, nml=parameters)
close(700)
!
allocate(parameterset(14))
!
parameterset(1) = TBASE
parameterset(2) = ABASE
parameterset(3) = C
parameterset(4) = B0
parameterset(5) = B1
parameterset(6) = B2
parameterset(7) = RAIN_SCALAR
parameterset(8) = DIF_ALB
parameterset(9) = SC_INT
parameterset(10) = SC_SLOPE
parameterset(11) = SNOW_TCRIT
parameterset(12) = SNOW_TRATE
parameterset(13) = TDAYCOEF
parameterset(14) = LWAVE_COR
!
!----------------------------------------
! inital progress indicator
print *, "Starting MTCLIM version 4.3\n"
!
!----------------------------------------
! allocate space in the data arrays for input and output data
!
!allocate(year(ndays))
allocate(yday(ndays))
allocate(tmax(ndays))
allocate(tmin(ndays))
allocate(prcp(ndays))
!
!allocate(swrad_obs(ndays))
!allocate(lwrad_obs(ndays))
!allocate(netlwrad_obs(ndays))
!
allocate(s_swe(ndays))
!
allocate(s_tmin(ndays))
allocate(s_tmax(ndays))
allocate(s_tday(ndays))
!
allocate(s_hum(ndays))
allocate(s_srad(ndays))
allocate(s_srad_dif(ndays))
allocate(s_lrad(ndays))
allocate(s_lradnet(ndays))
allocate(s_dayl(ndays))
!
!---------------------------------------
! read meteorological data from input file into data arrays
!---------------------------------------
! read meteorological data from input file into data arrays
open(unit=1234, file=in_meteo, status='old', action='read')
do i=1, ndays
   read(1234, '(i3,f9.2,f9.2,f9.2)') yday(i), tmax(i), tmin(i), prcp(i)
enddo
close(1234)
!
!----------------------------------------
! estimate daily air temperatures
call calc_tair(tmin, tmax, ndays, site_elev, base_elev, TDAYCOEF, tmin_lr, tmax_lr, s_tmax, s_tmin, s_tday)
print *, "Completed calc_tair()"
!
!----------------------------------------
! estimate daily precipitation
s_prcp = calc_prcp(prcp, ndays, site_isoh, base_isoh)
print *, "Completed calc_prcp()"
!----------------------------------------
! estimate daily snowpack
s_swe = snowpack(s_tmin, s_prcp, yday, ndays, SNOW_TCRIT, SNOW_TRATE)
print *, "Completed snowpack()"
!
!----------------------------------------
! estimate srad and humidity with iterative algorithm
! without Tdew input data, an iterative estimation of shortwave radiation and humidity is required
call calc_srad_humidity_iterative( &
     yday, s_tmax, s_tmin, s_tday, s_prcp, s_swe, &
     s_hum, s_srad, s_srad_dif ,s_lrad, s_lradnet, s_dayl, &
     ndays=ndays, outhum=outhum, lwrad=lwrad, netlwrad=netlwrad, &
     site_lat=site_lat, site_elev=site_elev, site_asp=site_asp, site_slp=site_slp, &
     site_ehoriz=site_ehoriz, site_whoriz=site_whoriz,&
     TBASE=TBASE, ABASE=ABASE, C=C, B0=B0, B1=B1, B2=B2, &
     RAIN_SCALAR=RAIN_SCALAR, DIF_ALB=DIF_ALB, LWAVE_COR=LWAVE_COR)

print *, "Completed radiation algorithm"
!----------------------------------------
! write output file
open(unit=1234, file=outprefix, status='unknown', action='write')
write(1234, '(A3,A9,A9,A9,A9,A9,A9,A9,A9,A9,A9)') "yd", "tmax", "tmin", "tday", "prcp", "VPD", &
     "swrad", "swrad_dif", "lwrad", "lwradnet", "daylen"
write(1234, '(A3,A9,A9,A9,A9,A9,A9,A9,A9,A9,A9)') " ","(deg C)", "(deg C)", "(deg C)", "(cm)", "(Pa)", &
     "(W m-2)", "(W m-2)", "(W m-2)", "(W m-2)", "(s)"
   do i=1, ndays
      write(1234, '(i3,f9.2,f9.2,f9.2,f9.2,f9.2,f9.2,f9.2,f9.2,f9.2,f9.2)') yday(i), s_tmax(i), s_tmin(i), s_tday(i),&
           s_prcp(i), s_hum(i), s_srad(i), s_srad_dif(i), s_lrad(i), s_lradnet(i), s_dayl(i)
   end do
close(unit=1234)
print *, "Results written in ", outprefix

!----------------------------------------
! test mtclim eval
!call eval_mtclim (parameterset, yday, tmin, tmax, prcp, ndays, &
!     outhum, lwrad, netlwrad &
!     s_tmax, s_tmin, s_tday, s_prcp, &
!     s_hum, s_srad, s_srad_dif, s_lrad, s_lradnet, s_dayl, &
!     site_elev, base_elev, tmin_lr, tmax_lr, site_isoh, base_isoh, &
!     site_lat, site_asp, site_slp, site_ehoriz, site_whoriz)

!----------------------------------------
! test mtclim objective swrad
!if (observ .eq. 1_i4) then
!   print *, "RMSE is ", objective_swrd(parameterset)
!end if
!
!----------------------------------------
! free previously allocated memory before returning
!
deallocate(yday)
deallocate(tmax)
deallocate(tmin)
deallocate(prcp)
!
deallocate(s_swe)
!
deallocate(s_tmin)
deallocate(s_tmax)
deallocate(s_tday)
!
deallocate(s_hum)
deallocate(s_srad)
deallocate(s_srad_dif)
deallocate(s_lrad)
deallocate(s_lradnet)
deallocate(s_dayl)
!
!------------------------------------------
isgood = .true.
! Finish
allgood = allgood .and. isgood
  if (allgood) then
     write(*,*) 'mo_mtclim o.k.'
  else
     write(*,*) 'mo_mtclim failed!'
  endif
!
end program main

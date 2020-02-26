PROGRAM mo_bootstrapping_test

  USE mo_bootstrapping_sensitivity_analysis,  only: bootstrap_si
  USE mo_kind,                                only: i4, i8, dp
  USE mo_model,                               only: model, par_range
  USE mo_sobol,                               only: sobol_array
  USE mo_sobol_index,                         only: sobol_index
  use mo_ansi_colors, only: color, c_red, c_green

  IMPLICIT NONE

!-----------------------------------------------------------

  integer(i4),   parameter                   :: npar = 4_i4             ! number of parameters
  integer(i4),   parameter                   :: nsets = 1500_i4         ! number of parameter sets
  integer(i4),   parameter                   :: nt = 10_i4              ! number of time points
  integer(i4),   parameter                   :: nboot = 19_i4           ! number of new datasets generated with bootstrapping

  integer(i4)                                :: skip = 30000_i4       
  integer(i8)                                :: seed = 123_i8

  real(dp),  dimension(nt)                   :: t                       ! i.e. time series

  real(dp),  dimension(npar,2)               :: range                   ! parameter range

  real(dp),  dimension(nsets,2*npar)         :: p                       ! sobol sequence
  real(dp),  dimension(nsets,npar)           :: pa                      ! parameter sample A
  real(dp),  dimension(nsets,npar)           :: pb                      ! parameter sample B
  real(dp),  dimension(nsets,npar)           :: pci                     ! parameter sample Ci

  real(dp),  dimension(nsets,nt)             :: ya                      ! model output A
  real(dp),  dimension(nsets,nt)             :: yb                      ! model output B
  real(dp),  dimension(nsets,npar,nt)        :: yci                     ! model output C(i), i=1,npar

  real(dp),  dimension(nboot+1,npar,nt)      :: si                      ! sobol index (main effect) 
                                                                        !    per parameter and time point
                                                                        !    for each bootstrapping sample
  real(dp),  dimension(nboot+1,npar,nt)      :: sti                     ! sobol index (total effect)
                                                                        !    per parameter and time point
                                                                        !    for each bootstrapping sample
  integer(i4)        :: i
  logical            :: isgood = .true.

!----------------------------------------------------------

  do i = 1, nt
     t(i) = real(i-1,dp)
  end do

  ! Generating parameter sets 
  call sobol_array(2*npar,nsets,skip,p)
  range = par_range()
  pa = p(:,:npar)
  pb = p(:,npar+1:)
  do i = 1, npar
     pa(:,i) = pa(:,i)*(range(i,2)-range(i,1)) + range(i,1)
     pb(:,i) = pb(:,i)*(range(i,2)-range(i,1)) + range(i,1)
  end do

  ! Calculating the model output
  ya = model(pa,t)
  yb = model(pb,t)

  do i = 1, npar
     pci = pa
     pci(:,i) = pb(:,i)                                           
     yci(:,i,:) = model(pci,t)
  end do

  ! Calculating the sobol index
  call sobol_index(ya, yb, yci, si(1,:,:), sti(1,:,:))

  ! Bootstrapping
  call bootstrap_si(ya, yb, yci, nboot, si(2:,:,:), sti(2:,:,:), seed = seed)


 if (abs(si(1,1,3) - 0.97236268) > 0.000001) isgood = .false.
 if (abs(si(1,2,3) - 0.01255484) > 0.000001) isgood = .false.
 if (abs(si(1,3,3) - 0.00000000) > 0.000001) isgood = .false.
 if (abs(si(1,4,3) - 0.00275693) > 0.000001) isgood = .false.
 if (abs(si(2,1,3) - 0.82035742) > 0.000001) isgood = .false.
 if (abs(si(2,2,3) - 0.00413255) > 0.000001) isgood = .false.
 if (abs(si(2,3,3) - 0.00000000) > 0.000001) isgood = .false.
 if (abs(si(2,4,3) - 0.01841593) > 0.000001) isgood = .false.
 if (abs(sti(1,1,3) - 0.97907962) > 0.000001) isgood = .false.
 if (abs(sti(1,2,3) - 0.01529379) > 0.000001) isgood = .false.
 if (abs(sti(1,3,3) - 0.00000000) > 0.000001) isgood = .false.
 if (abs(sti(1,4,3) - 0.00382390) > 0.000001) isgood = .false.
 if (abs(sti(2,1,3) - 1.01645796) > 0.000001) isgood = .false.
 if (abs(sti(2,2,3) - 0.01582186) > 0.000001) isgood = .false.
 if (abs(sti(2,3,3) - 0.00000000) > 0.000001) isgood = .false.
 if (abs(sti(2,4,3) - 0.00392133) > 0.000001) isgood = .false.

  if (isgood) then
     write(*,*) 'mo_bootstrapping_sensitivity_analysis ', color('o.k.', c_green)
  else
     write(*,*) 'mo_bootstrapping_sensitivity_analysis ', color('failed', c_red)
  end if

!-----------------------------------------------------------
END PROGRAM mo_bootstrapping_test

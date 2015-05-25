PROGRAM main

  ! \ls mo_*.f90 | sed -e 's/mo_/  USE mo_/' -e 's/.f90//' | \
  !    sed -e 's/USE mo_minpack/!USE mo_minpack/' -e 's/USE mo_template/!USE mo_template/' \
  !        -e 's/USE mo_nr$/!USE mo_nr/' -e 's/USE mo_pumpingtests/!USE mo_pumpingtests/'
  ! make EXCLUDE_FILES='mo_minpack.f90 mo_nr.f90 mo_pumpingtests.f90 mo_template.f90' 
  USE mo_anneal
  USE mo_append
  USE mo_boxcox
  USE mo_combinatorics
  USE mo_constants
  USE mo_corr
  USE mo_dds
  USE mo_delsa
  USE mo_elemeffects
  USE mo_errormeasures
  USE mo_file_utils
  USE mo_finish
  USE mo_fit
  USE mo_functions
  USE mo_histo
  USE mo_integrate
  USE mo_interpol
  USE mo_julian
  USE mo_kernel
  USE mo_kind
  USE mo_laplace_inversion
  USE mo_linear_algebra
  USE mo_linfit
  USE mo_mad
  USE mo_mcmc
  USE mo_message
  !USE mo_minpack
  USE mo_moment
  USE mo_ncread
  USE mo_ncwrite
  USE mo_nelmin
  USE mo_nml
  !USE mo_nr
  USE mo_nrutil
  USE mo_ode_generator
  USE mo_ode_solver
  USE mo_opt_functions
  USE mo_orderpack
  USE mo_percentile
  USE mo_pi_index
  USE mo_poly
  !USE mo_pumpingtests
  USE mo_quicksort
  USE mo_random_field
  USE mo_remap
  USE mo_sampling
  USE mo_sce
  USE mo_sobol
  USE mo_sobol_index
  USE mo_sort
  USE mo_spatialsimilarity
  USE mo_specan
  USE mo_spline
  USE mo_string_utils
  !USE mo_template
  USE mo_timer
  USE mo_utils
  USE mo_xor4096
  USE mo_xor4096_apps

  IMPLICIT NONE

  Write(*,*) 'I can compile.'

END PROGRAM main

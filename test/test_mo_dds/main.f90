program main

  use mo_kind, only: i4, i8, dp
  use mo_dds,  only: dds
  use mo_opt_functions, only: griewank

  implicit none

  integer(i4), parameter :: np = 10
  real(dp),    dimension(np) :: dv_ini  ! inital value of decision variables 
  real(dp),    dimension(np) :: dv_min  ! Min value of decision variables 
  real(dp),    dimension(np) :: dv_max  ! Max value of decision variables 
  real(dp)                   :: r_val   ! DDS perturbation parameter (-> 0.2 by default)
  integer(i4) :: seed     ! User seed to initialise the random number generator
  integer(i8) :: nIterMax ! Maximum number of iteration or function evaluation
  real(dp) :: to_max      ! Maximization or minimization of function
  integer(i8) :: nobj     ! Counter for Obj. function evaluations
  real(dp), dimension(np) :: dv_opt ! Best value of decision variables

  write(*,*) '**************************************************'
  write(*,*) '  Dynamically Dimensioned Search (DDS) Algorithm  '
  write(*,*) '         version 1.1 - by Bryan Tolson            '
  write(*,*) '**************************************************'

  nIterMax = 100000_i8
  r_val    = 0.20_dp
  to_max   = 1.0_dp
  seed     = 123456789_i4
  dv_min   = (/ -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0, -600.0 /)
  dv_max   = (/ 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0 /)
  dv_ini   = (/ -.226265E+01, -.130187E+01, -.151219E+01, 0.133983E+00, 0.988159E+00, &
       -.495074E+01, -.126574E+02, 0.572684E+00, 0.303864E+01, 0.343031E+01 /)
  dv_opt   = DDS(griewank, dv_ini, dv_min, dv_max, r_val, seed, nIterMax, to_max, nobj)

  write(*,*) ''
  write(*,'(A,10F9.4)') 'gFortran should be ', -3.1363960365479695, 4.4197167350424271, &
       4.09472260676029265E-002, -1.18042790202430650E-002, -6.79918707079235629E-003, &
       7.6358557562133633, 8.2796397347388364, -7.31881942986879475E-002, &
       3.52734183369751778E-002, -3.59188289890255974E-002
  write(*,'(A,10F9.4)') 'Output             ', dv_opt

  write(*,*) 'mo_dds double precision o.k.'

end program main

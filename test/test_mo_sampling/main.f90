! =============================================================================
! Name        : main_mo_sampling.f90
! Author      : Giovanni Dalmasso
! Contacts	  : giovanni.dalmasso@ufz.de - http://www.ufz.de/index.php?en=30665
!               Department of Computational Hydrosystems CHS - UZF Leipzig
! Date	      : Apr 16, 2013
! Version     : 1.0
! Copyright   : 2012-2013 Giovanni Dalmasso
! Description : main for testing mo_sampling.f90
! =============================================================================

program main_mo_sampling

    use mo_kind,            only    :   i4, sp
    use mo_ansi_colors, only: color, c_red, c_green
    use mo_sampling,        only    :   random, lhs, setupxor
    use mo_xor4096,         only    :   n_save_state

    implicit none

    integer(i4),    dimension(:),   allocatable ::  seed                ! seed
    integer(i4),    dimension(:,:), allocatable ::  save_the_setup_rnd  ! xor variable rnd
    integer(i4),    dimension(:,:), allocatable ::  save_the_setup_lhs  ! xor variable lhs
    real(sp)                                    ::  test                ! test variable
    real(sp),       dimension(:,:), allocatable ::  randomTest          ! random sampling
    real(sp),       dimension(:,:), allocatable ::  lhsTest             ! lhs sampling
    logical                                     ::  isgood              ! test variable
    logical,        dimension(:),   allocatable ::  dist                ! flag to chose between unif or gauss

    integer(i4),    parameter                   ::  Nsample = 10000_i4  ! number of samples
    integer(i4),    parameter                   ::  Npar = 3_i4         ! number of parameters

    if ( allocated(randomTest) ) deallocate(randomTest)
    if ( allocated(lhsTest) ) deallocate(lhsTest)
    if ( allocated(save_the_setup_rnd) ) deallocate(save_the_setup_rnd)
    if ( allocated(save_the_setup_lhs) ) deallocate(save_the_setup_lhs)
    if ( allocated(dist) ) deallocate(dist)
    if ( allocated(seed) ) deallocate(seed)
    allocate( randomTest(Nsample, Npar) )
    allocate( lhsTest(Nsample, Npar) )
    allocate( save_the_setup_rnd(Npar, n_save_state) )
    allocate( save_the_setup_lhs(Npar, n_save_state) )
    allocate( dist(Npar) )
    allocate( seed(Npar) )

    isgood = .true.  ! test variable
    test = 1.0e6_sp  ! test variable

    ! distributions
    dist = (/ .false., .true., .false. /)
    ! seed
    seed = (/ 100_i4, 200_i4, 300_i4 /)

    ! setup xor4096 with distributions
    call setupxor( Npar, dist, save_the_setup_rnd, seed=seed )
    ! test random sampling
    randomTest = random( Nsample, Npar, save_the_setup_rnd, dist=dist )
    ! isgood test
    if ( any(nint(randomTest(1,:)*test) .ne. (/ 641292_i4, 227264_i4, 161895_i4 /)) ) isgood = .false.    ! test

    ! setup xor4096
    call setupxor( Npar, save_the_setup_lhs, seed=seed )
    ! test Latin Hypercube Sampling
    lhsTest = lhs( Nsample, Npar, save_the_setup_lhs, dist=dist )
    print*, 'Ha01 ', lhsTest(1,:)*test
    if ( any(nint(lhsTest(1,:)*test) .ne. (/ 634054_i4, 736906_i4, 633938_i4 /)) ) isgood = .false.    ! test

    ! Is the program running properly?
    write(*,*) ''
    write(*,*) '------------------------------------------'
    if (isgood) then
        write(*,*) 'mo_sampling.f90 ', color('o.k.', c_green)
    else
        write(*,*) 'mo_sampling.f90 ', color('failed!', c_red)
    endif
    write(*,*) '-------------------END---------------------'

end program main_mo_sampling

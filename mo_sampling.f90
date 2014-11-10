!> \file mo_sampling.f90

!> \brief Random and Latin Hypercube Sampling.

!> \details Random and Latin Hypercube Sampling for a set of parameters with Uniform(0,1) or Gaussian(0,1) Distribution.

!> \authors Giovanni Dalmasso
!> \date Apr 2013

module mo_sampling

  ! This module is a template for the UFZ CHS Fortran library.

  ! Written  Giovanni Dalmasso, Apr 2013

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! The UFZ Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The UFZ Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the UFZ Fortran library (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2012-2013 Giovanni Dalmasso

  implicit none

  public :: random        ! Random Sampling
  public :: lhs           ! Latin Hypercube Sampling
  public :: setupxor      ! Initialize xor4096.f90 (uniform distribution)
  public :: setupxor_g    ! Initialize xor4096.f90 (gaussian distribution)

  ! ------------------------------------------------------------------

  !     NAME
  !         random

  !     PURPOSE
  !         Calculate the Random Sampling for a set of parameters with Uniform(0,1) or Gaussian(0,1) Distribution.
  !
  !>        \brief Random Sampling.
  !
  !>        \details Calculate the Random Sampling for a set of $N$ parameters with Uniform(0,1) or Gaussian(0,1) Distribution.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4)                                      :: Nsample"             number of samples
  !>        \param[in] "integer(i4)                                      :: Npar"                number of parameters
  !
  !     INTENT(INOUT)
  !>        \param[in] "integer(i4/i8),  dimension(Npar, n_save_state)   :: save_the_setup(:,:)"  xor4096 variable
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical,      dimension(Npar),    optional       :: dist(:)"             distributions
  !>                                                                                          if true = Uniform
  !>                                                                                          if false = Gaussian
  !>                                                                                          DEFAULT: true
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>       \return     real(sp/dp),       dimension(Nsample, Npar)       :: randomSamp(:,:)       Random Sampling
  !
  !     RESTRICTIONS
  !>       \note The random number generator has to be initialized before and therefore is called without a seed.
  !
  !     EXAMPLE
  !         -> see also example in test directory
  !
  !     LITERATURE
  !         none
  !
  !     HISTORY
  !>        \author Giovanni Dalmasso
  !>        \date Apr 2013
  !         Modified,
  !
  interface random
     module procedure random_sp, random_dp
  end interface random

  ! ------------------------------------------------------------------

  !     NAME
  !         lhs

  !     PURPOSE
  !         Calculate the Latin Hypercube Sampling for a set of parameters with Uniform(0,1) or Gaussian(0,1) Distribution.
  !
  !>        \brief Random Sampling.
  !
  !>        \details Calculate the Latin Hypercube Sampling for a set of $N$ parameters with Uniform(0,1)
  !>                  or Gaussian(0,1) Distribution.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4)                                      :: Nsample"             number of samples
  !>        \param[in] "integer(i4)                                      :: Npar"                number of parameters
  !
  !     INTENT(INOUT)
  !>        \param[in] "integer(i4/i8),  dimension(Npar, n_save_state)   :: save_the_setup(:,:)"  xor4096 variable
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical,      dimension(Npar),    optional       :: dist(:)"             distributions
  !>                                                                                          if true = Uniform
  !>                                                                                          if false = Gaussian
  !>                                                                                          DEFAULT: true
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>       \return     real(sp/dp),       dimension(Nsample, Npar)       :: lhsSamp(:,:)          Latin Hypercube Sampling
  !
  !     RESTRICTIONS
  !>       \note The random number generator has to be initialized before and therefore is called without a seed.
  !>             There is an i4 stream initialized inside the method.
  !>             Fortran90 implementation adapted from Chad Sprouse's C version
  !>             (http://home.online.no/~pjacklam/notes/invnorm/impl/sprouse/ltqnorm.c)
  !
  !     EXAMPLE
  !         -> see also example in test directory
  !
  !     LITERATURE
  !         Stein, M. 1987. Large Sample Properties of Simulations Using Latin Hypercube Sampling. Technometrics 29:143-151
  !
  !     HISTORY
  !>        \author Giovanni Dalmasso
  !>        \date Apr 2013
  !         Modified,
  !
  interface lhs
     module procedure lhs_sp, lhs_dp
  end interface lhs

  ! ------------------------------------------------------------------

  !     NAME
  !         setupxor

  !     PURPOSE
  !         Set up the module mo_xor4096.f90 with uniform distribution.
  !
  !>        \brief Random Sampling.
  !
  !>        \details The random number generator has to be initialized before using.
  !
  !     INTENT(IN)
  !>         \param[in] "integer(i4)                                              ::    sizeArray"            number of samples
  !>         \param[in] "logical,      dimension(sizeArray)                       ::    distribution"         distributions
  !>                                                                                                      if true = Uniform
  !>                                                                                                      if false = Gaussian
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !!>        \param[out] "integer(i4/i8),    dimension(sizeArray, n_save_state)    ::    save_the_setup(:,:)    
  !
  !     INTENT(IN), OPTIONAL
  !>          \param[in] "integer(i4/i8),    dimension(sizeArray),   optional      ::    seed(:)"              seed
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !
  !     RESTRICTIONS
  !>       \note The intent in variable sizeArray could be omitted if the problem is scalar (dim=0).
  !>              Therefore save_the_setup will have only one dimension equal to sizeArray and the seed will be a scalar.
  !>              If distribution is omitted the default is uniform.
  !
  !     EXAMPLE
  !         -> see also example in test directory
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Giovanni Dalmasso
  !>        \date Apr 2013
  !         Modified,
  !
  interface setupxor
     module procedure setupxor4096_i4_0d, &
          setupxor4096_sp_1d, setupxor4096_dp_1d, &
          setupxor4096_sp_1d_dist, setupxor4096_dp_1d_dist
  end interface setupxor

  ! ------------------------------------------------------------------

  !     NAME
  !         setupxor_g

  !     PURPOSE
  !         Set up the module mo_xor4096.f90 with gaussian distribution.
  !
  !>        \brief Random Sampling.
  !
  !>        \details The random number generator has to be initialized before using.
  !
  !     INTENT(IN)
  !>         \param[in] "integer(i4)                                              ::    sizeArray"            number of samples
  !>         \param[in] "logical,      dimension(sizeArray)                       ::    distribution"         distributions
  !>                                                                                                      if true = Uniform
  !>                                                                                                      if false = Gaussian
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !!>        \param[out] "integer(i4/i8),    dimension(sizeArray, n_save_state)    ::    save_the_setup(:,:)    
  !
  !     INTENT(IN), OPTIONAL
  !>          \param[in] "integer(i4/i8),    dimension(sizeArray),   optional      ::    seed(:)"              seed
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !
  !     RESTRICTIONS
  !>       \note The intent in variable sizeArray could be omitted if the problem is scalar (dim=0).
  !>              Therefore save_the_setup will have only one dimension equal to sizeArray and the seed will be a scalar.
  !>              If distribution is omitted the default is uniform.
  !
  !     EXAMPLE
  !         -> see also example in test directory
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Giovanni Dalmasso
  !>        \date Apr 2013
  !         Modified,
  !
  interface setupxor_g
     module procedure setupxor4096g_sp_1d, setupxor4096g_dp_1d
  end interface setupxor_g


  private

  ! ------------------------------------------------------------------
  ! Private Method
  ! Lower tail quantile for standard normal distribution function
  interface ltqnorm
     module procedure &
          ltqnorm_sp_d1, ltqnorm_dp_d1 
          ! ltqnorm_sp_d0, ltqnorm_dp_d0, 
  end interface ltqnorm

contains

  ! Random Sampling SINGLE PRECISION
  function random_sp( Nsample, Npar, save_the_setup, dist )    result(randomSamp)

    use mo_kind,            only    :   i4, sp
    use mo_xor4096,         only    :   xor4096, xor4096g, n_save_state

    implicit none

    ! Intent IN
    integer(i4),                                    intent(in)      ::  Nsample             ! number of samples
    integer(i4),                                    intent(in)      ::  Npar                ! number of parameters
    logical,        dimension(Npar),    optional,   intent(in)      ::  dist                ! distribution (unif, norm)

    ! Intent InOut
    integer(i4),    dimension(Npar, n_save_state),  intent(inout)   ::  save_the_setup      ! xor4096 variable

    ! Results
    real(sp),       dimension(Nsample, Npar)                        ::  randomSamp          ! random sampling

    ! Internal variables
    integer(i4)                         ::  j, jj   ! counters
    logical,        dimension(Npar)     ::  local_dist

    if ( Npar .eq. 0_i4 )   stop 'random_sp --> number of parameters could not be zero!!!'

    if( present(dist) ) then
       local_dist = dist
    else
       local_dist = .true.
    end if

    do j=1, Nsample
       do jj=1, Npar
          if ( local_dist(jj) )   then
             call xor4096( 0_i4, randomSamp(j,jj), save_state=save_the_setup(jj,:) )
          else
             call xor4096g( 0_i4, randomSamp(j,jj), save_state=save_the_setup(jj,:) )
          end if
       end do
    end do

  end function random_sp


  ! Random Sampling DOUBLE PRECISION
  function random_dp( Nsample, Npar, save_the_setup, dist )    result(randomSamp)

    use mo_kind,            only    :   i4, i8, dp
    use mo_xor4096,         only    :   xor4096, xor4096g, n_save_state

    implicit none

    ! Intent IN
    integer(i4),                                    intent(in)      ::  Nsample             ! number of samples
    integer(i4),                                    intent(in)      ::  Npar                ! number of parameters
    logical,        dimension(Npar),    optional,   intent(in)      ::  dist                ! distribution (unif, norm)

    ! Intent InOut
    integer(i8),    dimension(Npar, n_save_state),  intent(inout)   ::  save_the_setup      ! xor4096 variable

    ! Results
    real(dp),       dimension(Nsample, Npar)                        ::  randomSamp          ! random sampling

    ! Internal variables
    integer(i4)                         ::  j, jj   ! counters
    logical,        dimension(Npar)     ::  local_dist

    if ( Npar .eq. 0_i4 )   stop 'random_dp --> number of parameters could not be zero!!!'

    if( present(dist) ) then
       local_dist = dist
    else
       local_dist = .true.
    end if

    do j=1, Nsample
       do jj=1, Npar
          if ( local_dist(jj) )   then
             call xor4096( 0_i8, randomSamp(j,jj), save_state=save_the_setup(jj,:) )
          else
             call xor4096g( 0_i8, randomSamp(j,jj), save_state=save_the_setup(jj,:) )
          end if
       end do
    end do

  end function random_dp


  ! Latin Hypercube Sampling SINGLE PRECISION
  function lhs_sp( Nsample, Npar, save_the_setup, dist )    result(lhsSamp)

    use mo_kind,            only    :   i4, sp
    use mo_xor4096,         only    :   n_save_state
    use mo_xor4096_apps,    only    :   xor4096_array
    use mo_combinatorics,   only    :   random_index_permut

    implicit none

    ! Intent IN
    integer(i4),                                    intent(in)      ::  Nsample         ! number of samples
    integer(i4),                                    intent(in)      ::  Npar            ! number of parameters
    logical,        dimension(Npar),    optional,   intent(in)      ::  dist            ! distribution (unif, norm)

    ! Intent InOut
    integer(i4),    dimension(Npar, n_save_state),  intent(inout)   ::  save_the_setup  ! xor4096 variable

    ! Results
    real(sp),       dimension(Nsample, Npar)                        ::  lhsSamp         ! random sampling

    ! Internal variables
    integer(i4)                                 ::  j                   ! counter
    integer(i4)                                 ::  seed                ! seed
    integer(i4),    dimension(n_save_state)     ::  save_the_setup_i4   ! setup xor4096 i4 stream
    integer(i4),    dimension(Nsample)          ::  idx
    real(sp),       dimension(Nsample, Npar)    ::  ranSam              ! random sampling
    real(sp),       dimension(Nsample)          ::  p
    logical,        dimension(Npar)             ::  local_dist

    if ( Npar .eq. 0_i4 )   stop 'lhs_sp --> number of parameters could not be zero!!!'

    if( present(dist) ) then
       local_dist = dist
    else
       local_dist = .true.
    end if

    ! random sampling matrix
    call xor4096_array( ranSam, save_state=save_the_setup )

    ! seed for setup xor4096 i4 stream
    if ( save_the_setup(1,1) .gt. huge(1_i4)/save_the_setup(1,2) )  then
       seed = save_the_setup(3,1)
    else
       seed = save_the_setup(1,1)*save_the_setup(1,2) + save_the_setup(1,3)
    end if
    ! setup xor4096 i4 stream
    call setupxor( save_the_setup_i4, seed=seed )

    ! Latin Hypercube Sampling
    idx = random_index_permut( Nsample, save_the_setup_i4 )
    do j=1, Npar
       p = ( real(idx,sp) - ranSam(:,j) )/real(Nsample,sp)
       if( local_dist(j) )    then
          lhsSamp(:,j) = p                ! --> unif
       else
          lhsSamp(:,j) = ltqnorm(p(:))    ! --> gauss
       end if
    end do

  end function lhs_sp


  ! Latin Hypercube Sampling DOUBLE PRECISION
  function lhs_dp( Nsample, Npar, save_the_setup, dist )    result(lhsSamp)

    use mo_kind,            only    :   i4, i8, dp
    use mo_xor4096,         only    :   n_save_state
    use mo_xor4096_apps,    only    :   xor4096_array
    use mo_combinatorics,   only    :   random_index_permut

    implicit none

    ! Intent IN
    integer(i4),                                    intent(in)      ::  Nsample         ! number of samples
    integer(i4),                                    intent(in)      ::  Npar            ! number of parameters
    logical,        dimension(Npar),    optional,   intent(in)      ::  dist            ! distribution (unif, norm)

    ! Intent InOut
    integer(i8),    dimension(Npar, n_save_state),  intent(inout)   ::  save_the_setup  ! xor4096 variable

    ! Results
    real(dp),       dimension(Nsample, Npar)                        ::  lhsSamp         ! random sampling

    ! Internal variables
    integer(i4)                                 ::  j                   ! counter
    integer(i4)                                 ::  seed                ! seed
    integer(i4),    dimension(n_save_state)     ::  save_the_setup_i4   ! setup xor4096 i4 stream
    integer(i4),    dimension(Nsample)          ::  idx
    real(dp),       dimension(Nsample, Npar)    ::  ranSam              ! random sampling
    real(dp),       dimension(Nsample)          ::  p
    logical,        dimension(Npar)             ::  local_dist

    if ( Npar .eq. 0_i4 )   stop 'lhs_dp --> number of parameters could not be zero!!!'

    if ( present(dist) ) then
       local_dist = dist
    else
       local_dist = .true.
    end if

    ! random sampling matrix
    call xor4096_array( ranSam, save_state=save_the_setup )

    ! seed for setup xor4096 i4 stream
    if ( save_the_setup(1,1) .gt. huge(1_i4)/save_the_setup(1,2) )  then
       seed = int( save_the_setup(1,3),i4 )
    else
       seed = int( save_the_setup(1,1)*save_the_setup(1,2) + save_the_setup(1,3), i4 )
    end if
    ! setup xor4096 i4 stream
    call setupxor( save_the_setup_i4, seed=seed )

    ! Latin Hypercube Sampling
    idx = random_index_permut( Nsample, save_the_setup_i4 )
    do j=1, Npar
       p = ( real(idx,dp) - ranSam(:,j) )/real(Nsample,dp)
       if( local_dist(j) )    then
          lhsSamp(:,j) = p                ! --> unif
       else
          lhsSamp(:,j) = ltqnorm(p(:))    ! --> gauss
       end if
    end do

  end function lhs_dp

  ! ------------------------------------------------------------------

  ! SETUP XOR

  subroutine setupxor4096_i4_0d( save_the_setup, seed )

    use mo_kind,    only : i4
    use mo_xor4096, only : get_timeseed, xor4096, n_save_state

    implicit none

    ! Intent OUT
    integer(i4),    dimension(n_save_state),    intent(out) ::  save_the_setup

    ! Intent IN optional
    integer(i4),            optional,           intent(in)  ::  seed

    ! Internal variables
    integer(i4)             :: seedin
    integer(i4)             :: dummyRn

    if ( present(seed) ) then
       seedin = seed
    else
       call get_timeseed( seedin )
    end if

    call xor4096( seedin, dummyRn, save_state=save_the_setup )
    seedin = 0_i4

  end subroutine setupxor4096_i4_0d


  subroutine setupxor4096_sp_1d( sizeArray, save_the_setup, seed )

    use mo_kind,    only : i4, sp
    use mo_xor4096, only : get_timeseed, xor4096, n_save_state

    implicit none

    ! Intent IN
    integer(i4),                                        intent(in)  ::  sizeArray

    ! Intent OUT
    integer(i4),    dimension(sizeArray, n_save_state), intent(out) ::  save_the_setup

    ! Intent IN optional
    integer(i4),    dimension(sizeArray),   optional,   intent(in)  ::  seed

    ! Internal variables
    integer(i4),    dimension(sizeArray)  :: seedin
    real(sp),       dimension(sizeArray)  :: dummyRn

    if ( present(seed) ) then
       seedin(:) = seed(:)
    else
       call get_timeseed( seedin )
    end if

    call xor4096( seedin, dummyRn, save_state=save_the_setup )
    seedin(:) = 0_i4

  end subroutine setupxor4096_sp_1d


  subroutine setupxor4096_dp_1d( sizeArray, save_the_setup, seed )

    use mo_kind,    only : i4, i8, dp
    use mo_xor4096, only : get_timeseed, xor4096, n_save_state

    implicit none

    ! Intent IN
    integer(i4),                                        intent(in)  ::  sizeArray

    ! Intent OUT
    integer(i8),    dimension(sizeArray, n_save_state), intent(out) ::  save_the_setup

    ! Intent IN optional
    integer(i8),    dimension(sizeArray),   optional,   intent(in)  ::  seed

    ! Internal variables
    integer(i8),    dimension(sizeArray)  :: seedin
    real(dp),       dimension(sizeArray)  :: dummyRn

    if ( present(seed) ) then
       seedin(:) = seed(:)
    else
       call get_timeseed( seedin )
    end if

    call xor4096( seedin, dummyRn, save_state=save_the_setup )
    seedin(:) = 0_i8

  end subroutine setupxor4096_dp_1d


  subroutine setupxor4096g_sp_1d( sizeArray, save_the_setup_g, seed )

    use mo_kind,    only : i4, sp
    use mo_xor4096, only : get_timeseed, xor4096g, n_save_state

    implicit none

    ! Intent IN
    integer(i4),                                        intent(in)  :: sizeArray

    ! Intent OUT
    integer(i4),    dimension(sizeArray, n_save_state), intent(out) :: save_the_setup_g

    ! Intent IN optional
    integer(i4),    dimension(sizeArray),   optional,   intent(in)  ::  seed

    ! Internal variables
    integer(i4),    dimension(sizeArray)  :: seedin
    real(sp),       dimension(sizeArray)  :: dummyRn

    if ( present(seed) ) then
       seedin(:) = seed(:)
    else
       call get_timeseed( seedin )
    end if

    call xor4096g( seedin, dummyRn, save_state=save_the_setup_g )
    seedin(:) = 0_i4

  end subroutine setupxor4096g_sp_1d


  subroutine setupxor4096g_dp_1d( sizeArray, save_the_setup_g, seed )

    use mo_kind,    only : i4, i8, dp
    use mo_xor4096, only : get_timeseed, xor4096g, n_save_state

    implicit none

    ! Intent IN
    integer(i4),                                        intent(in)  :: sizeArray

    ! Intent OUT
    integer(i8),    dimension(sizeArray, n_save_state), intent(out) :: save_the_setup_g

    ! Intent IN optional
    integer(i8),    dimension(sizeArray),   optional,   intent(in)  ::  seed

    ! Internal variables
    integer(i8),    dimension(sizeArray)  :: seedin
    real(dp),       dimension(sizeArray)  :: dummyRn

    if ( present(seed) ) then
       seedin(:) = seed(:)
    else
       call get_timeseed( seedin )
    end if

    call xor4096g( seedin, dummyRn, save_state=save_the_setup_g )
    seedin(:) = 0_i8

  end subroutine setupxor4096g_dp_1d


  subroutine setupxor4096_sp_1d_dist( sizeArray, distribution, save_the_setup, seed )

    use mo_kind,    only : i4
    use mo_xor4096, only : n_save_state

    implicit none

    ! Intent IN
    integer(i4),                                        intent(in)  :: sizeArray
    logical,        dimension(sizeArray),               intent(in)  :: distribution

    ! Intent OUT
    integer(i4),    dimension(sizeArray, n_save_state), intent(out) :: save_the_setup

    ! Intent IN optional
    integer(i4),    dimension(sizeArray),   optional,   intent(in)  ::  seed

    ! Internal variables
    integer(i4)                                 ::  j                   ! counter
    integer(i4)                                 ::  Nunif, Ngauss       ! number of uniform and gaussian
    integer(i4),    dimension(:,:), allocatable ::  save_the_setup_u    ! xor variable unif
    integer(i4),    dimension(:,:), allocatable ::  save_the_setup_g    ! xor variable gauss


    Nunif = count(distribution)
    Ngauss = sizeArray - Nunif

    if ( allocated(save_the_setup_u) ) deallocate(save_the_setup_u)
    if ( allocated(save_the_setup_g) ) deallocate(save_the_setup_g)
    if ( Nunif .ne. 0_i4 )  allocate( save_the_setup_u(Nunif, n_save_state) )
    if ( Ngauss .ne. 0_i4 ) allocate( save_the_setup_g(Ngauss, n_save_state) )
    save_the_setup_u = 0_i4
    save_the_setup_g = 0_i4

    ! setup xor4096 - xor4096g
    if ( present(seed) ) then
       if ( Nunif .ne. 0_i4 )  call setupxor( Nunif, save_the_setup_u, seed=seed )
       if ( Ngauss .ne. 0_i4 ) call setupxor_g( Ngauss, save_the_setup_g, seed=seed )
    else
       if ( Nunif .ne. 0_i4 )  call setupxor( Nunif, save_the_setup_u )
       if ( Ngauss .ne. 0_i4 ) call setupxor_g( Ngauss, save_the_setup_g )
    end if

    do j=1,sizeArray
       if ( distribution(j) ) then
          save_the_setup(j,:) = save_the_setup_u(count(distribution(1:j)),:)
       else
          save_the_setup(j,:) = save_the_setup_g(count(.not.distribution(1:j)),:)
       end if
    end do

  end subroutine setupxor4096_sp_1d_dist


  subroutine setupxor4096_dp_1d_dist( sizeArray, distribution, save_the_setup, seed )

    use mo_kind,    only : i4, i8
    use mo_xor4096, only : n_save_state

    implicit none

    ! Intent IN
    integer(i4),                                        intent(in)  :: sizeArray
    logical,        dimension(sizeArray),               intent(in)  :: distribution

    ! Intent OUT
    integer(i8),    dimension(sizeArray, n_save_state), intent(out) :: save_the_setup

    ! Intent IN optional
    integer(i8),    dimension(sizeArray),   optional,   intent(in)  ::  seed

    ! Internal variables
    integer(i4)                                 ::  j, counter          ! counter
    integer(i4)                                 ::  Nunif, Ngauss       ! number of uniform and gaussian
    integer(i8),    dimension(:,:), allocatable ::  save_the_setup_u    ! xor variable unif
    integer(i8),    dimension(:,:), allocatable ::  save_the_setup_g    ! xor variable gauss
    integer(i8),    dimension(:),   allocatable ::  seedU
    integer(i8),    dimension(:),   allocatable ::  seedG


    Nunif = count(distribution)
    Ngauss = sizeArray - Nunif

    if ( allocated(save_the_setup_u) )  deallocate(save_the_setup_u)
    if ( allocated(save_the_setup_g) )  deallocate(save_the_setup_g)
    if ( allocated(seedU) )             deallocate(seedU)
    if ( allocated(seedG) )             deallocate(seedG)
    if ( Nunif .ne. 0_i4 )  then
       allocate( save_the_setup_u(Nunif, n_save_state) )
       allocate( seedU(Nunif) )
    end if
    if ( Ngauss .ne. 0_i4 ) then
       allocate( save_the_setup_g(Ngauss, n_save_state) )
       allocate( seedG(Ngauss) )
    end if

    ! setup xor4096 - xor4096g
    if ( present(seed) ) then
       ! seed for uniform
       if ( Nunif .ne. 0_i4 )  then
          counter = 1_i4
          do j=1, sizeArray
             if (distribution(j) )   then
                seedU(counter) = seed(j)
                counter = counter + 1_i4
             end if
          end do
          call setupxor( Nunif, save_the_setup_u, seed=seedU )
       end if
       ! seed for gaussian
       if ( Ngauss .ne. 0_i4 ) then
          counter = 1_i4
          do j=1, sizeArray
             if ( .not. distribution(j) )   then
                seedG(counter) = seed(j)
                counter = counter + 1_i4
             end if
          end do
          call setupxor_g( Ngauss, save_the_setup_g, seed=seedG )
       end if
    else
       if ( Nunif .ne. 0_i4 )  call setupxor( Nunif, save_the_setup_u )
       if ( Ngauss .ne. 0_i4 ) call setupxor_g( Ngauss, save_the_setup_g )
    end if

    do j=1,sizeArray
       if ( distribution(j) ) then
          save_the_setup(j,:) = save_the_setup_u(count(distribution(1:j)),:)
       else
          save_the_setup(j,:) = save_the_setup_g(count(.not.distribution(1:j)),:)
       end if
    end do

  end subroutine setupxor4096_dp_1d_dist



  ! ============================================================================
  !   PRIVATE METHODS
  ! ============================================================================


  ! ! Lower tail quantile for standard normal distribution function SINGLE PRECISION D0
  ! function ltqnorm_sp_d0( p )  result(ltqnorm)

  !   use mo_kind,            only    :   sp

  !   implicit none

  !   ! Intent IN
  !   real(sp),       intent(in)  ::  p

  !   ! Results
  !   real(sp)                    ::  ltqnorm

  !   ! Internal variables
  !   real(sp)                    ::  z,q,r

  !   ! Parameters
  !   real(sp),   parameter                   ::  p_low = 0.02425_sp
  !   real(sp),   parameter                   ::  p_high = 1.0_sp - p_low
  !   real(sp),   dimension(6),   parameter   ::  a = (/ -39.6968302866538_sp, 220.946098424521_sp, -275.928510446969_sp, &
  !        138.357751867269_sp, -30.6647980661472_sp, 2.50662827745924_sp /)
  !   real(sp),   dimension(5),   parameter   ::  b = (/ -54.4760987982241_sp, 161.585836858041_sp, -155.698979859887_sp, &
  !        66.8013118877197_sp, -13.2806815528857_sp /)
  !   real(sp),   dimension(6),   parameter   ::  c = (/ -0.00778489400243029_sp, -0.322396458041136_sp, -2.40075827716184_sp, &
  !        -2.54973253934373_sp, 4.37466414146497_sp, 2.93816398269878_sp /)
  !   real(sp),   dimension(4),   parameter   ::  d = (/ 0.00778469570904146_sp, 0.32246712907004_sp, &
  !        2.445134137143_sp, 3.75440866190742_sp /)

  !   if ( p .lt. 0.0_sp .or. p .gt. 1.0_sp )   then
  !      z = 0.0_sp          ! values out of domain
  !   else if ( abs(p) .lt. epsilon(1.0_sp) )  then
  !      z = -huge(1.0_sp)   ! out of range, minus "infinity"
  !   else if ( abs(p-1.0_sp) .lt. epsilon(1.0_sp) )  then
  !      z = huge(1.0_sp)    ! out of range, "infinity"
  !   else if ( p .lt. p_low ) then
  !      q = sqrt(-2.0_sp*log(p))
  !      z = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_sp)
  !   else if ( p .gt. p_high ) then
  !      q = sqrt(-2.0_sp*log(1.0_sp-p))
  !      z = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_sp)
  !   else
  !      q = p - 0.5_sp
  !      r = q*q
  !      z = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/(((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_sp)
  !   end if

  !   ltqnorm = z

  ! end function ltqnorm_sp_d0


  ! Lower tail quantile for standard normal distribution function SINGLE PRECISION D1
  function ltqnorm_sp_d1( p )  result(ltqnorm)

    use mo_kind,            only    :   i4, sp

    implicit none

    ! Intent IN
    real(sp),   dimension(:),   intent(in)  ::  p

    ! Results
    real(sp),   dimension(size(p))          ::  ltqnorm

    ! Internal variables
    integer(i4)                             ::  j       ! counter
    real(sp)                                ::  q,r
    real(sp),   dimension(size(p))          ::  z

    ! Parameters
    real(sp),   parameter                   ::  p_low = 0.02425_sp
    real(sp),   parameter                   ::  p_high = 1.0_sp - p_low
    real(sp),   dimension(6),   parameter   ::  a = (/ -39.6968302866538_sp, 220.946098424521_sp, -275.928510446969_sp, &
         138.357751867269_sp, -30.6647980661472_sp, 2.50662827745924_sp /)
    real(sp),   dimension(5),   parameter   ::  b = (/ -54.4760987982241_sp, 161.585836858041_sp, -155.698979859887_sp, &
         66.8013118877197_sp, -13.2806815528857_sp /)
    real(sp),   dimension(6),   parameter   ::  c = (/ -0.00778489400243029_sp, -0.322396458041136_sp, -2.40075827716184_sp, &
         -2.54973253934373_sp, 4.37466414146497_sp, 2.93816398269878_sp /)
    real(sp),   dimension(4),   parameter   ::  d = (/ 0.00778469570904146_sp, 0.32246712907004_sp, &
         2.445134137143_sp, 3.75440866190742_sp /)

    do j=1, size(p)
       if ( p(j) .lt. 0.0_sp .or. p(j) .gt. 1.0_sp )   then
          z(j) = 0.0_sp           ! values out of domain
       else if ( abs(p(j)) .lt. epsilon(1.0_sp) )  then
          z(j) = -huge(1.0_sp)    ! out of range, minus "infinity"
       else if ( abs(p(j)-1.0_sp) .lt. epsilon(1.0_sp) )  then
          z(j) = huge(1.0_sp)     ! out of range, "infinity"
       else if ( p(j) .lt. p_low ) then
          q = sqrt(-2.0_sp*log(p(j)))
          z(j) = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_sp)
       else if ( p(j) .gt. p_high ) then
          q = sqrt(-2.0_sp*log(1.0_sp-p(j)))
          z(j) = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_sp)
       else
          q = p(j) - 0.5_sp
          r = q*q
          z(j) = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/(((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_sp)
       end if
    end do

    ltqnorm(:) = z(:)

  end function ltqnorm_sp_d1


  ! ! Lower tail quantile for standard normal distribution function DOUBLE PRECISION D0
  ! function ltqnorm_dp_d0( p )  result(ltqnorm)

  !   use mo_kind,            only    :   dp

  !   implicit none

  !   ! Intent IN
  !   real(dp),       intent(in)  ::  p

  !   ! Results
  !   real(dp)                    ::  ltqnorm

  !   ! Internal variables
  !   real(dp)                    ::  z,q,r

  !   ! Parameters
  !   real(dp),   parameter                   ::  p_low = 0.02425_dp
  !   real(dp),   parameter                   ::  p_high = 1.0_dp - p_low
  !   real(dp),   dimension(6),   parameter   ::  a = (/ -39.6968302866538_dp, 220.946098424521_dp, -275.928510446969_dp, &
  !        138.357751867269_dp, -30.6647980661472_dp, 2.50662827745924_dp /)
  !   real(dp),   dimension(5),   parameter   ::  b = (/ -54.4760987982241_dp, 161.585836858041_dp, -155.698979859887_dp, &
  !        66.8013118877197_dp, -13.2806815528857_dp /)
  !   real(dp),   dimension(6),   parameter   ::  c = (/ -0.00778489400243029_dp, -0.322396458041136_dp, -2.40075827716184_dp, &
  !        -2.54973253934373_dp, 4.37466414146497_dp, 2.93816398269878_dp /)
  !   real(dp),   dimension(4),   parameter   ::  d = (/ 0.00778469570904146_dp, 0.32246712907004_dp, &
  !        2.445134137143_dp, 3.75440866190742_dp /)

  !   if ( p .lt. 0.0_dp .or. p .gt. 1.0_dp )   then
  !      z = 0.0_dp          ! values out of domain
  !   else if ( abs(p) .lt. epsilon(1.0_dp) )  then
  !      z = -huge(1.0_dp)   ! out of range, minus "infinity"
  !   else if ( abs(p-1.0_dp) .lt. epsilon(1.0_dp) )  then
  !      z = huge(1.0_dp)    ! out of range, "infinity"
  !   else if ( p .lt. p_low ) then
  !      q = sqrt(-2.0_dp*log(p))
  !      z = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
  !   else if ( p .gt. p_high ) then
  !      q = sqrt(-2.0_dp*log(1.0_dp-p))
  !      z = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
  !   else
  !      q = p - 0.5_dp
  !      r = q*q
  !      z = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/(((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
  !   end if

  !   ltqnorm = z

  ! end function ltqnorm_dp_d0


  ! Lower tail quantile for standard normal distribution function DOUBLE PRECISION D1
  function ltqnorm_dp_d1( p )  result(ltqnorm)

    use mo_kind,            only    :   i4, dp

    implicit none

    ! Intent IN
    real(dp),   dimension(:),   intent(in)  ::  p

    ! Results
    real(dp),   dimension(size(p))          ::  ltqnorm

    ! Internal variables
    integer(i4)                             ::  j       ! counter
    real(dp)                                ::  q,r
    real(dp),   dimension(size(p))          ::  z

    ! Parameters
    real(dp),   parameter                   ::  p_low = 0.02425_dp
    real(dp),   parameter                   ::  p_high = 1.0_dp - p_low
    real(dp),   dimension(6),   parameter   ::  a = (/ -39.6968302866538_dp, 220.946098424521_dp, -275.928510446969_dp, &
         138.357751867269_dp, -30.6647980661472_dp, 2.50662827745924_dp /)
    real(dp),   dimension(5),   parameter   ::  b = (/ -54.4760987982241_dp, 161.585836858041_dp, -155.698979859887_dp, &
         66.8013118877197_dp, -13.2806815528857_dp /)
    real(dp),   dimension(6),   parameter   ::  c = (/ -0.00778489400243029_dp, -0.322396458041136_dp, -2.40075827716184_dp, &
         -2.54973253934373_dp, 4.37466414146497_dp, 2.93816398269878_dp /)
    real(dp),   dimension(4),   parameter   ::  d = (/ 0.00778469570904146_dp, 0.32246712907004_dp, &
         2.445134137143_dp, 3.75440866190742_dp /)

    do j=1, size(p)
       if ( p(j) .lt. 0.0_dp .or. p(j) .gt. 1.0_dp )   then
          z(j) = 0.0_dp           ! values out of domain
       else if ( abs(p(j)) .lt. epsilon(1.0_dp) )  then
          z(j) = -huge(1.0_dp)    ! out of range, minus "infinity"
       else if ( abs(p(j)-1.0_dp) .lt. epsilon(1.0_dp) )  then
          z(j) = huge(1.0_dp)     ! out of range, "infinity"
       else if ( p(j) .lt. p_low ) then
          q = sqrt(-2.0_dp*log(p(j)))
          z(j) = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
       else if ( p(j) .gt. p_high ) then
          q = sqrt(-2.0_dp*log(1.0_dp-p(j)))
          z(j) = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
       else
          q = p(j) - 0.5_dp
          r = q*q
          z(j) = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/(((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
       end if
    end do

    ltqnorm(:) = z(:)

  end function ltqnorm_dp_d1

end module mo_sampling

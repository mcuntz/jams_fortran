module mo_mvn

  ! This module is a wrapper around xor4096g for sampling multivariate normal distribution and is
  ! part of the UFZ CHS Fortran library.

  ! Written  Stephan Thober, Oct 2012

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
  ! along with the UFZ Fortran library. If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2011-2012 Stephan Thober
  !
  implicit none
  !
  public :: mvn
  !
  ! ------------------------------------------------------------------

  !     NAME
  !         mvn

  !     PURPOSE
  !         returns an array of multivariate gaussian variates with
  !         specified covariance matrix (cov)
  !         This function is a wrapper around the random number generator
  !         xor4096g!

  !     CALLING SEQUENCE
  !         out = mvn( lcho [, init] )

  !     INTENT(IN)
  !         real(sp/dp), dimension(:,:) :: lcho ! lower cholesky factor of
  !                                               covariance matrix cov

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp), dimension(:)   :: mvn  ! array of multivariate gaussian
  !                                               variates

  !     INTENT(IN), OPTIONAL
  !         logical                     :: init ! if present, random number
  !                                               generator is initialized with time

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         one can use LAPACK for calculating the lower
  !         Cholesky Factor of the covariance matrix, which is
  !         required by MVN
  !
  !         external DPOTRF   ! Lapack routine for calculating the cholesky factor
  !
  !         integer(i4)              :: info
  !         real(sp), dimension(2,2) :: cov  ! covariance matrix
  !         real(sp), dimension(2,2) :: lcho ! lower cholesky factor
  !
  !         cov = (/ 1., 0.25; 0.25, 1./)
  !         lcho = cov
  !         call DPOTRF('L', size(lcho,1), lcho, size(lcho,1), info) ! changes lower triangular (incl. main diagonal) of Lcho
  !         if ( info .gt. 0_i4 ) then
  !             write(*,*)'GetCholeskyFactor_dp: The algorithm failed!'
  !             stop
  !         end if
  !
  !         ! DONT forget to set upper triangular matrix to zero
  !         lcho(1,2) = 0._sp
  !
  !         rand = mvn( lcho )
  !         CAUTION: one has to sample a time series to stabilize the covariances
  !         between the entries of mvn

  !     LITERATURE
  !         see affine transformation at
  !         http://en.wikipedia.org/wiki/Multivariate_normal_distribution

  !     HISTORY
  !         Written,  Stephan Thober, Oct 2012
  interface mvn
     module procedure mvn_SP, mvn_DP
  end interface mvn
  !
  private
  !
contains

  ! ------------------------------------------------------------------
  !
  function mvn_SP(lcho, init)
    !
    use mo_kind,    only: sp, i4
    use mo_xor4096, only: get_timeseed, xor4096g
    !
    implicit none
    !
    logical,     optional,                intent(in) :: init
    real(sp),    dimension(:,:),          intent(in) :: lcho
    real(sp),    dimension(size(lcho,1))             :: mvn_SP
    integer(i4), dimension(size(lcho,1))             :: seed
    !
    ! consistency check
    if ( size(lcho,1) /= size(lcho,2) ) &
         stop 'ERROR *** size mismatch in variable lcho. fun mvn_SP'
    !
    seed = 0_i4
    !
    if ( present(init) ) then
       !
       ! initialize random number generator with timeseed
       if ( init .and. .true. ) then
          call get_timeseed(seed)
       end if
       !
    end if
    !
    ! generate random
    call xor4096g(seed, mvn_SP)
    !
    mvn_SP = matmul(lcho,mvn_SP)
    !
  end function mvn_SP
  !
  function mvn_DP(lcho, init)
    !
    use mo_kind,    only: dp, i8
    use mo_xor4096, only: get_timeseed, xor4096g
    !
    implicit none
    !
    logical,     optional,       intent(in)    :: init
    real(dp),    dimension(:,:), intent(in)    :: lcho
    real(dp),    dimension(size(lcho,1))       :: mvn_DP
    integer(i8), dimension(size(lcho,1))       :: seed
    !
    ! consistency check
    if ( size(lcho,1) /= size(lcho,2) ) &
         stop 'ERROR *** size mismatch in variable lcho. fun mvn_DP'
    !
    seed = 0_i8
    !
    if ( present(init) ) then
       !
       ! initialize random number generator with timeseed
       call get_timeseed(seed)
       !
    end if
    !
    ! generate random
    call xor4096g(seed, mvn_DP)
    !
    mvn_DP = matmul(lcho,mvn_DP)
    !
  end function mvn_DP
  !
end module mo_mvn

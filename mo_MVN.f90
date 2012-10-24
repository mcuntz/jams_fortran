! ------------------------------------------------------------------------------
!
! module MVN
!     generates multivariate Gaussian Variates given a lower cholesky factor
!
! author: Stephan Thober
!
! created: 8 Aug 2012
!
! ------------------------------------------------------------------------------
module mo_MVN
  !
  implicit none
  !
  private
  !
  public :: MVN
  !
  interface MVN
     module procedure MVN_SP, MVN_DP
  end interface
  !
contains
  ! ------------------------------------------------------------------

  !     NAME
  !         MVN

  !     PURPOSE
  !         returns an array of multivariate gaussian variates with
  !         specified covariance matrix (cov)
  !         This function is a wrapper around the random number generator
  !         xor4096g!

  !     CALLING SEQUENCE
  !         out = MVN( lcho [, init] )
  
  !     INDENT(IN)
  !         real(sp/dp), dimension(:,:) :: lcho ! lower cholesky factor of 
  !                                               covariance matrix cov

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(sp/dp), dimension(:)   :: MVN  ! array of multivariate gaussian 
  !                                               variates

  !     INDENT(IN), OPTIONAL
  !         logical                     :: init ! if present, random number 
  !                                               generator is initialized with time

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         rand = MVN( lcho )
  !         CAUTION: one has to sample a time series to stabilize the covariances 
  !         between the entries

  !     LITERATURE
  !         see affine transformation at
  !         http://en.wikipedia.org/wiki/Multivariate_normal_distribution

  !     HISTORY
  !         Written,  Stephan Thober, Oct 2012
  !
  function MVN_SP(lcho, init)
    !
    use mo_kind,    only: sp, i4
    use mo_xor4096, only: get_timeseed, xor4096g
    !
    implicit none
    !
    logical,     optional,                intent(in) :: init
    real(sp),    dimension(:,:),          intent(in) :: lcho
    real(sp),    dimension(size(lcho,1))             :: MVN_SP
    integer(i4), dimension(size(lcho,1))             :: seed
    !
    ! consistency check
    if ( size(lcho,1) /= size(lcho,2) ) &
         stop 'ERROR *** size mismatch in variable lcho. fun MVN_SP'
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
    call xor4096g(seed, MVN_SP)
    !
    MVN_SP = matmul(lcho,MVN_SP)
    !
  end function MVN_SP
  !
  function MVN_DP(lcho, init)
    !
    use mo_kind,    only: dp, i8
    use mo_xor4096, only: get_timeseed, xor4096g
    !
    implicit none
    !
    logical,     optional,       intent(in)    :: init
    real(dp),    dimension(:,:), intent(in)    :: lcho
    real(dp),    dimension(size(lcho,1))       :: MVN_DP
    integer(i8), dimension(size(lcho,1))       :: seed
    !
    ! consistency check
    if ( size(lcho,1) /= size(lcho,2) ) &
         stop 'ERROR *** size mismatch in variable lcho. fun MVN_DP'
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
    call xor4096g(seed, MVN_DP)
    !
    MVN_DP = matmul(lcho,MVN_DP)
    !
  end function MVN_DP
  !
end module mo_MVN

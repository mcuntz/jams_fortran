!> \file mo_xor4096_apps.f90

!> \brief This module provides useful applications of the random number generator xor4096.

!> \details This module provides useful applications of the random number generator xor4096, e.g. 
!>          (1) generating several random numbers from one stream at once, 
!>          (2) generating multivariate normal distributed random numbers, and 
!>          (3) generation of random numbers within a given range.

!> \authors Juliane Mai, Stephan Thober
!> \date Feb 2013

MODULE mo_xor4096_apps

  ! This module is a template for the JAMS Fortran library.

  ! Written  Juliane Mai, Feb 2013
  ! Modified Stephan Thober, Feb 2013 - add multivariate normal

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2013 Juliane Mai, Stephan Thober
  !
  ! Permission is hereby granted, free of charge, to any person obtaining a copy
  ! of this software and associated documentation files (the "Software"), to deal
  ! in the Software without restriction, including without limitation the rights
  ! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  ! copies of the Software, and to permit persons to whom the Software is
  ! furnished to do so, subject to the following conditions:
  !
  ! The above copyright notice and this permission notice shall be included in all
  ! copies or substantial portions of the Software.
  !
  ! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  ! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  ! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  ! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  ! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  ! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  ! SOFTWARE.

  USE mo_kind,    ONLY: i4, i8, sp, dp
  USE mo_xor4096, ONLY: xor4096, xor4096g, n_save_state

  IMPLICIT NONE

  PUBLIC :: xor4096_array  ! Generation of (subsequent) UNIFORM  random numbers of stream(s)
  PUBLIC :: xor4096g_array ! Generation of (subsequent) GAUSSIAN random numbers of stream(s)
  PUBLIC :: xor4096g_mvn   ! Generation of multivariate NORMAL distributed random numbers
  PUBLIC :: xor4096_range  ! Generation of UNIFORM random numbers within a given range [a,b]

  ! ------------------------------------------------------------------

  !     NAME
  !         xor4096_array

  !     PURPOSE
  !>        \brief Generation of (subsequent) UNIFORM  random numbers of stream(s).
  !
  !>        \details Generation of (subsequent) UNIFORM  random numbers of stream(s).\n
  !>                 It is calling the xor4096 several times and returns a vector of random numbers.\n
  !>                 It can be called with an optional argument to save the state of the current stream 
  !>                 to allow for a later return to exactly this stream.\n
  !>                 The save_state variable has to be (I4) when random numbers are (I4) or (SP) and
  !>                 save_state variable has to be (I8) when random numbers are (I8) or (DP). 
  !>                 The first dimension of save_state can be allocated using n_save_state of the xor4096 module. 
  !>                 If there are N>1 streams used, the second dimension of save_state is N.
  !
  !     CALLING SEQUENCE
  !         call xor4096_array( rn, save_state=save_state )
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4/i8)/real(sp/dp)   :: RN(:)/RN(:,:)"
  !>                                            uniform distributed random numbers with
  !>                                            interval:\n
  !>                                                i4: (-2^31,2^31-1)\n
  !>                                                i8: (-2^63,2^63-1)\n
  !>                                                sp: (0.0_sp, 1.0_sp)\n
  !>                                                dp: (0.0_dp, 1.0_dp)\n
  !>                                            rn(i,j) is the i^{th} number in the j^{th} stream
  !
  !     INTENT(IN), OPTIONAL
  !          None
  !
  !     INTENT(INOUT), OPTIONAL
  !>        \param[in, out] "integer(i4/i8) :: save_state(size(rn,2),n_save_state)"
  !>                                            array carrying state of random number stream\n
  !>                                            this should be used if several streams are used, i.e.
  !>                                            each stream has its own save_state
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !>       \note The random number generator has to be initialized before and therefore is called without a seed.
  !
  !     EXAMPLE
  !         real(dp),    dimension(100)          :: rn
  !         integer(i8)                          :: seed
  !         integer(i8), dimension(n_save_state) :: save_state
  !
  !         ! initialize random number stream
  !         seed = 13546_i8
  !         call xor4096(seed, rn, save_state=save_state)
  !
  !         ! generating random numbers
  !         call xor4096_array(rn, save_state=save_state)
  !         rn --> 100 subsequent (dp) uniform random numbers from stream with seed 13546
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Feb 2013
  !         Modified,  - 
  !
  INTERFACE xor4096_array
     MODULE PROCEDURE &
          xor4096_array_i4_1d, xor4096_array_i4_2d, xor4096_array_i8_1d, xor4096_array_i8_2d, &
          xor4096_array_sp_1d, xor4096_array_sp_2d, xor4096_array_dp_1d, xor4096_array_dp_2d
  END INTERFACE xor4096_array

  ! ------------------------------------------------------------------

  !     NAME
  !         xor4096g_array

  !     PURPOSE
  !>        \brief Generation of (subsequent) GAUSSIAN random numbers of stream(s).
  !
  !>        \details Generation of (subsequent) GAUSSIAN random numbers of stream(s).\n
  !>                 It is calling the xor4096g several times and returns a vector of random numbers.\n
  !>                 It can be called with an optional argument to save the state of the current stream 
  !>                 to allow for a later return to exactly this stream.\n
  !>                 The save_state variable has to be (I4) when random numbers are (SP) and
  !>                 save_state variable has to be (I8) when random numbers are (DP). 
  !>                 The first dimension of save_state can be allocated using n_save_state of the xor4096 module. 
  !>                 If there are N>1 streams used, the second dimension of save_state is N.
  !
  !     CALLING SEQUENCE
  !         call xor4096g_array( rn, save_state=save_state )
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(sp/dp)   :: RN(:)/RN(:,:)"
  !>                                            gaussian distributed random numbers with
  !>                                            interval: mean 0.0 and variance 1.0
  !>                                            rn(i,j) is the i^{th} number in the j^{th} stream
  !
  !     INTENT(IN), OPTIONAL
  !          None
  !
  !     INTENT(INOUT), OPTIONAL
  !>        \param[in, out] "integer(i4/i8) :: save_state(size(rn,2),n_save_state)"
  !>                                            array carrying state of random number stream\n
  !>                                            this should be used if several streams are used, i.e.
  !>                                            each stream has its own save_state
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !>       \note The random number generator has to be initialized before and therefore is called without a seed.
  !
  !     EXAMPLE
  !         real(dp),    dimension(100)          :: rn
  !         integer(i8)                          :: seed
  !         integer(i8), dimension(n_save_state) :: save_state
  !
  !         ! initialize random number stream
  !         seed = 13546_i8
  !         call xor4096g(seed, rn, save_state=save_state)
  !
  !         ! generating random numbers
  !         call xor4096g_array(rn, save_state=save_state)
  !         rn --> 100 subsequent (dp) gaussian random numbers from stream with seed 13546
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Feb 2013
  !         Modified,  - 
  !
  INTERFACE xor4096g_array
     MODULE PROCEDURE &
          xor4096g_array_sp_1d, xor4096g_array_sp_2d, xor4096g_array_dp_1d, xor4096g_array_dp_2d
  END INTERFACE xor4096g_array

  ! ------------------------------------------------------------------

  !     NAME
  !         xor4096g_mvn

  !     PURPOSE
  !>        \brief Generation of multivariate NORMAL distributed random numbers.
  !
  !>        \details Generation of multivariate NORMAL distributed random numbers using a
  !>                 specified covariance matrix (cov).\n
  !>                 This function is a wrapper around the random number generator
  !>                 xor4096g!\n
  !>                 It can be called with an optional argument to save the state of the current stream 
  !>                 to allow for a later return to exactly this stream.\n
  !>                 The save_state variable has to be (I4) when random numbers are (SP) and
  !>                 save_state variable has to be (I8) when random numbers are (DP). 
  !>                 The first dimension of save_state can be allocated using n_save_state of the xor4096 module.
  !>                 The second dimension of save_state is according to the dimension of the covariance matrix.

  !     CALLING SEQUENCE
  !         call xor4096g_mvn( lcho, rn, save_state=save_state )

  !     INTENT(IN)
  !>        \param[in] "real(sp/dp), dimension(:,:) :: lcho"        lower cholesky factor of
  !>                                                                covariance matrix cov

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !>        \param[in, out] "integer(i4/i8) :: save_state( size(lcho,1),n_save_state)"
  !>                                                      array carrying state of random number stream\n
  !>                                                      this should be used if several streams are used, i.e.
  !>                                                      each stream has its own save_state

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(sp/dp), dimension(size(lcho,1))/dimension(:,size(lcho,1)) :: xor4096_mvn" 
  !>                                                      array of multivariate gaussian variates 

  !     RETURN
  !>       \return       

  !     RESTRICTIONS
  !>       \note The random number generator has to be initialized before and therefore is called without a seed.

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
  !         call xor4096g_mvn( lcho , mvn )
  !         CAUTION: one has to sample a time series to stabilize the covariances
  !         between the entries of mvn

  !     LITERATURE
  !         see affine transformation at
  !         http://en.wikipedia.org/wiki/Multivariate_normal_distribution

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Oct 2012
  interface xor4096g_mvn
     module procedure mvn_sp_1d, mvn_dp_1d, mvn_sp_2d, mvn_dp_2d
  end interface xor4096g_mvn
  
  ! ------------------------------------------------------------------

  !     NAME
  !         xor4096_range

  !     PURPOSE
  !>        \brief Generation of UNIFORM random numbers within a given range [a,b].
  !
  !>        \details Generation of UNIFORM random numbers within a given range [a,b].\n
  !>                 It is calling the xor4096 and scales the random numbers to the given range.\n
  !>                 It can be called with an optional argument to save the state of the current stream 
  !>                 to allow for a later return to exactly this stream.\n
  !>                 The save_state variable has to be (I4) when random numbers are (I4) or (SP) and
  !>                 save_state variable has to be (I8) when random numbers are (I8) or (DP). 
  !>                 The first dimension of save_state can be allocated using n_save_state of the xor4096 module. 
  !>                 If there are N>1 streams used, the second dimension of save_state is N.
  !
  !     CALLING SEQUENCE
  !         call xor4096_range( range, rn, save_state=save_state )
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8)/real(sp/dp) :: rn_range(2)" lower and upper bound of the random number 
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4/i8)/real(sp/dp)   :: RN/RN(:)/RN(:,:)"
  !>                                            uniform distributed random numbers with
  !>                                            interval:\n
  !>                                                i4: [a,b]\n
  !>                                                i8: [a,b]\n
  !>                                                sp: [a,b]\n
  !>                                                dp: [a,b]\n
  !>                                            rn(i,j) is the i^{th} number in the j^{th} stream
  !
  !     INTENT(IN), OPTIONAL
  !          None
  !
  !     INTENT(INOUT), OPTIONAL
  !>        \param[in, out] "integer(i4/i8) :: save_state(size(rn,2),n_save_state)"
  !>                                            array carrying state of random number stream\n
  !>                                            this should be used if several streams are used, i.e.
  !>                                            each stream has its own save_state
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !>       \note The random number generator has to be initialized before and therefore is called without a seed.
  !
  !     EXAMPLE
  !         real(dp),    dimension(100)          :: rn
  !         integer(i8)                          :: seed
  !         integer(i8), dimension(n_save_state) :: save_state
  !
  !         ! initialize random number stream
  !         seed = 13546_i8
  !         call xor4096(seed, rn, save_state=save_state)
  !
  !         ! generating random numbers
  !         call xor4096_range((/2.0_dp,4.0_dp/), rn, save_state=save_state)
  !         rn --> 100 subsequent (dp) uniform random numbers in interval [2.0,4.0] from stream with seed 13546
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Feb 2013
  !         Modified,  - 
  !
  INTERFACE xor4096_range
     MODULE PROCEDURE &
          xor4096_range_i4_0d, xor4096_range_i4_1d, xor4096_range_i4_2d, &
          xor4096_range_i8_0d, xor4096_range_i8_1d, xor4096_range_i8_2d, &
          xor4096_range_sp_0d, xor4096_range_sp_1d, xor4096_range_sp_2d, &
          xor4096_range_dp_0d, xor4096_range_dp_1d, xor4096_range_dp_2d
  END INTERFACE xor4096_range

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! *********************************************************************
  ! TEST XOR4096_ARRAY
  ! *********************************************************************

  SUBROUTINE xor4096_array_i4_1d(rn, save_state)

    IMPLICIT NONE

    INTEGER(I4), DIMENSION(:),                      INTENT(OUT)   :: rn
    INTEGER(I4), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)  :: i

    if (present(save_state)) then
       do i=1, size(rn,1)
          call xor4096(0_i4, rn(i), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(0_i4, rn(i))
       end do
    endif

  END SUBROUTINE xor4096_array_i4_1d

  SUBROUTINE xor4096_array_i4_2d(rn, save_state)

    IMPLICIT NONE

    INTEGER(I4), DIMENSION(:,:),           INTENT(OUT)   :: rn             ! dim2=number of streams
    INTEGER(I4), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)                        :: i
    INTEGER(I4), DIMENSION(size(rn,2)) :: seed_tmp

    seed_tmp = 0_i4

    if (present(save_state)) then
       if (size(save_state,1) .ne. size(rn,2) .or. size(save_state,2) .ne. n_save_state) then
          stop 'xor4096_array_i4_2d: dimensions of save_state matrix are not correct'
       end if
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:))
       end do
    endif

  END SUBROUTINE xor4096_array_i4_2d

  SUBROUTINE xor4096_array_i8_1d(rn, save_state)

    IMPLICIT NONE

    INTEGER(I8), DIMENSION(:),                      INTENT(OUT)   :: rn
    INTEGER(I8), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)  :: i

    if (present(save_state)) then
       do i=1, size(rn,1)
          call xor4096(0_i8, rn(i), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(0_i8, rn(i))
       end do
    endif

  END SUBROUTINE xor4096_array_i8_1d

  SUBROUTINE xor4096_array_i8_2d(rn, save_state)

    IMPLICIT NONE

    INTEGER(I8), DIMENSION(:,:),           INTENT(OUT)   :: rn             ! dim2=number of streams
    INTEGER(I8), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)                        :: i
    INTEGER(I8), DIMENSION(size(rn,2)) :: seed_tmp

    seed_tmp = 0_i8

    if (present(save_state)) then
       if (size(save_state,1) .ne. size(rn,2) .or. size(save_state,2) .ne. n_save_state) then
          stop 'xor4096_array_i8_2d: dimensions of save_state matrix are not correct'
       end if
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:))
       end do
    endif

  END SUBROUTINE xor4096_array_i8_2d

  SUBROUTINE xor4096_array_sp_1d(rn, save_state)

    IMPLICIT NONE

    REAL(SP),    DIMENSION(:),                      INTENT(OUT)   :: rn
    INTEGER(I4), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)  :: i

    if (present(save_state)) then
       do i=1, size(rn,1)
          call xor4096(0_i4, rn(i), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(0_i4, rn(i))
       end do
    endif

  END SUBROUTINE xor4096_array_sp_1d

  SUBROUTINE xor4096_array_sp_2d(rn, save_state)

    IMPLICIT NONE

    REAL(SP),    DIMENSION(:,:),           INTENT(OUT)   :: rn             ! dim2=number of streams
    INTEGER(I4), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)                        :: i
    INTEGER(I4), DIMENSION(size(rn,2)) :: seed_tmp

    seed_tmp = 0_i4

    if (present(save_state)) then
       if (size(save_state,1) .ne. size(rn,2) .or. size(save_state,2) .ne. n_save_state) then
          stop 'xor4096_array_sp_2d: dimensions of save_state matrix are not correct'
       end if
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:))
       end do
    endif

  END SUBROUTINE xor4096_array_sp_2d

  SUBROUTINE xor4096_array_dp_1d(rn, save_state)

    IMPLICIT NONE

    REAL(DP),    DIMENSION(:),                      INTENT(OUT)   :: rn
    INTEGER(I8), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)  :: i

    if (present(save_state)) then
       do i=1, size(rn,1)
          call xor4096(0_i8, rn(i), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(0_i8, rn(i))
       end do
    endif

  END SUBROUTINE xor4096_array_dp_1d

  SUBROUTINE xor4096_array_dp_2d(rn, save_state)

    IMPLICIT NONE

    REAL(DP),    DIMENSION(:,:),           INTENT(OUT)   :: rn             ! dim2=number of streams
    INTEGER(I8), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)                        :: i
    INTEGER(I8), DIMENSION(size(rn,2)) :: seed_tmp

    seed_tmp = 0_i8

    if (present(save_state)) then
       if (size(save_state,1) .ne. size(rn,2) .or. size(save_state,2) .ne. n_save_state) then
          stop 'xor4096_array_dp_2d: dimensions of save_state matrix are not correct'
       end if
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:))
       end do
    endif

  END SUBROUTINE xor4096_array_dp_2d

  ! *********************************************************************
  ! TEST XOR4096G_ARRAY
  ! *********************************************************************

  SUBROUTINE xor4096g_array_sp_1d(rn, save_state)

    IMPLICIT NONE

    REAL(SP),    DIMENSION(:),                      INTENT(OUT)   :: rn
    INTEGER(I4), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)  :: i

    if (present(save_state)) then
       do i=1, size(rn,1)
          call xor4096g(0_i4, rn(i), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096g(0_i4, rn(i))
       end do
    endif

  END SUBROUTINE xor4096g_array_sp_1d

  SUBROUTINE xor4096g_array_sp_2d(rn, save_state)

    IMPLICIT NONE

    REAL(SP),    DIMENSION(:,:),           INTENT(OUT)   :: rn             ! dim2=number of streams
    INTEGER(I4), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)                        :: i
    INTEGER(I4), DIMENSION(size(rn,2)) :: seed_tmp

    seed_tmp = 0_i4

    if (present(save_state)) then
       if (size(save_state,1) .ne. size(rn,2) .or. size(save_state,2) .ne. n_save_state) then
          stop 'xor4096g_array_sp_2d: dimensions of save_state matrix are not correct'
       end if
       do i=1, size(rn,1)
          call xor4096g(seed_tmp, rn(i,:), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096g(seed_tmp, rn(i,:))
       end do
    endif

  END SUBROUTINE xor4096g_array_sp_2d

  SUBROUTINE xor4096g_array_dp_1d(rn, save_state)

    IMPLICIT NONE

    REAL(DP),    DIMENSION(:),                      INTENT(OUT)   :: rn
    INTEGER(I8), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)  :: i

    if (present(save_state)) then
       do i=1, size(rn,1)
          call xor4096g(0_i8, rn(i), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096g(0_i8, rn(i))
       end do
    endif

  END SUBROUTINE xor4096g_array_dp_1d

  SUBROUTINE xor4096g_array_dp_2d(rn, save_state)

    IMPLICIT NONE

    REAL(DP),    DIMENSION(:,:),           INTENT(OUT)   :: rn             ! dim2=number of streams
    INTEGER(I8), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)                        :: i
    INTEGER(I8), DIMENSION(size(rn,2)) :: seed_tmp

    seed_tmp = 0_i8

    if (present(save_state)) then
       if (size(save_state,1) .ne. size(rn,2) .or. size(save_state,2) .ne. n_save_state) then
          stop 'xor4096g_array_dp_2d: dimensions of save_state matrix are not correct'
       end if
       do i=1, size(rn,1)
          call xor4096g(seed_tmp, rn(i,:), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096g(seed_tmp, rn(i,:))
       end do
    endif

  END SUBROUTINE xor4096g_array_dp_2d

  ! *********************************************************************
  ! XOR4096G_MVN
  ! *********************************************************************

  subroutine mvn_SP_1d(lcho, mvn, save_state)
    
    implicit none
                 
    real(sp),    dimension(:,:),           intent(in)    :: lcho
    real(sp),    dimension(size(lcho,1)),  intent(out)   :: mvn
    integer(i4), dimension(:,:), optional, intent(inout) :: save_state

    ! local variables
    integer(i4), dimension(size(lcho,1))  :: seed_tmp

    seed_tmp = 0_i4
    
    ! consistency check
    if ( size(lcho,1) /= size(lcho,2) ) then
       stop 'mvn_sp_1d: size mismatch in variable lcho.'
    end if
    
    ! ! initialize if non-zero
    ! if ( all(seed .ne. 0_i4) ) then
    !    if ( present(save_state) ) then
    !       call xor4096g(seed, mvn, save_state=save_state)
    !    else
    !       call xor4096g(seed, mvn)
    !    end if
    !    seed = 0_i4
    ! end if
    
    ! generate random
    if ( present(save_state) ) then
       if (size(save_state,1) .ne. size(lcho,1) .or. size(save_state,2) .ne. n_save_state) then
          stop 'mvn_sp_1d: dimensions of save_state matrix are not correct'
       end if
       call xor4096g(seed_tmp, mvn, save_state=save_state)
    else
       call xor4096g(seed_tmp, mvn)
    end if
    
    mvn = matmul(lcho,mvn)
    
  end subroutine mvn_SP_1d
  
  subroutine mvn_DP_1d(lcho, mvn, save_state)
    
    implicit none
               
    real(dp),    dimension(:,:),           intent(in)    :: lcho
    real(dp),    dimension(size(lcho,1)),  intent(out)   :: mvn
    integer(i8), dimension(:,:), optional, intent(inout) :: save_state

    ! local variables
    integer(i8), dimension(size(lcho,1))  :: seed_tmp

    seed_tmp = 0_i8

    ! consistency check
    if ( size(lcho,1) /= size(lcho,2) ) then
       stop 'mvn_dp_1d: size mismatch in variable lcho.'
    end if
    
    ! ! initialize if non-zero
    ! if ( all(seed .ne. 0_i8) ) then
    !    if ( present(save_state) ) then
    !       call xor4096g(seed, mvn, save_state=save_state)
    !    else
    !       call xor4096g(seed, mvn)
    !    end if
    !    seed = 0_i8
    ! end if
    
    ! generate random
    if ( present(save_state) ) then
       if (size(save_state,1) .ne. size(lcho,1) .or. size(save_state,2) .ne. n_save_state) then
          stop 'mvn_dp_1d: dimensions of save_state matrix are not correct'
       end if
       call xor4096g(seed_tmp, mvn, save_state=save_state)
    else
       call xor4096g(seed_tmp, mvn)
    end if
    
    mvn = matmul(lcho,mvn)
    
  end subroutine mvn_DP_1d

  subroutine mvn_SP_2d(lcho, mvn, save_state)
    
    implicit none
              
    real(sp),    dimension(:,:),           intent(in)    :: lcho
    real(sp),    dimension(:,:),           intent(out)   :: mvn         ! dim1=number of sets, 
    !                                                                                         ! dim2=size(lcho,1)
    integer(i4), dimension(:,:), optional, intent(inout) :: save_state

    ! local variables
    integer(i4)                           :: i
    integer(i4), dimension(size(lcho,1))  :: seed_tmp

    seed_tmp = 0_i4
    
    ! consistency check
    if ( size(lcho,1) /= size(lcho,2) ) then
       stop 'mvn_sp_2d: size mismatch in variable lcho.'
    end if
    
    ! ! initialize if non-zero
    ! if ( all(seed .ne. 0_i4) ) then
    !    if ( present(save_state) ) then
    !       call xor4096g(seed, mvn(1,:), save_state=save_state)
    !    else
    !       call xor4096g(seed, mvn(1,:))
    !    end if
    !    seed = 0_i4
    ! end if
    
    ! generate random
    if ( present(save_state) ) then
       if (size(save_state,1) .ne. size(lcho,1) .or. size(save_state,2) .ne. n_save_state) then
          stop 'mvn_sp_2d: dimensions of save_state matrix are not correct'
       end if
       do i=1, size(mvn,1)
          call xor4096g(seed_tmp, mvn(i,:), save_state=save_state)
          mvn(i,:) = matmul(lcho,mvn(i,:))
       end do
    else
       do i=1, size(mvn,1)
          call xor4096g(seed_tmp, mvn(i,:))
          mvn(i,:) = matmul(lcho,mvn(i,:))
       end do
    end if
    
  end subroutine mvn_SP_2d
  
  subroutine mvn_DP_2d(lcho, mvn, save_state)
    !
    implicit none
    !            
    real(dp),    dimension(:,:),           intent(in)    :: lcho
    real(dp),    dimension(:,:),           intent(out)   :: mvn         ! dim1=number of sets, 
    !                                                                                         ! dim2=size(lcho,1)
    integer(i8), dimension(:,:), optional, intent(inout) :: save_state

    ! local variables
    integer(i4)                           :: i
    integer(i8), dimension(size(lcho,1))  :: seed_tmp

    seed_tmp = 0_i8
    
    ! consistency check
    if ( size(lcho,1) /= size(lcho,2) ) then
       stop 'mvn_dp_2d: size mismatch in variable lcho.'
    end if
    
    ! ! initialize if non-zero
    ! if ( all(seed .ne. 0_i8) ) then
    !    if ( present(save_state) ) then
    !       call xor4096g(seed, mvn(1,:), save_state=save_state)
    !    else
    !       call xor4096g(seed, mvn(1,:))
    !    end if
    !    seed = 0_i8
    ! end if
    
    ! generate random
    if ( present(save_state) ) then
       if (size(save_state,1) .ne. size(lcho,1) .or. size(save_state,2) .ne. n_save_state) then
          stop 'mvn_dp_2d: dimensions of save_state matrix are not correct'
       end if
       do i=1, size(mvn,1)
          call xor4096g(seed_tmp, mvn(i,:), save_state=save_state)
          mvn(i,:) = matmul(lcho,mvn(i,:))
       end do
    else
       do i=1, size(mvn,1)
          call xor4096g(seed_tmp, mvn(i,:))
          mvn(i,:) = matmul(lcho,mvn(i,:))
       end do
    end if
    
  end subroutine mvn_DP_2d
  
  ! *********************************************************************
  ! XOR4096_RANGE
  ! *********************************************************************

  SUBROUTINE xor4096_range_i4_0d(rn_range, rn, save_state)

    IMPLICIT NONE

    INTEGER(I4), DIMENSION(2),                      INTENT(IN)    :: rn_range       ! lower and upper bound
    INTEGER(I4),                                    INTENT(OUT)   :: rn             
    INTEGER(I4), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    REAL(SP) :: rn_tmp

    if ( rn_range(1) .gt. rn_range(2) ) then
       stop 'xor4096_range_i4_0d: lower bound larger than upper bound'
    end if

    if (present(save_state)) then
       call xor4096(0_i4, rn, save_state=save_state)
    else
       call xor4096(0_i4, rn)
    endif

    ! scale [-2^31, 2^31) to [0.0, 1.0)
    rn_tmp = real(rn,sp) / 2.0_sp**32_i4 + 0.5_sp

    ! scale [0.0, 1.0) to [rn_range(1),rn_range(2)]
    rn = int( rn_tmp * ( real(rn_range(2) - rn_range(1) + 1_i4,sp) ) + real(rn_range(1),sp) , i4)

  END SUBROUTINE xor4096_range_i4_0d

  SUBROUTINE xor4096_range_i4_1d(rn_range, rn, save_state)

    IMPLICIT NONE
    
    INTEGER(I4), DIMENSION(2),                      INTENT(IN)    :: rn_range       ! lower and upper bound
    INTEGER(I4), DIMENSION(:),                      INTENT(OUT)   :: rn         
    INTEGER(I4), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)                       :: i
    REAL(SP),   DIMENSION(size(rn,1)) :: rn_tmp

    if ( rn_range(1) .gt. rn_range(2) ) then
       stop 'xor4096_range_i4_1d: lower bound larger than upper bound'
    end if

    if (present(save_state)) then
       do i=1, size(rn,1)
          call xor4096(0_i4, rn(i), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(0_i4, rn(i))
       end do
    endif

    ! scale [-2^31, 2^31) to [0.0, 1.0)
    rn_tmp = real(rn,sp) / 2.0_sp**32_i4 + 0.5_sp

    ! scale [0.0, 1.0) to [rn_range(1),rn_range(2)]
    rn = int( rn_tmp * ( real(rn_range(2) - rn_range(1) + 1_i4,sp) ) + real(rn_range(1),sp) , i4)

  END SUBROUTINE xor4096_range_i4_1d

  SUBROUTINE xor4096_range_i4_2d(rn_range, rn, save_state)

    IMPLICIT NONE

    INTEGER(I4), DIMENSION(2),             INTENT(IN)    :: rn_range       ! lower and upper bound
    INTEGER(I4), DIMENSION(:,:),           INTENT(OUT)   :: rn             ! dim2=number of streams
    INTEGER(I4), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)                                    :: i
    REAL(SP),    DIMENSION(size(rn,1), size(rn,2)) :: rn_tmp
    INTEGER(I4), DIMENSION(size(rn,2))             :: seed_tmp

    seed_tmp = 0_i4

    if ( rn_range(1) .gt. rn_range(2) ) then
       stop 'xor4096_range_i4_2d: lower bound larger than upper bound'
    end if

    if (present(save_state)) then
       if (size(save_state,1) .ne. size(rn,2) .or. size(save_state,2) .ne. n_save_state) then
          stop 'xor4096_range_i4_2d: dimensions of save_state matrix are not correct'
       end if
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:))
       end do
    endif

    ! scale [-2^31, 2^31) to [0.0, 1.0)
    rn_tmp = real(rn,sp) / 2.0_sp**32_i4 + 0.5_sp

    ! scale [0.0, 1.0) to [rn_range(1),rn_range(2)]
    rn = int( rn_tmp * ( real(rn_range(2) - rn_range(1) + 1_i4,sp) ) + real(rn_range(1),sp) , i4)

  END SUBROUTINE xor4096_range_i4_2d

  SUBROUTINE xor4096_range_i8_0d(rn_range, rn, save_state)

    IMPLICIT NONE

    INTEGER(I8), DIMENSION(2),                      INTENT(IN)    :: rn_range       ! lower and upper bound
    INTEGER(I8),                                    INTENT(OUT)   :: rn             
    INTEGER(I8), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    REAL(DP) :: rn_tmp

    if ( rn_range(1) .gt. rn_range(2) ) then
       stop 'xor4096_range_i8_0d: lower bound larger than upper bound'
    end if

    if (present(save_state)) then
       call xor4096(0_i8, rn, save_state=save_state)
    else
       call xor4096(0_i8, rn)
    endif

    ! scale [-2^63, 2^63) to [0.0, 1.0)
    rn_tmp = real(rn,dp) / 2.0_dp**64_i4 + 0.5_dp

    ! scale [0.0, 1.0) to [rn_range(1),rn_range(2)]
    rn = int( rn_tmp * ( real(rn_range(2) - rn_range(1) + 1_i8,dp) ) + real(rn_range(1),dp) , i8)

  END SUBROUTINE xor4096_range_i8_0d

  SUBROUTINE xor4096_range_i8_1d(rn_range, rn, save_state)

    IMPLICIT NONE

    INTEGER(I8), DIMENSION(2),                      INTENT(IN)    :: rn_range       ! lower and upper bound
    INTEGER(I8), DIMENSION(:),                      INTENT(OUT)   :: rn         
    INTEGER(I8), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)                       :: i
    REAL(DP),   DIMENSION(size(rn,1)) :: rn_tmp

    if ( rn_range(1) .gt. rn_range(2) ) then
       stop 'xor4096_range_i8_1d: lower bound larger than upper bound'
    end if

    if (present(save_state)) then
       do i=1, size(rn,1)
          call xor4096(0_i8, rn(i), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(0_i8, rn(i))
       end do
    endif

    ! scale [-2^63, 2^63) to [0.0, 1.0)
    rn_tmp = real(rn,dp) / 2.0_dp**64_i4 + 0.5_dp

    ! scale [0.0, 1.0) to [rn_range(1),rn_range(2)]
    rn = int( rn_tmp * ( real(rn_range(2) - rn_range(1) + 1_i8,dp) ) + real(rn_range(1),dp) , i8)

  END SUBROUTINE xor4096_range_i8_1d

  SUBROUTINE xor4096_range_i8_2d(rn_range, rn, save_state)

    IMPLICIT NONE

    INTEGER(I8), DIMENSION(2),             INTENT(IN)    :: rn_range       ! lower and upper bound
    INTEGER(I8), DIMENSION(:,:),           INTENT(OUT)   :: rn             ! dim2=number of streams
    INTEGER(I8), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)                                   :: i
    REAL(DP),   DIMENSION(size(rn,1), size(rn,2)) :: rn_tmp
    INTEGER(I8), DIMENSION(size(rn,2))  :: seed_tmp

    seed_tmp = 0_i8

    if ( rn_range(1) .gt. rn_range(2) ) then
       stop 'xor4096_range_i8_2d: lower bound larger than upper bound'
    end if

    if (present(save_state)) then
       if (size(save_state,1) .ne. size(rn,2) .or. size(save_state,2) .ne. n_save_state) then
          stop 'xor4096_range_i8_2d: dimensions of save_state matrix are not correct'
       end if
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:))
       end do
    endif

    ! scale [-2^63, 2^63) to [0.0, 1.0)
    rn_tmp = real(rn,dp) / 2.0_dp**64_i4 + 0.5_dp

    ! scale [0.0, 1.0) to [rn_range(1),rn_range(2)]
    rn = int( rn_tmp * ( real(rn_range(2) - rn_range(1) + 1_i8,dp) ) + real(rn_range(1),dp) , i8)

  END SUBROUTINE xor4096_range_i8_2d

  SUBROUTINE xor4096_range_sp_0d(rn_range, rn, save_state)

    IMPLICIT NONE

    REAL(SP),    DIMENSION(2),                      INTENT(IN)    :: rn_range       ! lower and upper bound
    REAL(SP),                                       INTENT(OUT)   :: rn             
    INTEGER(I4), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    if ( rn_range(1) .gt. rn_range(2) ) then
       stop 'xor4096_range_sp_0d: lower bound larger than upper bound'
    end if

    if (present(save_state)) then
       call xor4096(0_i4, rn, save_state=save_state)
    else
       call xor4096(0_i4, rn)
    endif
    rn = rn * ( rn_range(2) - rn_range(1) ) + rn_range(1)

  END SUBROUTINE xor4096_range_sp_0d

  SUBROUTINE xor4096_range_sp_1d(rn_range, rn, save_state)

    IMPLICIT NONE

    REAL(SP),    DIMENSION(2),                      INTENT(IN)    :: rn_range       ! lower and upper bound
    REAL(SP),    DIMENSION(:),                      INTENT(OUT)   :: rn         
    INTEGER(I4), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)  :: i

    if ( rn_range(1) .gt. rn_range(2) ) then
       stop 'xor4096_range_sp_1d: lower bound larger than upper bound'
    end if

    if (present(save_state)) then
       do i=1, size(rn,1)
          call xor4096(0_i4, rn(i), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(0_i4, rn(i))
       end do
    endif
    rn = rn * ( rn_range(2) - rn_range(1) ) + rn_range(1)

  END SUBROUTINE xor4096_range_sp_1d

  SUBROUTINE xor4096_range_sp_2d(rn_range, rn, save_state)

    IMPLICIT NONE

    REAL(SP),    DIMENSION(2),             INTENT(IN)    :: rn_range       ! lower and upper bound
    REAL(SP),    DIMENSION(:,:),           INTENT(OUT)   :: rn             ! dim2=number of streams
    INTEGER(I4), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)  :: i
    INTEGER(I4), DIMENSION(size(rn,2))  :: seed_tmp

    seed_tmp = 0_i4

    if ( rn_range(1) .gt. rn_range(2) ) then
       stop 'xor4096_range_sp_2d: lower bound larger than upper bound'
    end if

    if (present(save_state)) then
       if (size(save_state,1) .ne. size(rn,2) .or. size(save_state,2) .ne. n_save_state) then
          stop 'xor4096_range_sp_2d: dimensions of save_state matrix are not correct'
       end if
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:))
       end do
    endif
    rn = rn * ( rn_range(2) - rn_range(1) ) + rn_range(1)

  END SUBROUTINE xor4096_range_sp_2d

  SUBROUTINE xor4096_range_dp_0d(rn_range, rn, save_state)

    IMPLICIT NONE

    REAL(DP),    DIMENSION(2),                      INTENT(IN)    :: rn_range       ! lower and upper bound
    REAL(DP),                                       INTENT(OUT)   :: rn             
    INTEGER(I8), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    if ( rn_range(1) .gt. rn_range(2) ) then
       stop 'xor4096_range_dp_0d: lower bound larger than upper bound'
    end if

    if (present(save_state)) then
       call xor4096(0_i8, rn, save_state=save_state)
    else
       call xor4096(0_i8, rn)
    endif
    rn = rn * ( rn_range(2) - rn_range(1) ) + rn_range(1)

  END SUBROUTINE xor4096_range_dp_0d

  SUBROUTINE xor4096_range_dp_1d(rn_range, rn, save_state)

    IMPLICIT NONE

    REAL(DP),    DIMENSION(2),                      INTENT(IN)    :: rn_range       ! lower and upper bound
    REAL(DP),    DIMENSION(:),                      INTENT(OUT)   :: rn         
    INTEGER(I8), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)  :: i

    if ( rn_range(1) .gt. rn_range(2) ) then
       stop 'xor4096_range_dp_1d: lower bound larger than upper bound'
    end if

    if (present(save_state)) then
       do i=1, size(rn,1)
          call xor4096(0_i8, rn(i), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(0_i8, rn(i))
       end do
    endif
    rn = rn * ( rn_range(2) - rn_range(1) ) + rn_range(1)

  END SUBROUTINE xor4096_range_dp_1d

  SUBROUTINE xor4096_range_dp_2d(rn_range, rn, save_state)

    IMPLICIT NONE

    REAL(DP),    DIMENSION(2),             INTENT(IN)    :: rn_range       ! lower and upper bound
    REAL(DP),    DIMENSION(:,:),           INTENT(OUT)   :: rn             ! dim2=number of streams
    INTEGER(I8), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: save_state

    ! local variables
    INTEGER(I4)  :: i
    INTEGER(I8), DIMENSION(size(rn,2))  :: seed_tmp

    seed_tmp = 0_i8

    if ( rn_range(1) .gt. rn_range(2) ) then
       stop 'xor4096_range_dp_2d: lower bound larger than upper bound'
    end if

    if (present(save_state)) then
       if (size(save_state,1) .ne. size(rn,2) .or. size(save_state,2) .ne. n_save_state) then
          stop 'xor4096_range_dp_2d: dimensions of save_state matrix are not correct' 
       end if
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:), save_state=save_state)
       end do
    else
       do i=1, size(rn,1)
          call xor4096(seed_tmp, rn(i,:))
       end do
    endif
    rn = rn * ( rn_range(2) - rn_range(1) ) + rn_range(1)

  END SUBROUTINE xor4096_range_dp_2d


END MODULE mo_xor4096_apps

MODULE mo_specan

  ! This module provides routines for spectral analysis
  ! using routines from fftpack http://orion.math.iastate.edu/burkardt/f_src/fftpack/fftpack.html

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2013 Philipp Bosecker
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

  USE mo_kind,      ONLY: i4, sp, spc, dp, dpc
  USE mo_constants, ONLY: PI_sp, PI_dp
  USE mo_nelmin,    ONLY: nelminrange

  IMPLICIT NONE

  PUBLIC :: cost            ! Cost function for fitting linear reservoir model
  PUBLIC :: linres          ! Fitting linear reservoir model
  PUBLIC :: periodogram     ! Periodogram of given dataset using FFT
  PUBLIC :: setmeas         ! Set measurement as module variable

  ! ------------------------------------------------------------------

  !     NAME
  !         cost

  !     PURPOSE
  !         Define objective function to fit linear reservoir model linres to data.

  !     CALLING SEQUENCE
  !         func = cost(para)

  !     INTENT(IN)
  !         real(sp/dp) :: para(:)        parameter of function

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         see test example test/test_mo_specan

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written, December 2013, Philipp Bosecker

  ! ------------------------------------------------------------------

  INTERFACE cost
     MODULE PROCEDURE cost_sp, cost_dp
  END INTERFACE cost

  ! ------------------------------------------------------------------

  !     NAME
  !         linres

  !     PURPOSE
  !         Fitting a linear reservoir model f dependent on the frequency w
  !                  f(w)=1/(a(1+tc^2w^2)) 
  !         introduced by Gelhar et. al (1974) to the double logarithmic periodogram of a time series.
  !         The time series is for example groundwater fluctuations.
  !         The periodogram can be determined using the subroutine periodogram.
  !         The best fitting parameters a and tc are returned.
  !         In the linear reservoir model a is a discharge constant and
  !         tc is the characteristic respond time.

  !     CALLING SEQUENCE
  !         call linres(data,bestpara,a_min,a_max,tc_min,tc_max,initial_a,initial_tc)

  !     INTENT(IN)
  !         real(sp/dp) :: data(:)        Time series
  !         real(sp/dp) :: a_min          minimum of first parameter
  !         real(sp/dp) :: a_max          maximum of first parameter
  !         real(sp/dp) :: tc_min         minimum of second parameter
  !         real(sp/dp) :: tc_max         maximum of second parameter
  !         real(sp/dp) :: initial_a      first initial parameter
  !         real(sp/dp) :: initial_tc     second initial parameter


  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp(dp) :: bestpara(2)    best fitting parameter a and tc found for linear reservoir model

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         see test example test/test_mo_specan

  !     LITERATURE
  !       Lynn W. Gelhar and John L. Wilson,
  !           Ground-Water Quality Modeling,
  !           Ground Water, 1974

  !     HISTORY
  !         Written, December 2013, Philipp Bosecker

  ! ------------------------------------------------------------------

  INTERFACE linres
     MODULE PROCEDURE linres_sp, linres_dp
  END INTERFACE linres

  ! ------------------------------------------------------------------

  !     NAME
  !         periodogram

  !     PURPOSE
  !         Computes the periodogram of given dataset using FFT. The periodogram
  !         is calculated as the square of the absolut value of the discrete complex
  !         Fourier transformation. The output periodogram has the size of N/2.

  !     CALLING SEQUENCE
  !         call periodogram(data,output)

  !     INTENT(IN)
  !         real(sp/dp) :: data(:)        Time series

  !     INTENT(INOUT)
  !         NONE

  !     INTENT(OUT)
  !         real(sp/dp) :: output(:,2)    frequency, power

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         see test example test/test_mo_specan

  !     LITERATURE
  !       P.N. Swarztrauber,
  !           Vectorizing the FFTs, in Parallel Computations (G. Rodrigue, ed.),
  !           Academic Press, 1982, pp. 51--83

  !     HISTORY
  !         Written, December 2013, Philipp Bosecker

  ! ------------------------------------------------------------------

  INTERFACE periodogram
     MODULE PROCEDURE periodogram_sp, periodogram_dp
  END INTERFACE periodogram

  ! ------------------------------------------------------------------

  !     NAME
  !         setmeas

  !     PURPOSE
  !         Read in data/measurements and store them as global module variable meas.
  !         Global module variable is needed since optimization only allows for a single argument 
  !         (which is parameter set) of the objective function.

  !     CALLING SEQUENCE
  !         call setmeas(data)

  !     INTENT(IN)
  !         real(sp/dp) :: data(:)        data/measurements

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         see test example test/test_mo_specan

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written, December 2013, Philipp Bosecker

  ! ------------------------------------------------------------------

  INTERFACE setmeas
     MODULE PROCEDURE setmeas_dp, setmeas_sp
  END INTERFACE setmeas

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

  REAL(dp), DIMENSION(:), ALLOCATABLE :: meas_dp      ! data/measurement
  REAL(sp), DIMENSION(:), ALLOCATABLE :: meas_sp      ! data/measurement

  !**************************************************
  !
  ! Changes that have been done in FFTPACK routines:
  !
  ! - Declaration of data type
  ! - single and double precision versions
  ! - fixed that wrong data types were overgiven
  !
  !**************************************************

  ! Private routines, mostly from fftpack  
  INTERFACE cfftf
     MODULE PROCEDURE cfftf_sp, cfftf_dp
  END INTERFACE cfftf
  INTERFACE cfftf1
     MODULE PROCEDURE cfftf1_sp, cfftf1_dp
  END INTERFACE cfftf1
  INTERFACE cffti
     MODULE PROCEDURE cffti_sp, cffti_dp
  END INTERFACE cffti
  INTERFACE cffti1
     MODULE PROCEDURE cffti1_sp, cffti1_dp
  END INTERFACE cffti1
  INTERFACE i_factor
     MODULE PROCEDURE i_factor
  END INTERFACE i_factor
  !  INTERFACE passf
  !     MODULE PROCEDURE passf_sp, passf_dp
  !  END INTERFACE passf
  !  INTERFACE passf2
  !     MODULE PROCEDURE passf2_sp, passf2_dp
  !  END INTERFACE passf2
  !  INTERFACE passf3
  !     MODULE PROCEDURE passf3_sp, passf3_dp
  !  END INTERFACE passf3
  !  INTERFACE passf4
  !     MODULE PROCEDURE passf4_sp, passf4_dp
  !  END INTERFACE passf4
  !  INTERFACE passf5
  !     MODULE PROCEDURE passf5_sp, passf5_dp
  !  END INTERFACE passf5
  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION cost_sp(para)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN) :: para
    REAL(sp)                           :: cost_sp

    ! local variables
    INTEGER(i4)                         :: n
    INTEGER(i4)                         :: i
    REAL(sp), DIMENSION(:), ALLOCATABLE :: model


    n = size(meas_sp,1)
    allocate(model(n))

    do i = 1, n
       model(i) = 1._sp/(para(1)**2*(1._sp+para(2)**2*real(meas_sp(i),sp)))
    end do

    cost_sp = 0._sp
    do i = 1, n
       cost_sp = cost_sp + abs(log(model(i))-log(real(meas_sp(i),sp)))
    end do

  END FUNCTION cost_sp

  FUNCTION cost_dp(para)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: para
    REAL(dp)                           :: cost_dp

    ! local variables
    INTEGER(i4)                         :: n
    INTEGER(i4)                         :: i
    REAL(dp), DIMENSION(:), ALLOCATABLE :: model

    n = size(meas_dp,1)
    allocate(model(n))

    do i = 1, n
       model(i) = 1._dp/(para(1)**2*(1._dp+para(2)**2*real(meas_dp(i),dp)))
    end do

    cost_dp = 0._dp
    do i = 1, n
       cost_dp = cost_dp + abs(log(model(i))-log(real(meas_dp(i),dp)))
    end do

  END FUNCTION cost_dp

  ! ------------------------------------------------------------------

  SUBROUTINE linres_sp(data,a_min,a_max,tc_min,tc_max,initial_a,initial_tc, bestpara)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN)  :: data
    REAL(sp), DIMENSION(2), INTENT(OUT) :: bestpara
    REAL(sp),               INTENT(IN)  :: a_min
    REAL(sp),               INTENT(IN)  :: a_max
    REAL(sp),               INTENT(IN)  :: tc_min
    REAL(sp),               INTENT(IN)  :: tc_max
    REAL(sp),               INTENT(IN)  :: initial_a
    REAL(sp),               INTENT(IN)  :: initial_tc

    ! local variables
    REAL(sp), DIMENSION(2,2)            :: pararanges
    REAL(sp), DIMENSION(2)              :: initial
    ! REAL(sp), DIMENSION(2)              :: para

    initial(1) = initial_a
    initial(2) = initial_tc

    pararanges(1,1) = a_min
    pararanges(1,2) = a_max
    pararanges(2,1) = tc_min
    pararanges(2,2) = tc_max

    call setmeas(data)

    bestpara = nelminrange(cost_sp, initial, pararanges)

  END SUBROUTINE linres_sp

  SUBROUTINE linres_dp(data,a_min,a_max,tc_min,tc_max,initial_a,initial_tc, bestpara)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN)  :: data
    REAL(dp), DIMENSION(2), INTENT(OUT) :: bestpara
    REAL(dp),               INTENT(IN)  :: a_min
    REAL(dp),               INTENT(IN)  :: a_max
    REAL(dp),               INTENT(IN)  :: tc_min
    REAL(dp),               INTENT(IN)  :: tc_max
    REAL(dp),               INTENT(IN)  :: initial_a
    REAL(dp),               INTENT(IN)  :: initial_tc

    ! local variables
    REAL(dp), DIMENSION(2,2)            :: pararanges
    REAL(dp), DIMENSION(2)              :: initial
    ! REAL(dp), DIMENSION(2)              :: para

    initial(1) = initial_a
    initial(2) = initial_tc

    pararanges(1,1) = a_min
    pararanges(1,2) = a_max
    pararanges(2,1) = tc_min
    pararanges(2,2) = tc_max

    call setmeas(data)

    bestpara = nelminrange(cost_dp, initial, pararanges)

  END SUBROUTINE linres_dp

  ! ------------------------------------------------------------------

  SUBROUTINE periodogram_sp(data,output)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),   INTENT(IN)  :: data
    REAL(sp), DIMENSION(:,:), INTENT(OUT) :: output

    ! local variables
    INTEGER(i4)                             :: N, i
    COMPLEX(spc), DIMENSION(:), ALLOCATABLE :: c_data
    ! REAL(sp)                                :: mean_data
    REAL(sp),     DIMENSION(:), ALLOCATABLE :: wsave
    REAL(sp),     DIMENSION(:), ALLOCATABLE :: Iw
    REAL(sp),     DIMENSION(:), ALLOCATABLE :: w

    N = size(data)
    allocate(c_data(N))
    allocate(wsave(4*N+15))

    do i = 1, N
       c_data(i) = cmplx(data(i), 0.0_sp,kind=sp) !- mean_data
    end do

    call cffti(N,wsave)
    call cfftf(N,c_data,wsave)

    allocate(Iw(N))
    allocate(w(N))

    do i = 1, N/2
       Iw(i) = 2*(abs(c_data(i))**2.0_sp)/real(N,sp)
       output(i,2) = Iw(i)
    end do
    do i = 1, N/2
       w(i) = 2*PI_sp*real(i,sp)/real(N,sp)
       output(i,1) = w(i)
    end do

  END SUBROUTINE periodogram_sp

  SUBROUTINE periodogram_dp(data,output)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),   INTENT(IN)  :: data
    REAL(dp), DIMENSION(:,:), INTENT(OUT) :: output

    ! local variables
    INTEGER(i4)                             :: N, i
    COMPLEX(dpc), DIMENSION(:), ALLOCATABLE :: c_data
    ! REAL(dp)                                :: mean_data
    REAL(dp),     DIMENSION(:), ALLOCATABLE :: wsave
    REAL(dp),     DIMENSION(:), ALLOCATABLE :: Iw
    REAL(dp),     DIMENSION(:), ALLOCATABLE :: w

    N = size(data)
    allocate(c_data(N))
    allocate(wsave(4*N+15))

    do i = 1, N
       c_data(i) = cmplx(data(i), 0.0_dp,kind=dp)
    end do

    call cffti(N,wsave)
    call cfftf(N,c_data,wsave)


    allocate(Iw(N))
    allocate(w(N))
    do i = 1, N/2
       Iw(i) = (abs(c_data(i))**2.0_dp)/real(N,dp)
       output(i,2) = Iw(i)
    end do
    do i = 1, N/2
       w(i) = 2*PI_dp*real(i,dp)/real(N,dp)
       output(i,1) = w(i)
    end do

  END SUBROUTINE periodogram_dp

  ! ------------------------------------------------------------------

  SUBROUTINE setmeas_dp(data)

    REAL(dp), DIMENSION(:), INTENT(IN) :: data

    allocate(meas_dp(size(data)))
    meas_dp(:) = data(:)

  END SUBROUTINE setmeas_dp

  SUBROUTINE setmeas_sp(data)

    REAL(sp), DIMENSION(:), INTENT(IN) :: data

    allocate(meas_sp(size(data)))  
    meas_sp(:) = data(:)

  END SUBROUTINE setmeas_sp

  ! ------------------------------------------------------------------
  !                  PRIVATE ROUTINES
  ! ------------------------------------------------------------------

  SUBROUTINE cfftf_sp ( n, c, wsave )
    !
    !*******************************************************************************
    !
    !! CFFTF computes the forward complex discrete Fourier transform.
    !
    !
    !  Discussion:
    !
    !    This process is sometimes called Fourier analysis.
    !
    !    CFFTF computes the Fourier coefficients of a complex periodic sequence.
    !
    !    The transform is not normalized.  To obtain a normalized transform,
    !    the output must be divided by N.  Otherwise a call of CFFTF
    !    followed by a call of CFFTB will multiply the sequence by N.
    !
    !    The array WSAVE must be initialized by calling CFFTI.
    !
    !    The transform is defined by:
    !
    !      C_out(J) = sum ( 1 <= K <= N ) 
    !        C_in(K) * exp ( - sqrt ( -1 ) * ( J - 1 ) * ( K - 1 ) * 2 * PI / N )
    !
    !  Modified:
    !
    !    09 March 2001
    !
    !  Author:
    !
    !    Paul Swarztrauber,
    !    National Center for Atmospheric Research
    !
    !  Reference:
    !
    !    David Kahaner, Clever Moler, Steven Nash,
    !    Numerical Methods and Software,
    !    Prentice Hall, 1988.
    !
    !    P N Swarztrauber, 
    !    Vectorizing the FFT's, 
    !    in Parallel Computations,
    !    G. Rodrigue, editor, 
    !    Academic Press, 1982, pages 51-83.
    !
    !    B L Buzbee, 
    !    The SLATEC Common Math Library, 
    !    in Sources and Development of Mathematical Software,
    !    W. Cowell, editor,
    !    Prentice Hall, 1984, pages 302-318.
    !
    !  Parameters:
    !
    !    Input, integer N, the length of the sequence to be transformed.  
    !    The method is more efficient when N is the product of small primes.
    !
    !    Input/output, complex C(N).
    !    On input, the data sequence to be transformed.
    !    On output, the Fourier coefficients.
    !
    !    Input, real WSAVE(4*N+15).  The array must be initialized by calling 
    !    CFFTI.  A different WSAVE array must be used for each different
    !    value of N. 
    !
    implicit none
    !
    integer(i4), intent(in)     :: n
    !
    complex(spc), intent(inout) :: c(n)
    real(sp), intent(in)        :: wsave(4*n+15)
    !********************
    integer(i4)                 :: i
    integer(i4)                 :: wsave_int(4*n+15)
    complex(spc)                :: wsave_complex(4*n+15)

    do i = 1, 4*n+15
       wsave_complex(i) = cmplx(wsave(i),0.0_sp,kind=sp)
       wsave_int(i) = int(wsave(i))
    end do
    !********************
    !
    if ( n <= 1 ) then
       return
    end if

    !call cfftf1 ( n, c, wsave(1), wsave(2*n+1), wsave(4*n+1) )
    call cfftf1_sp ( n, c, wsave_complex(1), wsave(2*n+1), wsave_int(4*n+1) )

    return
  END SUBROUTINE cfftf_sp

  SUBROUTINE cfftf_dp ( n, c, wsave )
    !
    implicit none
    !
    integer(i4), intent(in)     :: n
    !
    complex(dpc), intent(inout) :: c(n)
    real(dp), intent(in)        :: wsave(4*n+15)
    !********************
    integer(i4)                 :: i
    integer(i4)                 :: wsave_int(4*n+15)
    complex(dpc)                :: wsave_complex(4*n+15)

    do i = 1, 4*n+15
       wsave_complex(i) = cmplx(wsave(i),0.0_dp,kind=dp)
       wsave_int(i) = int(wsave(i))
    end do
    !********************
    !
    if ( n <= 1 ) then
       return
    end if

    !call cfftf1 ( n, c, wsave(1), wsave(2*n+1), wsave(4*n+1) )
    call cfftf1_dp ( n, c, wsave_complex(1), wsave(2*n+1), wsave_int(4*n+1) )

    return
  END SUBROUTINE cfftf_dp

  ! ------------------------------------------------------------------

  SUBROUTINE cfftf1_sp ( n, c, ch, wa, ifac )
    !
    !*******************************************************************************
    !
    !! CFFTF1 is a lower level routine used by CFFTF.
    !
    !
    !  Modified:
    !
    !    09 March 2001
    !
    !  Author:
    !
    !    Paul Swarztrauber,
    !    National Center for Atmospheric Research
    !
    !  Parameters:
    !
    !    Input, integer N, the length of the sequence to be transformed.  
    !
    !    Input/output, complex C(N).
    !    On input, the data sequence to be transformed.
    !    On output, the Fourier coefficients.
    !
    !    Input, complex CH(N).
    !
    !    Input, real WA(2*N).
    !
    !    Input, integer IFAC(15).
    !    IFAC(1) = N, the number that was factored.
    !    IFAC(2) = NF, the number of factors.
    !    IFAC(3:2+NF), the factors.
    !
    implicit none
    !
    integer(i4)  :: n
    !
    complex(spc) :: c(n)
    complex(spc) :: ch(n)
    integer(i4)  :: idl1
    integer(i4)  :: ido
    integer(i4)  :: ifac(15)
    integer(i4)  :: ip
    integer(i4)  :: iw
    integer(i4)  :: ix2
    integer(i4)  :: ix3
    integer(i4)  :: ix4
    integer(i4)  :: k1
    integer(i4)  :: l1
    integer(i4)  :: l2
    integer(i4)  :: na
    integer(i4)  :: nac
    integer(i4)  :: nf
    real(sp)     :: wa(2*n)

    !***************
    real(sp)     :: c_r(n*2)
    real(sp)     :: ch_r(n*2)
    integer(i4)  :: i, ii

    ii = 1
    do i = 1, n
       c_r(ii) = real(real(c(i),sp),sp)
       ch_r(ii) = real(real(ch(i),sp),sp)
       ii = ii + 1
       c_r(ii) = real(aimag(c(i)),sp)
       ch_r(ii) = real(aimag(ch(i)),sp)
       ii = ii + 1
    end do
    !***************
    !
    nf = ifac(2)
    na = 0
    l1 = 1
    iw = 1

    do k1 = 1, nf

       ip = ifac(k1+2)
       l2 = ip * l1
       ido = n / l2
       idl1 = 2 * ido * l1

       if ( ip == 4 ) then

          ix2 = iw + 2 * ido
          ix3 = ix2 + 2 * ido

          if ( na == 0 ) then
             !call passf4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
             call passf4_sp ( 2*ido, l1, c_r, ch_r, wa(iw), wa(ix2), wa(ix3) )
          else
             !call passf4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
             call passf4_sp ( 2*ido, l1, ch_r, c_r, wa(iw), wa(ix2), wa(ix3) )
          end if

          na = 1 - na

       else if ( ip == 2 ) then

          if ( na == 0 ) then
             !call passf2 ( 2*ido, l1, c, ch, wa(iw) )
             call passf2_sp ( 2*ido, l1, c_r, ch_r, wa(iw) )
          else
             !call passf2 ( 2*ido, l1, ch, c, wa(iw) )
             call passf2_sp ( 2*ido, l1, ch_r, c_r, wa(iw) )
          end if

          na = 1 - na

       else if ( ip == 3 ) then

          ix2 = iw + 2 * ido

          if ( na == 0 ) then
             !call passf3 ( 2*ido, l1, c, ch, wa(iw), wa(ix2) )
             call passf3_sp ( 2*ido, l1, c_r, ch_r, wa(iw), wa(ix2) )
          else
             !call passf3 ( 2*ido, l1, ch, c, wa(iw), wa(ix2) )
             call passf3_sp ( 2*ido, l1, ch_r, c_r, wa(iw), wa(ix2) )
          end if

          na = 1 - na

       else if ( ip == 5 ) then

          ix2 = iw + 2 * ido
          ix3 = ix2 + 2 * ido
          ix4 = ix3 + 2 * ido

          if ( na == 0 ) then
             !call passf5 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
             call passf5_sp ( 2*ido, l1, c_r, ch_r, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
          else
             !call passf5 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
             call passf5_sp ( 2*ido, l1, ch_r, c_r, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
          end if

          na = 1 - na

       else

          if ( na == 0 ) then
             !call passf ( nac, 2*ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
             call passf_sp ( nac, 2*ido, ip, l1, idl1, c_r, c_r, c_r, ch_r, ch_r, wa(iw) )
          else
             !call passf ( nac, 2*ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
             call passf_sp ( nac, 2*ido, ip, l1, idl1, ch_r, ch_r, ch_r, c_r, c_r, wa(iw) )
          end if

          if ( nac /= 0 ) then
             na = 1 - na
          end if

       end if

       l1 = l2
       iw = iw + ( ip - 1 ) * 2 * ido

    end do

    !*********************
    ii = 1
    do i = 1, n
       c(i) = cmplx(c_r(ii), c_r(ii+1),kind=sp)
       ch(i) = cmplx(ch_r(ii), ch_r(ii+1),kind=sp)
       ii = ii + 2
    end do
    !*********************

    if ( na /= 0 ) then
       c(1:n) = ch(1:n)
    end if
    return
  END SUBROUTINE cfftf1_sp

  SUBROUTINE cfftf1_dp ( n, c, ch, wa, ifac )
    !
    implicit none
    !
    integer(i4)  :: n
    !
    complex(dpc) :: c(n)
    complex(dpc) :: ch(n)
    integer(i4)  :: idl1
    integer(i4)  :: ido
    integer(i4)  :: ifac(15)
    integer(i4)  :: ip
    integer(i4)  :: iw
    integer(i4)  :: ix2
    integer(i4)  :: ix3
    integer(i4)  :: ix4
    integer(i4)  :: k1
    integer(i4)  :: l1
    integer(i4)  :: l2
    integer(i4)  :: na
    integer(i4)  :: nac
    integer(i4)  :: nf
    real(dp)     :: wa(2*n)

    !***************
    real(dp)     :: c_r(n*2)
    real(dp)     :: ch_r(n*2)
    integer(i4)  :: i, ii

    ii = 1
    do i = 1, n
       c_r(ii) = real(real(c(i),dp),dp)
       ch_r(ii) = real(real(ch(i),dp),dp)
       ii = ii + 1
       c_r(ii) = real(aimag(c(i)),dp)
       ch_r(ii) = real(aimag(ch(i)),dp)
       ii = ii + 1
    end do
    !***************
    !
    nf = ifac(2)
    na = 0
    l1 = 1
    iw = 1

    do k1 = 1, nf

       ip = ifac(k1+2)
       l2 = ip * l1
       ido = n / l2
       idl1 = 2 * ido * l1

       if ( ip == 4 ) then

          ix2 = iw + 2 * ido
          ix3 = ix2 + 2 * ido

          if ( na == 0 ) then
             !call passf4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
             call passf4_dp ( 2*ido, l1, c_r, ch_r, wa(iw), wa(ix2), wa(ix3) )
          else
             !call passf4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
             call passf4_dp ( 2*ido, l1, ch_r, c_r, wa(iw), wa(ix2), wa(ix3) )
          end if

          na = 1 - na

       else if ( ip == 2 ) then

          if ( na == 0 ) then
             !call passf2 ( 2*ido, l1, c, ch, wa(iw) )
             call passf2_dp ( 2*ido, l1, c_r, ch_r, wa(iw) )
          else
             !call passf2 ( 2*ido, l1, ch, c, wa(iw) )
             call passf2_dp ( 2*ido, l1, ch_r, c_r, wa(iw) )
          end if

          na = 1 - na

       else if ( ip == 3 ) then

          ix2 = iw + 2 * ido

          if ( na == 0 ) then
             !call passf3 ( 2*ido, l1, c, ch, wa(iw), wa(ix2) )
             call passf3_dp ( 2*ido, l1, c_r, ch_r, wa(iw), wa(ix2) )
          else
             !call passf3 ( 2*ido, l1, ch, c, wa(iw), wa(ix2) )
             call passf3_dp ( 2*ido, l1, ch_r, c_r, wa(iw), wa(ix2) )
          end if

          na = 1 - na

       else if ( ip == 5 ) then

          ix2 = iw + 2 * ido
          ix3 = ix2 + 2 * ido
          ix4 = ix3 + 2 * ido

          if ( na == 0 ) then
             !call passf5 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
             call passf5_dp ( 2*ido, l1, c_r, ch_r, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
          else
             !call passf5 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
             call passf5_dp ( 2*ido, l1, ch_r, c_r, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
          end if

          na = 1 - na

       else

          if ( na == 0 ) then
             !call passf ( nac, 2*ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
             call passf_dp ( nac, 2*ido, ip, l1, idl1, c_r, c_r, c_r, ch_r, ch_r, wa(iw) )
          else
             !call passf ( nac, 2*ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
             call passf_dp ( nac, 2*ido, ip, l1, idl1, ch_r, ch_r, ch_r, c_r, c_r, wa(iw) )
          end if

          if ( nac /= 0 ) then
             na = 1 - na
          end if

       end if

       l1 = l2
       iw = iw + ( ip - 1 ) * 2 * ido

    end do

    !*********************
    ii = 1
    do i = 1, n
       c(i) = cmplx(c_r(ii), c_r(ii+1),kind=dp)
       ch(i) = cmplx(ch_r(ii), ch_r(ii+1),kind=dp)
       ii = ii + 2
    end do
    !*********************

    if ( na /= 0 ) then
       c(1:n) = ch(1:n)
    end if
    return
  END SUBROUTINE cfftf1_dp

  ! ------------------------------------------------------------------

  SUBROUTINE cffti_sp ( n, wsave )
    !
    !*******************************************************************************
    !
    !! CFFTI initializes WSAVE, used in CFFTF and CFFTB. 
    !
    !
    !  Discussion:
    !
    !    The prime factorization of N together with a tabulation of the 
    !    trigonometric functions are computed and stored in WSAVE.
    !
    !  Modified:
    !
    !    09 March 2001
    !
    !  Author:
    !
    !    Paul Swarztrauber,
    !    National Center for Atmospheric Research
    !
    !  Reference:
    !
    !    David Kahaner, Clever Moler, Steven Nash,
    !    Numerical Methods and Software,
    !    Prentice Hall, 1988.
    !
    !    P N Swarztrauber, 
    !    Vectorizing the FFT's, 
    !    in Parallel Computations,
    !    G. Rodrigue, editor, 
    !    Academic Press, 1982, pages 51-83.
    !
    !    B L Buzbee, 
    !    The SLATEC Common Math Library, 
    !    in Sources and Development of Mathematical Software,
    !    W. Cowell, editor,
    !    Prentice Hall, 1984, pages 302-318.
    !
    !  Parameters:
    !
    !    Input, integer N, the length of the sequence to be transformed.
    !
    !    Output, real WSAVE(4*N+15), contains data, dependent on the value
    !    of N, which is necessary for the CFFTF or CFFTB routines.  
    !
    implicit none
    !
    integer(i4), intent(in) :: n
    !
    real(sp), intent(out)   :: wsave(4*n+15)

    !************************
    integer(i4)             :: wsave_int(4*n+15)
    integer(i4)             :: i

    wsave(:)  = 0.0_sp
    wsave_int = int(wsave)

    !************************
    !
    if ( n <= 1 ) then
       return
    end if
    !call cffti1 ( n, wsave(2*n+1), wsave(4*n+1) )
    call cffti1_sp ( n, wsave(2*n+1), wsave_int(4*n+1) )

    !************************
    do i = 4*n+1, 4*n+15
       wsave(i) = real(wsave_int(i),sp)
    end do

    return
  END SUBROUTINE cffti_sp

  SUBROUTINE cffti_dp ( n, wsave )
    !
    implicit none
    !
    integer(i4), intent(in) :: n
    !
    real(dp), intent(out)   :: wsave(4*n+15)

    !************************
    integer(i4) :: wsave_int(4*n+15)
    integer(i4) :: i

    wsave(:)  = 0.0_dp
    wsave_int = int(wsave)
    !************************
    !
    if ( n <= 1 ) then
       return
    end if
    !call cffti1 ( n, wsave(2*n+1), wsave(4*n+1) )
    call cffti1_dp ( n, wsave(2*n+1), wsave_int(4*n+1) )

    !************************
    do i = 4*n+1, 4*n+15
       wsave(i) = real(wsave_int(i),dp)
    end do

    return
  END SUBROUTINE cffti_dp

  ! ------------------------------------------------------------------

  SUBROUTINE cffti1_sp ( n, wa, ifac )
    !
    !*******************************************************************************
    !
    !! CFFTI1 is a lower level routine used by CFFTI.
    !
    !
    !  Modified:
    !
    !    09 March 2001
    !
    !  Author:
    !
    !    Paul Swarztrauber,
    !    National Center for Atmospheric Research
    !
    !  Parameters:
    !
    !    Input, integer N, the length of the sequence to be transformed.
    !
    !    Input, real WA(2*N).
    !
    !    Input, integer IFAC(15).
    !    IFAC(1) = N, the number that was factored.
    !    IFAC(2) = NF, the number of factors.
    !    IFAC(3:2+NF), the factors.
    !
    implicit none
    !
    integer(i4), intent(in)    :: n
    integer(i4), intent(inout) :: ifac(15)
    real(sp), intent(inout)    :: wa(2*n)
    !
    real(sp)                   :: arg
    real(sp)                   :: argh
    real(sp)                   :: argld
    real(sp)                   :: fi
    integer(i4)                :: i
    integer(i4)                :: i1
    !integer(i4)               :: ib
    integer(i4)                :: ido
    integer(i4)                :: ii
    integer(i4)                :: ip
    integer(i4)                :: j
    integer(i4)                :: k1
    integer(i4)                :: l1
    integer(i4)                :: l2
    integer(i4)                :: ld
    integer(i4)                :: nf
    !
    call i_factor ( n, ifac )

    nf = ifac(2)

    argh = 2.0E+00 * PI_sp / real ( n,sp )
    i = 2
    l1 = 1

    do k1 = 1, nf

       ip = ifac(k1+2)
       ld = 0
       l2 = l1 * ip
       ido = n / l2

       do j = 1, ip-1

          i1 = i
          wa(i-1) = 1.0E+00
          wa(i) = 0.0E+00
          ld = ld + l1
          fi = 0.0E+00
          argld = real ( ld,sp ) * argh

          do ii = 4, 2*ido+2, 2
             i = i + 2
             fi = fi + 1.0E+00
             arg = fi * argld
             wa(i-1) = cos ( arg )
             wa(i) = sin ( arg )
          end do

          if ( ip > 5 ) then
             wa(i1-1) = wa(i-1)
             wa(i1) = wa(i)
          end if

       end do

       l1 = l2

    end do
    return
  END SUBROUTINE cffti1_sp

  SUBROUTINE cffti1_dp ( n, wa, ifac )
    !
    implicit none
    !
    integer(i4), intent(in)    :: n
    real(dp), intent(inout)    :: wa(2*n)
    integer(i4), intent(inout) :: ifac(15)
    !
    real(dp)                   :: arg
    real(dp)                   :: argh
    real(dp)                   :: argld
    real(dp)                   :: fi
    integer(i4)                :: i
    integer(i4)                :: i1
    !integer(i4)               :: ib
    integer(i4)                :: ido
    integer(i4)                :: ii
    integer(i4)                :: ip
    integer(i4)                :: j
    integer(i4)                :: k1
    integer(i4)                :: l1
    integer(i4)                :: l2
    integer(i4)                :: ld
    integer(i4)                :: nf
    !
    call i_factor ( n, ifac )

    nf = ifac(2)

    argh = 2.0E+00 * PI_dp / real ( n,dp )
    i = 2
    l1 = 1

    do k1 = 1, nf

       ip = ifac(k1+2)
       ld = 0
       l2 = l1 * ip
       ido = n / l2

       do j = 1, ip-1

          i1 = i
          wa(i-1) = 1.0E+00
          wa(i) = 0.0E+00
          ld = ld + l1
          fi = 0.0E+00
          argld = real ( ld,dp ) * argh

          do ii = 4, 2*ido+2, 2
             i = i + 2
             fi = fi + 1.0E+00
             arg = fi * argld
             wa(i-1) = cos ( arg )
             wa(i) = sin ( arg )
          end do

          if ( ip > 5 ) then
             wa(i1-1) = wa(i-1)
             wa(i1) = wa(i)
          end if

       end do

       l1 = l2

    end do
    return
  END SUBROUTINE cffti1_dp

  ! ------------------------------------------------------------------

  SUBROUTINE i_factor ( n, ifac )
    !
    !*******************************************************************************
    !
    !! I_FACTOR factors an integer.
    !
    !
    !  Modified:
    !
    !    14 March 2001
    !
    !  Author:
    !
    !    Paul Swarztrauber,
    !    National Center for Atmospheric Research
    !
    !  Parameters:
    !
    !    Input, integer N, the number to be factored.
    !
    !    Output, integer IFAC(15).
    !    IFAC(1) = N, the number that was factored.
    !    IFAC(2) = NF, the number of factors.
    !    IFAC(3:2+NF), the factors.
    !
    implicit none
    !
    integer(i4) :: i
    ! integer(i4) :: ib
    integer(i4) :: ifac(15)
    integer(i4) :: j
    integer(i4) :: n
    integer(i4) :: nf
    integer(i4) :: nl
    integer(i4) :: nq
    integer(i4) :: nr
    integer(i4) :: ntry
    !
    ifac(1) = n

    nf = 0
    nl = n

    if ( n == 0 ) then
       nf = 1
       ifac(2) = nf
       ifac(2+nf) = 0
       return
    end if

    if ( n < 1 ) then
       nf = nf + 1
       ifac(2+nf) = -1
       nl = - n
    end if

    if ( nl == 1 ) then
       nf = nf + 1
       ifac(2) = nf
       ifac(2+nf) = 1
       return
    end if

    j = 0

    do while ( nl > 1 )

       j = j + 1
       !
       !  Choose a trial divisor, NTRY.
       !
       if ( j == 1 ) then
          ntry = 4
       else if ( j == 2 ) then
          ntry = 2
       else if ( j == 3 ) then
          ntry = 3
       else if ( j == 4 ) then
          ntry = 5
       else
          ntry = ntry + 2
       end if
       !
       !  Divide by the divisor as many times as possible.
       !
       do

          nq = nl / ntry
          nr = nl - ntry * nq

          if ( nr /= 0 ) then
             exit
          end if

          nl = nq
          nf = nf + 1
          !
          !  Make sure factors of 2 appear in the front of the list.
          !
          if ( ntry /= 2 ) then

             ifac(2+nf) = ntry

          else

             do i = nf, 2, -1
                ifac(i+2) = ifac(i+1)
             end do
             ifac(3) = 2

          end if

       end do

    end do

    ifac(2) = nf

    return
  END SUBROUTINE i_factor

  SUBROUTINE passf_sp ( nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
    !
    !*******************************************************************************
    !
    !! PASSF is a lower level routine used by CFFTF1.
    !
    !
    !  Modified:
    !
    !    09 March 2001
    !
    !  Author:
    !
    !    Paul Swarztrauber,
    !    National Center for Atmospheric Research
    !
    !  Parameters:
    !
    implicit none
    !
    integer(i4)  :: idl1
    integer(i4)  :: ido
    integer(i4)  :: ip
    integer(i4)  :: l1
    !
    real(sp)     :: c1(ido,l1,ip)
    real(sp)     :: c2(idl1,ip)
    real(sp)     :: cc(ido,ip,l1)
    real(sp)     :: ch(ido,l1,ip)
    real(sp)     :: ch2(idl1,ip)
    integer(i4)  :: i
    integer(i4)  :: idij
    integer(i4)  :: idj
    integer(i4)  :: idl
    integer(i4)  :: idlj
    integer(i4)  :: idp
    integer(i4)  :: ik
    integer(i4)  :: inc
    integer(i4)  :: ipph
    integer(i4)  :: j
    integer(i4)  :: jc
    integer(i4)  :: k
    integer(i4)  :: l
    integer(i4)  :: lc
    integer(i4)  :: nac
    !integer(i4) :: nt
    real(sp)     :: wa(*)
    real(sp)     :: wai
    real(sp)     :: war
    !
    !nt = ip * idl1
    ipph = (ip+1) / 2
    idp = ip * ido

    if ( ido >= l1 ) then

       do j = 2, ipph
          jc = ip + 2 - j
          ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
          ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
       end do

       ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

    else

       do j = 2, ipph
          jc = ip + 2 - j
          ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
          ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
       end do

       ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

    end if

    idl = 2 - ido
    inc = 0

    do l = 2, ipph

       lc = ip + 2 - l
       idl = idl + ido

       do ik = 1, idl1
          c2(ik,l)  = ch2(ik,1) + wa(idl-1) * ch2(ik,2)
          c2(ik,lc) =           - wa(idl)   * ch2(ik,ip)
       end do

       idlj = idl
       inc = inc + ido

       do j = 3, ipph

          jc = ip + 2 - j

          idlj = idlj + inc
          if ( idlj > idp ) then
             idlj = idlj - idp
          end if

          war = wa(idlj-1)
          wai = wa(idlj)

          do ik = 1, idl1
             c2(ik,l)  = c2(ik,l)  + war * ch2(ik,j)
             c2(ik,lc) = c2(ik,lc) - wai * ch2(ik,jc)
          end do

       end do

    end do

    do j = 2, ipph
       ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
    end do

    do j = 2, ipph
       jc = ip + 2 - j
       do ik = 2, idl1, 2
          ch2(ik-1,j)  = c2(ik-1,j) - c2(ik,jc)
          ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
          ch2(ik,j)    = c2(ik,j)   + c2(ik-1,jc)
          ch2(ik,jc)   = c2(ik,j)   - c2(ik-1,jc)
       end do
    end do

    if ( ido == 2 ) then
       nac = 1
       return
    end if

    nac = 0

    c2(1:idl1,1)    = ch2(1:idl1,1)
    c1(1,1:l1,2:ip) = ch(1,1:l1,2:ip)
    c1(2,1:l1,2:ip) = ch(2,1:l1,2:ip)

    if ( ( ido / 2 ) <= l1 ) then

       idij = 0
       do j = 2, ip
          idij = idij + 2
          do i = 4, ido, 2
             idij = idij + 2
             c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) + wa(idij) * ch(i,1:l1,j)
             c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   - wa(idij) * ch(i-1,1:l1,j)
          end do
       end do

    else

       idj = 2 - ido

       do j = 2, ip
          idj = idj + ido
          do k = 1, l1
             idij = idj
             do i = 4, ido, 2
                idij = idij + 2
                c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) + wa(idij) * ch(i,k,j)
                c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   - wa(idij) * ch(i-1,k,j)
             end do
          end do
       end do

    end if

    return
  END SUBROUTINE passf_sp

  SUBROUTINE passf_dp ( nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
    !
    implicit none
    !
    integer(i4)  :: idl1
    integer(i4)  :: ido
    integer(i4)  :: ip
    integer(i4)  :: l1
    !
    real(dp)     :: c1(ido,l1,ip)
    real(dp)     :: c2(idl1,ip)
    real(dp)     :: cc(ido,ip,l1)
    real(dp)     :: ch(ido,l1,ip)
    real(dp)     :: ch2(idl1,ip)
    integer(i4)  :: i
    integer(i4)  :: idij
    integer(i4)  :: idj
    integer(i4)  :: idl
    integer(i4)  :: idlj
    integer(i4)  :: idp
    integer(i4)  :: ik
    integer(i4)  :: inc
    integer(i4)  :: ipph
    integer(i4)  :: j
    integer(i4)  :: jc
    integer(i4)  :: k
    integer(i4)  :: l
    integer(i4)  :: lc
    integer(i4)  :: nac
    !integer(i4) :: nt
    real(dp)     :: wa(*)
    real(dp)     :: wai
    real(dp)     :: war
    !
    !nt = ip * idl1
    ipph = (ip+1) / 2
    idp = ip * ido

    if ( ido >= l1 ) then

       do j = 2, ipph
          jc = ip + 2 - j
          ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
          ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
       end do

       ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

    else

       do j = 2, ipph
          jc = ip + 2 - j
          ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
          ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
       end do

       ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

    end if

    idl = 2 - ido
    inc = 0

    do l = 2, ipph

       lc = ip + 2 - l
       idl = idl + ido

       do ik = 1, idl1
          c2(ik,l)  = ch2(ik,1) + wa(idl-1) * ch2(ik,2)
          c2(ik,lc) =           - wa(idl)   * ch2(ik,ip)
       end do

       idlj = idl
       inc = inc + ido

       do j = 3, ipph

          jc = ip + 2 - j

          idlj = idlj + inc
          if ( idlj > idp ) then
             idlj = idlj - idp
          end if

          war = wa(idlj-1)
          wai = wa(idlj)

          do ik = 1, idl1
             c2(ik,l)  = c2(ik,l)  + war * ch2(ik,j)
             c2(ik,lc) = c2(ik,lc) - wai * ch2(ik,jc)
          end do

       end do

    end do

    do j = 2, ipph
       ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
    end do

    do j = 2, ipph
       jc = ip + 2 - j
       do ik = 2, idl1, 2
          ch2(ik-1,j)  = c2(ik-1,j) - c2(ik,jc)
          ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
          ch2(ik,j)    = c2(ik,j)   + c2(ik-1,jc)
          ch2(ik,jc)   = c2(ik,j)   - c2(ik-1,jc)
       end do
    end do

    if ( ido == 2 ) then
       nac = 1
       return
    end if

    nac = 0

    c2(1:idl1,1)    = ch2(1:idl1,1)
    c1(1,1:l1,2:ip) = ch(1,1:l1,2:ip)
    c1(2,1:l1,2:ip) = ch(2,1:l1,2:ip)

    if ( ( ido / 2 ) <= l1 ) then

       idij = 0
       do j = 2, ip
          idij = idij + 2
          do i = 4, ido, 2
             idij = idij + 2
             c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) + wa(idij) * ch(i,1:l1,j)
             c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   - wa(idij) * ch(i-1,1:l1,j)
          end do
       end do

    else

       idj = 2 - ido

       do j = 2, ip
          idj = idj + ido
          do k = 1, l1
             idij = idj
             do i = 4, ido, 2
                idij = idij + 2
                c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) + wa(idij) * ch(i,k,j)
                c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   - wa(idij) * ch(i-1,k,j)
             end do
          end do
       end do

    end if

    return
  END SUBROUTINE passf_dp

  ! ------------------------------------------------------------------

  SUBROUTINE passf2_sp ( ido, l1, cc, ch, wa1 )
    !
    !*******************************************************************************
    !
    !! PASSF2 is a lower level routine used by CFFTF1.
    !
    !
    !  Modified:
    !
    !    09 March 2001
    !
    !  Author:
    !
    !    Paul Swarztrauber,
    !    National Center for Atmospheric Research
    !
    !  Parameters:
    !
    implicit none
    !
    integer(i4) :: ido
    integer(i4) :: l1
    !
    real(sp)    :: cc(ido,2,l1)
    real(sp)    :: ch(ido,l1,2)
    integer(i4) :: i
    integer(i4) :: k
    real(sp)    :: ti2
    real(sp)    :: tr2
    real(sp)    :: wa1(ido)
    !
    if ( ido <= 2 ) then

       ch(1,1:l1,1) = cc(1,1,1:l1) + cc(1,2,1:l1)
       ch(1,1:l1,2) = cc(1,1,1:l1) - cc(1,2,1:l1)
       ch(2,1:l1,1) = cc(2,1,1:l1) + cc(2,2,1:l1)
       ch(2,1:l1,2) = cc(2,1,1:l1) - cc(2,2,1:l1)

    else

       do k = 1, l1
          do i = 2, ido, 2

             ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
             tr2         = cc(i-1,1,k) - cc(i-1,2,k)

             ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
             ti2       = cc(i,1,k) - cc(i,2,k)

             ch(i,k,2)   = wa1(i-1) * ti2 - wa1(i) * tr2
             ch(i-1,k,2) = wa1(i-1) * tr2 + wa1(i) * ti2

          end do
       end do

    end if
    return
  END SUBROUTINE passf2_sp

  SUBROUTINE passf2_dp ( ido, l1, cc, ch, wa1 )
    !
    implicit none
    !
    integer(i4) :: ido
    integer(i4) :: l1
    !
    real(dp)    :: cc(ido,2,l1)
    real(dp)    :: ch(ido,l1,2)
    integer(i4) :: i
    integer(i4) :: k
    real(dp)    :: ti2
    real(dp)    :: tr2
    real(dp)    :: wa1(ido)
    !
    if ( ido <= 2 ) then

       ch(1,1:l1,1) = cc(1,1,1:l1) + cc(1,2,1:l1)
       ch(1,1:l1,2) = cc(1,1,1:l1) - cc(1,2,1:l1)
       ch(2,1:l1,1) = cc(2,1,1:l1) + cc(2,2,1:l1)
       ch(2,1:l1,2) = cc(2,1,1:l1) - cc(2,2,1:l1)

    else

       do k = 1, l1
          do i = 2, ido, 2

             ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
             tr2         = cc(i-1,1,k) - cc(i-1,2,k)

             ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
             ti2       = cc(i,1,k) - cc(i,2,k)

             ch(i,k,2)   = wa1(i-1) * ti2 - wa1(i) * tr2
             ch(i-1,k,2) = wa1(i-1) * tr2 + wa1(i) * ti2

          end do
       end do

    end if
    return
  END SUBROUTINE passf2_dp

  ! ------------------------------------------------------------------

  SUBROUTINE passf3_sp ( ido, l1, cc, ch, wa1, wa2 )
    !
    !*******************************************************************************
    !
    !! PASSF3 is a lower level routine used by CFFTF1.
    !
    !
    !  Modified:
    !
    !    09 March 2001
    !
    !  Author:
    !
    !    Paul Swarztrauber,
    !    National Center for Atmospheric Research
    !
    !  Parameters:
    !
    implicit none
    !
    integer(i4)         :: ido
    integer(i4)         :: l1
    !
    real(sp)            :: cc(ido,3,l1)
    real(sp)            :: ch(ido,l1,3)
    real(sp)            :: ci2
    real(sp)            :: ci3
    real(sp)            :: cr2
    real(sp)            :: cr3
    real(sp)            :: di2
    real(sp)            :: di3
    real(sp)            :: dr2
    real(sp)            :: dr3
    integer(i4)         :: i
    integer(i4)         :: k
    real(sp)            :: taui
    real(sp), parameter :: taur = -0.5E+00
    real(sp)            :: ti2
    real(sp)            :: tr2
    real(sp)            :: wa1(ido)
    real(sp)            :: wa2(ido)
    !
    taui = - sqrt ( 3.0E+00 ) / 2.0E+00

    if ( ido == 2 ) then

       do k = 1, l1

          tr2 = cc(1,2,k) + cc(1,3,k)
          cr2 = cc(1,1,k) + taur * tr2
          ch(1,k,1) = cc(1,1,k) + tr2

          ti2 = cc(2,2,k) + cc(2,3,k)
          ci2 = cc(2,1,k) + taur * ti2
          ch(2,k,1) = cc(2,1,k) + ti2

          cr3 = taui * ( cc(1,2,k) - cc(1,3,k) )
          ci3 = taui * ( cc(2,2,k) - cc(2,3,k) )

          ch(1,k,2) = cr2 - ci3
          ch(1,k,3) = cr2 + ci3
          ch(2,k,2) = ci2 + cr3
          ch(2,k,3) = ci2 - cr3

       end do

    else

       do k = 1, l1
          do i = 2, ido, 2

             tr2 = cc(i-1,2,k) + cc(i-1,3,k)
             cr2 = cc(i-1,1,k) + taur * tr2
             ch(i-1,k,1) = cc(i-1,1,k) + tr2

             ti2 = cc(i,2,k) + cc(i,3,k)
             ci2 = cc(i,1,k) + taur * ti2
             ch(i,k,1) = cc(i,1,k) + ti2

             cr3 = taui * ( cc(i-1,2,k) - cc(i-1,3,k) )
             ci3 = taui * ( cc(i,2,k)   - cc(i,3,k) )

             dr2 = cr2 - ci3
             dr3 = cr2 + ci3
             di2 = ci2 + cr3
             di3 = ci2 - cr3

             ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
             ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
             ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
             ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3

          end do
       end do

    end if
    return
  END SUBROUTINE passf3_sp

  SUBROUTINE passf3_dp ( ido, l1, cc, ch, wa1, wa2 )
    !
    implicit none
    !
    integer(i4)         :: ido
    integer(i4)         :: l1
    !
    real(dp)            :: cc(ido,3,l1)
    real(dp)            :: ch(ido,l1,3)
    real(dp)            :: ci2
    real(dp)            :: ci3
    real(dp)            :: cr2
    real(dp)            :: cr3
    real(dp)            :: di2
    real(dp)            :: di3
    real(dp)            :: dr2
    real(dp)            :: dr3
    integer(i4)         :: i
    integer(i4)         :: k
    real(dp)            :: taui
    real(dp), parameter :: taur = -0.5E+00
    real(dp)            :: ti2
    real(dp)            :: tr2
    real(dp)            :: wa1(ido)
    real(dp)            :: wa2(ido)
    !
    taui = - sqrt ( 3.0E+00 ) / 2.0E+00

    if ( ido == 2 ) then

       do k = 1, l1

          tr2 = cc(1,2,k) + cc(1,3,k)
          cr2 = cc(1,1,k) + taur * tr2
          ch(1,k,1) = cc(1,1,k) + tr2

          ti2 = cc(2,2,k) + cc(2,3,k)
          ci2 = cc(2,1,k) + taur * ti2
          ch(2,k,1) = cc(2,1,k) + ti2

          cr3 = taui * ( cc(1,2,k) - cc(1,3,k) )
          ci3 = taui * ( cc(2,2,k) - cc(2,3,k) )

          ch(1,k,2) = cr2 - ci3
          ch(1,k,3) = cr2 + ci3
          ch(2,k,2) = ci2 + cr3
          ch(2,k,3) = ci2 - cr3

       end do

    else

       do k = 1, l1
          do i = 2, ido, 2

             tr2 = cc(i-1,2,k) + cc(i-1,3,k)
             cr2 = cc(i-1,1,k) + taur * tr2
             ch(i-1,k,1) = cc(i-1,1,k) + tr2

             ti2 = cc(i,2,k) + cc(i,3,k)
             ci2 = cc(i,1,k) + taur * ti2
             ch(i,k,1) = cc(i,1,k) + ti2

             cr3 = taui * ( cc(i-1,2,k) - cc(i-1,3,k) )
             ci3 = taui * ( cc(i,2,k)   - cc(i,3,k) )

             dr2 = cr2 - ci3
             dr3 = cr2 + ci3
             di2 = ci2 + cr3
             di3 = ci2 - cr3

             ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
             ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
             ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
             ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3

          end do
       end do

    end if
    return
  END SUBROUTINE passf3_dp

  ! ------------------------------------------------------------------

  SUBROUTINE passf4_sp ( ido, l1, cc, ch, wa1, wa2, wa3 )
    !
    !*******************************************************************************
    !
    !! PASSF4 is a lower level routine used by CFFTF1.
    !
    !
    !  Modified:
    !
    !    09 March 2001
    !
    !  Author:
    !
    !    Paul Swarztrauber,
    !    National Center for Atmospheric Research
    !
    !  Parameters:
    !
    implicit none
    !
    integer(i4) :: ido
    integer(i4) :: l1
    !
    real(sp)    :: cc(ido,4,l1)
    real(sp)    :: ch(ido,l1,4)
    !  real(sp) :: ci1
    real(sp)    :: ci2
    real(sp)    :: ci3
    real(sp)    :: ci4
    ! real(sp)  :: cr1
    real(sp)    :: cr2
    real(sp)    :: cr3
    real(sp)    :: cr4
    integer(i4) :: i
    integer(i4) :: k
    real(sp)    :: ti1
    real(sp)    :: ti2
    real(sp)    :: ti3
    real(sp)    :: ti4
    real(sp)    :: tr1
    real(sp)    :: tr2
    real(sp)    :: tr3
    real(sp)    :: tr4
    real(sp)    :: wa1(ido)
    real(sp)    :: wa2(ido)
    real(sp)    :: wa3(ido)
    !
    if ( ido == 2 ) then

       do k = 1, l1

          ti1 = cc(2,1,k) - cc(2,3,k)
          ti2 = cc(2,1,k) + cc(2,3,k)
          tr4 = cc(2,2,k) - cc(2,4,k)
          ti3 = cc(2,2,k) + cc(2,4,k)
          tr1 = cc(1,1,k) - cc(1,3,k)
          tr2 = cc(1,1,k) + cc(1,3,k)
          ti4 = cc(1,4,k) - cc(1,2,k)
          tr3 = cc(1,2,k) + cc(1,4,k)

          ch(1,k,1) = tr2 + tr3
          ch(1,k,3) = tr2 - tr3
          ch(2,k,1) = ti2 + ti3
          ch(2,k,3) = ti2 - ti3
          ch(1,k,2) = tr1 + tr4
          ch(1,k,4) = tr1 - tr4
          ch(2,k,2) = ti1 + ti4
          ch(2,k,4) = ti1 - ti4

       end do

    else

       do k = 1, l1
          do i = 2, ido, 2

             ti1 = cc(i,1,k)   - cc(i,3,k)
             ti2 = cc(i,1,k)   + cc(i,3,k)
             ti3 = cc(i,2,k)   + cc(i,4,k)
             tr4 = cc(i,2,k)   - cc(i,4,k)
             tr1 = cc(i-1,1,k) - cc(i-1,3,k)
             tr2 = cc(i-1,1,k) + cc(i-1,3,k)
             ti4 = cc(i-1,4,k) - cc(i-1,2,k)
             tr3 = cc(i-1,2,k) + cc(i-1,4,k)

             ch(i-1,k,1) = tr2 + tr3
             cr3         = tr2 - tr3
             ch(i,k,1)   = ti2 + ti3
             ci3         = ti2 - ti3

             cr2 = tr1 + tr4
             cr4 = tr1 - tr4
             ci2 = ti1 + ti4
             ci4 = ti1 - ti4

             ch(i-1,k,2) = wa1(i-1) * cr2 + wa1(i) * ci2
             ch(i,k,2)   = wa1(i-1) * ci2 - wa1(i) * cr2
             ch(i-1,k,3) = wa2(i-1) * cr3 + wa2(i) * ci3
             ch(i,k,3)   = wa2(i-1) * ci3 - wa2(i) * cr3
             ch(i-1,k,4) = wa3(i-1) * cr4 + wa3(i) * ci4
             ch(i,k,4)   = wa3(i-1) * ci4 - wa3(i) * cr4

          end do
       end do

    end if
    return
  END SUBROUTINE passf4_sp

  SUBROUTINE passf4_dp ( ido, l1, cc, ch, wa1, wa2, wa3 )
    !
    implicit none
    !
    integer(i4) :: ido
    integer(i4) :: l1
    !
    real(dp)    :: cc(ido,4,l1)
    real(dp)    :: ch(ido,l1,4)
    ! real(dp)  :: ci1
    real(dp)    :: ci2
    real(dp)    :: ci3
    real(dp)    :: ci4
    ! real(dp)  :: cr1
    real(dp)    :: cr2
    real(dp)    :: cr3
    real(dp)    :: cr4
    integer(i4) :: i
    integer(i4) :: k
    real(dp)    :: ti1
    real(dp)    :: ti2
    real(dp)    :: ti3
    real(dp)    :: ti4
    real(dp)    :: tr1
    real(dp)    :: tr2
    real(dp)    :: tr3
    real(dp)    :: tr4
    real(dp)    :: wa1(ido)
    real(dp)    :: wa2(ido)
    real(dp)    :: wa3(ido)
    !
    if ( ido == 2 ) then

       do k = 1, l1

          ti1 = cc(2,1,k) - cc(2,3,k)
          ti2 = cc(2,1,k) + cc(2,3,k)
          tr4 = cc(2,2,k) - cc(2,4,k)
          ti3 = cc(2,2,k) + cc(2,4,k)
          tr1 = cc(1,1,k) - cc(1,3,k)
          tr2 = cc(1,1,k) + cc(1,3,k)
          ti4 = cc(1,4,k) - cc(1,2,k)
          tr3 = cc(1,2,k) + cc(1,4,k)

          ch(1,k,1) = tr2 + tr3
          ch(1,k,3) = tr2 - tr3
          ch(2,k,1) = ti2 + ti3
          ch(2,k,3) = ti2 - ti3
          ch(1,k,2) = tr1 + tr4
          ch(1,k,4) = tr1 - tr4
          ch(2,k,2) = ti1 + ti4
          ch(2,k,4) = ti1 - ti4

       end do

    else

       do k = 1, l1
          do i = 2, ido, 2

             ti1 = cc(i,1,k)   - cc(i,3,k)
             ti2 = cc(i,1,k)   + cc(i,3,k)
             ti3 = cc(i,2,k)   + cc(i,4,k)
             tr4 = cc(i,2,k)   - cc(i,4,k)
             tr1 = cc(i-1,1,k) - cc(i-1,3,k)
             tr2 = cc(i-1,1,k) + cc(i-1,3,k)
             ti4 = cc(i-1,4,k) - cc(i-1,2,k)
             tr3 = cc(i-1,2,k) + cc(i-1,4,k)

             ch(i-1,k,1) = tr2 + tr3
             cr3         = tr2 - tr3
             ch(i,k,1)   = ti2 + ti3
             ci3         = ti2 - ti3

             cr2 = tr1 + tr4
             cr4 = tr1 - tr4
             ci2 = ti1 + ti4
             ci4 = ti1 - ti4

             ch(i-1,k,2) = wa1(i-1) * cr2 + wa1(i) * ci2
             ch(i,k,2)   = wa1(i-1) * ci2 - wa1(i) * cr2
             ch(i-1,k,3) = wa2(i-1) * cr3 + wa2(i) * ci3
             ch(i,k,3)   = wa2(i-1) * ci3 - wa2(i) * cr3
             ch(i-1,k,4) = wa3(i-1) * cr4 + wa3(i) * ci4
             ch(i,k,4)   = wa3(i-1) * ci4 - wa3(i) * cr4

          end do
       end do

    end if
    return
  END SUBROUTINE passf4_dp

  ! ------------------------------------------------------------------

  SUBROUTINE passf5_sp ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
    !
    !*******************************************************************************
    !
    !! PASSF5 is a lower level routine used by CFFTF1.
    !
    !
    !  Modified:
    !
    !    09 March 2001
    !
    !  Author:
    !
    !    Paul Swarztrauber,
    !    National Center for Atmospheric Research
    !
    !  Parameters:
    !
    implicit none
    !
    integer(i4) ido
    integer(i4) l1
    !
    real(sp) :: cc(ido,5,l1)
    real(sp) :: ch(ido,l1,5)
    real(sp) :: ci2
    real(sp) :: ci3
    real(sp) :: ci4
    real(sp) :: ci5
    real(sp) :: cr2
    real(sp) :: cr3
    real(sp) :: cr4
    real(sp) :: cr5
    real(sp) :: di2
    real(sp) :: di3
    real(sp) :: di4
    real(sp) :: di5
    real(sp) :: dr2
    real(sp) :: dr3
    real(sp) :: dr4
    real(sp) :: dr5
    integer(i4) i
    integer(i4) k
    real(sp), parameter :: ti11 = -0.951056516295154E+00
    real(sp), parameter :: ti12 = -0.587785252292473E+00
    real(sp) :: ti2
    real(sp) :: ti3
    real(sp) :: ti4
    real(sp) :: ti5
    real(sp) :: tr2
    real(sp) :: tr3
    real(sp) :: tr4
    real(sp) :: tr5
    !
    !  cos ( 72 ) = +0.3090
    !
    real(sp), parameter :: tr11 =  0.309016994374947E+00
    !
    !  cos ( 36 ) = +0.809016
    !
    real(sp), parameter :: tr12 = -0.809016994374947E+00
    real(sp) :: wa1(ido)
    real(sp) :: wa2(ido)
    real(sp) :: wa3(ido)
    real(sp) :: wa4(ido)
    !
    if ( ido == 2 ) then

       do k = 1, l1

          ti5 = cc(2,2,k) - cc(2,5,k)
          ti2 = cc(2,2,k) + cc(2,5,k)
          ti4 = cc(2,3,k) - cc(2,4,k)
          ti3 = cc(2,3,k) + cc(2,4,k)
          tr5 = cc(1,2,k) - cc(1,5,k)
          tr2 = cc(1,2,k) + cc(1,5,k)
          tr4 = cc(1,3,k) - cc(1,4,k)
          tr3 = cc(1,3,k) + cc(1,4,k)

          ch(1,k,1) = cc(1,1,k) + tr2 + tr3
          ch(2,k,1) = cc(2,1,k) + ti2 + ti3

          cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
          ci2 = cc(2,1,k) + tr11 * ti2 + tr12 * ti3
          cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
          ci3 = cc(2,1,k) + tr12 * ti2 + tr11 * ti3

          cr5 = ti11 * tr5 + ti12 * tr4
          ci5 = ti11 * ti5 + ti12 * ti4
          cr4 = ti12 * tr5 - ti11 * tr4
          ci4 = ti12 * ti5 - ti11 * ti4

          ch(1,k,2) = cr2 - ci5
          ch(1,k,5) = cr2 + ci5
          ch(2,k,2) = ci2 + cr5
          ch(2,k,3) = ci3 + cr4
          ch(1,k,3) = cr3 - ci4
          ch(1,k,4) = cr3 + ci4
          ch(2,k,4) = ci3 - cr4
          ch(2,k,5) = ci2 - cr5

       end do

    else

       do k = 1, l1
          do i = 2, ido, 2

             ti5 = cc(i,2,k) - cc(i,5,k)
             ti2 = cc(i,2,k) + cc(i,5,k)
             ti4 = cc(i,3,k) - cc(i,4,k)
             ti3 = cc(i,3,k) + cc(i,4,k)

             tr5 = cc(i-1,2,k) - cc(i-1,5,k)
             tr2 = cc(i-1,2,k) + cc(i-1,5,k)
             tr4 = cc(i-1,3,k) - cc(i-1,4,k)
             tr3 = cc(i-1,3,k) + cc(i-1,4,k)

             ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
             ch(i,k,1)   = cc(i,1,k)   + ti2 + ti3

             cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
             ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
             cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
             ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

             cr5 = ti11 * tr5 + ti12 * tr4
             ci5 = ti11 * ti5 + ti12 * ti4
             cr4 = ti12 * tr5 - ti11 * tr4
             ci4 = ti12 * ti5 - ti11 * ti4

             dr3 = cr3 - ci4
             dr4 = cr3 + ci4
             di3 = ci3 + cr4
             di4 = ci3 - cr4
             dr5 = cr2 + ci5
             dr2 = cr2 - ci5
             di5 = ci2 - cr5
             di2 = ci2 + cr5

             ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
             ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
             ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3
             ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
             ch(i-1,k,4) = wa3(i-1) * dr4 + wa3(i) * di4
             ch(i,k,4)   = wa3(i-1) * di4 - wa3(i) * dr4
             ch(i-1,k,5) = wa4(i-1) * dr5 + wa4(i) * di5
             ch(i,k,5)   = wa4(i-1) * di5 - wa4(i) * dr5

          end do
       end do

    end if
    return
  END SUBROUTINE passf5_sp

  SUBROUTINE passf5_dp ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
    !
    implicit none
    !
    integer(i4) ido
    integer(i4) l1
    !
    real(dp) :: cc(ido,5,l1)
    real(dp) :: ch(ido,l1,5)
    real(dp) :: ci2
    real(dp) :: ci3
    real(dp) :: ci4
    real(dp) :: ci5
    real(dp) :: cr2
    real(dp) :: cr3
    real(dp) :: cr4
    real(dp) :: cr5
    real(dp) :: di2
    real(dp) :: di3
    real(dp) :: di4
    real(dp) :: di5
    real(dp) :: dr2
    real(dp) :: dr3
    real(dp) :: dr4
    real(dp) :: dr5
    integer(i4) i
    integer(i4) k
    real(dp), parameter :: ti11 = -0.951056516295154E+00
    real(dp), parameter :: ti12 = -0.587785252292473E+00
    real(dp) :: ti2
    real(dp) :: ti3
    real(dp) :: ti4
    real(dp) :: ti5
    real(dp) :: tr2
    real(dp) :: tr3
    real(dp) :: tr4
    real(dp) :: tr5
    !
    !  cos ( 72 ) = +0.3090
    !
    real(dp), parameter :: tr11 =  0.309016994374947E+00
    !
    !  cos ( 36 ) = +0.809016
    !
    real(dp), parameter :: tr12 = -0.809016994374947E+00
    real(dp) :: wa1(ido)
    real(dp) :: wa2(ido)
    real(dp) :: wa3(ido)
    real(dp) :: wa4(ido)
    !
    if ( ido == 2 ) then

       do k = 1, l1

          ti5 = cc(2,2,k) - cc(2,5,k)
          ti2 = cc(2,2,k) + cc(2,5,k)
          ti4 = cc(2,3,k) - cc(2,4,k)
          ti3 = cc(2,3,k) + cc(2,4,k)
          tr5 = cc(1,2,k) - cc(1,5,k)
          tr2 = cc(1,2,k) + cc(1,5,k)
          tr4 = cc(1,3,k) - cc(1,4,k)
          tr3 = cc(1,3,k) + cc(1,4,k)

          ch(1,k,1) = cc(1,1,k) + tr2 + tr3
          ch(2,k,1) = cc(2,1,k) + ti2 + ti3

          cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
          ci2 = cc(2,1,k) + tr11 * ti2 + tr12 * ti3
          cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
          ci3 = cc(2,1,k) + tr12 * ti2 + tr11 * ti3

          cr5 = ti11 * tr5 + ti12 * tr4
          ci5 = ti11 * ti5 + ti12 * ti4
          cr4 = ti12 * tr5 - ti11 * tr4
          ci4 = ti12 * ti5 - ti11 * ti4

          ch(1,k,2) = cr2 - ci5
          ch(1,k,5) = cr2 + ci5
          ch(2,k,2) = ci2 + cr5
          ch(2,k,3) = ci3 + cr4
          ch(1,k,3) = cr3 - ci4
          ch(1,k,4) = cr3 + ci4
          ch(2,k,4) = ci3 - cr4
          ch(2,k,5) = ci2 - cr5

       end do

    else

       do k = 1, l1
          do i = 2, ido, 2

             ti5 = cc(i,2,k) - cc(i,5,k)
             ti2 = cc(i,2,k) + cc(i,5,k)
             ti4 = cc(i,3,k) - cc(i,4,k)
             ti3 = cc(i,3,k) + cc(i,4,k)

             tr5 = cc(i-1,2,k) - cc(i-1,5,k)
             tr2 = cc(i-1,2,k) + cc(i-1,5,k)
             tr4 = cc(i-1,3,k) - cc(i-1,4,k)
             tr3 = cc(i-1,3,k) + cc(i-1,4,k)

             ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
             ch(i,k,1)   = cc(i,1,k)   + ti2 + ti3

             cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
             ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
             cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
             ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

             cr5 = ti11 * tr5 + ti12 * tr4
             ci5 = ti11 * ti5 + ti12 * ti4
             cr4 = ti12 * tr5 - ti11 * tr4
             ci4 = ti12 * ti5 - ti11 * ti4

             dr3 = cr3 - ci4
             dr4 = cr3 + ci4
             di3 = ci3 + cr4
             di4 = ci3 - cr4
             dr5 = cr2 + ci5
             dr2 = cr2 - ci5
             di5 = ci2 - cr5
             di2 = ci2 + cr5

             ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
             ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
             ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3
             ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
             ch(i-1,k,4) = wa3(i-1) * dr4 + wa3(i) * di4
             ch(i,k,4)   = wa3(i-1) * di4 - wa3(i) * dr4
             ch(i-1,k,5) = wa4(i-1) * dr5 + wa4(i) * di5
             ch(i,k,5)   = wa4(i-1) * di5 - wa4(i) * dr5

          end do
       end do

    end if
    return
  END SUBROUTINE passf5_dp

END MODULE mo_specan

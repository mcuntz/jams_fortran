MODULE mo_percentile

  ! This module provides routines for median and percentiles.

  ! Literature
  !   WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !       Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !       Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996

  ! Written March 2011, Matthias Cuntz

  USE mo_kind, ONLY: i4, sp, dp, spc, dpc

  Implicit NONE

  PRIVATE

  PUBLIC :: ksmallest       ! Returns the kth smallest value in an array (median if k=N/2)
  PUBLIC :: median          ! Returns the median
  PUBLIC :: percentile      ! Returns the value below which a certain percent of array values fall.
                            ! (median if k=50)
  PUBLIC :: qmedian         ! Quick median calculation (k=N/2)

  ! Public
  INTERFACE ksmallest
     MODULE PROCEDURE ksmallest_sp, ksmallest_dp
  END INTERFACE ksmallest
  INTERFACE median
     MODULE PROCEDURE median_sp, median_dp
  END INTERFACE median
  INTERFACE percentile
     MODULE PROCEDURE percentile_0d_sp, percentile_0d_dp
     MODULE PROCEDURE percentile_1d_sp, percentile_1d_dp
  END INTERFACE percentile
  INTERFACE qmedian
     MODULE PROCEDURE qmedian_sp, qmedian_dp
  END INTERFACE qmedian

  ! Private
  INTERFACE swap
     MODULE PROCEDURE swap_i4, &
          swap_sp, swap_1d_sp, &
          swap_dp, swap_1d_dp, &
          swap_spc, swap_1d_spc, &
          swap_dpc, swap_1d_dpc, &
          masked_swap_sp, masked_swap_1d_sp, masked_swap_2d_sp, &
          masked_swap_dp, masked_swap_1d_dp, masked_swap_2d_dp, &
          masked_swap_spc, masked_swap_1d_spc, masked_swap_2d_spc, &
          masked_swap_dpc, masked_swap_1d_dpc, masked_swap_2d_dpc
  END INTERFACE swap

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         ksmallest

  !     PURPOSE
  !         Returns the kth smallest value in an array.
  !         If n>100, k=N/2 is the median for N even or odd.
  !         We do not rearrange the input array as in numerical recipes routine select.
  !
  !         If an optinal mask is given, values only on those locations that correspond
  !         to true values in the mask are used.

  !     CALLING SEQUENCE
  !         out = ksmallest(vec,k,mask=mask,before=before)
  
  !     INTENT(IN)
  !         real(sp/dp) :: vec(:)     1D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: out        k-th smallest value in input array

  !     INTENT(IN), OPTIONAL
  !         integer(i4) :: k          index of sorted array
  !         logical     :: mask(:)    1D-array of logical values with size(vec).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         real(sp/dp) :: before      (k-1)-th smallest value in input array for median/percentile calculations

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. /)
  !         ! Returns 4
  !         out = ksmallest(vec,4)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Mar 2011

  FUNCTION ksmallest_dp(arrin,k,mask,before)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN)  :: arrin
    INTEGER(i4),                         INTENT(IN)  :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp),                  OPTIONAL, INTENT(OUT) :: before
    REAL(dp)                                         :: ksmallest_dp

    INTEGER(i4) :: i, r, j, l, n
    REAL(dp)    :: a
    REAL(dp), DIMENSION(:), ALLOCATABLE :: arr

    if (present(mask)) then
       n = count(mask)
       allocate(arr(n))
       arr = pack(arrin,mask)
    else
       n = size(arrin)
       allocate(arr(n))
       arr = arrin
    endif

    if (k < 1) stop 'ksmallest_dp: k < 1'
    if (k > n) stop 'ksmallest_dp: k > n'

    l = 1
    r = n
    do
       if (r-l <= 1) then
          if (r-l == 1) call swap(arr(l), arr(r), arr(l)>arr(r))
          ksmallest_dp = arr(k)
          if (present(before)) before = maxval(arr(:k-1))
          deallocate(arr)
          RETURN
       else
          i = (l+r)/2
          call swap(arr(i),   arr(l+1))
          call swap(arr(l),   arr(r),   arr(l)>arr(r))
          call swap(arr(l+1), arr(r),   arr(l+1)>arr(r))
          call swap(arr(l),   arr(l+1), arr(l)>arr(l+1))
          i = l+1
          j = r
          a = arr(l+1)
          do
             do
                i = i+1
                if (arr(i) >= a) exit
             end do
             do
                j = j-1
                if (arr(j) <= a) exit
             end do
             if (j < i) exit
             call swap(arr(i), arr(j))
          end do
          arr(l+1) = arr(j)
          arr(j)   = a
          if (j >= k) r = j-1
          if (j <= k) l = i
       end if
    end do
    ksmallest_dp = arr(k)
    if (present(before)) before = maxval(arr(:k-1))

    deallocate(arr)

  END FUNCTION ksmallest_dp


  FUNCTION ksmallest_sp(arrin,k,mask,before)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN)  :: arrin
    INTEGER(i4),                         INTENT(IN)  :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp),                  OPTIONAL, INTENT(OUT) :: before
    REAL(sp)                                         :: ksmallest_sp

    INTEGER(i4) :: i, r, j, l, n
    REAL(sp)    :: a
    REAL(sp), DIMENSION(:), ALLOCATABLE :: arr

    if (present(mask)) then
       n = count(mask)
       allocate(arr(n))
       arr = pack(arrin,mask)
    else
       n = size(arrin)
       allocate(arr(n))
       arr = arrin
    endif

    if (k < 1) stop 'ksmallest_sp: k < 1'
    if (k > n) stop 'ksmallest_sp: k > n'

    l = 1
    r = n
    do
       if (r-l <= 1) then
          if (r-l == 1) call swap(arr(l), arr(r), arr(l)>arr(r))
          ksmallest_sp = arr(k)
          if (present(before)) before = maxval(arr(:k-1))
          deallocate(arr)
          RETURN
       else
          i = (l+r)/2
          call swap(arr(i),   arr(l+1))
          call swap(arr(l),   arr(r),   arr(l)>arr(r))
          call swap(arr(l+1), arr(r),   arr(l+1)>arr(r))
          call swap(arr(l),   arr(l+1), arr(l)>arr(l+1))
          i = l+1
          j = r
          a = arr(l+1)
          do
             do
                i = i+1
                if (arr(i) >= a) exit
             end do
             do
                j = j-1
                if (arr(j) <= a) exit
             end do
             if (j < i) exit
             call swap(arr(i), arr(j))
          end do
          arr(l+1) = arr(j)
          arr(j)   = a
          if (j >= k) r = j-1
          if (j <= k) l = i
       end if
    end do
    ksmallest_sp = arr(k)
    if (present(before)) before = maxval(arr(:k-1))

    deallocate(arr)

  END FUNCTION ksmallest_sp

  ! ------------------------------------------------------------------

  !     NAME
  !         median

  !     PURPOSE
  !         Returns the median of the values in an array.
  !         If size is even, then the mean of the size/2 and size/2+1 element is the median.
  !
  !         If an optinal mask is given, values only on those locations that correspond
  !         to true values in the mask are used.

  !     CALLING SEQUENCE
  !         out = median(vec,mask=mask)
  
  !     INTENT(IN)
  !         real(sp/dp) :: vec(:)     1D-array with input numbers

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: out        median of values in input array

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)    1D-array of logical values with size(vec).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. /)
  !         ! Returns 5.5
  !         out = median(vec)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Mar 2011

  FUNCTION median_dp(arrin,mask)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN) :: arrin
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(dp)                                        :: median_dp

    INTEGER(i4) :: n
    REAL(dp), DIMENSION(:), ALLOCATABLE :: arr
    REAL(dp) :: tmp

    if (present(mask)) then
       n = count(mask)
       allocate(arr(n))
       arr = pack(arrin,mask)

       if (n < 2) stop 'median_dp: n < 2'
    
       if (mod(n,2) == 0) then ! Even
          median_dp = ksmallest(arr,n/2+1,before=tmp)
          median_dp = 0.5_dp*(median_dp+tmp)
       else ! Odd
          median_dp = ksmallest(arr,(n+1)/2)
       endif

       deallocate(arr)
    else
       n = size(arrin)
       if (n < 2) stop 'median_dp: n < 2'
    
       if (mod(n,2) == 0) then ! Even
          median_dp = ksmallest(arrin,n/2+1,before=tmp)
          median_dp = 0.5_dp*(median_dp+tmp)
       else ! Odd
          median_dp = ksmallest(arrin,(n+1)/2)
       endif
    endif

  END FUNCTION median_dp


  FUNCTION median_sp(arrin,mask)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN) :: arrin
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(sp)                                        :: median_sp

    INTEGER(i4) :: n
    REAL(sp), DIMENSION(:), ALLOCATABLE :: arr
    REAL(sp) :: tmp

    if (present(mask)) then
       n = count(mask)
       allocate(arr(n))
       arr = pack(arrin,mask)

       if (n < 2) stop 'median_sp: n < 2'
    
       if (mod(n,2) == 0) then ! Even
          median_sp = ksmallest(arr,n/2+1,before=tmp)
          median_sp = 0.5_sp*(median_sp+tmp)
       else ! Odd
          median_sp = ksmallest(arr,(n+1)/2)
       endif

       deallocate(arr)
    else
       n = size(arrin)
       if (n < 2) stop 'median_sp: n < 2'
    
       if (mod(n,2) == 0) then ! Even
          median_sp = ksmallest(arrin,n/2+1,before=tmp)
          median_sp = 0.5_sp*(median_sp+tmp)
       else ! Odd
          median_sp = ksmallest(arrin,(n+1)/2)
       endif
    endif

  END FUNCTION median_sp

  ! ------------------------------------------------------------------

  !     NAME
  !         percentile

  !     PURPOSE
  !         Returns the value below which a certain percent of array values fall.
  !
  !         If an optinal mask is given, values only on those locations that correspond
  !         to true values in the mask are used.
  !
  !         Different definitions can be applied to interpolate the stepwise CDF of the given data.
  !         (1) Inverse empirical CDF (no interpolation, default MATHEMATICA)
  !         (2) Linear interpolation (California method)
  !         (3) Element numbered closest
  !         (4) Linear interpolation (hydrologist method)
  !         (5) Mean-based estimate (Weibull method, default IMSL)
  !         (6) Mode-based estimate
  !         (7) Median-based estimate
  !         (8) normal distribution estimate
  !         
  !         See: http://reference.wolfram.com/mathematica/tutorial/BasicStatistics.html

  !     CALLING SEQUENCE
  !         out = percentile(vec,k,mask=mask,mode_in=mode)
  
  !     INTENT(IN)
  !         real(sp/dp) :: vec(:)     1D-array with input numbers
  !         real(sp/dp) :: k[(:)]     Percentage of percentile, can be 1 dimensional 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         real(sp/dp) :: out[(size(k))]   k-th percentile of values in input array, can be 
  !                                         1 dimensional corresponding to k

  !     INTENT(IN), OPTIONAL
  !         logical     :: mask(:)    1D-array of logical values with size(vec).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.
  !         integer(i4) :: mode_in    Specifies the interpolation scheme applied.
  !                                   Default: 
  !                                       Inverse empirical CDF (no interpolation, default Mathematica)
  !                                       mode_in = 1_i4

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. /)
  !         ! Returns 10.
  !         out = percentile(vec,95.)
  !         ! Returns (10.,8)
  !         out = percentile(vec,(/95.,80./))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Mar  2011
  !         Modified, Stephan Thober, Dec  2011 - added 1 dimensional version
  !                   Juliane Mai,    July 2012 - different interpolation schemes

  FUNCTION percentile_0d_dp(arrin,k,mask,mode_in)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN) :: arrin
    REAL(dp),                            INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    INTEGER(i4),               OPTIONAL, INTENT(IN) :: mode_in                
    REAL(dp)                                        :: percentile_0d_dp

    INTEGER(i4)                         :: n, nn1, nn2
    INTEGER(i4)                         :: mode
    REAL(dp)                            :: kk, ks1, ks2
    REAL(dp), DIMENSION(:), ALLOCATABLE :: arr

    if (present(mask)) then
       n = count(mask)
    else
       n = size(arrin)
    endif

    if (present(mode_in)) then
       mode = mode_in
    else
       ! Default : Inverse empirical CDF
       mode = 1_i4
    end if

    if (n < 2) stop 'percentile_0d_dp: n < 2'
    
    select case (mode)
       ! Inverse empirical CDF: Mathematica default
       case(1_i4) 
          kk = k/100._dp*real(n,dp)
          nn1 = min(n, max(1_i4,ceiling(kk,kind=i4)))
          nn2 = nn1

       ! Linear interpolation (California method)
       case(2_i4)
          kk  = k/100._dp*real(n,dp)
          nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
          nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Element numbered closest
       case(3_i4)
          kk = 0.5_dp+k/100._dp*real(n,dp)
          nn1 = min(n, max(1_i4,floor(kk,kind=i4)))
          nn2 = nn1

       ! Linear interpolation (hydrologist method)
       case(4_i4)
          kk  = 0.5_dp+k/100._dp*(real(n,dp))
          nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
          nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Mean-based estimate (Weibull method): IMSL default
       case(5_i4)
          kk  = k/100._dp*(real(n,dp)+1._dp)
          nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
          nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Mode-based estimate
       case(6_i4)
          kk  = 1.0_dp+k/100._dp*(real(n,dp)-1._dp)
          nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
          nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Median-based estimate
       case(7_i4)
          kk  = 1.0_dp/3.0_dp+k/100._dp*(real(n,dp)+1.0_dp/3.0_dp)
          nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
          nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Normal distribution estimate
       case(8_i4)
          kk  = 3.0_dp/8.0_dp+k/100._dp*(real(n,dp)+1.0_dp/4.0_dp)
          nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
          nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! No valid mode
       case default
          stop 'percentile_0d_dp: mode > 8 not implemented'

    end select

    if (present(mask)) then
       allocate(arr(n))
       arr = pack(arrin,mask)
       if (nn1 .eq. nn2) then
          ! no interpolation
          percentile_0d_dp = ksmallest(arr,nn1)
       else
          ! interpolation
          ks2 = ksmallest(arr,nn2,before=ks1)
          percentile_0d_dp = ks1 + (ks2-ks1)*(kk-real(nn1,dp))
       end if
       deallocate(arr)
    else
       if (nn1 .eq. nn2) then
          ! no interpolation
          percentile_0d_dp = ksmallest(arrin,nn1)
       else
          ! interpolation
          ks2 = ksmallest(arrin,nn2,before=ks1)
          percentile_0d_dp = ks1 + (ks2-ks1)*(kk-real(nn1,dp))
       end if
    endif

  END FUNCTION percentile_0d_dp


  FUNCTION percentile_0d_sp(arrin,k,mask,mode_in)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN) :: arrin
    REAL(sp),                            INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    INTEGER(i4),               OPTIONAL, INTENT(IN) :: mode_in                
    REAL(sp)                                        :: percentile_0d_sp

    INTEGER(i4)                         :: n, nn1, nn2
    INTEGER(i4)                         :: mode
    REAL(sp)                            :: kk, ks1, ks2
    REAL(sp), DIMENSION(:), ALLOCATABLE :: arr

    if (present(mask)) then
       n = count(mask)
    else
       n = size(arrin)
    endif

    if (present(mode_in)) then
       mode = mode_in
    else
       ! Default : Inverse empirical CDF
       mode = 1_i4
    end if

    if (n < 2) stop 'percentile_0d_sp: n < 2'
    
    select case (mode)
       ! Inverse empirical CDF: Mathematica default
       case(1_i4) 
          kk = k/100._sp*real(n,sp)
          nn1 = min(n, max(1_i4,ceiling(kk,kind=i4)))
          nn2 = nn1

       ! Linear interpolation (California method)
       case(2_i4)
          kk  = k/100._sp*real(n,sp)
          nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
          nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Element numbered closest
       case(3_i4)
          kk = 0.5_sp+k/100._sp*real(n,sp)
          nn1 = min(n, max(1_i4,floor(kk,kind=i4)))
          nn2 = nn1

       ! Linear interpolation (hydrologist method)
       case(4_i4)
          kk  = 0.5_sp+k/100._sp*(real(n,sp))
          nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
          nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Mean-based estimate (Weibull method): IMSL default
       case(5_i4)
          kk  = k/100._sp*(real(n,sp)+1._sp)
          nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
          nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Mode-based estimate
       case(6_i4)
          kk  = 1.0_sp+k/100._sp*(real(n,sp)-1._sp)
          nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
          nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Median-based estimate
       case(7_i4)
          kk  = 1.0_sp/3.0_sp+k/100._sp*(real(n,sp)+1.0_sp/3.0_sp)
          nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
          nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! Normal distribution estimate
       case(8_i4)
          kk  = 3.0_sp/8.0_sp+k/100._sp*(real(n,sp)+1.0_sp/4.0_sp)
          nn1 = min(n, max(1_i4,  floor(kk,kind=i4)))
          nn2 = min(n, max(1_i4,ceiling(kk,kind=i4)))

       ! No valid mode
       case default
          stop 'percentile_0d_sp: mode > 8 not implemented'

    end select

    if (present(mask)) then
       allocate(arr(n))
       arr = pack(arrin,mask)
       if (nn1 .eq. nn2) then
          ! no interpolation
          percentile_0d_sp = ksmallest(arr,nn1)
       else
          ! interpolation
          ks2 = ksmallest(arr,nn2,before=ks1)
          percentile_0d_sp = ks1 + (ks2-ks1)*(kk-real(nn1,sp))
       end if
       deallocate(arr)
    else
       if (nn1 .eq. nn2) then
          ! no interpolation
          percentile_0d_sp = ksmallest(arrin,nn1)
       else
          ! interpolation
          ks2 = ksmallest(arrin,nn2,before=ks1)
          percentile_0d_sp = ks1 + (ks2-ks1)*(kk-real(nn1,sp))
       end if
    endif

  END FUNCTION percentile_0d_sp


  function percentile_1d_dp(arrin,k,mask,mode_in)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN) :: arrin
    REAL(dp),    DIMENSION(:),           INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    INTEGER(i4),               OPTIONAL, INTENT(IN) :: mode_in

    REAL(dp),    DIMENSION(size(k))                 :: percentile_1d_dp

    INTEGER(i4)                            :: i, n
    INTEGER(i4)                            :: mode
    INTEGER(i4), DIMENSION(size(k))        :: nn1, nn2
    REAL(dp),    DIMENSION(size(k))        :: kk
    REAL(dp)                               :: ks1, ks2
    REAL(dp),    DIMENSION(:), ALLOCATABLE :: arr

    if (present(mask)) then
       n = count(mask)
    else
       n = size(arrin)
    endif

    if (present(mode_in)) then
       mode = mode_in
    else
       ! Default : Inverse empirical CDF
       mode = 1_i4
    end if

    ! check consistency
    !if (size(k) > size(arr)) stop 'percentile_1d_dp: more Quantiles than data: size(k) > size(arr)'
    if (n < 2) stop 'percentile_1d_dp: n < 2'

    select case (mode)
       ! Inverse empirical CDF: Mathematica default
       case(1_i4)   
          kk(:) = k(:)/100._dp*real(n,dp)
          nn1(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))
          nn2 = nn1

      ! Linear interpolation (California method)
       case(2_i4)
          kk(:)  = k(:)/100._dp*real(n,dp)
          nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
          nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Element numbered closest
       case(3_i4)
          kk(:) = 0.5_dp+k(:)/100._dp*real(n,dp)
          nn1(:) = min(n, max(1_i4,floor(kk(:),kind=i4)))
          nn2 = nn1

       ! Linear interpolation (hydrologist method)
       case(4_i4)
          kk(:)  = 0.5_dp+k(:)/100._dp*(real(n,dp))
          nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
          nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Mean-based estimate (Weibull method): IMSL default
       case(5_i4)
          kk(:)  = k(:)/100._dp*(real(n,dp)+1._dp)
          nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
          nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Mode-based estimate
       case(6_i4)
          kk(:)  = 1.0_dp+k(:)/100._dp*(real(n,dp)-1._dp)
          nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
          nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Median-based estimate
       case(7_i4)
          kk(:)  = 1.0_dp/3.0_dp+k(:)/100._dp*(real(n,dp)+1.0_dp/3.0_dp)
          nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
          nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Normal distribution estimate
       case(8_i4)
          kk(:)  = 3.0_dp/8.0_dp+k(:)/100._dp*(real(n,dp)+1.0_dp/4.0_dp)
          nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
          nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! No valid mode
       case default
          stop 'percentile_1d_dp: mode > 8 not implemented'

    end select

    if (present(mask)) then
       allocate(arr(n))
       arr = pack(arrin,mask)
       do i=1, size(k)
          if (nn1(i) .eq. nn2(i)) then
             ! no interpolation
             percentile_1d_dp(i) = ksmallest(arr,nn1(i))
          else
             ! interpolation
             ks2 = ksmallest(arr,nn2(i),before=ks1)
             percentile_1d_dp(i) = ks1 + (ks2-ks1)*(kk(i)-real(nn1(i),dp))
          end if
       end do
       deallocate(arr)
    else
       do i=1, size(k)
          if (nn1(i) .eq. nn2(i)) then
             ! no interpolation
             percentile_1d_dp(i) = ksmallest(arrin,nn1(i))
          else
             ! interpolation
             ks2 = ksmallest(arrin,nn2(i),before=ks1)
             percentile_1d_dp(i) = ks1 + (ks2-ks1)*(kk(i)-real(nn1(i),dp))
          end if
       end do
    endif

  END function percentile_1d_dp


  function percentile_1d_sp(arrin,k,mask,mode_in)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN) :: arrin
    REAL(sp),    DIMENSION(:),           INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    INTEGER(i4),               OPTIONAL, INTENT(IN) :: mode_in

    REAL(sp),    DIMENSION(size(k))                 :: percentile_1d_sp

    INTEGER(i4)                            :: i, n
    INTEGER(i4)                            :: mode
    INTEGER(i4), DIMENSION(size(k))        :: nn1, nn2
    REAL(sp),    DIMENSION(size(k))        :: kk
    REAL(sp)                               :: ks1, ks2
    REAL(sp),    DIMENSION(:), ALLOCATABLE :: arr

    if (present(mask)) then
       n = count(mask)
    else
       n = size(arrin)
    endif

    if (present(mode_in)) then
       mode = mode_in
    else
       ! Default : Inverse empirical CDF
       mode = 1_i4
    end if

    ! check consistency
    !if (size(k) > size(arr)) stop 'percentile_1d_sp: more Quantiles than data: size(k) > size(arr)'
    if (n < 2) stop 'percentile_1d_sp: n < 2'

    select case (mode)
       ! Inverse empirical CDF: Mathematica default
       case(1_i4)   
          kk(:) = k(:)/100._sp*real(n,sp)
          nn1(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))
          nn2 = nn1

      ! Linear interpolation (California method)
       case(2_i4)
          kk(:)  = k(:)/100._sp*real(n,sp)
          nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
          nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Element numbered closest
       case(3_i4)
          kk(:) = 0.5_sp+k(:)/100._sp*real(n,sp)
          nn1(:) = min(n, max(1_i4,floor(kk(:),kind=i4)))
          nn2 = nn1

       ! Linear interpolation (hydrologist method)
       case(4_i4)
          kk(:)  = 0.5_sp+k(:)/100._sp*(real(n,sp))
          nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
          nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Mean-based estimate (Weibull method): IMSL default
       case(5_i4)
          kk(:)  = k(:)/100._sp*(real(n,sp)+1._sp)
          nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
          nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Mode-based estimate
       case(6_i4)
          kk(:)  = 1.0_sp+k(:)/100._sp*(real(n,sp)-1._sp)
          nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
          nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Median-based estimate
       case(7_i4)
          kk(:)  = 1.0_sp/3.0_sp+k(:)/100._sp*(real(n,sp)+1.0_sp/3.0_sp)
          nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
          nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! Normal distribution estimate
       case(8_i4)
          kk(:)  = 3.0_sp/8.0_sp+k(:)/100._sp*(real(n,sp)+1.0_sp/4.0_sp)
          nn1(:) = min(n, max(1_i4,  floor(kk(:),kind=i4)))
          nn2(:) = min(n, max(1_i4,ceiling(kk(:),kind=i4)))

       ! No valid mode
       case default
          stop 'percentile_1d_sp: mode > 8 not implemented'

    end select

    if (present(mask)) then
       allocate(arr(n))
       arr = pack(arrin,mask)
       do i=1, size(k)
          if (nn1(i) .eq. nn2(i)) then
             ! no interpolation
             percentile_1d_sp(i) = ksmallest(arr,nn1(i))
          else
             ! interpolation
             ks2 = ksmallest(arr,nn2(i),before=ks1)
             percentile_1d_sp(i) = ks1 + (ks2-ks1)*(kk(i)-real(nn1(i),sp))
          end if
       end do
       deallocate(arr)
    else
       do i=1, size(k)
          if (nn1(i) .eq. nn2(i)) then
             ! no interpolation
             percentile_1d_sp(i) = ksmallest(arrin,nn1(i))
          else
             ! interpolation
             ks2 = ksmallest(arrin,nn2(i),before=ks1)
             percentile_1d_sp(i) = ks1 + (ks2-ks1)*(kk(i)-real(nn1(i),sp))
          end if
       end do
    endif

  END function percentile_1d_sp

  ! ------------------------------------------------------------------

  !     NAME
  !         qmedian

  !     PURPOSE
  !         Returns the quickly calculated median by overwriting the input array.
  !         It takes the size/2-th element as the median.

  !     CALLING SEQUENCE
  !         out = qmedian(vec)
  
  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         real(sp/dp) :: vec(:)     1D-array with input numbers

  !     INTENT(OUT)
  !         real(sp/dp) :: out        median of values in input array

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. /)
  !         ! Returns 5.
  !         out = qmedian(vec)
  !         -> see also example in test directory

  !     LITERATURE
  !         Niklaus Wirth. "Algorithms and Data Structures". Prentice-Hall, Inc., 1985. ISBN 0-13-022005-1.

  !     HISTORY
  !         Written, Filip Hroch as part of Munipack: http://munipack.physics.muni.cz
  !         Modified, Matthias Cuntz, Jul 2012 - function, k=n/2

  function qmedian_dp(dat)

    IMPLICIT NONE

    real(dp), dimension(:), intent(inout) :: dat
    real(dp) :: qmedian_dp

    real(dp) :: w
    integer(i4) :: n,k
    integer(i4) :: l,r,i,j

    logical :: iseven

    n = size(dat)
    k = n/2_i4+1_i4
    l = 1_i4
    r = n
    do while( l < r )
       qmedian_dp = dat(k)
       i = l
       j = r
       do
          do while( dat(i) < qmedian_dp )
             i = i + 1_i4
          enddo
          do while( qmedian_dp < dat(j) )
             j = j - 1_i4
          enddo
          if ( i <= j ) then
             w      = dat(i)
             dat(i) = dat(j)
             dat(j) = w
             i = i + 1_i4
             j = j - 1_i4
          endif
          if ( i > j ) exit
       enddo
       if ( j < k ) l = i
       if ( k < i ) r = j
    enddo
    if (mod(n,2_i4) == 0_i4) then
       qmedian_dp = 0.5_dp*(dat(k) + maxval(dat(:k-1)))
    else
       qmedian_dp = dat(k)
    endif

  end function qmedian_dp

  function qmedian_sp(dat)

    IMPLICIT NONE

    real(sp), dimension(:), intent(inout) :: dat
    real(sp) :: qmedian_sp

    real(sp) :: w
    integer(i4) :: n,k
    integer(i4) :: l,r,i,j

    n = size(dat)
    k = n/2
    l = 1
    r = n
    do while( l < r )
       qmedian_sp = dat(k)
       i = l
       j = r
       do
          do while( dat(i) < qmedian_sp )
             i = i + 1
          enddo
          do while( qmedian_sp < dat(j) )
             j = j - 1
          enddo
          if( i <= j ) then
             w = dat(i)
             dat(i) = dat(j)
             dat(j) = w
             i = i + 1
             j = j - 1
          endif
          if( i > j ) exit
       enddo
       if( j < k ) l = i
       if( k < i ) r = j
    enddo

  end function qmedian_sp

  ! ------------------------------------------------------------------

  SUBROUTINE swap_i4(a,b)
    INTEGER(i4), INTENT(INOUT) :: a,b
    INTEGER(i4) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i4


  SUBROUTINE swap_sp(a,b)
    REAL(sp), INTENT(INOUT) :: a,b
    REAL(sp) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_sp


  SUBROUTINE masked_swap_sp(a,b,mask)
    REAL(sp), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    REAL(sp) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_sp


  SUBROUTINE masked_swap_1d_sp(a,b,mask)
    REAL(sp), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    REAL(sp), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_1d_sp


  SUBROUTINE masked_swap_2d_sp(a,b,mask)
    REAL(sp), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
    REAL(sp), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_2d_sp


  SUBROUTINE swap_1d_sp(a,b)
    REAL(sp), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(sp), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_1d_sp


  SUBROUTINE swap_dp(a,b)
    REAL(dp), INTENT(INOUT) :: a,b
    REAL(dp) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_dp


  SUBROUTINE swap_1d_dp(a,b)
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(dp), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_1d_dp


  SUBROUTINE masked_swap_dp(a,b,mask)
    REAL(dp), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    REAL(dp) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_dp


  SUBROUTINE masked_swap_1d_dp(a,b,mask)
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    REAL(dp), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_1d_dp


  SUBROUTINE masked_swap_2d_dp(a,b,mask)
    REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
    REAL(dp), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_2d_dp


  SUBROUTINE swap_spc(a,b)
    COMPLEX(spc), INTENT(INOUT) :: a,b
    COMPLEX(spc) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_spc


  SUBROUTINE masked_swap_spc(a,b,mask)
    COMPLEX(spc), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    COMPLEX(spc) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_spc


  SUBROUTINE masked_swap_1d_spc(a,b,mask)
    COMPLEX(spc), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    COMPLEX(spc), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_1d_spc


  SUBROUTINE masked_swap_2d_spc(a,b,mask)
    COMPLEX(spc), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
    COMPLEX(spc), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_2d_spc


  SUBROUTINE swap_1d_spc(a,b)
    COMPLEX(spc), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(spc), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_1d_spc


  SUBROUTINE swap_dpc(a,b)
    COMPLEX(dpc), INTENT(INOUT) :: a,b
    COMPLEX(dpc) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_dpc


  SUBROUTINE masked_swap_dpc(a,b,mask)
    COMPLEX(dpc), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    COMPLEX(dpc) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_dpc


  SUBROUTINE masked_swap_1d_dpc(a,b,mask)
    COMPLEX(dpc), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    COMPLEX(dpc), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_1d_dpc


  SUBROUTINE masked_swap_2d_dpc(a,b,mask)
    COMPLEX(dpc), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
    COMPLEX(dpc), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_2d_dpc


  SUBROUTINE swap_1d_dpc(a,b)
    COMPLEX(dpc), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(dpc), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_1d_dpc

  ! ------------------------------------------------------------------

END MODULE mo_percentile

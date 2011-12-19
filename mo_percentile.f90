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
  PUBLIC :: percentile      ! Returns the percent smallest value in an array (median if k=50)

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
  !         out = ksmallest(vec,k,mask=mask)
  
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
  !         None

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

  FUNCTION ksmallest_dp(arrin,k,mask)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN) :: arrin
    INTEGER(i4),                         INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(dp)                                        :: ksmallest_dp

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

    deallocate(arr)

  END FUNCTION ksmallest_dp


  FUNCTION ksmallest_sp(arrin,k,mask)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN) :: arrin
    INTEGER(i4),                         INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(sp)                                        :: ksmallest_sp

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

    if (present(mask)) then
       n = count(mask)
       allocate(arr(n))
       arr = pack(arrin,mask)
    else
       n = size(arrin)
       allocate(arr(n))
       arr = arrin
    endif

    if (n < 2) stop 'median_dp: n < 2'
    
    if (mod(n,2) == 0) then ! Even
       median_dp = 0.5_dp * (ksmallest(arr,n/2) + ksmallest(arr,n/2+1))
    else ! Odd
       median_dp = ksmallest(arr,(n+1)/2)
    endif

    deallocate(arr)

  END FUNCTION median_dp


  FUNCTION median_sp(arrin,mask)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN) :: arrin
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(sp)                                        :: median_sp

    INTEGER(i4) :: n
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

    if (n < 2) stop 'median_sp: n < 2'
    
    if (mod(n,2) == 0) then ! Even
       median_sp = 0.5_sp * (ksmallest(arr,n/2) + ksmallest(arr,n/2+1))
    else ! Odd
       median_sp = ksmallest(arr,(n+1)/2)
    endif

    deallocate(arr)

  END FUNCTION median_sp

  ! ------------------------------------------------------------------

  !     NAME
  !         percentile

  !     PURPOSE
  !         Returns the percent smallest value in an array.
  !
  !         If an optinal mask is given, values only on those locations that correspond
  !         to true values in the mask are used.

  !     CALLING SEQUENCE
  !         out = percentile(vec,k,mask=mask)
  
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

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         vec = (/ 1.,2.,3.,4.,5.,6.,7.,8.,9.,10. /)
  !         ! Returns 9
  !         out = percentile(vec,95.)
  !         ! Returns (9,8)
  !         out = percentile(vec,(/95.,80./))
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Mar 2011
  !         Modified, Stephan Thober, Dec 2011 - added 1 dimensional version

  FUNCTION percentile_0d_dp(arrin,k,mask)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN) :: arrin
    REAL(dp),                            INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(dp)                                        :: percentile_0d_dp

    INTEGER(i4) :: n, nn
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

    if (n < 2) stop 'percentile_0d_dp: n < 2'
    
    nn = floor(k/100._dp*real(n,dp),kind=i4)
    percentile_0d_dp = ksmallest(arr,nn)

    deallocate(arr)

  END FUNCTION percentile_0d_dp


  FUNCTION percentile_0d_sp(arrin,k,mask)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN) :: arrin
    REAL(sp),                            INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(sp)                                        :: percentile_0d_sp

    INTEGER(i4) :: n, nn
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

    if (n < 2) stop 'percentile_0d_sp: n < 2'
    
    nn = floor(k/100._sp*real(n,sp),kind=i4)
    percentile_0d_sp = ksmallest(arr,nn)

    deallocate(arr)

  END FUNCTION percentile_0d_sp

  function percentile_1d_dp(arrin,k,mask)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:),           INTENT(IN) :: arrin
    REAL(dp),    Dimension(:),           INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask

    REAL(dp),    Dimension(size(k))                 :: percentile_1d_dp

    INTEGER(i4) :: n, nn
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

    ! check consistency
    !if (size(k) > size(arr)) stop 'percentile_1d_dp: more Quantiles than data: size(k) > size(arr)'
    if (n < 2) stop 'percentile_1d_dp: n < 2'
    
    do n = 1, size(k)
       nn = floor(k(n)/100._dp*real(size(arr),dp),kind=i4)
       percentile_1d_dp(n) = ksmallest(arr,nn)
    end do

    deallocate(arr)

  END function percentile_1d_dp

  function percentile_1d_sp(arrin,k,mask)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:),           INTENT(IN) :: arrin
    REAL(sp),    Dimension(:),           INTENT(IN) :: k
    LOGICAL,     DIMENSION(:), OPTIONAL, INTENT(IN) :: mask

    REAL(sp),    Dimension(size(k))                 :: percentile_1d_sp

    INTEGER(i4) :: n, nn
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

    ! check consistency
    !if (size(k) > size(arr)) stop 'percentile_1d_sp: more Quantiles than data: size(k) > size(arr)'
    if (n < 2) stop 'percentile_1d_sp: n < 2'
    
    do n = 1, size(k)
       nn = floor(k(n)/100._sp*real(size(arr),sp),kind=i4)
       percentile_1d_sp(n) = ksmallest(arr,nn)
    end do

    deallocate(arr)

  END function percentile_1d_sp

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

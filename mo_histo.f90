
MODULE mo_histo

  ! This module calculates the Histogram of data (x,y) and is 
  ! part of the UFZ CHS Fortran library.

  ! Written  Juliane Mai, Feb 2012
  ! Modified 

  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: histo, histoplot           ! Histogram of data

  INTERFACE histo
     MODULE PROCEDURE histo_sp_1d, histo_dp_1d, histo_sp_2d, histo_dp_2d
  END INTERFACE histo

  INTERFACE histoplot
     MODULE PROCEDURE GnuPlot_Histo_sp_bin, GnuPlot_Histo_dp_bin, GnuPlot_Histo_sp_mean, GnuPlot_Histo_dp_mean
  END INTERFACE histoplot

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         histo

  !     PURPOSE
  !         Calculates the histogram of n data 1D points (Xi),i=1,n , 
  !         i.e. the n data points are sorted into k categories (bins). 
  !         In case of multidimensional data points (Xij),j=1,n , 
  !         only the first coordinate (X1j) is taken for determination of the category,
  !         while the other coordinates are averged per bin.  
  !
  !         With no optional arguments the number of bins k is the integer part of Sqrt(n) 
  !         where n is the number of data points. E.g. n=15 will create k=3 bins.
  !         The bin width w is calculated depending on the minimal and maximal value of 
  !         the first coordinate):
  !                     w = (maxval(X) - minval(X)) / k
  !
  !         Optionally, one can fix the number of bins k (integer i4) 
  !         while the width of the bins is determined. 
  !         Or, on can fix the binwidth w (float sp/dp) while k is determined automatically.
  !         If an optional mask is given, the histogram is taken into account only data points (Xij) 
  !         which have a true value in the mask.
  !         (Xij) can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         call hist(x, binx, bincount, width)                          ! with no optional arguments
  !         call hist(x, binx, bincount, width, bins=k,     mask=maske)  ! to fix the number of bins k
  !         call hist(x, binx, bincount, width, binwidth=w, mask=maske)  ! to fix the bin width w
  
  !     INDENT(IN)
  !         real(sp/dp) :: x(:)     1D array of x values
  !       or
  !         real(sp/dp) :: x(:,:)   2D array of x values

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(sp/dp) :: binx(:)       x-values of the histogram (center of the bin)
  !       or
  !         real(sp/dp) :: binx(:,:)     center of the bins and the averaged 2nd, 3rd... coordinates
  !         integer(i4) :: bincount(:)   number of values within each bin
  !         real(sp/dp) :: width         width of a bin

  !     INDENT(IN), OPTIONAL
  !         logical     :: mask(:)    1D-array of logical values with size(x,1).
  !                                   If present, only the data points in (Xij) corresponding 
  !                                   to the true values in mask are used.
  !         integer(i4) :: bins       Number of bins to be generated.
  !                                   If present, the number of bins will be fixed to this value.
  !                                   Otherwise, is set to integer part of sqrt(size(x)) or 
  !                                   using the optionally set binwidth.
  !         real(sp/dp) :: binwidth   Width of the bins to be generated.
  !                                   If present, the width of bins will be fixed to this value and
  !                                   the INDENT(OUT) width will be set to this value.
  !                                   Otherwise, will be determined using the number of bins.  

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.
  !         Use more than one data point: Size(x,1) > 1.
  !         If you fix the number of bins, it has to be larger than 1: bins > 1.
  !         If you fix the width of the bins, it has to be larger than 0: binwidth > 0.
  !         Either the number of bins or the width of the bins can be fixed, not both.
  !         Not all values of the first coordinate should be equal.
  !         Dimension of x and mask have to match: Size(x,1) = Size(mask). 

  !     EXAMPLE
  !         x = (/3.0, 4.0, 6.0, 7.0/)
  !
  !         binx = (/4.0, 6.0/)
  !         bincount = (/2, 2/)
  !         width = 2.0
  !
  !         x = (/ (/3.0, 4.0, 6.0, 7.0/) , (/1.0, 2.0, 5.0, 2.0/) /)
  !
  !         binx = (/ (/4.0, 6.0/) , (/1.5, 3.5/) /)
  !         bincount = (/2, 2/)
  !         width = 2.0
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         

  !     HISTORY
  !         Written,  Juliane Mai, Feb 2012
  !         Modified, 

  SUBROUTINE histo_dp_1d(x, binx, bincount, width, mask, bins, binwidth)

    IMPLICIT NONE

    REAL(DP),    DIMENSION(:),                     INTENT(IN)  :: x
    REAL(DP),    DIMENSION(:),        ALLOCATABLE, INTENT(OUT) :: binx
    INTEGER(I4), DIMENSION(:),        ALLOCATABLE, INTENT(OUT) :: bincount
    REAL(DP),                                      INTENT(OUT) :: width

    INTEGER(I4),            OPTIONAL,              INTENT(IN)  :: bins
    REAL(DP),               OPTIONAL,              INTENT(IN)  :: binwidth
    LOGICAL,  DIMENSION(:), OPTIONAL,              INTENT(IN)  :: mask

    ! Local variables
    LOGICAL, DIMENSION(size(x))              :: maske
    INTEGER(I4)                              :: n        ! number of data points
    INTEGER(I4)                              :: k        ! number of bins
    REAL(DP)                                 :: w        ! binwidth
    INTEGER(I4)                              :: i, binnr
    REAL(DP), DIMENSION(:),ALLOCATABLE       :: helpbinx
    INTEGER(I4), DIMENSION(:),ALLOCATABLE    :: helpbincount

    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error histo_dp: size(mask) /= size(x)'
       maske = mask
       n = count(maske)
    else
       maske(:) = .true.
       n = size(x)
    endif

    if (n .le. (1.0_dp+tiny(1.0_dp))) then
       stop 'Error histo_dp: size(x) must be at least 2'
    end if
    if (maxval(x(:), mask=maske)-minval(x(:), mask=maske) <= tiny(1.0_dp)) then
       stop 'Error histo_dp: all entries of x(:) are equal'
    end if

    if (present(bins) .and. present(binwidth)) then
       stop 'Error histo_dp: Either fix number of bins or binwidth, not both.' 
    end if 

    if (present(bins)) then
       if (bins .le. 1_i4) stop 'Error histo_dp: number of bins <= 1'
       k = bins
       w = (maxval(x(:), mask=maske)-minval(x(:), mask=maske)) / real(k,dp) 
    else
       if (present(binwidth)) then
          if (binwidth .le. tiny(1.0_dp)) stop 'Error histo_dp: width of bins too small'
          w = binwidth
          k = Ceiling( (maxval(x(:), mask=maske)-minval(x(:), mask=maske)) / w )
       else
          k = Floor(Sqrt(real(n,dp)))
          w = (maxval(x(:), mask=maske)-minval(x(:), mask=maske)) / real(k,dp)
       end if
    endif

    ! Histogram

    allocate (binx(k),bincount(k))

    do i=1,k
       binx(i) = minval(x(:), mask=maske) + real(2*i-1,dp)/2.0_dp*w
       bincount(i) = 0_i4
    end do

    do i=1,size(x)
       If (maske(i)) then
          binnr           = Floor((x(i)-minval(x(:), mask=maske))/w )
          ! maxval(x) has to be assigned to last bin k not k+1
          binnr           = binnr - Floor(real(binnr,dp)/real(k,dp)) + 1_i4
          bincount(binnr) = bincount(binnr) + 1_i4
       end if
    end do

    width = w

    if (count(bincount .gt. 0_i4) .lt. k) then
       print *, ' '
       print *, 'WARNING histo_dp: Empty bins have been deleted. size(binx) < bins'
       print *, ' '

       allocate(helpbinx(count(bincount .gt. 0_i4)))
       allocate(helpbincount(count(bincount .gt. 0_i4)))
       binnr=1_i4
       do i=1,k
          if (bincount(i) .gt. 0_i4) then
             helpbinx(binnr) = binx(i)
             helpbincount(binnr) = bincount(i)
             binnr=binnr+1_i4
          end if
       end do
       deallocate(binx,bincount)
       allocate(binx(size(helpbinx)))
       allocate(bincount(size(helpbinx)))
       binx=helpbinx
       bincount=helpbincount
       deallocate(helpbinx,helpbincount)
    end if

  END SUBROUTINE histo_dp_1d

  
  SUBROUTINE histo_sp_1d(x, binx, bincount, width, mask, bins, binwidth)

    IMPLICIT NONE

    REAL(SP),    DIMENSION(:),                     INTENT(IN)  :: x
    REAL(SP),    DIMENSION(:),        ALLOCATABLE, INTENT(OUT) :: binx
    INTEGER(I4), DIMENSION(:),        ALLOCATABLE, INTENT(OUT) :: bincount
    REAL(SP),                                      INTENT(OUT) :: width

    INTEGER(I4),            OPTIONAL,              INTENT(IN)  :: bins
    REAL(SP),               OPTIONAL,              INTENT(IN)  :: binwidth
    LOGICAL,  DIMENSION(:), OPTIONAL,              INTENT(IN)  :: mask

    ! Local variables
    LOGICAL, DIMENSION(size(x))              :: maske
    INTEGER(I4)                              :: n        ! number of data points
    INTEGER(I4)                              :: k        ! number of bins
    REAL(SP)                                 :: w        ! binwidth
    INTEGER(I4)                              :: i, binnr
    REAL(SP), DIMENSION(:),ALLOCATABLE       :: helpbinx
    INTEGER(I4), DIMENSION(:),ALLOCATABLE    :: helpbincount

    if (present(mask)) then
       if (size(mask) /= size(x)) stop 'Error histo_sp: size(mask) /= size(x)'
       maske = mask
       n = count(maske)
    else
       maske(:) = .true.
       n = size(x)
    endif

    if (n .le. (1.0_sp+tiny(1.0_sp))) then
       stop 'Error histo_sp: size(x) must be at least 2'
    end if
    if (maxval(x(:), mask=maske)-minval(x(:), mask=maske) <= tiny(1.0_sp)) then
       stop 'Error histo_sp: all entries of x(:) are equal'
    end if

    if (present(bins) .and. present(binwidth)) then
       stop 'Error histo_sp: Either fix number of bins or binwidth, not both.' 
    end if 

    if (present(bins)) then
       if (bins .le. 1_i4) stop 'Error histo_sp: number of bins <= 1'
       k = bins
       w = (maxval(x(:), mask=maske)-minval(x(:), mask=maske)) / real(k,sp) 
    else
       if (present(binwidth)) then
          if (binwidth .le. tiny(1.0_sp)) stop 'Error histo_sp: width of bins too small'
          w = binwidth
          k = Ceiling( (maxval(x(:), mask=maske)-minval(x(:), mask=maske)) / w )
       else
          k = Floor(Sqrt(real(n,sp)))
          w = (maxval(x(:), mask=maske)-minval(x(:), mask=maske)) / real(k,sp)
       end if
    endif

    ! Histogram

    allocate (binx(k),bincount(k))

    do i=1,k
       binx(i) = minval(x(:), mask=maske) + real(2*i-1,sp)/2.0_sp*w
       bincount(i) = 0_i4
    end do

    do i=1,size(x)
       If (maske(i)) then
          binnr           = Floor((x(i)-minval(x(:), mask=maske))/w )
          ! maxval(x) has to be assigned to last bin k not k+1
          binnr           = binnr - Floor(real(binnr,sp)/real(k,sp)) + 1_i4
          bincount(binnr) = bincount(binnr) + 1_i4
       end if
    end do

    width = w

    if (count(bincount .gt. 0_i4) .lt. k) then
       print *, ' '
       print *, 'WARNING histo_sp: Empty bins have been deleted. size(binx) < bins'
       print *, ' '

       allocate(helpbinx(count(bincount .gt. 0_i4)))
       allocate(helpbincount(count(bincount .gt. 0_i4)))
       binnr=1_i4
       do i=1,k
          if (bincount(i) .gt. 0_i4) then
             helpbinx(binnr) = binx(i)
             helpbincount(binnr) = bincount(i)
             binnr=binnr+1_i4
          end if
       end do
       deallocate(binx,bincount)
       allocate(binx(size(helpbinx)))
       allocate(bincount(size(helpbinx)))
       binx=helpbinx
       bincount=helpbincount
       deallocate(helpbinx,helpbincount)
    end if

  END SUBROUTINE histo_sp_1d

  SUBROUTINE histo_dp_2d(x, binx, bincount, width, mask, bins, binwidth)

    IMPLICIT NONE

    REAL(DP),    DIMENSION(:,:),                   INTENT(IN)  :: x
    REAL(DP),    DIMENSION(:,:),      ALLOCATABLE, INTENT(OUT) :: binx
    INTEGER(I4), DIMENSION(:),        ALLOCATABLE, INTENT(OUT) :: bincount
    REAL(DP),                                      INTENT(OUT) :: width

    INTEGER(I4),            OPTIONAL,              INTENT(IN)  :: bins
    REAL(DP),               OPTIONAL,              INTENT(IN)  :: binwidth
    LOGICAL,  DIMENSION(:), OPTIONAL,              INTENT(IN)  :: mask

    ! Local variables
    LOGICAL, DIMENSION(size(x,1))            :: maske
    INTEGER(I4)                              :: n        ! number of data points
    INTEGER(I4)                              :: k        ! number of bins
    REAL(DP)                                 :: w        ! binwidth
    INTEGER(I4)                              :: i, j, binnr
    REAL(DP), DIMENSION(:,:),ALLOCATABLE     :: helpbinx
    INTEGER(I4), DIMENSION(:),ALLOCATABLE    :: helpbincount

    if (present(mask)) then
       if (size(mask) /= size(x,1)) stop 'Error histo_dp: size(mask) /= size(x,1)'
       maske = mask
       n = count(maske)
    else
       maske(:) = .true.
       n = size(x,1)
    endif

    if (n .le. (1.0_dp+tiny(1.0_dp))) then
       stop 'Error histo_dp: size(x,1) must be at least 2'
    end if
    if (maxval(x(:,1), mask=maske)-minval(x(:,1), mask=maske) <= tiny(1.0_dp)) then
       stop 'Error histo_dp: all entries of x(:,1) are equal'
    end if

    if (present(bins) .and. present(binwidth)) then
       stop 'Error histo_dp: Either fix number of bins or binwidth, not both.' 
    end if 

    if (present(bins)) then
       if (bins .le. 1_i4) stop 'Error histo_dp: number of bins <= 1'
       k = bins
       w = (maxval(x(:,1), mask=maske)-minval(x(:,1), mask=maske)) / real(k,dp) 
    else
       if (present(binwidth)) then
          if (binwidth .le. tiny(1.0_dp)) stop 'Error histo_dp: width of bins too small'
          w = binwidth
          k = Ceiling( (maxval(x(:,1), mask=maske)-minval(x(:,1), mask=maske)) / w )
       else
          k = Floor(Sqrt(real(n,dp)))
          w = (maxval(x(:,1), mask=maske)-minval(x(:,1), mask=maske)) / real(k,dp)
       end if
    endif

    ! Histogram

    allocate (binx(k,size(x,2)),bincount(k))

    binx = 0.0_dp
    bincount = 0_i4
    do i=1,k
       binx(i,1) = minval(x(:,1), mask=maske) + real(2*i-1,dp)/2.0_dp*w
    end do

    do i=1,size(x,1)
       If (maske(i)) then
          binnr           = Floor((x(i,1)-minval(x(:,1), mask=maske))/w )
          ! maxval(x) has to be assigned to last bin k not k+1
          binnr           = binnr - Floor(real(binnr,dp)/real(k,dp)) + 1_i4
          bincount(binnr) = bincount(binnr) + 1_i4
          do j=2,size(x,2)
             binx(binnr,j)     = binx(binnr,j) + x(i,j)
          end do
       end if
    end do

    width = w

    if (count(bincount .gt. 0_i4) .lt. k) then
       print *, ' '
       print *, 'WARNING histo_dp: Empty bins have been deleted. size(binx) < bins'
       print *, ' '

       allocate(helpbinx(count(bincount .gt. 0_i4),size(x,2)))
       allocate(helpbincount(count(bincount .gt. 0_i4)))
       binnr=1_i4
       do i=1,k
          if (bincount(i) .gt. 0_i4) then
             helpbinx(binnr,1) = binx(i,1)
             do j=2,size(x,2)
                helpbinx(binnr,j) = binx(i,j)/real(bincount(i),dp)
             end do
             helpbincount(binnr) = bincount(i)
             binnr=binnr+1_i4
          !else
          !   print*, 'Delete binx = ', binx(i,1)
          end if
       end do
       deallocate(binx,bincount)
       allocate(binx(size(helpbinx,1),size(helpbinx,2)))
       allocate(bincount(size(helpbinx)))
       binx=helpbinx
       bincount=helpbincount
       deallocate(helpbinx,helpbincount)
    else
       do i=1,k
          do j=2,size(x,2)
             binx(i,j)     = binx(i,j) / real(bincount(i),dp)
          end do
       end do
    end if

  END SUBROUTINE histo_dp_2d

  
  SUBROUTINE histo_sp_2d(x, binx, bincount, width, mask, bins, binwidth)

    IMPLICIT NONE

    REAL(SP),    DIMENSION(:,:),                   INTENT(IN)  :: x
    REAL(SP),    DIMENSION(:,:),      ALLOCATABLE, INTENT(OUT) :: binx
    INTEGER(I4), DIMENSION(:),        ALLOCATABLE, INTENT(OUT) :: bincount
    REAL(SP),                                      INTENT(OUT) :: width

    INTEGER(I4),            OPTIONAL,              INTENT(IN)  :: bins
    REAL(SP),               OPTIONAL,              INTENT(IN)  :: binwidth
    LOGICAL,  DIMENSION(:), OPTIONAL,              INTENT(IN)  :: mask

    ! Local variables
    LOGICAL, DIMENSION(size(x,1))            :: maske
    INTEGER(I4)                              :: n        ! number of data points
    INTEGER(I4)                              :: k        ! number of bins
    REAL(SP)                                 :: w        ! binwidth
    INTEGER(I4)                              :: i, j, binnr
    REAL(SP), DIMENSION(:,:),ALLOCATABLE     :: helpbinx
    INTEGER(I4), DIMENSION(:),ALLOCATABLE    :: helpbincount

    if (present(mask)) then
       if (size(mask) /= size(x,1)) stop 'Error histo_sp: size(mask) /= size(x,1)'
       maske = mask
       n = count(maske)
    else
       maske(:) = .true.
       n = size(x,1)
    endif

    if (n .le. (1.0_sp+tiny(1.0_sp))) then
       stop 'Error histo_sp: size(x,1) must be at least 2'
    end if
    if (maxval(x(:,1), mask=maske)-minval(x(:,1), mask=maske) <= tiny(1.0_sp)) then
       stop 'Error histo_sp: all entries of x(:,1) are equal'
    end if

    if (present(bins) .and. present(binwidth)) then
       stop 'Error histo_sp: Either fix number of bins or binwidth, not both.' 
    end if 

    if (present(bins)) then
       if (bins .le. 1_i4) stop 'Error histo_sp: number of bins <= 1'
       k = bins
       w = (maxval(x(:,1), mask=maske)-minval(x(:,1), mask=maske)) / real(k,sp) 
    else
       if (present(binwidth)) then
          if (binwidth .le. tiny(1.0_sp)) stop 'Error histo_sp: width of bins too small'
          w = binwidth
          k = Ceiling( (maxval(x(:,1), mask=maske)-minval(x(:,1), mask=maske)) / w )
       else
          k = Floor(Sqrt(real(n,sp)))
          w = (maxval(x(:,1), mask=maske)-minval(x(:,1), mask=maske)) / real(k,sp)
       end if
    endif

    ! Histogram

    allocate (binx(k,size(x,2)),bincount(k))

    binx = 0.0_sp
    bincount = 0_i4
    do i=1,k
       binx(i,1) = minval(x(:,1), mask=maske) + real(2*i-1,sp)/2.0_sp*w
    end do

    do i=1,size(x,1)
       If (maske(i)) then
          binnr           = Floor((x(i,1)-minval(x(:,1), mask=maske))/w )
          ! maxval(x) has to be assigned to last bin k not k+1
          binnr           = binnr - Floor(real(binnr,sp)/real(k,sp)) + 1_i4
          bincount(binnr) = bincount(binnr) + 1_i4
          do j=2,size(x,2)
             binx(binnr,j)     = binx(binnr,j) + x(i,j)
          end do
       end if
    end do

    width = w

    if (count(bincount .gt. 0_i4) .lt. k) then
       print *, ' '
       print *, 'WARNING histo_sp: Empty bins have been deleted. size(binx) < bins'
       print *, ' '

       allocate(helpbinx(count(bincount .gt. 0_i4),size(x,2)))
       allocate(helpbincount(count(bincount .gt. 0_i4)))
       binnr=1_i4
       do i=1,k
          if (bincount(i) .gt. 0_i4) then
             helpbinx(binnr,1) = binx(i,1)
             do j=2,size(x,2)
                helpbinx(binnr,j) = binx(i,j)/real(bincount(i),sp)
             end do
             helpbincount(binnr) = bincount(i)
             binnr=binnr+1_i4
          !else
          !   print*, 'Delete binx = ', binx(i,1)
          end if
       end do
       deallocate(binx,bincount)
       allocate(binx(size(helpbinx,1),size(helpbinx,2)))
       allocate(bincount(size(helpbinx)))
       binx=helpbinx
       bincount=helpbincount
       deallocate(helpbinx,helpbincount)
    else
       do i=1,k
          do j=2,size(x,2)
             binx(i,j)     = binx(i,j) / real(bincount(i),sp)
          end do
       end do
    end if

  END SUBROUTINE histo_sp_2d

  ! ------------------------------------------------------------------

  !     NAME
  !         histoplot

  !     PURPOSE
  !         Plots a histogram using Gnuplot.
  !
  !         Therefore, the x-values (center of the bin) and the y-values 
  !         (either number of data points in each bin or an average y-value of each bin) 
  !         have to be given. Further, the plotter needs the specific width of the bins.
  !         Optionally, a filename can be given.

  !     CALLING SEQUENCE
  !         call histoplot(binx, bincount, width)                       ! traditional histogram without options
  !         call histoplot(binx, biny, width)                           ! histogram with averaged second coordinate
  !         call histoplot(binx, bincount, width, filename='file.eps')  ! traditional histogram with optional filename
  
  !     INDENT(IN)
  !         real(sp/dp) :: binx(:)     1D array of x values
  !         integer(i4) :: bincount(:) 1D array of data points fallen in each bin
  !       or
  !         real(sp/dp) :: biny(:)     1D array of averaged second coordinate
  !         real(sp/dp) :: width       width of a bin

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         character   :: filename   Name of the eps-file.
  !                                   Default, 'Histogram.eps'. 

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Only working on EVE cluster system.
  !          

  !     EXAMPLE
  !
  !         binx = (/4.0, 6.0/)
  !         bincount = (/2, 2/)
  !         width = 2.0
  !
  !         histoplot(binx, bincount, width)
  !         --> Histogram.eps
  !
  !         binx = (/ (/4.0, 6.0/) , (/1.5, 3.5/) /)
  !         bincount = (/2, 2/)
  !         width = 2.0
  !
  !         histoplot(binx(:,1), binx(:,2), width, filename='Histogram_averaged.eps')
  !         --> Histogram_averaged.eps
  !
  !         -> see also example in test directory

  !     LITERATURE
  !         

  !     HISTORY
  !         Written,  Juliane Mai, Feb 2012
  !         Modified, 

  SUBROUTINE GnuPlot_Histo_sp_bin(binx,bincount,width,filename)

  REAL(SP),    DIMENSION(:), INTENT(IN) :: binx
  INTEGER(I4), DIMENSION(:), INTENT(IN) :: bincount
  REAL(SP),                  INTENT(IN) :: width
  CHARACTER(*), OPTIONAL,    INTENT(IN) :: filename

  INTEGER(I4)       :: i
  CHARACTER(256)    :: gnudat, gnucmd, epsname
  CHARACTER(256)     :: str1, str2

  open (unit=10, file='data.dat', status='unknown')
  do i=1,size(binx)
     write(10,*) binx(i),bincount(i)
  end do
  close(unit=10)

  if (present(filename)) then
     epsname = filename
  else
     epsname = "Histogram.eps"
  end if

  call system('if [[ $(echo $LOADEDMODULES | grep -o gnuplot) != "gnuplot" ]]; ' // &
              'then source ~/.bashrc; module load gnuplot/4.4.4 ; fi')

  gnudat='data.dat'
  gnucmd='gnuplot.cmd'
  
   open(unit=11,file=trim(gnucmd),status='unknown')
   write(11,'(a)') 'set terminal postscript eps enhanced "Helvetica" 20 '
   write(11,'(a)') 'set output "'// trim(adjustl(epsname)) //'" '
   write(11,'(a)') 'set title "Histogram" '
   write(11,'(a)') 'set nokey '
   write (str1,*)  minval(binx)-2*(maxval(binx)-minval(binx))/size(binx)
   write (str2,*)  maxval(binx)+2*(maxval(binx)-minval(binx))/size(binx)
   write(11,'(a)') 'set xrange ['//trim(str1)//':'//trim(str2)//']'
   write (str1,*)  0
   write (str2,*)  maxval(bincount)+1
   write(11,'(a)') 'set yrange ['//trim(str1)//':'//trim(str2)//']'
   write(11,'(a)') 'set style fill solid 0.25 border'
   write (str1,*)  width   
   write(11,'(a)') 'set boxwidth '//trim(str1)
   write(11,'(a)') 'plot "'//trim(gnudat)//'" with boxes'
   write(11,'(a)') 'quit'
   close(11)

   call system(trim('gnuplot '//trim(gnucmd)))

   call system('rm -f ' // trim(gnudat) )
   call system('rm -f ' // trim(gnucmd) )

  END SUBROUTINE GnuPlot_Histo_sp_bin

  SUBROUTINE GnuPlot_Histo_dp_bin(binx,bincount,width,filename)

  REAL(DP),    DIMENSION(:), INTENT(IN) :: binx
  INTEGER(I4), DIMENSION(:), INTENT(IN) :: bincount
  REAL(DP),                  INTENT(IN) :: width
  CHARACTER(*), OPTIONAL,    INTENT(IN) :: filename

  INTEGER(I4)       :: i
  CHARACTER(256)    :: gnudat, gnucmd, epsname
  CHARACTER(256)     :: str1, str2

  open (unit=10, file='data.dat', status='unknown')
  do i=1,size(binx)
     write(10,*) binx(i),bincount(i)
  end do
  close(unit=10)

  if (present(filename)) then
     epsname = filename
  else
     epsname = "Histogram.eps"
  end if

  call system('if [[ $(echo $LOADEDMODULES | grep -o gnuplot) != "gnuplot" ]]; ' // &
              'then source ~/.bashrc; module load gnuplot/4.4.4 ; fi')

  gnudat='data.dat'
  gnucmd='gnuplot.cmd'
  
   open(unit=11,file=trim(gnucmd),status='unknown')
   write(11,'(a)') 'set terminal postscript eps enhanced "Helvetica" 20 '
   write(11,'(a)') 'set output "'// trim(adjustl(epsname)) //'" '
   write(11,'(a)') 'set title "Histogram" '
   write(11,'(a)') 'set nokey '
   write (str1,*) minval(binx)-2*(maxval(binx)-minval(binx))/size(binx)
   write (str2,*) maxval(binx)+2*(maxval(binx)-minval(binx))/size(binx)
   write(11,'(a)') 'set xrange ['//trim(str1)//':'//trim(str2)//']'
   write (str1,*) 0
   write (str2,*) maxval(bincount)+1
   write(11,'(a)') 'set yrange ['//trim(str1)//':'//trim(str2)//']'
   write(11,'(a)') 'set style fill solid 0.25 border'
   write (str1,*)  width   
   write(11,'(a)') 'set boxwidth '//trim(str1)
   write(11,'(a)') 'plot "'//trim(gnudat)//'" with boxes'
   write(11,'(a)') 'quit'
   close(11)

   call system(trim('gnuplot '//trim(gnucmd)))

   call system('rm -f ' // trim(gnudat) )
   call system('rm -f ' // trim(gnucmd) )

  END SUBROUTINE GnuPlot_Histo_dp_bin

  SUBROUTINE GnuPlot_Histo_sp_mean(binx,biny,width,filename)

  REAL(SP),    DIMENSION(:),                   INTENT(IN) :: binx
  REAL(SP),    DIMENSION(:),                   INTENT(IN) :: biny      ! Mean of y-Values in bin x
  REAL(SP),                                    INTENT(IN) :: width
  CHARACTER(*),                   OPTIONAL,    INTENT(IN) :: filename

  INTEGER(I4)       :: i
  CHARACTER(256)    :: gnudat, gnucmd, epsname
  CHARACTER(256)    :: str1, str2

  open (unit=10, file='data.dat', status='unknown')
  do i=1,size(binx)
     write(10,*) binx(i),biny(i)
  end do
  close(unit=10)

  if (present(filename)) then
     epsname = filename
  else
     epsname = "Histogram.eps"
  end if

  call system('if [[ $(echo $LOADEDMODULES | grep -o gnuplot) != "gnuplot" ]]; ' // &
              'then source ~/.bashrc ; module load gnuplot/4.4.4 ; fi')
  gnudat='data.dat'
  gnucmd='gnuplot.cmd'
  
   open(unit=11,file=trim(gnucmd),status='unknown')
   write(11,'(a)') 'set terminal postscript eps enhanced "Helvetica" 20 '
   write(11,'(a)') 'set output "'// trim(adjustl(epsname)) //'" '
   write(11,'(a)') 'set title "Histogram" '
   write(11,'(a)') 'set nokey '
   write (str1,*)  minval(binx)-2*(maxval(binx)-minval(binx))/size(binx)
   write (str2,*)  maxval(binx)+2*(maxval(binx)-minval(binx))/size(binx)
   write(11,'(a)') 'set xrange ['//trim(str1)//':'//trim(str2)//']'
   write (str1,*)  minval(biny)-2*(maxval(biny)-minval(biny))/size(biny)
   write (str2,*)  maxval(biny)+2*(maxval(biny)-minval(biny))/size(biny)
   write(11,'(a)') 'set yrange ['//trim(str1)//':'//trim(str2)//']'
   write(11,'(a)') 'set style fill solid 0.25 border'
   write (str1,*)  width   
   write(11,'(a)') 'set boxwidth '//trim(str1)
   write(11,'(a)') 'plot "'//trim(gnudat)//'" with boxes'
   write(11,'(a)') 'quit'
   close(11)

   call system(trim('gnuplot '//trim(gnucmd)))

   call system('rm -f ' // trim(gnudat) )
   call system('rm -f ' // trim(gnucmd) )

  END SUBROUTINE GnuPlot_Histo_sp_mean

  SUBROUTINE GnuPlot_Histo_dp_mean(binx,biny,width,filename)

  REAL(DP),    DIMENSION(:),                   INTENT(IN) :: binx
  REAL(DP),    DIMENSION(:),                   INTENT(IN) :: biny      ! Mean of y-Values in bin x
  REAL(DP),                                    INTENT(IN) :: width
  CHARACTER(*),                   OPTIONAL,    INTENT(IN) :: filename

  INTEGER(I4)       :: i
  CHARACTER(256)    :: gnudat, gnucmd, epsname
  CHARACTER(256)    :: str1, str2

  open (unit=10, file='data.dat', status='unknown')
  do i=1,size(binx,1)
     write(10,*) binx(i),biny(i)
  end do
  close(unit=10)

  if (present(filename)) then
     epsname = filename
  else
     epsname = "Histogram.eps"
  end if

  call system('if [[ $(echo $LOADEDMODULES | grep -o gnuplot) != "gnuplot" ]];' // &
              'then source ~/.bashrc; module load gnuplot/4.4.4 ; fi')

  gnudat='data.dat'
  gnucmd='gnuplot.cmd'
  
   open(unit=11,file=trim(gnucmd),status='unknown')
   write(11,'(a)') 'set terminal postscript eps enhanced "Helvetica" 20 '
   write(11,'(a)') 'set output "'// trim(adjustl(epsname)) //'" '
   write(11,'(a)') 'set title "Histogram" '
   write(11,'(a)') 'set nokey '
   write (str1,*)  minval(binx)-2*(maxval(binx)-minval(binx))/size(binx)
   write (str2,*)  maxval(binx)+2*(maxval(binx)-minval(binx))/size(binx)
   write(11,'(a)') 'set xrange ['//trim(str1)//':'//trim(str2)//']'
   write (str1,*)  minval(biny)-2*(maxval(biny)-minval(biny))/size(biny)
   write (str2,*)  maxval(biny)+2*(maxval(biny)-minval(biny))/size(biny)
   write(11,'(a)') 'set yrange ['//trim(str1)//':'//trim(str2)//']'
   write(11,'(a)') 'set style fill solid 0.25 border'
   write (str1,*)  width   
   write(11,'(a)') 'set boxwidth '//trim(str1)
   write(11,'(a)') 'plot "'//trim(gnudat)//'" with boxes'
   write(11,'(a)') 'quit'
   close(11)

   call system(trim('gnuplot '//trim(gnucmd)))

   call system('rm -f ' // trim(gnudat) )
   call system('rm -f ' // trim(gnucmd) )

  END SUBROUTINE GnuPlot_Histo_dp_mean

  ! ------------------------------------------------------------------

END MODULE mo_histo

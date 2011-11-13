MODULE mo_template

  ! This module is a template for the UFZ CHS Fortran library.

  ! The module provides an example procedures and demonstrates the agreed coding standard:
  !     - Use mo_kind, only: i4, sp, dp, lgt, etc.
  !     - Always make procedures for single and double precision (sp, dp), i.e. include interface to module procedures.
  !     - Include optional mask argument, if possible.
  !     - Make the module private by default.
  !     - Make routines available explicitly, i.e. public.
  !     - Break lines at column 130 at most.
  !     - Do not use tabs in files.
  !     - Give 1-line descriptions after the public definition.
  !     - Documentation:
  !       * document before the individual routines;
  !       * do one documentation per interface, i.e. no separate docu for sp and dp;
  !       * follow the documentation structure before the function mean_sp in mo_template.f90;
  !       * break comment lines at column 130 at most as well.
  !     - Sort routines alphabetically in the file and in the public definitions.
  !     - The modules should be tested with at least two different compilers (of different vendors)

  ! Written  Matthias Cuntz, Nov 2011
  ! Modified Matthias Cuntz, Nov 2011 - add private
  !                          Nov 2011 - add public

  ! Always use the number precisions of mo_kind
  USE mo_kind, ONLY: i4, sp, dp, lgt

  ! Of course
  IMPLICIT NONE

  ! Make everything private by default
  PRIVATE

  ! Explicitly make public only the routines, parameters, etc. that shall be provided
  ! Sort alphabetically and give 1-line descriptions
  PUBLIC :: mean            ! 1st moment of an array, i.e. the mean
  PUBLIC :: PI_dp           ! Constant Pi in double precision
  PUBLIC :: PI_sp           ! Constant Pi in single precision

  ! Interfaces for single and double precision routines; sort alphabetically
  INTERFACE mean
     MODULE PROCEDURE mean_sp, mean_dp
  END INTERFACE mean

  ! Public parameters
  REAL(dp), PARAMETER :: PI_dp = 3.141592653589793238462643383279502884197_dp
  REAL(sp), PARAMETER :: PI_sp = 3.141592653589793238462643383279502884197_sp

  ! Private global parameters (not used, only for demonstration)
  INTEGER(i4), PARAMETER :: iTest=1

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         mean

  !     PURPOSE
  !         Calculates the average value of a vector, i.e. the first moment of a series of numbers:
  !             mean = sum(x)/n
  !         If an optinal mask is given, the mean is only over those locations that correspond to true values in the mask:
  !             mean = sum(x, mask)/count(mask)
  !         x can be single or double precision. The mean will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = mean(vec, mask=mask)
  
  !     INDENT(IN)
  !         real(sp/dp) :: vec(:)     1D-array with input numbers

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         real(sp/dp) :: mean       average of all elements in dat

  !     INDENT(IN), OPTIONAL
  !         logical(lpt) :: mask(:)   1D-array of logical values with size(vec).
  !                                   If present, only those locations in vec corresponding to the true values in mask are used.

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         vec = (/ 1., 2, 3., -999., 5., 6. /)
  !         m   = mean(vec, mask=(vec >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !         Sokal RR & Rohlf FJ - Biometry: the principle and practice of statistics in biological research,
  !             Freeman & Co., ISBN 0-7167-2411-1
  !         Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !             The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !             Cambridge University Press, UK, 1996

  !     HISTORY
  !         Written,  Matthias Cuntz, Nov 2011
  !         Modified, Matthias Cuntz, Nov 2011 - include mask
  !                   Matthias Cuntz, Nov 2011 - test size(mask) == size(dat)

  FUNCTION mean_dp(dat, mask)

    IMPLICIT NONE

    REAL(dp),     DIMENSION(:),           INTENT(IN)  :: dat
    LOGICAL(lgt), DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(dp)                                          :: mean_dp

    REAL(dp) :: n

    LOGICAL(lgt), DIMENSION(size(dat)) :: maske

    if (present(mask)) then
       if (size(mask) /= size(dat)) stop 'Error mean_dp: size(mask) /= size(dat)'
       maske = mask
       n = real(count(maske),dp)
    else
       maske(:) = .true.
       n = real(size(dat),dp)
    endif
    if (n <= (1.0_dp+tiny(1.0_dp))) stop 'mean_dp: n must be at least 2'

    ! Mean
    mean_dp  = sum(dat(:), mask=maske)/n

  END FUNCTION mean_dp

  
  FUNCTION mean_sp(dat, mask)

    IMPLICIT NONE

    REAL(sp),     DIMENSION(:),           INTENT(IN)  :: dat
    LOGICAL(lgt), DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(sp)                                          :: mean_sp

    REAL(sp) :: n

    LOGICAL(lgt), DIMENSION(size(dat)) :: maske

    if (present(mask)) then
       if (size(mask) /= size(dat)) stop 'Error mean_sp: size(mask) /= size(dat)'
       maske = mask
       n = real(count(maske),sp)
    else
       maske(:) = .true.
       n = real(size(dat),sp)
    endif
    if (n <= (1.0_sp+tiny(1.0_sp))) stop 'mean_sp: n must be at least 2'

    ! Mean
    mean_sp  = sum(dat(:), mask=maske)/n

  END FUNCTION mean_sp

  ! ------------------------------------------------------------------

END MODULE mo_template

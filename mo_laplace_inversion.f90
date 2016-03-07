MODULE mo_laplace_inversion

  ! This module serves numerical laplace-inversion

  ! ?
  !

  ! Written Sebastian Mueller, June 2014

  ! License
  ! -------
  ! This file is part of the JAMS Fortran library.

  ! The JAMS Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The JAMS Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the JAMS Fortran library (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2014 Sebastian Mueller

  use mo_kind, only: i4, i8, SP, DP
  use mo_combinatorics,   only: binomcoeffi
  use mo_functions,       only: factorial

  Implicit NONE

  PUBLIC :: NLInvSteh            ! Laplace back-transformed function by Stehfest algorithm

  ! ------------------------------------------------------------------

  !     NAME
  !         NLInvSteh

  !     PURPOSE
  !         Calculates numerically the laplace-inversion of a given function for given time-values (and optionaly space-values)
  !         that are given in an array
  !     CALLING SEQUENCE
  !         out = NLInvStehvec(func, para, [r,] t, n)

  !     INTENT(IN)
  !>        REAL(sp/dp)                 :: func         function to be inverted (f([r,]s,para))
  !>        REAL(sp/dp), DIMENSION(:)   :: para         Parameters for the given function
  !>       [REAL(sp/dp), DIMENSION(:)   :: r            input radius]
  !>        REAL(sp/dp), DIMENSION(:)   :: t            input time
  !>        integer(i4/i8)              :: n            boundary for the stehfest-algorithm (recommended: 6/8/10/12)

  !     INTENT(INOUT)
  !>        None

  !     INTENT(OUT)
  !>        None

  !     INTENT(IN), OPTIONAL
  !>        None

  !     INTENT(INOUT), OPTIONAL
  !>        None

  !     INTENT(OUT), OPTIONAL
  !>        None

  !     RETURNS
  !>        real(dp)                    :: NLInvSteh

  !     RESTRICTIONS
  !         1) Stehfest-boundary n must be positiv

  !     EXAMPLE
  !         func(s) = 1/(s^2)

  !         NLInvSteh(func,t,12) ~ t

  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written, Sebastian Mueller, June 2014

  INTERFACE NLInvSteh
     MODULE PROCEDURE NLInvSteh_dp,     NLInvSteh_sp,       NLInvStehspace_dp,      NLInvStehspace_sp,&
          NLInvStehvec_dp,  NLInvStehvec_sp,    NLInvStehvecspace_dp,   NLInvStehvecspace_sp
  END INTERFACE NLInvSteh

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  ! c_n=csteh(n, m) computed from the Stehfestalg. m=N/2 n=n k=k (the values are all written in an array)
  function csteh_dp(m)

    implicit none

    integer(i4), intent(in)  :: m
    real(dp), dimension(2*m) :: csteh_dp

    ! local variables
    integer(i4)              :: k, n

    !Check if the given boundary is not to big
    if (m>=11_i4)            stop 'The Stehfest-boundary must be less than 22 for double-precision.'

    csteh_dp            =   0.0_dp

    do n=1_i4, 2*m
       do k=FLOOR((real(n,dp)+1.0_dp)/2.0_dp), min(n,m)
          csteh_dp(n)  =   csteh_dp(n) +   real(k,dp)**real(m+1_i4,dp)*real(binomcoeffi(int(2_i4*k,i8),int(k,i8)),dp)/&
               real(factorial(int(m-k,i8))*factorial(int(n-k,i8))*factorial(int(2_i4*k-n,i8)),dp)

!          csteh_dp(n)  =   csteh_dp(n) +   real(k**(m+1_i4)*binomcoeffi(2_i4*k,k),dp)/&
!               real(factorial(m-k)*factorial(n-k)*factorial(2_i4*k-n),dp)

       end do
       csteh_dp(n)     =   (-1.0_dp)**(n+m)*csteh_dp(n)
!write(*,*) "csteh_dp(",n,")= ", csteh_dp(n)
    end do

  end function csteh_dp

  function csteh_sp(m)

    implicit none

    integer(i4), intent(in)  :: m
    real(sp), dimension(2*m) :: csteh_sp

    ! local variables
    integer(i4)              :: k, n

    !Check if the given boundary is not to big
    if (m>=7_i4)            stop 'The Stehfest-boundary must be less than 14 for single-precision.'

    csteh_sp            =   0.0_sp

    do n=1_i4, 2*m
       do k=FLOOR((real(n,sp)+1.0_sp)/2.0_sp), min(n,m)
          csteh_sp(n)  =   csteh_sp(n) +   real(k,sp)**real(m+1_i4,sp)*real(binomcoeffi(int(2_i4*k,i8),int(k,i8)),sp)/&
               real(factorial(int(m-k,i8))*factorial(int(n-k,i8))*factorial(int(2_i4*k-n,i8)),sp)
       end do
       csteh_sp(n)     =   (-1.0_sp)**(n+m)*csteh_sp(n)
!write(*,*) "csteh_sp(",n,")= ", csteh_sp(n)
    end do

  end function csteh_sp

  !Stehfestalg for only-time functions with array-input

  function NLInvStehvec_dp( func, para, t, n)

    !Interface for the given function
    INTERFACE
       function func(s, para)
         use mo_kind, only: dp
         implicit none
         real(dp), dimension(:), intent(in) :: s
         real(dp), dimension(:), intent(in) :: para
         real(dp), dimension(size(s))       :: func
       end function func
    END INTERFACE

    real(dp),    dimension(:), intent(in)   :: para
    real(dp),    dimension(:), intent(in)   :: t
    integer(i4),               intent(in)   :: n

    real(dp),    dimension(size(t))         :: NLInvStehvec_dp

    !intern variables
    integer(i4)                                                     :: m,i
    real(dp),   dimension(size(t))                                  :: a
    real(dp),   dimension(2*Ceiling(real(n,dp)/2.0_dp))             :: csteh
    real(dp),   dimension(size(t), 2*Ceiling(real(n,dp)/2.0_dp))    :: funcpoints

    !Check if the given boundary is a positiv number
    if (n<=0_i4)            stop 'The Stehfest-boundary n must be positiv.'
    !Check if there are time-values
    if (size(t)==0_i4)      stop 'There are no given time-values.'

    !If n is odd take the next bigger even number N=2m
    m                       =   Ceiling(real(n,dp)/2.0_dp)

    !Define the coefficients c_n
    csteh                   =   csteh_dp(m)

    !Evaluate all necessary function points to call the function only one time
    !set all necessary s-values into a matrix
    do i=1_i4, 2*m
       funcpoints(:,i)     =   real(i,dp)/t*log(2.0_dp)
    end do

    !get the function-values for all needed s-values
    funcpoints              =   reshape(func(pack(funcpoints,.true.), para), (/size(t),2*m/))

    !Calculate the Lapalce-inversion
    a                       =   log(2.0_dp)/t

    do i = 1, size(t)
       NLInvStehvec_dp(i)  =   dot_product(csteh, funcpoints(i,:))
    end do

    NLInvStehvec_dp         =   a*NLInvStehvec_dp

  end function NLInvStehvec_dp

  !-----------------------------------------------------------------

  function NLInvStehvec_sp( func, para, t, n)

    !Interface for the given function
    INTERFACE
       function func(s, para)
         use mo_kind, only: sp
         implicit none
         real(sp), dimension(:), intent(in) :: s
         real(sp), dimension(:), intent(in) :: para
         real(sp), dimension(size(s))       :: func
       end function func
    END INTERFACE

    real(sp),    dimension(:), intent(in)   :: para
    real(sp),    dimension(:), intent(in)   :: t
    integer(i4),               intent(in)   :: n

    real(sp),    dimension(size(t))         :: NLInvStehvec_sp

    !intern variables
    integer(i4)                                                     :: m,i
    real(sp),   dimension(size(t))                                  :: a
    real(sp),   dimension(2*Ceiling(real(n,sp)/2.0_sp))             :: csteh
    real(sp),   dimension(size(t),2*Ceiling(real(n,sp)/2.0_sp))     :: funcpoints

    !Check if the given boundary is a positiv number
    if (n<=0_i4)            stop 'The Stehfest-boundary n must be positiv.'

    !Check if there are time-values
    if (size(t)==0_i4)      stop 'There are no given time-values.'

    !If n is odd take the next bigger even number N=2m
    m                       =   Ceiling(real(n,sp)/2.0_sp)

    !Define the coefficients c_n
    csteh                   =   csteh_sp(m)

    !Evaluate all necessary function points to call the function only one time
    !set all necessary s-values into a matrix
    do i=1_i4, 2*m
       funcpoints(:,i)     =   real(i,sp)/t*log(2.0_sp)
    end do

    !get the function-values for all needed s-values
    funcpoints              =   reshape(func(pack(funcpoints,.true.), para), (/size(t),2*m/))

    !Calculate the Lapalce-inversion
    a                       =   log(2.0_sp)/t

    do i = 1, size(t)
       NLInvStehvec_sp(i)  =   dot_product(csteh, funcpoints(i,:))
    end do

    NLInvStehvec_sp         =   a*NLInvStehvec_sp

  end function NLInvStehvec_sp

  !------------------------------------------------------------------
  !Stehfestalg for space and time functions with array-input

  function NLInvStehvecspace_dp( func, para, r, t, n)

    !Interface for the given function with additional radial-values r
    INTERFACE
       function func(r, s, para)
         use mo_kind, only: dp
         implicit none
         real(dp), dimension(:), intent(in)     :: r
         real(dp), dimension(:), intent(in)     :: s
         real(dp), dimension(:), intent(in)     :: para
         real(dp), dimension(size(r),size(s))   :: func
       end function func
    END INTERFACE

    real(dp),    dimension(:), intent(in)       :: para
    real(dp),    dimension(:), intent(in)       :: r
    real(dp),    dimension(:), intent(in)       :: t
    integer(i4),               intent(in)       :: n


    real(dp),    dimension(size(r),size(t))     :: NLInvStehvecspace_dp

    !intern variables
    integer(i4)                                                             :: m,i,j
    real(dp),   dimension(size(t))                                          :: a
    real(dp),   dimension(2*Ceiling(real(n,dp)/2.0_dp))                     :: csteh
    real(dp),   dimension(size(t), 2*Ceiling(real(n,dp)/2.0_dp))            :: timepoints
    real(dp),   dimension(size(r), size(t), 2*Ceiling(real(n,dp)/2.0_dp))   :: funcpoints

    !Check if the given boundary is a positiv number
    if (n<=0_i4)            stop 'The Stehfest-boundary n must be positiv.'
    !Check if there are time-values and if they are valid
    if (size(t)==0_i4)      stop 'There are no given time-values.'
    !Check if there are space-values and if they are valid
    if (size(r)==0_i4)      stop 'There are no given space-values.'

    !If n is odd take the next bigger even number N=2m
    m                       =   Ceiling(real(n,dp)/2.0_dp)

    !Define the coefficients c_n
    csteh                   =   csteh_dp(m)

    !Evaluate all necessary function points to call the function only one time
    !set all necessary s-values into a matrix
    do i=1_i4, 2*m
       timepoints(:,i)     =   real(i,dp)/t*log(2.0_dp)
    end do

    !get the function-values for all needed s-values
    funcpoints              =   reshape(func(r,pack(timepoints,.true.), para), (/size(r),size(t),2*m/))

    !Calculate the Lapalce-inversion
    a                       =   log(2.0_dp)/t

    do i = 1, size(r)
       do j = 1, size(t)
          NLInvStehvecspace_dp(i,j)   =   dot_product(csteh, funcpoints(i,j,:))
       end do
       NLInvStehvecspace_dp(i,:)       =   a*NLInvStehvecspace_dp(i,:)
    end do

  end function NLInvStehvecspace_dp

  function NLInvStehvecspace_sp( func, para, r, t, n)

    !Interface for the given function with additional radial-values r
    INTERFACE
       function func(r, s, para)
         use mo_kind, only: sp
         implicit none
         real(sp), dimension(:), intent(in)     :: r
         real(sp), dimension(:), intent(in)     :: s
         real(sp), dimension(:), intent(in)     :: para
         real(sp), dimension(size(r),size(s))   :: func
       end function func
    END INTERFACE

    real(sp),    dimension(:), intent(in)       :: para
    real(sp),    dimension(:), intent(in)       :: r
    real(sp),    dimension(:), intent(in)       :: t
    integer(i4),               intent(in)       :: n

    real(sp),    dimension(size(r),size(t))     :: NLInvStehvecspace_sp

    !intern variables
    integer(i4)                                                             :: m,i,j
    real(sp),   dimension(size(t))                                          :: a
    real(sp),   dimension(2*Ceiling(real(n,sp)/2.0_sp))                     :: csteh
    real(sp),   dimension(size(t), 2*Ceiling(real(n,sp)/2.0_sp))            :: timepoints
    real(sp),   dimension(size(r), size(t), 2*Ceiling(real(n,sp)/2.0_sp))   :: funcpoints

    !Check if the given boundary is a positiv number
    if (n<=0_i4)            stop 'The Stehfest-boundary n must be positiv.'
    !Check if there are time-values and if they are valid
    if (size(t)==0_i4)      stop 'There are no given time-values.'
    !Check if there are space-values and if they are valid
    if (size(r)==0_i4)      stop 'There are no given space-values.'

    !If n is odd take the next bigger even number N=2m
    m                       =   Ceiling(real(n,sp)/2.0_sp)

    !Define the coefficients c_n
    csteh                   =   csteh_sp(m)

    !Evaluate all necessary time points to call the function only ONE time
    !set all necessary s-values into a matrix
    do i=1_i4, 2*m
       timepoints(:,i)     =   real(i,sp)/t*log(2.0_sp)
    end do

    !get the function-values for all needed s-values
    funcpoints              =   reshape(func(r,pack(timepoints,.true.), para), (/size(r),size(t),2*m/))

    !Calculate the Lapalce-inversion
    a                       =   log(2.0_sp)/t

    do i = 1, size(r)
       do j = 1, size(t)
          NLInvStehvecspace_sp(i,j)   =   dot_product(csteh, funcpoints(i,j,:))
       end do
       NLInvStehvecspace_sp(i,:)       =   a*NLInvStehvecspace_sp(i,:)
    end do

  end function NLInvStehvecspace_sp

  !stehfest for single input value

  function NLInvSteh_dp( func, para, t, n)

    !Interface for the given function
    INTERFACE
       function func(s, para)
         use mo_kind, only: dp
         implicit none
         real(dp),               intent(in) :: s
         real(dp), dimension(:), intent(in) :: para
         real(dp)                           :: func
       end function func
    END INTERFACE

    real(dp),    dimension(:), intent(in)   :: para
    real(dp),                  intent(in)   :: t
    integer(i4),               intent(in)   :: n
    real(dp)                                :: NLInvSteh_dp

    !intern variables
    integer(i4)                                         :: m,i
    real(dp),   dimension(2*Ceiling(real(n,dp)/2.0_dp)) :: csteh, funcpoints

    !Check if the given boundary is a positiv number
    if (n<=0_i4)            stop 'The Stehfest-boundary n must be positiv.'

    !If n is odd take the next bigger even number N=2m
    m                       =   Ceiling(real(n,dp)/2.0_dp)

    !Define the coefficients c_n
    csteh                   =   csteh_dp(m)

    !Evaluate all necessary function points
    !get the function-values for all needed s-values
    do i=1_i4, 2*m
       funcpoints(i)       =   func(real(i,dp)/t*log(2.0_dp), para)
    end do

    !Calculate the Lapalce-inversion
    NLInvSteh_dp            =   (log(2.0_dp)/t)*dot_product(csteh, funcpoints)

  end function NLInvSteh_dp

  !------------------------------------------------------------------
  !Stehfestalg for only-time functions with single-input

  function NLInvSteh_sp( func, para, t, n)

    !Interface for the given function
    INTERFACE
       function func(s, para)
         use mo_kind, only: sp
         implicit none
         real(sp),               intent(in) :: s
         real(sp), dimension(:), intent(in) :: para
         real(sp)                           :: func
       end function func
    END INTERFACE

    real(sp),    dimension(:), intent(in)   :: para
    real(sp),                  intent(in)   :: t
    integer(i4),               intent(in)   :: n
    real(sp)                                :: NLInvSteh_sp

    !intern variables
    integer(i4)                                         :: m,i
    real(sp),   dimension(2*Ceiling(real(n,sp)/2.0_sp)) :: csteh, funcpoints

    !Check if the given boundary is a positiv number
    if (n<=0_i4)            stop 'The Stehfest-boundary n must be positiv.'

    !If n is odd take the next bigger even number N=2m
    m                       =   Ceiling(real(n,sp)/2.0_sp)

    !Define the coefficients c_n
    csteh                   =   csteh_sp(m)

    !Evaluate all necessary function points
    !get the function-values for all needed s-values
    do i=1_i4, 2*m
       funcpoints(i)       =   func(real(i,sp)/t*log(2.0_sp), para)
    end do

    !Calculate the Lapalce-inversion
    NLInvSteh_sp        =   (log(2.0_sp)/t)*dot_product(csteh, funcpoints)

  end function NLInvSteh_sp

  !------------------------------------------------------------------
  !Stehfestalg for space and time functions with single-input

  function NLInvStehspace_dp( func, para, r, t, n)

    !Interface for the given function with additional radial-values r
    INTERFACE
       function func(r, s, para)
         use mo_kind, only: dp
         implicit none
         real(dp),               intent(in)     :: r
         real(dp),               intent(in)     :: s
         real(dp), dimension(:), intent(in)     :: para
         real(dp)                               :: func
       end function func
    END INTERFACE

    real(dp),    dimension(:), intent(in)       :: para
    real(dp),                  intent(in)       :: r,t
    integer(i4),               intent(in)       :: n

    real(dp)                                    :: NLInvStehspace_dp

    !intern variables
    integer(i4)                                         :: m,i
    real(dp),   dimension(2*Ceiling(real(n,dp)/2.0_dp)) :: csteh, funcpoints

    !Check if the given boundary is a positiv number
    if (n<=0_i4)            stop 'The Stehfest-boundary n must be positiv.'

    !If n is odd take the next bigger even number N=2m
    m                       =   Ceiling(real(n,dp)/2.0_dp)

    !Define the coefficients c_n
    csteh                   =   csteh_dp(m)

    !get the function-values for all needed s-values
    do i=1_i4, 2*m
       funcpoints(i)       =   func(r, real(i,dp)/t*log(2.0_dp), para)
    end do

    !Calculate the Lapalce-inversion
    NLInvStehspace_dp       =   (log(2.0_dp)/t)*dot_product(csteh, funcpoints)

  end function NLInvStehspace_dp


  function NLInvStehspace_sp( func, para, r, t, n)

    !Interface for the given function with additional radial-values r
    INTERFACE
       function func(r, s, para)
         use mo_kind, only: sp
         implicit none
         real(sp),               intent(in)     :: r
         real(sp),               intent(in)     :: s
         real(sp), dimension(:), intent(in)     :: para
         real(sp)                               :: func
       end function func
    END INTERFACE

    real(sp),    dimension(:), intent(in)       :: para
    real(sp),                  intent(in)       :: r,t
    integer(i4),               intent(in)       :: n

    real(sp)                                    :: NLInvStehspace_sp

    !intern variables
    integer(i4)                                         :: m,i
    real(sp),   dimension(2*Ceiling(real(n,sp)/2.0_sp)) :: csteh, funcpoints

    !Check if the given boundary is a positiv number
    if (n<=0_i4)            stop 'The Stehfest-boundary n must be positiv.'

    !If n is odd take the next bigger even number N=2m
    m                       =   Ceiling(real(n,sp)/2.0_sp)

    !Define the coefficients c_n
    csteh                   =   csteh_sp(m)

    !get the function-values for all needed s-values
    do i=1_i4, 2*m
       funcpoints(i)       =   func(r, real(i,sp)/t*log(2.0_sp), para)
    end do

    !Calculate the Lapalce-inversion
    NLInvStehspace_sp       =   (log(2.0_sp)/t)*dot_product(csteh, funcpoints)

  end function NLInvStehspace_sp

end module mo_laplace_inversion

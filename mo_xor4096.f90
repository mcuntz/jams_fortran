module mo_xor4096

  ! -----------------------------------------------------------------------------
  ! The original version of this source code  (without multiple streams, optional
  ! arguments and gaussian distributed RN) is under GNU General Public Licence
  !      xorgens.c
  ! Copyright (C) 2004 R. P. Brent.                                       
  ! This program is free software; you can redistribute it and/or         
  ! modify it under the terms of the GNU General Public License,       
  ! version 2, June 1991, as published by the Free Software Foundation.   
  ! For details see http://www.gnu.org/copyleft/gpl.html .                
  !
  ! Author: Richard P. Brent (random@rpbrent.co.uk)
  ! -----------------------------------------------------------------------------

  ! Written  Juliane Mai, Nov 2011

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

  ! Copyright 2011 Juliane Mai

  use mo_kind, only: i4, i8, sp, dp

  Implicit NONE

  PRIVATE

  PUBLIC :: get_timeseed    ! Returns a seed dependend on time
  PUBLIC :: xor4096         ! Generates uniform distributed random number
  PUBLIC :: xor4096g        ! Generates gaussian distributed random number

  ! Interfaces for single and double precision routines
  INTERFACE get_timeseed
     MODULE PROCEDURE   get_timeseed_i4_0d, get_timeseed_i4_1d, &
                        get_timeseed_i8_0d, get_timeseed_i8_1d
  END INTERFACE get_timeseed

  INTERFACE xor4096
     MODULE PROCEDURE   xor4096s_0d, xor4096s_1d, xor4096f_0d, xor4096f_1d, &
                        xor4096l_0d, xor4096l_1d, xor4096d_0d, xor4096d_1d
  END INTERFACE xor4096

  INTERFACE xor4096g
     MODULE PROCEDURE xor4096gf_0d, xor4096gf_1d, xor4096gd_0d, xor4096gd_1d
  END INTERFACE xor4096g

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         get_timeseed

  !     PURPOSE
  !         Function which returns a scalar or a vector of integers
  !         dependent on time.
  !         This returned values can be used a seed for initializing
  !         random number generators like xor4096.
  !
  !         If the return variable is a vector, only the first entry
  !         depends on time while the others are just the entry before
  !         increased by 1000.

  !     CALLING SEQUENCE
  !         call get_timeseed(seed)
  
  !     INDENT(IN)
  !         None

  !     INDENT(INOUT)
  !         integer(i4/i8)                     :: seed
  !           OR
  !         integer(i4/i8), dimension(:)       :: seed

  !     INDENT(OUT)
  !         None                                                       

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         integer(i4) :: seed
  !         call get_timeseed(seed)
  !         --> e.g. seed = 327_i4
  !
  !         integer(i8), dimension(3) :: seed
  !         call get_timeseed(seed)
  !         --> e.g. seed = (/ 327_i8, 1327_i8, 2327_i8 /)

  !     LITERATURE
  !         

  !     HISTORY
  !         Written,  Juliane Mai, Aug 2012

  ! ------------------------------------------------------------------

  subroutine get_timeseed_i4_0d(seed)
    
    implicit none
    integer(i4), intent(inout)  :: seed

    ! local variables
    integer(i4), dimension(8)  :: time_array

    call date_and_time(values=time_array)
    seed = &
         time_array(5) * 3600000_i4  + &   ! hour
         time_array(6) * 60000_i4    + &   ! minutes
         time_array(7) * 1000_i4     + &   ! seconds
         time_array(8) * 1_i4              ! milliseconds
    
  end subroutine get_timeseed_i4_0d

  subroutine get_timeseed_i4_1d(seed)
    
    implicit none
    integer(i4), dimension(:), intent(inout)  :: seed

    ! local variables
    integer(i4), dimension(8)  :: time_array
    integer(i4)                :: i

    call date_and_time(values=time_array)
    seed(1) = &
         time_array(5) * 3600000_i4  + &   ! hour
         time_array(6) * 60000_i4    + &   ! minutes
         time_array(7) * 1000_i4     + &   ! seconds
         time_array(8) * 1_i4              ! milliseconds
    do i=2,size(seed)
       seed(i) = seed(i-1) + 1000_i4
    end do
    
  end subroutine get_timeseed_i4_1d

  subroutine get_timeseed_i8_0d(seed)
    
    implicit none
    integer(i8), intent(inout)  :: seed

    ! local variables
    integer(i4), dimension(8)  :: time_array

    call date_and_time(values=time_array)
    seed = &
         int(time_array(5),i8) * 3600000_i8  + &   ! hour
         int(time_array(6),i8) * 60000_i8    + &   ! minutes
         int(time_array(7),i8) * 1000_i8     + &   ! seconds
         int(time_array(8),i8) * 1_i8              ! milliseconds
    
  end subroutine get_timeseed_i8_0d

  subroutine get_timeseed_i8_1d(seed)
    
    implicit none
    integer(i8), dimension(:), intent(inout)  :: seed

    ! local variables
    integer(i4), dimension(8)  :: time_array
    integer(i4)                :: i

    call date_and_time(values=time_array)
    seed(1) = &
         int(time_array(5),i8) * 3600000_i8  + &   ! hour
         int(time_array(6),i8) * 60000_i8    + &   ! minutes
         int(time_array(7),i8) * 1000_i8     + &   ! seconds
         int(time_array(8),i8) * 1_i8              ! milliseconds
    do i=2,size(seed)
       seed(i) = seed(i-1) + 1000_i8
    end do
    
  end subroutine get_timeseed_i8_1d

  ! ------------------------------------------------------------------

  !     NAME
  !         xor4096

  !     PURPOSE
  !         Generates a uniform distributed random number based on xor4096 algorithm proposed by 
  !         Brent et.al (2006).
  !         ****************************************************************************************
  !         The original version of this source code (without multiple streams and 
  !         optional arguments) is under GNU General Public Licence
  !              xorgens.c
  !         Copyright (C) 2004 R. P. Brent.                                       
  !         This program is free software; you can redistribute it and/or         
  !         modify it under the terms of the GNU General Public License,       
  !         version 2, June 1991, as published by the Free Software Foundation.   
  !         For details see http://www.gnu.org/copyleft/gpl.html .                
  !
  !         Author: Richard P. Brent (random@rpbrent.co.uk)
  !         ****************************************************************************************
  !         The period of the generator is 
  !              (2^4096 - 1)*2^32 for single precision version and 
  !              (2^4096 - 1)*2^64 for double precision version.
  !
  !         The generator is based on bitwise XOR compositions of left- and right-shifted numbers.
  !
  !         The generator has to be called once with a non-zero seed value which initializes a new 
  !         random number stream. The subsequent calls are with seed=0 which returns the following 
  !         numbers within the beforehand initialized stream. Since the random numbers are based on 
  !         the initial seed, the precison of the seed (sp/dp) determines the precision of the 
  !         returned random number (sp/dp).
  !  
  !         If one initialize the generator with an array of seeds, one initializes n independent 
  !         streams of random numbers.
  !         Lets assume that the streams with seed 1_sp is 10, 20, 30 ...
  !                                                2_sp is 40, 50, 60 ...
  !                                                3_sp is 70, 80, 90 ...
  !         What you get is:
  !         1st call
  !             call xor( (/ 1_SP, 2_SP, 3_SP /), RN )   --> RN = (/ 10, 40, 70 /)
  !         2nd call
  !             call xor( (/ 0_SP, 0_SP, 0_SP /), RN )   --> RN = (/ 20, 50, 80 /)
  !         3rd call
  !             call xor( (/ 0_SP, 0_SP, 0_SP /), RN )   --> RN = (/ 30, 60, 90 /)
  !
  !         Since after every initialization one looses the old stream of random numbers,
  !         it might be necessary to switch between different streams. Therefore, one needs ALL
  !         of the three optional arguments.
  !         What you get is:
  !         1st call of 1st stream
  !                                   call xor( 1_SP, RN, opt1_1, opt2_1, opt3_1 )   RN = 10
  !         2nd call of 1st stream
  !                                   call xor( 0_SP, RN, opt1_1, opt2_1, opt3_1 )   RN = 20
  !         1st call of 2nd stream
  !                                   call xor( 2_SP, RN, opt1_2, opt2_2, opt3_2 )   RN = 40
  !         3rd call of 1st stream
  !                                   call xor( 0_SP, RN, opt1_1, opt2_1, opt3_1 )   RN = 30
  !         Note: If you would have called 4 times without optional arguments, 
  !               you would have get 50 in the 4th call.

  !     CALLING SEQUENCE
  !         call xor4096(seed, rn) or
  !         call xor4096(seed, rn, i, w, x)
  
  !     INDENT(IN)
  !         integer(i4/i8) :: seed/seed(:)     value or 1D-array with non-zero seeds for 
  !                                            initialization or zero for subsequent calls

  !     INDENT(INOUT)
  !             none

  !     INDENT(OUT)
  !         integer(i4/i8)/real(sp/dp)   :: RN/RN(size(seed))       
  !                                            uniform distributed random number with
  !                                            interval:
  !                                                i4: (-2^31,2^31-1)  
  !                                                i8: (-2^63,2^63-1) 
  !                                                sp: (0.0_sp, 1.0_sp) 
  !                                                dp: (0.0_dp, 1.0_dp)                                                       

  !     INDENT(IN), OPTIONAL
  !         none

  !     INDENT(INOUT), OPTIONAL
  !         integer(i4/i8), dimension(size(seed))              :: i
  !         integer(i4/i8), dimension(size(seed))              :: w
  !         integer(i4/i8), dimension(size(seed),0:127/0:63)   :: x

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         In case of optional arguments all three (i,w,x) have to be given.
  !         If random numbers are in single precision (sp), one needs a seed, i, w, x in (i4).
  !         If random numbers are in double precision (dp), one needs a seed, i, w, x in (i8).
  !         The size of optional x array depends on precision.

  !     EXAMPLE
  !         seed = (/ 1_SP, 100_SP, 2_SP /)
  !         call xor4096(seed,RN)
  !         print*, RN --> (/ 10, 20, 30 /)
  !
  !         seed = (/ 0_SP, 0_SP, 0_SP /)
  !         call xor4096(seed,RN)
  !         print*, RN --> (/ 20, 10, 50 /)
  !         -> see also example in test_mo_xor4096 directory

  !     LITERATURE
  !         Brent RP - Some long-period random number generators using shifts and xors, 2010
  !         Brent RP - From Mersenne Primes to Rndom Number Generators, 2006
  !         L''Ecuyer P & Simard R - ACM: TestU01: A C Library for Empirical Testing of 
  !                   Random Number Generators, 2007

  !     HISTORY
  !         Written,  Juliane Mai, Nov 2011

  subroutine xor4096s_0d(seed,SingleIntegerRN,iin,win,xin)
    implicit none

    integer(i4), intent(in)              :: seed
    integer(i4), intent(out)             :: SingleIntegerRN
    integer(i4), optional, intent(inout) :: iin
    integer(i4), optional, intent(inout) :: win
    integer(i4), optional, intent(inout) :: xin(0:127)

    integer(i4)        :: wlen, r, s, a, b, c, d

    integer(i4), save  :: w
    integer(i4), save  :: x(0:127)                 ! x(0) ... x(r-1)
    integer(i4)        :: weyl = 1640531527_i4    !Z'61C88647'       ! Hexadecimal notation
    integer(i4)        :: t,v
    integer(i4), save  :: i = -1                   ! i<0 indicates first call
    integer(i4)        :: k

!$omp   threadprivate(x,i,w) 

    wlen = 32
    r = 128
    s = 95
    a = 17
    b = 12
    c = 13
    d = 15

    if ( present(iin) .and. (seed .eq. 0) ) i = iin
    if ( present(win) .and. (seed .eq. 0) ) w = win
    if ( present(xin) .and. (seed .eq. 0) ) x = xin

    If ((i .lt. 0) .or. (seed .ne. 0)) then     ! Initialization necessary
        If (seed .ne. 0) then                   ! v must be nonzero
            v = seed
        else
            v = NOT(seed)
        end if

        do k=wlen,1,-1                          ! Avoid correlations for close seeds
            ! This recurrence has period of 2^32-1
            v = IEOR(v,ISHFT(v,13))
            v = IEOR(v,ISHFT(v,-17))
            v = IEOR(v,ISHFT(v, 5))
        end do

        ! Initialize circular array
        w = v
        do k=0,r-1
            w = w + weyl
            v = IEOR(v,ISHFT(v,13))
            v = IEOR(v,ISHFT(v,-17))
            v = IEOR(v,ISHFT(v, 5))
            x(k) = v + w
        end do

        ! Discard first 4*r results (Gimeno)
        i = r-1
        do k = 4*r,1,-1
            i = IAND(i+1,r-1)
            t = x(i)
            v = x(IAND(i+(r-s),r-1))
            t = IEOR(t,ISHFT(t,a))
            t = IEOR(t,ISHFT(t,-b))
            v = IEOR(v,ISHFT(v,c))
            v = IEOR(v,IEOR(t,ISHFT(v,-d)))
            x(i) = v
        end do
    end if ! end of initialization

    ! Apart from initialization (above), this is the generator
    i = IAND(i+1,r-1)
    t = x(i)
    v = x(IAND(i+(r-s),r-1))
    t = IEOR(t,ISHFT(t,a))
    t = IEOR(t,ISHFT(t,-b))
    v = IEOR(v,ISHFT(v,c))
    v = IEOR(v,IEOR(t,ISHFT(v,-d)))
    x(i) = v

    w = w + weyl

    SingleIntegerRN = v+w

    if ( present(iin) ) then
        iin=i
    End if

    if ( present(win) ) then
        win=w
    End if

    If ( present(xin) ) then
        xin=x
    End if
  end subroutine xor4096s_0d

!******************************************************************************************

  subroutine xor4096s_1d(seed,SingleIntegerRN,iin,win,xin)
    implicit none

    integer(i4), dimension(:),          intent(in)     :: seed
    integer(i4), dimension(size(seed,1)), intent(out)    :: SingleIntegerRN
    integer(i4), optional, dimension(size(seed,1)),       intent(inout) :: iin
    integer(i4), optional, dimension(size(seed,1)),       intent(inout) :: win
    integer(i4), optional, dimension(size(seed,1),0:127), intent(inout) :: xin

    integer(i4)                        :: m
    integer(i4)                        :: wlen, r, s, a, b, c, d
    integer(i4)                        :: weyl = 1640531527_i4    !Z'61C88647'       ! Hexadecimal notation
    integer(i4)                        :: k, j
    integer(i4), dimension(size(seed,1)) :: t,v
    integer(i4), dimension(:,:), allocatable, save   :: x               ! x(0) ... x(r-1)
    integer(i4), dimension(:),   allocatable, save   :: i,w             ! i<0 indicates first call

!$omp   threadprivate(x,i,w) 

    if ( present(iin) .and. (Any(seed .eq. 0)) ) i = iin
    if ( present(win) .and. (Any(seed .eq. 0)) ) w = win
    if ( present(xin) .and. (Any(seed .eq. 0)) ) x = xin

    m = size(seed,1)
    if (.not. allocated(i)) then
       allocate(i(m))
       i = -1
    end if
    if (.not. allocated(x)) then
       allocate(x(m,0:127))
    end if
    if (.not. allocated(w)) then
       allocate(w(m))
    end if

    wlen = 32
    r = 128
    s = 95
    a = 17
    b = 12
    c = 13
    d = 15

    Do j = 1, m
       If ((i(j) .lt. 0) .or. (seed(j) .ne. 0)) then     ! Initialization necessary
          If (seed(j) .ne. 0) then                   ! v must be nonzero
             v(j) = seed(j)
          else
             v(j) = NOT(seed(j))
          end if

          do k=wlen,1,-1                          ! Avoid correlations for close seeds
             ! This recurrence has period of 2^32-1
             v(j) = IEOR(v(j),ISHFT(v(j),13))
             v(j) = IEOR(v(j),ISHFT(v(j),-17))
             v(j) = IEOR(v(j),ISHFT(v(j), 5))
          end do

          ! Initialize circular array
          w(j) = v(j)
          do k=0,r-1
             w(j) = w(j) + weyl
             v(j) = IEOR(v(j),ISHFT(v(j),13))
             v(j) = IEOR(v(j),ISHFT(v(j),-17))
             v(j) = IEOR(v(j),ISHFT(v(j), 5))
             x(j,k) = v(j) + w(j)
          end do

          ! Discard first 4*r results (Gimeno)
          i(j) = r-1
          do k = 4*r,1,-1
             i(j) = IAND(i(j)+1,r-1)
             t(j) = x(j,i(j))
             v(j) = x(j,IAND(i(j)+(r-s),r-1))
             t(j) = IEOR(t(j),ISHFT(t(j),a))
             t(j) = IEOR(t(j),ISHFT(t(j),-b))
             v(j) = IEOR(v(j),ISHFT(v(j),c))
             v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
             x(j,i(j)) = v(j)
          end do
       end if ! end of initialization
    end do

    ! Apart from initialization (above), this is the generator

    do j=1,m
       i(j) = IAND(i(j)+1,r-1)
       t(j) = x(j,i(j))
       v(j) = x(j,IAND(i(j)+(r-s),r-1))
       t(j) = IEOR(t(j),ISHFT(t(j),a))
       t(j) = IEOR(t(j),ISHFT(t(j),-b))
       v(j) = IEOR(v(j),ISHFT(v(j),c))
       v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
       x(j,i(j)) = v(j)

       w(j) = w(j) + weyl
    end do

    SingleIntegerRN = v+w

    if ( present(iin) ) then
        iin=i
    End if

    if ( present(win) ) then
        win=w
    End if

    If ( present(xin) ) then
        xin=x
    End if

  end subroutine xor4096s_1d

  !******************************************************************************************

  subroutine xor4096f_0d(seed,SingleRealRN,iin,win,xin)

    implicit none

    integer(i4),  intent(in)  :: seed
    real(SP),      intent(out) :: SingleRealRN
    integer(i4), optional, intent(inout) :: iin
    integer(i4), optional, intent(inout) :: win
    integer(i4), optional, intent(inout) :: xin(0:127)

    integer(i4)        :: wlen, r, s, a, b, c, d
    integer(i4), save  :: w
    integer(i4), save  :: x(0:127)                 ! x(0) ... x(r-1)
    integer(i4)        :: weyl = 1640531527_i4    !Z'61C88647'       ! Hexadecimal notation
    integer(i4)        :: t,v
    integer(i4), save  :: i = -1                   ! i<0 indicates first call
    integer(i4)        :: k

    real(SP)            :: t24 = 1.0_SP/16777216.0_SP     ! = 0.5^24 = 1/2^24

!$omp   threadprivate(x,i,w) 

    ! produces a 24bit Integer Random Number (0...16777216) and
    ! scales it afterwards to (0.0,1.0)

    wlen = 32
    r = 128
    s = 95
    a = 17
    b = 12
    c = 13
    d = 15

    if ( present(iin) .and. (seed .eq. 0) ) i = iin
    if ( present(win) .and. (seed .eq. 0) ) w = win
    if ( present(xin) .and. (seed .eq. 0) ) x = xin

    If ((i .lt. 0) .or. (seed .ne. 0)) then     ! Initialization necessary
        If (seed .ne. 0) then                   ! v must be nonzero
            v = seed
        else
            v = NOT(seed)
        end if

        do k=wlen,1,-1                          ! Avoid correlations for close seeds
            ! This recurrence has period of 2^32-1
            v = IEOR(v,ISHFT(v,13))
            v = IEOR(v,ISHFT(v,-17))
            v = IEOR(v,ISHFT(v, 5))
        end do

        ! Initialize circular array
        w = v
        do k=0,r-1
            w = w + weyl
            v = IEOR(v,ISHFT(v,13))
            v = IEOR(v,ISHFT(v,-17))
            v = IEOR(v,ISHFT(v, 5))
            x(k) = v + w
        end do

        ! Discard first 4*r results (Gimeno)
        i = r-1
        do k = 4*r,1,-1
            i = IAND(i+1,r-1)
            t = x(i)
            v = x(IAND(i+(r-s),r-1))
            t = IEOR(t,ISHFT(t,a))
            t = IEOR(t,ISHFT(t,-b))
            v = IEOR(v,ISHFT(v,c))
            v = IEOR(v,IEOR(t,ISHFT(v,-d)))
            x(i) = v
        end do
    end if ! end of initialization

    ! Apart from initialization (above), this is the generator
    v = 0_i4
    Do While (v .eq. 0_i4)
        i = IAND(i+1,r-1)
        t = x(i)
        v = x(IAND(i+(r-s),r-1))
        t = IEOR(t,ISHFT(t,a))
        t = IEOR(t,ISHFT(t,-b))
        v = IEOR(v,ISHFT(v,c))
        v = IEOR(v,IEOR(t,ISHFT(v,-d)))
        x(i) = v
        w = w + weyl
        v = v + w
        v = ISHFT(v,-8)
    End Do

    SingleRealRN = t24*v

    if ( present(iin) ) then
        iin=i
    End if

    if ( present(win) ) then
        win=w
    End if

    If ( present(xin) ) then
        xin=x
    End if

  end subroutine xor4096f_0d

  !******************************************************************************************

  subroutine xor4096f_1d(seed,SingleRealRN,iin,win,xin)

    implicit none

    integer(i4), dimension(:),          intent(in)  :: seed
    real(SP),     dimension(size(seed)), intent(out) :: SingleRealRN
    integer(i4), optional, dimension(size(seed)),       intent(inout) :: iin
    integer(i4), optional, dimension(size(seed)),       intent(inout) :: win
    integer(i4), optional, dimension(size(seed),0:127), intent(inout) :: xin

    integer(i4)                         :: m
    integer(i4)                         :: wlen, r, s, a, b, c, d
    integer(i4)                         :: weyl =  1640531527_i4              !Z'61C88647' = Hexadecimal notation
    integer(i4)                         :: k, j
    real(SP), save                       :: t24 = 1.0_SP/16777216.0_SP      ! = 0.5^24 = 1/2^24
    integer(i4), dimension(size(seed))  :: t,v
    integer(i4), dimension(:,:), allocatable, save  :: x                   ! x(0) ... x(r-1)
    integer(i4), dimension(:),   allocatable, save  :: i,w                 ! i<0 indicates first call

!$omp   threadprivate(x,i,w) 

    m= size(seed)

    if ( present(iin) .and. (Any(seed .eq. 0)) ) i = iin
    if ( present(win) .and. (Any(seed .eq. 0)) ) w = win
    if ( present(xin) .and. (Any(seed .eq. 0)) ) x = xin

    ! produces a 24bit Integer Random Number (0...16777216) and
    ! scales it afterwards to (0.0,1.0)

    if (.not. allocated(i)) then
       allocate(i(m))
       i = -1
    end if
    if (.not. allocated(x)) then
       allocate(x(m,0:127))
    end if
    if (.not. allocated(w)) then
       allocate(w(m))
    end if

    wlen = 32
    r = 128
    s = 95
    a = 17
    b = 12
    c = 13
    d = 15

    Do j = 1,m !Loop over every stream
       If ((i(j) .lt. 0) .or. (seed(j) .ne. 0)) then     ! Initialization necessary
          If (seed(j) .ne. 0) then                   ! v must be nonzero
             v(j) = seed(j)
          else
             v(j) = NOT(seed(j))
          end if

          do k=wlen,1,-1                          ! Avoid correlations for close seeds
             ! This recurrence has period of 2^32-1
             v(j) = IEOR(v(j),ISHFT(v(j),13))
             v(j) = IEOR(v(j),ISHFT(v(j),-17))
             v(j) = IEOR(v(j),ISHFT(v(j), 5))
          end do

          ! Initialize circular array
          w(j) = v(j)
          do k=0,r-1
             w(j) = w(j) + weyl
             v(j) = IEOR(v(j),ISHFT(v(j),13))
             v(j) = IEOR(v(j),ISHFT(v(j),-17))
             v(j) = IEOR(v(j),ISHFT(v(j), 5))
             x(j,k) = v(j) + w(j)
          end do

          ! Discard first 4*r results (Gimeno)
          i(j) = r-1
          do k = 4*r,1,-1
             i(j) = IAND(i(j)+1,r-1)
             t(j) = x(j,i(j))
             v(j) = x(j,IAND(i(j)+(r-s),r-1))
             t(j) = IEOR(t(j),ISHFT(t(j),a))
             t(j) = IEOR(t(j),ISHFT(t(j),-b))
             v(j) = IEOR(v(j),ISHFT(v(j),c))
             v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
             x(j,i(j)) = v(j)
          end do
       end if ! end of initialization
    end do

    ! Apart from initialization (above), this is the generator
    v = 0_i4
    Do j=1,m
       Do While (v(j) .eq. 0_i4)
          i(j) = IAND(i(j)+1,r-1)
          t(j) = x(j,i(j))
          v(j) = x(j,IAND(i(j)+(r-s),r-1))
          t(j) = IEOR(t(j),ISHFT(t(j),a))
          t(j) = IEOR(t(j),ISHFT(t(j),-b))
          v(j) = IEOR(v(j),ISHFT(v(j),c))
          v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
          x(j,i(j)) = v(j)
          w(j) = w(j) + weyl
          v(j) = v(j) + w(j)
          v(j) = ISHFT(v(j),-8)
       End Do
    End Do

    SingleRealRN = t24*v

    if ( present(iin) ) then
        iin=i
    End if

    if ( present(win) ) then
        win=w
    End if

    If ( present(xin) ) then
        xin=x
    End if

  end subroutine xor4096f_1d

!******************************************************************************************

  subroutine xor4096l_0d(seed,DoubleIntegerRN,iin,win,xin)

    implicit none

    integer(i8), intent(in)  :: seed
    integer(i8), intent(out) :: DoubleIntegerRN
    integer(i8), optional, intent(inout) :: iin
    integer(i8), optional, intent(inout) :: win
    integer(i8), optional, intent(inout) :: xin(0:63)

    integer(i8)        :: wlen, r, s, a, b, c, d
    integer(i8), save  :: w
    integer(i8), save  :: x(0:63)                  ! x(0) ... x(r-1)
    integer(i8)        :: weyl = 7046029254386353131_i8
    integer(i8)        :: t,v
    integer(i8), save  :: i = -1                   ! i<0 indicates first call
    integer(i8)        :: k

!$omp   threadprivate(x,i,w) 

    if ( present(iin) .and. (seed .eq. 0) ) i = iin
    if ( present(win) .and. (seed .eq. 0) ) w = win
    if ( present(xin) .and. (seed .eq. 0) ) x = xin

    wlen = 64_i8
    r = 64_i8
    s = 53_i8
    a = 33_i8
    b = 26_i8
    c = 27_i8
    d = 29_i8

    If ((i .lt. 0) .or. (seed .ne. 0)) then     ! Initialization necessary
        If (seed .ne. 0) then                   ! v must be nonzero
            v = seed
        else
            v = NOT(seed)
        end if

        do k=wlen,1,-1                          ! Avoid correlations for close seeds
            ! This recurrence has period of 2^64-1
            v = IEOR(v,ISHFT(v,7))
            v = IEOR(v,ISHFT(v,-9))
        end do

        ! Initialize circular array
        w = v
        do k=0,r-1
            w = w + weyl
            v = IEOR(v,ISHFT(v,7))
            v = IEOR(v,ISHFT(v,-9))
            x(k) = v + w
        end do

        ! Discard first 4*r results (Gimeno)
        i = r-1
        do k = 4*r,1,-1
            i = IAND(i+1,r-1)
            t = x(i)
            v = x(IAND(i+(r-s),r-1))
            t = IEOR(t,ISHFT(t,a))
            t = IEOR(t,ISHFT(t,-b))
            v = IEOR(v,ISHFT(v,c))
            v = IEOR(v,IEOR(t,ISHFT(v,-d)))
            x(i) = v
        end do
    end if ! end of initialization

    ! Apart from initialization (above), this is the generator
    i = IAND(i+1,r-1)
    t = x(i)
    v = x(IAND(i+(r-s),r-1))
    t = IEOR(t,ISHFT(t,a))
    t = IEOR(t,ISHFT(t,-b))
    v = IEOR(v,ISHFT(v,c))
    v = IEOR(v,IEOR(t,ISHFT(v,-d)))
    x(i) = v

    w = w + weyl

    DoubleIntegerRN = v+w

    if ( present(iin) ) then
        iin=i
    End if

    if ( present(win) ) then
        win=w
    End if

    If ( present(xin) ) then
        xin=x
    End if

  end subroutine xor4096l_0d

!******************************************************************************************

  subroutine xor4096l_1d(seed, DoubleIntegerRN, iin, win, xin)

    implicit none

    integer(i8), dimension(:), intent(in)  :: seed
    integer(i8), dimension(size(seed)), intent(out) :: DoubleIntegerRN
    integer(i8), optional, dimension(size(seed)),      intent(inout) :: iin
    integer(i8), optional, dimension(size(seed)),      intent(inout) :: win
    integer(i8), optional, dimension(size(seed),0:63), intent(inout) :: xin

    integer(i4)        :: m
    integer(i8)        :: wlen, r, s, a, b, c, d
    integer(i8)        :: weyl = 7046029254386353131_i8 
    integer(i8)        :: k, j
    integer(i8), dimension(size(seed))              :: t,v
    integer(i8), dimension(:,:), allocatable, save  :: x                  ! x(0) ... x(r-1)
    integer(i8), dimension(:),   allocatable, save  :: i,w                ! i<0 indicates first call

!$omp   threadprivate(x,i,w) 

    if ( present(iin) .and. (Any(seed .eq. 0)) ) i = iin
    if ( present(win) .and. (Any(seed .eq. 0)) ) w = win
    if ( present(xin) .and. (Any(seed .eq. 0)) ) x = xin

    m = size(seed)
    wlen = 64_i8
    r = 64_i8
    s = 53_i8
    a = 33_i8
    b = 26_i8
    c = 27_i8
    d = 29_i8

    if (.not. allocated(i)) then
       allocate(i(m))
       i = -1
    end if
    if (.not. allocated(x)) then
       allocate(x(m,0:63))
    end if
    if (.not. allocated(w)) then
       allocate(w(m))
    end if

    Do j=1,m
       If ((i(j) .lt. 0) .or. (seed(j) .ne. 0)) then     ! Initialization necessary
          If (seed(j) .ne. 0) then                   ! v must be nonzero
             v(j) = seed(j)
          else
             v(j) = NOT(seed(j))
          end if

          do k=wlen,1,-1                          ! Avoid correlations for close seeds
             ! This recurrence has period of 2^64-1
             v(j) = IEOR(v(j),ISHFT(v(j),7))
             v(j) = IEOR(v(j),ISHFT(v(j),-9))
          end do

          ! Initialize circular array
          w(j) = v(j)
          do k=0,r-1
             w(j) = w(j) + weyl
             v(j) = IEOR(v(j),ISHFT(v(j),7))
             v(j) = IEOR(v(j),ISHFT(v(j),-9))
             x(j,k) = v(j) + w(j)
          end do

          ! Discard first 4*r results (Gimeno)
          i(j) = r-1
          do k = 4*r,1,-1
             i(j) = IAND(i(j)+1,r-1)
             t(j) = x(j,i(j))
             v(j) = x(j,IAND(i(j)+(r-s),r-1))
             t(j) = IEOR(t(j),ISHFT(t(j),a))
             t(j) = IEOR(t(j),ISHFT(t(j),-b))
             v(j) = IEOR(v(j),ISHFT(v(j),c))
             v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
             x(j,i(j)) = v(j)
          end do
       end if ! end of initialization
    end do

    ! Apart from initialization (above), this is the generator
    do j=1,m
       i(j) = IAND(i(j)+1,r-1)
       t(j) = x(j,i(j))
       v(j) = x(j,IAND(i(j)+(r-s),r-1))
       t(j) = IEOR(t(j),ISHFT(t(j),a))
       t(j) = IEOR(t(j),ISHFT(t(j),-b))
       v(j) = IEOR(v(j),ISHFT(v(j),c))
       v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
       x(j,i(j)) = v(j)

       w(j) = w(j) + weyl
    end do

    DoubleIntegerRN = v+w

    if ( present(iin) ) then
        iin=i
    End if

    if ( present(win) ) then
        win=w
    End if

    If ( present(xin) ) then
        xin=x
    End if

  end subroutine xor4096l_1d

!******************************************************************************************

  subroutine xor4096d_0d(seed,DoubleRealRN,iin,win,xin)

    implicit none

    integer(i8), intent(in)  :: seed
    real(DP),     intent(out) :: DoubleRealRN
    integer(i8), optional, intent(inout) :: iin
    integer(i8), optional, intent(inout) :: win
    integer(i8), optional, intent(inout) :: xin(0:63)

    integer(i8)        :: wlen, r, s, a, b, c, d

    integer(i8), save  :: w
    integer(i8), save  :: x(0:63)                  ! x(0) ... x(r-1)
    integer(i8)        :: weyl = 7046029254386353131_i8
    integer(i8)        :: t,v
    integer(i8), save  :: i = -1                   ! i<0 indicates first call
    integer(i8)        :: k

    real(DP)            :: t53 = 1.0_DP/9007199254740992.0_DP                     ! = 0.5^53 = 1/2^53

!$omp   threadprivate(x,i,w) 

    ! produces a 53bit Integer Random Number (0...9 007 199 254 740 992) and
    ! scales it afterwards to (0.0,1.0)

    if ( present(iin) .and. (seed .eq. 0) ) i = iin
    if ( present(win) .and. (seed .eq. 0) ) w = win
    if ( present(xin) .and. (seed .eq. 0) ) x = xin

    wlen = 64_i8
    r = 64_i8
    s = 53_i8
    a = 33_i8
    b = 26_i8
    c = 27_i8
    d = 29_i8

    If ((i .lt. 0) .or. (seed .ne. 0)) then     ! Initialization necessary
        If (seed .ne. 0) then                   ! v must be nonzero
            v = seed
        else
            v = NOT(seed)
        end if

        do k=wlen,1,-1                          ! Avoid correlations for close seeds
            ! This recurrence has period of 2^64-1
            v = IEOR(v,ISHFT(v,7))
            v = IEOR(v,ISHFT(v,-9))
        end do

        ! Initialize circular array
        w = v
        do k=0,r-1
            w = w + weyl
            v = IEOR(v,ISHFT(v,7))
            v = IEOR(v,ISHFT(v,-9))
            x(k) = v + w
        end do

        ! Discard first 4*r results (Gimeno)
        i = r-1
        do k = 4*r,1,-1
            i = IAND(i+1,r-1)
            t = x(i)
            v = x(IAND(i+(r-s),r-1))
            t = IEOR(t,ISHFT(t,a))
            t = IEOR(t,ISHFT(t,-b))
            v = IEOR(v,ISHFT(v,c))
            v = IEOR(v,IEOR(t,ISHFT(v,-d)))
            x(i) = v
        end do
    end if ! end of initialization

    ! Apart from initialization (above), this is the generator
    v = 0_i8
    Do While (v .eq. 0_i8)
        i = IAND(i+1,r-1)
        t = x(i)
        v = x(IAND(i+(r-s),r-1))
        t = IEOR(t,ISHFT(t,a))
        t = IEOR(t,ISHFT(t,-b))
        v = IEOR(v,ISHFT(v,c))
        v = IEOR(v,IEOR(t,ISHFT(v,-d)))
        x(i) = v
        w = w + weyl
        v = v + w
        v = ISHFT(v,-11)
    End Do

    DoubleRealRN = t53*v

    if ( present(iin) ) then
        iin=i
    End if

    if ( present(win) ) then
        win=w
    End if

    If ( present(xin) ) then
        xin=x
    End if

  end subroutine xor4096d_0d

!******************************************************************************************

  subroutine xor4096d_1d(seed,DoubleRealRN,iin,win,xin)

    implicit none

    integer(i8), dimension(:),          intent(in)  :: seed
    real(DP),     dimension(size(seed)), intent(out) :: DoubleRealRN
    integer(i8), optional, dimension(size(seed)),      intent(inout) :: iin
    integer(i8), optional, dimension(size(seed)),      intent(inout) :: win
    integer(i8), optional, dimension(size(seed),0:63), intent(inout) :: xin

    integer(i4)                       :: m
    integer(i8)                       :: wlen, r, s, a, b, c, d
    integer(i8)                       :: weyl = 7046029254386353131_i8
    real(DP)                           :: t53  = 1.0_DP/9007199254740992.0_DP  ! = 0.5^53 = 1/2^53
    integer(i8)                       :: k,j
    integer(i8), dimension(size(seed))              :: t,v
    integer(i8), dimension(:,:), allocatable, save  :: x       ! x(0) ... x(r-1)
    integer(i8), dimension(:),   allocatable, save  :: w
    integer(i8), dimension(:),   allocatable, save  :: i       ! i<0 indicates first call

!$omp   threadprivate(x,i,w) 

    ! produces a 53bit Integer Random Number (0...9 007 199 254 740 992) and
    ! scales it afterwards to (0.0,1.0)

    if ( present(iin) .and. (Any(seed .eq. 0)) ) i = iin
    if ( present(win) .and. (Any(seed .eq. 0)) ) w = win
    if ( present(xin) .and. (Any(seed .eq. 0)) ) x = xin

    m = size(seed)

    if (.not. allocated(i)) then
       allocate(i(m))
       i = -1
    end if
    if (.not. allocated(x)) then
       allocate(x(m,0:63))
    end if
    if (.not. allocated(w)) then
       allocate(w(m))
    end if

    wlen = 64_i8
    r = 64_i8
    s = 53_i8
    a = 33_i8
    b = 26_i8
    c = 27_i8
    d = 29_i8

    Do j=1,m
       If ((i(j) .lt. 0) .or. (seed(j) .ne. 0)) then     ! Initialization necessary
          If (seed(j) .ne. 0) then                   ! v must be nonzero
             v(j) = seed(j)
          else
             v(j) = NOT(seed(j))
          end if

          do k=wlen,1,-1                          ! Avoid correlations for close seeds
             ! This recurrence has period of 2^64-1
             v(j) = IEOR(v(j),ISHFT(v(j),7))
             v(j) = IEOR(v(j),ISHFT(v(j),-9))
          end do

          ! Initialize circular array
          w(j) = v(j)
          do k=0,r-1
             w(j) = w(j) + weyl
             v(j) = IEOR(v(j),ISHFT(v(j),7))
             v(j) = IEOR(v(j),ISHFT(v(j),-9))
             x(j,k) = v(j) + w(j)
          end do

          ! Discard first 4*r results (Gimeno)
          i(j) = r-1
          do k = 4*r,1,-1
             i(j) = IAND(i(j)+1,r-1)
             t(j) = x(j,i(j))
             v(j) = x(j,IAND(i(j)+(r-s),r-1))
             t(j) = IEOR(t(j),ISHFT(t(j),a))
             t(j) = IEOR(t(j),ISHFT(t(j),-b))
             v(j) = IEOR(v(j),ISHFT(v(j),c))
             v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
             x(j,i(j)) = v(j)
          end do
       end if ! end of initialization
    end do

    ! Apart from initialization (above), this is the generator
    v = 0_i8
    Do j=1,m
       Do While (v(j) .eq. 0_i8)
          i(j) = IAND(i(j)+1,r-1)
          t(j) = x(j,i(j))
          v(j) = x(j,IAND(i(j)+(r-s),r-1))
          t(j) = IEOR(t(j),ISHFT(t(j),a))
          t(j) = IEOR(t(j),ISHFT(t(j),-b))
          v(j) = IEOR(v(j),ISHFT(v(j),c))
          v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
          x(j,i(j)) = v(j)
          w(j) = w(j) + weyl
          v(j) = v(j) + w(j)
          v(j) = ISHFT(v(j),-11)
       End Do
    End Do

    DoubleRealRN = t53*v

    if ( present(iin) ) then
        iin=i
    End if

    if ( present(win) ) then
        win=w
    End if

    If ( present(xin) ) then
        xin=x
    End if

  end subroutine xor4096d_1d

!******************************************************************************************

  ! ------------------------------------------------------------------

  !     NAME
  !         xor4096g

  !     PURPOSE
  !         Generates a gaussian distributed random number. First, a uniform distributed random 
  !         number based on xor4096 algorithm is generated. Second, this number is transformed
  !         using the Polar method of Box-Mueller-transform to calculate the gaussian distributed
  !         numbers.
  !
  !         ****************************************************************************************
  !         The original version of this source code (without multiple streams, 
  !         optional arguments and gaussian distributed RN) is under GNU General Public Licence
  !              xorgens.c
  !         Copyright (C) 2004 R. P. Brent.                                       
  !         This program is free software; you can redistribute it and/or         
  !         modify it under the terms of the GNU General Public License,       
  !         version 2, June 1991, as published by the Free Software Foundation.   
  !         For details see http://www.gnu.org/copyleft/gpl.html .                
  !
  !         Author: Richard P. Brent (random@rpbrent.co.uk)
  !         ****************************************************************************************
  !
  !         The generator has to be called once with a non-zero seed value which initializes a new 
  !         random number stream. The subsequent calls are with seed=0 which returns the following 
  !         numbers within the beforehand initialized stream. Since the random numbers are based on 
  !         the initial seed, the precison of the seed (sp/dp) determines the precision of the 
  !         returned random number (sp/dp).
  !
  !         The returned values are gaussian distributed with mean 0 and variance 1.
  !  
  !         The polar method of the Box-Mueller transform transforms two uniform distributed
  !         random numbers u1 and u2 into two gaussian distributed random numbers x1 and x2.
  !         First, two uniform numbers a1, a2 distributed between (-1,1) are calculated:
  !                      a1 = 2u1 - 1
  !                      a2 = 2u2 - 1
  !         Second, q = a1^2 + a2^2 is calculated. If (q = 0) or (q > 1) one has to generate 
  !         a new x1 and x2, since a divison by q is needed afterwards and q needs to be distributed
  !         within the unit circle.
  !         Third,
  !                      p = Sqrt( -2 LN(q) / q )
  !         is calculated.
  !         The numbers z1 and z2 with
  !                      z1 = a1 * p
  !                      z2 = a2 * p
  !         are then gaussian distributed with mean 0 and variance 1.

  !     CALLING SEQUENCE
  !         call xor4096g(seed, rn) or
  !         call xor4096g(seed, rn, i, w, x, Flag, y)
  
  !     INDENT(IN)
  !         integer(i4/i8) :: seed/seed(:)     value or 1D-array with non-zero seeds for 
  !                                            initialization or zero for subsequent calls

  !     INDENT(INOUT)
  !             none

  !     INDENT(OUT)
  !         real(sp/dp)   :: RN/RN(size(seed))       
  !                                            gaussian distributed random number with
  !                                            interval:
  !                                                sp: RN ~ N(0.0_sp, 1.0_sp) 
  !                                                dp: RN ~ N(0.0_dp, 1.0_dp)                                                       

  !     INDENT(IN), OPTIONAL
  !         none

  !     INDENT(INOUT), OPTIONAL
  !         integer(i4/i8), dimension(size(seed))              :: i
  !         integer(i4/i8), dimension(size(seed))              :: w
  !         integer(i4/i8), dimension(size(seed),0:127/0:63)   :: x
  !         integer(i4/i8), dimension(size(seed))              :: Flag
  !         real(sp/dp),    dimension(size(seed))              :: y
  !         

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         In case of optional arguments all five (i,w,x,Flag,y) have to be given.
  !         If random numbers are in single precision (sp), one needs a seed, i, w, x, Flag in (i4)
  !         and y in (sp).
  !         If random numbers are in double precision (dp), one needs a seed, i, w, x, Flag in (i8)
  !         and y in (dp).
  !         The size of optional x array depends on precision.

  !     EXAMPLE
  !         seed = (/ 1_SP, 100_SP, 2_SP /)
  !         call xor4096g(seed,RN)
  !         print*, RN --> (/ 0.1, 0.05, 0.3 /)
  !
  !         seed = (/ 0_SP, 0_SP, 0_SP /)
  !         call xor4096(seed,RN)
  !         print*, RN --> (/ 0.2, 0.9, 0.4 /)
  !         -> see also example in test_mo_xor4096 directory

  !     LITERATURE
  !         Brent RP - Some long-period random number generators using shifts and xors, 2010
  !         Brent RP - From Mersenne Primes to Rndom Number Generators, 2006
  !         L''Ecuyer P & Simard R - ACM: TestU01: A C Library for Empirical Testing of 
  !                   Random Number Generators, 2007
  !         http://en.wikipedia.org/wiki/Marsaglia_polar_method
  !         http://de.wikipedia.org/wiki/Polar-Methode

  !     HISTORY
  !         Written,  Juliane Mai, Nov 2011
  !

subroutine xor4096gf_0d(seed,SingleRealRN,iIn,wIn,xIn,FlagIn,y2In)

    implicit none

    integer(i4),                intent(in)        :: seed
    real(SP),                    intent(out)      :: SingleRealRN
    integer(i4), optional,      intent(inout)     :: iin
    integer(i4), optional,      intent(inout)     :: win
    integer(i4), optional,      intent(inout)     :: xin(0:127)
    integer(i4), optional,      intent(inout)     :: FlagIn
    real(SP),     optional,      intent(inout)    :: y2In

    integer(i4)        :: wlen, r, s, a, b, c, d
    integer(i4), save  :: w
    integer(i4), save  :: x(0:127)                 ! x(0) ... x(r-1)
    integer(i4)        :: weyl = 1640531527_i4    
    integer(i4)        :: t,v
    integer(i4), save  :: i = -1                   ! i<0 indicates first call
    integer(i4)        :: k
    real(SP)           :: t24 = 1.0_SP/16777216.0_SP     ! = 0.5^24 = 1/2^24

    real(SP)            :: rn1, rn2               ! uniform random numbers
    real(SP)            :: x1,x2,y1,ww            ! for Box-Mueller transform
    integer(i4),save    :: Flag = 1               ! if Flag = 1 return y1 else return y2
    real(SP),save       :: y2

!$omp   threadprivate(x,i,w,y2,flag) 

    ! produces a 24bit Integer Random Number (0...16777216) and
    ! scales it afterwards to (0.0,1.0)
    ! transform using polar-method

    if ( present(iin) .and. seed .eq. 0 )    i    = iin
    if ( present(win) .and. seed .eq. 0 )    w    = win
    if ( present(xin) .and. seed .eq. 0 )    x    = xin
    if ( present(Flagin) .and. seed .eq. 0 ) Flag = Flagin
    if ( present(y2in) .and. seed .eq. 0 )   y2   = y2in

    wlen = 32
    r = 128
    s = 95
    a = 17
    b = 12
    c = 13
    d = 15

    If ((i .lt. 0) .or. (seed .ne. 0)) then     ! Initialization necessary
       If (seed .ne. 0) then                   ! v must be nonzero
          v = seed
       else
          v = NOT(seed)
       end if

       do k=wlen,1,-1                          ! Avoid correlations for close seeds
          ! This recurrence has period of 2^32-1
          v = IEOR(v,ISHFT(v,13))
          v = IEOR(v,ISHFT(v,-17))
          v = IEOR(v,ISHFT(v, 5))
       end do

       ! Initialize circular array
       w = v
       do k=0,r-1
          w = w + weyl
          v = IEOR(v,ISHFT(v,13))
          v = IEOR(v,ISHFT(v,-17))
          v = IEOR(v,ISHFT(v, 5))
          x(k) = v + w
       end do

       ! Discard first 4*r results (Gimeno)
       i = r-1
       do k = 4*r,1,-1
          i = IAND(i+1,r-1)
          t = x(i)
          v = x(IAND(i+(r-s),r-1))
          t = IEOR(t,ISHFT(t,a))
          t = IEOR(t,ISHFT(t,-b))
          v = IEOR(v,ISHFT(v,c))
          v = IEOR(v,IEOR(t,ISHFT(v,-d)))
          x(i) = v
       end do
       Flag = 1
    end if ! end of initialization

    If (Flag .eq. 1) then
    ! Polar method of Box-Mueller-transform to generate Gaussian distributed random number
    ww = 1.0_SP
    do while (ww .ge. 1.0_SP)

       ! Apart from initialization (above), this is the generator
       v = 0_i4
       Do While (v .eq. 0_i4)
          i = IAND(i+1,r-1)
          t = x(i)
          v = x(IAND(i+(r-s),r-1))
          t = IEOR(t,ISHFT(t,a))
          t = IEOR(t,ISHFT(t,-b))
          v = IEOR(v,ISHFT(v,c))
          v = IEOR(v,IEOR(t,ISHFT(v,-d)))
          x(i) = v
          w = w + weyl
          v = v + w
          v = ISHFT(v,-8)
       End Do

       rn1 = t24*v

       v = 0_i4
       Do While (v .eq. 0_i4)
          i = IAND(i+1,r-1)
          t = x(i)
          v = x(IAND(i+(r-s),r-1))
          t = IEOR(t,ISHFT(t,a))
          t = IEOR(t,ISHFT(t,-b))
          v = IEOR(v,ISHFT(v,c))
          v = IEOR(v,IEOR(t,ISHFT(v,-d)))
          x(i) = v
          w = w + weyl
          v = v + w
          v = ISHFT(v,-8)
       End Do

       rn2 = t24*v

       x1 = 2.0_SP * rn1 -1.0_SP
       x2 = 2.0_SP * rn2 -1.0_SP

       ww = x1*x1 + x2*x2
    end do ! end of polar method

    ww = Sqrt( (-2.0_SP * Log(ww)) / ww)
    y1 = x1 * ww
    y2 = x2 * ww

    end if  ! Only if Flag = 1

    If (Flag .eq. 1) then
        Flag = 2
        SingleRealRN = y1
    else
        Flag = 1
        SingleRealRN = y2
    end if

    If ( present(iin) )    iin=i
    If ( present(win) )    win=w
    If ( present(xin) )    xin=x
    If ( present(Flagin) ) Flagin=Flag
    If ( present(y2in) )   y2in=y2

end subroutine xor4096gf_0d

!******************************************************************************************

subroutine xor4096gf_1d(seed,SingleRealRN,iin,win,xin,FlagIn,y2In)

    implicit none

    integer(i4), dimension(:),                          intent(in)        :: seed
    real(SP),    dimension(size(seed)),                 intent(out)       :: SingleRealRN
    integer(i4), dimension(size(seed)),       optional, intent(inout)     :: iin
    integer(i4), dimension(size(seed)),       optional, intent(inout)     :: win
    integer(i4), dimension(size(seed),0:127), optional, intent(inout)     :: xin
    integer(i4), dimension(size(seed)),       optional, intent(inout)     :: FlagIn
    real(SP),    dimension(size(seed)),       optional, intent(inout)     :: y2In

    integer(i4)                         :: m
    integer(i4)                         :: wlen, r, s, a, b, c, d
    integer(i4)                         :: weyl =  1640531527_i4        
    integer(i4)                         :: k, j
    real(SP)                            :: t24 = 1.0_SP/16777216.0_SP      ! = 0.5^24 = 1/2^24
    integer(i4), dimension(size(seed))  :: t,v
    integer(i4), dimension(:,:), allocatable, save  :: x                   ! x(0) ... x(r-1)
    integer(i4), dimension(:),   allocatable, save  :: i,w                 ! i<0 indicates first call

    real(SP),    dimension(size(seed))              :: rn1, rn2     ! uniform random numbers
    real(SP),    dimension(size(seed))              :: x1,x2,y1,ww  ! for Box-Mueller transform
    real(SP),    dimension(:), allocatable, save    :: y2
    integer(i4), dimension(:), allocatable, save    :: Flag         ! if Flag = 1 return y1 else return y2

!$omp   threadprivate(x,i,w,y2,flag) 

    m= size(seed)

    if ( present(iin) .and. (Any(seed .eq. 0)) )    i    = iin
    if ( present(win) .and. (Any(seed .eq. 0)) )    w    = win
    if ( present(xin) .and. (Any(seed .eq. 0)) )    x    = xin
    if ( present(Flagin) .and. (Any(seed .eq. 0)) ) Flag = Flagin
    if ( present(y2in) .and. (Any(seed .eq. 0)) )   y2   = y2in

    ! produces a 24bit Integer Random Number (0...16777216) and
    ! scales it afterwards to (0.0,1.0)
    ! transform using polar-method

    wlen = 32
    r = 128
    s = 95
    a = 17
    b = 12
    c = 13
    d = 15

    if (.not. allocated(i)) then
       allocate(i(m))
       i = -1
    end if
    if (.not. allocated(x)) then
       allocate(x(m,0:127))
    end if
    if (.not. allocated(w)) then
       allocate(w(m))
    end if
    if (.not. allocated(Flag)) then
       allocate(Flag(m))
       Flag = 1
    end if
    if (.not. allocated(y2)) then
       allocate(y2(m))
    end if

    Do j=1,m !Loop over every stream
       If ((i(j) .lt. 0) .or. (seed(j) .ne. 0)) then     ! Initialization necessary
          If (seed(j) .ne. 0) then                   ! v must be nonzero
             v(j) = seed(j)
          else
             v(j) = NOT(seed(j))
          end if

          do k=wlen,1,-1                          ! Avoid correlations for close seeds
             ! This recurrence has period of 2^32-1
             v(j) = IEOR(v(j),ISHFT(v(j),13))
             v(j) = IEOR(v(j),ISHFT(v(j),-17))
             v(j) = IEOR(v(j),ISHFT(v(j), 5))
          end do

          ! Initialize circular array
          w(j) = v(j)
          do k=0,r-1
             w(j) = w(j) + weyl
             v(j) = IEOR(v(j),ISHFT(v(j),13))
             v(j) = IEOR(v(j),ISHFT(v(j),-17))
             v(j) = IEOR(v(j),ISHFT(v(j), 5))
             x(j,k) = v(j) + w(j)
          end do

          ! Discard first 4*r results (Gimeno)
          i(j) = r-1
          do k = 4*r,1,-1
             i(j) = IAND(i(j)+1,r-1)
             t(j) = x(j,i(j))
             v(j) = x(j,IAND(i(j)+(r-s),r-1))
             t(j) = IEOR(t(j),ISHFT(t(j),a))
             t(j) = IEOR(t(j),ISHFT(t(j),-b))
             v(j) = IEOR(v(j),ISHFT(v(j),c))
             v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
             x(j,i(j)) = v(j)
          end do
          Flag(j) = 1
       end if ! end of initialization
    end do

    Do j=1,m        !Loop over every stream
        If (Flag(j) .eq. 1) then
        ! Polar method of Box-Mueller-transform to generate Gaussian distributed random number
            ww(j) = 1.0_SP
            do while (ww(j) .ge. 1.0_SP)

                ! Apart from initialization (above), this is the generator
                v(j) = 0_i4
                Do While (v(j) .eq. 0_i4)
                    i(j) = IAND(i(j)+1,r-1)
                    t(j) = x(j,i(j))
                    v(j) = x(j,IAND(i(j)+(r-s),r-1))
                    t(j) = IEOR(t(j),ISHFT(t(j),a))
                    t(j) = IEOR(t(j),ISHFT(t(j),-b))
                    v(j) = IEOR(v(j),ISHFT(v(j),c))
                    v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
                    x(j,i(j)) = v(j)
                    w(j) = w(j) + weyl
                    v(j) = v(j) + w(j)
                    v(j) = ISHFT(v(j),-8)
                End Do

                rn1(j) = t24*v(j)

                v(j) = 0_i4
                Do While (v(j) .eq. 0_i4)
                    i(j) = IAND(i(j)+1,r-1)
                    t(j) = x(j,i(j))
                    v(j) = x(j,IAND(i(j)+(r-s),r-1))
                    t(j) = IEOR(t(j),ISHFT(t(j),a))
                    t(j) = IEOR(t(j),ISHFT(t(j),-b))
                    v(j) = IEOR(v(j),ISHFT(v(j),c))
                    v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
                    x(j,i(j)) = v(j)
                    w(j) = w(j) + weyl
                    v(j) = v(j) + w(j)
                    v(j) = ISHFT(v(j),-8)
                End Do

                rn2(j) = t24*v(j)

                x1(j) = 2.0_SP * rn1(j) -1.0_SP
                x2(j) = 2.0_SP * rn2(j) -1.0_SP
                ww(j) = x1(j)*x1(j) + x2(j)*x2(j)
            end do ! end of polar method

            ww(j) = Sqrt( (-2.0_SP * Log(ww(j))) / ww(j))
            y1(j) = x1(j) * ww(j)
            y2(j) = x2(j) * ww(j)
        end if  ! Only if Flag = 1
    end do ! Loop over each stream


    Do j=1,m
        If (Flag(j) .eq. 1) then
            Flag(j) = 2
            SingleRealRN(j) = y1(j)
        else
            Flag(j) = 1
            SingleRealRN(j) = y2(j)
        end if
    end Do

    If ( present(iin) )    iin=i
    If ( present(win) )    win=w
    If ( present(xin) )    xin=x
    If ( present(Flagin) ) Flagin=Flag
    If ( present(y2in) )   y2in=y2

end subroutine xor4096gf_1d

!******************************************************************************************

subroutine xor4096gd_0d(seed,DoubleRealRN,iin,win,xin,FlagIn,y2In)

    implicit none

    integer(i8),                intent(in)        :: seed
    real(DP),                   intent(out)       :: DoubleRealRN
    integer(i8), optional,      intent(inout)     :: iin
    integer(i8), optional,      intent(inout)     :: win
    integer(i8), optional,      intent(inout)     :: xin(0:63)
    integer(i8), optional,      intent(inout)     :: FlagIn
    real(DP),    optional,      intent(inout)     :: y2In

    integer(i8)        :: wlen, r, s, a, b, c, d

    integer(i8), save  :: w
    integer(i8), save  :: x(0:63)                  ! x(0) ... x(r-1)
    integer(i8)        :: weyl = 7046029254386353131_i8
    integer(i8)        :: t,v
    integer(i8), save  :: i = -1_i8                ! i<0 indicates first call
    integer(i8)        :: k

    real(DP)            :: t53 = 1.0_DP/9007199254740992.0_DP     ! = 0.5^53 = 1/2^53

    real(DP)            :: rn1, rn2                 ! uniform random numbers
    real(DP)            :: x1,x2,y1,ww              ! for Box-Mueller transform
    real(DP), save      :: y2
    integer(i8), save   :: Flag = 1_i8             ! if Flag = 1 return y1 else return y2

!$omp   threadprivate(x,i,w,y2,flag) 

    ! produces a 53bit Integer Random Number (0...9 007 199 254 740 992) and
    ! scales it afterwards to (0.0,1.0)
    ! transform using polar-method

    if ( present(iin) .and. seed .eq. 0_i8 )    i    = iin
    if ( present(win) .and. seed .eq. 0_i8 )    w    = win
    if ( present(xin) .and. seed .eq. 0_i8 )    x    = xin
    if ( present(Flagin) .and. seed .eq. 0_i8 ) Flag = Flagin
    if ( present(y2in) .and. seed .eq. 0_i8 )   y2   = y2in

    wlen = 64_i8
    r = 64_i8
    s = 53_i8
    a = 33_i8
    b = 26_i8
    c = 27_i8
    d = 29_i8

    If ((i .lt. 0_i8) .or. (seed .ne. 0_i8)) then     ! Initialization necessary
       If (seed .ne. 0) then                   ! v must be nonzero
          v = seed
       else
          v = NOT(seed)
       end if

       do k=wlen,1,-1                          ! Avoid correlations for close seeds
          ! This recurrence has period of 2^64-1
          v = IEOR(v,ISHFT(v,7))
          v = IEOR(v,ISHFT(v,-9))
       end do

       ! Initialize circular array
       w = v
       do k=0,r-1
          w = w + weyl
          v = IEOR(v,ISHFT(v,7))
          v = IEOR(v,ISHFT(v,-9))
          x(k) = v + w
       end do

       ! Discard first 4*r results (Gimeno)
       i = r-1
       do k = 4*r,1,-1
          i = IAND(i+1,r-1)
          t = x(i)
          v = x(IAND(i+(r-s),r-1))
          t = IEOR(t,ISHFT(t,a))
          t = IEOR(t,ISHFT(t,-b))
          v = IEOR(v,ISHFT(v,c))
          v = IEOR(v,IEOR(t,ISHFT(v,-d)))
          x(i) = v
       end do
       Flag = 1_i8
    end if ! end of initialization

    If (Flag .eq. 1_i8) then
    ! Polar method of Box-Mueller-transform to generate Gaussian distributed random number
    ww = 1.0_DP
    do while (ww .ge. 1.0_DP)

       ! Apart from initialization (above), this is the generator
       v = 0_i8
       Do While (v .eq. 0_i8)
          i = IAND(i+1,r-1)
          t = x(i)
          v = x(IAND(i+(r-s),r-1))
          t = IEOR(t,ISHFT(t,a))
          t = IEOR(t,ISHFT(t,-b))
          v = IEOR(v,ISHFT(v,c))
          v = IEOR(v,IEOR(t,ISHFT(v,-d)))
          x(i) = v
          w = w + weyl
          v = v + w
          v = ISHFT(v,-11)
       End Do

       rn1 = t53*v

       v = 0_i8
       Do While (v .eq. 0_i8)
          i = IAND(i+1,r-1)
          t = x(i)
          v = x(IAND(i+(r-s),r-1))
          t = IEOR(t,ISHFT(t,a))
          t = IEOR(t,ISHFT(t,-b))
          v = IEOR(v,ISHFT(v,c))
          v = IEOR(v,IEOR(t,ISHFT(v,-d)))
          x(i) = v
          w = w + weyl
          v = v + w
          v = ISHFT(v,-11)
       End Do

       rn2 = t53*v

       x1 = 2.0_DP * rn1 -1.0_DP
       x2 = 2.0_DP * rn2 -1.0_DP
       ww = x1*x1 + x2*x2
    end do ! end of polar method

    ww = Sqrt( (-2.0_DP * Log(ww)) / ww)
    y1 = x1 * ww
    y2 = x2 * ww

    end if ! Only if Flag = 1

    If (Flag .eq. 1) then
        Flag = 2
        DoubleRealRN = y1
    else
        Flag = 1
        DoubleRealRN = y2
    end if

    If ( present(iin) )    iin=i
    If ( present(win) )    win=w
    If ( present(xin) )    xin=x
    If ( present(Flagin) ) Flagin=Flag
    If ( present(y2in) )   y2in=y2

end subroutine xor4096gd_0d

!******************************************************************************************

subroutine xor4096gd_1d(seed,DoubleRealRN,iin,win,xin,FlagIn,y2In)

    implicit none

    integer(i8), dimension(:),                          intent(in)        :: seed
    real(DP),    dimension(size(seed)),                 intent(out)       :: DoubleRealRN
    integer(i8), dimension(size(seed)),       optional, intent(inout)     :: iin
    integer(i8), dimension(size(seed)),       optional, intent(inout)     :: win
    integer(i8), dimension(size(seed),0:63),  optional, intent(inout)     :: xin
    integer(i8), dimension(size(seed)),       optional, intent(inout)     :: FlagIn
    real(DP),    dimension(size(seed)),       optional, intent(inout)     :: y2In

    integer(i4)                         :: m
    integer(i8)                         :: wlen, r, s, a, b, c, d
    integer(i8)                         :: weyl =  7046029254386353131_i8              !Z'61C88647' = Hexadecimal notation
    integer(i8)                         :: k, j
    real(DP)                            :: t53 = 1.0_DP/9007199254740992.0_DP      ! = 0.5^24 = 1/2^24
    integer(i8), dimension(size(seed))  :: t,v
    integer(i8), dimension(:,:), allocatable, save  :: x                   ! x(0) ... x(r-1)
    integer(i8), dimension(:),   allocatable, save  :: i,w                   ! i<0 indicates first call

    real(DP),    dimension(size(seed))              :: rn1, rn2     ! uniform random numbers
    real(DP),    dimension(size(seed))              :: x1,x2,y1,ww  ! for Box-Mueller transform
    real(DP),    dimension(:),   allocatable, save  :: y2
    integer(i8), dimension(:),   allocatable, save  :: Flag         ! if Flag = 1 return y1 else return y2

!$omp   threadprivate(x,i,w,y2,flag) 

    ! produces a 53bit Integer Random Number (0...9 007 199 254 740 992) and
    ! scales it afterwards to (0.0,1.0)
    ! transform using polar-method

    m= size(seed)

    if ( present(iin) .and. (Any(seed .eq. 0)) )    i    = iin
    if ( present(win) .and. (Any(seed .eq. 0)) )    w    = win
    if ( present(xin) .and. (Any(seed .eq. 0)) )    x    = xin
    if ( present(Flagin) .and. (Any(seed .eq. 0)) ) Flag = Flagin
    if ( present(y2in) .and. (Any(seed .eq. 0)) )   y2   = y2in

    wlen = 64_i8
    r = 64_i8
    s = 53_i8
    a = 33_i8
    b = 26_i8
    c = 27_i8
    d = 29_i8

    if (.not. allocated(i)) then
       allocate(i(m))
       i = -1_i8
    end if
    if (.not. allocated(x)) then
       allocate(x(m,0:63))
    end if
    if (.not. allocated(w)) then
       allocate(w(m))
    end if
    if (.not. allocated(Flag)) then
       allocate(Flag(m))
       Flag = 1
    end if
    if (.not. allocated(y2)) then
       allocate(y2(m))
    end if

    Do j=1,m !Loop over every stream
       If ((i(j) .lt. 0) .or. (seed(j) .ne. 0)) then     ! Initialization necessary
          If (seed(j) .ne. 0) then                   ! v must be nonzero
             v(j) = seed(j)
          else
             v(j) = NOT(seed(j))
          end if

          do k=wlen,1,-1                          ! Avoid correlations for close seeds
             ! This recurrence has period of 2^64-1
             v(j) = IEOR(v(j),ISHFT(v(j),7))
             v(j) = IEOR(v(j),ISHFT(v(j),-9))
          end do

          ! Initialize circular array
          w(j) = v(j)
          do k=0,r-1
             w(j) = w(j) + weyl
             v(j) = IEOR(v(j),ISHFT(v(j),7))
             v(j) = IEOR(v(j),ISHFT(v(j),-9))
             x(j,k) = v(j) + w(j)
          end do

          ! Discard first 4*r results (Gimeno)
          i(j) = r-1
          do k = 4*r,1,-1
             i(j) = IAND(i(j)+1,r-1)
             t(j) = x(j,i(j))
             v(j) = x(j,IAND(i(j)+(r-s),r-1))
             t(j) = IEOR(t(j),ISHFT(t(j),a))
             t(j) = IEOR(t(j),ISHFT(t(j),-b))
             v(j) = IEOR(v(j),ISHFT(v(j),c))
             v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
             x(j,i(j)) = v(j)
          end do
          Flag(j) = 1
       end if ! end of initialization
    end do

    Do j=1,m        !Loop over every stream
        If (Flag(j) .eq. 1) then
        ! Polar method of Box-Mueller-transform to generate Gaussian distributed random number
            ww(j) = 1.0_DP
            do while (ww(j) .ge. 1.0_DP)

                ! Apart from initialization (above), this is the generator
                v(j) = 0_i8
                Do While (v(j) .eq. 0_i8)
                    i(j) = IAND(i(j)+1,r-1)
                    t(j) = x(j,i(j))
                    v(j) = x(j,IAND(i(j)+(r-s),r-1))
                    t(j) = IEOR(t(j),ISHFT(t(j),a))
                    t(j) = IEOR(t(j),ISHFT(t(j),-b))
                    v(j) = IEOR(v(j),ISHFT(v(j),c))
                    v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
                    x(j,i(j)) = v(j)
                    w(j) = w(j) + weyl
                    v(j) = v(j) + w(j)
                    v(j) = ISHFT(v(j),-11)
                End Do

                rn1(j) = t53*v(j)

                v(j) = 0_i8
                Do While (v(j) .eq. 0_i8)
                    i(j) = IAND(i(j)+1,r-1)
                    t(j) = x(j,i(j))
                    v(j) = x(j,IAND(i(j)+(r-s),r-1))
                    t(j) = IEOR(t(j),ISHFT(t(j),a))
                    t(j) = IEOR(t(j),ISHFT(t(j),-b))
                    v(j) = IEOR(v(j),ISHFT(v(j),c))
                    v(j) = IEOR(v(j),IEOR(t(j),ISHFT(v(j),-d)))
                    x(j,i(j)) = v(j)
                    w(j) = w(j) + weyl
                    v(j) = v(j) + w(j)
                    v(j) = ISHFT(v(j),-11)
                End Do

                rn2(j) = t53*v(j)

                x1(j) = 2.0_DP * rn1(j) -1.0_DP
                x2(j) = 2.0_DP * rn2(j) -1.0_DP
                ww(j) = x1(j)*x1(j) + x2(j)*x2(j)
            end do ! end of polar method

            ww(j) = Sqrt( (-2.0_DP * Log(ww(j))) / ww(j))
            y1(j) = x1(j) * ww(j)
            y2(j) = x2(j) * ww(j)

        end if  ! Only if Flag = 1
    end do ! Loop over each stream


    Do j=1,m
        If (Flag(j) .eq. 1) then
            Flag(j) = 2
            DoubleRealRN(j) = y1(j)
        else
            Flag(j) = 1
            DoubleRealRN(j) = y2(j)
        end if
    End do

    If ( present(iin) )    iin=i
    If ( present(win) )    win=w
    If ( present(xin) )    xin=x
    If ( present(Flagin) ) Flagin=Flag
    If ( present(y2in) )   y2in=y2

end subroutine xor4096gd_1d

end module mo_xor4096

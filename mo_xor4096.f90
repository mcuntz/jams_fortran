module mo_xor4096

  ! Written  Juliane Mai, Nov 2011

  use mo_kind, only: i4, i8, SP, DP

  Implicit NONE

  PRIVATE

  PUBLIC :: xor4096         ! Generates uniform distributed random number
  PUBLIC :: xor4096g        ! Generates gaussian distributed random number

  ! Interfaces for single and double precision routines
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
  !         xor4096

  !     PURPOSE
  !         Generates a uniform distributed random number based on xor4096 algorithm proposed by 
  !         Brent et.al (2006).
  !         ****************************************************************************************
  !         The original version of this source code is under GNU General Public Licence
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
  !         What you get is:
  !         1st call
  !             call xor( (/ 1_SP, 2_SP, 3_SP /), RN )   --> e.g. 10, 20, 30
  !         2nd call
  !             call xor( (/ 0_SP, 0_SP, 0_SP /), RN )   --> e.g. 20, 10, 50
  !         3rd call
  !             call xor( (/ 0_SP, 0_SP, 0_SP /), RN )   --> e.g. 30, 70, 10
  !
  !         Since after every initialization one looses the old stream of random numbers,
  !         it might be necessary to switch between different streams. Therefore, one needs all
  !         of the three optional arguments.
  !         What you get is:
  !         1st call of 1st stream
  !                                   call xor( 1_SP, RN, opt1_1, opt2_1, opt3_1 )   --> e.g. 10
  !         2nd call of 1st stream
  !                                   call xor( 0_SP, RN, opt1_1, opt2_1, opt3_1 )   --> e.g. 20
  !         1st call of 2nd stream
  !                                   call xor( 2_SP, RN, opt1_2, opt2_2, opt3_2 )   --> e.g. 20
  !         3rd call of 1st stream
  !                                   call xor( 0_SP, RN, opt1_1, opt2_1, opt3_1 )   --> e.g. 30
  !         Note: If you would have called 4 times without optional arguments, 
  !               you would have get 10 in the 4th call.

  !     CALLING SEQUENCE
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
  !         L'Ecuyer P & Simard R - ACM: TestU01: A C Library for Empirical Testing of 
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
    integer(i4), dimension(size(seed)), intent(out)    :: SingleIntegerRN
    integer(i4), optional, dimension(size(seed)),       intent(inout) :: iin
    integer(i4), optional, dimension(size(seed)),       intent(inout) :: win
    integer(i4), optional, dimension(size(seed),0:127), intent(inout) :: xin

    integer(i4)                        :: m
    integer(i4)                        :: wlen, r, s, a, b, c, d
    integer(i4)                        :: weyl = 1640531527_i4    !Z'61C88647'       ! Hexadecimal notation
    integer(i4)                        :: k, j
    integer(i4), dimension(size(seed)) :: t,v
    integer(i4), dimension(:,:), allocatable, save   :: x               ! x(0) ... x(r-1)
    integer(i4), dimension(:),   allocatable, save   :: i,w             ! i<0 indicates first call

    if ( present(iin) .and. (Any(seed .eq. 0)) ) i = iin
    if ( present(win) .and. (Any(seed .eq. 0)) ) w = win
    if ( present(xin) .and. (Any(seed .eq. 0)) ) x = xin

    m = size(seed)
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

    if ( present(iin) .and. (seed .eq. 0) ) i = iin
    if ( present(win) .and. (seed .eq. 0) ) w = win
    if ( present(xin) .and. (seed .eq. 0) ) x = xin

    wlen = 64
    r = 64
    s = 53
    a = 33
    b = 26
    c = 27
    d = 29

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

  subroutine xor4096l_1d(seed,DoubleIntegerRN,iin,win,xin)

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

    if ( present(iin) .and. (Any(seed .eq. 0)) ) i = iin
    if ( present(win) .and. (Any(seed .eq. 0)) ) w = win
    if ( present(xin) .and. (Any(seed .eq. 0)) ) x = xin

    m = size(seed)
    wlen = 64
    r = 64
    s = 53
    a = 33
    b = 26
    c = 27
    d = 29

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

    ! produces a 53bit Integer Random Number (0...9 007 199 254 740 992) and
    ! scales it afterwards to (0.0,1.0)

    if ( present(iin) .and. (seed .eq. 0) ) i = iin
    if ( present(win) .and. (seed .eq. 0) ) w = win
    if ( present(xin) .and. (seed .eq. 0) ) x = xin

    wlen = 64
    r = 64
    s = 53
    a = 33
    b = 26
    c = 27
    d = 29

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

    wlen = 64
    r = 64
    s = 53
    a = 33
    b = 26
    c = 27
    d = 29

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


  !
  ! ******************************************************************************************
  ! CALL XOR4096G (SEED,RN)
  ! ******************************************************************************************
  !      generates GAUSSIAN distributed Random Numbers
  !      seed is (scalar or vector) and (Single Integer or Double Integer)
  !      RN   is (scalar or vector) and (Single Real or Double Real)
  !
  !      NOTE:
  !             Every XOR4096 has to be called once with a non-zero seed for initialization and
  !             afterwards every time with Seed = 0 to guarantee that the random numbers are
  !             from the same stream.
  !             If several independent streams are necessary use VECTOR versions
  !             to initialize several streams at once.
  ! INTERFACE for subroutines:
  !     xor4096g_0d  ... 1 single precision real RN gaussian distributed N~[0,1] with one seed
  !     xor4096g_1d  ... N single precision real RN gaussian distributed N~[0,1] with n seeds
  !     xor4096gd_0d ... 1 double precision real RN gaussian distributed N~[0,1] with one seed
  !     xor4096gd_1d ... N double precision real RN gaussian distributed N~[0,1] with n seeds

subroutine xor4096gf_0d(seed,SingleRealRN,iIn,wIn,xIn,FlagIn,y2In)

    implicit none

    integer(i4),                intent(in)        :: seed
    real(SP),                    intent(out)       :: SingleRealRN
    integer(i4), optional,      intent(inout)     :: iin
    integer(i4), optional,      intent(inout)     :: win
    integer(i4), optional,      intent(inout)     :: xin(0:127)
    integer(i4), optional,      intent(inout)     :: FlagIn
    real(SP),     optional,      intent(inout)     :: y2In

    integer(i4)        :: wlen, r, s, a, b, c, d
    integer(i4), save  :: w
    integer(i4), save  :: x(0:127)                 ! x(0) ... x(r-1)
    integer(i4)        :: weyl = 1640531527_i4    !Z'61C88647'       ! Hexadecimal notation
    integer(i4)        :: t,v
    integer(i4), save  :: i = -1                   ! i<0 indicates first call
    integer(i4)        :: k
    real(SP)            :: t24 = 1.0_SP/16777216.0_SP     ! = 0.5^24 = 1/2^24

    real(SP)            :: rn1, rn2               ! uniform random numbers
    real(SP)            :: x1,x2,y1,ww            ! for Box-Mueller transform
    integer(i4),save   :: Flag = 1               ! if Flag = 1 return y1 else return y2
    real(SP),save       :: y2

    ! produces a 24bit Integer Random Number (0...16777216) and
    ! scales it afterwards to (0.0,1.0)

    ! SEED = 400
    !xor4096gf = -0.712919831275940
    !xor4096gf = -0.891133427619934
    !xor4096gf =  0.449485123157501
    !xor4096gf =  0.488143056631088
    !xor4096gf = -0.329735398292542


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
        !print '(A14,1X,F10.7,A9,B32.32)', 'xor4096gf_0d = ',y1,'   (bin) ',y1
        SingleRealRN = y1
    else
        Flag = 1
        !print '(A14,1X,F10.7,A9,B32.32)', 'xor4096gf_0d = ',y2,'   (bin) ',y1
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
    real(SP),     dimension(size(seed)),                 intent(out)       :: SingleRealRN
    integer(i4), dimension(size(seed)),       optional, intent(inout)     :: iin
    integer(i4), dimension(size(seed)),       optional, intent(inout)     :: win
    integer(i4), dimension(size(seed),0:127), optional, intent(inout)     :: xin
    integer(i4), dimension(size(seed)),       optional, intent(inout)     :: FlagIn
    real(SP),     dimension(size(seed)),       optional, intent(inout)     :: y2In

    integer(i4)                         :: m
    integer(i4)                         :: wlen, r, s, a, b, c, d
    integer(i4)                         :: weyl =  1640531527_i4              !Z'61C88647' = Hexadecimal notation
    integer(i4)                         :: k, j
    real(SP)                             :: t24 = 1.0_SP/16777216.0_SP      ! = 0.5^24 = 1/2^24
    integer(i4), dimension(size(seed))  :: t,v
    integer(i4), dimension(:,:), allocatable, save  :: x                   ! x(0) ... x(r-1)
    integer(i4), dimension(:),   allocatable, save  :: i,w                 ! i<0 indicates first call

    real(SP),     dimension(size(seed))              :: rn1, rn2     ! uniform random numbers
    real(SP),     dimension(size(seed))              :: x1,x2,y1,ww  ! for Box-Mueller transform
    real(SP),     dimension(:), allocatable, save    :: y2
    integer(i4), dimension(:), allocatable, save    :: Flag         ! if Flag = 1 return y1 else return y2

    m= size(seed)

    if ( present(iin) .and. (Any(seed .eq. 0)) )    i    = iin
    if ( present(win) .and. (Any(seed .eq. 0)) )    w    = win
    if ( present(xin) .and. (Any(seed .eq. 0)) )    x    = xin
    if ( present(Flagin) .and. (Any(seed .eq. 0)) ) Flag = Flagin
    if ( present(y2in) .and. (Any(seed .eq. 0)) )   y2   = y2in

    ! produces a 24bit Integer Random Number (0...16777216) and
    ! scales it afterwards to (0.0,1.0)

    ! SEED = 400
    !xor4096gf = -0.712919831275940
    !xor4096gf = -0.891133427619934
    !xor4096gf =  0.449485123157501
    !xor4096gf =  0.488143056631088
    !xor4096gf = -0.329735398292542


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
        !print '(A13,I6,A4)', 'xor4096gf_1d(',i(1)-5+Flag,') = '
        If (Flag(j) .eq. 1) then
            Flag(j) = 2
            !print '(A11,I4,A4,F10.7,A9,B32.32)', '   stream #',j,'  = ',y1(j),'   (bin) ',y1(j)
            SingleRealRN(j) = y1(j)
        else
            Flag(j) = 1
            !print '(A11,I4,A4,F10.7,A9,B32.32)', '   stream #',j,'  = ',y2(j),'   (bin) ',y2(j)
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
    real(DP),                    intent(out)       :: DoubleRealRN
    integer(i8), optional,      intent(inout)     :: iin
    integer(i8), optional,      intent(inout)     :: win
    integer(i8), optional,      intent(inout)     :: xin(0:63)
    integer(i8), optional,      intent(inout)     :: FlagIn
    real(DP),     optional,      intent(inout)     :: y2In

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
    integer(i8),save   :: Flag = 1_i8             ! if Flag = 1 return y1 else return y2

    ! produces a 53bit Integer Random Number (0...9 007 199 254 740 992) and
    ! scales it afterwards to (0.0,1.0)

    !Flag = 1_i8

    if ( present(iin) .and. seed .eq. 0_i8 )    i    = iin
    if ( present(win) .and. seed .eq. 0_i8 )    w    = win
    if ( present(xin) .and. seed .eq. 0_i8 )    x    = xin
    if ( present(Flagin) .and. seed .eq. 0_i8 ) Flag = Flagin
    if ( present(y2in) .and. seed .eq. 0_i8 )   y2   = y2in

    ! SEED = 400
    !xor4096gd = -0.034398645363107
    !xor4096gd =  1.870749670146693
    !xor4096gd = -1.673243084001851
    !xor4096gd =  0.071335181927899
    !xor4096gd = -0.824476113514766

    wlen = 64
    r = 64
    s = 53
    a = 33
    b = 26
    c = 27
    d = 29

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
        !print '(A14,1X,F18.15,A9,B64.64)', 'xor4096gd_0d = ',y1
        DoubleRealRN = y1
    else
        Flag = 1
        !print '(A14,1X,F18.15,A9,B64.64)', 'xor4096gd_0d = ',y2
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
    real(DP),     dimension(size(seed)),                 intent(out)       :: DoubleRealRN
    integer(i8), dimension(size(seed)),       optional, intent(inout)     :: iin
    integer(i8), dimension(size(seed)),       optional, intent(inout)     :: win
    integer(i8), dimension(size(seed),0:63),  optional, intent(inout)     :: xin
    integer(i8), dimension(size(seed)),       optional, intent(inout)     :: FlagIn
    real(DP),     dimension(size(seed)),       optional, intent(inout)     :: y2In

    integer(i4)                         :: m
    integer(i8)                         :: wlen, r, s, a, b, c, d
    integer(i8)                         :: weyl =  7046029254386353131_i8              !Z'61C88647' = Hexadecimal notation
    integer(i8)                         :: k, j
    real(DP)                             :: t53 = 1.0_DP/9007199254740992.0_DP      ! = 0.5^24 = 1/2^24
    integer(i8), dimension(size(seed))  :: t,v
    integer(i8), dimension(:,:), allocatable, save  :: x                   ! x(0) ... x(r-1)
    integer(i8), dimension(:),   allocatable, save  :: i,w                   ! i<0 indicates first call

    real(DP),     dimension(size(seed))              :: rn1, rn2     ! uniform random numbers
    real(DP),     dimension(size(seed))              :: x1,x2,y1,ww  ! for Box-Mueller transform
    real(DP),     dimension(:),   allocatable, save  :: y2
    integer(i8), dimension(:),   allocatable, save  :: Flag         ! if Flag = 1 return y1 else return y2

    ! produces a 53bit Integer Random Number (0...9 007 199 254 740 992) and
    ! scales it afterwards to (0.0,1.0)

    m= size(seed)

    if ( present(iin) .and. (Any(seed .eq. 0)) )    i    = iin
    if ( present(win) .and. (Any(seed .eq. 0)) )    w    = win
    if ( present(xin) .and. (Any(seed .eq. 0)) )    x    = xin
    if ( present(Flagin) .and. (Any(seed .eq. 0)) ) Flag = Flagin
    if ( present(y2in) .and. (Any(seed .eq. 0)) )   y2   = y2in

    ! SEED = 400
    !xor4096gd = -0.034398645363107
    !xor4096gd =  1.870749670146693
    !xor4096gd = -1.673243084001851
    !xor4096gd =  0.071335181927899
    !xor4096gd = -0.824476113514766

    wlen = 64
    r = 64
    s = 53
    a = 33
    b = 26
    c = 27
    d = 29

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
                !print*, 'ww = ',ww(j)
            end do ! end of polar method

            ww(j) = Sqrt( (-2.0_DP * Log(ww(j))) / ww(j))
            y1(j) = x1(j) * ww(j)
            y2(j) = x2(j) * ww(j)

        end if  ! Only if Flag = 1
    end do ! Loop over each stream


    Do j=1,m
        !print '(A13,I6,A4)', 'xor4096gf_1d(',i(1)-5+Flag,') = '
        If (Flag(j) .eq. 1) then
            Flag(j) = 2
            !print '(A11,I4,A4,F18.15,A9,B64.64)', '   stream #',j,'  = ',y1(j),'   (bin) ',y1(j)
            DoubleRealRN(j) = y1(j)
        else
            Flag(j) = 1
            !print '(A11,I4,A4,F18.15,A9,B64.64)', '   stream #',j,'  = ',y2(j),'   (bin) ',y2(j)
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

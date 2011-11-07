module mo_xor4096

  ! ******************************************************************************************
  ! The original version of this source code is under GNU General Public Licence
  !        xorgens.c
  ! Author: Richard P. Brent (random@rpbrent.co.uk)
  ! ******************************************************************************************
    
  ! ******************************************************************************************
  ! CALL XOR4096 (SEED,RN)
  ! ******************************************************************************************
  !      generates UNIFORM distributed Random Numbers
  !      seed is (scalar or vector) and (Single Integer or Double Integer)
  !      RN   is (scalar or vector) and (Single Integer/Real or Double Integer/Real)
  !
  !      NOTE:
  !             Every XOR4096 has to be called once with a non-zero seed for initialization and
  !             afterwards every time with Seed = 0 to guarantee that the random numbers are
  !             from the same stream.
  !             If several independent streams are necessary use VECTOR versions
  !             to initialize several streams at once.
  ! INTERFACE for subroutines:
  !     xor4096s_0d ... 1 single precision integer RN uniform distributed [-2^31,2^31-1] with one seed
  !     xor4096s_1d ... N single precision integer RN uniform distributed [-2^31,2^31-1] with n seeds
  !     xor4096f_0d ... 1 single precision real    RN uniform distributed (0.0,1.0)      with one seed
  !     xor4096f_1d ... N single precision real    RN uniform distributed (0.0,1.0)      with n seeds
  !     xor4096l_0d ... 1 double precision integer RN uniform distributed [-2^63,2^63-1] with one seed
  !     xor4096l_1d ... N double precision integer RN uniform distributed [-2^63,2^63-1] with n seeds
  !     xor4096d_0d ... 1 double precision real    RN uniform distributed (0.0,1.0)      with one seed
  !     xor4096d_1d ... N double precision real    RN uniform distributed (0.0,1.0)      with n seeds
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
  
  use mo_nrtype, only: I4B, I8B, SP, DP

  Implicit NONE

  PRIVATE

  PUBLIC :: xor4096
  PUBLIC :: xor4096g

  INTERFACE xor4096
     MODULE PROCEDURE   xor4096s_0d, xor4096s_1d, xor4096f_0d, xor4096f_1d, &
                        xor4096l_0d, xor4096l_1d, xor4096d_0d, xor4096d_1d
  END INTERFACE xor4096

  INTERFACE xor4096g
     MODULE PROCEDURE xor4096gf_0d, xor4096gf_1d, xor4096gd_0d, xor4096gd_1d
  END INTERFACE xor4096g

CONTAINS

  subroutine xor4096s_0d(seed,SingleIntegerRN,iin,win,xin)
    implicit none

    integer(I4B), intent(in)              :: seed
    integer(I4B), intent(out)             :: SingleIntegerRN
    integer(I4B), optional, intent(inout) :: iin
    integer(I4B), optional, intent(inout) :: win
    integer(I4B), optional, intent(inout) :: xin(0:127)

    integer(I4B)        :: wlen, r, s, a, b, c, d

    integer(I4B), save  :: w
    integer(I4B), save  :: x(0:127)                 ! x(0) ... x(r-1)
    integer(I4B)        :: weyl = 1640531527_I4B    !Z'61C88647'       ! Hexadecimal notation
    integer(I4B)        :: t,v
    integer(I4B), save  :: i = -1                   ! i<0 indicates first call
    integer(I4B)        :: k

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
    !if ( present(iin) .and. (seed .eq. 0) ) print*, 'iin = ',iin
    !if ( present(win) .and. (seed .eq. 0) ) print*, 'win = ',win
    !if ( present(xin) .and. (seed .eq. 0) ) print*, 'xin = ',xin(2)

    i = IAND(i+1,r-1)
    t = x(i)
    v = x(IAND(i+(r-s),r-1))
    t = IEOR(t,ISHFT(t,a))
    t = IEOR(t,ISHFT(t,-b))
    v = IEOR(v,ISHFT(v,c))
    v = IEOR(v,IEOR(t,ISHFT(v,-d)))
    x(i) = v

    w = w + weyl
    
    !print '(A12,I6,A4,I12,A9,B32.32)', 'xor4096s_0d(',i,') = ',v+w,'   (bin) ',v+w

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

    integer(I4B), dimension(:),          intent(in)     :: seed
    integer(I4B), dimension(size(seed)), intent(out)    :: SingleIntegerRN
    integer(I4B), optional, dimension(size(seed)),       intent(inout) :: iin
    integer(I4B), optional, dimension(size(seed)),       intent(inout) :: win
    integer(I4B), optional, dimension(size(seed),0:127), intent(inout) :: xin

    integer(I4B)                        :: m
    integer(I4B)                        :: wlen, r, s, a, b, c, d
    integer(I4B)                        :: weyl = 1640531527_I4B    !Z'61C88647'       ! Hexadecimal notation
    integer(I4B)                        :: k, j
    integer(I4B), dimension(size(seed)) :: t,v
    integer(I4B), dimension(:,:), allocatable, save   :: x               ! x(0) ... x(r-1)
    integer(I4B), dimension(:),   allocatable, save   :: i,w             ! i<0 indicates first call

    if ( present(iin) .and. (Any(seed .eq. 0)) ) i = iin
    if ( present(win) .and. (Any(seed .eq. 0)) ) w = win
    if ( present(xin) .and. (Any(seed .eq. 0)) ) x = xin

    m = size(seed)
    if (.not. allocated(i)) then
       allocate(i(m))
       i = -1
    else
       !write(*,*) 'I am here'
       !i = 0
       !if (present(iin)) i = iin
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

    ! SEED = 3
    !xor4096s =   -1381736336   (bin) 10101101101001000110000001110000
    !xor4096s =    1670483215   (bin) 01100011100100011000110100001111
    !xor4096s =     -70622511   (bin) 11111011110010100110001011010001
    !xor4096s =    1980498182   (bin) 01110110000011000000000100000110
    !xor4096s =    -909439004   (bin) 11001001110010110000111111100100
    ! -----------------------------------------
    ! SEED = 5
    !xor4096s =   -1632616910   (bin) 10011110101100000011111000110010
    !xor4096s =      28429767   (bin) 00000001101100011100110111000111
    !xor4096s =     692184865   (bin) 00101001010000011110011100100001
    !xor4096s =    -964703497   (bin) 11000110011111111100101011110111
    !xor4096s =      94001477   (bin) 00000101100110100101100101000101

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

    !print '(A12,I6,A4)', 'xor4096s_1d(',i(1),') = '
    !do j=1,m
    !   print '(A11,I4,A4,I12,A9,B32.32)', '   stream #',j,'  = ',v(j)+w(j),'   (bin) ',v(j)+w(j)
    !end do
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

    integer(I4B),  intent(in)  :: seed
    real(SP),      intent(out) :: SingleRealRN
    integer(I4B), optional, intent(inout) :: iin
    integer(I4B), optional, intent(inout) :: win
    integer(I4B), optional, intent(inout) :: xin(0:127)

    integer(I4B)        :: wlen, r, s, a, b, c, d
    integer(I4B), save  :: w
    integer(I4B), save  :: x(0:127)                 ! x(0) ... x(r-1)
    integer(I4B)        :: weyl = 1640531527_I4B    !Z'61C88647'       ! Hexadecimal notation
    integer(I4B)        :: t,v
    integer(I4B), save  :: i = -1                   ! i<0 indicates first call
    integer(I4B)        :: k

    real(SP)            :: t24 = 1.0_SP/16777216.0_SP     ! = 0.5^24 = 1/2^24

    ! produces a 24bit Integer Random Number (0...16777216) and
    ! scales it afterwards to (0.0,1.0)

    ! Seed = 3
    !xor4096f =   0.6782894   (bin) 00111111001011011010010001100000
    !xor4096f =   0.3889397   (bin) 00111110110001110010001100011010
    !xor4096f =   0.9835569   (bin) 00111111011110111100101001100010
    !xor4096f =   0.4611207   (bin) 00111110111011000001100000000010
    !xor4096f =   0.7882547   (bin) 00111111010010011100101100001111

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
    v = 0_I4B
    Do While (v .eq. 0_I4B)
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
    
    !print '(A12,I6,A4,F10.7,A9,B32.32)', 'xor4096f_0d(',i,') = ',t24*v,'   (bin) ',t24*v

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

    integer(I4B), dimension(:),          intent(in)  :: seed
    real(SP),     dimension(size(seed)), intent(out) :: SingleRealRN
    integer(I4B), optional, dimension(size(seed)),       intent(inout) :: iin
    integer(I4B), optional, dimension(size(seed)),       intent(inout) :: win
    integer(I4B), optional, dimension(size(seed),0:127), intent(inout) :: xin

    integer(I4B)                         :: m
    integer(I4B)                         :: wlen, r, s, a, b, c, d
    integer(I4B)                         :: weyl =  1640531527_I4B              !Z'61C88647' = Hexadecimal notation
    integer(I4B)                         :: k, j
    real(SP), save                       :: t24 = 1.0_SP/16777216.0_SP      ! = 0.5^24 = 1/2^24
    integer(I4B), dimension(size(seed))  :: t,v
    integer(I4B), dimension(:,:), allocatable, save  :: x                   ! x(0) ... x(r-1)
    integer(I4B), dimension(:),   allocatable, save  :: i,w                 ! i<0 indicates first call

    m= size(seed)

    if ( present(iin) .and. (Any(seed .eq. 0)) ) i = iin
    if ( present(win) .and. (Any(seed .eq. 0)) ) w = win
    if ( present(xin) .and. (Any(seed .eq. 0)) ) x = xin

    ! produces a 24bit Integer Random Number (0...16777216) and
    ! scales it afterwards to (0.0,1.0)

    ! Seed = 3
    !xor4096f =   0.6782894   (bin) 00111111001011011010010001100000
    !xor4096f =   0.3889397   (bin) 00111110110001110010001100011010
    !xor4096f =   0.9835569   (bin) 00111111011110111100101001100010
    !xor4096f =   0.4611207   (bin) 00111110111011000001100000000010
    !xor4096f =   0.7882547   (bin) 00111111010010011100101100001111

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
    v = 0_I4B
    Do j=1,m
       Do While (v(j) .eq. 0_I4B)
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

    !print '(A12,I6,A4)', 'xor4096f_1d(',i(1),') = '
    !do j=1,m
    !   print '(A11,I4,A4,F10.7,A9,B32.32)', '   stream #',j,'  = ',t24*v(j),'   (bin) ',t24*v(j)
    !end do
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

    integer(I8B), intent(in)  :: seed
    integer(I8B), intent(out) :: DoubleIntegerRN
    integer(I8B), optional, intent(inout) :: iin
    integer(I8B), optional, intent(inout) :: win
    integer(I8B), optional, intent(inout) :: xin(0:63)

    integer(I8B)        :: wlen, r, s, a, b, c, d
    integer(I8B), save  :: w
    integer(I8B), save  :: x(0:63)                  ! x(0) ... x(r-1)
    integer(I8B)        :: weyl = 7046029254386353131_I8B
    integer(I8B)        :: t,v
    integer(I8B), save  :: i = -1                   ! i<0 indicates first call
    integer(I8B)        :: k

    ! SEED = 8974719815813083380_I8B
    !xor4096l =    865112168962364054   (bin) 0000110000000001011111101110111100100001110111101001001010010110
    !xor4096l =   5875351260125253729   (bin) 0101000110001001011100001101011110010110011010001101000001100001
    !xor4096l =  -1631495050823355867   (bin) 1110100101011011110000111110010100110010110111011011001000100101
    !xor4096l =   7150529841878842630   (bin) 0110001100111011110010010000000110001101000001001011110100000110
    !xor4096l =   8167262367244560672   (bin) 0111000101010111111100011101100111011011000001110010100100100000
    ! -----------------------------------------
    ! SEED = 599111145707029239_I8B
    !xor4096l =   7707611802973968921   (bin) 0110101011110110111100000001100100011110110111110001101000011001
    !xor4096l =   8714808719791153023   (bin) 0111100011110001001110000101101110001100101000001101111101111111
    !xor4096l =   8607982351506415915   (bin) 0111011101110101101100100101011101110110010111110010000100101011
    !xor4096l =  -8342501239476650915   (bin) 1000110000111001011110110101001010111010101111001111110001011101
    !xor4096l =  -3070962797875769889   (bin) 1101010101100001101111111001101100100010001011001100110111011111

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
    
    !print '(A12,I6,A4,I20,A9,B64.64)', 'xor4096l_0d(',i,') = ',v+w,'   (bin) ',v+w
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

    integer(I8B), dimension(:), intent(in)  :: seed
    integer(I8B), dimension(size(seed)), intent(out) :: DoubleIntegerRN
    integer(I8B), optional, dimension(size(seed)),      intent(inout) :: iin
    integer(I8B), optional, dimension(size(seed)),      intent(inout) :: win
    integer(I8B), optional, dimension(size(seed),0:63), intent(inout) :: xin

    integer(I4B)        :: m
    integer(I8B)        :: wlen, r, s, a, b, c, d
    integer(I8B)        :: weyl = 7046029254386353131_I8B !B'0110000111001000100001100100011010000000101101011000001111101011'
    integer(I8B)        :: k, j
    integer(I8B), dimension(size(seed))              :: t,v
    integer(I8B), dimension(:,:), allocatable, save  :: x                  ! x(0) ... x(r-1)
    integer(I8B), dimension(:),   allocatable, save  :: i,w                ! i<0 indicates first call

    ! SEED = 8974719815813083380_I8B
    !xor4096l =    865112168962364054   (bin) 0000110000000001011111101110111100100001110111101001001010010110
    !xor4096l =   5875351260125253729   (bin) 0101000110001001011100001101011110010110011010001101000001100001
    !xor4096l =  -1631495050823355867   (bin) 1110100101011011110000111110010100110010110111011011001000100101
    !xor4096l =   7150529841878842630   (bin) 0110001100111011110010010000000110001101000001001011110100000110
    !xor4096l =   8167262367244560672   (bin) 0111000101010111111100011101100111011011000001110010100100100000
    ! -----------------------------------------
    ! SEED = 599111145707029239_I8B
    !xor4096l =   7707611802973968921   (bin) 0110101011110110111100000001100100011110110111110001101000011001
    !xor4096l =   8714808719791153023   (bin) 0111100011110001001110000101101110001100101000001101111101111111
    !xor4096l =   8607982351506415915   (bin) 0111011101110101101100100101011101110110010111110010000100101011
    !xor4096l =  -8342501239476650915   (bin) 1000110000111001011110110101001010111010101111001111110001011101
    !xor4096l =  -3070962797875769889   (bin) 1101010101100001101111111001101100100010001011001100110111011111

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

    !print '(A12,I6,A4)', 'xor4096l_1d(',i(1),') = '
    !do j=1,m
    !   print '(A11,I4,A4,I20,A9,B64.64)', '   stream #',j,'  = ',v(j)+w(j),'   (bin) ',v(j)+w(j)
    !end do
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

    integer(I8B), intent(in)  :: seed
    real(DP),     intent(out) :: DoubleRealRN
    integer(I8B), optional, intent(inout) :: iin
    integer(I8B), optional, intent(inout) :: win
    integer(I8B), optional, intent(inout) :: xin(0:63)

    integer(I8B)        :: wlen, r, s, a, b, c, d

    integer(I8B), save  :: w
    integer(I8B), save  :: x(0:63)                  ! x(0) ... x(r-1)
    integer(I8B)        :: weyl = 7046029254386353131_I8B
    integer(I8B)        :: t,v
    integer(I8B), save  :: i = -1                   ! i<0 indicates first call
    integer(I8B)        :: k

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
    v = 0_I8B
    Do While (v .eq. 0_I8B)
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

    !print '(A12,I6,A4,F18.15,A9,B64.64)', 'xor4096d_0d(',i,') = ',t53*v,'   (bin) ',t53*v
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

    integer(I8B), dimension(:),          intent(in)  :: seed
    real(DP),     dimension(size(seed)), intent(out) :: DoubleRealRN
    integer(I8B), optional, dimension(size(seed)),      intent(inout) :: iin
    integer(I8B), optional, dimension(size(seed)),      intent(inout) :: win
    integer(I8B), optional, dimension(size(seed),0:63), intent(inout) :: xin

    integer(I4B)                       :: m
    integer(I8B)                       :: wlen, r, s, a, b, c, d
    integer(I8B)                       :: weyl = 7046029254386353131_I8B
    real(DP)                           :: t53  = 1.0_DP/9007199254740992.0_DP  ! = 0.5^53 = 1/2^53
    integer(I8B)                       :: k,j
    integer(I8B), dimension(size(seed))              :: t,v
    integer(I8B), dimension(:,:), allocatable, save  :: x       ! x(0) ... x(r-1)
    integer(I8B), dimension(:),   allocatable, save  :: w
    integer(I8B), dimension(:),   allocatable, save  :: i       ! i<0 indicates first call             

    ! produces a 53bit Integer Random Number (0...9 007 199 254 740 992) and
    ! scales it afterwards to (0.0,1.0)

    ! SEED = 3
    !xor4096d =   0.721549534183639   (bin) 0011111111100111000101101110111100001100011110000110100010010001
    !xor4096d =   0.620277262221687   (bin) 0011111111100011110110010100111110110011011101100011100111011001
    !xor4096d =   0.486520536087663   (bin) 0011111111011111001000110010011100000111110101010000110010101100
    !xor4096d =   0.032477880286792   (bin) 0011111110100000101000001111000011010011010110011011000010000000
    !xor4096d =   0.139152986577077   (bin) 0011111111000001110011111100001111011011001111101010000001101100

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
    v = 0_I8B
    Do j=1,m
       Do While (v(j) .eq. 0_I8B)
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

    !print '(A12,I6,A4)', 'xor4096d_1d(',i(1),') = '
    !do j=1,m
    !   print '(A11,I4,A4,F18.15,A9,B64.64)', '   stream #',j,'  = ',t53*v(j),'   (bin) ',t53*v(j)
    !end do
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

subroutine xor4096gf_0d(seed,SingleRealRN,iIn,wIn,xIn,FlagIn,y2In)

    implicit none

    integer(I4B),                intent(in)        :: seed
    real(SP),                    intent(out)       :: SingleRealRN
    integer(I4B), optional,      intent(inout)     :: iin
    integer(I4B), optional,      intent(inout)     :: win
    integer(I4B), optional,      intent(inout)     :: xin(0:127)
    integer(I4B), optional,      intent(inout)     :: FlagIn
    real(SP),     optional,      intent(inout)     :: y2In

    integer(I4B)        :: wlen, r, s, a, b, c, d
    integer(I4B), save  :: w
    integer(I4B), save  :: x(0:127)                 ! x(0) ... x(r-1)
    integer(I4B)        :: weyl = 1640531527_I4B    !Z'61C88647'       ! Hexadecimal notation
    integer(I4B)        :: t,v
    integer(I4B), save  :: i = -1                   ! i<0 indicates first call
    integer(I4B)        :: k
    real(SP)            :: t24 = 1.0_SP/16777216.0_SP     ! = 0.5^24 = 1/2^24

    real(SP)            :: rn1, rn2               ! uniform random numbers
    real(SP)            :: x1,x2,y1,ww            ! for Box-Mueller transform
    integer(I4B),save   :: Flag = 1               ! if Flag = 1 return y1 else return y2
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
       v = 0_I4B
       Do While (v .eq. 0_I4B)
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

       v = 0_I4B
       Do While (v .eq. 0_I4B)
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

    integer(I4B), dimension(:),                          intent(in)        :: seed
    real(SP),     dimension(size(seed)),                 intent(out)       :: SingleRealRN
    integer(I4B), dimension(size(seed)),       optional, intent(inout)     :: iin
    integer(I4B), dimension(size(seed)),       optional, intent(inout)     :: win
    integer(I4B), dimension(size(seed),0:127), optional, intent(inout)     :: xin
    integer(I4B), dimension(size(seed)),       optional, intent(inout)     :: FlagIn
    real(SP),     dimension(size(seed)),       optional, intent(inout)     :: y2In

    integer(I4B)                         :: m
    integer(I4B)                         :: wlen, r, s, a, b, c, d
    integer(I4B)                         :: weyl =  1640531527_I4B              !Z'61C88647' = Hexadecimal notation
    integer(I4B)                         :: k, j
    real(SP)                             :: t24 = 1.0_SP/16777216.0_SP      ! = 0.5^24 = 1/2^24
    integer(I4B), dimension(size(seed))  :: t,v
    integer(I4B), dimension(:,:), allocatable, save  :: x                   ! x(0) ... x(r-1)
    integer(I4B), dimension(:),   allocatable, save  :: i,w                 ! i<0 indicates first call
    
    real(SP),     dimension(size(seed))              :: rn1, rn2     ! uniform random numbers
    real(SP),     dimension(size(seed))              :: x1,x2,y1,ww  ! for Box-Mueller transform
    real(SP),     dimension(:), allocatable, save    :: y2
    integer(I4B), dimension(:), allocatable, save    :: Flag         ! if Flag = 1 return y1 else return y2
    
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
                v(j) = 0_I4B
                Do While (v(j) .eq. 0_I4B)
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

                v(j) = 0_I4B
                Do While (v(j) .eq. 0_I4B)
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

    integer(I8B),                intent(in)        :: seed
    real(DP),                    intent(out)       :: DoubleRealRN
    integer(I8B), optional,      intent(inout)     :: iin
    integer(I8B), optional,      intent(inout)     :: win
    integer(I8B), optional,      intent(inout)     :: xin(0:63)
    integer(I8B), optional,      intent(inout)     :: FlagIn
    real(DP),     optional,      intent(inout)     :: y2In

    integer(I8B)        :: wlen, r, s, a, b, c, d

    integer(I8B), save  :: w
    integer(I8B), save  :: x(0:63)                  ! x(0) ... x(r-1)
    integer(I8B)        :: weyl = 7046029254386353131_I8B
    integer(I8B)        :: t,v
    integer(I8B), save  :: i = -1_I8B                ! i<0 indicates first call
    integer(I8B)        :: k

    real(DP)            :: t53 = 1.0_DP/9007199254740992.0_DP     ! = 0.5^53 = 1/2^53

    real(DP)            :: rn1, rn2                 ! uniform random numbers
    real(DP)            :: x1,x2,y1,ww              ! for Box-Mueller transform
    real(DP), save      :: y2
    integer(I8B),save   :: Flag = 1_I8B             ! if Flag = 1 return y1 else return y2

    ! produces a 53bit Integer Random Number (0...9 007 199 254 740 992) and
    ! scales it afterwards to (0.0,1.0)

    !Flag = 1_I8B

    if ( present(iin) .and. seed .eq. 0_I8B )    i    = iin
    if ( present(win) .and. seed .eq. 0_I8B )    w    = win
    if ( present(xin) .and. seed .eq. 0_I8B )    x    = xin
    if ( present(Flagin) .and. seed .eq. 0_I8B ) Flag = Flagin
    if ( present(y2in) .and. seed .eq. 0_I8B )   y2   = y2in

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

    If ((i .lt. 0_I8B) .or. (seed .ne. 0_I8B)) then     ! Initialization necessary
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
       Flag = 1_I8B
    end if ! end of initialization

    If (Flag .eq. 1_I8B) then
    ! Polar method of Box-Mueller-transform to generate Gaussian distributed random number
    ww = 1.0_DP
    do while (ww .ge. 1.0_DP)

       ! Apart from initialization (above), this is the generator
       v = 0_I8B
       Do While (v .eq. 0_I8B)
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

       v = 0_I8B
       Do While (v .eq. 0_I8B)
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

    integer(I8B), dimension(:),                          intent(in)        :: seed
    real(DP),     dimension(size(seed)),                 intent(out)       :: DoubleRealRN
    integer(I8B), dimension(size(seed)),       optional, intent(inout)     :: iin
    integer(I8B), dimension(size(seed)),       optional, intent(inout)     :: win
    integer(I8B), dimension(size(seed),0:63),  optional, intent(inout)     :: xin
    integer(I8B), dimension(size(seed)),       optional, intent(inout)     :: FlagIn
    real(DP),     dimension(size(seed)),       optional, intent(inout)     :: y2In

    integer(I4B)                         :: m
    integer(I8B)                         :: wlen, r, s, a, b, c, d
    integer(I8B)                         :: weyl =  7046029254386353131_I8B              !Z'61C88647' = Hexadecimal notation
    integer(I8B)                         :: k, j
    real(DP)                             :: t53 = 1.0_DP/9007199254740992.0_DP      ! = 0.5^24 = 1/2^24
    integer(I8B), dimension(size(seed))  :: t,v
    integer(I8B), dimension(:,:), allocatable, save  :: x                   ! x(0) ... x(r-1)
    integer(I8B), dimension(:),   allocatable, save  :: i,w                   ! i<0 indicates first call

    real(DP),     dimension(size(seed))              :: rn1, rn2     ! uniform random numbers
    real(DP),     dimension(size(seed))              :: x1,x2,y1,ww  ! for Box-Mueller transform
    real(DP),     dimension(:),   allocatable, save  :: y2
    integer(I8B), dimension(:),   allocatable, save  :: Flag         ! if Flag = 1 return y1 else return y2

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
       i = -1_I8B
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
                v(j) = 0_I8B
                Do While (v(j) .eq. 0_I8B)
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

                v(j) = 0_I8B
                Do While (v(j) .eq. 0_I8B)
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

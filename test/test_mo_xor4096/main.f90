!*******************************************************
!
!   TO TEST RNGs
!
!*******************************************************

program RNG

    use mo_kind,  only: i4, i8, SP, DP
    use mo_xor4096, only: xor4096, xor4096g

    implicit none

    integer(i4)                    :: ISeedSP_D0
    integer(i8)                    :: ISeedDP_D0
    integer(i4), allocatable       :: ISeedSP(:)
    integer(i8), allocatable       :: ISeedDP(:)

    integer(i4)                    :: SingleIntegerRN_D0
    integer(i4), allocatable       :: SingleIntegerRN(:)
    integer(i8)                    :: DoubleIntegerRN_D0
    integer(i8), allocatable       :: DoubleIntegerRN(:)
    real(SP)                       :: SingleRealRN_D0
    real(SP),    allocatable       :: SingleRealRN(:)
    real(DP)                       :: DoubleRealRN_D0
    real(DP),    allocatable       :: DoubleRealRN(:)

    integer(i4)                    :: NumberOfStreams = 3

    character(10)                  :: CheckString
    integer(i4), allocatable       :: CheckISP_D0(:)
    integer(i4), allocatable       :: CheckISP_D1(:,:)
    integer(i8), allocatable       :: CheckIDP_D0(:)
    integer(i8), allocatable       :: CheckIDP_D1(:,:)
    real(SP),    allocatable       :: CheckRSP_D0(:)
    real(SP),    allocatable       :: CheckRSP_D1(:,:)
    real(DP),    allocatable       :: CheckRDP_D0(:)
    real(DP),    allocatable       :: CheckRDP_D1(:,:)

    integer(i8)                    :: i,j
    
    ! Needed for optional versions: Single
    integer(i4)                    :: SingleIin_D0, SingleWin_D0
    integer(i4)                    :: SingleFlagIn_D0
    real(SP)                       :: SingleY2In_D0
    integer(i4),dimension(0:127)   :: SingleXin_D0
    ! Needed for optional versions: Double
    integer(i8)                    :: DoubleIin_D0, DoubleWin_D0
    integer(i8)                    :: DoubleFlagIn_D0
    real(DP)                       :: DoubleY2In_D0
    integer(i8),dimension(0:63)    :: DoubleXin_D0

    allocate(ISeedSP(NumberOfStreams))
    allocate(ISeedDP(NumberOfStreams))
    allocate(SingleIntegerRN(NumberOfStreams))
    allocate(DoubleIntegerRN(NumberOfStreams))
    allocate(SingleRealRN(NumberOfStreams))
    allocate(DoubleRealRN(NumberOfStreams))

    print*,'----------------------------------------'
    print*,'UNIFORM Single Precision (32 bit) Integer'
    print*,'          1 STREAM                       '
    print*,'----------------------------------------'
    ISeedSP_D0 = 3_i4
    allocate (CheckISP_D0(5))
    CheckISP_D0 = (/-1381736336,1670483215,-70622511,1980498182,-909439004/)
    do i=1,5
        call xor4096(ISeedSP_D0, SingleIntegerRN_D0)
        if (SingleIntegerRN_D0 .eq. CheckISP_D0(i)) then
            CheckString = 'Ok'
        else
            CheckString = 'Failed'
        end if
        print '(A8,I6,A4,I12,A4,A10)', 'xor4096(',i,') = ',SingleIntegerRN_D0,'    ',CheckString
        ISeedSP_D0 = 0_i4
    end do
    deallocate (CheckISP_D0)

    print*,'----------------------------------------'
    print*,'UNIFORM Single Precision (32 bit) Integer'
    print*,'          3 STREAMS                      '
    print*,'----------------------------------------'
    ISeedSP = (/ 3_i4, 5_i4, 7_i4 /)
    allocate (CheckISP_D1(5,NumberOfStreams))
    CheckISP_D1(1,:) = (/-1381736336,-1632616910,739198204/)
    CheckISP_D1(2,:) = (/1670483215,28429767,1655944655/)
    CheckISP_D1(3,:) = (/-70622511,692184865,1023508511/)
    CheckISP_D1(4,:) = (/1980498182,-964703497,290140738/)
    CheckISP_D1(5,:) = (/-909439004,94001477,-1671307991/)
    do i=1,5
        call xor4096(ISeedSP,SingleIntegerRN)
        print '(A8,I6,A4)', 'xor4096(',i,') = '
        do j=1,NumberOfStreams
            if (SingleIntegerRN(j) .eq. CheckISP_D1(i,j)) then
                CheckString = 'Ok'
            else
                CheckString = 'Failed'
            end if
            print '(A11,I4,A4,I12,A4,A10)', '   stream #',j,'  = ',SingleIntegerRN(j),'    ',CheckString
        end do
        ISeedSP = 0_i4
    end do
    deallocate (CheckISP_D1)

    print*,'-----------------------------------------'
    print*,'UNIFORM Single Precision (32 bit) Integer'
    print*,'          WITH OPTIONAL                  '
    print*,'-----------------------------------------'
    ISeedSP_D0 = 3_i4
    allocate (CheckISP_D0(6))
    CheckISP_D0 = (/-1381736336,1670483215,-70622511,1980498182,-909439004,-1632616910/)
    do i=1,3
        call xor4096(ISeedSP_D0,SingleIntegerRN_D0,SingleIin_D0,SingleWin_D0,SingleXin_D0)
        if (SingleIntegerRN_D0 .eq. CheckISP_D0(i)) then
            CheckString = 'Ok'
        else
            CheckString = 'Failed'
        end if
        print '(A8,I6,A4,I12,A4,A10)', 'xor4096(',i,') = ',SingleIntegerRN_D0,'    ',CheckString
        ISeedSP_D0 = 0_i4
    end do
    ! Start a new Stream
    call xor4096(5_i4,SingleIntegerRN_D0)
    if (SingleIntegerRN_D0 .eq. CheckISP_D0(6)) then
        CheckString = 'Ok'
    else
        CheckString = 'Failed'
    end if
    print '(A8,I6,A4,I12,A4,A10)', 'xor4096(',1,') = ',SingleIntegerRN_D0,'    ',CheckString
    ISeedSP_D0 = 0_i4
    ! Go back to the old stream
    do i=4,5
        call xor4096(ISeedSP_D0,SingleIntegerRN_D0,SingleIin_D0,SingleWin_D0,SingleXin_D0)
        if (SingleIntegerRN_D0 .eq. CheckISP_D0(i)) then
            CheckString = 'Ok'
        else
            CheckString = 'Failed'
        end if
        print '(A8,I6,A4,I12,A4,A10)', 'xor4096(',i,') = ',SingleIntegerRN_D0,'    ',CheckString
        ISeedSP_D0 = 0_i4
    end do
    deallocate (CheckISP_D0)

    deallocate(ISeedSP)
    deallocate(ISeedDP)
    deallocate(SingleIntegerRN)
    deallocate(DoubleIntegerRN)
    deallocate(SingleRealRN)
    deallocate(DoubleRealRN)


!!$ -----------------------------------------
!!$ UNIFORM Single Precision (32 bit) Integer
!!$           WITH OPTIONAL                  
!!$ -----------------------------------------
!!$xor4096s_0d(     1) =  -1381736336   (bin) 10101101101001000110000001110000
!!$xor4096s_0d(     2) =   1670483215   (bin) 01100011100100011000110100001111
!!$xor4096s_0d(     3) =    -70622511   (bin) 11111011110010100110001011010001
!!$xor4096s_0d(     1) =  -1632616910   (bin) 10011110101100000011111000110010
!!$xor4096s_0d(     4) =   1980498182   (bin) 01110110000011000000000100000110
!!$xor4096s_0d(     5) =   -909439004   (bin) 11001001110010110000111111100100
!!$ ---------------------------------------- 

end program





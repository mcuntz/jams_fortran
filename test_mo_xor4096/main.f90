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
    integer(i4),allocatable        :: ISeedSP(:)
    integer(i8),allocatable        :: ISeedDP(:)

    integer(i4)                    :: SingleIntegerRN_D0
    integer(i4),allocatable        :: SingleIntegerRN(:)
    integer(i8)                    :: DoubleIntegerRN_D0
    integer(i8),allocatable        :: DoubleIntegerRN(:)
    real(SP)                        :: SingleRealRN_D0
    real(SP),allocatable            :: SingleRealRN(:)
    real(DP)                        :: DoubleRealRN_D0
    real(DP),allocatable            :: DoubleRealRN(:)

    integer(i4)                    :: NumberOfStreams = 3

    integer(i8)                    :: i,j,k
    
    ! Needed for optinal versions: Single
    integer(i4)                    :: SingleIin_D0, SingleWin_D0
    integer(i4)                    :: SingleFlagIn_D0
    real(SP)                        :: SingleY2In_D0
    integer(i4),dimension(0:127)   :: SingleXin_D0
    ! Needed for optinal versions: Double
    integer(i8)                    :: DoubleIin_D0, DoubleWin_D0
    integer(i8)                    :: DoubleFlagIn_D0
    real(DP)                        :: DoubleY2In_D0
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
    do i=1,5
        call xor4096(ISeedSP_D0, SingleIntegerRN_D0)
        print '(A12,I6,A4,I12,A9,B32.32)', 'xor4096s_0d(',i,') = ',SingleIntegerRN_D0, &
                                           '   (bin) ',SingleIntegerRN_D0
        ISeedSP_D0 = 0_i4
    end do

    print*,'----------------------------------------'
    print*,'UNIFORM Single Precision (32 bit) Integer'
    print*,'          3 STREAMS                      '
    print*,'----------------------------------------'
    ISeedSP = (/ 3_i4, 5_i4, 7_i4 /)
    do i=1,5
        call xor4096(ISeedSP,SingleIntegerRN)
        print '(A12,I6,A4)', 'xor4096s_1d(',i,') = '
        do j=1,NumberOfStreams
            print '(A11,I4,A4,I12,A9,B32.32)', '   stream #',j,'  = ',SingleIntegerRN(j), &
                                               '   (bin) ',SingleIntegerRN(j)
        end do
        ISeedSP = 0_i4
    end do

    print*,'-----------------------------------------'
    print*,'UNIFORM Single Precision (32 bit) Integer'
    print*,'          WITH OPTIONAL                  '
    print*,'-----------------------------------------'
    ISeedSP_D0 = 3_i4
    do i=1,3
        call xor4096(ISeedSP_D0,SingleIntegerRN_D0,SingleIin_D0,SingleWin_D0,SingleXin_D0)
        print '(A12,I6,A4,I12,A9,B32.32)', 'xor4096s_0d(',i,') = ',SingleIntegerRN_D0, &
                                           '   (bin) ',SingleIntegerRN_D0
        ISeedSP_D0 = 0_i4
    end do
    ! Start a new Stream
    call xor4096(5_i4,SingleIntegerRN_D0)
    print '(A12,I6,A4,I12,A9,B32.32)', 'xor4096s_0d(',1,') = ',SingleIntegerRN_D0, &
                                       '   (bin) ',SingleIntegerRN_D0
    ISeedSP_D0 = 0_i4
    ! Go back to the old stream
    do i=4,5
        call xor4096(ISeedSP_D0,SingleIntegerRN_D0,SingleIin_D0,SingleWin_D0,SingleXin_D0)
        print '(A12,I6,A4,I12,A9,B32.32)', 'xor4096s_0d(',i,') = ',SingleIntegerRN_D0, &
                                           '   (bin) ',SingleIntegerRN_D0
        ISeedSP_D0 = 0_i4
    end do

    print*,'---------------------------------------- '
    print*,'UNIFORM Single Precision (32 bit) Real'
    print*,'          1 STREAM                       '
    print*,'---------------------------------------- '
    ISeedSP_D0 = 3_i4
    do i=1,5
        call xor4096(ISeedSP_D0, SingleRealRN_D0)
        print '(A12,I6,A4,F10.7,A9,B32.32)', 'xor4096f_0d(',i,') = ',SingleRealRN_D0, &
                                             '   (bin) ',SingleRealRN_D0
        ISeedSP_D0 = 0_i4
    end do

    print*,'-----------------------------------------'
    print*,'UNIFORM Single Precision (32 bit) Real   '
    print*,'          3 STREAMS                      '
    print*,'-----------------------------------------'
    ISeedSP = (/ 3_i4, 5_i4, 7_i4 /)
    do i=1,5
        call xor4096(ISeedSP,SingleRealRN)
        print '(A12,I6,A4)', 'xor4096f_1d(',i,') = '
        do j=1,NumberOfStreams
            print '(A11,I4,A4,F10.7,A9,B32.32)', '   stream #',j,'  = ',SingleRealRN(j), &
                                                 '   (bin) ',SingleRealRN(j)
        end do
        ISeedSP = 0_i4
    end do

    print*,'---------------------------------------- '
    print*,'UNIFORM Double Precision (64 bit) Integer'
    print*,'          1 STREAM                       '
    print*,'---------------------------------------- '
    ISeedDP_D0 = 3_i8
    do i=1,5
        call xor4096(ISeedDP_D0, DoubleIntegerRN_D0)
        print '(A12,I6,A4,I20,A9,B64.64)', 'xor4096l_0d(',i,') = ',DoubleIntegerRN_D0, &
                                           '   (bin) ',DoubleIntegerRN_D0
        ISeedDP_D0 = 0_i8
    end do

    print*,'-----------------------------------------'
    print*,'UNIFORM Double Precision (64 bit) Integer'
    print*,'          3 STREAMS                      '
    print*,'---------------------------------------- '
    ISeedDP = (/ 8974719815813083380_i8, 599111145707029239_i8 , 3_i8/)
    do i=1,5
        call xor4096(ISeedDP, DoubleIntegerRN)
        print '(A12,I6,A4)', 'xor4096l_1d(',i,') = '
        do j=1,NumberOfStreams
            print '(A11,I4,A4,I20,A9,B64.64)', '   stream #',j,'  = ',DoubleIntegerRN(j), &
                                               '   (bin) ',DoubleIntegerRN(j)
        end do
        ISeedDP = 0_i8
    end do

    print*,'----------------------------------------'
    print*,'UNIFORM Double Precision (64 bit) Real  '
    print*,'          1 STREAM                      '
    print*,'----------------------------------------'
    ISeedDP_D0 = 3_i8
    do i=1,5
        call xor4096(ISeedDP_D0, DoubleRealRN_D0)
        print '(A12,I6,A4,F18.15,A9,B64.64)', 'xor4096d_0d(',i,') = ',DoubleRealRN_D0, &
                                              '   (bin) ',DoubleRealRN_D0
        ISeedDP_D0 = 0_i8
    end do

    print*,'----------------------------------------'
    print*,'UNIFORM Double Precision (64 bit) Real  '
    print*,'          3 STREAMS                      '
    print*,'----------------------------------------'
    ISeedDP = (/ 3_i8, 8974719815813083380_i8, 599111145707029239_i8 /)
    do i=1,5
        call xor4096(ISeedDP, DoubleRealRN)
        print '(A12,I6,A4)', 'xor4096d_1d(',i,') = '
        do j=1,NumberOfStreams
            print '(A11,I4,A4,F18.15,A9,B64.64)', '   stream #',j,'  = ',DoubleRealRN(j), &
                                                  '   (bin) ',DoubleRealRN(j)
        end do
        ISeedDP = 0_i8
    end do

    print*,'-----------------------------------------'
    print*,'GAUSSIAN Single Precision (32 bit) Real  '
    print*,'          1 STREAM                       '
    print*,'-----------------------------------------'

    ISeedSP_D0 = 400_i4
    do i=1,5
        call xor4096g(ISeedSP_D0, SingleRealRN_D0)
        print '(A14,1X,F10.7,A9,B32.32)', 'xor4096gf_0d = ',SingleRealRN_D0, &
                                          '   (bin) ',SingleRealRN_D0
        ISeedSP_D0 = 0_i4
    end do

    print*,'-----------------------------------------'
    print*,'GAUSSIAN Single Precision (32 bit) Real  '
    print*,'          3 STREAMS                       '
    print*,'-----------------------------------------'

    ISeedSP = (/ 400_i4, 5_i4, 7_i4 /)
    do i=1,5
        call xor4096g(ISeedSP, SingleRealRN)
        print '(A13,I6,A4)', 'xor4096gf_1d(',i,') = '
        do j=1,NumberOfStreams
            print '(A11,I4,A4,F10.7,A9,B32.32)', '   stream #',j,'  = ',SingleRealRN(j), &
                                                 '   (bin) ',SingleRealRN(j)
        end do
        ISeedSP = 0_i4
    end do

    print*,'-----------------------------------------'
    print*,'Gaussian Single Precision (32 bit) Real  '
    print*,'          WITH OPTIONAL                  '
    print*,'-----------------------------------------'
    ISeedSP_D0 = 400_i4
    do i=1,3
        call xor4096g(ISeedSP_D0, SingleRealRN_D0, SingleIin_D0,SingleWin_D0,SingleXin_D0,SingleFlagin_D0,SingleY2In_D0)
        print '(A12,I6,A4,F10.7,A9,B32.32)', 'xor4096gf_0d(',i,') = ',SingleRealRN_D0, &
                                           '   (bin) ',SingleRealRN_D0
        ISeedSP_D0 = 0_i4
    end do
    ! Start a new Stream
    call xor4096g(5_i4, SingleRealRN_D0)
    print '(A12,I6,A4,F10.7,A9,B32.32)', 'xor4096gf_0d(',1,') = ',SingleRealRN_D0, &
                                       '   (bin) ',SingleRealRN_D0
    ISeedSP_D0 = 0_i4
    ! Go back to the old stream
    do i=4,5
        call xor4096g(ISeedSP_D0, SingleRealRN_D0, SingleIin_D0,SingleWin_D0,SingleXin_D0,SingleFlagin_D0,SingleY2In_D0)
        print '(A12,I6,A4,F10.7,A9,B32.32)', 'xor4096gf_0d(',i,') = ',SingleRealRN_D0, &
                                           '   (bin) ',SingleRealRN_D0
        ISeedSP_D0 = 0_i4
    end do

    print*,'-----------------------------------------'
    print*,'GAUSSIAN Double Precision (64 bit) Real  '
    print*,'          1 STREAM                       '
    print*,'-----------------------------------------'

    ISeedDP_D0 = 400_i8
    do i=1,5
        call xor4096g(ISeedDP_D0, DoubleRealRN_D0)
        print '(A14,1X,F18.15,A9,B64.64)', 'xor4096gd_0d = ',DoubleRealRN_D0
        ISeedDP_D0 = 0_i8
    end do

    print*,'-----------------------------------------'
    print*,'GAUSSIAN Double Precision (64 bit) Real  '
    print*,'          3 STREAMS                      '
    print*,'-----------------------------------------'

    ISeedDP = (/ 400_i8, 5_i8, 7_i8 /)
    do i=1,5
        call xor4096g(ISeedDP, DoubleRealRN)
        print '(A13,I6,A4)', 'xor4096gf_1d(',i,') = '
        do j=1,NumberOfStreams
            print '(A11,I4,A4,F18.15,A9,B64.64)', '   stream #',j,'  = ',DoubleRealRN(j), &
                                                  '   (bin) ',DoubleRealRN(j)
        end do
        ISeedDP = 0_i8
    end do

    print*,'-----------------------------------------'
    print*,'Gaussian Double Precision (64 bit) Real  '
    print*,'          WITH OPTIONAL                  '
    print*,'-----------------------------------------'
    ISeedDP_D0 = 400_i8
    do i=1,3
        call xor4096g(ISeedDP_D0, DoubleRealRN_D0, DoubleIin_D0,DoubleWin_D0,DoubleXin_D0,DoubleFlagin_D0,DoubleY2In_D0)
        print '(A12,I6,A4,F18.15,A9,B64.64)', 'xor4096gd_0d(',i,') = ',DoubleRealRN_D0, &
                                           '   (bin) ',DoubleRealRN_D0
        ISeedDP_D0 = 0_i8
    end do
    ! Start a new Stream
    call xor4096g(5_i8, DoubleRealRN_D0)
    print '(A12,I6,A4,F18.15,A9,B64.64)', 'xor4096gd_0d(',1,') = ',DoubleRealRN_D0, &
                                       '   (bin) ',DoubleRealRN_D0
    ISeedDP_D0 = 0_i8
    ! Go back to the old stream
    do i=4,5
        call xor4096g(ISeedDP_D0, DoubleRealRN_D0, DoubleIin_D0,DoubleWin_D0,DoubleXin_D0,DoubleFlagin_D0,DoubleY2In_D0)
        print '(A12,I6,A4,F18.15,A9,B64.64)', 'xor4096gd_0d(',i,') = ',DoubleRealRN_D0, &
                                           '   (bin) ',DoubleRealRN_D0
        ISeedDP_D0 = 0_i8
    end do

    deallocate(ISeedSP)
    deallocate(ISeedDP)
    deallocate(SingleIntegerRN)
    deallocate(DoubleIntegerRN)
    deallocate(SingleRealRN)
    deallocate(DoubleRealRN)
    
end program

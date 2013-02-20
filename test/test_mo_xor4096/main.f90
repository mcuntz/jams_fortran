!*******************************************************
!
!   TO TEST RNGs
!
!*******************************************************

program RNG

    use mo_kind,  only: i4, i8, SP, DP
    use mo_xor4096, only: xor4096, xor4096g, n_save_state

    implicit none

    integer(i4)                    :: ISeedSP_D0
    integer(i4), allocatable       :: ISeedSP(:)
    integer(i8), allocatable       :: ISeedDP(:)

    integer(i4)                    :: SingleIntegerRN_D0
    integer(i4), allocatable       :: SingleIntegerRN(:)
    integer(i8), allocatable       :: DoubleIntegerRN(:)
    real(SP),    allocatable       :: SingleRealRN(:)
    real(DP),    allocatable       :: DoubleRealRN(:)

    integer(i4)                    :: NumberOfStreams = 3

    character(80)                  :: CheckString, CheckStringShort
    integer(i4), allocatable       :: CheckISP_D0(:)
    integer(i4), allocatable       :: CheckISP_D1(:,:)

    integer(i8)                    :: i,j
    
    ! Needed for optional versions: Single
    integer(i4), dimension(n_save_state)   :: save_state_d0
    integer(i4), dimension(3, n_save_state) :: save_state_d1_3
    integer(i4), dimension(2, n_save_state) :: save_state_d1_2

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
    CheckString = 'mo_xor4096: Uniform Single Integer (1 stream) o.k.'
    do i=1,5
        call xor4096(ISeedSP_D0, SingleIntegerRN_D0)
        if (SingleIntegerRN_D0 .eq. CheckISP_D0(i)) then
            CheckStringShort = 'ok'
        else
            CheckString = 'mo_xor4096: Uniform Single Integer (1 stream) failed'
            CheckStringShort = 'fail'
        end if
        print '(A8,I6,A4,I12,A4,A10)', 'xor4096(',i,') = ',SingleIntegerRN_D0,'    ',CheckStringShort
        ISeedSP_D0 = 0_i4
    end do
    print *, CheckString
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
    CheckString = 'mo_xor4096: Uniform Single Integer (3 streams) o.k.'
    CheckStringShort = 'ok'
    do i=1,5
        call xor4096(ISeedSP,SingleIntegerRN)
        print '(A8,I6,A4)', 'xor4096(',i,') = '
        do j=1,NumberOfStreams
            if (SingleIntegerRN(j) .eq. CheckISP_D1(i,j)) then
                CheckStringShort = 'ok'
            else
                CheckString = 'mo_xor4096: Uniform Single Integer (3 streams) failed'
                CheckStringShort = 'fail'
            end if
            print '(A11,I4,A4,I12,A4,A10)', '   stream #',j,'  = ',SingleIntegerRN(j),'    ',CheckStringShort
        end do
        ISeedSP = 0_i4
    end do
    print *, CheckString
    deallocate (CheckISP_D1)

    print*,'-----------------------------------------'
    print*,'UNIFORM Single Precision (32 bit) Integer'
    print*,'          WITH OPTIONAL                  '
    print*,'-----------------------------------------'
    ISeedSP_D0 = 3_i4
    allocate (CheckISP_D0(6))
    CheckISP_D0 = (/-1381736336,1670483215,-70622511,1980498182,-909439004,-1632616910/)
    CheckString = 'mo_xor4096: Uniform Single Integer (with optional) o.k.'
    CheckStringShort = 'ok'
    do i=1,3
        call xor4096(ISeedSP_D0,SingleIntegerRN_D0,save_state=save_state_d0)
        if (SingleIntegerRN_D0 .eq. CheckISP_D0(i)) then
            CheckStringShort = 'ok'
        else
            CheckString = 'mo_xor4096: Uniform Single Integer (1st stream) failed'
            CheckStringShort = 'fail'
        end if
        print '(A8,I6,A4,I12,A4,A10)', 'xor4096(',i,') = ',SingleIntegerRN_D0,'    ',CheckStringShort
        ISeedSP_D0 = 0_i4
    end do
    ! Start a new Stream
    CheckStringShort = 'ok'
    call xor4096(5_i4,SingleIntegerRN_D0)
    if (SingleIntegerRN_D0 .eq. CheckISP_D0(6)) then
        CheckStringShort = 'ok'
    else
        CheckString = 'mo_xor4096: Uniform Single Integer (2nd stream) failed'
        CheckStringShort = 'fail'
    end if
    print '(A8,I6,A4,I12,A4,A10)', 'xor4096(',1,') = ',SingleIntegerRN_D0,'    ',CheckStringShort
    ISeedSP_D0 = 0_i4
    ! Go back to the old stream
    CheckStringShort = 'ok'
    do i=4,5
        call xor4096(ISeedSP_D0,SingleIntegerRN_D0,save_state=save_state_d0)
        if (SingleIntegerRN_D0 .eq. CheckISP_D0(i)) then
            CheckStringShort = 'ok'
        else
            CheckString = 'mo_xor4096: Uniform Single Integer (restart 1st stream) failed'
            CheckStringShort = 'fail'
        end if
        print '(A8,I6,A4,I12,A4,A10)', 'xor4096(',i,') = ',SingleIntegerRN_D0,'    ',CheckStringShort
        ISeedSP_D0 = 0_i4
    end do
    print *, CheckString
    deallocate (CheckISP_D0)

    print*,'----------------------------------------'
    print*,'GAUSSIAN Single Precision               '
    print*,'          3 STREAMS and 2 STREAMS       '
    print*,'----------------------------------------'

    print*, '3 streams'
    ISeedSP = (/ 3_i4, 5_i4, 7_i4 /)
    call xor4096g(ISeedSP,SingleRealRN, save_state=save_state_d1_3)
    print*, SingleRealRN
    ISeedSP = 0_i4
    call xor4096g(ISeedSP,SingleRealRN, save_state=save_state_d1_3)
    print*, SingleRealRN
    call xor4096g(ISeedSP,SingleRealRN, save_state=save_state_d1_3)
    print*, SingleRealRN

    print*, '2 streams'
    call xor4096g( (/ 2_i4, 37_i4 /), SingleRealRN(1:2), save_state=save_state_d1_2)
    print*, SingleRealRN(1:2)
    ISeedSP = 0_i4
    call xor4096g((/ 0_i4, 0_i4 /), SingleRealRN(1:2), save_state=save_state_d1_2)
    print*, SingleRealRN(1:2)
    call xor4096g((/ 0_i4, 0_i4 /), SingleRealRN(1:2), save_state=save_state_d1_2)
    print*, SingleRealRN(1:2)

    print*, 'again 3 streams'
    ISeedSP = (/ 3_i4, 5_i4, 7_i4 /)
    call xor4096g(ISeedSP,SingleRealRN, save_state=save_state_d1_3)
    print*, SingleRealRN
    ISeedSP = 0_i4
    call xor4096g(ISeedSP,SingleRealRN, save_state=save_state_d1_3)
    print*, SingleRealRN
    call xor4096g(ISeedSP,SingleRealRN, save_state=save_state_d1_3)
    print*, SingleRealRN

    deallocate(ISeedSP)
    deallocate(ISeedDP)
    deallocate(SingleIntegerRN)
    deallocate(DoubleIntegerRN)
    deallocate(SingleRealRN)
    deallocate(DoubleRealRN)

end program RNG





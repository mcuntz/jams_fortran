program main
  
  use mo_kind,   only: i4, i8, dp, sp
  use mo_string_utils, only: tolower, toupper, separator, num2str, nonull, compress
  use mo_string_utils, only: startsWith, equalStrings, splitString, countsubstring
  use mo_string_utils, only: count_split, split
#ifndef ABSOFT
  use mo_string_utils, only: divide_string
#endif

  implicit none
  
  character(len=100)                        :: sout, sundef
#ifndef ABSOFT
  character(256), dimension(:), allocatable :: strArr
#endif

  logical :: isgood

  write(*,*) ''
  write(*,*) 'Test mo_string_utils.f90'

  isgood = .true.
  ! tolower
  sout = tolower('Hallo')
  if (.not.(lle(trim(sout),'hallo') .and. lge(trim(sout),'hallo'))) isgood =.false.
  ! toupper
  sout = toupper('Hallo')
  if (.not.(lle(trim(sout),'HALLO') .and. lge(trim(sout),'HALLO'))) isgood =.false.
  ! num2str
  sout = separator
  if (.not.(lle(sout(1:3),'---') .and. lge(sout(1:3),'---'))) isgood =.false.
  sout = num2str(3.141592653589793238462643383279502884197_dp, '(F5.3)')
  if (.not.(lle(trim(sout),'3.142') .and. lge(trim(sout),'3.142'))) isgood =.false.
  sout = num2str(3.141592653589793238462643383279502884197_sp, '(F5.3)')
  if (.not.(lle(trim(sout),'3.142') .and. lge(trim(sout),'3.142'))) isgood =.false.
  sout = num2str(101_i4, '(I3)')
  if (.not.(lle(trim(sout),'101') .and. lge(trim(sout),'101'))) isgood =.false.
  sout = num2str(101_i8, '(I3)')
  if (.not.(lle(trim(sout),'101') .and. lge(trim(sout),'101'))) isgood =.false.
  sout = num2str(.true., '(L1)')
  if (.not.(lle(trim(sout),'T') .and. lge(trim(sout),'T'))) isgood =.false.
  ! nonull
  if (.not. nonull(sout)) isgood =.false.
  if (nonull(sundef)) isgood =.false.
  ! compress
  sout = compress('H a      l l o       ')
  if (.not.(lle(trim(sout),'Hallo') .and. lge(trim(sout),'Hallo'))) isgood =.false.
  ! startsWith
  if (.not. (startsWith("Thisisatest","T") .and. startsWith("Thisisatest","This"))) isgood = .false.
  if ((startsWith("Thisisatest","t") .or. startsWith("Thisisatest","test_*"))) isgood = .false.
  ! equalStrings
  if (.not. (equalStrings("Thisis","Thisis") .and. equalStrings("test*_<","test*_<"))) isgood = .false.
  if ((equalStrings("Thisis","THISIS") .or. equalStrings("test*_<","est*_<"))) isgood = .false.
  if (.not. isgood) write(*,*) 'mo_string_utils failed after equalstrings'
  ! splitString
  if (allocated(strArr)) deallocate(strArr)
  allocate(strArr(countsubstring('I want to test this routine!', ' ')))
  strArr = splitString('I want to test this routine!', ' ')
  isgood = isgood .and. (strArr(1) .eq. 'I')
  isgood = isgood .and. (strArr(2) .eq. 'want')
  isgood = isgood .and. (strArr(3) .eq. 'to')
  isgood = isgood .and. (strArr(4) .eq. 'test')
  isgood = isgood .and. (strArr(5) .eq. 'this')
  isgood = isgood .and. (strArr(6) .eq. 'routine!')
  if (.not. isgood) write(*,*) 'mo_string_utils splitstring 1 failed'
  if (allocated(strArr)) deallocate(strArr)
  allocate(strArr(countsubstring('I,want,to,test,this,routine! ', ',')))
  strArr = splitString('I,want,to,test,this,routine! ', ',')
  isgood = isgood .and. (strArr(1) .EQ. 'I')
  isgood = isgood .and. (strArr(2) .EQ. 'want')
  isgood = isgood .and. (strArr(3) .EQ. 'to')
  isgood = isgood .and. (strArr(4) .EQ. 'test')
  isgood = isgood .and. (strArr(5) .EQ. 'this')
  isgood = isgood .and. (strArr(6) .EQ. 'routine! ')
  if (.not. isgood) write(*,*) 'mo_string_utils splitstring 2 failed'
  strArr = splitString('I,want,to,test,this,routine! ', ',')
  isgood = isgood .and. (strArr(6) .EQ. 'routine!')
  if (.not. isgood) write(*,*) 'mo_string_utils splitstring 2.1 failed'
  if (allocated(strArr)) deallocate(strArr)
  allocate(strArr(countsubstring('w!hat_s-a+bout=-sp.eci,al-chara<cte>rs?', '-')))
  strArr = splitString('w!hat_s-a+bout=-sp.eci,al-chara<cte>rs?', '-')
  isgood = isgood .and. (strArr(1) .EQ. 'w!hat_s')
  isgood = isgood .and. (strArr(2) .EQ. 'a+bout=')
  isgood = isgood .and. (strArr(3) .EQ. 'sp.eci,al')
  isgood = isgood .and. (strArr(4) .EQ. 'chara<cte>rs?')
  if (.not. isgood) write(*,*) 'mo_string_utils splitstring 3 failed'
  ! divide_string does not allow multi-character splits, so pgi does not work
  if (allocated(strArr)) deallocate(strArr)
  allocate(strArr(countsubstring('multi_+character_*splits_+should work_+', '_+')))
  strArr = splitString('multi_+character_*splits_+should work_+', '_+')
  isgood = isgood .and. (strArr(1) .EQ. 'multi')
  isgood = isgood .and. (strArr(2) .EQ. 'character_*splits')
  isgood = isgood .and. (strArr(3) .EQ. 'should work')
  isgood = isgood .and. (strArr(4) .EQ. '')
  if (.not. isgood) write(*,*) 'mo_string_utils splitstring 4 failed'
  ! split
  if (allocated(strArr)) deallocate(strArr)
  allocate(strArr(count_split('I want to test this routine!', ' ')))
  strArr = split('I want to test this routine!', ' ')
  isgood = isgood .and. (strArr(1) .eq. 'I')
  isgood = isgood .and. (strArr(2) .eq. 'want')
  isgood = isgood .and. (strArr(3) .eq. 'to')
  isgood = isgood .and. (strArr(4) .eq. 'test')
  isgood = isgood .and. (strArr(5) .eq. 'this')
  isgood = isgood .and. (strArr(6) .eq. 'routine!')
  if (.not. isgood) write(*,*) 'mo_string_utils split 5 failed'
  if (allocated(strArr)) deallocate(strArr)
  allocate(strArr(count_split('I,want,to,test,this,routine! ', ',')))
  strArr = split('I,want,to,test,this,routine! ', ',')
  isgood = isgood .and. (strArr(1) .EQ. 'I')
  isgood = isgood .and. (strArr(2) .EQ. 'want')
  isgood = isgood .and. (strArr(3) .EQ. 'to')
  isgood = isgood .and. (strArr(4) .EQ. 'test')
  isgood = isgood .and. (strArr(5) .EQ. 'this')
  isgood = isgood .and. (strArr(6) .EQ. 'routine! ')
  if (.not. isgood) write(*,*) 'mo_string_utils split 6 failed'
  strArr = split('I,want,to,test,this,routine! ', ',')
  isgood = isgood .and. (strArr(6) .EQ. 'routine!')
  if (.not. isgood) write(*,*) 'mo_string_utils split 6.1 failed'
  if (allocated(strArr)) deallocate(strArr)
  allocate(strArr(count_split('w!hat_s-a+bout=-sp.eci,al-chara<cte>rs?', '-')))
  strArr = split('w!hat_s-a+bout=-sp.eci,al-chara<cte>rs?', '-')
  isgood = isgood .and. (strArr(1) .EQ. 'w!hat_s')
  isgood = isgood .and. (strArr(2) .EQ. 'a+bout=')
  isgood = isgood .and. (strArr(3) .EQ. 'sp.eci,al')
  isgood = isgood .and. (strArr(4) .EQ. 'chara<cte>rs?')
  if (.not. isgood) write(*,*) 'mo_string_utils split 7 failed'
  if (allocated(strArr)) deallocate(strArr)
  allocate(strArr(count_split('multi_+character_*splits_+should work_+', '_+')))
  strArr = split('multi_+character_*splits_+should work_+', '_+')
  isgood = isgood .and. (strArr(1) .EQ. 'multi')
  isgood = isgood .and. (strArr(2) .EQ. 'character_*splits')
  isgood = isgood .and. (strArr(3) .EQ. 'should work')
  isgood = isgood .and. (strArr(4) .EQ. '')
  if (.not. isgood) write(*,*) 'mo_string_utils split 8 failed'

#ifndef ABSOFT
  call DIVIDE_STRING('I want to test this routine!', ' ', strArr)
  isgood = isgood .and. (strArr(1) .EQ. 'I')
  isgood = isgood .and. (strArr(2) .EQ. 'want')
  isgood = isgood .and. (strArr(3) .EQ. 'to')
  isgood = isgood .and. (strArr(4) .EQ. 'test')
  isgood = isgood .and. (strArr(5) .EQ. 'this')
  isgood = isgood .and. (strArr(6) .EQ. 'routine!')
  call DIVIDE_STRING('I,want,to,test,this,routine!', ',', strArr)
  isgood = isgood .and. (strArr(1) .EQ. 'I')
  isgood = isgood .and. (strArr(2) .EQ. 'want')
  isgood = isgood .and. (strArr(3) .EQ. 'to')
  isgood = isgood .and. (strArr(4) .EQ. 'test')
  isgood = isgood .and. (strArr(5) .EQ. 'this')
  isgood = isgood .and. (strArr(6) .EQ. 'routine!')
  call DIVIDE_STRING("w!hat_s-a+bout=-sp.eci,al-chara<cte>rs?", '-', strArr)
  isgood = isgood .and. (strArr(1) .EQ. 'w!hat_s')
  isgood = isgood .and. (strArr(2) .EQ. 'a+bout=')
  isgood = isgood .and. (strArr(3) .EQ. 'sp.eci,al')
  isgood = isgood .and. (strArr(4) .EQ. 'chara<cte>rs?')
#endif

  if (allocated(strArr)) deallocate(strArr)

  Write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_string_utils o.k.'
  else
     write(*,*) 'mo_string_utils failed'
  endif

END PROGRAM main

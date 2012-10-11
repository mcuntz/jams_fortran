PROGRAM main
  
  USE mo_kind,   ONLY: i4, i8, dp, sp
  USE mo_string_utils, ONLY: tolower, toupper, separator, num2str, nonull
#ifndef ABSOFT
  USE mo_string_utils, ONLY: DIVIDE_STRING
#endif

  IMPLICIT NONE
  
  CHARACTER(len=100)                        :: sout, sundef
#ifndef ABSOFT
  CHARACTER(256), dimension(:), allocatable :: strArr
#endif

  LOGICAL :: isgood

  Write(*,*) ''
  Write(*,*) 'Test mo_string_utils.f90'

  isgood = .true.
  sout = tolower('Hallo')
  if (.not.(lle(trim(sout),'hallo') .and. lge(trim(sout),'hallo'))) isgood =.false.
  sout = toupper('Hallo')
  if (.not.(lle(trim(sout),'HALLO') .and. lge(trim(sout),'HALLO'))) isgood =.false.
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
  if (.not. nonull(sout)) isgood =.false.
  if (nonull(sundef)) isgood =.false.
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

  Write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_string_utils o.k.'
  else
     write(*,*) 'mo_string_utils failed'
  endif

END PROGRAM main

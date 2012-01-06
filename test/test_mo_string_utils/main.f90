PROGRAM main
  
  USE mo_kind,   ONLY: i4, i8, dp, sp
  USE mo_string_utils, ONLY: tolower, toupper, separator, num2str, nullstring

  IMPLICIT NONE
  
  CHARACTER(len=100) :: sout, sundef

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
  if (.not. nullstring(sout)) isgood =.false.
  if (nullstring(sundef)) isgood =.false.

  Write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_string_utils o.k.'
  else
     write(*,*) 'mo_string_utils failed'
  endif

END PROGRAM main

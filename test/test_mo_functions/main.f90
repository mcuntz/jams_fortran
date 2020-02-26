program main

  use mo_kind,      only: i4, i8, sp, dp
  use mo_functions, only: factorial, factln, gamm, gammln
  use mo_ansi_colors, only: color, c_red, c_green

  implicit none

  logical                              :: isgood

  isgood = .true.
  write(*,*) ''
  write(*,*) 'Test mo_functions'

  ! Factorial
  if ( factorial(5_i4) .ne. 120_i4 ) isgood = .false.
  if ( factorial(5_i8) .ne. 120_i8 ) isgood = .false.

  ! factln
  if ( nint(exp(factln(5_i4)),i4) .ne. 120_i4 ) isgood = .false.
  if ( nint(exp(factln(5_i8)),i8) .ne. 120_i8 ) isgood = .false.

  ! gamm
  if ( nint(gamm(6._sp),i4) .ne. 120_i4 ) isgood = .false.
  if ( nint(gamm(6._dp),i8) .ne. 120_i8 ) isgood = .false.

  ! gammln
  if ( nint(exp(gammln(6._sp)),i4) .ne. 120_i4 ) isgood = .false.
  if ( nint(exp(gammln(6._dp)),i8) .ne. 120_i8 ) isgood = .false.

  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_functions ', color('o.k.', c_green)
  else
     write(*,*) 'mo_functions ', color('failed', c_red)
  end if

end program main

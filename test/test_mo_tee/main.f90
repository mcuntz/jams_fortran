program main

  use mo_tee,  only: tee

  implicit none

  character(len=*), parameter :: nfile1 = '../fortran/test/test_mo_tee/messages.txt'
  character(len=*), parameter :: nfile2 = 'tee_make_check_test_file'

  character(len=100) :: iread1, iread2
  integer :: ierr1, ierr2

  logical :: isgood, allgood

  write(*,*) ''
  write(*,*) 'Test mo_tee.f90'

  allgood = .true.

  ! tee_filename
  isgood = .true.
  write(*,*) ''
  call tee(nfile2, 'Message 1')
  call tee(nfile2, 'Message', ' ', '1', overwrite=.true.)
  call tee(nfile2, 'Another ', 'message')

  open(98, file=nfile1, status='old', action='read', form="formatted")
  open(99, file=nfile2, status='old', action='read', form="formatted")
  ierr1 = 0
  ierr2 = 0
  do while ((ierr1==0) .and. (ierr2==0))
     read(98, '(a)', iostat=ierr1) iread1
     read(99, '(a)', iostat=ierr2) iread2
     if (trim(iread1) /= trim(iread2)) isgood = .false.
  end do
  close(98)
  close(99)

  allgood = allgood .and. isgood

  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_tee filename o.k.'
  else
     write(*,*) 'mo_nml filename failed!'
  endif


  ! tee_unit
  isgood = .true.
  write(*,*) ''
  open(97, file=nfile2, action='write', form="formatted", status='replace')
  call tee(97, 'Message', advance='no')
  call tee(97, ' ', '1')
  call tee(97, 'Another ', 'message')
  close(97)

  open(98, file=nfile1, status='old', action='read', form="formatted")
  open(99, file=nfile2, status='old', action='read', form="formatted")
  ierr1 = 0
  ierr2 = 0
  do while ((ierr1==0) .and. (ierr2==0))
     read(98, '(a)', iostat=ierr1) iread1
     read(99, '(a)', iostat=ierr2) iread2
     if (trim(iread1) /= trim(iread2)) isgood = .false.
  end do
  close(98)
  close(99)

  allgood = allgood .and. isgood

  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_tee unit o.k.'
  else
     write(*,*) 'mo_nml unit failed!'
  endif

  write(*,*) ''
  if (allgood) then
     write(*,*) 'mo_tee o.k.'
  else
     write(*,*) 'mo_nml failed!'
  endif

end program main

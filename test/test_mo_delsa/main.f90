program main

  ! This programs produces the 2-parameter toy model presented in Rakovec et al (2014, WRR)
  ! which calculates first-order DELSA method shown in bottom Figure 5.

  use mo_kind,        only: dp
  use mo_delsa,       only: delsa
  use mo_ansi_colors, only: color, c_red, c_green

  implicit none

  integer :: ii, jj, set

  integer, parameter :: npara = 2    ! number of parameter
  integer, parameter :: nsets = 100  ! number of parametersets

  real(dp), dimension(nsets,npara) :: parbase            ! array of base model parameters
  real(dp), dimension(nsets,npara) :: parpert            ! array of perturbed model parameters (alltogether)
  real(dp), dimension(nsets,npara) :: parpert1, parpert2 ! array of perturbed model parameters (individual param)
  real(dp), dimension(nsets)       :: outbase            ! model output using the base runs
  real(dp), dimension(nsets,npara) :: outpert            ! model output using the perturbed runs
  real(dp), dimension(npara)       :: varprior           ! prior parameter variance

  real(dp), dimension(nsets,npara) :: delsafirst, delsatest ! first-order DELSA results,delsa test results for a check

  logical :: isgood = .true.

  ! Data Files for delsa(parbase,parpert,outbase,outpert,varprior,delsafirst)
  ! read text files which were created externally:
  character(len=*), parameter :: fileparbase    = 'dummy_data/parbase.txt'
  character(len=*), parameter :: fileparpert1   = 'dummy_data/parpert1.txt'
  character(len=*), parameter :: fileparpert2   = 'dummy_data/parpert2.txt'
  character(len=*), parameter :: fileoutbase    = 'dummy_data/outbase.txt'
  character(len=*), parameter :: fileoutpert    = 'dummy_data/outpert.txt'
  character(len=*), parameter :: filetestoutput = 'dummy_data/test_output.txt'

  ! ************************************************************************
  ! read parameter sets and corresponding model outputs
  ! ************************************************************************

  ! read base parameter set:
  open(unit=50, file=fileparbase, action='read')
  do set=1, size(parbase, 1)
     read(50,*) parbase(set,:)
  end do
  close(50)

  ! read parameter set, in which is perturbed parameter 1:
  open(unit=50, file=fileparpert1, action='read')
  do set=1, size(parbase, 1)
     read(50,*) parpert1(set, :)
  end do
  close(50)

  ! read parameter set, in which is perturbed parameter 2:
  open(unit=50, file=fileparpert2, action='read')
  do set=1, size(parbase, 1)
     read(50,*) parpert2(set, :)
  end do
  close(50)

  ! compose an array with only perturbed parameters
  parpert(:, 1) = parpert1(:, 1)
  parpert(:, 2) = parpert2(:, 2)

  ! read model output using the base parameter set
  open(unit=50, file=fileoutbase, action='read')
  do set=1, size(outbase, 1)
     read(50,*) outbase(set)
  end do
  close(50)

  ! read model output using the perturbed parameter sets
  open(unit=50, file=fileoutpert, action='read')
  do set=1, size(outpert, 1)
     read(50,*) outpert(set, :)
  end do
  close(50)

  ! prior variance (this comes from last sentence in Rakovec et al, 2014, WRR of Sect. 3.2)
  ! note, there is a typo in the paper, those values are standard deviations and not variances!
  varprior = (/28.83865_dp**2, 0.08660254_dp**2/)
  write(*,*) 'varprior: '
  write(*,*) varprior
  write(*,*) ' '

  call delsa(parbase, parpert, outbase, outpert, varprior, delsafirst)
  write(*,*) 'delsa first, param 1:'
  write(*,*) delsafirst(:, 1)
  write(*,*) ' '
  write(*,*) 'delsa first, param 2:'
  write(*,*) delsafirst(:, 2)
  write(*,*) ' '

  ! save results in a text file:
  open(unit=999, file='delsafirst_fortran.txt.make_check_test_file', action='write')

  do set=1, size(outbase, 1)
     write(999,'(f12.6)', advance='no') delsafirst(set, 1)
     write(999,'(f12.6)') delsafirst(set, 2)
  end do
  close(999)
  write(*,*) 'results saved in delsafirst_fortran.txt.make_check_test_file'


  ! ************************************************************************
  !       Verify that module runs correctly
  ! ************************************************************************
  open(unit=50, file=filetestoutput, action='read')
  do set=1, size(delsatest, 1)
     read(50,*) delsatest(set, :)
  end do
  close(50)

  do ii=1, size(delsatest, 1)
     do jj=1, size(delsatest, 2)
     if (nint(delsafirst(ii, jj)*1000000._dp) .ne. nint(delsatest(ii, jj)*1000000._dp)) isgood = .false.
     end do
  end do

  if (isgood) then
     write(*,*) 'mo_delsa ', color('o.k.', c_green)
  else
     write(*,*) 'mo_delsa ', color('failed', c_red)
  end if

end program main

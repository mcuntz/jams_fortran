program pi_index_test

  use mo_kind,        only: dp,i4, i8
  use mo_ansi_colors, only: color, c_red, c_green
  use mo_sobol,       only: sobol
  use mo_model,       only: getrange, model
  use mo_pi_index,    only: pi_index

  implicit none

  integer(i4), parameter                 :: npara = 4                    ! number of parameter
  integer(i4), parameter                 :: nsets = 400                  ! number of parametersets 
  integer(i4), parameter                 :: nx = 10                      ! number of variable values
  real(dp),    parameter                 :: delta = 0.01_dp              ! delta used for derivative approximation
  integer(i8)                            :: skip = 30000_i8              ! used for generating parameter sets via sobol sequences

  ! Samples and Model Output
  real(dp),    dimension(npara,2)        :: ParaRange                    ! array: Parameter ranges (Min, Max)
  real(dp),    dimension(nx)             :: x                            ! array: variable values, e.g. time points

  ! PI index
  real(dp),    dimension(npara)          :: sample, sample_delta         ! parameter sets: nominal and 1% changed set
  real(dp),    dimension(nx)             :: output_ref, output_delta     ! model output (time series) 
  !                                                                      ! for nominal and 1% changed parameter set
  real(dp),    dimension(nsets*nx,npara) :: s_matrix                     ! sensitivity matrix S
  real(dp),    dimension(npara)          :: pi                           ! Parameter importance index 
  real(dp),    dimension(:), allocatable :: bi                           ! B index
  integer(i4), dimension(:), allocatable :: counter_pi                   ! number of valid entries in ith column of matrix S 
  integer(i8)                            :: seed                         ! seed sobol sequence
 
  ! Dummy variables
  integer(i4)                            :: set, para, i
  logical                                :: isgood = .true.

  ParaRange = GetRange()

  ! ------------------------------------------------------------------
  ! variable values, e.g. time points
  ! ------------------------------------------------------------------
  do i=1,nx
     x(i) = real(i,dp) * 0.1_dp
  end do

  ! ------------------------------------------------------------------
  ! Generating sensitivity matrix S using Sobol sequences
  ! ------------------------------------------------------------------
  seed = 0
  call sobol(int(npara,i8), seed, sample)
  do set = 1,nsets
     ! sobol sequence
     call sobol(int(npara,i8), skip, sample)
     ! scaling
     sample(:)  = ParaRange(:,1) + sample(:) * (ParaRange(:,2) - ParaRange(:,1)) 
     output_ref = model(sample,x)
     do para=1,npara
        sample_delta(:)    = sample(:)
        sample_delta(para) = (1.+delta) * sample_delta(para)
        output_delta       = model(sample_delta,x)
        s_matrix((set-1)*nx+1:set*nx,para) = ( output_delta - output_ref ) / ( delta ) * sample(para) / output_ref
     end do
  end do

  ! ------------------------------------------------------------------
  ! Calculate Parameter Importance (PI) index
  ! ------------------------------------------------------------------
  call pi_index(pi, s=s_matrix, norm=1_i4, domad=.true., counter=counter_pi, b_index=bi)
  
  ! ------------------------------------------------------------------
  ! Printing and checking
  ! ------------------------------------------------------------------
  write(*,'(A48)') '------------------------------------------------'
  write(*,'(A48)') '        PI index    &    B index                '
  write(*,'(A48)') '------------------------------------------------'
  write(*,'(A48)') 'para  PI       B index     (number valid values)'
  do para=1,npara
     write(*,'(I4,A2,F7.4,A2,F7.4,A6,I5,A15)') para, '  ', pi(para), '  ', bi(para), '     (', counter_pi(para),'              )'
  end do
  write(*,'(A48)') '------------------------------------------------'
  write(*,*) ' '

  if (nint(pi(1)*10000._dp)        .ne.  702_i4)  isgood = .false.
  if (nint(pi(2)*10000._dp)        .ne. 1415_i4)  isgood = .false.
  if (nint(pi(3)*10000._dp)        .ne. 3059_i4)  isgood = .false.
  if (nint(pi(4)*10000._dp)        .ne. 4824_i4)  isgood = .false.
  if (nint(bi(1)*10000._dp)        .ne.  271_i4)  isgood = .false.
  if (nint(bi(2)*10000._dp)        .ne.  788_i4)  isgood = .false.
  if (nint(bi(3)*10000._dp)        .ne. 2659_i4)  isgood = .false.
  if (nint(bi(4)*10000._dp)        .ne. 6282_i4)  isgood = .false.

  if (isgood) then
     write(*,*) 'mo_pi_index ', color('o.k.', c_green)
  else
     write(*,*) 'mo_pi_index ', color('failed', c_red)
  end if

END PROGRAM pi_index_test

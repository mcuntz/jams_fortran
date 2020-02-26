program sobol_index_test

  use mo_kind,        only: dp,i4, i8
  use mo_ansi_colors, only: color, c_red, c_green
  use mo_sobol,       only: sobol
  use mo_model,       only: getrange, model
  use mo_sobol_index, only: sobol_index

  implicit none

  integer(i4), parameter                 :: npara = 3                       ! number of parameter
  integer(i4), parameter                 :: nsets = 1000                    ! number of parametersets 
  integer(i4), parameter                 :: nx = 100                        ! number of variable values
  integer(i8)                            :: skip = 30000_i8                 ! used for generating parameter sets 
  !                                                                         ! via sobol sequences

  ! Samples and Model Output
  real(dp),    dimension(npara,2)        :: ParaRange                       ! array: Parameter ranges (Min, Max)
  real(dp),    dimension(nx)             :: x                               ! array: variable values

  ! Sobol index
  real(dp),    dimension(2*npara)        :: sample_sobol                    ! sobol sequence
  real(dp),    dimension(npara)          :: sample_a                        ! parameter sample A
  real(dp),    dimension(npara)          :: sample_b                        ! parameter sample B
  real(dp),    dimension(npara)          :: sample_ci                       ! parameter sample Ci
  real(dp),    dimension(nsets,nx)       :: ya, yb                          ! model output A and B, 
  !                                                                         !    i.e. time series with nx time points
  real(dp),    dimension(nsets,npara,nx) :: yc                              ! model output C(i), i=1,npara, 
  !                                                                         !    i.e. time series with nx time points
  real(dp),    dimension(npara,nx)       :: si                              ! Sobol index (main effect)  per time point and param.
  !                                                                         !    based on ya, yb, and yc
  real(dp),    dimension(npara,nx)       :: sti                             ! Sobol index (total effect) per time point and param.
  !                                                                         !    based on ya, yb, and yc
  real(dp),    dimension(npara,2)        :: smean                           ! Mean SI and STI
  real(dp),    dimension(npara,2)        :: wmean                           ! Variance weighted mean SI and STI

  integer(i8) :: seed

  ! Dummy variables
  integer(i4)                            :: set, para, i
  logical                                :: isgood = .true.

  ParaRange = GetRange()

  ! ------------------------------------------------------------------
  ! variable values, e.g. time points
  ! ------------------------------------------------------------------
  do i=1,nx
     x(i) = real(i-1,dp) * 1.0_dp + 50.0_dp
  end do

  ! ------------------------------------------------------------------
  ! Generating model outputs for parameter sets A, B, and C(i)
  ! ------------------------------------------------------------------
  seed = 0_i8
  call sobol(int(2*npara,i8), seed, sample_sobol)
  do set = 1,nsets
     ! sobol sequence
     call sobol(int(2*npara,i8), skip, sample_sobol)
     sample_a = sample_sobol(1:npara)
     sample_b = sample_sobol(npara+1:2*npara)
     ! scaling
     sample_a = ParaRange(:,1) + sample_a * (ParaRange(:,2) - ParaRange(:,1)) 
     sample_b = ParaRange(:,1) + sample_b * (ParaRange(:,2) - ParaRange(:,1)) 

     ya(set,:) = model(sample_a,x)
     yb(set,:) = model(sample_b,x)
     
     do para=1,npara
        sample_ci       = sample_b
        sample_ci(para) = sample_a(para)
        yc(set,para,:)  = model(sample_ci,x)
     end do
  end do

  ! ------------------------------------------------------------------
  ! Calculate Sobol index: Main Effect (Si) and Total Effect (STi)
  ! ------------------------------------------------------------------
  call sobol_index(ya, yb, yc, si, sti, method=6, smean=smean, wmean=wmean)
  
  ! ------------------------------------------------------------------
  ! Printing and checking
  ! ------------------------------------------------------------------
  write(*,'(A40)')         '----------------------------------------'
  write(*,'(A24,F5.1,A8)') '     Sobol index   (x = ',x(3),')       '
  write(*,'(A40)')         '----------------------------------------'
  write(*,'(A40)')         'para  SI         STI                    '
  do para=1,npara
     write(*,'(I4,A2,F7.4,A4,F7.4)') para, '   ', si(para,3),'  ', sti(para,3)
  end do
  write(*,'(A40)') '----------------------------------------'
  write(*,*) ' '

  if (nint(si(1,3)*10000._dp)      .ne.  3226_i4)  isgood = .false.
  if (nint(si(2,3)*10000._dp)      .ne.  1404_i4)  isgood = .false.
  if (nint(si(3,3)*10000._dp)      .ne. -  58_i4)  isgood = .false.
  if (nint(sti(1,3)*10000._dp)     .ne.  7826_i4)  isgood = .false.
  if (nint(sti(2,3)*10000._dp)     .ne.  1490_i4)  isgood = .false.
  if (nint(sti(3,3)*10000._dp)     .ne.  5348_i4)  isgood = .false.

  write(*,'(A40)')         '----------------------------------------'
  write(*,'(A40)')         '     Mean Sobol index                   '
  write(*,'(A40)')         '----------------------------------------'
  write(*,'(A40)')         'para  SI         STI                    '
  do para=1,npara
     write(*,'(I4,A2,F7.4,A4,F7.4)') para, '   ', smean(para,1),'  ', smean(para,2)
  end do
  write(*,'(A40)') '----------------------------------------'
  write(*,*) ' '

  if (nint(smean(1,1)*10000._dp)     .ne.  2417_i4)  isgood = .false.
  if (nint(smean(2,1)*10000._dp)     .ne.  3689_i4)  isgood = .false.
  if (nint(smean(3,1)*10000._dp)     .ne.    13_i4)  isgood = .false.
  if (nint(smean(1,2)*10000._dp)     .ne.  5784_i4)  isgood = .false.
  if (nint(smean(2,2)*10000._dp)     .ne.  3820_i4)  isgood = .false.
  if (nint(smean(3,2)*10000._dp)     .ne.  3980_i4)  isgood = .false.

  write(*,'(A40)')         '----------------------------------------'
  write(*,'(A40)')         '     Var. weighted mean Sobol index     '
  write(*,'(A40)')         '----------------------------------------'
  write(*,'(A40)')         'para  SI         STI                    '
  do para=1,npara
     write(*,'(I4,A2,F7.4,A4,F7.4)') para, '   ', wmean(para,1),'  ', wmean(para,2)
  end do
  write(*,'(A40)') '----------------------------------------'
  write(*,*) ' '

  if (nint(wmean(1,1)*10000._dp)     .ne.  2319_i4)  isgood = .false.
  if (nint(wmean(2,1)*10000._dp)     .ne.  3963_i4)  isgood = .false.
  if (nint(wmean(3,1)*10000._dp)     .ne.    20_i4)  isgood = .false.
  if (nint(wmean(1,2)*10000._dp)     .ne.  5535_i4)  isgood = .false.
  if (nint(wmean(2,2)*10000._dp)     .ne.  4102_i4)  isgood = .false.
  if (nint(wmean(3,2)*10000._dp)     .ne.  3813_i4)  isgood = .false.

  if (isgood) then
     write(*,*) 'mo_sobol_index ', color('o.k.', c_green)
  else
     write(*,*) 'mo_sobol_index ', color('failed', c_red)
  end if

end program sobol_index_test

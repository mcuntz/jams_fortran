PROGRAM elemeffects_test

  use mo_kind,        only: dp,i4
  use mo_elemeffects, only: elemeffects
  use mo_model,       only: GetRange, Model

  IMPLICIT NONE

  integer(i4), parameter              :: npara = 4                       ! number of parameter
  integer(i4), parameter              :: nsets = 400                     ! number of parametersets 
                                                                         ! = 20 x (number of parameter) 
                                                                         !      x (number of parameter + 1)
  integer(i4), parameter              :: nx = 10                         ! number of variable values

  ! Data Files
  ! file: Morris samples
  character(len=*), parameter :: fileParasets    = '../FORTRAN_chs_lib/test/test_mo_elemeffects/data/Morris_parasets.dat'
  ! file: changed parameter in Morris samples
  character(len=*), parameter :: fileChangedPara = '../FORTRAN_chs_lib/test/test_mo_elemeffects/data/Morris_changedpara.dat'

  ! Samples and Model Output
  real(dp),    dimension(nsets,npara) :: MorrisSample                    ! array: Morris samples
  integer(i4), dimension(nsets)       :: MorrisChangedPara               ! array: changed parameter in Morris samples
  real(dp),    dimension(npara,2)     :: ParaRange                       ! array: Parameter ranges (Min, Max)
  real(dp),    dimension(nsets,npara) :: MorrisSampleScaled              ! array: Morris samples scaled to ranges
  real(dp),    dimension(nx)          :: x                               ! array: variable values
  real(dp),    dimension(nsets,nx)    :: ModelOutput                     ! array: model output for each parameter set

  ! Elementary Effects
  real(dp),    dimension(npara)       :: elemEffect                      ! Elementary Effects
  integer(i4), dimension(npara)       :: counter                         ! number of model runs the elementary effect is based on

  ! Dummy variables
  integer(i4)                         :: set, para, i
  real(dp)                            :: fdummy
  logical                             :: isgood = .true.

  ! ************************************************************************
  !       Reading Morris samples
  !       generated with e.g. MATLAB
  ! ************************************************************************

  open(unit=50, file=fileParasets, action='read')
  do set=1,nsets
     read(50,*) MorrisSample(set,:)
  end do
  close(50)

  open(unit=50, file=fileChangedPara, action='read')
  do set=1,nsets
     read(50,*) fdummy
     MorrisChangedPara(set) = int(fdummy,i4)
  end do
  close(50)

  ! ************************************************************************
  !       Running the model with Morris sets
  !       --> Morris sets have to be scaled to parameter ranges first
  ! ************************************************************************
  
  ! scaling
  ParaRange = GetRange()
  do para=1,npara
     MorrisSampleScaled(:,para) = ParaRange(para,1) + MorrisSample(:,para) * (ParaRange(para,2) - ParaRange(para,1)) 
  end do

  ! variable values
  do i=1,nx
     x(i) = real(i,dp) * 0.1_dp
  end do

  ! running model
  do set=1,nsets
     ModelOutput(set,:) = model(MorrisSampleScaled(set,:),x)
  end do

  ! ************************************************************************
  !       Calculating Elementary Effects
  ! ************************************************************************
  
  call elemeffects(ModelOutput,MorrisSample,MorrisChangedPara,elemeffect,counter)
  
  write(*,'(A40)') '----------------------------------------'
  write(*,'(A40)') '        Elementary Effects              '
  write(*,'(A40)') '----------------------------------------'
  write(*,'(A40)') 'para  ElemEffect  (number model runs)   '
  do para=1,npara
     write(*,'(I4,A2,F7.4,A6,I5,A13)') para, '  ', elemEffect(para), '     (', counter(para),'            )'
  end do
  write(*,'(A40)') '----------------------------------------'
  write(*,*) ' '

  ! ************************************************************************
  !       Verify that module runs correctly
  ! ************************************************************************
  if (any(counter .ne. 80_i4))                    isgood = .false.
  if (nint(elemEffect(1)*1000._dp) .ne.  605_i4)  isgood = .false.
  if (nint(elemEffect(2)*1000._dp) .ne.  770_i4)  isgood = .false.
  if (nint(elemEffect(3)*1000._dp) .ne. 1100_i4)  isgood = .false.
  if (nint(elemEffect(4)*1000._dp) .ne. 2000_i4)  isgood = .false.

  if (isgood) then
     write(*,*) 'mo_elemeffects o.k.'
  else
     write(*,*) 'mo_elemeffects failed'
  end if

END PROGRAM elemeffects_test

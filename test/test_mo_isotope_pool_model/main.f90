program main

  use mo_kind,   only: dp, i4
  use mo_utils,  only: ne
  use mo_isotope_pool_model, only: isotope_pool_model

  implicit none

  integer(i4), parameter :: nn = 3 ! # of pools
  real(dp),    parameter :: vpdb = 0.0112372_dp

  ! input / output
  real(dp) :: dt
  real(dp), dimension(nn)    :: Ci, C, Ct, S, Rs, Si, beta, trash
  real(dp), dimension(nn,nn) :: F, alpha

  integer(i4) :: i, j, nstep

  logical :: isgood, allgood

  write(*,*) ''
  write(*,*) 'Test mo_isotope_pool_model'
  
  allgood = .true.
  
  ! -----------------------
  ! No fractionation - Sink

  isgood = .true.

  dt    = 0.023_dp
  C     = 1.2345_dp
  S     = 0.004_dp
  Rs    = vpdb ! ~alpha=1
  Si    = 0.003_dp
  beta  = 0._dp
  trash = 0._dp

  do i=1, nn
     do j=1, nn
        if (i == j) then
           F(i,j)     = 0._dp
           alpha(i,j) = 1.0_dp   ! sinks
        else if  (i < j) then    ! outgoing fluxes
           F(i,j)     = 0.001_dp
           alpha(i,j) = 1._dp
        else                     ! incoming fluxes
           F(i,j)     = 0.002_dp
           alpha(i,j) = 1._dp
        endif
     end do
  end do

  Ci = C * vpdb
  nstep = 20
  do i=1, nstep
     Ct = C
     C = C + sum(F, dim=1)*dt - sum(F, dim=2)*dt + S*dt - Si*dt
     ! subroutine isotope_pool_model(dt, Ci, C, F, S, Rs, Si, alpha, beta, trash)
     call isotope_pool_model(dt, Ci, Ct, F, S=S, Rs=Rs, Si=Si, alpha=alpha, trash=trash)
     ! print*, 'd ', (Ci/C/vpdb-1.)*1000.
  end do
  ! print*, 'Trash ', trash

  if (any(abs(Ci/C/vpdb-1.)*1000.>1e-9)) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model no frac sink o.k.'
  else
     write(*,*) 'mo_isotope_pool_model no frac sink failed!'
  endif
  
  ! -----------------------
  ! No fractionation - beta

  isgood = .true.

  dt    = 0.023_dp
  C     = 1.2345_dp
  S     = 0.004_dp
  Rs    = vpdb ! ~alpha=1
  Si    = 0.003_dp
  beta  = Si / C
  trash = 0._dp

  do i=1, nn
     do j=1, nn
        if (i == j) then
           F(i,j)     = 0._dp
           alpha(i,j) = 1.0_dp   ! sinks
        else if  (i < j) then    ! outgoing fluxes
           F(i,j)     = 0.001_dp
           alpha(i,j) = 1._dp
        else                     ! incoming fluxes
           F(i,j)     = 0.002_dp
           alpha(i,j) = 1._dp
        endif
     end do
  end do

  Ci = C * vpdb
  nstep = 20
  do i=1, nstep
     Ct = C
     C = C + sum(F, dim=1)*dt - sum(F, dim=2)*dt + S*dt - beta*C*dt
     ! subroutine isotope_pool_model(dt, Ci, C, F, S, Rs, Si, alpha, beta, trash)
     call isotope_pool_model(dt, Ci, Ct, F, S=S, Rs=Rs, beta=beta, alpha=alpha, trash=trash)
     ! print*, 'd ', (Ci/C/vpdb-1.)*1000.
  end do
  ! print*, 'Trash ', trash

  if (any(abs(Ci/C/vpdb-1.)*1000.>1e-9)) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model no frac beta o.k.'
  else
     write(*,*) 'mo_isotope_pool_model no frac beta failed!'
  endif
  
  ! -----------------------
  ! One pool, different source composition

  isgood = .true.

  dt     = 23._dp
  C(1)   = 1.2345_dp
  C(2:)  = 0._dp
  S(1)   = 0.004_dp
  S(2:)  = 0._dp
  Rs(1)  = vpdb * 1.1_dp
  Rs(2:) = 0._dp
  Si(1)  = 0.003_dp
  Si(2:) = 0._dp
  beta   = 1._dp
  where (C > 0.0_dp) beta  = Si / C
  trash  = 0._dp
  F      = 0._dp
  alpha  = 1._dp

  Ci = C * vpdb
  nstep = 2000
  do i=1, nstep
     Ct = C
     C = C + sum(F, dim=1)*dt - sum(F, dim=2)*dt + S*dt - beta*C*dt
     ! subroutine isotope_pool_model(dt, Ci, C, F, S, Rs, Si, alpha, beta, trash)
     call isotope_pool_model(dt, Ci, Ct, F, S=S, Rs=Rs, beta=beta, alpha=alpha, trash=trash)
     !print*, 'd ', (Ci(1)/C(1)/vpdb-1.)*1000.
  end do
  ! print*, 'Trash ', trash
  ! write(*,*) 'delta_source ', (Rs(1)/vpdb-1.)*1000.
  ! write(*,*) 'delta end    ', (Ci(1)/C(1)/vpdb-1.)*1000.
  
  if (abs(Ci(1)/C(1)/vpdb - Rs(1)/vpdb)*1000.>1e-9) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model one pool source o.k.'
  else
     write(*,*) 'mo_isotope_pool_model one pool source failed!'
  endif
  
  ! -----------------------
  ! Two pool, source in 1st pool, everything shuffled to 2nd pool

  isgood = .true.

  dt      = 23._dp
  C(1:2)  = 1.2345_dp
  C(3:)   = 0._dp
  S(1)    = 0.04_dp
  S(2:)   = 0._dp
  Rs(1)   = vpdb * 1.1_dp
  Rs(2:)  = 0._dp ! only source in 1st pool
  Si(1)   = 0._dp
  Si(2)   = S(1)  ! remove in 2nd as much as is coming from 1st pool
  Si(3:)  = 0._dp
  beta    = 1._dp
  where (C > 0.0_dp) beta  = Si / C
  trash   = 0._dp
  F       = 0._dp
  F(1,2)  = S(1)  ! shuffle to 2nd pool
  alpha   = 1._dp

  Ci = C * vpdb
  nstep = 200
  do i=1, nstep
     Ct = C
     C = C + sum(F, dim=1)*dt - sum(F, dim=2)*dt + S*dt - beta*C*dt
     ! subroutine isotope_pool_model(dt, Ci, C, F, S, Rs, Si, alpha, beta, trash)
     call isotope_pool_model(dt, Ci, Ct, F, S=S, Rs=Rs, beta=beta, alpha=alpha, trash=trash)
     ! print*, 'd ', (Ci(1:2)/C(1:2)/vpdb-1.)*1000.
  end do
  ! print*, 'Trash ', trash
  ! write(*,*) 'delta_source ', (Rs(1)/vpdb-1.)*1000.
  ! write(*,*) 'delta end    ', (Ci(2)/C(2)/vpdb-1.)*1000.
  
  if (abs(Ci(2)/C(2)/vpdb - Rs(1)/vpdb)*1000.>1e-9) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model two pool shuffle o.k.'
  else
     write(*,*) 'mo_isotope_pool_model two pool shuffle failed!'
  endif

  ! ----------------------
  
  if (allgood) then
     write(*,*) 'mo_isotope_pool_model o.k.'
  else
     write(*,*) 'mo_isotope_pool_model failed!'
  endif

END PROGRAM main

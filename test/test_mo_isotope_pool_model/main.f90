program main

  use mo_kind,   only: dp, i4
  use mo_isotope_pool_model, only: isotope_pool_model, isotope_luc_model

  implicit none

  integer(i4), parameter :: nland = 10 ! # of pools
  integer(i4), parameter :: nn    = 3 ! # of pools
  real(dp),    parameter :: vpdb  = 0.0112372_dp

  ! input / output
  real(dp) :: dt, Citot
  real(dp), dimension(nn)    :: Ci, C, Ct, S, Rs, Si, beta, trash, A, At
  real(dp), dimension(nn,nn) :: F, alpha, dA
  real(dp), dimension(nn,nland)    :: Ci2, C2, Ct2, S2, Rs2, Si2, beta2, trash2, A2, At2
  real(dp), dimension(nn,nn,nland) :: F2, alpha2, dA2

  integer(i4) :: i, j, nstep

  logical :: isgood, allgood

  write(*,*) ''
  write(*,*) 'Test mo_isotope_pool_model'
  
  allgood = .true.
  
  ! -----------------------
  ! isotope_pool_model - No fractionation - Sink

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
     write(*,*) 'mo_isotope_pool_model pool no frac sink o.k.'
  else
     write(*,*) 'mo_isotope_pool_model pool no frac sink failed!'
  endif
  
  ! -----------------------
  ! isotope_pool_model - No fractionation - beta

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
     write(*,*) 'mo_isotope_pool_model pool no frac beta o.k.'
  else
     write(*,*) 'mo_isotope_pool_model pool no frac beta failed!'
  endif
  
  ! -----------------------
  ! isotope_pool_model - One pool, different source composition

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
     write(*,*) 'mo_isotope_pool_model pool one pool source o.k.'
  else
     write(*,*) 'mo_isotope_pool_model pool one pool source failed!'
  endif
  
  ! -----------------------
  ! isotope_pool_model - Two pool, source in 1st pool, everything shuffled to 2nd pool

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
     write(*,*) 'mo_isotope_pool_model pool two pool shuffle o.k.'
  else
     write(*,*) 'mo_isotope_pool_model pool two pool shuffle failed!'
  endif
  
  ! -----------------------
  ! isotope_pool_model_2d - No fractionation - Sink

  isgood = .true.

  dt     = 0.023_dp
  C2     = 1.2345_dp
  S2     = 0.004_dp
  Rs2    = vpdb ! ~alpha=1
  Si2    = 0.003_dp
  beta2  = 0._dp
  trash2 = 0._dp

  do i=1, nn
     do j=1, nn
        if (i == j) then
           F2(i,j,:)     = 0._dp
           alpha2(i,j,:) = 1.0_dp   ! sinks
        else if  (i < j) then    ! outgoing fluxes
           F2(i,j,:)     = 0.001_dp
           alpha2(i,j,:) = 1._dp
        else                     ! incoming fluxes
           F2(i,j,:)     = 0.002_dp
           alpha2(i,j,:) = 1._dp
        endif
     end do
  end do

  Ci2 = C2 * vpdb
  nstep = 20
  do i=1, nstep
     Ct2 = C2
     C2 = C2 + sum(F2, dim=1)*dt - sum(F2, dim=2)*dt + S2*dt - Si2*dt
     call isotope_pool_model(dt, Ci2, Ct2, F2, S=S2, Rs=Rs2, Si=Si2, alpha=alpha2, trash=trash2)
     ! print*, 'd ', (Ci/C/vpdb-1.)*1000.
  end do
  ! print*, 'Trash ', trash

  if (any(abs(Ci2/C2/vpdb-1.)*1000.>1e-9)) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model 2d pool no frac sink o.k.'
  else
     write(*,*) 'mo_isotope_pool_model 2d pool no frac sink failed!'
  endif
  
  ! -----------------------
  ! isotope_pool_model2d - Two pool, source in 1st pool, everything shuffled to 2nd pool

  isgood = .true.

  dt      = 23._dp
  C2(1:2,:nland/2)   = 1.2345_dp
  C2(1:2,nland/2+1:) = 2.2345_dp
  C2(3:,:)   = 0._dp
  S2(1,:)    = 0.04_dp
  S2(2:,:)   = 0._dp
  Rs2(1,:)   = vpdb * 1.1_dp
  Rs2(2:,:)  = 0._dp ! only source in 1st pool
  Si2(1,:)   = 0._dp
  Si2(2,:)   = S(1)  ! remove in 2nd as much as is coming from 1st pool
  Si2(3:,:)  = 0._dp
  beta2    = 1._dp
  where (C2 > 0.0_dp) beta2  = Si2 / C2
  trash2   = 0._dp
  F2       = 0._dp
  F2(1,2,:) = S2(1,:)  ! shuffle to 2nd pool
  alpha2   = 1._dp

  Ci2 = C2 * vpdb
  nstep = 200
  do i=1, nstep
     Ct2 = C2
     C2 = C2 + sum(F2, dim=1)*dt - sum(F2, dim=2)*dt + S2*dt - beta2*C2*dt
     call isotope_pool_model(dt, Ci2, Ct2, F2, S=S2, Rs=Rs2, beta=beta2, alpha=alpha2, trash=trash2)
     ! print*, 'd ', (Ci(1:2)/C(1:2)/vpdb-1.)*1000.
  end do
  ! print*, 'Trash ', trash
  ! write(*,*) 'delta_source ', (Rs(1)/vpdb-1.)*1000.
  ! write(*,*) 'delta end    ', (Ci(2)/C(2)/vpdb-1.)*1000.
  
  if (any(abs(Ci2(2,:)/C2(2,:)/vpdb - Rs2(1,:)/vpdb)*1000.>1e-9)) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model 2d pool two pool shuffle o.k.'
  else
     write(*,*) 'mo_isotope_pool_model 2d pool two pool shuffle failed!'
  endif
  
  ! -----------------------
  ! isotope_luc_model - Same isotopic composition

  isgood = .true.

  dt    = 0.023_dp
  C     = 1.2345_dp
  Rs    = vpdb ! ~alpha=1
  Ci    = Rs*C
  A     = 0.1_dp
  trash = 0._dp

  do i=1, nn
     do j=1, nn
        if (i == j) then
           dA(i,j) = 0._dp
        else if  (i < j) then    ! outgoing fluxes
           if (i==1) then
              dA(i,j) = 0.01_dp
           else
              dA(i,j) = 0.02_dp
           endif
        else                     ! incoming fluxes
           if (j==1) then
              dA(i,j) = 0.01_dp
           else
              dA(i,j) = 0.02_dp
           endif
        endif
     end do
  end do

  nstep = 20
  do i=1, nstep
     Ct = C
     At = A
     A = At - sum(dA,2) + sum(dA,1)
     C = (Ct*(At-sum(dA,2)) + sum(dA*spread(Ct,dim=2,ncopies=nn),1)) / A
     call isotope_luc_model(Ci, At, dA, trash=trash)
  end do

  if (any(abs(Ci/C/vpdb-1.)*1000.>1e-9)) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model luc Ci all same o.k.'
  else
     write(*,*) 'mo_isotope_pool_model luc Ci all same failed!'
  endif
  
  ! -----------------------
  ! isotope_luc_model - Same isotopic composition - empty land-use class

  isgood = .true.

  dt    = 0.023_dp
  C     = 1.2345_dp
  Rs    = vpdb ! R
  Ci    = Rs*C
  A     = 0.1_dp
  trash = 0._dp

  dA(1,:) = (/ 0._dp, 0.01_dp, 0.02_dp /)
  dA(2,:) = (/ 0.02_dp, 0._dp, 0.01_dp /)
  dA(3,:) = (/ 0.01_dp, 0.01_dp, 0._dp /)

  nstep = 20
  do i=1, nstep
     Ct = C
     At = A
     A = At - sum(dA,2) + sum(dA,1)
     if (A(1) < 0._dp) then
        dA(1,2) = 0._dp ! set outgoing fluxes = 0.
        dA(1,3) = 0._dp
        A = At - sum(dA,2) + sum(dA,1)
     endif
     if (A(2) < 0._dp) then
        dA(2,1) = 0._dp 
        dA(2,3) = 0._dp
        A = At - sum(dA,2) + sum(dA,1)
     endif
     if (A(3) < 0._dp) then
        dA(3,1) = 0._dp 
        dA(3,2) = 0._dp
        A = At - sum(dA,2) + sum(dA,1)
     endif
     C = (Ct*(At-sum(dA,2)) + sum(dA*spread(Ct,dim=2,ncopies=nn),1)) / A
     call isotope_luc_model(Ci, At, dA, trash=trash)
  end do

  if (any(abs(Ci/C/vpdb-1.)*1000.>1e-9)) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model luc empty class o.k.'
  else
     write(*,*) 'mo_isotope_pool_model luc empty class failed!'
  endif
  
  ! -----------------------
  ! isotope_luc_model - Shuffle to last pool

  isgood = .true.

  dt    = 0.023_dp
  C     = 1.2345_dp
  forall(i=1:nn) Rs(i) = (1._dp+real(i,dp)/10._dp) * vpdb ! R
  Ci    = Rs*C
  A     = 0.1_dp
  trash = 0._dp
  Citot = sum(A*Ci)/sum(A)

  dA(1,:) = (/ 0._dp, 0.0_dp, 0.02_dp /) ! 1->3
  dA(2,:) = (/ 0.0_dp, 0._dp, 0.01_dp /) ! 2->3
  dA(3,:) = 0._dp

  nstep = 20
  do i=1, nstep
     Ct = C
     At = A
     A = At - sum(dA,2) + sum(dA,1)
     if (A(1) < 0._dp) then
        dA(1,2) = 0._dp ! set outgoing fluxes = 0.
        dA(1,3) = At(1)
        A = At - sum(dA,2) + sum(dA,1)
     endif
     if (A(2) < 0._dp) then
        dA(2,1) = 0._dp 
        dA(2,3) = At(2)
        A = At - sum(dA,2) + sum(dA,1)
     endif
     C = (Ct*(At-sum(dA,2)) + sum(dA*spread(Ct,dim=2,ncopies=nn),1)) / A
     call isotope_luc_model(Ci, At, dA, trash=trash)
  end do

  if (abs(Ci(3) - Citot) > epsilon(1.0_dp)) isgood = .false.
  if (abs(sum(A*Ci)/sum(A) - Citot) > epsilon(1.0_dp)) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model luc shuffle o.k.'
  else
     write(*,*) 'mo_isotope_pool_model luc shuffle failed!'
  endif
  
  ! -----------------------
  ! isotope_luc_model - Same isotopic composition, passing C

  isgood = .true.

  dt    = 0.023_dp
  C     = 1.2345_dp
  Rs    = vpdb ! ~alpha=1
  Ci    = Rs*C
  A     = 0.1_dp
  trash = 0._dp

  do i=1, nn
     do j=1, nn
        if (i == j) then
           dA(i,j) = 0._dp
        else if  (i < j) then    ! outgoing fluxes
           if (i==1) then
              dA(i,j) = 0.01_dp
           else
              dA(i,j) = 0.02_dp
           endif
        else                     ! incoming fluxes
           if (j==1) then
              dA(i,j) = 0.01_dp
           else
              dA(i,j) = 0.02_dp
           endif
        endif
     end do
  end do

  nstep = 20
  do i=1, nstep
     Ct = C
     At = A
     A = At - sum(dA,2) + sum(dA,1)
     C = (Ct*(At-sum(dA,2)) + sum(dA*spread(Ct,dim=2,ncopies=nn),1)) / A
     call isotope_luc_model(Ci, At, dA, Ct, trash=trash)
  end do

  if (any(abs(Ci/C/vpdb-1.)*1000.>1e-9)) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model luc Ci all same C o.k.'
  else
     write(*,*) 'mo_isotope_pool_model luc Ci all same C failed!'
  endif
  
  ! -----------------------
  ! isotope_luc_model - Same isotopic composition - empty land-use class, passing C

  isgood = .true.

  dt    = 0.023_dp
  C     = 1.2345_dp
  Rs    = vpdb ! R
  Ci    = Rs*C
  A     = 0.1_dp
  trash = 0._dp

  dA(1,:) = (/ 0._dp, 0.01_dp, 0.02_dp /)
  dA(2,:) = (/ 0.02_dp, 0._dp, 0.01_dp /)
  dA(3,:) = (/ 0.01_dp, 0.01_dp, 0._dp /)

  nstep = 20
  do i=1, nstep
     Ct = C
     At = A
     A = At - sum(dA,2) + sum(dA,1)
     if (A(1) < 0._dp) then
        dA(1,2) = 0._dp ! set outgoing fluxes = 0.
        dA(1,3) = 0._dp
        A = At - sum(dA,2) + sum(dA,1)
     endif
     if (A(2) < 0._dp) then
        dA(2,1) = 0._dp 
        dA(2,3) = 0._dp
        A = At - sum(dA,2) + sum(dA,1)
     endif
     if (A(3) < 0._dp) then
        dA(3,1) = 0._dp 
        dA(3,2) = 0._dp
        A = At - sum(dA,2) + sum(dA,1)
     endif
     C = (Ct*(At-sum(dA,2)) + sum(dA*spread(Ct,dim=2,ncopies=nn),1)) / A
     call isotope_luc_model(Ci, At, dA, Ct, trash=trash)
  end do

  if (any(abs(Ci/C/vpdb-1.)*1000.>1e-9)) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model luc empty class C o.k.'
  else
     write(*,*) 'mo_isotope_pool_model luc empty class C failed!'
  endif
  
  ! -----------------------
  ! isotope_luc_model - Shuffle to last pool, passing C

  isgood = .true.

  dt    = 0.023_dp
  C     = 1.2345_dp
  forall(i=1:nn) Rs(i) = (1._dp+real(i,dp)/10._dp) * vpdb ! R
  Ci    = Rs*C
  A     = 0.1_dp
  trash = 0._dp
  Citot = sum(A*Ci)/sum(A)

  dA(1,:) = (/ 0._dp, 0.0_dp, 0.02_dp /) ! 1->3
  dA(2,:) = (/ 0.0_dp, 0._dp, 0.01_dp /) ! 2->3
  dA(3,:) = 0._dp

  nstep = 20
  do i=1, nstep
     Ct = C
     At = A
     A = At - sum(dA,2) + sum(dA,1)
     if (A(1) < 0._dp) then
        dA(1,2) = 0._dp ! set outgoing fluxes = 0.
        dA(1,3) = At(1)
        A = At - sum(dA,2) + sum(dA,1)
     endif
     if (A(2) < 0._dp) then
        dA(2,1) = 0._dp 
        dA(2,3) = At(2)
        A = At - sum(dA,2) + sum(dA,1)
     endif
     where (A > 0._dp)
        C = (Ct*(At-sum(dA,2)) + sum(dA*spread(Ct,dim=2,ncopies=nn),1)) / A
     elsewhere
        C = 0._dp
     endwhere
     call isotope_luc_model(Ci, At, dA, Ct, trash=trash)
  end do

  if (abs(Ci(3) - Citot) > epsilon(1.0_dp)) isgood = .false.
  if (abs(sum(A*Ci)/sum(A) - Citot) > epsilon(1.0_dp)) isgood = .false.
  
  allgood = allgood .and. isgood
  
  if (isgood) then
     write(*,*) 'mo_isotope_pool_model luc shuffle C o.k.'
  else
     write(*,*) 'mo_isotope_pool_model luc shuffle C failed!'
  endif

  ! ----------------------
  
  if (allgood) then
     write(*,*) 'mo_isotope_pool_model o.k.'
  else
     write(*,*) 'mo_isotope_pool_model failed!'
  endif

END PROGRAM main

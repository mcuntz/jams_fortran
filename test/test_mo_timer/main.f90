program main

  use mo_kind,  only: i4, dp
  use mo_timer, only: max_timers, timers_init, timer_clear, timer_start, timer_stop, timer_get, &
       max_timers, cycles_max, clock_rate, cycles1, cycles2, cputime, status, timer_check, timer_print

  implicit none

  integer(i4) :: i
  real(dp), dimension(1000000) :: r

  ! Initialise timers
  call timers_init
  do i=1, max_timers
     call timer_clear(i)
  end do
  !
  ! Do something
  call timer_start(1)
  ! r = huge(1.0_dp)
  r = 1.0_dp/epsilon(1.0_dp)
  do i=1, 10
     r = r**(1.0_dp/real(i,dp))
  end do
  call timer_check(1)
  call timer_print(1)
  call timer_stop(1)
  !
  ! Do something else
  call timer_start(2)
  do i=1, 10
     r = exp(sin(r))
  end do
  call timer_check(2)
  call timer_print(2)
  call timer_stop(2)
  !
  ! Wrte out timings
  write(*,*) ''
  write(*,*) 'Timings [s]'
  write(*,*) 'Something done: ', timer_get(1)
  write(*,*) 'Something else: ', timer_get(2)
  !
  ! test a bit more
  write(*,*) 'max_timers ', max_timers
  write(*,*) 'cycles_max ', cycles_max
  write(*,*) 'clock_rate ', clock_rate
  write(*,*) 'cycles1    ', cycles1   
  write(*,*) 'cycles2    ', cycles2   
  write(*,*) 'cputime    ', cputime   
  write(*,*) 'status     ', status
  do i=1, 2
     call timer_check(i)
  end do
  do i=1, 2
     call timer_print(i)
  end do
  !
  write(*,*) 'mo_timer single precision o.k.'
  !
end program main

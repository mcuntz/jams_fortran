program main
  ! ------------------------------------------------------------------------------
  !
  ! program remap-con
  !     performs conservative remaping based on given weights
  !
  ! author: Stephan Thober, Matthias Cuntz
  !
  ! created: 10 Aug 2012
  !
  ! ------------------------------------------------------------------------------

  use mo_kind,      only: i4, i8, sp, dp
  use mo_ansi_colors, only: color, c_green
  use mo_timer,     only: max_timers, timers_init, timer_clear, timer_start, timer_stop, timer_get
  use mo_readdata,  only: readdata
  use mo_mainvar,   only: src_data, dst_data, src_add, dst_add, wts, nodata
  use mo_remap,     only: remap
  use mo_setnc,     only: setnc
  use mo_ncwrite,   only: create_netcdf, write_static_netcdf, write_dynamic_netcdf, close_netcdf, V
  use mo_mainvar,   only: time, outfile

  implicit none

  integer(i4) :: hnc, i
  real(dp), dimension(:), allocatable :: src_1d
  real(dp), dimension(:), allocatable :: dst_1d

  write(*,*) 'Start'
  !
  ! Ini timers
  call timers_init
  do i=1, max_timers
     call timer_clear(int(i, i8))
  end do

  !
  ! Read data
  call timer_start(1_i8)
  call readdata
  call timer_stop(1_i8)
  write(*,*) 'readdata...OK'

  !
  ! Prepare output
  call timer_start(3_i8)
  call setnc
  call timer_stop(3_i8)
  write(*,*) 'setnc...OK'

  call timer_start(4_i8)
  call create_netcdf(trim(outfile), hnc)
  call timer_stop(4_i8)
  write(*,*) 'create_netcdf...OK'

  call timer_start(5_i8)
  call write_static_netcdf(hnc)
  call timer_stop(5_i8)
  write(*,*) 'write_static_netcdf...OK'

  !
  ! Remap
  call timer_start(2_i8)
  allocate(src_1d(size(src_data,1)*size(src_data,2)))
  allocate(dst_1d(size(dst_data,1)*size(dst_data,2)))
  do i=1, size(src_data,3)
     src_1d = reshape(src_data(:,:,i),shape(src_1d))                                   ! 2d to 1d
     call remap(dst_1d, wts, dst_add, src_add, src_1d, miss=real(nodata,dp))           ! remap (in dp)
     dst_data(:,:,i) = reshape(real(dst_1d,sp), (/size(dst_data,1),size(dst_data,2)/)) ! 1d to 2d
  end do
  deallocate(src_1d)
  deallocate(dst_1d)
  call timer_stop(2_i8)
  write(*,*) 'remap_con...OK'

  !
  ! Write out results
  call timer_start(6_i8)
  do i=1, size(src_data,3)
     V(3)%G0_d => time(i)
     V(4)%G2_f => dst_data(:,:,i)
     call write_dynamic_netcdf(hnc,i)
  enddo
  call timer_stop(6_i8)
  write(*,*) 'write_dynamic_netcdf...OK'

  call close_netcdf(hnc)

  write(*,*) ''
  write(*,*) 'Timings [s]'
  write(*,*) 'Read data:    ', timer_get(1_i8)
  write(*,*) 'Remap:        ', timer_get(2_i8)
  write(*,*) 'SetNC:        ', timer_get(3_i8)
  write(*,*) 'CreateNC:     ', timer_get(4_i8)
  write(*,*) 'WriteStatic:  ', timer_get(5_i8)
  write(*,*) 'WriteDynamic: ', timer_get(6_i8)

  write(*,*) ''
  write(*,*) 'Done!'

  write(*,*) 'mo_remap double precision ', color('o.k.', c_green)


end program main

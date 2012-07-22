
Module ModLine
  !
  ! robust fit of the line
  !
  use mo_kind, only: dp

  implicit none

  ! debug mode
  logical, parameter, private :: debug = .false.

  integer, private :: n
  real(dp), private, pointer :: x(:), y(:)
  real(dp), private, allocatable, dimension(:,:) :: jacmin

contains

  subroutine line(ndata,xdata,ydata,t,dt,s0)

    use mo_minpack, only: hybrj, qrfac, covar, qrfac, qform, qrinv

    implicit none

    integer, parameter :: npar = 2
    integer, parameter :: lwa = (npar*(3*npar+13))/2
    integer, parameter :: lr = (npar*(npar+1))/2

    integer :: ndata
    real(dp), target :: xdata(:),ydata(:)
    real(dp) :: t(2),dt(2),s0

    integer :: i, j, info, ipvt(2), iflag, nfev, njev
    real(dp) :: tol, fvec(2), fjac(2,2),h(2,2),cov(2,2)
    real(dp) :: fjacdiag(2),acnorm(2),sum1,sum2,temp,temp2
    real(dp) :: v(2,2),w(2)
    real(dp), dimension(npar) :: diag,wa1,wa2,wa3,wa4,qtf
    real(dp), dimension(lr) :: r

    LOGICAL :: isgood

! initialization
    n = ndata
    x => xdata
    y => ydata
    allocate(jacmin(2,2))

    
! 0th estimator, it get value of parameters by minimization of the absolute 
! deviations, uses simplex method
    t = 0
    dt = 1.0
    t = (/29.223, -1.913 /)
    t = 0
    dt = 0

! iterate to get better solution
    tol = sqrt(epsilon(tol))
    info = 0
!    call hybrj1(minfun,2,t,fvec,fjac,2,tol,info,wa,(2*(2+13))/2)
!    write(*,*) 'hybrj1 finished with code',info
    diag = 1.0
    call hybrj(minfun,2,t,fvec,fjac,2,tol,200*(2+1),diag,2,100.0_dp,1, &
         info,nfev,njev,r,lr,qtf,wa1,wa2,wa3,wa4)
    if (debug) write(*,*) 'hybrj finished with code',info

! form inverse matrix
    if (debug) write(*,*) "hybr:"
    if (debug) write(*,'(2F15.5)') fjac

    fjac = jacmin
    call qrfac(2,2,fjac,2,.true.,ipvt,2,fjacdiag,acnorm,w)
    do i = 1,2
       fjac(i,i) = fjacdiag(i)
    end do

    call covar(2,fjac,2,ipvt,epsilon(fjac),w)

    if (debug) write(*,*) " covariance matrix = "
    if (debug) write(*,'(2f15.5)') fjac


!    call covar(2,fjac,2,ipvt,epsilon(fjac),w)


    iflag = 2
    call minfun(2,t,fvec,fjac,2,iflag)
    call qrfac(2,2,fjac,2,.false.,ipvt,2,fjacdiag,acnorm,w) 
    if (debug) write(*,*) acnorm
    if (debug) write(*,'(2F15.5)') fjac
    call qform(2,2,fjac,2,w)
    v = fjac
    w = fjacdiag
    if (debug) write(*,*) ipvt
    if (debug) write(*,"(A,2(2F15.5/))") "eigean. vectors", v
       
    if (debug) write(*,*) "eigrn values"
    if (debug) write(*,'(2F15.5)') w

    call qrinv(2,2,v,w,2,fjac)
    if (debug) write(*,*) "qrinv:"
    if (debug) write(*,'(2F15.5)') fjac

!    goto 666
! inverse matrix
    iflag = 2
    call minfun(2,t,fvec,fjac,2,iflag)
    call qrfac(2,2,fjac,2,.true.,ipvt,2,fjacdiag,acnorm,w) 
    do i = 1,2
      h(i,i) = 1.0/fjacdiag(i)
    enddo
    if (debug) write(*,*) h
    do j = 2,1,-1
       do i = j-1,1,-1
          h(i,j) = -sum(fjac(i,i+1:j)*h(i+1:j,j))/fjacdiag(i)
       enddo
    enddo
    if (debug) write(*,*) "inverse = "
    if (debug) write(*,'(2F15.5)') h

! compute inverse matrix
!    fjac(1,1:2) = (/ 1.0, -1.0/)
!    fjac(2,1:2) = (/ -1.0, 1.0 /)
    iflag = 2
    call minfun(2,t,fvec,fjac,2,iflag)
    if (debug) write(*,*) " input jacobian = "
    if (debug) write(*,'(2F15.5)') fjac
    call qrfac(2,2,fjac,2,.true.,ipvt,2,fjacdiag,acnorm,w)

    if (debug) write(*,*) fjacdiag
    if (debug) write(*,*) acnorm
    if (debug) write(*,*) ipvt
!    ipvt = (/ 1,2/)
    if (debug) write(*,'(a,2F15.5)') " inverse jacobian = ",fjac(1,1:2)
    if (debug) write(*,'(a,2F15.5)') "                    ",fjac(2,1:2)

    do i = 1,2
       fjac(i,i) = fjacdiag(i)
    end do

    call covar(2,fjac,2,ipvt,epsilon(fjac),w)

    if (debug) write(*,*) " covariance matrix = "
    if (debug) write(*,'(2f15.5)') fjac


666 continue

    sum1 = 0.0
    sum2 = 0.0
    s0 = 0
    do i = 1, n
       temp = (y(i) - (t(1) + t(2)*x(i)))
       temp2 = temp**2
       sum1 = sum1 + temp2
       sum2 = sum2 + 1.0
       s0 = s0 + temp**2
       if (debug) write(*,'(5f15.5)') x(i), y(i),temp,sum1, s0
    enddo
    if (debug) write(*,*) sum1,sum2**2
!    s0 = sum1
    cov = s0*fjac/(n - 2)!*(sum2/n)**2
    write(*,*) "covariance:"
    write(*,'(2f15.5)') cov
    do i = 1,2
       dt(i) = sqrt(cov(i,i))
    enddo

   
    if( debug ) then
       write(*,*) "Exit for "
       if( info == 0 ) then
          write(*,*) " improper input parameters."
       elseif( info == 1 )then
          write(*,*) " algorithm estimates that the relative error between x and the solution is at most tol."
       elseif( info == 2 )then
          write(*,*) " number of calls to fcn with iflag = 1 has reached 100*(n+1)."
       elseif( info == 3 ) then
          write(*,*) " tol is too small. no further improvement in the approximate solution x is possible."
       elseif( info == 4 )then
          write(*,*) " iteration is not making good progress."
       endif

       write(*,*) "Solution"
       write(*,*) " info = ", info
       write(*,*) " solution = ",t," with dev = ",dt
       write(*,*) " residual sum s0 = ",s0
       write(*,*) " standard deviation = ",sqrt(s0/(n - 2))
       write(*,*) " No. of data = ",n
       
    else
       write(*,*) "Solution"
       write(*,*) t(1)," +- ",dt(1)
       write(*,*) t(2)," +- ",dt(2)
       write(*,*) " residual sum: ",s0
       write(*,*) " No. of data:  ",n
    endif
    isgood = .True.
    isgood = isgood .and. (anint(1000._dp*t(1)) == 29223._dp)
    isgood = isgood .and. (anint(1000._dp*t(2)) == -1913._dp)
    isgood = isgood .and. (anint(1000._dp*dt(1)) == 1278._dp)
    isgood = isgood .and. (anint(1000._dp*dt(2)) == 122._dp)

    s0 = 0.0
    do i = 1, n
       sum1 = y(i) - (t(1) + t(2)*x(i))
       s0 =  s0 + sum1**2
    enddo
    if (debug) write(*,*) s0

    if (isgood) then
       write(*,*) 'mo_minpack double precision o.k.'
    else
       write(*,*) 'mo_minpack double precision failed!'
    endif

  end subroutine line
!
!--------------------------------------------------------------------------
!
  subroutine minfun(nt,t,fvec,fjac,ldjac,iflag)
    
    implicit none

    ! integer, intent(in) :: nt, ldjac
    ! integer, intent(inout) :: iflag
    ! real(dp), intent(in) :: t(nt)
    ! real(dp), intent(inout) :: fvec(nt), fjac(ldjac,nt)
    integer :: nt, ldjac
    integer :: iflag
    real(dp) :: t(nt)
    real(dp) :: fvec(nt), fjac(ldjac,nt)

    integer :: i
    real(dp) :: r, sa,sb

    if( ldjac /= 2 .or. nt /= 2 ) stop

    if( iflag == 0 .and. debug )then
       write(*,'(a,2F15.5,a,2f15.5)') "minpack: par = ",t," sum = ", &
            sum(abs(y(1:n) - (t(1) + t(2)*x(1:n)))**2)
    endif

    if( iflag == 1 )then
!       fvec = 0.0_dp
       sa = 0.0; sb = 0.0;
       do i = 1, n
          r = y(i) - (t(1) + t(2)*x(i))
!          write(*,*) r,t(1) + t(2)*x(i),y(i),x(i),t
          sa = sa - r
          sb = sb - r*x(i)
!          rp = r
!          fvec(1) = fvec(1) - rp*1.0_dp
!          fvec(2) = fvec(2) - rp*x(i)
!          if(debug) write(*,'(5f15.5)') x(i),y(i),r,sa,sb
       enddo
!       fvec = fvec
       fvec = (/sa,sb/)
!       fvec = 0.0
       if( debug ) write(*,'(a,2f15.5)') " vector of residuals = ",fvec
    endif

    if( iflag == 2 )then
!       fjac = 0.0_dp
       sa = 0.0; sb = 0.0;
       do i = 1, n
!          r = y(i) - (t(1) + t(2)*x(i))
          sa = sa - x(i)
          sb = sb - x(i)**2
!          rp = 1.0_dp
!          fjac(1,1) = fjac(1,1) - rp*1.0_dp
!          fjac(1,2) = fjac(1,2) - rp*x(i)
!          fjac(2,2) = fjac(2,2) - rp*x(i)**2
       enddo
!       fjac(2,1) = fjac(1,2)
!       fjac = - fjac
       fjac(1,:) = (/real(n,kind=dp),-sa/)
       fjac(2,:) = (/-sa,-sb/)
       if( debug ) write(*,'(a,4f15.5)') " jacobian = ",fjac
       jacmin = fjac
    endif

  end subroutine minfun
!
!--------------------------------------------------------------------------
!
end Module ModLine

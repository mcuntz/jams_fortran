module mo_minpack

  ! F90 interface to minpack from munipack http://munipack.physics.muni.cz/

  ! License
  ! -------
  ! This file is part of the JAMS Fortran library.

  ! The JAMS Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The JAMS Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the JAMS Fortran library (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2012 Matthias Cuntz

  use mo_kind, only: dp

  implicit none

  ! set to .true. for debug prints
  logical, parameter, private :: dbg = .false.

  interface

     subroutine chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err)
       use mo_kind, only: dp
       implicit none
       integer, intent(in) :: m,n,ldfjac,mode
       real(dp) :: x(n),fvec(m),fjac(ldfjac,n),xp(n),fvecp(m),err(m)
     end subroutine chkder

     subroutine covar(n,r,ldr,ipvt,tol,wa)
       use mo_kind, only: dp
       implicit none
       integer, intent(in) ::  n,ldr
       integer::  ipvt(n)
       real(dp) :: tol
       real(dp) :: r(ldr,n),wa(n)
     end subroutine covar

     subroutine dmchar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp, &
          maxexp,eps,epsneg,xmin,xmax)
       use mo_kind, only: dp
       implicit none
       integer :: ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,maxexp
       real(dp) :: eps,epsneg,xmin,xmax
     end subroutine dmchar

     subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
       use mo_kind, only: dp
       implicit none
       integer :: n,lr
       real(dp) :: delta
       real(dp) :: r(lr),diag(n),qtb(n),x(n),wa1(n),wa2(n)
     end subroutine dogleg

     function dpmpar(i)
       use mo_kind, only: dp
       implicit none
       real(dp) :: dpmpar
       integer, intent(in) :: i
     end function dpmpar

     function enorm(n,x)
       use mo_kind, only: dp
       implicit none
       integer, intent(in) :: n
       real(dp) :: enorm,x(n)
     end function enorm

     subroutine errjac(n,x,fjac,ldfjac,nprob)
       use mo_kind, only: dp
       implicit none
       integer :: n,ldfjac,nprob
       real(dp) :: x(n),fjac(ldfjac,n)
     end subroutine errjac

     subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,wa2)
       use mo_kind, only: dp
       implicit none
       integer :: n,ldfjac,iflag,ml,mu
       real(dp) :: epsfcn
       real(dp) :: x(n),fvec(n),fjac(ldfjac,n),wa1(n),wa2(n)
       interface
          subroutine fcn(n,x,fvec,iflag)
            use mo_kind, only: dp
            implicit none
            integer :: n,iflag
            real(dp) :: x(n),fvec(n)
          end subroutine fcn
       end interface
     end subroutine fdjac1

     subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,ldfjac,iflag
       real(dp) :: epsfcn
       real(dp) :: x(n),fvec(m),fjac(ldfjac,n),wa(m)
       interface
          subroutine fcn(m,n,x,fvec,iflag)
            use mo_kind, only: dp
            implicit none
            integer :: m,n,iflag
            real(dp) :: x(n),fvec(m)
          end subroutine fcn
       end interface
     end subroutine fdjac2

     subroutine grdfcn(n,x,g,nprob)
       use mo_kind, only: dp
       implicit none
       integer :: n,nprob
       real(dp) :: x(n),g(n)
     end subroutine grdfcn

     subroutine hesfcn(n,x,h,ldh,nprob)
       use mo_kind, only: dp
       implicit none
       integer :: n,ldh,nprob
       real(dp) :: x(n),h(ldh,n)
     end subroutine hesfcn

     subroutine initpt(n,x,nprob,factor)
       use mo_kind, only: dp
       implicit none
       integer :: n,nprob
       real(dp) :: factor
       real(dp) :: x(n)
     end subroutine initpt

     subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag, &
          mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr, &
          qtf,wa1,wa2,wa3,wa4)
       use mo_kind, only: dp
       implicit none
       integer :: n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr
       real(dp) :: xtol,epsfcn,factor
       real(dp) :: x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr), &
            qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
       interface
          subroutine fcn(n,x,fvec,iflag)
            use mo_kind, only: dp
            implicit none
            integer,intent(in) :: n
            integer,intent(inout) :: iflag
            real(dp),intent(in) :: x(n)
            real(dp),intent(out) :: fvec(n)
          end subroutine fcn
       end interface
     end subroutine hybrd

     subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
       use mo_kind, only: dp
       implicit none
       integer :: n,info,lwa
       real(dp) :: tol
       real(dp) :: x(n),fvec(n),wa(lwa)
       interface
          subroutine fcn(n,x,fvec,iflag)
            use mo_kind, only: dp
            implicit none
            integer :: n,iflag
            real(dp) :: x(n),fvec(n)
          end subroutine fcn
       end interface
     end subroutine hybrd1

     subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,mode, &
          factor,nprint,info,nfev,njev,r,lr,qtf,wa1,wa2, &
          wa3,wa4)
       use mo_kind, only: dp
       implicit none
       integer :: n,ldfjac,maxfev,mode,nprint,info,nfev,njev,lr
       real(dp) :: xtol,factor
       real(dp) :: x(n),fvec(n),fjac(ldfjac,n),diag(n),r(lr), &
            qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
       interface
          subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
            use mo_kind, only: dp
            implicit none
            integer :: n,ldfjac,iflag
            real(dp) :: x(n),fvec(n),fjac(ldfjac,n)
          end subroutine fcn
       end interface
     end subroutine hybrj

     subroutine hybrj1(fcn,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa)
       use mo_kind, only: dp
       implicit none
       integer :: n,ldfjac,info,lwa
       real(dp) :: tol
       real(dp) :: x(n),fvec(n),fjac(ldfjac,n),wa(lwa)
       interface
          subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
            use mo_kind, only: dp
            implicit none
            integer :: n,ldfjac,iflag
            real(dp) :: x(n),fvec(n),fjac(ldfjac,n)
          end subroutine fcn
       end interface
     end subroutine hybrj1

     subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol, &
          maxfev,diag,mode,factor,nprint,info,nfev,njev, &
          ipvt,qtf,wa1,wa2,wa3,wa4)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,ldfjac,maxfev,mode,nprint,info,nfev,njev
       integer :: ipvt(n)
       real(dp) :: ftol,xtol,gtol,factor
       real(dp) :: x(n),fvec(m),fjac(ldfjac,n),diag(n),qtf(n), &
            wa1(n),wa2(n),wa3(n),wa4(m)
       interface
          subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
            use mo_kind, only: dp
            implicit none
            integer :: m,n,ldfjac,iflag
            real(dp) :: x(n),fvec(m),fjac(ldfjac,n)
          end subroutine fcn
       end interface
     end subroutine lmder

     subroutine lmder1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,lwa)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,ldfjac,info,lwa
       integer :: ipvt(n)
       real(dp) :: tol
       real(dp) :: x(n),fvec(m),fjac(ldfjac,n),wa(lwa)
       interface
          subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
            use mo_kind, only: dp
            implicit none
            integer :: m,n,ldfjac,iflag
            real(dp) :: x(n),fvec(m),fjac(ldfjac,n)
          end subroutine fcn
       end interface
     end subroutine lmder1

     subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn, &
          diag,mode,factor,nprint,info,nfev,fjac,ldfjac, &
          ipvt,qtf,wa1,wa2,wa3,wa4)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,maxfev,mode,nprint,info,nfev,ldfjac
       integer :: ipvt(n)
       real(dp) :: ftol,xtol,gtol,epsfcn,factor
       real(dp) :: x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n), &
            wa1(n),wa2(n),wa3(n),wa4(m)
       interface
          subroutine fcn(m,n,x,fvec,iflag)
            use mo_kind, only: dp
            implicit none
            integer :: m,n,iflag
            real(dp) :: x(n),fvec(m)
          end subroutine fcn
       end interface
     end subroutine lmdif

     subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,info,lwa
       integer :: iwa(n)
       real(dp) :: tol
       real(dp) :: x(n),fvec(m),wa(lwa)
       interface
          subroutine fcn(m,n,x,fvec,iflag)
            use mo_kind, only: dp
            implicit none
            integer :: m,n,iflag
            real(dp) :: x(n),fvec(m)
          end subroutine fcn
       end interface
     end subroutine lmdif1

     subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,wa2)
       use mo_kind, only: dp
       implicit none
       integer :: n,ldr
       integer :: ipvt(n)
       real(dp) :: delta,par
       real(dp) :: r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa1(n),wa2(n)
     end subroutine lmpar

     subroutine lmstr(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol, &
          maxfev,diag,mode,factor,nprint,info,nfev,njev, &
          ipvt,qtf,wa1,wa2,wa3,wa4)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,ldfjac,maxfev,mode,nprint,info,nfev,njev
       integer :: ipvt(n)
       real(dp) :: ftol,xtol,gtol,factor
       real(dp) :: x(n),fvec(m),fjac(ldfjac,n),diag(n),qtf(n), &
            wa1(n),wa2(n),wa3(n),wa4(m)
       interface
          subroutine fcn(m,n,x,fvec,fjrow,iflag)
            use mo_kind, only: dp
            implicit none
            integer :: m,n,iflag
            real(dp) :: x(n),fvec(m),fjrow(n)
          end subroutine fcn
       end interface
     end subroutine lmstr

     subroutine lmstr1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,lwa)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,ldfjac,info,lwa
       integer :: ipvt(n)
       real(dp) :: tol
       real(dp) :: x(n),fvec(m),fjac(ldfjac,n),wa(lwa)
       interface
          subroutine fcn(m,n,x,fvec,fjrow,iflag)
            use mo_kind, only: dp
            implicit none
            integer m,n,iflag
            real(dp) x(n),fvec(m),fjrow(n)
          end subroutine fcn
       end interface
     end subroutine lmstr1

     subroutine objfcn(n,x,f,nprob)
       use mo_kind, only: dp
       implicit none
       integer :: n,nprob
       real(dp) :: f
       real(dp) :: x(n)
     end subroutine objfcn

     subroutine qform(m,n,q,ldq,wa)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,ldq
       real(dp) :: q(ldq,m),wa(m)
     end subroutine qform

     subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,lda,lipvt
       integer :: ipvt(lipvt)
       logical :: pivot
       real(dp) :: a(lda,n),rdiag(n),acnorm(n),wa(n)
     end subroutine qrfac

     subroutine qrsolv2(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
       use mo_kind, only: dp
       implicit none
       integer :: n,ldr
       integer :: ipvt(n)
       real(dp) :: r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)
     end subroutine qrsolv2

     subroutine r1mpyq(m,n,a,lda,v,w)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,lda
       real(dp) :: a(lda,n),v(n),w(n)
     end subroutine r1mpyq

     subroutine r1updt(m,n,s,ls,u,v,w,sing)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,ls
       logical :: sing
       real(dp) :: s(ls),u(m),v(n),w(m)
     end subroutine r1updt

     subroutine rwupdt(n,r,ldr,w,b,alpha,cos,sin)
       use mo_kind, only: dp
       implicit none
       integer :: n,ldr
       real(dp) :: alpha
       real(dp) :: r(ldr,n),w(n),b(n),cos(n),sin(n)
     end subroutine rwupdt

     subroutine ssqfcn(m,n,x,fvec,nprob)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,nprob
       real(dp) :: x(n),fvec(m)
     end subroutine ssqfcn

     subroutine ssqjac(m,n,x,fjac,ldfjac,nprob)
       use mo_kind, only: dp
       implicit none
       integer :: m,n,ldfjac,nprob
       real(dp) :: x(n),fjac(ldfjac,n)
     end subroutine ssqjac

     subroutine vecfcn(n,x,fvec,nprob)
       use mo_kind, only: dp
       implicit none
       integer :: n,nprob
       real(dp) :: x(n),fvec(n)
     end subroutine vecfcn

     subroutine vecjac(n,x,fjac,ldfjac,nprob)
       use mo_kind, only: dp
       implicit none
       integer :: n,ldfjac,nprob
       real(dp) :: x(n),fjac(ldfjac,n)
     end subroutine vecjac

  end interface

contains

  !----------------------------------------------------------------------
  !
  ! add-ons
  !
  !----------------------------------------------------------------------

  subroutine qrinv(m,n,a,a1,lda,diag)

    ! compute inverse matrix in least-quare sense by the use
    ! of the QR factorisation

    ! implementation is reference, no effective, test for different
    ! m and n are sparse

    use mo_kind, only: dp
    implicit none

    integer :: m,n,lda
    real(dp) :: a(lda,n),a1(lda,n), diag(n)
    real(dp) :: r(m,n), q(m,n), qtb(m,n)
    integer :: i,ipvt(n)
    real(dp) :: rdiag(n),acnorm(n),wa(n),x(n),b(n),sdiag(n)
    character(len=10) :: fmt = '(xxf17.7)'

    if(dbg) then
       write(*,*) 'qrinv:'
       write(fmt(2:3),'(i2.2)') n
       write(*,fmt) a
    endif

    ! form the r matrix, r is upper trinagle (without diagonal)
    ! of the factorized a, diagonal is presented in rdiag

    call qrfac(m,n,a,lda,.true.,ipvt,n,rdiag,acnorm,wa)
    if(dbg) then
       write(*,*) 'qrfac:'
       write(*,fmt) a,rdiag,acnorm
    end if

    r = a
    do i = 1, n
       r(i,i) = rdiag(i)
    end do

    if(dbg) then
       write(*,*) 'r:'
       write(*,fmt) r
    endif

    ! form the q matrix

    call qform(m,n,a,lda,wa)

    if(dbg) then
       write(*,*) 'qform:'
       write(*,fmt) a,rdiag,acnorm
       write(*,*) 'ipvt:'
       write(*,*) ipvt
    end if

    q = a
    qtb = transpose(q)
    do i = 1, n
       b = 0
       b(i) = 1.0
       b = matmul(qtb,b)
       call qrsolv2(n,r,m,ipvt,diag,b,x,sdiag,wa)
       a1(i,:) = x
    enddo

  end subroutine qrinv

  !---------------------------------------------------------------------
  !
  ! Fortran 90 interfaces
  !
  !---------------------------------------------------------------------

  !    subroutine chkder_8(x,fvec,fjac,xp,fvecp,mode,err)
  !      integer, intent(in) :: mode
  !      real(dp) :: x(:),fvec(:),fjac(:,:),xp(:),fvecp(:),err(:)
  !      integer :: m,n,ldfjac
  !
  !      m = size(fvec)
  !      n = size(x)
  !      ldfjac = size(fjac,1)
  !
  !      call chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err)
  !    end subroutine chkder_8

end module mo_minpack

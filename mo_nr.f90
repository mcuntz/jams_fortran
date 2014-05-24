MODULE mo_nr

  ! This module provides the interfaces to the numerical recipes f90 routines.
  ! The routines itself are in individual files.
  !
  ! WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !    Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !    Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996

  ! USAGE
  !    All routines that are already brought to the FORTRAN_chs_lib structure are in the directory nr_chs.
  !    To use them, copy mo_nrutil, mo_rn and the routines from nr_chs to your directory and 'use mo_nr'
  !    in your code.
  !    Documentation is provided in the Fortran77 and Fortran90 books as pdf files in the directory nr_ori.

  ! EXTENTION
  !    The original numerical recipes routines are given in the directory nr_ori.
  !    There is also the documentation in form of the two Fortran books in PDF format and a Windows help file.
  !    
  !    If you want to use a routine that is not yet in nr_chs:
  !    - Copy routine from nr_ori to nr_chs.
  !    - Change nrtype to mo_kind.
  !    - Change nrutil to mo_nrutil.
  !    - Change nr to mo_nr.
  !    - Change I?B to I?.
  !    - Add mo_constants if constants such as PI are used.
  !    - Make sp and dp versions.
  !    - Look for the interface in mo_nr.f90 and add change it accordingly, for example add the dp version.
  !    - If you change input/output or similar of the routine the original documentation in the PDF files
  !      is not valid anymore.
  !      Please add the documentation structure from mo_template.f90 or consider making completely new module.

  ! Written Numerical Recipes, 1996
  ! Modified Matthias Cuntz, Nov 2011 - adapted to FORTRAN_chs_lib structure
  !                                   - added some double precision and 2d array routines

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! The UFZ Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The UFZ Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2011 Matthias Cuntz


  ! Note on Numerical Recipes License
  ! ---------------------------------
  ! Be aware that some code is under the Numerical Recipes License 3rd
  ! edition <http://www.nr.com/aboutNR3license.html>

  ! The Numerical Recipes Personal Single-User License lets you personally
  ! use Numerical Recipes code ("the code") on any number of computers,
  ! but only one computer at a time. You are not permitted to allow anyone
  ! else to access or use the code. You may, under this license, transfer
  ! precompiled, executable applications incorporating the code to other,
  ! unlicensed, persons, providing that (i) the application is
  ! noncommercial (i.e., does not involve the selling or licensing of the
  ! application for a fee), and (ii) the application was first developed,
  ! compiled, and successfully run by you, and (iii) the code is bound
  ! into the application in such a manner that it cannot be accessed as
  ! individual routines and cannot practicably be unbound and used in
  ! other programs. That is, under this license, your application user
  ! must not be able to use Numerical Recipes code as part of a program
  ! library or "mix and match" workbench.

  ! Businesses and organizations that purchase the disk or code download,
  ! and that thus acquire one or more Numerical Recipes Personal
  ! Single-User Licenses, may permanently assign those licenses, in the
  ! number acquired, to individual employees. Such an assignment must be
  ! made before the code is first used and, once made, it is irrevocable
  ! and can not be transferred. 

  ! If you do not hold a Numerical Recipes License, this code is only for
  ! informational and educational purposes but cannot be used.

  INTERFACE
     SUBROUTINE airy(x,ai,bi,aip,bip)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP), INTENT(OUT) :: ai,bi,aip,bip
     END SUBROUTINE airy
  END INTERFACE
  INTERFACE
     SUBROUTINE amebsa(p,y,pb,yb,ftol,func,iter,temptr)
       USE mo_kind
       INTEGER(I4), INTENT(INOUT) :: iter
       REAL(SP), INTENT(INOUT) :: yb
       REAL(SP), INTENT(IN) :: ftol,temptr
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: y,pb
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE amebsa
  END INTERFACE
  INTERFACE
     SUBROUTINE amoeba(p,y,ftol,func,iter)
       USE mo_kind
       INTEGER(I4), INTENT(OUT) :: iter
       REAL(SP), INTENT(IN) :: ftol
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE amoeba
  END INTERFACE
  INTERFACE
     SUBROUTINE anneal(x,y,iorder)
       USE mo_kind
       INTEGER(I4), DIMENSION(:), INTENT(INOUT) :: iorder
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
     END SUBROUTINE anneal
  END INTERFACE
  INTERFACE asolve
     SUBROUTINE asolve_sp(b,x,itrnsp)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: b
       REAL(SP), DIMENSION(:), INTENT(OUT) :: x
       INTEGER(I4), INTENT(IN) :: itrnsp
     END SUBROUTINE asolve_sp
     !MC
     SUBROUTINE asolve_dp(b,x,itrnsp)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: b
       REAL(DP), DIMENSION(:), INTENT(OUT) :: x
       INTEGER(I4), INTENT(IN) :: itrnsp
     END SUBROUTINE asolve_dp
  END INTERFACE
  INTERFACE atimes
     SUBROUTINE atimes_sp(x,r,itrnsp)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(:), INTENT(OUT) :: r
       INTEGER(I4), INTENT(IN) :: itrnsp
     END SUBROUTINE atimes_sp
     !MC
     SUBROUTINE atimes_dp(x,r,itrnsp)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(:), INTENT(OUT) :: r
       INTEGER(I4), INTENT(IN) :: itrnsp
     END SUBROUTINE atimes_dp
  END INTERFACE
  INTERFACE
     SUBROUTINE avevar(data,ave,var)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: data
       REAL(SP), INTENT(OUT) :: ave,var
     END SUBROUTINE avevar
  END INTERFACE
  INTERFACE
     SUBROUTINE balanc(a)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
     END SUBROUTINE balanc
  END INTERFACE
  INTERFACE
     SUBROUTINE banbks(a,m1,m2,al,indx,b)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: m1,m2
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: indx
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,al
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
     END SUBROUTINE banbks
  END INTERFACE
  INTERFACE
     SUBROUTINE bandec(a,m1,m2,al,indx,d)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: m1,m2
       INTEGER(I4), DIMENSION(:), INTENT(OUT) :: indx
       REAL(SP), INTENT(OUT) :: d
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: al
     END SUBROUTINE bandec
  END INTERFACE
  INTERFACE
     SUBROUTINE banmul(a,m1,m2,x,b)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: m1,m2
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(:), INTENT(OUT) :: b
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
     END SUBROUTINE banmul
  END INTERFACE
  INTERFACE
     SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)
       USE mo_kind
       REAL(SP), INTENT(IN) :: d1,d2
       REAL(SP), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
       REAL(SP), DIMENSION(4,4), INTENT(OUT) :: c
     END SUBROUTINE bcucof
  END INTERFACE
  INTERFACE
     SUBROUTINE bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,&
          ansy1,ansy2)
       USE mo_kind
       REAL(SP), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
       REAL(SP), INTENT(IN) :: x1l,x1u,x2l,x2u,x1,x2
       REAL(SP), INTENT(OUT) :: ansy,ansy1,ansy2
     END SUBROUTINE bcuint
  END INTERFACE
  INTERFACE beschb
     SUBROUTINE beschb_s(x,gam1,gam2,gampl,gammi)
       USE mo_kind
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(OUT) :: gam1,gam2,gampl,gammi
     END SUBROUTINE beschb_s
     !BL
     SUBROUTINE beschb_v(x,gam1,gam2,gampl,gammi)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(:), INTENT(OUT) :: gam1,gam2,gampl,gammi
     END SUBROUTINE beschb_v
  END INTERFACE
  INTERFACE bessi
     FUNCTION bessi_s(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: bessi_s
     END FUNCTION bessi_s
     !JM
     FUNCTION bessi_v(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: bessi_v
     END FUNCTION bessi_v
     FUNCTION dbessi_s(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: dbessi_s
     END FUNCTION dbessi_s
     !JM
     FUNCTION dbessi_v(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(size(x)) :: dbessi_v
     END FUNCTION dbessi_v
  END INTERFACE
  INTERFACE bessi0
     FUNCTION bessi0_s(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: bessi0_s
     END FUNCTION bessi0_s
     !JM
     FUNCTION bessi0_v(x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: bessi0_v
     END FUNCTION bessi0_v
     FUNCTION dbessi0_s(x)
       USE mo_kind
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: dbessi0_s
     END FUNCTION dbessi0_s
     !JM
     FUNCTION dbessi0_v(x)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(size(x)) :: dbessi0_v
     END FUNCTION dbessi0_v
  END INTERFACE
  INTERFACE bessi1
     FUNCTION bessi1_s(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: bessi1_s
     END FUNCTION bessi1_s
     !JM
     FUNCTION bessi1_v(x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: bessi1_v
     END FUNCTION bessi1_v
     FUNCTION dbessi1_s(x)
       USE mo_kind
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: dbessi1_s
     END FUNCTION dbessi1_s
     !JM
     FUNCTION dbessi1_v(x)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(size(x)) :: dbessi1_v
     END FUNCTION dbessi1_v
  END INTERFACE
  INTERFACE
     SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x,xnu
       REAL(SP), INTENT(OUT) :: ri,rk,rip,rkp
     END SUBROUTINE bessik
  END INTERFACE
  INTERFACE bessj
     FUNCTION bessj_s(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: bessj_s
     END FUNCTION bessj_s
     !JM
     FUNCTION bessj_v(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: bessj_v
     END FUNCTION bessj_v
  END INTERFACE
  INTERFACE bessj0
     FUNCTION bessj0_s(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: bessj0_s
     END FUNCTION bessj0_s
     !JM
     FUNCTION bessj0_v(x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: bessj0_v
     END FUNCTION bessj0_v
  END INTERFACE
  INTERFACE bessj1
     FUNCTION bessj1_s(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: bessj1_s
     END FUNCTION bessj1_s
     !JM
     FUNCTION bessj1_v(x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: bessj1_v
     END FUNCTION bessj1_v
  END INTERFACE
  INTERFACE bessjy
     SUBROUTINE bessjy_s(x,xnu,rj,ry,rjp,ryp)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x,xnu
       REAL(SP), INTENT(OUT) :: rj,ry,rjp,ryp
     END SUBROUTINE bessjy_s
     !JM
     SUBROUTINE bessjy_v(x,xnu,rj,ry,rjp,ryp)
       USE mo_kind
       REAL(SP), INTENT(IN) :: xnu
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
     END SUBROUTINE bessjy_v
  END INTERFACE
  INTERFACE bessk
     FUNCTION bessk_s(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: bessk_s
     END FUNCTION bessk_s
     !JM
     FUNCTION bessk_v(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: bessk_v
     END FUNCTION bessk_v
     FUNCTION dbessk_s(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: dbessk_s
     END FUNCTION dbessk_s
     !JM
     FUNCTION dbessk_v(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(size(x)) :: dbessk_v
     END FUNCTION dbessk_v
  END INTERFACE
  INTERFACE bessk0
     FUNCTION bessk0_s(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: bessk0_s
     END FUNCTION bessk0_s
     !JM
     FUNCTION bessk0_v(x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: bessk0_v
     END FUNCTION bessk0_v
     FUNCTION dbessk0_s(x)
       USE mo_kind
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: dbessk0_s
     END FUNCTION dbessk0_s
     !JM
     FUNCTION dbessk0_v(x)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(size(x)) :: dbessk0_v
     END FUNCTION dbessk0_v
  END INTERFACE
  INTERFACE bessk1
     FUNCTION bessk1_s(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: bessk1_s
     END FUNCTION bessk1_s
     !JM
     FUNCTION bessk1_v(x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: bessk1_v
     END FUNCTION bessk1_v
     FUNCTION dbessk1_s(x)
       USE mo_kind
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: dbessk1_s
     END FUNCTION dbessk1_s
     !JM
     FUNCTION dbessk1_v(x)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       REAL(DP), DIMENSION(size(x)) :: dbessk1_v
     END FUNCTION dbessk1_v
  END INTERFACE
  INTERFACE bessy
     FUNCTION bessy_s(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: bessy_s
     END FUNCTION bessy_s
     !BL
     FUNCTION bessy_v(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: bessy_v
     END FUNCTION bessy_v
  END INTERFACE
  INTERFACE bessy0
     FUNCTION bessy0_s(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: bessy0_s
     END FUNCTION bessy0_s
     !BL
     FUNCTION bessy0_v(x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: bessy0_v
     END FUNCTION bessy0_v
  END INTERFACE
  INTERFACE bessy1
     FUNCTION bessy1_s(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: bessy1_s
     END FUNCTION bessy1_s
     !BL
     FUNCTION bessy1_v(x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: bessy1_v
     END FUNCTION bessy1_v
  END INTERFACE
  INTERFACE beta
     FUNCTION beta_s(z,w)
       USE mo_kind
       REAL(SP), INTENT(IN) :: z,w
       REAL(SP) :: beta_s
     END FUNCTION beta_s
     !BL
     FUNCTION beta_v(z,w)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: z,w
       REAL(SP), DIMENSION(size(z)) :: beta_v
     END FUNCTION beta_v
  END INTERFACE
  INTERFACE betacf
     FUNCTION betacf_s(a,b,x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b,x
       REAL(SP) :: betacf_s
     END FUNCTION betacf_s
     !BL
     FUNCTION betacf_v(a,b,x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,x
       REAL(SP), DIMENSION(size(x)) :: betacf_v
     END FUNCTION betacf_v
  END INTERFACE
  INTERFACE betai
     FUNCTION betai_s(a,b,x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b,x
       REAL(SP) :: betai_s
     END FUNCTION betai_s
     !BL
     FUNCTION betai_v(a,b,x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,x
       REAL(SP), DIMENSION(size(a)) :: betai_v
     END FUNCTION betai_v
  END INTERFACE
  INTERFACE bico
     FUNCTION bico_s(n,k)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n,k
       REAL(SP) :: bico_s
     END FUNCTION bico_s
     !BL
     FUNCTION bico_v(n,k)
       USE mo_kind
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: n,k
       REAL(SP), DIMENSION(size(n)) :: bico_v
     END FUNCTION bico_v
  END INTERFACE
  INTERFACE
     FUNCTION bnldev(pp,n)
       USE mo_kind
       REAL(SP), INTENT(IN) :: pp
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP) :: bnldev
     END FUNCTION bnldev
  END INTERFACE
  INTERFACE
     FUNCTION brent(ax,bx,cx,func,tol,xmin)
       USE mo_kind
       REAL(SP), INTENT(IN) :: ax,bx,cx,tol
       REAL(SP), INTENT(OUT) :: xmin
       REAL(SP) :: brent
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION brent
  END INTERFACE
  INTERFACE
     SUBROUTINE broydn(x,check)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
       LOGICAL(LGT), INTENT(OUT) :: check
     END SUBROUTINE broydn
  END INTERFACE
  INTERFACE
     SUBROUTINE bsstep(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
       REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
       REAL(SP), INTENT(INOUT) :: x
       REAL(SP), INTENT(IN) :: htry,eps
       REAL(SP), INTENT(OUT) :: hdid,hnext
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
       END INTERFACE
     END SUBROUTINE bsstep
  END INTERFACE
  INTERFACE
     SUBROUTINE caldat(julian,mm,id,iyyy)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: julian
       INTEGER(I4), INTENT(OUT) :: mm,id,iyyy
     END SUBROUTINE caldat
  END INTERFACE
  INTERFACE
     FUNCTION chder(a,b,c)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), DIMENSION(:), INTENT(IN) :: c
       REAL(SP), DIMENSION(size(c)) :: chder
     END FUNCTION chder
  END INTERFACE
  INTERFACE chebev
     FUNCTION chebev_s(a,b,c,x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b,x
       REAL(SP), DIMENSION(:), INTENT(IN) :: c
       REAL(SP) :: chebev_s
     END FUNCTION chebev_s
     !BL
     FUNCTION chebev_v(a,b,c,x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), DIMENSION(:), INTENT(IN) :: c,x
       REAL(SP), DIMENSION(size(x)) :: chebev_v
     END FUNCTION chebev_v
  END INTERFACE
  INTERFACE
     FUNCTION chebft(a,b,n,func)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), DIMENSION(n) :: chebft
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION chebft
  END INTERFACE
  INTERFACE
     FUNCTION chebpc(c)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: c
       REAL(SP), DIMENSION(size(c)) :: chebpc
     END FUNCTION chebpc
  END INTERFACE
  INTERFACE
     FUNCTION chint(a,b,c)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), DIMENSION(:), INTENT(IN) :: c
       REAL(SP), DIMENSION(size(c)) :: chint
     END FUNCTION chint
  END INTERFACE
  INTERFACE choldc
     SUBROUTINE choldc_sp(a,p)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
       REAL(SP), DIMENSION(:), INTENT(OUT) :: p
     END SUBROUTINE choldc_sp
     SUBROUTINE choldc_dp(a,p)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
       REAL(DP), DIMENSION(:), INTENT(OUT) :: p
     END SUBROUTINE choldc_dp
  END INTERFACE
  INTERFACE
     SUBROUTINE cholsl(a,p,b,x)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
       REAL(SP), DIMENSION(:), INTENT(IN) :: p,b
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
     END SUBROUTINE cholsl
  END INTERFACE
  INTERFACE
     SUBROUTINE chsone(bins,ebins,knstrn,df,chsq,prob)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: knstrn
       REAL(SP), INTENT(OUT) :: df,chsq,prob
       REAL(SP), DIMENSION(:), INTENT(IN) :: bins,ebins
     END SUBROUTINE chsone
  END INTERFACE
  INTERFACE
     SUBROUTINE chstwo(bins1,bins2,knstrn,df,chsq,prob)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: knstrn
       REAL(SP), INTENT(OUT) :: df,chsq,prob
       REAL(SP), DIMENSION(:), INTENT(IN) :: bins1,bins2
     END SUBROUTINE chstwo
  END INTERFACE
  INTERFACE
     SUBROUTINE cisi(x,ci,si)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP), INTENT(OUT) :: ci,si
     END SUBROUTINE cisi
  END INTERFACE
  INTERFACE
     SUBROUTINE cntab1(nn,chisq,df,prob,cramrv,ccc)
       USE mo_kind
       INTEGER(I4), DIMENSION(:,:), INTENT(IN) :: nn
       REAL(SP), INTENT(OUT) :: chisq,df,prob,cramrv,ccc
     END SUBROUTINE cntab1
  END INTERFACE
  INTERFACE
     SUBROUTINE cntab2(nn,h,hx,hy,hygx,hxgy,uygx,uxgy,uxy)
       USE mo_kind
       INTEGER(I4), DIMENSION(:,:), INTENT(IN) :: nn
       REAL(SP), INTENT(OUT) :: h,hx,hy,hygx,hxgy,uygx,uxgy,uxy
     END SUBROUTINE cntab2
  END INTERFACE
  INTERFACE
     FUNCTION convlv(data,respns,isign)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: data
       REAL(SP), DIMENSION(:), INTENT(IN) :: respns
       INTEGER(I4), INTENT(IN) :: isign
       REAL(SP), DIMENSION(size(data)) :: convlv
     END FUNCTION convlv
  END INTERFACE
  INTERFACE
     FUNCTION correl(data1,data2)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
       REAL(SP), DIMENSION(size(data1)) :: correl
     END FUNCTION correl
  END INTERFACE
  INTERFACE
     SUBROUTINE cosft1(y)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
     END SUBROUTINE cosft1
  END INTERFACE
  INTERFACE
     SUBROUTINE cosft2(y,isign)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE cosft2
  END INTERFACE
  INTERFACE covsrt
     SUBROUTINE covsrt_sp(covar,maska)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
       LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
     END SUBROUTINE covsrt_sp
     SUBROUTINE covsrt_dp(covar,maska)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: covar
       LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
     END SUBROUTINE covsrt_dp
  END INTERFACE covsrt
  INTERFACE
     SUBROUTINE cyclic(a,b,c,alpha,beta,r,x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN):: a,b,c,r
       REAL(SP), INTENT(IN) :: alpha,beta
       REAL(SP), DIMENSION(:), INTENT(OUT):: x
     END SUBROUTINE cyclic
  END INTERFACE
  INTERFACE
     SUBROUTINE daub4(a,isign)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE daub4
  END INTERFACE
  INTERFACE dawson
     FUNCTION dawson_s(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: dawson_s
     END FUNCTION dawson_s
     !BL
     FUNCTION dawson_v(x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: dawson_v
     END FUNCTION dawson_v
  END INTERFACE
  INTERFACE
     FUNCTION dbrent(ax,bx,cx,func,dbrent_dfunc,tol,xmin)
       USE mo_kind
       REAL(SP), INTENT(IN) :: ax,bx,cx,tol
       REAL(SP), INTENT(OUT) :: xmin
       REAL(SP) :: dbrent
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
          !BL
          FUNCTION dbrent_dfunc(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: dbrent_dfunc
          END FUNCTION dbrent_dfunc
       END INTERFACE
     END FUNCTION dbrent
  END INTERFACE
  INTERFACE
     SUBROUTINE ddpoly(c,x,pd)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP), DIMENSION(:), INTENT(IN) :: c
       REAL(SP), DIMENSION(:), INTENT(OUT) :: pd
     END SUBROUTINE ddpoly
  END INTERFACE
  INTERFACE
     FUNCTION decchk(string,ch)
       USE mo_kind
       CHARACTER(1), DIMENSION(:), INTENT(IN) :: string
       CHARACTER(1), INTENT(OUT) :: ch
       LOGICAL(LGT) :: decchk
     END FUNCTION decchk
  END INTERFACE
  INTERFACE dfpmin
     SUBROUTINE dfpmin_sp(p,gtol,iter,fret,func,dfunc)
       USE mo_kind
       INTEGER(I4), INTENT(OUT) :: iter
       REAL(SP), INTENT(IN) :: gtol
       REAL(SP), INTENT(OUT) :: fret
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
       INTERFACE
          FUNCTION func(p)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: p
            REAL(SP) :: func
          END FUNCTION func
          !BL
          FUNCTION dfunc(p)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: p
            REAL(SP), DIMENSION(size(p)) :: dfunc
          END FUNCTION dfunc
       END INTERFACE
     END SUBROUTINE dfpmin_sp
     SUBROUTINE dfpmin_dp(p,gtol,iter,fret,func,dfunc)
       USE mo_kind
       INTEGER(I4), INTENT(OUT) :: iter
       REAL(DP), INTENT(IN) :: gtol
       REAL(DP), INTENT(OUT) :: fret
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
       INTERFACE
          FUNCTION func(p)
            USE mo_kind
            REAL(DP), DIMENSION(:), INTENT(IN) :: p
            REAL(DP) :: func
          END FUNCTION func
          !BL
          FUNCTION dfunc(p)
            USE mo_kind
            REAL(DP), DIMENSION(:), INTENT(IN) :: p
            REAL(DP), DIMENSION(size(p)) :: dfunc
          END FUNCTION dfunc
       END INTERFACE
     END SUBROUTINE dfpmin_dp
  END INTERFACE dfpmin
  INTERFACE
     FUNCTION dfridr(func,x,h,err)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x,h
       REAL(SP), INTENT(OUT) :: err
       REAL(SP) :: dfridr
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION dfridr
  END INTERFACE
  INTERFACE
     SUBROUTINE dftcor(w,delta,a,b,endpts,corre,corim,corfac)
       USE mo_kind
       REAL(SP), INTENT(IN) :: w,delta,a,b
       REAL(SP), INTENT(OUT) :: corre,corim,corfac
       REAL(SP), DIMENSION(:), INTENT(IN) :: endpts
     END SUBROUTINE dftcor
  END INTERFACE
  INTERFACE
     SUBROUTINE dftint(func,a,b,w,cosint,sinint)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b,w
       REAL(SP), INTENT(OUT) :: cosint,sinint
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE dftint
  END INTERFACE
  INTERFACE
     SUBROUTINE difeq(k,k1,k2,jsf,is1,isf,indexv,s,y)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: is1,isf,jsf,k,k1,k2
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: indexv
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: s
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: y
     END SUBROUTINE difeq
  END INTERFACE
  INTERFACE
     FUNCTION eclass(lista,listb,n)
       USE mo_kind
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: lista,listb
       INTEGER(I4), INTENT(IN) :: n
       INTEGER(I4), DIMENSION(n) :: eclass
     END FUNCTION eclass
  END INTERFACE
  INTERFACE
     FUNCTION eclazz(equiv,n)
       USE mo_kind
       INTERFACE
          FUNCTION equiv(i,j)
            USE mo_kind
            LOGICAL(LGT) :: equiv
            INTEGER(I4), INTENT(IN) :: i,j
          END FUNCTION equiv
       END INTERFACE
       INTEGER(I4), INTENT(IN) :: n
       INTEGER(I4), DIMENSION(n) :: eclazz
     END FUNCTION eclazz
  END INTERFACE
  INTERFACE ei
     FUNCTION ei_sp(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: ei_sp
     END FUNCTION ei_sp
     FUNCTION ei_dp(x)
       USE mo_kind
       REAL(dp), INTENT(IN) :: x
       REAL(dp) :: ei_dp
     END FUNCTION ei_dp
  END INTERFACE
  INTERFACE
     SUBROUTINE eigsrt(d,v)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: v
     END SUBROUTINE eigsrt
  END INTERFACE
  INTERFACE elle
     FUNCTION elle_s(phi,ak)
       USE mo_kind
       REAL(SP), INTENT(IN) :: phi,ak
       REAL(SP) :: elle_s
     END FUNCTION elle_s
     !BL
     FUNCTION elle_v(phi,ak)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: phi,ak
       REAL(SP), DIMENSION(size(phi)) :: elle_v
     END FUNCTION elle_v
  END INTERFACE
  INTERFACE ellf
     FUNCTION ellf_s(phi,ak)
       USE mo_kind
       REAL(SP), INTENT(IN) :: phi,ak
       REAL(SP) :: ellf_s
     END FUNCTION ellf_s
     !BL
     FUNCTION ellf_v(phi,ak)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: phi,ak
       REAL(SP), DIMENSION(size(phi)) :: ellf_v
     END FUNCTION ellf_v
  END INTERFACE
  INTERFACE ellpi
     FUNCTION ellpi_s(phi,en,ak)
       USE mo_kind
       REAL(SP), INTENT(IN) :: phi,en,ak
       REAL(SP) :: ellpi_s
     END FUNCTION ellpi_s
     !BL
     FUNCTION ellpi_v(phi,en,ak)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: phi,en,ak
       REAL(SP), DIMENSION(size(phi)) :: ellpi_v
     END FUNCTION ellpi_v
  END INTERFACE
  INTERFACE
     SUBROUTINE elmhes(a)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
     END SUBROUTINE elmhes
  END INTERFACE
  INTERFACE erf
     FUNCTION erf_s(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: erf_s
     END FUNCTION erf_s
     !BL
     FUNCTION erf_v(x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: erf_v
     END FUNCTION erf_v
  END INTERFACE
  INTERFACE erfc
     FUNCTION erfc_s(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: erfc_s
     END FUNCTION erfc_s
     !BL
     FUNCTION erfc_v(x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: erfc_v
     END FUNCTION erfc_v
  END INTERFACE
  INTERFACE erfcc
     FUNCTION erfcc_s(x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: erfcc_s
     END FUNCTION erfcc_s
     !BL
     FUNCTION erfcc_v(x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: erfcc_v
     END FUNCTION erfcc_v
  END INTERFACE
  INTERFACE
     SUBROUTINE eulsum(sum,term,jterm)
       USE mo_kind
       REAL(SP), INTENT(INOUT) :: sum
       REAL(SP), INTENT(IN) :: term
       INTEGER(I4), INTENT(IN) :: jterm
     END SUBROUTINE eulsum
  END INTERFACE
  INTERFACE
     FUNCTION evlmem(fdt,d,xms)
       USE mo_kind
       REAL(SP), INTENT(IN) :: fdt,xms
       REAL(SP), DIMENSION(:), INTENT(IN) :: d
       REAL(SP) :: evlmem
     END FUNCTION evlmem
  END INTERFACE
  INTERFACE expdev
     SUBROUTINE expdev_s(harvest)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: harvest
     END SUBROUTINE expdev_s
     !BL
     SUBROUTINE expdev_v(harvest)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
     END SUBROUTINE expdev_v
  END INTERFACE
  INTERFACE expint
     FUNCTION expint_sp(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: expint_sp
     END FUNCTION expint_sp
     !MC
     FUNCTION expint_dp(n,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: expint_dp
     END FUNCTION expint_dp
  END INTERFACE
  INTERFACE factln
     FUNCTION factln_s(n)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP) :: factln_s
     END FUNCTION factln_s
     !BL
     FUNCTION factln_v(n)
       USE mo_kind
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: n
       REAL(SP), DIMENSION(size(n)) :: factln_v
     END FUNCTION factln_v
  END INTERFACE
  INTERFACE factrl
     FUNCTION factrl_s(n)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP) :: factrl_s
     END FUNCTION factrl_s
     !BL
     FUNCTION factrl_v(n)
       USE mo_kind
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: n
       REAL(SP), DIMENSION(size(n)) :: factrl_v
     END FUNCTION factrl_v
  END INTERFACE
  INTERFACE
     SUBROUTINE fasper(x,y,ofac,hifac,px,py,jmax,prob)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
       REAL(SP), INTENT(IN) :: ofac,hifac
       INTEGER(I4), INTENT(OUT) :: jmax
       REAL(SP), INTENT(OUT) :: prob
       REAL(SP), DIMENSION(:), POINTER :: px,py
     END SUBROUTINE fasper
  END INTERFACE
  INTERFACE
     SUBROUTINE fdjac(x,fvec,df)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: fvec
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: df
     END SUBROUTINE fdjac
  END INTERFACE
  INTERFACE
     SUBROUTINE fgauss(x,a,y,dyda)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
       REAL(SP), DIMENSION(:), INTENT(OUT) :: y
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
     END SUBROUTINE fgauss
  END INTERFACE
  INTERFACE
     SUBROUTINE fit(x,y,a,b,siga,sigb,chi2,q,sig)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
       REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
       REAL(SP), DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
     END SUBROUTINE fit
  END INTERFACE
  INTERFACE
     SUBROUTINE fitexy(x,y,sigx,sigy,a,b,siga,sigb,chi2,q)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sigx,sigy
       REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
     END SUBROUTINE fitexy
  END INTERFACE
  INTERFACE
     SUBROUTINE fixrts(d)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
     END SUBROUTINE fixrts
  END INTERFACE
  INTERFACE
     FUNCTION fleg(x,n)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), DIMENSION(n) :: fleg
     END FUNCTION fleg
  END INTERFACE
  INTERFACE
     SUBROUTINE flmoon(n,nph,jd,frac)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n,nph
       INTEGER(I4), INTENT(OUT) :: jd
       REAL(SP), INTENT(OUT) :: frac
     END SUBROUTINE flmoon
  END INTERFACE
  INTERFACE four1
     SUBROUTINE four1_dp(data,isign)
       USE mo_kind
       COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE four1_dp
     !BL
     SUBROUTINE four1_sp(data,isign)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE four1_sp
  END INTERFACE
  INTERFACE
     SUBROUTINE four1_alt(data,isign)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE four1_alt
  END INTERFACE
  INTERFACE
     SUBROUTINE four1_gather(data,isign)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE four1_gather
  END INTERFACE
  INTERFACE
     SUBROUTINE four2(data,isign)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
       INTEGER(I4),INTENT(IN) :: isign
     END SUBROUTINE four2
  END INTERFACE
  INTERFACE
     SUBROUTINE four2_alt(data,isign)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE four2_alt
  END INTERFACE
  INTERFACE
     SUBROUTINE four3(data,isign)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
       INTEGER(I4),INTENT(IN) :: isign
     END SUBROUTINE four3
  END INTERFACE
  INTERFACE
     SUBROUTINE four3_alt(data,isign)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE four3_alt
  END INTERFACE
  INTERFACE
     SUBROUTINE fourcol(data,isign)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE fourcol
  END INTERFACE
  INTERFACE
     SUBROUTINE fourcol_3d(data,isign)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE fourcol_3d
  END INTERFACE
  INTERFACE
     SUBROUTINE fourn_gather(data,nn,isign)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: nn
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE fourn_gather
  END INTERFACE
  INTERFACE fourrow
     SUBROUTINE fourrow_dp(data,isign)
       USE mo_kind
       COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE fourrow_dp
     !BL
     SUBROUTINE fourrow_sp(data,isign)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE fourrow_sp
  END INTERFACE
  INTERFACE
     SUBROUTINE fourrow_3d(data,isign)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE fourrow_3d
  END INTERFACE
  INTERFACE
     FUNCTION fpoly(x,n)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), DIMENSION(n) :: fpoly
     END FUNCTION fpoly
  END INTERFACE
  INTERFACE
     SUBROUTINE fred2(a,b,t,f,w,g,ak)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), DIMENSION(:), INTENT(OUT) :: t,f,w
       INTERFACE
          FUNCTION g(t)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: t
            REAL(SP), DIMENSION(size(t)) :: g
          END FUNCTION g
          !BL
          FUNCTION ak(t,s)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
            REAL(SP), DIMENSION(size(t),size(s)) :: ak
          END FUNCTION ak
       END INTERFACE
     END SUBROUTINE fred2
  END INTERFACE
  INTERFACE
     FUNCTION fredin(x,a,b,t,f,w,g,ak)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,t,f,w
       REAL(SP), DIMENSION(size(x)) :: fredin
       INTERFACE
          FUNCTION g(t)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: t
            REAL(SP), DIMENSION(size(t)) :: g
          END FUNCTION g
          !BL
          FUNCTION ak(t,s)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
            REAL(SP), DIMENSION(size(t),size(s)) :: ak
          END FUNCTION ak
       END INTERFACE
     END FUNCTION fredin
  END INTERFACE
  INTERFACE
     SUBROUTINE frenel(x,s,c)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP), INTENT(OUT) :: s,c
     END SUBROUTINE frenel
  END INTERFACE
  INTERFACE
     SUBROUTINE frprmn(p,ftol,iter,fret)
       USE mo_kind
       INTEGER(I4), INTENT(OUT) :: iter
       REAL(SP), INTENT(IN) :: ftol
       REAL(SP), INTENT(OUT) :: fret
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
     END SUBROUTINE frprmn
  END INTERFACE
  INTERFACE
     SUBROUTINE ftest(data1,data2,f,prob)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: f,prob
       REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
     END SUBROUTINE ftest
  END INTERFACE
  INTERFACE
     FUNCTION gamdev(ia)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: ia
       REAL(SP) :: gamdev
     END FUNCTION gamdev
  END INTERFACE
  INTERFACE gammln
     FUNCTION gammln_s(xx)
       USE mo_kind
       REAL(SP), INTENT(IN) :: xx
       REAL(SP) :: gammln_s
     END FUNCTION gammln_s
     !BL
     FUNCTION gammln_v(xx)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: xx
       REAL(SP), DIMENSION(size(xx)) :: gammln_v
     END FUNCTION gammln_v
     !MC
     FUNCTION dgammln_s(z)
       USE mo_kind
       REAL(DP), INTENT(IN) :: z
       REAL(DP) :: dgammln_s
     END FUNCTION dgammln_s
     !MC
     FUNCTION dgammln_v(z)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: z
       REAL(DP), DIMENSION(size(z)) :: dgammln_v
     END FUNCTION dgammln_v
  END INTERFACE
  INTERFACE gamma
     FUNCTION gamma_s(xx)
       USE mo_kind
       REAL(SP), INTENT(IN) :: xx
       REAL(SP) :: gamma_s
     END FUNCTION gamma_s
     !MC
     FUNCTION gamma_v(xx)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: xx
       REAL(SP), DIMENSION(size(xx)) :: gamma_v
     END FUNCTION gamma_v
     !MC
     FUNCTION dgamma_s(z)
       USE mo_kind
       REAL(DP), INTENT(IN) :: z
       REAL(DP) :: dgamma_s
     END FUNCTION dgamma_s
     !MC
     FUNCTION dgamma_v(z)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: z
       REAL(DP), DIMENSION(size(z)) :: dgamma_v
     END FUNCTION dgamma_v
  END INTERFACE
  INTERFACE gammp
     FUNCTION gammp_s(a,x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,x
       REAL(SP) :: gammp_s
     END FUNCTION gammp_s
     !BL
     FUNCTION gammp_v(a,x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
       REAL(SP), DIMENSION(size(a)) :: gammp_v
     END FUNCTION gammp_v
  END INTERFACE
  INTERFACE gammq
     FUNCTION gammq_s(a,x,gln)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,x
       REAL(SP), OPTIONAL, INTENT(OUT) :: gln
       REAL(SP) :: gammq_s
     END FUNCTION gammq_s
     !MC
     FUNCTION gammq_v(a,x,gln)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
       REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
       REAL(SP), DIMENSION(size(a)) :: gammq_v
     END FUNCTION gammq_v
     !MC
     FUNCTION dgammq_s(a,x,gln)
       USE mo_kind
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: a,x
       REAL(DP), OPTIONAL, INTENT(OUT) :: gln
       REAL(DP) :: dgammq_s
     END FUNCTION dgammq_s
     !MC
     FUNCTION dgammq_v(a,x,gln)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
       REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
       REAL(DP), DIMENSION(size(a)) :: dgammq_v
     END FUNCTION dgammq_v
  END INTERFACE
  INTERFACE igamma
     FUNCTION igamma_s(a,x)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,x
       REAL(SP) :: igamma_s
     END FUNCTION igamma_s
     !MC
     FUNCTION igamma_v(a,x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
       REAL(SP), DIMENSION(size(a)) :: igamma_v
     END FUNCTION igamma_v
     !MC
     FUNCTION digamma_s(a,x)
       USE mo_kind
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: a,x
       REAL(DP) :: digamma_s
     END FUNCTION digamma_s
     !MC
     FUNCTION digamma_v(a,x)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
       REAL(DP), DIMENSION(size(a)) :: digamma_v
     END FUNCTION digamma_v
  END INTERFACE
  INTERFACE gasdev
     SUBROUTINE gasdev_s(harvest)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: harvest
     END SUBROUTINE gasdev_s
     !BL
     SUBROUTINE gasdev_v(harvest)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
     END SUBROUTINE gasdev_v
  END INTERFACE
  INTERFACE
     SUBROUTINE gaucof(a,b,amu0,x,w)
       USE mo_kind
       REAL(SP), INTENT(IN) :: amu0
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
       REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
     END SUBROUTINE gaucof
  END INTERFACE
  INTERFACE
     SUBROUTINE gauher(x,w)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
     END SUBROUTINE gauher
  END INTERFACE
  INTERFACE
     SUBROUTINE gaujac(x,w,alf,bet)
       USE mo_kind
       REAL(SP), INTENT(IN) :: alf,bet
       REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
     END SUBROUTINE gaujac
  END INTERFACE
  INTERFACE
     SUBROUTINE gaulag(x,w,alf)
       USE mo_kind
       REAL(SP), INTENT(IN) :: alf
       REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
     END SUBROUTINE gaulag
  END INTERFACE
  INTERFACE
     SUBROUTINE gauleg(x1,x2,x,w)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x1,x2
       REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
     END SUBROUTINE gauleg
  END INTERFACE
  INTERFACE gaussj
     SUBROUTINE gaussj_sp(a,b)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
     END SUBROUTINE gaussj_sp
     SUBROUTINE gaussj_dp(a,b)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b
     END SUBROUTINE gaussj_dp
  END INTERFACE
  INTERFACE gcf
     FUNCTION gcf_s(a,x,gln)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,x
       REAL(SP), OPTIONAL, INTENT(OUT) :: gln
       REAL(SP) :: gcf_s
     END FUNCTION gcf_s
     !BL
     FUNCTION gcf_v(a,x,gln)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
       REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
       REAL(SP), DIMENSION(size(a)) :: gcf_v
     END FUNCTION gcf_v
     !MC
     FUNCTION dgcf_s(a,x,gln)
       USE mo_kind
       REAL(DP), INTENT(IN) :: a,x
       REAL(DP), OPTIONAL, INTENT(OUT) :: gln
       REAL(DP) :: dgcf_s
     END FUNCTION dgcf_s
     !MC
     FUNCTION dgcf_v(a,x,gln)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
       REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
       REAL(DP), DIMENSION(size(a)) :: dgcf_v
     END FUNCTION dgcf_v
  END INTERFACE
  INTERFACE golden
     FUNCTION golden_sp(ax,bx,cx,func,tol,xmin)
       USE mo_kind
       REAL(SP), INTENT(IN) :: ax,bx,cx,tol
       REAL(SP), INTENT(OUT) :: xmin
       REAL(SP) :: golden_sp
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION golden_sp
     FUNCTION golden_dp(ax,bx,cx,func,tol,xmin)
       USE mo_kind
       REAL(DP), INTENT(IN) :: ax,bx,cx,tol
       REAL(DP), INTENT(OUT) :: xmin
       REAL(DP) :: golden_dp
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(DP), INTENT(IN) :: x
            REAL(DP) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION golden_dp
  END INTERFACE golden
  INTERFACE gser
     FUNCTION gser_s(a,x,gln)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,x
       REAL(SP), OPTIONAL, INTENT(OUT) :: gln
       REAL(SP) :: gser_s
     END FUNCTION gser_s
     !BL
     FUNCTION gser_v(a,x,gln)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
       REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
       REAL(SP), DIMENSION(size(a)) :: gser_v
     END FUNCTION gser_v
     !MC
     FUNCTION dgser_s(a,x,gln)
       USE mo_kind
       REAL(DP), INTENT(IN) :: a,x
       REAL(DP), OPTIONAL, INTENT(OUT) :: gln
       REAL(DP) :: dgser_s
     END FUNCTION dgser_s
     !MC
     FUNCTION dgser_v(a,x,gln)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
       REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
       REAL(DP), DIMENSION(size(a)) :: dgser_v
     END FUNCTION dgser_v
  END INTERFACE
  INTERFACE
     SUBROUTINE hqr(a,wr,wi)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(OUT) :: wr,wi
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
     END SUBROUTINE hqr
  END INTERFACE
  INTERFACE
     SUBROUTINE hunt(xx,x,jlo)
       USE mo_kind
       INTEGER(I4), INTENT(INOUT) :: jlo
       REAL(SP), INTENT(IN) :: x
       REAL(SP), DIMENSION(:), INTENT(IN) :: xx
     END SUBROUTINE hunt
  END INTERFACE
  INTERFACE
     SUBROUTINE hypdrv(s,ry,rdyds)
       USE mo_kind
       REAL(SP), INTENT(IN) :: s
       REAL(SP), DIMENSION(:), INTENT(IN) :: ry
       REAL(SP), DIMENSION(:), INTENT(OUT) :: rdyds
     END SUBROUTINE hypdrv
  END INTERFACE
  INTERFACE
     FUNCTION hypgeo(a,b,c,z)
       USE mo_kind
       COMPLEX(SPC), INTENT(IN) :: a,b,c,z
       COMPLEX(SPC) :: hypgeo
     END FUNCTION hypgeo
  END INTERFACE
  INTERFACE
     SUBROUTINE hypser(a,b,c,z,series,deriv)
       USE mo_kind
       COMPLEX(SPC), INTENT(IN) :: a,b,c,z
       COMPLEX(SPC), INTENT(OUT) :: series,deriv
     END SUBROUTINE hypser
  END INTERFACE
  INTERFACE
     FUNCTION icrc(crc,buf,jinit,jrev)
       USE mo_kind
       CHARACTER(1), DIMENSION(:), INTENT(IN) :: buf
       INTEGER(I2), INTENT(IN) :: crc,jinit
       INTEGER(I4), INTENT(IN) :: jrev
       INTEGER(I2) :: icrc
     END FUNCTION icrc
  END INTERFACE
  INTERFACE
     FUNCTION igray(n,is)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n,is
       INTEGER(I4) :: igray
     END FUNCTION igray
  END INTERFACE
  INTERFACE
     RECURSIVE SUBROUTINE index_bypack(arr,index,partial)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: arr
       INTEGER(I4), DIMENSION(:), INTENT(INOUT) :: index
       INTEGER, OPTIONAL, INTENT(IN) :: partial
     END SUBROUTINE index_bypack
  END INTERFACE
  INTERFACE indexx
     SUBROUTINE indexx_sp(arr,index)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: arr
       INTEGER(I4), DIMENSION(:), INTENT(OUT) :: index
     END SUBROUTINE indexx_sp
     SUBROUTINE indexx_i4(iarr,index)
       USE mo_kind
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: iarr
       INTEGER(I4), DIMENSION(:), INTENT(OUT) :: index
     END SUBROUTINE indexx_i4
  END INTERFACE
  INTERFACE
     FUNCTION interp(uc)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: uc
       REAL(DP), DIMENSION(2*size(uc,1)-1,2*size(uc,1)-1) :: interp
     END FUNCTION interp
  END INTERFACE
  INTERFACE
     ! original from NR --> "rank" überdeckt eingebaute Funktion mit gleichem Namen
     ! FUNCTION rank(indx)
     ! renamed
     FUNCTION rank_nr(indx)
       USE mo_kind
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: indx
       INTEGER(I4), DIMENSION(size(indx)) :: rank_nr
     END FUNCTION rank_nr
  END INTERFACE
  INTERFACE
     FUNCTION irbit1(iseed)
       USE mo_kind
       INTEGER(I4), INTENT(INOUT) :: iseed
       INTEGER(I4) :: irbit1
     END FUNCTION irbit1
  END INTERFACE
  INTERFACE
     FUNCTION irbit2(iseed)
       USE mo_kind
       INTEGER(I4), INTENT(INOUT) :: iseed
       INTEGER(I4) :: irbit2
     END FUNCTION irbit2
  END INTERFACE
  INTERFACE
     SUBROUTINE jacobi(a,d,v,nrot)
       USE mo_kind
       INTEGER(I4), INTENT(OUT) :: nrot
       REAL(SP), DIMENSION(:), INTENT(OUT) :: d
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
     END SUBROUTINE jacobi
  END INTERFACE
  INTERFACE
     SUBROUTINE jacobn(x,y,dfdx,dfdy)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x
       REAL(SP), DIMENSION(:), INTENT(IN) :: y
       REAL(SP), DIMENSION(:), INTENT(OUT) :: dfdx
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dfdy
     END SUBROUTINE jacobn
  END INTERFACE
  INTERFACE
     FUNCTION julday(mm,id,iyyy)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: mm,id,iyyy
       INTEGER(I4) :: julday
     END FUNCTION julday
  END INTERFACE
  INTERFACE
     SUBROUTINE kendl1(data1,data2,tau,z,prob)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: tau,z,prob
       REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
     END SUBROUTINE kendl1
  END INTERFACE
  INTERFACE
     SUBROUTINE kendl2(tab,tau,z,prob)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: tab
       REAL(SP), INTENT(OUT) :: tau,z,prob
     END SUBROUTINE kendl2
  END INTERFACE
  INTERFACE
     FUNCTION kermom(y,m)
       USE mo_kind
       REAL(DP), INTENT(IN) :: y
       INTEGER(I4), INTENT(IN) :: m
       REAL(DP), DIMENSION(m) :: kermom
     END FUNCTION kermom
  END INTERFACE
  INTERFACE
     SUBROUTINE ks2d1s(x1,y1,quadvl,d1,prob)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x1,y1
       REAL(SP), INTENT(OUT) :: d1,prob
       INTERFACE
          SUBROUTINE quadvl(x,y,fa,fb,fc,fd)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x,y
            REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
          END SUBROUTINE quadvl
       END INTERFACE
     END SUBROUTINE ks2d1s
  END INTERFACE
  INTERFACE
     SUBROUTINE ks2d2s(x1,y1,x2,y2,d,prob)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x1,y1,x2,y2
       REAL(SP), INTENT(OUT) :: d,prob
     END SUBROUTINE ks2d2s
  END INTERFACE
  INTERFACE
     SUBROUTINE ksone(data,func,d,prob)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: d,prob
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE ksone
  END INTERFACE
  INTERFACE
     SUBROUTINE kstwo(data1,data2,d,prob)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: d,prob
       REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
     END SUBROUTINE kstwo
  END INTERFACE
  INTERFACE laguer
     SUBROUTINE laguer_sp(a,x,its)
       USE mo_kind
       INTEGER(I4), INTENT(OUT) :: its
       COMPLEX(SPC), INTENT(INOUT) :: x
       COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
     END SUBROUTINE laguer_sp
     SUBROUTINE laguer_dp(a,x,its)
       USE mo_kind
       INTEGER(I4), INTENT(OUT) :: its
       COMPLEX(DPC), INTENT(INOUT) :: x
       COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: a
     END SUBROUTINE laguer_dp
  END INTERFACE laguer
  INTERFACE
     SUBROUTINE lfit(x,y,sig,a,maska,covar,chisq,funcs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
       LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
       REAL(SP), INTENT(OUT) :: chisq
       INTERFACE
          SUBROUTINE funcs(x,arr)
            USE mo_kind
            REAL(SP),INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(OUT) :: arr
          END SUBROUTINE funcs
       END INTERFACE
     END SUBROUTINE lfit
  END INTERFACE
  INTERFACE linbcg
     SUBROUTINE linbcg_sp(b,x,itol,tol,itmax,iter,err)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: b
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
       INTEGER(I4), INTENT(IN) :: itol,itmax
       REAL(SP), INTENT(IN) :: tol
       INTEGER(I4), INTENT(OUT) :: iter
       REAL(SP), INTENT(OUT) :: err
     END SUBROUTINE linbcg_sp
     !MC
     SUBROUTINE linbcg_dp(b,x,itol,tol,itmax,iter,err)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: b
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
       INTEGER(I4), INTENT(IN) :: itol,itmax
       REAL(DP), INTENT(IN) :: tol
       INTEGER(I4), INTENT(OUT) :: iter
       REAL(DP), INTENT(OUT) :: err
     END SUBROUTINE linbcg_dp
  END INTERFACE
  INTERFACE
     SUBROUTINE linmin(p,xi,fret)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: fret
       REAL(SP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
     END SUBROUTINE linmin
  END INTERFACE
  INTERFACE lnsrch
     SUBROUTINE lnsrch_sp(xold,fold,g,p,x,f,stpmax,check,func)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: xold,g
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
       REAL(SP), INTENT(IN) :: fold,stpmax
       REAL(SP), DIMENSION(:), INTENT(OUT) :: x
       REAL(SP), INTENT(OUT) :: f
       LOGICAL(LGT), INTENT(OUT) :: check
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP) :: func
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE lnsrch_sp
     SUBROUTINE lnsrch_dp(xold,fold,g,p,x,f,stpmax,check,func)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: xold,g
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
       REAL(DP), INTENT(IN) :: fold,stpmax
       REAL(DP), DIMENSION(:), INTENT(OUT) :: x
       REAL(DP), INTENT(OUT) :: f
       LOGICAL(LGT), INTENT(OUT) :: check
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(DP) :: func
            REAL(DP), DIMENSION(:), INTENT(IN) :: x
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE lnsrch_dp
  END INTERFACE lnsrch
  INTERFACE
     FUNCTION locate(xx,x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: xx
       REAL(SP), INTENT(IN) :: x
       INTEGER(I4) :: locate
     END FUNCTION locate
  END INTERFACE
  INTERFACE
     FUNCTION lop(u)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: u
       REAL(DP), DIMENSION(size(u,1),size(u,1)) :: lop
     END FUNCTION lop
  END INTERFACE
  INTERFACE
     SUBROUTINE lubksb(a,indx,b)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: indx
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
     END SUBROUTINE lubksb
  END INTERFACE
  INTERFACE ludcmp
     SUBROUTINE ludcmp_sp(a,indx,d)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
       INTEGER(I4), DIMENSION(:), INTENT(OUT) :: indx
       REAL(SP), INTENT(OUT) :: d
     END SUBROUTINE ludcmp_sp
     SUBROUTINE ludcmp_dp(a,indx,d)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
       INTEGER(I4), DIMENSION(:), INTENT(OUT) :: indx
       REAL(DP), INTENT(OUT) :: d
     END SUBROUTINE ludcmp_dp
  END INTERFACE ludcmp
  INTERFACE
     SUBROUTINE machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,&
          maxexp,eps,epsneg,xmin,xmax)
       USE mo_kind
       INTEGER(I4), INTENT(OUT) :: ibeta,iexp,irnd,it,machep,maxexp,&
            minexp,negep,ngrd
       REAL(SP), INTENT(OUT) :: eps,epsneg,xmax,xmin
     END SUBROUTINE machar
  END INTERFACE
  INTERFACE
     SUBROUTINE medfit(x,y,a,b,abdev)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
       REAL(SP), INTENT(OUT) :: a,b,abdev
     END SUBROUTINE medfit
  END INTERFACE
  INTERFACE
     SUBROUTINE memcof(data,xms,d)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: xms
       REAL(SP), DIMENSION(:), INTENT(IN) :: data
       REAL(SP), DIMENSION(:), INTENT(OUT) :: d
     END SUBROUTINE memcof
  END INTERFACE
  INTERFACE
     SUBROUTINE mgfas(u,maxcyc)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
       INTEGER(I4), INTENT(IN) :: maxcyc
     END SUBROUTINE mgfas
  END INTERFACE
  INTERFACE
     SUBROUTINE mglin(u,ncycle)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
       INTEGER(I4), INTENT(IN) :: ncycle
     END SUBROUTINE mglin
  END INTERFACE
  INTERFACE
     SUBROUTINE midexp(funk,aa,bb,s,n)
       USE mo_kind
       REAL(SP), INTENT(IN) :: aa,bb
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4), INTENT(IN) :: n
       INTERFACE
          FUNCTION funk(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: funk
          END FUNCTION funk
       END INTERFACE
     END SUBROUTINE midexp
  END INTERFACE
  INTERFACE
     SUBROUTINE midinf(funk,aa,bb,s,n)
       USE mo_kind
       REAL(SP), INTENT(IN) :: aa,bb
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4), INTENT(IN) :: n
       INTERFACE
          FUNCTION funk(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: funk
          END FUNCTION funk
       END INTERFACE
     END SUBROUTINE midinf
  END INTERFACE
  INTERFACE
     SUBROUTINE midpnt(func,a,b,s,n)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4), INTENT(IN) :: n
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE midpnt
  END INTERFACE
  INTERFACE
     SUBROUTINE midsql(funk,aa,bb,s,n)
       USE mo_kind
       REAL(SP), INTENT(IN) :: aa,bb
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4), INTENT(IN) :: n
       INTERFACE
          FUNCTION funk(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: funk
          END FUNCTION funk
       END INTERFACE
     END SUBROUTINE midsql
  END INTERFACE
  INTERFACE
     SUBROUTINE midsqu(funk,aa,bb,s,n)
       USE mo_kind
       REAL(SP), INTENT(IN) :: aa,bb
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4), INTENT(IN) :: n
       INTERFACE
          FUNCTION funk(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: funk
          END FUNCTION funk
       END INTERFACE
     END SUBROUTINE midsqu
  END INTERFACE
  INTERFACE
     RECURSIVE SUBROUTINE miser(func,regn,ndim,npts,dith,ave,var)
       USE mo_kind
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP) :: func
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
          END FUNCTION func
       END INTERFACE
       REAL(SP), DIMENSION(:), INTENT(IN) :: regn
       INTEGER(I4), INTENT(IN) :: ndim,npts
       REAL(SP), INTENT(IN) :: dith
       REAL(SP), INTENT(OUT) :: ave,var
     END SUBROUTINE miser
  END INTERFACE
  INTERFACE
     SUBROUTINE mmid(y,dydx,xs,htot,nstep,yout,derivs)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: nstep
       REAL(SP), INTENT(IN) :: xs,htot
       REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
       REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
       END INTERFACE
     END SUBROUTINE mmid
  END INTERFACE
  INTERFACE mnbrak
     SUBROUTINE mnbrak_sp(ax,bx,cx,fa,fb,fc,func)
       USE mo_kind
       REAL(SP), INTENT(INOUT) :: ax,bx
       REAL(SP), INTENT(OUT) :: cx,fa,fb,fc
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE mnbrak_sp
     SUBROUTINE mnbrak_dp(ax,bx,cx,fa,fb,fc,func)
       USE mo_kind
       REAL(DP), INTENT(INOUT) :: ax,bx
       REAL(DP), INTENT(OUT) :: cx,fa,fb,fc
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(DP), INTENT(IN) :: x
            REAL(DP) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE mnbrak_dp
  END INTERFACE mnbrak
  INTERFACE
     SUBROUTINE mnewt(ntrial,x,tolx,tolf,usrfun)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: ntrial
       REAL(SP), INTENT(IN) :: tolx,tolf
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
       INTERFACE
          SUBROUTINE usrfun(x,fvec,fjac)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(OUT) :: fvec
            REAL(SP), DIMENSION(:,:), INTENT(OUT) :: fjac
          END SUBROUTINE usrfun
       END INTERFACE
     END SUBROUTINE mnewt
  END INTERFACE
  INTERFACE
     SUBROUTINE moment(data,ave,adev,sdev,var,skew,curt)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: ave,adev,sdev,var,skew,curt
       REAL(SP), DIMENSION(:), INTENT(IN) :: data
     END SUBROUTINE moment
  END INTERFACE
  INTERFACE
     SUBROUTINE mp2dfr(a,s,n,m)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       INTEGER(I4), INTENT(OUT) :: m
       CHARACTER(1), DIMENSION(:), INTENT(INOUT) :: a
       CHARACTER(1), DIMENSION(:), INTENT(OUT) :: s
     END SUBROUTINE mp2dfr
  END INTERFACE
  INTERFACE
     SUBROUTINE mpdiv(q,r,u,v,n,m)
       USE mo_kind
       CHARACTER(1), DIMENSION(:), INTENT(OUT) :: q,r
       CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
       INTEGER(I4), INTENT(IN) :: n,m
     END SUBROUTINE mpdiv
  END INTERFACE
  INTERFACE
     SUBROUTINE mpinv(u,v,n,m)
       USE mo_kind
       CHARACTER(1), DIMENSION(:), INTENT(OUT) :: u
       CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
       INTEGER(I4), INTENT(IN) :: n,m
     END SUBROUTINE mpinv
  END INTERFACE
  INTERFACE
     SUBROUTINE mpmul(w,u,v,n,m)
       USE mo_kind
       CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
       CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
       INTEGER(I4), INTENT(IN) :: n,m
     END SUBROUTINE mpmul
  END INTERFACE
  INTERFACE
     SUBROUTINE mppi(n)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
     END SUBROUTINE mppi
  END INTERFACE
  INTERFACE
     SUBROUTINE mprove(a,alud,indx,b,x)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,alud
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: indx
       REAL(SP), DIMENSION(:), INTENT(IN) :: b
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
     END SUBROUTINE mprove
  END INTERFACE
  INTERFACE
     SUBROUTINE mpsqrt(w,u,v,n,m)
       USE mo_kind
       CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w,u
       CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
       INTEGER(I4), INTENT(IN) :: n,m
     END SUBROUTINE mpsqrt
  END INTERFACE
  INTERFACE
     SUBROUTINE mrqcof(x,y,sig,a,maska,alpha,beta,chisq,funcs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,a,sig
       REAL(SP), DIMENSION(:), INTENT(OUT) :: beta
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: alpha
       REAL(SP), INTENT(OUT) :: chisq
       LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
       INTERFACE
          SUBROUTINE funcs(x,a,yfit,dyda)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
            REAL(SP), DIMENSION(:), INTENT(OUT) :: yfit
            REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
          END SUBROUTINE funcs
       END INTERFACE
     END SUBROUTINE mrqcof
  END INTERFACE
  INTERFACE mrqmin
     SUBROUTINE mrqmin_sp(x,y,sig,a,maska,covar,alpha,chisq,funcs,alamda)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: covar,alpha
       REAL(SP), INTENT(OUT) :: chisq
       REAL(SP), INTENT(INOUT) :: alamda
       LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
       INTERFACE
          SUBROUTINE funcs(x,a,yfit,dyda)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
            REAL(SP), DIMENSION(:), INTENT(OUT) :: yfit
            REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
          END SUBROUTINE funcs
       END INTERFACE
     END SUBROUTINE mrqmin_sp
     SUBROUTINE mrqmin_dp(x,y,sig,a,maska,covar,alpha,chisq,funcs,alamda)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: x,y,sig
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: a
       REAL(DP), DIMENSION(:,:), INTENT(OUT) :: covar,alpha
       REAL(DP), INTENT(OUT) :: chisq
       REAL(DP), INTENT(INOUT) :: alamda
       LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
       INTERFACE
          SUBROUTINE funcs(x,a,yfit,dyda)
            USE mo_kind
            REAL(DP), DIMENSION(:), INTENT(IN) :: x,a
            REAL(DP), DIMENSION(:), INTENT(OUT) :: yfit
            REAL(DP), DIMENSION(:,:), INTENT(OUT) :: dyda
          END SUBROUTINE funcs
       END INTERFACE
     END SUBROUTINE mrqmin_dp
  END INTERFACE mrqmin
  INTERFACE
     SUBROUTINE newt(x,check)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
       LOGICAL(LGT), INTENT(OUT) :: check
     END SUBROUTINE newt
  END INTERFACE
  INTERFACE
     SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: ystart
       REAL(SP), INTENT(IN) :: x1,x2,eps,h1,hmin
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
          !BL
          SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
            REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
            REAL(SP), INTENT(INOUT) :: x
            REAL(SP), INTENT(IN) :: htry,eps
            REAL(SP), INTENT(OUT) :: hdid,hnext
            INTERFACE
               SUBROUTINE derivs(x,y,dydx)
                 USE mo_kind
                 REAL(SP), INTENT(IN) :: x
                 REAL(SP), DIMENSION(:), INTENT(IN) :: y
                 REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
               END SUBROUTINE derivs
            END INTERFACE
          END SUBROUTINE rkqs
       END INTERFACE
     END SUBROUTINE odeint
  END INTERFACE
  INTERFACE
     SUBROUTINE orthog(anu,alpha,beta,a,b)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: anu,alpha,beta
       REAL(SP), DIMENSION(:), INTENT(OUT) :: a,b
     END SUBROUTINE orthog
  END INTERFACE
  INTERFACE
     SUBROUTINE pade(cof,resid)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: cof
       REAL(SP), INTENT(OUT) :: resid
     END SUBROUTINE pade
  END INTERFACE
  INTERFACE
     FUNCTION pccheb(d)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: d
       REAL(SP), DIMENSION(size(d)) :: pccheb
     END FUNCTION pccheb
  END INTERFACE
  INTERFACE
     SUBROUTINE pcshft(a,b,d)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
     END SUBROUTINE pcshft
  END INTERFACE
  INTERFACE
     SUBROUTINE pearsn(x,y,r,prob,z)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: r,prob,z
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
     END SUBROUTINE pearsn
  END INTERFACE
  INTERFACE
     SUBROUTINE period(x,y,ofac,hifac,px,py,jmax,prob)
       USE mo_kind
       INTEGER(I4), INTENT(OUT) :: jmax
       REAL(SP), INTENT(IN) :: ofac,hifac
       REAL(SP), INTENT(OUT) :: prob
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
       REAL(SP), DIMENSION(:), POINTER :: px,py
     END SUBROUTINE period
  END INTERFACE
  INTERFACE plgndr
     FUNCTION plgndr_s(l,m,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: l,m
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: plgndr_s
     END FUNCTION plgndr_s
     !BL
     FUNCTION plgndr_v(l,m,x)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: l,m
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: plgndr_v
     END FUNCTION plgndr_v
  END INTERFACE
  INTERFACE
     FUNCTION poidev(xm)
       USE mo_kind
       REAL(SP), INTENT(IN) :: xm
       REAL(SP) :: poidev
     END FUNCTION poidev
  END INTERFACE
  INTERFACE
     FUNCTION polcoe(x,y)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
       REAL(SP), DIMENSION(size(x)) :: polcoe
     END FUNCTION polcoe
  END INTERFACE
  INTERFACE
     FUNCTION polcof(xa,ya)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
       REAL(SP), DIMENSION(size(xa)) :: polcof
     END FUNCTION polcof
  END INTERFACE
  INTERFACE
     SUBROUTINE poldiv(u,v,q,r)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: u,v
       REAL(SP), DIMENSION(:), INTENT(OUT) :: q,r
     END SUBROUTINE poldiv
  END INTERFACE
  INTERFACE
     SUBROUTINE polin2(x1a,x2a,ya,x1,x2,y,dy)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
       REAL(SP), INTENT(IN) :: x1,x2
       REAL(SP), INTENT(OUT) :: y,dy
     END SUBROUTINE polin2
  END INTERFACE
  INTERFACE polint
     SUBROUTINE spolint(xa,ya,x,y,dy)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
       REAL(SP), INTENT(IN) :: x
       REAL(SP), INTENT(OUT) :: y,dy
     END SUBROUTINE spolint
     SUBROUTINE dpolint(xa,ya,x,y,dy)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
       REAL(DP), INTENT(IN) :: x
       REAL(DP), INTENT(OUT) :: y,dy
     END SUBROUTINE dpolint
  END INTERFACE
  INTERFACE
     SUBROUTINE powell(p,xi,ftol,iter,fret)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: xi
       INTEGER(I4), INTENT(OUT) :: iter
       REAL(SP), INTENT(IN) :: ftol
       REAL(SP), INTENT(OUT) :: fret
     END SUBROUTINE powell
  END INTERFACE
  INTERFACE
     FUNCTION predic(data,d,nfut)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: data,d
       INTEGER(I4), INTENT(IN) :: nfut
       REAL(SP), DIMENSION(nfut) :: predic
     END FUNCTION predic
  END INTERFACE
  INTERFACE
     FUNCTION probks(alam)
       USE mo_kind
       REAL(SP), INTENT(IN) :: alam
       REAL(SP) :: probks
     END FUNCTION probks
  END INTERFACE
  INTERFACE psdes
     SUBROUTINE psdes_s(lword,rword)
       USE mo_kind
       INTEGER(I4), INTENT(INOUT) :: lword,rword
     END SUBROUTINE psdes_s
     !BL
     SUBROUTINE psdes_v(lword,rword)
       USE mo_kind
       INTEGER(I4), DIMENSION(:), INTENT(INOUT) :: lword,rword
     END SUBROUTINE psdes_v
  END INTERFACE
  INTERFACE
     SUBROUTINE pwt(a,isign)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE pwt
  END INTERFACE
  INTERFACE
     SUBROUTINE pwtset(n)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
     END SUBROUTINE pwtset
  END INTERFACE
  INTERFACE pythag
     FUNCTION pythag_dp(a,b)
       USE mo_kind
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP) :: pythag_dp
     END FUNCTION pythag_dp
     !BL
     FUNCTION pythag_sp(a,b)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP) :: pythag_sp
     END FUNCTION pythag_sp
  END INTERFACE
  INTERFACE
     SUBROUTINE pzextr(iest,xest,yest,yz,dy)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: iest
       REAL(SP), INTENT(IN) :: xest
       REAL(SP), DIMENSION(:), INTENT(IN) :: yest
       REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
     END SUBROUTINE pzextr
  END INTERFACE
  INTERFACE
     SUBROUTINE qrdcmp(a,c,d,sing)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
       REAL(SP), DIMENSION(:), INTENT(OUT) :: c,d
       LOGICAL(LGT), INTENT(OUT) :: sing
     END SUBROUTINE qrdcmp
  END INTERFACE
  INTERFACE qromb
     FUNCTION qromb_sp(func,a,b,prec)
       USE mo_kind
       INTEGER(I4), INTENT(IN), optional :: prec
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP) :: qromb_sp
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION qromb_sp
     ! JM
     FUNCTION qromb_dp(func,a,b, prec)
       USE mo_kind
       INTEGER(i4), INTENT(in), optional :: prec
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP) :: qromb_dp
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(DP), DIMENSION(:), INTENT(IN) :: x
            REAL(DP), DIMENSION(size(x)) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION qromb_dp
  END INTERFACE
  INTERFACE
     FUNCTION qromo(func,a,b,choose)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP) :: qromo
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
          END FUNCTION func
       END INTERFACE
       INTERFACE
          SUBROUTINE choose(funk,aa,bb,s,n)
            USE mo_kind
            REAL(SP), INTENT(IN) :: aa,bb
            REAL(SP), INTENT(INOUT) :: s
            INTEGER(I4), INTENT(IN) :: n
            INTERFACE
               FUNCTION funk(x)
                 USE mo_kind
                 REAL(SP), DIMENSION(:), INTENT(IN) :: x
                 REAL(SP), DIMENSION(size(x)) :: funk
               END FUNCTION funk
            END INTERFACE
          END SUBROUTINE choose
       END INTERFACE
     END FUNCTION qromo
  END INTERFACE
  INTERFACE
     SUBROUTINE qroot(p,b,c,eps)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: p
       REAL(SP), INTENT(INOUT) :: b,c
       REAL(SP), INTENT(IN) :: eps
     END SUBROUTINE qroot
  END INTERFACE
  INTERFACE
     SUBROUTINE qrsolv(a,c,d,b)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
       REAL(SP), DIMENSION(:), INTENT(IN) :: c,d
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
     END SUBROUTINE qrsolv
  END INTERFACE
  INTERFACE
     SUBROUTINE qrupdt(r,qt,u,v)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: r,qt
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: u
       REAL(SP), DIMENSION(:), INTENT(IN) :: v
     END SUBROUTINE qrupdt
  END INTERFACE
  INTERFACE
     FUNCTION qsimp(func,a,b)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP) :: qsimp
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION qsimp
  END INTERFACE
  INTERFACE
     FUNCTION qtrap(func,a,b)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP) :: qtrap
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION qtrap
  END INTERFACE
  INTERFACE
     SUBROUTINE quadct(x,y,xx,yy,fa,fb,fc,fd)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x,y
       REAL(SP), DIMENSION(:), INTENT(IN) :: xx,yy
       REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
     END SUBROUTINE quadct
  END INTERFACE
  INTERFACE
     SUBROUTINE quadmx(a)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: a
     END SUBROUTINE quadmx
  END INTERFACE
  INTERFACE
     SUBROUTINE quadvl(x,y,fa,fb,fc,fd)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x,y
       REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
     END SUBROUTINE quadvl
  END INTERFACE
  INTERFACE
     ! original from NR --> "ran" überdeckt eingebaute Funktion mit gleichem Namen
     ! FUNCTION ran(idum)
     ! renamed
     FUNCTION ran_nr(idum)
       INTEGER(selected_int_kind(9)), INTENT(INOUT) :: idum
       REAL :: ran_nr
     END FUNCTION ran_nr
  END INTERFACE
  INTERFACE ran0
     SUBROUTINE ran0_s(harvest)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: harvest
     END SUBROUTINE ran0_s
     !BL
     SUBROUTINE ran0_v(harvest)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
     END SUBROUTINE ran0_v
  END INTERFACE
  INTERFACE ran1
     SUBROUTINE ran1_s(harvest)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: harvest
     END SUBROUTINE ran1_s
     !BL
     SUBROUTINE ran1_v(harvest)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
     END SUBROUTINE ran1_v
  END INTERFACE
  INTERFACE ran2
     SUBROUTINE ran2_s(harvest)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: harvest
     END SUBROUTINE ran2_s
     !BL
     SUBROUTINE ran2_v(harvest)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
     END SUBROUTINE ran2_v
  END INTERFACE
  INTERFACE ran3
     SUBROUTINE ran3_s(harvest)
       USE mo_kind
       REAL(SP), INTENT(OUT) :: harvest
     END SUBROUTINE ran3_s
     !BL
     SUBROUTINE ran3_v(harvest)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
     END SUBROUTINE ran3_v
  END INTERFACE
  INTERFACE
     SUBROUTINE ratint(xa,ya,x,y,dy)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
       REAL(SP), INTENT(IN) :: x
       REAL(SP), INTENT(OUT) :: y,dy
     END SUBROUTINE ratint
  END INTERFACE
  INTERFACE
     SUBROUTINE ratlsq(func,a,b,mm,kk,cof,dev)
       USE mo_kind
       REAL(DP), INTENT(IN) :: a,b
       INTEGER(I4), INTENT(IN) :: mm,kk
       REAL(DP), DIMENSION(:), INTENT(OUT) :: cof
       REAL(DP), INTENT(OUT) :: dev
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(DP), DIMENSION(:), INTENT(IN) :: x
            REAL(DP), DIMENSION(size(x)) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE ratlsq
  END INTERFACE
  INTERFACE ratval
     FUNCTION ratval_s(x,cof,mm,kk)
       USE mo_kind
       REAL(DP), INTENT(IN) :: x
       INTEGER(I4), INTENT(IN) :: mm,kk
       REAL(DP), DIMENSION(mm+kk+1), INTENT(IN) :: cof
       REAL(DP) :: ratval_s
     END FUNCTION ratval_s
     !BL
     FUNCTION ratval_v(x,cof,mm,kk)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: x
       INTEGER(I4), INTENT(IN) :: mm,kk
       REAL(DP), DIMENSION(mm+kk+1), INTENT(IN) :: cof
       REAL(DP), DIMENSION(size(x)) :: ratval_v
     END FUNCTION ratval_v
  END INTERFACE
  INTERFACE rc
     FUNCTION rc_s(x,y)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x,y
       REAL(SP) :: rc_s
     END FUNCTION rc_s
     !BL
     FUNCTION rc_v(x,y)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
       REAL(SP), DIMENSION(size(x)) :: rc_v
     END FUNCTION rc_v
  END INTERFACE
  INTERFACE rd
     FUNCTION rd_s(x,y,z)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x,y,z
       REAL(SP) :: rd_s
     END FUNCTION rd_s
     !BL
     FUNCTION rd_v(x,y,z)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z
       REAL(SP), DIMENSION(size(x)) :: rd_v
     END FUNCTION rd_v
  END INTERFACE
  INTERFACE realft
     SUBROUTINE realft_dp(data,isign,zdata)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
       COMPLEX(DPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
     END SUBROUTINE realft_dp
     !BL
     SUBROUTINE realft_sp(data,isign,zdata)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
       INTEGER(I4), INTENT(IN) :: isign
       COMPLEX(SPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
     END SUBROUTINE realft_sp
  END INTERFACE
  INTERFACE
     RECURSIVE FUNCTION recur1(a,b) RESULT(u)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
       REAL(SP), DIMENSION(size(a)) :: u
     END FUNCTION recur1
  END INTERFACE
  INTERFACE
     FUNCTION recur2(a,b,c)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c
       REAL(SP), DIMENSION(size(a)) :: recur2
     END FUNCTION recur2
  END INTERFACE
  INTERFACE
     SUBROUTINE relax(u,rhs)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
     END SUBROUTINE relax
  END INTERFACE
  INTERFACE
     SUBROUTINE relax2(u,rhs)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
     END SUBROUTINE relax2
  END INTERFACE
  INTERFACE
     FUNCTION resid(u,rhs)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,rhs
       REAL(DP), DIMENSION(size(u,1),size(u,1)) :: resid
     END FUNCTION resid
  END INTERFACE
  INTERFACE rf
     FUNCTION rf_s(x,y,z)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x,y,z
       REAL(SP) :: rf_s
     END FUNCTION rf_s
     !BL
     FUNCTION rf_v(x,y,z)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z
       REAL(SP), DIMENSION(size(x)) :: rf_v
     END FUNCTION rf_v
  END INTERFACE
  INTERFACE rj
     FUNCTION rj_s(x,y,z,p)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x,y,z,p
       REAL(SP) :: rj_s
     END FUNCTION rj_s
     !BL
     FUNCTION rj_v(x,y,z,p)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z,p
       REAL(SP), DIMENSION(size(x)) :: rj_v
     END FUNCTION rj_v
  END INTERFACE
  INTERFACE
     SUBROUTINE rk4(y,dydx,x,h,yout,derivs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
       REAL(SP), INTENT(IN) :: x,h
       REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
       END INTERFACE
     END SUBROUTINE rk4
  END INTERFACE
  INTERFACE
     SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
       REAL(SP), INTENT(IN) :: x,h
       REAL(SP), DIMENSION(:), INTENT(OUT) :: yout,yerr
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
       END INTERFACE
     END SUBROUTINE rkck
  END INTERFACE
  INTERFACE
     SUBROUTINE rkdumb(vstart,x1,x2,nstep,derivs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: vstart
       REAL(SP), INTENT(IN) :: x1,x2
       INTEGER(I4), INTENT(IN) :: nstep
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
       END INTERFACE
     END SUBROUTINE rkdumb
  END INTERFACE
  INTERFACE
     SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
       REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
       REAL(SP), INTENT(INOUT) :: x
       REAL(SP), INTENT(IN) :: htry,eps
       REAL(SP), INTENT(OUT) :: hdid,hnext
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
       END INTERFACE
     END SUBROUTINE rkqs
  END INTERFACE
  INTERFACE
     SUBROUTINE rlft2(data,spec,speq,isign)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: data
       COMPLEX(SPC), DIMENSION(:,:), INTENT(OUT) :: spec
       COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: speq
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE rlft2
  END INTERFACE
  INTERFACE
     SUBROUTINE rlft3(data,spec,speq,isign)
       USE mo_kind
       REAL(SP), DIMENSION(:,:,:), INTENT(INOUT) :: data
       COMPLEX(SPC), DIMENSION(:,:,:), INTENT(OUT) :: spec
       COMPLEX(SPC), DIMENSION(:,:), INTENT(OUT) :: speq
       INTEGER(I4), INTENT(IN) :: isign
     END SUBROUTINE rlft3
  END INTERFACE
  INTERFACE
     SUBROUTINE rotate(r,qt,i,a,b)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), TARGET, INTENT(INOUT) :: r,qt
       INTEGER(I4), INTENT(IN) :: i
       REAL(SP), INTENT(IN) :: a,b
     END SUBROUTINE rotate
  END INTERFACE
  INTERFACE rsolv
     SUBROUTINE rsolv_sp(a,d,b)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
       REAL(SP), DIMENSION(:), INTENT(IN) :: d
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
     END SUBROUTINE rsolv_sp
     SUBROUTINE rsolv_dp(a,d,b)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
       REAL(DP), DIMENSION(:), INTENT(IN) :: d
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
     END SUBROUTINE rsolv_dp
  END INTERFACE rsolv
  INTERFACE
     FUNCTION rstrct(uf)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: uf
       REAL(DP), DIMENSION((size(uf,1)+1)/2,(size(uf,1)+1)/2) :: rstrct
     END FUNCTION rstrct
  END INTERFACE
  INTERFACE rtbis
     FUNCTION rtbis_s(func,x1,x2,xacc)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x1,x2,xacc
       REAL(SP) :: rtbis_s
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION rtbis_s
     FUNCTION rtbis_d(func,x1,x2,xacc)
       USE mo_kind
       REAL(DP), INTENT(IN) :: x1,x2,xacc
       REAL(DP) :: rtbis_d
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(DP), INTENT(IN) :: x
            REAL(DP) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION rtbis_d
  END INTERFACE
  INTERFACE
     FUNCTION rtflsp(func,x1,x2,xacc)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x1,x2,xacc
       REAL(SP) :: rtflsp
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION rtflsp
  END INTERFACE
  INTERFACE
     FUNCTION rtnewt(funcd,x1,x2,xacc)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x1,x2,xacc
       REAL(SP) :: rtnewt
       INTERFACE
          SUBROUTINE funcd(x,fval,fderiv)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), INTENT(OUT) :: fval,fderiv
          END SUBROUTINE funcd
       END INTERFACE
     END FUNCTION rtnewt
  END INTERFACE
  INTERFACE
     FUNCTION rtsafe(funcd,x1,x2,xacc)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x1,x2,xacc
       REAL(SP) :: rtsafe
       INTERFACE
          SUBROUTINE funcd(x,fval,fderiv)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), INTENT(OUT) :: fval,fderiv
          END SUBROUTINE funcd
       END INTERFACE
     END FUNCTION rtsafe
  END INTERFACE
  INTERFACE
     FUNCTION rtsec(func,x1,x2,xacc)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x1,x2,xacc
       REAL(SP) :: rtsec
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION rtsec
  END INTERFACE
  INTERFACE
     SUBROUTINE rzextr(iest,xest,yest,yz,dy)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: iest
       REAL(SP), INTENT(IN) :: xest
       REAL(SP), DIMENSION(:), INTENT(IN) :: yest
       REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
     END SUBROUTINE rzextr
  END INTERFACE
  INTERFACE
     FUNCTION savgol(nl,nrr,ld,m)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: nl,nrr,ld,m
       REAL(SP), DIMENSION(nl+nrr+1) :: savgol
     END FUNCTION savgol
  END INTERFACE
  INTERFACE
     SUBROUTINE scrsho(func)
       USE mo_kind
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE scrsho
  END INTERFACE
  INTERFACE
     FUNCTION select(k,arr)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: k
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
       REAL(SP) :: select
     END FUNCTION select
  END INTERFACE
  INTERFACE
     FUNCTION select_bypack(k,arr)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: k
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
       REAL(SP) :: select_bypack
     END FUNCTION select_bypack
  END INTERFACE
  INTERFACE
     SUBROUTINE select_heap(arr,heap)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: arr
       REAL(SP), DIMENSION(:), INTENT(OUT) :: heap
     END SUBROUTINE select_heap
  END INTERFACE
  INTERFACE
     FUNCTION select_inplace(k,arr)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: k
       REAL(SP), DIMENSION(:), INTENT(IN) :: arr
       REAL(SP) :: select_inplace
     END FUNCTION select_inplace
  END INTERFACE
  INTERFACE
     SUBROUTINE simplx(a,m1,m2,m3,icase,izrov,iposv)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
       INTEGER(I4), INTENT(IN) :: m1,m2,m3
       INTEGER(I4), INTENT(OUT) :: icase
       INTEGER(I4), DIMENSION(:), INTENT(OUT) :: izrov,iposv
     END SUBROUTINE simplx
  END INTERFACE
  INTERFACE
     SUBROUTINE simpr(y,dydx,dfdx,dfdy,xs,htot,nstep,yout,derivs)
       USE mo_kind
       REAL(SP), INTENT(IN) :: xs,htot
       REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx,dfdx
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: dfdy
       INTEGER(I4), INTENT(IN) :: nstep
       REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
       END INTERFACE
     END SUBROUTINE simpr
  END INTERFACE
  INTERFACE
     SUBROUTINE sinft(y)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
     END SUBROUTINE sinft
  END INTERFACE
  INTERFACE
     SUBROUTINE slvsm2(u,rhs)
       USE mo_kind
       REAL(DP), DIMENSION(3,3), INTENT(OUT) :: u
       REAL(DP), DIMENSION(3,3), INTENT(IN) :: rhs
     END SUBROUTINE slvsm2
  END INTERFACE
  INTERFACE
     SUBROUTINE slvsml(u,rhs)
       USE mo_kind
       REAL(DP), DIMENSION(3,3), INTENT(OUT) :: u
       REAL(DP), DIMENSION(3,3), INTENT(IN) :: rhs
     END SUBROUTINE slvsml
  END INTERFACE
  INTERFACE
     SUBROUTINE sncndn(uu,emmc,sn,cn,dn)
       USE mo_kind
       REAL(SP), INTENT(IN) :: uu,emmc
       REAL(SP), INTENT(OUT) :: sn,cn,dn
     END SUBROUTINE sncndn
  END INTERFACE
  INTERFACE snrm
     FUNCTION snrm_sp(sx,itol)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: sx
       INTEGER(I4), INTENT(IN) :: itol
       REAL(SP) :: snrm_sp
     END FUNCTION snrm_sp
     !MC
     FUNCTION snrm_dp(sx,itol)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: sx
       INTEGER(I4), INTENT(IN) :: itol
       REAL(DP) :: snrm_dp
     END FUNCTION snrm_dp
  END INTERFACE
  INTERFACE
     SUBROUTINE sobseq(x,init)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(OUT) :: x
       INTEGER(I4), OPTIONAL, INTENT(IN) :: init
     END SUBROUTINE sobseq
  END INTERFACE
  INTERFACE
     SUBROUTINE solvde(itmax,conv,slowc,scalv,indexv,nb,y)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: itmax,nb
       REAL(SP), INTENT(IN) :: conv,slowc
       REAL(SP), DIMENSION(:), INTENT(IN) :: scalv
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: indexv
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: y
     END SUBROUTINE solvde
  END INTERFACE
  INTERFACE
     SUBROUTINE sor(a,b,c,d,e,f,u,rjac)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: a,b,c,d,e,f
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
       REAL(DP), INTENT(IN) :: rjac
     END SUBROUTINE sor
  END INTERFACE
  INTERFACE
     SUBROUTINE sort(arr)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
     END SUBROUTINE sort
  END INTERFACE
  INTERFACE
     SUBROUTINE sort2(arr,slave)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave
     END SUBROUTINE sort2
  END INTERFACE
  INTERFACE
     SUBROUTINE sort3(arr,slave1,slave2)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave1,slave2
     END SUBROUTINE sort3
  END INTERFACE
  INTERFACE
     SUBROUTINE sort_bypack(arr)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
     END SUBROUTINE sort_bypack
  END INTERFACE
  INTERFACE
     SUBROUTINE sort_byreshape(arr)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
     END SUBROUTINE sort_byreshape
  END INTERFACE
  INTERFACE
     SUBROUTINE sort_heap(arr)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
     END SUBROUTINE sort_heap
  END INTERFACE
  INTERFACE
     SUBROUTINE sort_pick(arr)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
     END SUBROUTINE sort_pick
  END INTERFACE
  INTERFACE
     SUBROUTINE sort_radix(arr)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
     END SUBROUTINE sort_radix
  END INTERFACE
  INTERFACE
     SUBROUTINE sort_shell(arr)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
     END SUBROUTINE sort_shell
  END INTERFACE
  INTERFACE
     SUBROUTINE spctrm(p,k,ovrlap,unit,n_window)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(OUT) :: p
       INTEGER(I4), INTENT(IN) :: k
       LOGICAL(LGT), INTENT(IN) :: ovrlap
       INTEGER(I4), OPTIONAL, INTENT(IN) :: n_window,unit
     END SUBROUTINE spctrm
  END INTERFACE
  INTERFACE
     SUBROUTINE spear(data1,data2,d,zd,probd,rs,probrs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
       REAL(SP), INTENT(OUT) :: d,zd,probd,rs,probrs
     END SUBROUTINE spear
  END INTERFACE
  INTERFACE sphbes
     SUBROUTINE sphbes_s(n,x,sj,sy,sjp,syp)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), INTENT(IN) :: x
       REAL(SP), INTENT(OUT) :: sj,sy,sjp,syp
     END SUBROUTINE sphbes_s
     !BL
     SUBROUTINE sphbes_v(n,x,sj,sy,sjp,syp)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(:), INTENT(OUT) :: sj,sy,sjp,syp
     END SUBROUTINE sphbes_v
  END INTERFACE
  INTERFACE
     SUBROUTINE splie2(x1a,x2a,ya,y2a)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: y2a
     END SUBROUTINE splie2
  END INTERFACE
  INTERFACE
     FUNCTION splin2(x1a,x2a,ya,y2a,x1,x2)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya,y2a
       REAL(SP), INTENT(IN) :: x1,x2
       REAL(SP) :: splin2
     END FUNCTION splin2
  END INTERFACE
  INTERFACE
     SUBROUTINE spline(x,y,yp1,ypn,y2)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
       REAL(SP), INTENT(IN) :: yp1,ypn
       REAL(SP), DIMENSION(:), INTENT(OUT) :: y2
     END SUBROUTINE spline
  END INTERFACE
  INTERFACE
     FUNCTION splint(xa,ya,y2a,x)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: splint
     END FUNCTION splint
  END INTERFACE
  INTERFACE sprsax
     SUBROUTINE sprsax_dp(sa,x,b)
       USE mo_kind
       TYPE(sprs2_dp), INTENT(IN) :: sa
       REAL(DP), DIMENSION (:), INTENT(IN) :: x
       REAL(DP), DIMENSION (:), INTENT(OUT) :: b
     END SUBROUTINE sprsax_dp
     !BL
     SUBROUTINE sprsax_sp(sa,x,b)
       USE mo_kind
       TYPE(sprs2_sp), INTENT(IN) :: sa
       REAL(SP), DIMENSION (:), INTENT(IN) :: x
       REAL(SP), DIMENSION (:), INTENT(OUT) :: b
     END SUBROUTINE sprsax_sp
  END INTERFACE
  INTERFACE sprsdiag
     SUBROUTINE sprsdiag_dp(sa,b)
       USE mo_kind
       TYPE(sprs2_dp), INTENT(IN) :: sa
       REAL(DP), DIMENSION(:), INTENT(OUT) :: b
     END SUBROUTINE sprsdiag_dp
     !BL
     SUBROUTINE sprsdiag_sp(sa,b)
       USE mo_kind
       TYPE(sprs2_sp), INTENT(IN) :: sa
       REAL(SP), DIMENSION(:), INTENT(OUT) :: b
     END SUBROUTINE sprsdiag_sp
  END INTERFACE
  INTERFACE sprsin
     SUBROUTINE sprsin_sp(a,thresh,sa)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
       REAL(SP), INTENT(IN) :: thresh
       TYPE(sprs2_sp), INTENT(OUT) :: sa
     END SUBROUTINE sprsin_sp
     !BL
     SUBROUTINE sprsin_dp(a,thresh,sa)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
       REAL(DP), INTENT(IN) :: thresh
       TYPE(sprs2_dp), INTENT(OUT) :: sa
     END SUBROUTINE sprsin_dp
  END INTERFACE
  INTERFACE
     SUBROUTINE sprstp(sa)
       USE mo_kind
       TYPE(sprs2_sp), INTENT(INOUT) :: sa
     END SUBROUTINE sprstp
  END INTERFACE
  INTERFACE sprstx
     SUBROUTINE sprstx_dp(sa,x,b)
       USE mo_kind
       TYPE(sprs2_dp), INTENT(IN) :: sa
       REAL(DP), DIMENSION (:), INTENT(IN) :: x
       REAL(DP), DIMENSION (:), INTENT(OUT) :: b
     END SUBROUTINE sprstx_dp
     !BL
     SUBROUTINE sprstx_sp(sa,x,b)
       USE mo_kind
       TYPE(sprs2_sp), INTENT(IN) :: sa
       REAL(SP), DIMENSION (:), INTENT(IN) :: x
       REAL(SP), DIMENSION (:), INTENT(OUT) :: b
     END SUBROUTINE sprstx_sp
  END INTERFACE
  INTERFACE
     SUBROUTINE stifbs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
       REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
       REAL(SP), INTENT(IN) :: htry,eps
       REAL(SP), INTENT(INOUT) :: x
       REAL(SP), INTENT(OUT) :: hdid,hnext
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
       END INTERFACE
     END SUBROUTINE stifbs
  END INTERFACE
  INTERFACE
     SUBROUTINE stiff(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
       REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
       REAL(SP), INTENT(INOUT) :: x
       REAL(SP), INTENT(IN) :: htry,eps
       REAL(SP), INTENT(OUT) :: hdid,hnext
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
       END INTERFACE
     END SUBROUTINE stiff
  END INTERFACE
  INTERFACE
     SUBROUTINE stoerm(y,d2y,xs,htot,nstep,yout,derivs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: y,d2y
       REAL(SP), INTENT(IN) :: xs,htot
       INTEGER(I4), INTENT(IN) :: nstep
       REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
       INTERFACE
          SUBROUTINE derivs(x,y,dydx)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP), DIMENSION(:), INTENT(IN) :: y
            REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
          END SUBROUTINE derivs
       END INTERFACE
     END SUBROUTINE stoerm
  END INTERFACE
  INTERFACE svbksb
     SUBROUTINE svbksb_dp(u,w,v,b,x)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,v
       REAL(DP), DIMENSION(:), INTENT(IN) :: w,b
       REAL(DP), DIMENSION(:), INTENT(OUT) :: x
     END SUBROUTINE svbksb_dp
     !BL
     SUBROUTINE svbksb_sp(u,w,v,b,x)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: u,v
       REAL(SP), DIMENSION(:), INTENT(IN) :: w,b
       REAL(SP), DIMENSION(:), INTENT(OUT) :: x
     END SUBROUTINE svbksb_sp
  END INTERFACE
  INTERFACE svdcmp
     SUBROUTINE svdcmp_dp(a,w,v)
       USE mo_kind
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
       REAL(DP), DIMENSION(:), INTENT(OUT) :: w
       REAL(DP), DIMENSION(:,:), INTENT(OUT) :: v
     END SUBROUTINE svdcmp_dp
     !BL
     SUBROUTINE svdcmp_sp(a,w,v)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
       REAL(SP), DIMENSION(:), INTENT(OUT) :: w
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
     END SUBROUTINE svdcmp_sp
  END INTERFACE
  INTERFACE
     SUBROUTINE svdfit(x,y,sig,a,v,w,chisq,funcs)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
       REAL(SP), DIMENSION(:), INTENT(OUT) :: a,w
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
       REAL(SP), INTENT(OUT) :: chisq
       INTERFACE
          FUNCTION funcs(x,n)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            INTEGER(I4), INTENT(IN) :: n
            REAL(SP), DIMENSION(n) :: funcs
          END FUNCTION funcs
       END INTERFACE
     END SUBROUTINE svdfit
  END INTERFACE
  INTERFACE
     SUBROUTINE svdvar(v,w,cvm)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(IN) :: v
       REAL(SP), DIMENSION(:), INTENT(IN) :: w
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: cvm
     END SUBROUTINE svdvar
  END INTERFACE
  INTERFACE
     FUNCTION toeplz(r,y)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: r,y
       REAL(SP), DIMENSION(size(y)) :: toeplz
     END FUNCTION toeplz
  END INTERFACE
  INTERFACE
     SUBROUTINE tptest(data1,data2,t,prob)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
       REAL(SP), INTENT(OUT) :: t,prob
     END SUBROUTINE tptest
  END INTERFACE
  INTERFACE
     SUBROUTINE tqli(d,e,z)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: d,e
       REAL(SP), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
     END SUBROUTINE tqli
  END INTERFACE
  INTERFACE trapzd
     SUBROUTINE strapzd(func,a,b,s,n)
       USE mo_kind
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP), INTENT(INOUT) :: s
       INTEGER(I4), INTENT(IN) :: n
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(size(x)) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE strapzd
     SUBROUTINE dtrapzd(func,a,b,s,n)
       USE mo_kind
       REAL(DP), INTENT(IN) :: a,b
       REAL(DP), INTENT(INOUT) :: s
       INTEGER(I4), INTENT(IN) :: n
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(DP), DIMENSION(:), INTENT(IN) :: x
            REAL(DP), DIMENSION(size(x)) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE dtrapzd
  END INTERFACE
  INTERFACE
     SUBROUTINE tred2(a,d,e,novectors)
       USE mo_kind
       REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
       REAL(SP), DIMENSION(:), INTENT(OUT) :: d,e
       LOGICAL(LGT), OPTIONAL, INTENT(IN) :: novectors
     END SUBROUTINE tred2
  END INTERFACE
  !	On a purely serial machine, for greater efficiency, remove
  !	the generic name tridag from the following interface,
  !	and put it on the next one after that.
  INTERFACE 
     RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
       REAL(SP), DIMENSION(:), INTENT(OUT) :: u
     END SUBROUTINE tridag_par
     !MC
     RECURSIVE SUBROUTINE dtridag_par(a,b,c,r,u)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
       REAL(DP), DIMENSION(:), INTENT(OUT) :: u
     END SUBROUTINE dtridag_par
  END INTERFACE
  INTERFACE tridag
     SUBROUTINE tridag_ser(a,b,c,r,u)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
       REAL(SP), DIMENSION(:), INTENT(OUT) :: u
     END SUBROUTINE tridag_ser
     !MC
     SUBROUTINE dtridag_ser(a,b,c,r,u)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
       REAL(DP), DIMENSION(:), INTENT(OUT) :: u
     END SUBROUTINE dtridag_ser
  END INTERFACE
  INTERFACE
     SUBROUTINE ttest(data1,data2,t,prob)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
       REAL(SP), INTENT(OUT) :: t,prob
     END SUBROUTINE ttest
  END INTERFACE
  INTERFACE
     SUBROUTINE tutest(data1,data2,t,prob)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
       REAL(SP), INTENT(OUT) :: t,prob
     END SUBROUTINE tutest
  END INTERFACE
  INTERFACE
     SUBROUTINE twofft(data1,data2,fft1,fft2)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
       COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: fft1,fft2
     END SUBROUTINE twofft
  END INTERFACE
  INTERFACE
     FUNCTION vander(x,q)
       USE mo_kind
       REAL(DP), DIMENSION(:), INTENT(IN) :: x,q
       REAL(DP), DIMENSION(size(x)) :: vander
     END FUNCTION vander
  END INTERFACE
  INTERFACE
     SUBROUTINE vegas(region,func,init,ncall,itmx,nprn,tgral,sd,chi2a)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: region
       INTEGER(I4), INTENT(IN) :: init,ncall,itmx,nprn
       REAL(SP), INTENT(OUT) :: tgral,sd,chi2a
       INTERFACE
          FUNCTION func(pt,wgt)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(IN) :: pt
            REAL(SP), INTENT(IN) :: wgt
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE vegas
  END INTERFACE
  INTERFACE
     SUBROUTINE voltra(t0,h,t,f,g,ak)
       USE mo_kind
       REAL(SP), INTENT(IN) :: t0,h
       REAL(SP), DIMENSION(:), INTENT(OUT) :: t
       REAL(SP), DIMENSION(:,:), INTENT(OUT) :: f
       INTERFACE
          FUNCTION g(t)
            USE mo_kind
            REAL(SP), INTENT(IN) :: t
            REAL(SP), DIMENSION(:), POINTER :: g
          END FUNCTION g
          !BL
          FUNCTION ak(t,s)
            USE mo_kind
            REAL(SP), INTENT(IN) :: t,s
            REAL(SP), DIMENSION(:,:), POINTER :: ak
          END FUNCTION ak
       END INTERFACE
     END SUBROUTINE voltra
  END INTERFACE
  INTERFACE
     SUBROUTINE wt1(a,isign,wtstep)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
       INTEGER(I4), INTENT(IN) :: isign
       INTERFACE
          SUBROUTINE wtstep(a,isign)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
            INTEGER(I4), INTENT(IN) :: isign
          END SUBROUTINE wtstep
       END INTERFACE
     END SUBROUTINE wt1
  END INTERFACE
  INTERFACE
     SUBROUTINE wtn(a,nn,isign,wtstep)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
       INTEGER(I4), DIMENSION(:), INTENT(IN) :: nn
       INTEGER(I4), INTENT(IN) :: isign
       INTERFACE
          SUBROUTINE wtstep(a,isign)
            USE mo_kind
            REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
            INTEGER(I4), INTENT(IN) :: isign
          END SUBROUTINE wtstep
       END INTERFACE
     END SUBROUTINE wtn
  END INTERFACE
  INTERFACE
     FUNCTION wwghts(n,h,kermom)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       REAL(SP), INTENT(IN) :: h
       REAL(SP), DIMENSION(n) :: wwghts
       INTERFACE
          FUNCTION kermom(y,m)
            USE mo_kind
            REAL(DP), INTENT(IN) :: y
            INTEGER(I4), INTENT(IN) :: m
            REAL(DP), DIMENSION(m) :: kermom
          END FUNCTION kermom
       END INTERFACE
     END FUNCTION wwghts
  END INTERFACE
  INTERFACE
     SUBROUTINE zbrac(func,x1,x2,succes)
       USE mo_kind
       REAL(SP), INTENT(INOUT) :: x1,x2
       LOGICAL(LGT), INTENT(OUT) :: succes
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE zbrac
  END INTERFACE
  INTERFACE
     SUBROUTINE zbrak(func,x1,x2,n,xb1,xb2,nb)
       USE mo_kind
       INTEGER(I4), INTENT(IN) :: n
       INTEGER(I4), INTENT(OUT) :: nb
       REAL(SP), INTENT(IN) :: x1,x2
       REAL(SP), DIMENSION(:), POINTER :: xb1,xb2
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END SUBROUTINE zbrak
  END INTERFACE
  INTERFACE
     FUNCTION zbrent(func,x1,x2,tol)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x1,x2,tol
       REAL(SP) :: zbrent
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION zbrent
  END INTERFACE
  INTERFACE
     SUBROUTINE zrhqr(a,rtr,rti)
       USE mo_kind
       REAL(SP), DIMENSION(:), INTENT(IN) :: a
       REAL(SP), DIMENSION(:), INTENT(OUT) :: rtr,rti
     END SUBROUTINE zrhqr
  END INTERFACE
  INTERFACE
     FUNCTION zriddr(func,x1,x2,xacc)
       USE mo_kind
       REAL(SP), INTENT(IN) :: x1,x2,xacc
       REAL(SP) :: zriddr
       INTERFACE
          FUNCTION func(x)
            USE mo_kind
            REAL(SP), INTENT(IN) :: x
            REAL(SP) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION zriddr
  END INTERFACE
  INTERFACE zroots
     SUBROUTINE zroots_sp(a,roots,polish)
       USE mo_kind
       COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
       COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: roots
       LOGICAL(lgt), INTENT(IN) :: polish
     END SUBROUTINE zroots_sp
     SUBROUTINE zroots_dp(a,roots,polish)
       USE mo_kind
       COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: a
       COMPLEX(DPC), DIMENSION(:), INTENT(OUT) :: roots
       LOGICAL(lgt), INTENT(IN) :: polish
     END SUBROUTINE zroots_dp
  END INTERFACE zroots
END MODULE mo_nr

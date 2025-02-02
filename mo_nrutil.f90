MODULE mo_nrutil

  ! This module provides common utilities for the numerical recipes f90 routines.

  ! WH Press, SA Teukolsky, WT Vetterling, BP Flannery,
  !    Numerical Recipes in Fortran 90 - The Art of Parallel Scientific Computing, 2nd Edition
  !    Volume 2 of Fortran Numerical Recipes, Cambridge University Press, UK, 1996

  ! Written Numerical Recipes, 1996
  ! Modified Matthias Cuntz, Nov 2011 - adapted to JAMS Fortran structure
  !                                   - added some double precision and 2d array routines

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 1996-2011 Numerical Recipes, Matthias Cuntz - mc (at) macu (dot) de
  !
  ! Permission is hereby granted, free of charge, to any person obtaining a copy
  ! of this software and associated documentation files (the "Software"), to deal
  ! in the Software without restriction, including without limitation the rights
  ! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  ! copies of the Software, and to permit persons to whom the Software is
  ! furnished to do so, subject to the following conditions:
  !
  ! The above copyright notice and this permission notice shall be included in all
  ! copies or substantial portions of the Software.
  !
  ! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  ! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  ! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  ! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  ! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  ! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  ! SOFTWARE.

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

  USE mo_kind

  IMPLICIT NONE

  INTEGER(I4), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  INTEGER(I4), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
  INTEGER(I4), PARAMETER :: NPAR_CUMSUM=16
  INTEGER(I4), PARAMETER :: NPAR_CUMPROD=8
  INTEGER(I4), PARAMETER :: NPAR_POLY=8
  INTEGER(I4), PARAMETER :: NPAR_POLYTERM=8

  ! -------------------------------------------------------------
  ! interfaces sorted alphabetically
  ! -------------------------------------------------------------
  INTERFACE array_copy
     MODULE PROCEDURE array_copy_r, array_copy_d, array_copy_i
  END INTERFACE array_copy
  INTERFACE arth
     MODULE PROCEDURE arth_r, arth_d, arth_i
  END INTERFACE arth
  INTERFACE assert
     MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
  END INTERFACE assert
  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE assert_eq

  INTERFACE cumsum
     MODULE PROCEDURE cumsum_r,cumsum_i
  END INTERFACE cumsum

  INTERFACE diagadd
     MODULE PROCEDURE diagadd_rvsp,diagadd_rsp,diagadd_rvdp,diagadd_rdp
  END INTERFACE diagadd
  INTERFACE diagmult
     MODULE PROCEDURE diagmult_rvsp,diagmult_rsp, diagmult_rvdp,diagmult_rdp
  END INTERFACE diagmult

  INTERFACE geop
     MODULE PROCEDURE geop_r, geop_d, geop_i, geop_c, geop_dv
  END INTERFACE geop
  INTERFACE get_diag
     MODULE PROCEDURE get_diag_rv, get_diag_dv
  END INTERFACE get_diag

  INTERFACE imaxloc
     MODULE PROCEDURE imaxloc_r,imaxloc_i
  END INTERFACE imaxloc
  INTERFACE iminloc
     MODULE PROCEDURE iminloc_s,iminloc_d
  END INTERFACE iminloc

  INTERFACE outerdiff
     MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
  END INTERFACE outerdiff
  INTERFACE outerprod
     MODULE PROCEDURE outerprod_r,outerprod_d
  END INTERFACE outerprod

  INTERFACE poly
     MODULE PROCEDURE poly_rr,poly_rrv,poly_dd,poly_ddv,&
          poly_rc,poly_ddc,poly_cc,poly_dcdc,poly_msk_rrv,poly_msk_ddv
  END INTERFACE poly
  INTERFACE poly_term
     MODULE PROCEDURE poly_term_rr, poly_term_dd, poly_term_cc, poly_term_dcdc
  END INTERFACE poly_term
  INTERFACE put_diag
     MODULE PROCEDURE put_diag_rv, put_diag_r
  END INTERFACE put_diag

  INTERFACE reallocate
     MODULE PROCEDURE reallocate_rv,reallocate_rm,&
          reallocate_iv,reallocate_im,reallocate_hv, &
          reallocate_dv, reallocate_dm
  END INTERFACE reallocate

  INTERFACE scatter_add
     MODULE PROCEDURE scatter_add_r,scatter_add_d
  END INTERFACE scatter_add
  INTERFACE scatter_max
     MODULE PROCEDURE scatter_max_r,scatter_max_d
  END INTERFACE scatter_max
  INTERFACE swap
     MODULE PROCEDURE swap_i,  swap_r,  swap_d,  swap_c,  swap_z, &
          swap_iv, swap_rv, swap_dv, swap_cv, swap_zv, &
          swap_im, swap_rm, swap_dm, swap_cm, swap_zm, &
          masked_swap_is, masked_swap_iv, masked_swap_im, &
          masked_swap_rs, masked_swap_rv, masked_swap_rm, &
          masked_swap_ds, masked_swap_dv, masked_swap_dm, &
          masked_swap_cs, masked_swap_cv, masked_swap_cm, &
          masked_swap_zs, masked_swap_zv, masked_swap_zm
  END INTERFACE swap

  INTERFACE unit_matrix
     MODULE PROCEDURE unit_matrix_sp, unit_matrix_dp
  END INTERFACE unit_matrix

  INTERFACE vabs
     MODULE PROCEDURE vabs_sp, vabs_dp
  END INTERFACE vabs

CONTAINS
  !BL
  SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
    REAL(SP), DIMENSION(:), INTENT(IN) :: src
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER(I4), INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_r
  !BL
  SUBROUTINE array_copy_d(src,dest,n_copied,n_not_copied)
    REAL(DP), DIMENSION(:), INTENT(IN) :: src
    REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER(I4), INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_d
  !BL
  SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
    INTEGER(I4), DIMENSION(:), INTENT(IN) :: src
    INTEGER(I4), DIMENSION(:), INTENT(OUT) :: dest
    INTEGER(I4), INTENT(OUT) :: n_copied, n_not_copied
    n_copied=min(size(src),size(dest))
    n_not_copied=size(src)-n_copied
    dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_i
  !BL
  !BL
  SUBROUTINE swap_i(a,b)
    INTEGER(I4), INTENT(INOUT) :: a,b
    INTEGER(I4) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i
  !BL
  SUBROUTINE swap_r(a,b)
    REAL(SP), INTENT(INOUT) :: a,b
    REAL(SP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r
  !MC
  SUBROUTINE swap_d(a,b)
    REAL(DP), INTENT(INOUT) :: a,b
    REAL(DP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_d
  !BL
  SUBROUTINE swap_c(a,b)
    COMPLEX(SPC), INTENT(INOUT) :: a,b
    COMPLEX(SPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_c
  !BL
  SUBROUTINE swap_z(a,b)
    COMPLEX(DPC), INTENT(INOUT) :: a,b
    COMPLEX(DPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_z
  !MC
  SUBROUTINE swap_iv(a,b)
    INTEGER(I4), DIMENSION(:), INTENT(INOUT) :: a,b
    INTEGER(I4), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_iv
  !BL
  SUBROUTINE swap_rv(a,b)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv
  !MC
  SUBROUTINE swap_dv(a,b)
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(DP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_dv
  !BL
  SUBROUTINE swap_cv(a,b)
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cv
  !BL
  SUBROUTINE swap_zv(a,b)
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zv
  !MC
  SUBROUTINE swap_im(a,b)
    INTEGER(I4), DIMENSION(:,:), INTENT(INOUT) :: a,b
    INTEGER(I4), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_im
  !MC
  SUBROUTINE swap_rm(a,b)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rm
  !MC
  SUBROUTINE swap_dm(a,b)
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    REAL(DP), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_dm
  !BL
  SUBROUTINE swap_cm(a,b)
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cm
  !BL
  SUBROUTINE swap_zm(a,b)
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zm
  !MC
  SUBROUTINE masked_swap_is(a,b,mask)
    INTEGER(I4), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    INTEGER(I4) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_is
  !MC
  SUBROUTINE masked_swap_iv(a,b,mask)
    INTEGER(I4), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    INTEGER(I4), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_iv
  !MC
  SUBROUTINE masked_swap_im(a,b,mask)
    INTEGER(I4), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
    INTEGER(I4), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_im
  !BL
  SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(SP), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    REAL(SP) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_rs
  !BL
  SUBROUTINE masked_swap_rv(a,b,mask)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rv
  !BL
  SUBROUTINE masked_swap_rm(a,b,mask)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rm
  !MC
  SUBROUTINE masked_swap_ds(a,b,mask)
    REAL(DP), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    REAL(DP) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_ds
  !MC
  SUBROUTINE masked_swap_dv(a,b,mask)
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    REAL(DP), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_dv
  !MC
  SUBROUTINE masked_swap_dm(a,b,mask)
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
    REAL(DP), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_dm
  !MC
  SUBROUTINE masked_swap_cs(a,b,mask)
    COMPLEX(SPC), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    COMPLEX(SPC) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_cs
  !MC
  SUBROUTINE masked_swap_cv(a,b,mask)
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    COMPLEX(SPC), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_cv
  !MC
  SUBROUTINE masked_swap_cm(a,b,mask)
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
    COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_cm
  !MC
  SUBROUTINE masked_swap_zs(a,b,mask)
    COMPLEX(DPC), INTENT(INOUT) :: a,b
    LOGICAL, INTENT(IN) :: mask
    COMPLEX(DPC) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_zs
  !MC
  SUBROUTINE masked_swap_zv(a,b,mask)
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    COMPLEX(DPC), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_zv
  !MC
  SUBROUTINE masked_swap_zm(a,b,mask)
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
    COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_zm
  !BL
  !BL
  !GD
  FUNCTION reallocate_dv(p,n)
    REAL(DP), DIMENSION(:), POINTER :: p, reallocate_dv
    INTEGER(I4), INTENT(IN) :: n
    INTEGER(I4) :: nold,ierr
    allocate(reallocate_dv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_dv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_dv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_dv
  !GD
  FUNCTION reallocate_dm(p,n,m)
    REAL(DP), DIMENSION(:,:), POINTER :: p, reallocate_dm
    INTEGER(I4), INTENT(IN) :: n,m
    INTEGER(I4) :: nold,mold,ierr
    allocate(reallocate_dm(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_dm: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_dm(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  END FUNCTION reallocate_dm
  !BL
  FUNCTION reallocate_rv(p,n)
    REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
    INTEGER(I4), INTENT(IN) :: n
    INTEGER(I4) :: nold,ierr
    allocate(reallocate_rv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_rv
  !BL
  FUNCTION reallocate_iv(p,n)
    INTEGER(I4), DIMENSION(:), POINTER :: p, reallocate_iv
    INTEGER(I4), INTENT(IN) :: n
    INTEGER(I4) :: nold,ierr
    allocate(reallocate_iv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_iv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_iv
  !BL
  FUNCTION reallocate_hv(p,n)
    CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
    INTEGER(I4), INTENT(IN) :: n
    INTEGER(I4) :: nold,ierr
    allocate(reallocate_hv(n),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_hv: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p)
    reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
    deallocate(p)
  END FUNCTION reallocate_hv
  !BL
  FUNCTION reallocate_rm(p,n,m)
    REAL(SP), DIMENSION(:,:), POINTER :: p, reallocate_rm
    INTEGER(I4), INTENT(IN) :: n,m
    INTEGER(I4) :: nold,mold,ierr
    allocate(reallocate_rm(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_rm: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_rm(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  END FUNCTION reallocate_rm
  !BL
  FUNCTION reallocate_im(p,n,m)
    INTEGER(I4), DIMENSION(:,:), POINTER :: p, reallocate_im
    INTEGER(I4), INTENT(IN) :: n,m
    INTEGER(I4) :: nold,mold,ierr
    allocate(reallocate_im(n,m),stat=ierr)
    if (ierr /= 0) call &
         nrerror('reallocate_im: problem in attempt to allocate memory')
    if (.not. associated(p)) RETURN
    nold=size(p,1)
    mold=size(p,2)
    reallocate_im(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
  END FUNCTION reallocate_im
  !BL
  FUNCTION ifirstloc(mask)
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    INTEGER(I4) :: ifirstloc
    INTEGER(I4), DIMENSION(1) :: loc
    loc=maxloc(merge(1,0,mask))
    ifirstloc=loc(1)
    if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
  END FUNCTION ifirstloc
  !BL
  FUNCTION imaxloc_r(arr)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4) :: imaxloc_r
    INTEGER(I4), DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc_r=imax(1)
  END FUNCTION imaxloc_r
  !BL
  FUNCTION imaxloc_i(iarr)
    INTEGER(I4), DIMENSION(:), INTENT(IN) :: iarr
    INTEGER(I4), DIMENSION(1) :: imax
    INTEGER(I4) :: imaxloc_i
    imax=maxloc(iarr(:))
    imaxloc_i=imax(1)
  END FUNCTION imaxloc_i
  !BL
  FUNCTION iminloc_s(arr)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4), DIMENSION(1) :: imin
    INTEGER(I4) :: iminloc_s
    imin=minloc(arr(:))
    iminloc_s=imin(1)
  END FUNCTION iminloc_s
  !JM
  FUNCTION iminloc_d(arr)
    REAL(DP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4), DIMENSION(1) :: imin
    INTEGER(I4) :: iminloc_d
    imin=minloc(arr(:))
    iminloc_d=imin(1)
  END FUNCTION iminloc_d
  !BL
  SUBROUTINE assert1(n1,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1
    if (.not. n1) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert1'
    end if
  END SUBROUTINE assert1
  !BL
  SUBROUTINE assert2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2
    if (.not. (n1 .and. n2)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert2'
    end if
  END SUBROUTINE assert2
  !BL
  SUBROUTINE assert3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3
    if (.not. (n1 .and. n2 .and. n3)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert3'
    end if
  END SUBROUTINE assert3
  !BL
  SUBROUTINE assert4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1,n2,n3,n4
    if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert4'
    end if
  END SUBROUTINE assert4
  !BL
  SUBROUTINE assert_v(n,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, DIMENSION(:), INTENT(IN) :: n
    if (.not. all(n)) then
       write (*,*) 'nrerror: an assertion failed with this tag:', &
            string
       STOP 'program terminated by assert_v'
    end if
  END SUBROUTINE assert_v
  !BL
  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2
  !BL
  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3
  !BL
  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4
  !BL
  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn
  !BL
  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror
  !BL
  FUNCTION arth_r(first,increment,n)
    REAL(SP), INTENT(IN) :: first,increment
    INTEGER(I4), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: arth_r
    INTEGER(I4) :: k,k2
    REAL(SP) :: temp
    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_r(k)=arth_r(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_r
  !BL
  FUNCTION arth_d(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER(I4) :: k,k2
    REAL(DP) :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_d(k)=arth_d(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_d(k)=arth_d(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_d
  !BL
  FUNCTION arth_i(first,increment,n)
    INTEGER(I4), INTENT(IN) :: first,increment,n
    INTEGER(I4), DIMENSION(n) :: arth_i
    INTEGER(I4) :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_i
  !BL
  !BL
  FUNCTION geop_r(first,factor,n)
    REAL(SP), INTENT(IN) :: first,factor
    INTEGER(I4), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: geop_r
    INTEGER(I4) :: k,k2
    REAL(SP) :: temp
    if (n > 0) geop_r(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_r(k)=geop_r(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_r(k)=geop_r(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_r(k+1:min(k2,n))=temp*geop_r(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_r
  !BL
  FUNCTION geop_d(first,factor,n)
    REAL(DP), INTENT(IN) :: first,factor
    INTEGER(I4), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: geop_d
    INTEGER(I4) :: k,k2
    REAL(DP) :: temp
    if (n > 0) geop_d(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_d(k)=geop_d(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_d(k)=geop_d(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_d(k+1:min(k2,n))=temp*geop_d(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_d
  !BL
  FUNCTION geop_i(first,factor,n)
    INTEGER(I4), INTENT(IN) :: first,factor,n
    INTEGER(I4), DIMENSION(n) :: geop_i
    INTEGER(I4) :: k,k2,temp
    if (n > 0) geop_i(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_i(k)=geop_i(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_i(k)=geop_i(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_i(k+1:min(k2,n))=temp*geop_i(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_i
  !BL
  FUNCTION geop_c(first,factor,n)
    COMPLEX(SP), INTENT(IN) :: first,factor
    INTEGER(I4), INTENT(IN) :: n
    COMPLEX(SP), DIMENSION(n) :: geop_c
    INTEGER(I4) :: k,k2
    COMPLEX(SP) :: temp
    if (n > 0) geop_c(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_c(k)=geop_c(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop_c(k)=geop_c(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_c(k+1:min(k2,n))=temp*geop_c(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_c
  !BL
  FUNCTION geop_dv(first,factor,n)
    REAL(DP), DIMENSION(:), INTENT(IN) :: first,factor
    INTEGER(I4), INTENT(IN) :: n
    REAL(DP), DIMENSION(size(first),n) :: geop_dv
    INTEGER(I4) :: k,k2
    REAL(DP), DIMENSION(size(first)) :: temp
    if (n > 0) geop_dv(:,1)=first(:)
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
       end do
    else
       do k=2,NPAR2_GEOP
          geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop_dv(:,k+1:min(k2,n))=geop_dv(:,1:min(k,n-k))*&
               spread(temp,2,size(geop_dv(:,1:min(k,n-k)),2))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop_dv
  !BL
  !BL
  RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    REAL(SP), OPTIONAL, INTENT(IN) :: seed
    REAL(SP), DIMENSION(size(arr)) :: ans
    INTEGER(I4) :: n,j
    REAL(SP) :: sd
    n=size(arr)
    if (n == 0_i4) RETURN
    sd=0.0_sp
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  END FUNCTION cumsum_r
  !BL
  RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
    INTEGER(I4), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4), OPTIONAL, INTENT(IN) :: seed
    INTEGER(I4), DIMENSION(size(arr)) :: ans
    INTEGER(I4) :: n,j,sd
    n=size(arr)
    if (n == 0_i4) RETURN
    sd=0_i4
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  END FUNCTION cumsum_i
  !BL
  !BL
  RECURSIVE FUNCTION cumprod(arr,seed) RESULT(ans)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    REAL(SP), OPTIONAL, INTENT(IN) :: seed
    REAL(SP), DIMENSION(size(arr)) :: ans
    INTEGER(I4) :: n,j
    REAL(SP) :: sd
    n=size(arr)
    if (n == 0_i4) RETURN
    sd=1.0_sp
    if (present(seed)) sd=seed
    ans(1)=arr(1)*sd
    if (n < NPAR_CUMPROD) then
       do j=2,n
          ans(j)=ans(j-1)*arr(j)
       end do
    else
       ans(2:n:2)=cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)*arr(3:n:2)
    end if
  END FUNCTION cumprod
  !BL
  !BL
  FUNCTION poly_rr(x,coeffs)
    REAL(SP), INTENT(IN) :: x
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
    REAL(SP) :: poly_rr
    REAL(SP) :: pow
    REAL(SP), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_rr=0.0_sp
    else if (n < NPAR_POLY) then
       poly_rr=coeffs(n)
       do i=n-1,1,-1
          poly_rr=x*poly_rr+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_sp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_rr=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_rr
  !BL
  FUNCTION poly_dd(x,coeffs)
    REAL(DP), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs
    REAL(DP) :: poly_dd
    REAL(DP) :: pow
    REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_dd=0.0_dp
    else if (n < NPAR_POLY) then
       poly_dd=coeffs(n)
       do i=n-1,1,-1
          poly_dd=x*poly_dd+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_dp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_dd=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_dd
  !BL
  FUNCTION poly_rc(x,coeffs)
    COMPLEX(SPC), INTENT(IN) :: x
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
    COMPLEX(SPC) :: poly_rc
    COMPLEX(SPC) :: pow
    COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_rc=0.0_sp
    else if (n < NPAR_POLY) then
       poly_rc=coeffs(n)
       do i=n-1,1,-1
          poly_rc=x*poly_rc+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_sp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_rc=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_rc
  !
  FUNCTION poly_ddc(x,coeffs)
    COMPLEX(DPC), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs
    COMPLEX(DPC) :: poly_ddc
    COMPLEX(DPC) :: pow
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_ddc=0.0_dp
    else if (n < NPAR_POLY) then
       poly_ddc=coeffs(n)
       do i=n-1,1,-1
          poly_ddc=x*poly_ddc+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_dp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_ddc=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_ddc
  !BL
  FUNCTION poly_cc(x,coeffs)
    COMPLEX(SPC), INTENT(IN) :: x
    COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs
    COMPLEX(SPC) :: poly_cc
    COMPLEX(SPC) :: pow
    COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_cc=0.0_sp
    else if (n < NPAR_POLY) then
       poly_cc=coeffs(n)
       do i=n-1,1,-1
          poly_cc=x*poly_cc+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_sp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_cc=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_cc
  !
  FUNCTION poly_dcdc(x,coeffs)
    COMPLEX(DPC), INTENT(IN) :: x
    COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: coeffs
    COMPLEX(DPC) :: poly_dcdc
    COMPLEX(DPC) :: pow
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER(I4) :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
       poly_dcdc=0.0_dp
    else if (n < NPAR_POLY) then
       poly_dcdc=coeffs(n)
       do i=n-1,1,-1
          poly_dcdc=x*poly_dcdc+coeffs(i)
       end do
    else
       allocate(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       do
          vec(n+1)=0.0_dp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          if (nn == 1) exit
          pow=pow*pow
          n=nn
       end do
       poly_dcdc=vec(1)
       deallocate(vec)
    end if
  END FUNCTION poly_dcdc
  !BL
  FUNCTION poly_rrv(x,coeffs)
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
    REAL(SP), DIMENSION(size(x)) :: poly_rrv
    INTEGER(I4) :: i,n,m
    m=size(coeffs)
    n=size(x)
    if (m <= 0) then
       poly_rrv=0.0_sp
    else if (m < n .or. m < NPAR_POLY) then
       poly_rrv=coeffs(m)
       do i=m-1,1,-1
          poly_rrv=x*poly_rrv+coeffs(i)
       end do
    else
       do i=1,n
          poly_rrv(i)=poly_rr(x(i),coeffs)
       end do
    end if
  END FUNCTION poly_rrv
  !BL
  FUNCTION poly_ddv(x,coeffs)
    REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
    REAL(DP), DIMENSION(size(x)) :: poly_ddv
    INTEGER(I4) :: i,n,m
    m=size(coeffs)
    n=size(x)
    if (m <= 0) then
       poly_ddv=0.0_dp
    else if (m < n .or. m < NPAR_POLY) then
       poly_ddv=coeffs(m)
       do i=m-1,1,-1
          poly_ddv=x*poly_ddv+coeffs(i)
       end do
    else
       do i=1,n
          poly_ddv(i)=poly_dd(x(i),coeffs)
       end do
    end if
  END FUNCTION poly_ddv
  !BL
  FUNCTION poly_msk_rrv(x,coeffs,mask)
    REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(x)) :: poly_msk_rrv
    poly_msk_rrv=unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0_sp)
  END FUNCTION poly_msk_rrv
  !BL
  FUNCTION poly_msk_ddv(x,coeffs,mask)
    REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    REAL(DP), DIMENSION(size(x)) :: poly_msk_ddv
    poly_msk_ddv=unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_dp)
  END FUNCTION poly_msk_ddv
  !BL
  !BL
  RECURSIVE FUNCTION poly_term_rr(a,b) RESULT(u)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a
    REAL(SP), INTENT(IN) :: b
    REAL(SP), DIMENSION(size(a)) :: u
    INTEGER(I4) :: n,j
    n=size(a)
    if (n <= 0) RETURN
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  END FUNCTION poly_term_rr
  !
  RECURSIVE FUNCTION poly_term_dd(a,b) RESULT(u)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a
    REAL(DP), INTENT(IN) :: b
    REAL(DP), DIMENSION(size(a)) :: u
    INTEGER(I4) :: n,j
    n=size(a)
    if (n <= 0) RETURN
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term_dd(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  END FUNCTION poly_term_dd
  !BL
  RECURSIVE FUNCTION poly_term_cc(a,b) RESULT(u)
    COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
    COMPLEX(SPC), INTENT(IN) :: b
    COMPLEX(SPC), DIMENSION(size(a)) :: u
    INTEGER(I4) :: n,j
    n=size(a)
    if (n <= 0) RETURN
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  END FUNCTION poly_term_cc
  !
  RECURSIVE FUNCTION poly_term_dcdc(a,b) RESULT(u)
    COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: a
    COMPLEX(DPC), INTENT(IN) :: b
    COMPLEX(DPC), DIMENSION(size(a)) :: u
    INTEGER(I4) :: n,j
    n=size(a)
    if (n <= 0) RETURN
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
       do j=2,n
          u(j)=a(j)+b*u(j-1)
       end do
    else
       u(2:n:2)=poly_term_dcdc(a(2:n:2)+a(1:n-1:2)*b,b*b)
       u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  END FUNCTION poly_term_dcdc
  !BL
  !BL
  FUNCTION zroots_unity(n,nn)
    USE mo_constants, only: TWOPI
    INTEGER(I4), INTENT(IN) :: n,nn
    COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
    INTEGER(I4) :: k
    REAL(SP) :: theta
    zroots_unity(1)=1.0
    theta=TWOPI/n
    k=1
    do
       if (k >= nn) exit
       zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
       zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
            zroots_unity(2:min(k,nn-k))
       k=2*k
    end do
  END FUNCTION zroots_unity
  !BL
  FUNCTION outerprod_r(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outerprod_r
    outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_r
  !BL
  FUNCTION outerprod_d(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_d
    outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_d
  !BL
  FUNCTION outerdiv(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outerdiv
    outerdiv = spread(a,dim=2,ncopies=size(b)) / &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiv
  !BL
  FUNCTION outersum(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outersum
    outersum = spread(a,dim=2,ncopies=size(b)) + &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outersum
  !BL
  FUNCTION outerdiff_r(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outerdiff_r
    outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_r
  !BL
  FUNCTION outerdiff_d(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff_d
    outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_d
  !BL
  FUNCTION outerdiff_i(a,b)
    INTEGER(I4), DIMENSION(:), INTENT(IN) :: a,b
    INTEGER(I4), DIMENSION(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_i
  !BL
  FUNCTION outerand(a,b)
    LOGICAL, DIMENSION(:), INTENT(IN) :: a,b
    LOGICAL, DIMENSION(size(a),size(b)) :: outerand
    outerand = spread(a,dim=2,ncopies=size(b)) .and. &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerand
  !BL
  SUBROUTINE scatter_add_r(dest,source,dest_index)
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
    REAL(SP), DIMENSION(:), INTENT(IN) :: source
    INTEGER(I4), DIMENSION(:), INTENT(IN) :: dest_index
    INTEGER(I4) :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_add_r')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
    end do
  END SUBROUTINE scatter_add_r
  !BL
  SUBROUTINE scatter_add_d(dest,source,dest_index)
    REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
    REAL(DP), DIMENSION(:), INTENT(IN) :: source
    INTEGER(I4), DIMENSION(:), INTENT(IN) :: dest_index
    INTEGER(I4) :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_add_d')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
    end do
  END SUBROUTINE scatter_add_d
  !BL
  SUBROUTINE scatter_max_r(dest,source,dest_index)
    REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
    REAL(SP), DIMENSION(:), INTENT(IN) :: source
    INTEGER(I4), DIMENSION(:), INTENT(IN) :: dest_index
    INTEGER(I4) :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_max_r')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
    end do
  END SUBROUTINE scatter_max_r
  !BL
  SUBROUTINE scatter_max_d(dest,source,dest_index)
    REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
    REAL(DP), DIMENSION(:), INTENT(IN) :: source
    INTEGER(I4), DIMENSION(:), INTENT(IN) :: dest_index
    INTEGER(I4) :: m,n,j,i
    n=assert_eq2(size(source),size(dest_index),'scatter_max_d')
    m=size(dest)
    do j=1,n
       i=dest_index(j)
       if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
    end do
  END SUBROUTINE scatter_max_d
  !BL
  SUBROUTINE diagadd_rvsp(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), DIMENSION(:), INTENT(IN) :: diag
    INTEGER(I4) :: j,n
    n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rvsp')
    do j=1,n
       mat(j,j)=mat(j,j)+diag(j)
    end do
  END SUBROUTINE diagadd_rvsp
  !BL
  SUBROUTINE diagadd_rsp(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), INTENT(IN) :: diag
    INTEGER(I4) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)+diag
    end do
  END SUBROUTINE diagadd_rsp
  !BL
  SUBROUTINE diagadd_rvdp(mat,diag)
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(DP), DIMENSION(:), INTENT(IN) :: diag
    INTEGER(I4) :: j,n
    n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rvdp')
    do j=1,n
       mat(j,j)=mat(j,j)+diag(j)
    end do
  END SUBROUTINE diagadd_rvdp
  !BL
  SUBROUTINE diagadd_rdp(mat,diag)
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(DP), INTENT(IN) :: diag
    INTEGER(I4) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)+diag
    end do
  END SUBROUTINE diagadd_rdp
  !BL
  SUBROUTINE diagmult_rvsp(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), DIMENSION(:), INTENT(IN) :: diag
    INTEGER(I4) :: j,n
    n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
    do j=1,n
       mat(j,j)=mat(j,j)*diag(j)
    end do
  END SUBROUTINE diagmult_rvsp
  !BL
  SUBROUTINE diagmult_rsp(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), INTENT(IN) :: diag
    INTEGER(I4) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)*diag
    end do
  END SUBROUTINE diagmult_rsp
  !BL
  SUBROUTINE diagmult_rvdp(mat,diag)
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(DP), DIMENSION(:), INTENT(IN) :: diag
    INTEGER(I4) :: j,n
    n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
    do j=1,n
       mat(j,j)=mat(j,j)*diag(j)
    end do
  END SUBROUTINE diagmult_rvdp
  !BL
  SUBROUTINE diagmult_rdp(mat,diag)
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(DP), INTENT(IN) :: diag
    INTEGER(I4) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)*diag
    end do
  END SUBROUTINE diagmult_rdp
  !BL
  FUNCTION get_diag_rv(mat)
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: mat
    REAL(SP), DIMENSION(size(mat,1)) :: get_diag_rv
    INTEGER(I4) :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
    do j=1,size(mat,1)
       get_diag_rv(j)=mat(j,j)
    end do
  END FUNCTION get_diag_rv
  !BL
  FUNCTION get_diag_dv(mat)
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: mat
    REAL(DP), DIMENSION(size(mat,1)) :: get_diag_dv
    INTEGER(I4) :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_dv')
    do j=1,size(mat,1)
       get_diag_dv(j)=mat(j,j)
    end do
  END FUNCTION get_diag_dv
  !BL
  SUBROUTINE put_diag_rv(diagv,mat)
    REAL(SP), DIMENSION(:), INTENT(IN) :: diagv
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER(I4) :: j,n
    n=assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
    do j=1,n
       mat(j,j)=diagv(j)
    end do
  END SUBROUTINE put_diag_rv
  !BL
  SUBROUTINE put_diag_r(scal,mat)
    REAL(SP), INTENT(IN) :: scal
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER(I4) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=scal
    end do
  END SUBROUTINE put_diag_r
  !BL
  SUBROUTINE unit_matrix_sp(mat)
    REAL(SP), DIMENSION(:,:), INTENT(OUT) :: mat
    INTEGER(I4) :: i,n
    n=min(size(mat,1),size(mat,2))
    mat(:,:)=0.0_sp
    do i=1,n
       mat(i,i)=1.0_sp
    end do
  END SUBROUTINE unit_matrix_sp
  SUBROUTINE unit_matrix_dp(mat)
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: mat
    INTEGER(I4) :: i,n
    n=min(size(mat,1),size(mat,2))
    mat(:,:)=0.0_dp
    do i=1,n
       mat(i,i)=1.0_dp
    end do
  END SUBROUTINE unit_matrix_dp
  !BL
  FUNCTION upper_triangle(j,k,extra)
    INTEGER(I4), INTENT(IN) :: j,k
    INTEGER(I4), OPTIONAL, INTENT(IN) :: extra
    LOGICAL, DIMENSION(j,k) :: upper_triangle
    INTEGER(I4) :: n
    n=0
    if (present(extra)) n=extra
    upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
  END FUNCTION upper_triangle
  !BL
  FUNCTION lower_triangle(j,k,extra)
    INTEGER(I4), INTENT(IN) :: j,k
    INTEGER(I4), OPTIONAL, INTENT(IN) :: extra
    LOGICAL, DIMENSION(j,k) :: lower_triangle
    INTEGER(I4) :: n
    n=0
    if (present(extra)) n=extra
    lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
  END FUNCTION lower_triangle
  !BL
  FUNCTION vabs_sp(v)
    REAL(SP), DIMENSION(:), INTENT(IN) :: v
    REAL(SP) :: vabs_sp
    vabs_sp=sqrt(dot_product(v,v))
  END FUNCTION vabs_sp
  FUNCTION vabs_dp(v)
    REAL(DP), DIMENSION(:), INTENT(IN) :: v
    REAL(DP) :: vabs_dp
    vabs_dp=sqrt(dot_product(v,v))
  END FUNCTION vabs_dp
  !BL
END MODULE mo_nrutil

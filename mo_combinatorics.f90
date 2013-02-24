!> \file mo_combinatorics.f90

!> \brief Package for combinatorial calculations.

!> \details This package provides routines and functions for combinatorial calculations.

!> \authors Matthias Cuntz, Giovanni Dalmasso, Juliane Mai, Stephan Thober
!> \date Feb 2013

MODULE mo_combinatorics

  ! Written  Matthias Cuntz, Giovanni Dalmasso, Juliane Mai, Stephan Thober Feb 2013
  ! Modified 

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
  ! along with the UFZ Fortran library. If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2011-2013 Matthias Cuntz, Giovanni Dalmasso, Juliane Mai, Stephan Thober


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

  USE mo_kind, ONLY: i4, i8, sp, dp

  IMPLICIT NONE

  PUBLIC :: binomcoeffi     ! Binomial coefficient (n choose k)
  PUBLIC :: random_kofn    ! Random selection of k of n
  PUBLIC :: next_kofn      ! Next selection of k of n to a given one
  PUBLIC :: all_kofn       ! All selections of k of n
  ! PUBLIC :: random_permut  ! Random permutation to a given one
  ! PUBLIC :: next_permut    ! Next permutation to a given one
  ! PUBLIC :: all_permut     ! All permutations of n

  ! ------------------------------------------------------------------

  !     NAME
  !         binomcoeffi

  !     PURPOSE
  !
  !>        \brief The binomial coefficient.
  !
  !>        \details Calculates the binomial coefficient (n choose k):
  !>                     \f[ C(n,k) = \frac{n!}{k! (n-k)!} \f], 
  !>                 i.e. the number of possibilities to select from n numbers k
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8) :: n/n(:)"       from n numbers ...
  !>        \param[in] "integer(i4/i8) :: k/k(:)"       ... k will be selected
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     integer(i4/i8) :: binomcoeffi/binomcoeffi(:) &mdash; Binomial coefficient (n choose k)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         bico = binomcoeffi(5,2)
  !         bico --> 10
  !         -> see also example in test directory

  !     LITERATURE
  !         

  !     HISTORY
  !>        \author Matthias Cuntz, Juliane Mai
  !>        \date Feb 2013
  !         Modified, 

  INTERFACE binomcoeffi
     MODULE PROCEDURE binomcoeffi_i4_d0, binomcoeffi_i8_d0, binomcoeffi_i4_d1, binomcoeffi_i8_d1
  END INTERFACE binomcoeffi

  ! ------------------------------------------------------------------

  !     NAME
  !         all_kofn

  !     PURPOSE
  !
  !>        \brief All possible selections of k numbers from n.
  !
  !>        \details Returns all possible selections of a k-subset from n numbers.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8) :: n"                 select from n numbers ...
  !>        \param[in] "integer(i4/i8) :: k"                 ... k numbers
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     integer(i4/i8), allocatable, dimension(:,:) :: all &mdash; All selections. Dim_2=k.
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         all = all_kofn(3, 2)
  !         all --> (/ (/ 1,2 /), (/ 1,3 /), (/ 2,3 /) /) 
  !         -> see also example in test directory

  !     LITERATURE
  !         Nijenhuis, A., & Wilf, H. S. (1978). 
  !         Combinatorial algorithms for computers and calculators (2nd ed.). 
  !         Academic Pr.

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Feb 2013
  !         Modified, 

  INTERFACE all_kofn
     MODULE PROCEDURE all_kofn_i4, all_kofn_i8
  END INTERFACE all_kofn

  ! ------------------------------------------------------------------

  !     NAME
  !         next_kofn

  !     PURPOSE
  !
  !>        \brief The next selection of k numbers from n to a given one.
  !
  !>        \details Determines the next selection of k numbers from k to a given one./n
  !>                 If one has n=5 numbers and want to pick k=3, the first possible selection is (1,2,3).
  !>                 The next selection will be (1,2,4). The subsequent selection to (1,3,4) will be (1,3,5)./n
  !>                 The given subset (previous) has to be sorted.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8) :: n"                 select from n numbers ...
  !>        \param[in] "integer(i4/i8) :: k"                 ... k numbers
  !>        \param[in] "integer(i4/i8) :: previous(:)"       previous selection (dim_1 = k)
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     integer(i4/i8) :: next(size(previous)) &mdash; Next selection of k from n numbers
  !
  !     RESTRICTIONS
  !>        \note The given subset (previous) has to be sorted.
  !
  !     EXAMPLE
  !         next = next_kofn(5, 3, (/1,3,4/))
  !         next --> (/ 1,3,5 /)
  !         -> see also example in test directory

  !     LITERATURE
  !         Nijenhuis, A., & Wilf, H. S. (1978). 
  !         Combinatorial algorithms for computers and calculators (2nd ed.). 
  !         Academic Pr.

  !     HISTORY
  !>        \author Giovanni Dalmasso, Juliane Mai
  !>        \date Feb 2013
  !         Modified, 

  INTERFACE next_kofn
     MODULE PROCEDURE next_kofn_i4, next_kofn_i8
  END INTERFACE next_kofn

  ! ------------------------------------------------------------------

  !     NAME
  !         random_kofn

  !     PURPOSE
  !
  !>        \brief A random selection of a k-subset of n numbers (1..n).
  !
  !>        \details Returns a random k-subset of n numbers (1..n)/n
  !>                 The returned subset will be sorted./n
  !>                 The code is adapted from Nijenhuis (1978) - routine RANKSB.
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4/i8) :: n"                 select from n numbers ...
  !>        \param[in] "integer(i4/i8) :: k"                 ... k numbers
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !>        \param[in,out] "integer(i4/i8), optional :: save_state(n_save_state)"  an array for saving the state 
  !>                                                                               of an uniform random number stream
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !>        \return     integer(i4/i8) :: set(k) &mdash; Random selection of k from n numbers
  !
  !     RESTRICTIONS
  !>        \note A random number stream for generating uniform random numbers has to be initialized before running this function. 
  !
  !     EXAMPLE
  !         set = random_kofn(5, 3)
  !         set --> (/ 2,3,5 /)
  !         -> see also example in test directory

  !     LITERATURE
  !         Nijenhuis, A., & Wilf, H. S. (1978). 
  !         Combinatorial algorithms for computers and calculators (2nd ed.). 
  !         Academic Pr.

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Feb 2013
  !         Modified, 

  INTERFACE random_kofn
     MODULE PROCEDURE random_kofn_i4 , random_kofn_i8
  END INTERFACE random_kofn

  ! ------------------------------------------------------------------

  PRIVATE

  INTERFACE factln
     MODULE PROCEDURE factln_i4_s_dp, factln_i4_v_dp, factln_i8_s_dp, factln_i8_v_dp
  END INTERFACE factln

  INTERFACE gammln
     MODULE PROCEDURE gammln_s_sp, gammln_s_dp, gammln_v_sp, gammln_v_dp
  END INTERFACE gammln

  INTERFACE arth
     MODULE PROCEDURE arth_sp_i4, arth_dp_i4, arth_sp_i8, arth_dp_i8, arth_i4, arth_i8
  END INTERFACE arth

  INTEGER(i4), PARAMETER :: NPAR_ARTH=16, NPAR2_ARTH=8

  ! ------------------------------------------------------------------

CONTAINS

  FUNCTION binomcoeffi_i4_d0(n,k)
    ! Returns the binomial coefficient as an integer number.
    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: n, k
    INTEGER(i4)             :: binomcoeffi_i4_d0

    binomcoeffi_i4_d0 = nint(exp(factln(n)-factln(k)-factln(n-k)))

  END FUNCTION binomcoeffi_i4_d0

  FUNCTION binomcoeffi_i8_d0(n,k)
    ! Returns the binomial coefficient as an integer number.
    IMPLICIT NONE

    INTEGER(i8), INTENT(IN) :: n, k
    INTEGER(i8)             :: binomcoeffi_i8_d0

    binomcoeffi_i8_d0 = nint(exp(factln(n)-factln(k)-factln(n-k)))

  END FUNCTION binomcoeffi_i8_d0

  FUNCTION binomcoeffi_i4_d1(n,k)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:), INTENT(IN) :: n, k
    INTEGER(i4), DIMENSION(size(n))       :: binomcoeffi_i4_d1

    if (size(n) /= size(k)) stop 'Error binomcoeffi: size(n) /= size(k)'

    binomcoeffi_i4_d1 = nint(exp(factln(n)-factln(k)-factln(n-k)))

  END FUNCTION binomcoeffi_i4_d1

  FUNCTION binomcoeffi_i8_d1(n,k)

    IMPLICIT NONE

    INTEGER(i8), DIMENSION(:), INTENT(IN) :: n, k
    INTEGER(i8), DIMENSION(size(n))       :: binomcoeffi_i8_d1

    if (size(n) /= size(k)) stop 'Error binomcoeffi: size(n) /= size(k)'

    binomcoeffi_i8_d1 = nint(exp(factln(n)-factln(k)-factln(n-k)))

  END FUNCTION binomcoeffi_i8_d1

  FUNCTION all_kofn_i4(n, k) result(all)

    IMPLICIT NONE

    INTEGER(I4),                         INTENT(IN)  :: n             ! from n numbers will be selected
    INTEGER(I4),                         INTENT(IN)  :: k             ! k numbers will be selected
    INTEGER(I4), DIMENSION(:,:), allocatable         :: all           ! all subsets

    ! local variables
    integer(i4)         :: bico
    integer(i4)         :: i

    bico = binomcoeffi(n,k)
    allocate(all(bico,k))

    do i=1,k
       all(1,i) = i
    end do

    do i=2,bico
       all(i,:) = next_kofn(n, k,  all(i-1,:))
    end do

  END FUNCTION all_kofn_i4

  FUNCTION all_kofn_i8(n, k) result(all)

    IMPLICIT NONE

    INTEGER(I8),                         INTENT(IN)  :: k             ! k numbers will be selected
    INTEGER(I8),                         INTENT(IN)  :: n             ! from n numbers will be selected
    INTEGER(I8), DIMENSION(:,:), allocatable         :: all           ! all subsets

    ! local variables
    integer(i8)         :: bico
    integer(i8)         :: i

    bico = binomcoeffi(n,k)
    allocate(all(bico,k))

    do i=1,k
       all(1,i) = i
    end do

    do i=2,bico
       all(i,:) = next_kofn(n, k,  all(i-1,:))
    end do

  END FUNCTION all_kofn_i8

  FUNCTION next_kofn_i4(n, k, previous) result(next)

    IMPLICIT NONE

    INTEGER(I4),                         INTENT(IN)  :: k             ! k numbers will be selected
    INTEGER(I4),                         INTENT(IN)  :: n             ! from n numbers will be selected
    INTEGER(I4), DIMENSION(k),           INTENT(IN)  :: previous      ! previous subset
    INTEGER(I4), DIMENSION(k)                        :: next          ! next subset

    ! local variables
    integer(i4)         :: indx
    logical             :: found
    integer(i4)         :: i

    ! find index from the back to start changing the values
    indx = k
    found = .false.
    do while ( (.not. found) .and. (indx .gt. 1) ) 
       if ( previous(indx) .eq. n-k+indx ) then
          indx = indx-1
       else
          found = .true.
       end if
    end do
    if( (indx .gt. 1) .or. ( previous(1) .ne. (n-k+1)) ) then
       ! it is a subset in between
       next = previous
       next(indx) = previous(indx) + 1
       do i=indx+1,k
          next(i) = next(i-1) + 1
       end do
    else
       ! there do not exist a next subset and first one (/1, ..., k/) is returned
       do i=1, k
          next(i) = i
       end do
    end if        

  END FUNCTION next_kofn_i4

  FUNCTION next_kofn_i8(n, k, previous) result(next)

    IMPLICIT NONE

    INTEGER(I8),                         INTENT(IN)  :: k             ! k numbers will be selected
    INTEGER(I8),                         INTENT(IN)  :: n             ! from n numbers will be selected
    INTEGER(I8), DIMENSION(k),           INTENT(IN)  :: previous      ! previous subset
    INTEGER(I8), DIMENSION(k)                        :: next          ! next subset

    ! local variables
    integer(i8)         :: indx
    logical             :: found
    integer(i8)         :: i

    ! find index from the back to start changing the values
    indx = k
    found = .false.
    do while ( (.not. found) .and. (indx .gt. 1) ) 
       if ( previous(indx) .eq. n-k+indx ) then
          indx = indx-1
       else
          found = .true.
       end if
    end do

     if( (indx .gt. 1) .or. ( previous(1) .ne. (n-k+1)) ) then
       ! it is a subset in between
       next = previous
       next(indx) = previous(indx) + 1
       do i=indx+1,k
          next(i) = next(i-1) + 1
       end do
    else
       ! there do not exist a next subset and first one (/1, ..., k/) is returned
       do i=1, k
          next(i) = i
       end do
    end if        

  END FUNCTION next_kofn_i8

  FUNCTION random_kofn_i4(n, k, save_state) result(set)
    
    use mo_xor4096, only: xor4096, n_save_state

    implicit none

    INTEGER(I4),                                    INTENT(IN)    :: k             ! k numbers will be selected
    INTEGER(I4),                                    INTENT(IN)    :: n             ! from n numbers will be selected
    INTEGER(I4), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state    ! for saving state of a uniform, sp 
    !                                                                              ! random number stream
    INTEGER(I4), DIMENSION(k)                                     :: set           ! random subset

    ! local variables
    logical     :: gout
    integer(i4) :: x, r, ds, p, s, c, i, l, m,m0, iseed
    real(sp)    :: rn

    iseed = 0_i4

    ! (A)
    c=k
    do  i=1,k
       set(i)=(i-1)*n/k
    end do
    
    ! (B)
    do while (c .gt. 0)
       if (present(save_state)) then
          call xor4096(iSeed, rn, save_state=save_state)
       else
          call xor4096(iSeed, rn)
       end if
       x = 1 + int( real(n,sp)*rn, i4 )
       l = 1 + int( real(x*k-1,sp)/real(n,sp), i4 )
       if(x.le.set(l)) cycle
       set(l) = set(l) + 1
       c = c - 1
    end do
    i = 0
    p = 0
    s = k

    ! (C)
    gout = .false.
    do while (.not. gout)
       i = i+1
       if (i .gt. k) then
          gout = .true.
       else
          if (set(i) .eq. int(real((i-1)*n,sp)/real(k),i4)) then
             set(i) = 0
          else
             p=p+1
             m=set(i)
             set(i) = 0
             set(p) = m
          end if
       end if
    end do

    ! (D)
    gout = .false.
    do while (.not. gout)
       l  = 1+int(real(set(p)*k-1,sp)/real(n,sp),i4)
       ds = set(p) - int( real((l-1)*n,sp)/real(k,sp) ,i4)
       set(p) = 0
       set(s) = l
       s=s-ds
       p=p-1
       if (p .le. 0) gout = .true.
    end do
    l = k
     
    do while( l .ne. 0)
       
       ! (E)
       if (set(l) .ne. 0) then
          r = l
          m0 = 1+int(real((set(l)-1)*n,sp)/real(k,sp),i4)
          m = int(real(set(l)* n, sp) / real(k,sp),i4) - m0 + 1
       end if
       
       ! (F)
       if (present(save_state)) then
          call xor4096(iSeed, rn, save_state=save_state)
       else
          call xor4096(iSeed, rn)
       end if
       x = m0 + int(real(m,sp) * rn,i4)
       i = l
       
       ! (G)
       gout = .false.
       do while (.not. gout)
          i = i + 1
          if (i .le. r) then
             if( x .ge. set(i)) then
                set(i-1) = set(i)
                x = x+1
             else 
                gout = .true.
             end if
          else
             gout = .true.
          end if
       end do
       
       ! (H)
       set(i-1) = x
       m = m-1
       l=l-1

    end do

  END FUNCTION random_kofn_i4

  FUNCTION random_kofn_i8(n, k, save_state) result(set)
    
    use mo_xor4096, only: xor4096, n_save_state

    implicit none

    INTEGER(I8),                                    INTENT(IN)    :: k             ! k numbers will be selected
    INTEGER(I8),                                    INTENT(IN)    :: n             ! from n numbers will be selected
    INTEGER(I8), DIMENSION(n_save_state), OPTIONAL, INTENT(INOUT) :: save_state    ! for saving state of a uniform, dp 
    !                                                                              ! random number stream
    INTEGER(I8), DIMENSION(k)                                     :: set           ! random subset

    ! local variables
    logical     :: gout
    integer(i8) :: x, r, ds, p, s, c, i, l, m,m0, iseed
    real(dp)    :: rn

    iseed = 0_i8

    ! (A)
    c=k
    do  i=1,k
       set(i)=(i-1)*n/k
    end do
    
    ! (B)
    do while (c .gt. 0)
       if (present(save_state)) then
          call xor4096(iSeed, rn, save_state=save_state)
       else
          call xor4096(iSeed, rn)
       end if
       x = 1 + int( real(n,dp)*rn, i8 )
       l = 1 + int( real(x*k-1,dp)/real(n,dp), i8 )
       if(x.le.set(l)) cycle
       set(l) = set(l) + 1
       c = c - 1
    end do
    i = 0
    p = 0
    s = k

    ! (C)
    gout = .false.
    do while (.not. gout)
       i = i+1
       if (i .gt. k) then
          gout = .true.
       else
          if (set(i) .eq. int(real((i-1)*n,dp)/real(k),i8)) then
             set(i) = 0
          else
             p=p+1
             m=set(i)
             set(i) = 0
             set(p) = m
          end if
       end if
    end do

    ! (D)
    gout = .false.
    do while (.not. gout)
       l  = 1+int(real(set(p)*k-1,dp)/real(n,dp),i8)
       ds = set(p) - int( real((l-1)*n,dp)/real(k,dp) ,i8)
       set(p) = 0
       set(s) = l
       s=s-ds
       p=p-1
       if (p .le. 0) gout = .true.
    end do
    l = k
     
    do while( l .ne. 0)
       
       ! (E)
       if (set(l) .ne. 0) then
          r = l
          m0 = 1+int(real((set(l)-1)*n,dp)/real(k,dp),i8)
          m = int(real(set(l)* n, dp) / real(k,dp),i8) - m0 + 1
       end if
       
       ! (F)
       if (present(save_state)) then
          call xor4096(iSeed, rn, save_state=save_state)
       else
          call xor4096(iSeed, rn)
       end if
       x = m0 + int(real(m,dp) * rn,i8)
       i = l
       
       ! (G)
       gout = .false.
       do while (.not. gout)
          i = i + 1
          if (i .le. r) then
             if( x .ge. set(i)) then
                set(i-1) = set(i)
                x = x+1
             else 
                gout = .true.
             end if
          else
             gout = .true.
          end if
       end do
       
       ! (H)
       set(i-1) = x
       m = m-1
       l=l-1

    end do

  END FUNCTION random_kofn_i8


  ! ------------------------------------------------------------------
  ! PRIVATE
  ! ------------------------------------------------------------------

  FUNCTION factln_i4_s_dp(n)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: n
    REAL(dp)                :: factln_i4_s_dp

    INTEGER(i4), PARAMETER :: TMAX=100
    REAL(dp), DIMENSION(TMAX), SAVE :: a
    LOGICAL, SAVE :: init=.true.

    !$omp   threadprivate(a, init)

    if (init) then
       a(1:TMAX) = gammln(arth(1.0_dp,1.0_dp,TMAX))
       init=.false.
    end if
    if (n < 0) stop 'Error factln_i4_s_dp: n < 0'
    if (n < TMAX) then
       factln_i4_s_dp = a(n+1)
    else
       factln_i4_s_dp = gammln(real(n,dp)+1.0_dp)
    end if
  END FUNCTION factln_i4_s_dp

  FUNCTION factln_i4_v_dp(n)

    IMPLICIT NONE

    INTEGER(i4), DIMENSION(:),       INTENT(IN) :: n
    REAL(dp),    DIMENSION(size(n))             :: factln_i4_v_dp

    LOGICAL,     DIMENSION(size(n))             :: mask
    INTEGER(i4), PARAMETER                      :: TMAX=100
    REAL(dp),    DIMENSION(TMAX), SAVE          :: a
    LOGICAL, SAVE                               :: init=.true.

    !$omp   threadprivate(a, init)

    if (init) then
       a(1:TMAX) = gammln(arth(1.0_dp,1.0_dp,TMAX))
       init = .false.
    end if
    if (any(n < 0)) stop 'Error factln_i4_v_dp: some n < 0'
    mask = (n >= TMAX)
    factln_i4_v_dp = unpack(gammln(real(pack(n,mask),dp)+1.0_dp),mask,0.0_dp)
    where (.not. mask) factln_i4_v_dp = a(n+1)

  END FUNCTION factln_i4_v_dp

  FUNCTION factln_i8_s_dp(n)

    IMPLICIT NONE

    INTEGER(i8), INTENT(IN) :: n
    REAL(dp)                :: factln_i8_s_dp

    INTEGER(i8), PARAMETER :: TMAX=100
    REAL(dp), DIMENSION(TMAX), SAVE :: a
    LOGICAL, SAVE :: init=.true.

    !$omp   threadprivate(a, init)

    if (init) then
       a(1:TMAX) = gammln(arth(1.0_dp,1.0_dp,TMAX))
       init=.false.
    end if
    if (n < 0) stop 'Error factln_i8_s_dp: n < 0'
    if (n < TMAX) then
       factln_i8_s_dp = a(n+1)
    else
       factln_i8_s_dp = gammln(real(n,dp)+1.0_dp)
    end if
  END FUNCTION factln_i8_s_dp

  FUNCTION factln_i8_v_dp(n)

    IMPLICIT NONE

    INTEGER(i8), DIMENSION(:),       INTENT(IN) :: n
    REAL(dp),    DIMENSION(size(n))             :: factln_i8_v_dp

    LOGICAL,     DIMENSION(size(n))             :: mask
    INTEGER(i8), PARAMETER                      :: TMAX=100
    REAL(dp),    DIMENSION(TMAX), SAVE          :: a
    LOGICAL, SAVE                               :: init=.true.

    !$omp   threadprivate(a, init)

    if (init) then
       a(1:TMAX) = gammln(arth(1.0_dp,1.0_dp,TMAX))
       init = .false.
    end if
    if (any(n < 0)) stop 'Error factln_i8_v_dp: some n < 0'
    mask = (n >= TMAX)
    factln_i8_v_dp = unpack(gammln(real(pack(n,mask),dp)+1.0_dp),mask,0.0_dp)
    where (.not. mask) factln_i8_v_dp = a(n+1)

  END FUNCTION factln_i8_v_dp

  FUNCTION gammln_s_sp(z)
    !  Uses Lanczos-type approximation to ln(gamma) for z > 0.
    !  Reference:
    !       Lanczos, C. ''A precision approximation of the gamma
    !               function'', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
    !  Accuracy: About 14 significant digits except for small regions
    !            in the vicinity of 1 and 2.
    !  Programmer: Alan Miller
    !              1 Creswick Street, Brighton, Vic. 3187, Australia
    !  e-mail: amiller @ bigpond.net.au
    !  Latest revision - 14 October 1996

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: z
    REAL(sp)             :: gammln_s_sp

    ! Local variables
    REAL(sp)  :: a(9) = (/ 0.9999999999995183_sp, 676.5203681218835_sp, &
         -1259.139216722289_sp, 771.3234287757674_sp, &
         -176.6150291498386_sp, 12.50734324009056_sp, &
         -0.1385710331296526_sp, 0.9934937113930748E-05_sp, &
         0.1659470187408462E-06_sp /)
    REAL(sp)  :: zero = 0.0_sp,   &
         one = 1.0_sp, &
         lnsqrt2pi =  0.9189385332046727_sp, &
         half = 0.5_sp, &
         sixpt5 = 6.5_sp, &
         seven = 7.0_sp
    REAL(sp)  :: tmp
    INTEGER   :: j

    if (z <= 0.0_sp) stop 'Error gammln_s_sp: z <= 0'
    gammln_s_sp = zero
    tmp = z + seven
    DO j = 9, 2, -1
       gammln_s_sp = gammln_s_sp + a(j)/tmp
       tmp = tmp - one
    END DO
    gammln_s_sp = gammln_s_sp + a(1)
    gammln_s_sp = LOG(gammln_s_sp) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)

  END FUNCTION gammln_s_sp

  FUNCTION gammln_v_sp(z)
    !  Uses Lanczos-type approximation to ln(gamma) for z > 0.
    !  Reference:
    !       Lanczos, C. ''A precision approximation of the gamma
    !               function'', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
    !  Accuracy: About 14 significant digits except for small regions
    !            in the vicinity of 1 and 2.
    !  Programmer: Alan Miller
    !              1 Creswick Street, Brighton, Vic. 3187, Australia
    !  e-mail: amiller @ bigpond.net.au
    !  Latest revision - 14 October 1996

    IMPLICIT NONE

    REAL(sp), DIMENSION(:), INTENT(IN) :: z
    REAL(sp), DIMENSION(size(z))       :: gammln_v_sp

    ! Local variables
    REAL(sp)  :: a(9) = (/ 0.9999999999995183_sp, 676.5203681218835_sp, &
         -1259.139216722289_sp, 771.3234287757674_sp, &
         -176.6150291498386_sp, 12.50734324009056_sp, &
         -0.1385710331296526_sp, 0.9934937113930748E-05_sp, &
         0.1659470187408462E-06_sp /)
    REAL(sp)  :: zero = 0.0_sp,   &
         one = 1.0_sp, &
         lnsqrt2pi =  0.9189385332046727_sp, &
         half = 0.5_sp, &
         sixpt5 = 6.5_sp, &
         seven = 7.0_sp
    REAL(sp), DIMENSION(size(z)) :: tmp
    INTEGER   :: j

    if (any(z <= 0.0_sp)) stop 'Error gammln_v_sp: some z <= 0'
    gammln_v_sp = zero
    tmp = z + seven
    DO j = 9, 2, -1
       gammln_v_sp = gammln_v_sp + a(j)/tmp
       tmp = tmp - one
    END DO
    gammln_v_sp = gammln_v_sp  + a(1)
    gammln_v_sp = LOG(gammln_v_sp) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)

  END FUNCTION gammln_v_sp

  FUNCTION gammln_s_dp(z)
    !  Uses Lanczos-type approximation to ln(gamma) for z > 0.
    !  Reference:
    !       Lanczos, C. ''A precision approximation of the gamma
    !               function'', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
    !  Accuracy: About 14 significant digits except for small regions
    !            in the vicinity of 1 and 2.
    !  Programmer: Alan Miller
    !              1 Creswick Street, Brighton, Vic. 3187, Australia
    !  e-mail: amiller @ bigpond.net.au
    !  Latest revision - 14 October 1996

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: z
    REAL(dp)             :: gammln_s_dp

    ! Local variables
    REAL(dp)  :: a(9) = (/ 0.9999999999995183_dp, 676.5203681218835_dp, &
         -1259.139216722289_dp, 771.3234287757674_dp, &
         -176.6150291498386_dp, 12.50734324009056_dp, &
         -0.1385710331296526_dp, 0.9934937113930748E-05_dp, &
         0.1659470187408462E-06_dp /)
    REAL(dp)  :: zero = 0.0_dp,   &
         one = 1.0_dp, &
         lnsqrt2pi =  0.9189385332046727_dp, &
         half = 0.5_dp, &
         sixpt5 = 6.5_dp, &
         seven = 7.0_dp
    REAL(dp)  :: tmp
    INTEGER   :: j

    if (z <= 0.0_dp) stop 'Error gammln_s_dp: z <= 0'
    gammln_s_dp = zero
    tmp = z + seven
    DO j = 9, 2, -1
       gammln_s_dp = gammln_s_dp + a(j)/tmp
       tmp = tmp - one
    END DO
    gammln_s_dp = gammln_s_dp + a(1)
    gammln_s_dp = LOG(gammln_s_dp) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)

  END FUNCTION gammln_s_dp

  FUNCTION gammln_v_dp(z)
    !  Uses Lanczos-type approximation to ln(gamma) for z > 0.
    !  Reference:
    !       Lanczos, C. ''A precision approximation of the gamma
    !               function'', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
    !  Accuracy: About 14 significant digits except for small regions
    !            in the vicinity of 1 and 2.
    !  Programmer: Alan Miller
    !              1 Creswick Street, Brighton, Vic. 3187, Australia
    !  e-mail: amiller @ bigpond.net.au
    !  Latest revision - 14 October 1996

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN) :: z
    REAL(dp), DIMENSION(size(z))       :: gammln_v_dp

    ! Local variables
    REAL(dp)  :: a(9) = (/ 0.9999999999995183_dp, 676.5203681218835_dp, &
         -1259.139216722289_dp, 771.3234287757674_dp, &
         -176.6150291498386_dp, 12.50734324009056_dp, &
         -0.1385710331296526_dp, 0.9934937113930748E-05_dp, &
         0.1659470187408462E-06_dp /)
    REAL(dp)  :: zero = 0.0_dp,   &
         one = 1.0_dp, &
         lnsqrt2pi =  0.9189385332046727_dp, &
         half = 0.5_dp, &
         sixpt5 = 6.5_dp, &
         seven = 7.0_dp
    REAL(dp), DIMENSION(size(z)) :: tmp
    INTEGER   :: j

    if (any(z <= 0.0_dp)) stop 'Error gammln_v_dp: some z <= 0'
    gammln_v_dp = zero
    tmp = z + seven
    DO j = 9, 2, -1
       gammln_v_dp = gammln_v_dp + a(j)/tmp
       tmp = tmp - one
    END DO
    gammln_v_dp = gammln_v_dp  + a(1)
    gammln_v_dp = LOG(gammln_v_dp) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)

  END FUNCTION gammln_v_dp

  ! Return an arithmetic progression as an array.
  FUNCTION arth_sp_i4(first,increment,n)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: first,increment
    INTEGER(i4), INTENT(IN) :: n
    REAL(sp), DIMENSION(n) :: arth_sp_i4

    INTEGER(i4) :: k,k2
    REAL(sp) :: temp

    if (n > 0) arth_sp_i4(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_sp_i4(k)=arth_sp_i4(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_sp_i4(k)=arth_sp_i4(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_sp_i4(k+1:min(k2,n))=temp+arth_sp_i4(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if

  END FUNCTION arth_sp_i4

  ! Return an arithmetic progression as an array.
  FUNCTION arth_dp_i4(first,increment,n)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: first,increment
    INTEGER(i4), INTENT(IN) :: n
    REAL(dp), DIMENSION(n) :: arth_dp_i4

    INTEGER(i4) :: k,k2
    REAL(dp) :: temp

    if (n > 0) arth_dp_i4(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_dp_i4(k)=arth_dp_i4(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_dp_i4(k)=arth_dp_i4(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_dp_i4(k+1:min(k2,n))=temp+arth_dp_i4(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if

  END FUNCTION arth_dp_i4

  ! Return an arithmetic progression as an array.
  FUNCTION arth_sp_i8(first,increment,n)

    IMPLICIT NONE

    REAL(sp), INTENT(IN) :: first,increment
    INTEGER(i8), INTENT(IN) :: n
    REAL(sp), DIMENSION(n) :: arth_sp_i8

    INTEGER(i8) :: k,k2
    REAL(sp) :: temp

    if (n > 0) arth_sp_i8(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_sp_i8(k)=arth_sp_i8(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_sp_i8(k)=arth_sp_i8(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_sp_i8(k+1:min(k2,n))=temp+arth_sp_i8(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if

  END FUNCTION arth_sp_i8

  ! Return an arithmetic progression as an array.
  FUNCTION arth_dp_i8(first,increment,n)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: first,increment
    INTEGER(i8), INTENT(IN) :: n
    REAL(dp), DIMENSION(n) :: arth_dp_i8

    INTEGER(i8) :: k,k2
    REAL(dp) :: temp

    if (n > 0) arth_dp_i8(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_dp_i8(k)=arth_dp_i8(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_dp_i8(k)=arth_dp_i8(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_dp_i8(k+1:min(k2,n))=temp+arth_dp_i8(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if

  END FUNCTION arth_dp_i8

  FUNCTION arth_i4(first,increment,n)

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: first,increment,n
    INTEGER(i4), DIMENSION(n) :: arth_i4

    INTEGER(i4) :: k,k2,temp

    if (n > 0) arth_i4(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i4(k)=arth_i4(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i4(k)=arth_i4(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i4(k+1:min(k2,n))=temp+arth_i4(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if

  END FUNCTION arth_i4

  FUNCTION arth_i8(first,increment,n)

    IMPLICIT NONE

    INTEGER(i8), INTENT(IN) :: first,increment,n
    INTEGER(i8), DIMENSION(n) :: arth_i8

    INTEGER(i8) :: k,k2,temp

    if (n > 0) arth_i8(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i8(k)=arth_i8(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i8(k)=arth_i8(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i8(k+1:min(k2,n))=temp+arth_i8(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if

  END FUNCTION arth_i8

END MODULE mo_combinatorics

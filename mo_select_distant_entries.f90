!> \file mo_select_distant_entries.f90

!> \brief Selects elements of distance matrix with most largest distance to each other

!> \details Selects from a distance matrix (symmetric and positive) a specified 
!>          number (N) of elements which have large distances to each other.

!> \authors Juliane Mai
!> \date November 2014

module mo_select_distant_entries

  ! Written November 2014, Juliane Mai

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
  ! along with the JAMS makefile project (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2014 Juliane Mai

  use mo_kind,         only: i4, sp, dp
  use mo_sort,         only: sort
  use mo_xor4096,      only: get_timeseed, n_save_state, xor4096
  use mo_xor4096_apps, only: xor4096_range

  implicit none

  public :: select_distant_entries   ! selects elements of distance matrix with most largest distance to each other

  ! ------------------------------------------------------------------

  !     NAME
  !         select_distant_entries

  !     PURPOSE
  !         Selects from a distance matrix (symmetric and positive) a specified 
  !         number (N) of elements which have large distances to each other.
  !
  !         Indexes are chosen iteratively. Therefore, it might be possible to 
  !         find selections which have an overall larger distance. This is only
  !         a good GUESS of such a selection.
  !
  !         There are two ways implemented how the first pair is chosen.
  !         (1) Default
  !             The first pair is the one with the overall largest distance,
  !             i.e. largest entry in distance matrix
  !             Therfore, the selection is deterministic.
  !         (2) The first pair is randomly chosen.
  !             Therefore, the returened selection is random.

  !     CALLING SEQUENCE
  !         indexes = select_distant_entries( distance_matrix, N, first=0 )

  !     INTENT(IN)
  !         real(sp/dp) :: distance_matrix   2D array of distances between element i and j
  !         integer(i4) :: n_elements        number of elements to return 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         integer(i4) :: first             method how first pair is chosen
  !                                          0 : first pair is the one with the overall 
  !                                              largest distance,
  !                                              i.e. result is deterministic
  !                                          1 : first pair is randomly chosen,
  !                                              i.e. result is random
  !                                          DEFAULT: 0 
  !         integer(i4) :: seed              seed used for random number generator
  !                                          (only considered if first=1)
  !                                          DEFAULT: timeseed

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         integer(i4) :: select_distant_entries(n_elements)   indexes of elements with largest distance to each other

  !     RESTRICTIONS
  !         There is no initial mask of distance matrix allowed.

  !     EXAMPLE
  !         -> see also example in test directory
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !         Written,  Juliane Mai, Nov 2014

  interface select_distant_entries
     module procedure select_distant_entries_dp, select_distant_entries_sp
  end interface select_distant_entries

  private

contains

  function select_distant_entries_dp( distance_matrix, n_elements, first, seed )
    
    implicit none

    real(dp), dimension(:,:), intent(in) :: distance_matrix           ! square matrix with distance between element i and j 
    integer(i4),              intent(in) :: n_elements                ! number of elements to return 
    integer(i4),    optional, intent(in) :: first                     ! method how first pair is chosen
    !                                                                 !    0 : first pair is the one with the overall 
    !                                                                 !        largest distance,
    !                                                                 !        i.e. result is deterministic
    !                                                                 !    1 : first pair is randomly chosen,
    !                                                                 !        i.e. result is random
    !                                                                 ! DEFAULT: 0 
    integer(i4),    optional, intent(in) :: seed                      ! seed used for random number generator
    !                                                                 ! (only considered if first=1)
    !                                                                 ! DEFAULT: timeseed
    integer(i4), dimension(n_elements)   :: select_distant_entries_dp ! indexes of elements with largest distance to each other

    ! local variables
    integer(i4)                             :: ii, jj, kk
    integer(i4)                             :: ifirst      ! method how first pair is chosen
    integer(i4)                             :: iseed       ! seed for random number generator
    integer(i4), dimension(n_save_state)    :: save_state
    integer(i4)                             :: n_total     ! number of elements --> number of rows/cols of distance_matrix
    integer(i4)                             :: maxpos(2)
    logical,    dimension(:,:), allocatable :: mask

    ! handling of optionals
    if ( present(first) ) then
       ifirst = first
    else
       ifirst = 0
    end if

    ! handling of optionals
    if ( present(seed) ) then
       iseed = seed
    else
       call get_timeseed(iseed)
    end if

    n_total = size(distance_matrix,1)

    ! initialization
    select_distant_entries_dp(:) = 0

    if (ifirst .eq. 0) then
       ! first pair is the one with overall maximal distance
       maxpos  = maxloc(distance_matrix)
       select_distant_entries_dp(1) = maxpos(1)
       select_distant_entries_dp(2) = maxpos(2)
    else
       ! first pair is chosen randomly
       ! (a) initialize random stream
       call xor4096(iseed, maxpos(1), save_state=save_state)
       iseed = 0
       ! (b) generate first index
       call xor4096_range((/1_i4, n_total/), maxpos(1), save_state=save_state)
       ! (c) generate second index
       call xor4096_range((/1_i4, n_total/), maxpos(2), save_state=save_state)
       do while (maxpos(1) .eq. maxpos(2))
          call xor4096_range((/1_i4, n_total/), maxpos(2), save_state=save_state)
       end do
       select_distant_entries_dp(1) = maxpos(1)
       select_distant_entries_dp(2) = maxpos(2)
    end if

    do ii = 3, n_elements
       allocate( mask(n_total,n_total) )
       mask = .false.
       ! select rows and columns
       do jj = 1, ii-1
          ! select row
          mask(select_distant_entries_dp(jj),:) = .true. 
          ! select column
          mask(:,select_distant_entries_dp(jj)) = .true. 
          ! de-select diagonal elements and already chosen pairs
          do kk = 1, jj
             mask(select_distant_entries_dp(jj),select_distant_entries_dp(kk)) = .false.
             mask(select_distant_entries_dp(kk),select_distant_entries_dp(jj)) = .false.
          end do
       end do
       ! find largest distance of selected entries
       maxpos  = maxloc(distance_matrix, mask=mask)
       ! either row or column index is already part of the list
       ! --> add the other one
       if ( count(select_distant_entries_dp .eq. maxpos(1)) .eq. 0 ) then
          select_distant_entries_dp(ii) = maxpos(1)
       else
          select_distant_entries_dp(ii) = maxpos(2)
       end if
       deallocate( mask )
    end do

    call sort(select_distant_entries_dp)
    
  end function select_distant_entries_dp

  function select_distant_entries_sp( distance_matrix, n_elements, first, seed )
    
    implicit none

    real(sp), dimension(:,:), intent(in) :: distance_matrix           ! square matrix with distance between element i and j 
    integer(i4),              intent(in) :: n_elements                ! number of elements to return 
    integer(i4),    optional, intent(in) :: first                     ! method how first pair is chosen
    !                                                                 !    0 : first pair is the one with the overall 
    !                                                                 !        largest distance,
    !                                                                 !        i.e. result is deterministic
    !                                                                 !    1 : first pair is randomly chosen,
    !                                                                 !        i.e. result is random
    !                                                                 ! DEFAULT: 0 
    integer(i4),    optional, intent(in) :: seed                      ! seed used for random number generator
    !                                                                 ! (only considered if first=1)
    !                                                                 ! DEFAULT: timeseed
    integer(i4), dimension(n_elements)   :: select_distant_entries_sp ! indexes of elements with largest distance to each other

    ! local variables
    integer(i4)                             :: ii, jj, kk
    integer(i4)                             :: ifirst      ! method how first pair is chosen
    integer(i4)                             :: iseed       ! seed for random number generator
    integer(i4), dimension(n_save_state)    :: save_state
    integer(i4)                             :: n_total     ! number of elements --> number of rows/cols of distance_matrix
    integer(i4)                             :: maxpos(2)
    logical,    dimension(:,:), allocatable :: mask

    ! handling of optionals
    if ( present(first) ) then
       ifirst = first
    else
       ifirst = 0
    end if

    ! handling of optionals
    if ( present(seed) ) then
       iseed = seed
    else
       call get_timeseed(iseed)
    end if

    n_total = size(distance_matrix,1)

    ! initialization
    select_distant_entries_sp(:) = 0

    if (ifirst .eq. 0) then
       ! first pair is the one with overall maximal distance
       maxpos  = maxloc(distance_matrix)
       select_distant_entries_sp(1) = maxpos(1)
       select_distant_entries_sp(2) = maxpos(2)
    else
       ! first pair is chosen randomly
       ! (a) initialize random stream
       call xor4096(iseed, maxpos(1), save_state=save_state)
       iseed = 0
       ! (b) generate first index
       call xor4096_range((/1_i4, n_total/), maxpos(1), save_state=save_state)
       ! (c) generate second index
       call xor4096_range((/1_i4, n_total/), maxpos(2), save_state=save_state)
       do while (maxpos(1) .eq. maxpos(2))
          call xor4096_range((/1_i4, n_total/), maxpos(2), save_state=save_state)
       end do
       select_distant_entries_sp(1) = maxpos(1)
       select_distant_entries_sp(2) = maxpos(2)
    end if

    do ii = 3, n_elements
       allocate( mask(n_total,n_total) )
       mask = .false.
       ! select rows and columns
       do jj = 1, ii-1
          ! select row
          mask(select_distant_entries_sp(jj),:) = .true. 
          ! select column
          mask(:,select_distant_entries_sp(jj)) = .true. 
          ! de-select diagonal elements and already chosen pairs
          do kk = 1, jj
             mask(select_distant_entries_sp(jj),select_distant_entries_sp(kk)) = .false.
             mask(select_distant_entries_sp(kk),select_distant_entries_sp(jj)) = .false.
          end do
       end do
       ! find largest distance of selected entries
       maxpos  = maxloc(distance_matrix, mask=mask)
       ! either row or column index is already part of the list
       ! --> add the other one
       if ( count(select_distant_entries_sp .eq. maxpos(1)) .eq. 0 ) then
          select_distant_entries_sp(ii) = maxpos(1)
       else
          select_distant_entries_sp(ii) = maxpos(2)
       end if
       deallocate( mask )
    end do

    call sort(select_distant_entries_sp)
    
  end function select_distant_entries_sp

end module mo_select_distant_entries

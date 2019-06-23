module mo_remap

  ! -----------------------------------------------------------------------
  !
  !     this routine performs a remapping based on addresses and weights
  !     computed in a setup phase, e.g. with: cdo -f nc -gencon,grid.txt infile wts_file.nc
  !
  ! -----------------------------------------------------------------------
  !
  !     CVS:$Id: remap.f,v 1.5 2000/04/19 21:56:25 pwjones Exp $
  !
  !     Copyright (c) 1997, 1998 the Regents of the University of
  !       California.
  !
  !     This software and ancillary information (herein called software)
  !     called SCRIP is made available under the terms described here.
  !     The software has been approved for release with associated
  !     LA-CC Number 98-45.
  !
  !     Unless otherwise indicated, this software has been authored
  !     by an employee or employees of the University of California,
  !     operator of the Los Alamos National Laboratory under Contract
  !     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
  !     Government has rights to use, reproduce, and distribute this
  !     software.  The public may copy and use this software without
  !     charge, provided that this Notice and any statement of authorship
  !     are reproduced on all copies.  Neither the Government nor the
  !     University makes any warranty, express or implied, or assumes
  !     any liability or responsibility for the use of this software.
  !
  !     If software is modified to produce derivative works, such modified
  !     software should be clearly marked, so as not to confuse it with
  !     the version available from Los Alamos National Laboratory.
  !
  ! -----------------------------------------------------------------------

  ! Modified, Matthias Cuntz & Stephan Thober, Aug. 2012
  !                    - adapted to JAMS library, called mo_remap.f90
  !                    - miss

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 1997-2012 The Regents University of California, Matthias Cuntz, Stephan Thober
  !                         mc (at) macu (dot) de
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

  use mo_kind, only: i4, dp

  implicit none

  private

  public :: remap ! Performs remapping based on weights computed elsewhere

  ! -----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  !
  !     performs the remapping based on weights computed elsewhere
  !
  !-----------------------------------------------------------------------

  subroutine remap(dst_array, map_wts, dst_add, src_add, src_array, src_grad1, src_grad2, src_grad3, miss)

    implicit none

    integer(i4), dimension(:),   intent(in)           :: dst_add   ! destination address for each link
    integer(i4), dimension(:),   intent(in)           :: src_add   ! source      address for each link
    real(dp),    dimension(:,:), intent(in)           :: map_wts   ! remapping weights for each link
    real(dp),    dimension(:),   intent(in)           :: src_array ! array with source field to be remapped

    real(dp),    dimension(:),   intent(in), optional :: src_grad1 ! gradient arrays on source grid necessary for
    real(dp),    dimension(:),   intent(in), optional :: src_grad2 ! higher-order remappings
    real(dp),    dimension(:),   intent(in), optional :: src_grad3

    real(dp),                    intent(in), optional :: miss      ! Set missing cells to miss value (default: 0)

    real(dp),    dimension(:),   intent(inout)        :: dst_array ! array for remapped field on destination grid

    ! Local
    integer(i4) :: n, iorder
    real(dp) :: imiss
    logical, dimension(:), allocatable :: mask

    !
    ! check the order of the interpolation
    if (present(src_grad1)) then
       iorder = 2
    else
       iorder = 1
    endif
    !
    !
    if (present(miss)) then
       imiss = miss
    else
       imiss = 0.0_dp
    endif
    !
    ! remap
    dst_array = 0.0_dp
    allocate(mask(size(dst_array)))
    mask = .false.
    select case (iorder)
    case(1)
       !
       ! first order remapping
       do n=1,size(dst_add)
          dst_array(dst_add(n)) = dst_array(dst_add(n)) + &
               src_array(src_add(n))*map_wts(1,n)
          mask(dst_add(n)) = .true.
       end do
    case(2)
       !
       ! second order remapping
       if (size(map_wts,DIM=1) == 3) then
          do n=1, size(dst_add)
             dst_array(dst_add(n)) = dst_array(dst_add(n)) + &
                  src_array(src_add(n))*map_wts(1,n) + &
                  src_grad1(src_add(n))*map_wts(2,n) + &
                  src_grad2(src_add(n))*map_wts(3,n)
             mask(dst_add(n)) = .true.
          end do
       else if (size(map_wts,DIM=1) == 4) then
          do n=1,size(dst_add)
             dst_array(dst_add(n)) = dst_array(dst_add(n)) + &
                  src_array(src_add(n))*map_wts(1,n) + &
                  src_grad1(src_add(n))*map_wts(2,n) + &
                  src_grad2(src_add(n))*map_wts(3,n) + &
                  src_grad3(src_add(n))*map_wts(4,n)
             mask(dst_add(n)) = .true.
          end do
       endif
    end select

    dst_array = merge(dst_array, imiss, mask)
    deallocate(mask)

  end subroutine remap

end module mo_remap

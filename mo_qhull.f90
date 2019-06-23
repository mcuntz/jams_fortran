! License
! -------
! This file is part of the JAMS Fortran package, distributed under the MIT License.
!
! Copyright (c) 2016 Matthias Cuntz - mc (at) macu (dot) de
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
module mo_qhull

  use mo_kind, only: i4, dp

  implicit none

  PUBLIC :: qhull ! convex hull calculation from number of points

  integer(i4) :: qhull_f
  external :: qhull_f

contains

  function qhull(points, flags, outfile)

    implicit none

    real(dp), dimension(:,:), intent(in) :: points
    character(len=*), optional           :: flags
    character(len=*), optional           :: outfile
    integer(i4)                          :: qhull

    character(len=250) :: iflags
#ifdef pgiFortran
    character(len=299)  :: ioutfile
#else
    character(len=1024) :: ioutfile
#endif

    if (present(flags)) then
       iflags = trim(flags)
    else
       iflags = "qhull FS Fv n"
    endif

    if (present(outfile)) then
       ioutfile = trim(outfile)
    else
       ioutfile = "qhull.out"
    endif

    qhull = qhull_f(trim(iflags), trim(ioutfile), size(points,1), size(points,2), points)

  end function qhull

end module mo_qhull

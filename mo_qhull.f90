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

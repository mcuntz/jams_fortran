module mo_map_functions

  use mo_kind, only: int8=>i1, int16=>i2, int32=>i4, int64=>i8, real32=>sp, real64=>dp
#if defined (__pgiFortran__) || defined (__NAGf90Fortran__) || defined (__NAG__)
  use mo_kind, only: real128=>dp
#else
  use mo_kind, only: real128=>qp
#endif

  implicit none

  public

contains

  pure integer(kind=int8) function xpowx_i1(x) result(res)
    integer(kind=int8),intent(in) :: x
    res = x**x
  end function xpowx_i1

  pure integer(kind=int16) function xpowx_i2(x) result(res)
    integer(kind=int16),intent(in) :: x
    res = x**x
  end function xpowx_i2

  pure integer(kind=int32) function xpowx_i4(x) result(res)
    integer(kind=int32),intent(in) :: x
    res = x**x
  end function xpowx_i4

  pure integer(kind=int64) function xpowx_i8(x) result(res)
    integer(kind=int64),intent(in) :: x
    res = x**x
  end function xpowx_i8

  pure real(kind=real32) function xpowx_r4(x) result(res)
    real(kind=real32),intent(in) :: x
    res = x**x
  end function xpowx_r4

  pure real(kind=real64) function xpowx_r8(x) result(res)
    real(kind=real64),intent(in) :: x
    res = x**x
  end function xpowx_r8

  pure real(kind=real128) function xpowx_r16(x) result(res)
    real(kind=real128),intent(in) :: x
    res = x**x
  end function xpowx_r16

  pure complex(kind=real32) function xpowx_c4(x) result(res)
    complex(kind=real32),intent(in) :: x
    res = x**x
  end function xpowx_c4

  pure complex(kind=real64) function xpowx_c8(x) result(res)
    complex(kind=real64),intent(in) :: x
    res = x**x
  end function xpowx_c8

  pure complex(kind=real128) function xpowx_c16(x) result(res)
    complex(kind=real128),intent(in) :: x
    res = x**x
  end function xpowx_c16

end module mo_map_functions

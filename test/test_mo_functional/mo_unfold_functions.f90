module mo_unfold_functions

  use mo_kind, only: int8=>i1, int16=>i2, int32=>i4, int64=>i8, real32=>sp, real64=>dp
#if defined (pgiFortran) || defined (NAGf90Fortran) || defined (NAG)
  use mo_kind, only: real128=>dp
#else
  use mo_kind, only: real128=>qp
#endif

  implicit none

  public

contains

  pure integer(kind=int8) function addone_i1(x) result(res)
    use mo_kind, only: int8=>i1
    integer(kind=int8),intent(in) :: x
    res = x+1
  end function addone_i1

  pure integer(kind=int16) function addone_i2(x) result(res)
    use mo_kind, only: int16=>i2
    integer(kind=int16),intent(in) :: x
    res = x+1
  end function addone_i2

  pure integer(kind=int32) function addone_i4(x) result(res)
    use mo_kind, only: int32=>i4
    integer(kind=int32),intent(in) :: x
    res = x+1
  end function addone_i4

  pure integer(kind=int64) function addone_i8(x) result(res)
    use mo_kind, only: int64=>i8
    integer(kind=int64),intent(in) :: x
    res = x+1
  end function addone_i8

  pure real(kind=real32) function addone_r4(x) result(res)
    use mo_kind, only: real32=>sp
    real(kind=real32),intent(in) :: x
    res = x+1
  end function addone_r4

  pure real(kind=real64) function addone_r8(x) result(res)
    use mo_kind, only: real64=>dp
    real(kind=real64),intent(in) :: x
    res = x+1
  end function addone_r8

  pure real(kind=real128) function addone_r16(x) result(res)
#if defined (pgiFortran) || defined (NAGf90Fortran) || defined (NAG)
    use mo_kind, only: real128=>dp
#else
    use mo_kind, only: real128=>qp
#endif
    real(kind=real128),intent(in) :: x
    res = x+1
  end function addone_r16

  pure complex(kind=real32) function addone_c4(x) result(res)
    use mo_kind, only: real32=>sp
    complex(kind=real32),intent(in) :: x
    res = x+1
  end function addone_c4

  pure complex(kind=real64) function addone_c8(x) result(res)
    use mo_kind, only: real64=>dp
    complex(kind=real64),intent(in) :: x
    res = x+1
  end function addone_c8

  pure complex(kind=real128) function addone_c16(x) result(res)
#if defined (pgiFortran) || defined (NAGf90Fortran) || defined (NAG)
    use mo_kind, only: real128=>dp
#else
    use mo_kind, only: real128=>qp
#endif
    complex(kind=real128),intent(in) :: x
    res = x+1
  end function addone_c16

end module mo_unfold_functions

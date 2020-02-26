module mo_filter_functions

  use mo_kind, only: int8=>i1, int16=>i2, int32=>i4, int64=>i8, real32=>sp, real64=>dp
#if defined (__pgiFortran__) || defined (__NAGf90Fortran__) || defined (__NAG__)
  use mo_kind, only: real128=>dp
#else
  use mo_kind, only: real128=>qp
#endif

  implicit none

  public

contains

  pure logical function gt3lt5_i1(x) result(res)
    use mo_kind, only: int8=>i1
    integer(kind=int8),intent(in) :: x
    res = x > 3 .and. x < 5
  end function gt3lt5_i1

  pure logical function gt3lt5_i2(x) result(res)
    use mo_kind, only: int16=>i2
    integer(kind=int16),intent(in) :: x
    res = x > 3 .and. x < 5
  end function gt3lt5_i2

  pure logical function gt3lt5_i4(x) result(res)
    use mo_kind, only: int32=>i4
    integer(kind=int32),intent(in) :: x
    res = x > 3 .and. x < 5
  end function gt3lt5_i4

  pure logical function gt3lt5_i8(x) result(res)
    use mo_kind, only: int64=>i8
    integer(kind=int64),intent(in) :: x
    res = x > 3 .and. x < 5
  end function gt3lt5_i8

  pure logical function gt3lt5_r4(x) result(res)
    use mo_kind, only: real32=>sp
    real(kind=real32),intent(in) :: x
    res = x > 3 .and. x < 5
  end function gt3lt5_r4

  pure logical function gt3lt5_r8(x) result(res)
    use mo_kind, only: real64=>dp
    real(kind=real64),intent(in) :: x
    res = x > 3 .and. x < 5
  end function gt3lt5_r8

  pure logical function gt3lt5_r16(x) result(res)
#if defined (__pgiFortran__) || defined (__NAGf90Fortran__) || defined (__NAG__)
    use mo_kind, only: real128=>dp
#else
    use mo_kind, only: real128=>qp
#endif
    real(kind=real128),intent(in) :: x
    res = x > 3 .and. x < 5
  end function gt3lt5_r16

  pure logical function gt3lt5_c4(x) result(res)
    use mo_kind, only: real32=>sp
    complex(kind=real32),intent(in) :: x
    res = real(x) > 3 .and. real(x) < 5
  end function gt3lt5_c4

  pure logical function gt3lt5_c8(x) result(res)
    use mo_kind, only: real64=>dp
    complex(kind=real64),intent(in) :: x
    res = real(x) > 3 .and. real(x) < 5
  end function gt3lt5_c8

  pure logical function gt3lt5_c16(x) result(res)
#if defined (__pgiFortran__) || defined (__NAGf90Fortran__) || defined (__NAG__)
    use mo_kind, only: real128=>dp
#else
    use mo_kind, only: real128=>qp
#endif
    complex(kind=real128),intent(in) :: x
    res = real(x) > 3 .and. real(x) < 5
  end function gt3lt5_c16

end module mo_filter_functions

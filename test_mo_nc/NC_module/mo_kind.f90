module mo_kind
  ! Integer and floating point section
  ! see:   http://fortranwiki.org/fortran/show/Real+precision
  integer, parameter                        :: i4 = selected_int_kind(9)
  integer, parameter                        :: sp = selected_real_kind(6,37)
  integer, parameter                        :: dp = selected_real_kind(15,307)
end module mo_kind

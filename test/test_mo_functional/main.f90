program main

  use mo_kind,       only: int8=>i1, int16=>i2, int32=>i4, int64=>i8, real32=>sp, real64=>dp
#if defined (__pgiFortran__) || defined (__NAGf90Fortran__) || defined (__NAG__)
  use mo_kind,       only: real128=>dp
#else
  use mo_kind,       only: real128=>qp
#endif
  use mo_testing,    only: assert, initialize_tests, report_tests
  use mo_ansi_colors, only: color, c_red, c_green

  use mo_functional

  use mo_filter_functions
  use mo_fold_functions
  use mo_iterfold_functions
  use mo_map_functions
  use mo_unfold_functions

  implicit none

  logical, dimension(:), allocatable :: tests
  logical :: test_failed
  integer :: n, ntests
  integer, parameter :: stdout = 6

  real(kind=real32),dimension(1000) :: x

  complex(kind=real32),dimension(:),allocatable :: c4,c4_res
  complex(kind=real64),dimension(:),allocatable :: c8,c8_res
  complex(kind=real128),dimension(:),allocatable :: c16,c16_res
  complex(kind=real64) :: c8_start
  complex(kind=real128) :: c16_start
  complex(kind=real32),dimension(:),allocatable :: c_r4
  complex(kind=real64),dimension(:),allocatable :: c_r8
  complex(kind=real128),dimension(:),allocatable :: c_r16

  logical  :: isgood

  isgood = .true.

  
  n = 1
  ntests = 19
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(arange(1_int8,3_int8) == [1_int8,2_int8,3_int8]),&
       'arange, int8')
  n = n + 1

  tests(n) = assert(all(arange(1_int16,3_int16) == [1_int16,2_int16,3_int16]),&
       'arange, int16')
  n = n + 1

  tests(n) = assert(all(arange(1_int32,3_int32) == [1_int32,2_int32,3_int32]),&
       'arange, int32')
  n = n + 1

  tests(n) = assert(all(arange(1_int64,3_int64) == [1_int64,2_int64,3_int64]),&
       'arange, int64')
  n = n + 1

  tests(n) = assert(all(arange(1._real32,3._real32) == [1._real32,2._real32,3._real32]),&
       'arange, real32')
  n = n + 1

  tests(n) = assert(all(arange(1._real64,3._real64) == [1._real64,2._real64,3._real64]),&
       'arange, real32')
  n = n + 1

  tests(n) = assert(all(arange(1._real128,3._real128) == [1._real128,2._real128,3._real128]),&
       'arange, real128')
  n = n + 1

  tests(n) = assert(all(arange(1._real128,3._real128) == [1._real128,2._real128,3._real128]),&
       'arange, real128')
  n = n + 1

  tests(n) = assert(all(arange(cmplx(1._real32,0._real32, kind=real32),&
       cmplx(3._real32,0._real32, kind=real32),&
       cmplx(1._real32,0._real32, kind=real32))&
       == [cmplx(1._real32,0._real32, kind=real32),&
       cmplx(2._real32,0._real32, kind=real32),&
       cmplx(3._real32,0._real32, kind=real32)]),&
       'arange, complex real32')
  n = n + 1

  tests(n) = assert(all(arange(cmplx(1._real64,0._real64, kind=real64),&
       cmplx(3._real64,0._real64, kind=real64),&
       cmplx(1._real64,0._real64, kind=real64))&
       == [cmplx(1._real64,0._real64, kind=real64),&
       cmplx(2._real64,0._real64, kind=real64),&
       cmplx(3._real64,0._real64, kind=real64)]),&
       'arange, complex real64')
  n = n + 1

  tests(n) = assert(all(arange(cmplx(1._real128,0._real128, kind=real128),&
       cmplx(3._real128,0._real128, kind=real128),&
       cmplx(1._real128,0._real128, kind=real128))&
       == [cmplx(1._real128,0._real128, kind=real128),&
       cmplx(2._real128,0._real128, kind=real128),&
       cmplx(3._real128,0._real128, kind=real128)]),&
       'arange, complex real128')
  n = n + 1

  tests(n) = assert(all(arange(cmplx(1,1, kind=real32),cmplx(3,5, kind=real32),cmplx(1,2, kind=real32))&
       == [cmplx(1,1, kind=real32),cmplx(2,3, kind=real32),cmplx(3,5, kind=real32)]),&
       'arange, incrementing both parts of complex numbers')
  n = n + 1

  tests(n) = assert(all(arange(1,10) == arange(1,10,1)),&
       'arange increment equals 1 when ommited')
  n = n + 1

  tests(n) = assert(all(arange(1.,10.) == real(arange(1,10,1))),&
       'integer and real arange variants produce same values')
  n = n + 1

  tests(n) = assert(all(arange(0.,2.4,0.8) == [0.,0.8,1.6,2.4]),&
       'custom increment value')
  n = n + 1

  tests(n) = assert(all(arange(3,-1,-1) == [3,2,1,0,-1]),&
       'negative increment value')
  n = n + 1

  tests(n) = assert(size(arange(1.0,1.4,0.1)) == 5,&
       'real32-typed arange returns array of expected size')
  n = n + 1

  tests(n) = assert(size(arange(1.0_real64,1.4_real64,0.1_real64)) == 5,&
       'real64-typed arange returns array of expected size')
  n = n + 1

  tests(n) = assert(size(arange(1.0_real128,1.4_real128,0.1_real128)) == 5,&
       'real128-typed arange returns array of expected size')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)

  n = 1
  ntests = 11
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(complement([1_int8,2_int8],[2_int8]) == [1]),&
       'complement, int8')
  n = n + 1

  tests(n) = assert(all(complement([1_int16,2_int16],[2_int16]) == [1]),&
       'complement, int16')
  n = n + 1

  tests(n) = assert(all(complement([1_int32,2_int32],[2_int32]) == [1]),&
       'complement, int32')
  n = n + 1

  tests(n) = assert(all(complement([1_int64,2_int64],[2_int64]) == [1]),&
       'complement, int64')
  n = n + 1

  tests(n) = assert(all(complement([1._real32,2._real32],[2._real32]) == [1]),&
       'complement, real32')
  n = n + 1

  tests(n) = assert(all(complement([1._real64,2._real64],[2._real64]) == [1]),&
       'complement, real64')
  n = n + 1

  tests(n) = assert(all(complement([1._real128,2._real128],[2._real128]) == [1]),&
       'complement, real128')
  n = n + 1

  tests(n) = assert(all(complement([cmplx(1._real32,0._real32, kind=real32),&
       cmplx(2._real32,0._real32, kind=real32)],&
       [cmplx(2._real32,0._real32, kind=real32)])&
       == [cmplx(1._real32,0._real32, kind=real32)]),'complement, complex real32')
  n = n + 1

  tests(n) = assert(all(complement([cmplx(1._real64,0._real64, kind=real64),&
       cmplx(2._real64,0._real64, kind=real64)],&
       [cmplx(2._real64,0._real64, kind=real64)])&
       == [cmplx(1._real64,0._real64, kind=real64)]),'complement, complex real64')
  n = n + 1

  tests(n) = assert(all(complement([cmplx(1._real64,0._real64, kind=real64),&
       cmplx(2._real64,0._real64, kind=real64)],&
       [cmplx(2._real64,0._real64, kind=real64)])&
       == [cmplx(1._real64,0._real64, kind=real64)]),'complement, complex real64')
  n = n + 1

  tests(n) = assert(all(complement([1,2],[2]) == ([1,2].complement.[2])),&
       'complement operator, x.complement.y')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 11
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(filter(gt3lt5_i1,[3_int8,4_int8,5_int8]) == [4]),&
       'filter, int8')
  n = n + 1

  tests(n) = assert(all(filter(gt3lt5_i2,[3_int16,4_int16,5_int16]) == [4]),&
       'filter, int16')
  n = n + 1

  tests(n) = assert(all(filter(gt3lt5_i4,[3,4,5]) == [4]),&
       'filter, int32')
  n = n + 1

  tests(n) = assert(all(filter(gt3lt5_i8,[3_int64,4_int64,5_int64]) == [4]),&
       'filter, int64')
  n = n + 1

  tests(n) = assert(all(filter(gt3lt5_r4,[3.,4.,5.]) == [4]),&
       'filter, real32')
  n = n + 1

  tests(n) = assert(all(filter(gt3lt5_r8,[3._real64,4._real64,5._real64]) == [4]),&
       'filter, real64')
  n = n + 1

  tests(n) = assert(all(filter(gt3lt5_r16,[3._real128,4._real128,5._real128]) == [4]),&
       'filter, real128')
  n = n + 1

  tests(n) = assert(all(filter(gt3lt5_c4,&
       [cmplx(3.,0., kind=real32),cmplx(4.,0., kind=real32),cmplx(5.,0., kind=real32)]) == [cmplx(4.,0., kind=real32)]),&
       'filter, complex real32')
  n = n + 1

  ! Need to assign to a variable first because cmplx() by default
  ! returns single-precision complex number which breaks the generic
  ! interface
  c8 = [cmplx(3._real64,0._real64, kind=real64),&
       cmplx(4._real64,0._real64, kind=real64),&
       cmplx(5._real64,0._real64, kind=real64)]
  c16 = [cmplx(3._real128,0._real128, kind=real128),&
       cmplx(4._real128,0._real128, kind=real128),&
       cmplx(5._real128,0._real128, kind=real128)]

  tests(n) = assert(all(filter(gt3lt5_c8,c8) == [cmplx(4.,0., kind=real64)]),&
       'filter, complex real64')
  n = n + 1

  tests(n) = assert(all(filter(gt3lt5_c16,c16) == [cmplx(4.,0., kind=real128)]),&
       'filter, complex real128')
  n = n + 1

  tests(n) = assert(size(filter(gt3lt5_i4,[1,2,3,5,6])) == 0,&
       'filter returns empty array')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 10
  call initialize_tests(tests,ntests)

  tests(n) = assert(foldl(sum_i1,0_int8,[1_int8,2_int8,3_int8,4_int8,5_int8]) == 15,&
       'foldl, int8')
  n = n + 1

  tests(n) = assert(foldl(sum_i2,0_int16,[1_int16,2_int16,3_int16,4_int16,5_int16]) == 15,&
       'foldl, int16')
  n = n + 1

  tests(n) = assert(foldl(sum_i4,0_int32,[1_int32,2_int32,3_int32,4_int32,5_int32]) == 15,&
       'foldl, int32')
  n = n + 1

  tests(n) = assert(foldl(sum_i8,0_int64,[1_int64,2_int64,3_int64,4_int64,5_int64]) == 15,&
       'foldl, int64')
  n = n + 1

  tests(n) = assert(foldl(sum_r4,0._real32,[1._real32,2._real32,3._real32,4._real32,5._real32]) == 15,&
       'foldl, real32')
  n = n + 1

  tests(n) = assert(foldl(sum_r8,0._real64,[1._real64,2._real64,3._real64,4._real64,5._real64]) == 15,&
       'foldl, real64')
  n = n + 1

  tests(n) = assert(foldl(sum_r16,0._real128,[1._real128,2._real128,3._real128,4._real128,5._real128]) == 15,&
       'foldl, real128')
  n = n + 1

  c4 = arange(cmplx(1,0, kind=real32),cmplx(5,0, kind=real32))
  c8 = c4
  c16 = c4

  c8_start = cmplx(0,0)
  c16_start = c8_start

  tests(n) = assert(foldl(sum_c4,cmplx(0.,0., kind=real32),c4) == cmplx(15,0, kind=real32),&
       'foldl, complex real32')
  n = n + 1

  tests(n) = assert(foldl(sum_c8,c8_start,c8) == cmplx(15._real64,0._real64, kind=real64),&
       'foldl, complex real64')
  n = n + 1

  tests(n) = assert(foldl(sum_c16,c16_start,c16) == cmplx(15._real128,0._real128, kind=real128),&
       'foldl, complex real128')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 10
  call initialize_tests(tests,ntests)

  tests(n) = assert(foldr(sum_i1,0_int8,[1_int8,2_int8,3_int8,4_int8,5_int8]) == 15,&
       'foldr, int8')
  n = n + 1

  tests(n) = assert(foldr(sum_i2,0_int16,[1_int16,2_int16,3_int16,4_int16,5_int16]) == 15,&
       'foldr, int16')
  n = n + 1

  tests(n) = assert(foldr(sum_i4,0_int32,[1_int32,2_int32,3_int32,4_int32,5_int32]) == 15,&
       'foldr, int32')
  n = n + 1

  tests(n) = assert(foldr(sum_i8,0_int64,[1_int64,2_int64,3_int64,4_int64,5_int64]) == 15,&
       'foldr, int64')
  n = n + 1

  tests(n) = assert(foldr(sum_r4,0._real32,[1._real32,2._real32,3._real32,4._real32,5._real32]) == 15,&
       'foldr, real32')
  n = n + 1

  tests(n) = assert(foldr(sum_r8,0._real64,[1._real64,2._real64,3._real64,4._real64,5._real64]) == 15,&
       'foldr, real64')
  n = n + 1

  tests(n) = assert(foldr(sum_r16,0._real128,[1._real128,2._real128,3._real128,4._real128,5._real128]) == 15,&
       'foldr, real128')
  n = n + 1

  c4 = arange(cmplx(1,0, kind=real32),cmplx(5,0, kind=real32))
  c8 = c4
  c16 = c4

  c8_start = cmplx(0,0, kind=real32)
  c16_start = c8_start

  tests(n) = assert(foldr(sum_c4,cmplx(0.,0., kind=real32),c4) == cmplx(15,0, kind=real32),&
       'foldr, complex real32')
  n = n + 1

  tests(n) = assert(foldr(sum_c8,c8_start,c8) == cmplx(15._real64,0._real64, kind=real64),&
       'foldr, complex real64')
  n = n + 1

  tests(n) = assert(foldr(sum_c16,c16_start,c16) == cmplx(15._real128,0._real128, kind=real128),&
       'foldr, complex real128')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 10
  call initialize_tests(tests,ntests)

  tests(n) = assert(foldt(sum_i1,0_int8,[1_int8,2_int8,3_int8,4_int8,5_int8]) == 15,&
       'foldt, int8')
  n = n + 1

  tests(n) = assert(foldt(sum_i2,0_int16,[1_int16,2_int16,3_int16,4_int16,5_int16]) == 15,&
       'foldt, int16')
  n = n + 1

  tests(n) = assert(foldt(sum_i4,0_int32,[1_int32,2_int32,3_int32,4_int32,5_int32]) == 15,&
       'foldt, int32')
  n = n + 1

  tests(n) = assert(foldt(sum_i8,0_int64,[1_int64,2_int64,3_int64,4_int64,5_int64]) == 15,&
       'foldt, int64')
  n = n + 1

  tests(n) = assert(foldt(sum_r4,0._real32,[1._real32,2._real32,3._real32,4._real32,5._real32]) == 15,&
       'foldt, real32')
  n = n + 1

  tests(n) = assert(foldt(sum_r8,0._real64,[1._real64,2._real64,3._real64,4._real64,5._real64]) == 15,&
       'foldt, real64')
  n = n + 1

  tests(n) = assert(foldt(sum_r16,0._real128,[1._real128,2._real128,3._real128,4._real128,5._real128]) == 15,&
       'foldt, real128')
  n = n + 1

  c4 = arange(cmplx(1,0, kind=real32),cmplx(5,0, kind=real32))
  c8 = c4
  c16 = c4

  c8_start = cmplx(0,0, kind=real32)
  c16_start = c8_start

  tests(n) = assert(foldt(sum_c4,cmplx(0.,0., kind=real32),c4) == cmplx(15,0, kind=real32),&
       'foldt, complex real32')
  n = n + 1

  tests(n) = assert(foldt(sum_c8,c8_start,c8) == cmplx(15._real64,0._real64, kind=real64),&
       'foldt, complex real64')
  n = n + 1

  tests(n) = assert(foldt(sum_c16,c16_start,c16) == cmplx(15._real128,0._real128, kind=real128),&
       'foldt, complex real128')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)



  c_r4 = [(1,2),(2,4)]
  c_r8 = c_r4
  c_r16 = c_r4

  n = 1
  ntests = 11
  call initialize_tests(tests,ntests)

  tests(n) = assert(head([1_int8,2_int8]) == 1_int8,'head, int8')
  n = n + 1

  tests(n) = assert(head([1_int16,2_int16]) == 1_int16,'head, int16')
  n = n + 1

  tests(n) = assert(head([1_int32,2_int32]) == 1_int32,'head, int32')
  n = n + 1

  tests(n) = assert(head([1_int64,2_int64]) == 1_int64,'head, int64')
  n = n + 1

  tests(n) = assert(head([1._real32,2._real32]) == 1._real32,'head, real32')
  n = n + 1

  tests(n) = assert(head([1._real64,2._real64]) == 1._real64,'head, real64')
  n = n + 1

  tests(n) = assert(head([1._real128,2._real128]) == 1._real128,'head, real128')
  n = n + 1

  tests(n) = assert(head(c_r4) == c_r4(1),'head, complex real32')
  n = n + 1

  tests(n) = assert(head(c_r8) == c_r8(1),'head, complex real64')
  n = n + 1

  tests(n) = assert(head(c_r16) == c_r16(1),'head, complex real128')
  n = n + 1

#ifdef __pgiFortran__
  tests(n) = assert(.true.,'skip head operator, .head.x')
#else
  tests(n) = assert(head([1,2]) == .head.[1,2],'head operator, .head.x')
#endif
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  c_r4 = [(1,2),(2,4)]
  c_r8 = c_r4
  c_r16 = c_r4

  n = 1
  ntests = 13
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(init([1_int8,2_int8]) == [1_int8]),'init, int8')
  n = n + 1

  tests(n) = assert(all(init([1_int16,2_int16]) == [1_int16]),'init, int16')
  n = n + 1

  tests(n) = assert(all(init([1_int32,2_int32]) == [1_int32]),'init, int32')
  n = n + 1

  tests(n) = assert(all(init([1_int64,2_int64]) == [1_int64]),'init, int64')
  n = n + 1

  tests(n) = assert(all(init([1._real32,2._real32]) == [1._real32]),'init, real32')
  n = n + 1

  tests(n) = assert(all(init([1._real64,2._real64]) == [1._real64]),'init, real64')
  n = n + 1

  tests(n) = assert(all(init([1._real128,2._real128]) == [1._real128]),'init, real128')
  n = n + 1

  tests(n) = assert(all(init(c_r4) == [c_r4(1)]),'init, complex real32')
  n = n + 1

  tests(n) = assert(all(init(c_r8) == [c_r8(1)]),'init, complex real64')
  n = n + 1

  tests(n) = assert(all(init(c_r16) == [c_r16(1)]),'init, complex real128')
  n = n + 1

  tests(n) = assert(size(init([1])) == 0,'size(init([1])) == 0')
  n = n + 1

  tests(n) = assert(size(init(init([1]))) == 0,'size(init(init([1]))) == 0')
  n = n + 1

#ifdef __pgiFortran__
  tests(n) = assert(.true.,'skip init operator, .init.x')
#else
  tests(n) = assert(all(init([1,2]) == .init.[1,2]),'init operator, .init.x')
#endif
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 12
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(insert(2_int8,2,[1_int8,3_int8]) == [1,2,3]),&
       'insert, int8')
  n = n + 1

  tests(n) = assert(all(insert(2_int16,2,[1_int16,3_int16]) == [1,2,3]),&
       'insert, int16')
  n = n + 1

  tests(n) = assert(all(insert(2_int32,2,[1_int32,3_int32]) == [1,2,3]),&
       'insert, int32')
  n = n + 1

  tests(n) = assert(all(insert(2_int64,2,[1_int64,3_int64]) == [1,2,3]),&
       'insert, int64')
  n = n + 1

  tests(n) = assert(all(insert(2._real32,2,[1._real32,3._real32]) == [1,2,3]),&
       'insert, real32')
  n = n + 1

  tests(n) = assert(all(insert(2._real64,2,[1._real64,3._real64]) == [1,2,3]),&
       'insert, real64')
  n = n + 1

  tests(n) = assert(all(insert(2._real128,2,[1._real128,3._real128]) == [1,2,3]),&
       'insert, real128')
  n = n + 1

  tests(n) = assert(all(insert(cmplx(2._real32,0._real32, kind=real32),2,&
       [cmplx(1._real32,0._real32, kind=real32),cmplx(3._real32,0._real32, kind=real32)])&
       == arange(cmplx(1._real32,0._real32, kind=real32),cmplx(3._real32,0._real32, kind=real32))),&
       'insert, real32')
  n = n + 1

  tests(n) = assert(all(insert(cmplx(2._real64,0._real64, kind=real64),2,&
       [cmplx(1._real64,0._real64, kind=real64),cmplx(3._real64,0._real64, kind=real64)])&
       == arange(cmplx(1._real64,0._real64, kind=real64),cmplx(3._real64,0._real64, kind=real64))),&
       'insert, real64')
  n = n + 1

  tests(n) = assert(all(insert(cmplx(2._real128,0._real128, kind=real128),2,&
       [cmplx(1._real128,0._real128, kind=real128),cmplx(3._real128,0._real128, kind=real128)])&
       == arange(cmplx(1._real128,0._real128, kind=real128),cmplx(3._real128,0._real128, kind=real128))),&
       'insert, real128')
  n = n + 1

  tests(n) = assert(all(insert(1,1,arange(1,0)) == [1]),&
       'insert into empty array')
  n = n + 1

  tests(n) = assert(all(insert(2,2,[1]) == [1,2]),&
       'insert out of bounds')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 11
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(intersection([1_int8,2_int8],[2_int8,3_int8]) == [2]),&
       'intersection, int8')
  n = n + 1

  tests(n) = assert(all(intersection([1_int16,2_int16],[2_int16,3_int16]) == [2]),&
       'intersection, int16')
  n = n + 1

  tests(n) = assert(all(intersection([1_int32,2_int32],[2_int32,3_int32]) == [2]),&
       'intersection, int32')
  n = n + 1

  tests(n) = assert(all(intersection([1_int64,2_int64],[2_int64,3_int64]) == [2]),&
       'intersection, int64')
  n = n + 1

  tests(n) = assert(all(intersection([1._real32,2._real32],[2._real32,3._real32]) == [2]),&
       'intersection, real32')
  n = n + 1

  tests(n) = assert(all(intersection([1._real64,2._real64],[2._real64,3._real64]) == [2]),&
       'intersection, real64')
  n = n + 1

  tests(n) = assert(all(intersection([1._real128,2._real128],[2._real128,3._real128]) == [2]),&
       'intersection, real128')
  n = n + 1

  tests(n) = assert(all(intersection([cmplx(1._real32,0._real32, kind=real32),cmplx(2._real32,0._real32, kind=real32)],&
       [cmplx(2._real32,0._real32, kind=real32),cmplx(3._real32,0._real32, kind=real32)])&
       == [cmplx(2._real32,0._real32, kind=real32)]),&
       'intersection, complex real32')
  n = n + 1

  tests(n) = assert(all(intersection([cmplx(1._real64,0._real64, kind=real64),cmplx(2._real64,0._real64, kind=real64)],&
       [cmplx(2._real64,0._real64, kind=real64),cmplx(3._real64,0._real64, kind=real64)])&
       == [cmplx(2._real64,0._real64, kind=real64)]),&
       'intersection, complex real64')
  n = n + 1

  tests(n) = assert(all(intersection([cmplx(1._real128,0._real128, kind=real128),cmplx(2._real128,0._real128, kind=real128)],&
       [cmplx(2._real128,0._real128, kind=real128),cmplx(3._real128,0._real128, kind=real128)])&
       == [cmplx(2._real128,0._real128, kind=real128)]),&
       'intersection, complex real128')
  n = n + 1

  tests(n) = assert(all(intersection([1,2],[2,3]) == ([1,2].intersection.[2,3])),&
       'intersection operator, x.intersection.y')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 10
  call initialize_tests(tests,ntests)

  tests(n) = assert(iterfold(isum_i1,0_int8,[1_int8,2_int8,3_int8,4_int8,5_int8]) == 15,&
       'iterfold, int8')
  n = n + 1

  tests(n) = assert(iterfold(isum_i2,0_int16,[1_int16,2_int16,3_int16,4_int16,5_int16]) == 15,&
       'iterfold, int16')
  n = n + 1

  tests(n) = assert(iterfold(isum_i4,0_int32,[1_int32,2_int32,3_int32,4_int32,5_int32]) == 15,&
       'iterfold, int32')
  n = n + 1

  tests(n) = assert(iterfold(isum_i8,0_int64,[1_int64,2_int64,3_int64,4_int64,5_int64]) == 15,&
       'iterfold, int64')
  n = n + 1

  tests(n) = assert(iterfold(isum_r4,0._real32,[1._real32,2._real32,3._real32,4._real32,5._real32]) == 15,&
       'iterfold, real32')
  n = n + 1

  tests(n) = assert(iterfold(isum_r8,0._real64,[1._real64,2._real64,3._real64,4._real64,5._real64]) == 15,&
       'iterfold, real64')
  n = n + 1

  tests(n) = assert(iterfold(isum_r16,0._real128,[1._real128,2._real128,3._real128,4._real128,5._real128]) == 15,&
       'iterfold, real128')
  n = n + 1

  c4 = arange(cmplx(1,0, kind=real32),cmplx(5,0, kind=real32))
  c8 = c4
  c16 = c4

  c8_start = cmplx(0,0, kind=real32)
  c16_start = c8_start

  tests(n) = assert(iterfold(isum_c4,cmplx(0.,0., kind=real32),c4) == cmplx(15,0, kind=real32),&
       'iterfold, complex real32')
  n = n + 1

  tests(n) = assert(iterfold(isum_c8,c8_start,c8) == cmplx(15._real64,0._real64, kind=real64),&
       'iterfold, complex real64')
  n = n + 1

  tests(n) = assert(iterfold(isum_c16,c16_start,c16) == cmplx(15._real128,0._real128, kind=real128),&
       'iterfold, complex real128')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  c_r4 = [(1,2),(2,4)]
  c_r8 = c_r4
  c_r16 = c_r4

  n = 1
  ntests = 11
  call initialize_tests(tests,ntests)

  tests(n) = assert(last([1_int8,2_int8]) == 2_int8,'last, int8')
  n = n + 1

  tests(n) = assert(last([1_int16,2_int16]) == 2_int16,'last, int16')
  n = n + 1

  tests(n) = assert(last([1_int32,2_int32]) == 2_int32,'last, int32')
  n = n + 1

  tests(n) = assert(last([1_int64,2_int64]) == 2_int64,'last, int64')
  n = n + 1

  tests(n) = assert(last([1._real32,2._real32]) == 2._real32,'last, real32')
  n = n + 1

  tests(n) = assert(last([1._real64,2._real64]) == 2._real64,'last, real64')
  n = n + 1

  tests(n) = assert(last([1._real128,2._real128]) == 2._real128,'last, real128')
  n = n + 1

  tests(n) = assert(last(c_r4) == c_r4(2),'last, complex real32')
  n = n + 1

  tests(n) = assert(last(c_r8) == c_r8(2),'last, complex real64')
  n = n + 1

  tests(n) = assert(last(c_r16) == c_r16(2),'last, complex real128')
  n = n + 1

#ifdef __pgiFortran__
  tests(n) = assert(.true.,'skip last operator, .last.x')
#else
  tests(n) = assert(last([1,2]) == .last.[1,2],'last operator, .last.x')
#endif
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 11
  call initialize_tests(tests,ntests)

  tests(n) = assert(limit(2_int8,1_int8,3_int8) == 2_int8,&
       'limit, int8')
  n = n + 1

  tests(n) = assert(limit(2_int16,1_int16,3_int16) == 2_int16,&
       'limit, int16')
  n = n + 1

  tests(n) = assert(limit(2_int32,1_int32,3_int32) == 2_int32,&
       'limit, int32')
  n = n + 1

  tests(n) = assert(limit(2_int64,1_int64,3_int64) == 2_int64,&
       'limit, int64')
  n = n + 1

  tests(n) = assert(limit(2._real32,1._real32,3._real32) == 2._real32,&
       'limit, real32')
  n = n + 1

  tests(n) = assert(limit(2._real64,1._real64,3._real64) == 2._real64,&
       'limit, real32')
  n = n + 1

  tests(n) = assert(limit(2._real128,1._real128,3._real128) == 2._real128,&
       'limit, real128')
  n = n + 1

  tests(n) = assert(limit(cmplx(-0.5,1.5, kind=real32),cmplx(0,0, kind=real32),cmplx(1,1, kind=real32)) == cmplx(0,1, kind=real32),&
       'limit, complex real32')
  n = n + 1

  tests(n) = assert(limit(cmplx(-0.5_real64,1.5_real64, kind=real64),cmplx(0._real64,0._real64, kind=real64),&
       cmplx(1._real64,1._real64, kind=real64)) == cmplx(0._real64,1._real64, kind=real64),&
       'limit, complex real64')
  n = n + 1

  tests(n) = assert(limit(cmplx(-0.5_real128,1.5_real128, kind=real128),cmplx(0._real128,0._real128, kind=real128),&
       cmplx(1._real128,1._real128, kind=real128)) == cmplx(0._real128,1._real128, kind=real128),&
       'limit, complex real128')
  n = n + 1

  tests(n) = assert(all(limit(arange(1,3),2,2) == [2,2,2]),&
       'limit works on arrays')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 10
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(map(xpowx_i1,[1_int8,2_int8,3_int8])&
       == [1_int8,4_int8,27_int8]),'map, int8')
  n = n + 1

  tests(n) = assert(all(map(xpowx_i2,[1_int16,2_int16,3_int16])&
       == [1_int16,4_int16,27_int16]),'map, int16')
  n = n + 1

  tests(n) = assert(all(map(xpowx_i4,[1_int32,2_int32,3_int32])&
       == [1_int32,4_int32,27_int32]),'map, int32')
  n = n + 1

  tests(n) = assert(all(map(xpowx_i8,[1_int64,2_int64,3_int64])&
       == [1_int64,4_int64,27_int64]),'map, int64')
  n = n + 1

  tests(n) = assert(all(map(xpowx_r4,[1._real32,2._real32,3._real32])&
       == [1._real32,4._real32,27._real32]),'map, real32')
  n = n + 1

  tests(n) = assert(all(map(xpowx_r8,[1._real64,2._real64,3._real64])&
       == [1._real64,4._real64,27._real64]),'map, real64')
  n = n + 1

  tests(n) = assert(all(map(xpowx_r16,[1._real128,2._real128,3._real128])&
       == [1._real128,4._real128,27._real128]),'map, real128')
  n = n + 1

  c4 = [cmplx(1.,0., kind=real32),cmplx(2.,0., kind=real32),cmplx(3.,0., kind=real32)]
  c8 = [cmplx(1._real64,0._real64, kind=real64),&
       cmplx(2._real64,0._real64, kind=real64),&
       cmplx(3._real64,0._real64, kind=real64)]
  c16 = [cmplx(1._real128,0._real128, kind=real128),&
       cmplx(2._real128,0._real128, kind=real128),&
       cmplx(3._real128,0._real128, kind=real128)]

  tests(n) = assert(all(map(xpowx_c4,c4) == c4**c4),'map, complex real32')
  n = n + 1

  tests(n) = assert(all(map(xpowx_c8,c8) == c8**c8),'map, complex real64')
  n = n + 1

  tests(n) = assert(all(map(xpowx_c16,c16) == c16**c16),'map, complex real128')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 11
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(reverse(arange(1_int8,3_int8)) == [3,2,1]),&
       'reverse, int8')
  n = n + 1

  tests(n) = assert(all(reverse(arange(1_int16,3_int16)) == [3,2,1]),&
       'reverse, int16')
  n = n + 1

  tests(n) = assert(all(reverse(arange(1_int32,3_int32)) == [3,2,1]),&
       'reverse, int32')
  n = n + 1

  tests(n) = assert(all(reverse(arange(1_int64,3_int64)) == [3,2,1]),&
       'reverse, int64')
  n = n + 1

  tests(n) = assert(all(reverse(arange(1._real32,3._real32)) == [3,2,1]),&
       'reverse, real32')
  n = n + 1

  tests(n) = assert(all(reverse(arange(1._real64,3._real64)) == [3,2,1]),&
       'reverse, real64')
  n = n + 1

  tests(n) = assert(all(reverse(arange(1._real128,3._real128)) == [3,2,1]),&
       'reverse, real128')
  n = n + 1

  tests(n) = assert(all(reverse(arange(cmplx(1._real32,0._real32, kind=real32),&
       cmplx(3._real32,0._real32, kind=real32)))&
       == arange(cmplx(3._real32,0._real32, kind=real32),&
       cmplx(1._real32,0._real32, kind=real32),&
       cmplx(-1._real32,0._real32, kind=real32))),&
       'reverse, complex real32')
  n = n + 1

  tests(n) = assert(all(reverse(arange(cmplx(1._real64,0._real64, kind=real64),&
       cmplx(3._real64,0._real64, kind=real64)))&
       == arange(cmplx(3._real64,0._real64, kind=real64),&
       cmplx(1._real64,0._real64, kind=real64),&
       cmplx(-1._real64,0._real64, kind=real64))),&
       'reverse, complex real64')
  n = n + 1

  tests(n) = assert(all(reverse(arange(cmplx(1._real128,0._real128, kind=real64),&
       cmplx(3._real128,0._real128, kind=real64)))&
       == arange(cmplx(3._real128,0._real128, kind=real64),&
       cmplx(1._real128,0._real128, kind=real64),&
       cmplx(-1._real128,0._real128, kind=real64))),&
       'reverse, complex real128')
  n = n + 1

#ifdef __pgiFortran__
  tests(n) = assert(.true.,'skip reverse operator, .reverse.x')
#else
  tests(n) = assert(all(reverse([1,2,3]) == .reverse.[1,2,3]),&
       'reverse operator, .reverse.x')
#endif
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 12
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(set([1_int8,2_int8,2_int8,3_int8]) == [1,2,3]),&
       'set, int8')
  n = n + 1

  tests(n) = assert(all(set([1_int16,2_int16,2_int16,3_int16]) == [1,2,3]),&
       'set, int16')
  n = n + 1

  tests(n) = assert(all(set([1_int32,2_int32,2_int32,3_int32]) == [1,2,3]),&
       'set, int32')
  n = n + 1

  tests(n) = assert(all(set([1_int64,2_int64,2_int64,3_int64]) == [1,2,3]),&
       'set, int64')
  n = n + 1

  tests(n) = assert(all(set([1._real32,2._real32,2._real32,3._real32]) == [1,2,3]),&
       'set, real32')
  n = n + 1

  tests(n) = assert(all(set([1._real64,2._real64,2._real64,3._real64]) == [1,2,3]),&
       'set, real64')
  n = n + 1

  tests(n) = assert(all(set([1._real128,2._real128,2._real128,3._real128]) == [1,2,3]),&
       'set, real128')
  n = n + 1

  c4 = [cmplx(1,0, kind=real32),cmplx(2,0, kind=real32),cmplx(2,0, kind=real32),cmplx(3,0, kind=real32)]
  c4_res = [cmplx(1,0, kind=real32),cmplx(2,0, kind=real32),cmplx(3,0, kind=real32)]
  tests(n) = assert(all(set(c4) == c4_res),'set, complex real32')
  n = n + 1

  c8 = c4
  c8_res = c4_res
  tests(n) = assert(all(set(c8) == c8_res),'set, complex real64')
  n = n + 1

  c16 = c4
  c16_res = c4_res
  tests(n) = assert(all(set(c16) == c16_res),'set, complex real128')
  n = n + 1

  tests(n) = assert(all(set(arange(1,0)) == arange(1,0)),&
       'set of empty array is an empty array')
  n = n + 1

#ifdef __pgiFortran__
  tests(n) = assert(.true.,'skip set operator, .set.x')
#else
  tests(n) = assert(all(set([1,2,2,3]) == .set.[1,2,2,3]),&
       'set operator, .set.x')
#endif
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 12
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(sort([3_int8,2_int8,1_int8]) == [1,2,3]),&
       'sort, int8')
  n = n + 1

  tests(n) = assert(all(sort([3_int16,2_int16,1_int16]) == [1,2,3]),&
       'sort, int16')
  n = n + 1

  tests(n) = assert(all(sort([3_int32,2_int32,1_int32]) == [1,2,3]),&
       'sort, int32')
  n = n + 1

  tests(n) = assert(all(sort([3_int64,2_int64,1_int64]) == [1,2,3]),&
       'sort, int8')
  n = n + 1

  tests(n) = assert(all(sort([3._real32,2._real32,1._real32]) == [1,2,3]),&
       'sort, real32')
  n = n + 1

  tests(n) = assert(all(sort([3._real64,2._real64,1._real64]) == [1,2,3]),&
       'sort, real64')
  n = n + 1

  tests(n) = assert(all(sort([3._real128,2._real128,1._real128]) == [1,2,3]),&
       'sort, real128')
  n = n + 1

  tests(n) = assert(all(sort(arange(cmplx(3._real32,0._real32, kind=real32),&
       cmplx(1._real32,0._real32, kind=real32),&
       cmplx(-1._real32,0._real32, kind=real32)))&
       == arange(cmplx(1._real32,0._real32, kind=real32),&
       cmplx(3._real32,0._real32, kind=real32))),&
       'sort, complex real32')
  n = n + 1

  tests(n) = assert(all(sort(arange(cmplx(3._real64,0._real64, kind=real64),&
       cmplx(1._real64,0._real64, kind=real64),&
       cmplx(-1._real64,0._real64, kind=real64)))&
       == arange(cmplx(1._real64,0._real64, kind=real64),&
       cmplx(3._real64,0._real64, kind=real64))),&
       'sort, complex real64')
  n = n + 1

  tests(n) = assert(all(sort(arange(cmplx(3._real128,0._real128, kind=real128),&
       cmplx(1._real128,0._real128, kind=real128),&
       cmplx(-1._real128,0._real128, kind=real128)))&
       == arange(cmplx(1._real128,0._real128, kind=real128),&
       cmplx(3._real128,0._real128, kind=real128))),&
       'sort, complex real128')
  n = n + 1

  call random_number(x)
  tests(n) = assert(all(tail(sort(x)) >= init(sort(x))),&
       'all(tail(sort(x)) >= init(sort(x))')
  n = n + 1

#ifdef __pgiFortran__
  tests(n) = assert(.true.,'skip sort operator, .sort.x')
#else
  tests(n) = assert(all(sort(x) == .sort.x),&
       'sort operator, .sort.x')
#endif
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 16
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(split(arange(1_int8,10_int8),1) == arange(1,5)),&
       'split(x,1), int8')
  n = n + 1

  tests(n) = assert(all(split(arange(1_int8,10_int8),2) == arange(6,10)),&
       'split(x,2), int8')
  n = n + 1

  tests(n) = assert(all(split(arange(1_int16,10_int16),1) == arange(1,5)),&
       'split(x,1), int16')
  n = n + 1

  tests(n) = assert(all(split(arange(1_int16,10_int16),2) == arange(6,10)),&
       'split(x,2), int16')
  n = n + 1

  tests(n) = assert(all(split(arange(1_int32,10_int32),1) == arange(1,5)),&
       'split(x,1), int32')
  n = n + 1

  tests(n) = assert(all(split(arange(1_int32,10_int32),2) == arange(6,10)),&
       'split(x,2), int32')
  n = n + 1

  tests(n) = assert(all(split(arange(1_int64,10_int64),1) == arange(1,5)),&
       'split(x,1), int64')
  n = n + 1

  tests(n) = assert(all(split(arange(1_int64,10_int64),2) == arange(6,10)),&
       'split(x,2), int64')
  n = n + 1

  tests(n) = assert(all(split(arange(1._real32,10._real32),1) == arange(1,5)),&
       'split(x,1), real32')
  n = n + 1

  tests(n) = assert(all(split(arange(1._real32,10._real32),2) == arange(6,10)),&
       'split(x,2), real32')
  n = n + 1

  tests(n) = assert(all(split(arange(1._real64,10._real64),1) == arange(1,5)),&
       'split(x,1), real64')
  n = n + 1

  tests(n) = assert(all(split(arange(1._real64,10._real64),2) == arange(6,10)),&
       'split(x,2), real64')
  n = n + 1

  tests(n) = assert(all(split(arange(1._real128,10._real128),1) == arange(1,5)),&
       'split(x,1), real128')
  n = n + 1

  tests(n) = assert(all(split(arange(1._real128,10._real128),2) == arange(6,10)),&
       'split(x,2), real128')
  n = n + 1

  tests(n) = assert(all(split([1],1) == arange(1,0)),'split([1],1) returns an empty array')
  n = n + 1

  tests(n) = assert(all(split([1],2) == [1]),'split([1],2) returns [1]')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 11
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(subscript([1_int8,2_int8,3_int8],[2_int8]) == [2_int8]),&
       'subscript, int8')
  n = n + 1

  tests(n) = assert(all(subscript([1_int16,2_int16,3_int16],[2_int16]) == [2_int16]),&
       'subscript, int16')
  n = n + 1

  tests(n) = assert(all(subscript([1_int32,2_int32,3_int32],[2_int32]) == [2_int32]),&
       'subscript, int32')
  n = n + 1

  tests(n) = assert(all(subscript([1_int64,2_int64,3_int64],[2_int64]) == [2_int64]),&
       'subscript, int64')
  n = n + 1

  tests(n) = assert(all(subscript([1._real32,2._real32,3._real32],[2]) == [2._real32]),&
       'subscript, real32')
  n = n + 1

  tests(n) = assert(all(subscript([1._real64,2._real64,3._real64],[2]) == [2._real64]),&
       'subscript, real64')
  n = n + 1

  tests(n) = assert(all(subscript([1._real128,2._real128,3._real128],[2]) == [2._real128]),&
       'subscript, real128')
  n = n + 1

  tests(n) = assert(size(subscript([1,2,3],[0])) == 0,&
       'subscript out of bounds returns empty array')
  n = n + 1

  tests(n) = assert(all(subscript(arange(cmplx(1._real32,0._real32, kind=real32),&
       cmplx(3._real32,0._real32, kind=real32)),[2])&
       == [cmplx(2._real32,0._real32, kind=real32)]),&
       'subscript, complex real32')
  n = n + 1

  tests(n) = assert(all(subscript(arange(cmplx(1._real64,0._real64, kind=real64),&
       cmplx(3._real64,0._real64, kind=real64)),[2])&
       == [cmplx(2._real64,0._real64, kind=real64)]),&
       'subscript, complex real64')
  n = n + 1

  tests(n) = assert(all(subscript(arange(cmplx(1._real128,0._real128, kind=real128),&
       cmplx(3._real128,0._real128, kind=real128)),[2])&
       == [cmplx(2._real128,0._real128, kind=real128)]),&
       'subscript, complex real128')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  c_r4 = [(1,2),(2,4)]
  c_r8 = c_r4
  c_r16 = c_r4

  n = 1
  ntests = 13
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(tail([1_int8,2_int8]) == [2_int8]),'tail, int8')
  n = n + 1

  tests(n) = assert(all(tail([1_int16,2_int16]) == [2_int16]),'tail, int16')
  n = n + 1

  tests(n) = assert(all(tail([1_int32,2_int32]) == [2_int32]),'tail, int32')
  n = n + 1

  tests(n) = assert(all(tail([1_int64,2_int64]) == [2_int64]),'tail, int64')
  n = n + 1

  tests(n) = assert(all(tail([1._real32,2._real32]) == [2._real32]),'tail, real32')
  n = n + 1

  tests(n) = assert(all(tail([1._real64,2._real64]) == [2._real64]),'tail, real64')
  n = n + 1

  tests(n) = assert(all(tail([1._real128,2._real128]) == [2._real128]),'tail, real128')
  n = n + 1

  tests(n) = assert(all(tail(c_r4) == [c_r4(2)]),'tail, complex real32')
  n = n + 1

  tests(n) = assert(all(tail(c_r8) == [c_r8(2)]),'tail, complex real64')
  n = n + 1

  tests(n) = assert(all(tail(c_r16) == [c_r16(2)]),'tail, complex real128')
  n = n + 1

  tests(n) = assert(size(tail([1._real32])) == 0,'size(tail([1])) == 0')
  n = n + 1

  tests(n) = assert(size(tail(tail([1._real32]))) == 0,'size(tail(tail([1]))) == 0')
  n = n + 1

#ifdef __pgiFortran__
  tests(n) = assert(.true.,'skip tail operator, .tail.x')
#else
  tests(n) = assert(all(tail([1,2]) == .tail.[1,2]),'tail operator, .tail.x')
#endif
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 10
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(unfold(addone_i1,[1_int8],3_int8) == [1,2,3]),&
       'unfold, int8')
  n = n + 1

  tests(n) = assert(all(unfold(addone_i2,[1_int16],3_int16) == [1,2,3]),&
       'unfold, int16')
  n = n + 1

  tests(n) = assert(all(unfold(addone_i4,[1_int32],3_int32) == [1,2,3]),&
       'unfold, int32')
  n = n + 1

  tests(n) = assert(all(unfold(addone_i8,[1_int64],3_int64) == [1,2,3]),&
       'unfold, int64')
  n = n + 1

  tests(n) = assert(all(unfold(addone_r4,[1._real32],3_int32) == [1,2,3]),&
       'unfold, real32')
  n = n + 1

  tests(n) = assert(all(unfold(addone_r8,[1._real64],3_int32) == [1,2,3]),&
       'unfold, real64')
  n = n + 1

  tests(n) = assert(all(unfold(addone_r16,[1._real128],3_int32) == [1,2,3]),&
       'unfold, real128')
  n = n + 1

  tests(n) = assert(all(unfold(addone_c4,[cmplx(1._real32,0._real32, kind=real32)],3)&
       == arange(cmplx(1._real32,0._real32, kind=real32),&
       cmplx(3._real32,0._real32, kind=real32))),&
       'unfold, complex real32')
  n = n + 1

  c8 = [cmplx(1._real64,0._real64, kind=real64)]
  tests(n) = assert(all(unfold(addone_c8,c8,3)&
       == arange(cmplx(1._real64,0._real64, kind=real64),&
       cmplx(3._real64,0._real64, kind=real64))),&
       'unfold, complex real64')
  n = n + 1

  c16 = [cmplx(1._real128,0._real128, kind=real128)]
  tests(n) = assert(all(unfold(addone_c16,c16,3)&
       == arange(cmplx(1._real128,0._real128, kind=real128),&
       cmplx(3._real128,0._real128, kind=real128))),&
       'unfold, complex real128')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)


  n = 1
  ntests = 14
  call initialize_tests(tests,ntests)

  tests(n) = assert(all(union([1_int8,2_int8],[2_int8,3_int8]) == [1,2,3]),&
       'union, int8')
  n = n + 1

  tests(n) = assert(all(union([1_int16,2_int16],[2_int16,3_int16]) == [1,2,3]),&
       'union, int16')
  n = n + 1

  tests(n) = assert(all(union([1_int32,2_int32],[2_int32,3_int32]) == [1,2,3]),&
       'union, int32')
  n = n + 1

  tests(n) = assert(all(union([1_int64,2_int64],[2_int64,3_int64]) == [1,2,3]),&
       'union, int64')
  n = n + 1

  tests(n) = assert(all(union([1._real32,2._real32],[2._real32,3._real32]) == [1,2,3]),&
       'union, real32')
  n = n + 1

  tests(n) = assert(all(union([1._real64,2._real64],[2._real64,3._real64]) == [1,2,3]),&
       'union, real64')
  n = n + 1

  tests(n) = assert(all(union([1._real128,2._real128],[2._real128,3._real128]) == [1,2,3]),&
       'union, real128')
  n = n + 1

  tests(n) = assert(all(union([cmplx(1._real32,0._real32, kind=real32),cmplx(2._real32,0._real32, kind=real32)],&
       [cmplx(2._real32,0._real32, kind=real32),cmplx(3._real32,0._real32, kind=real32)])&
       == [cmplx(1._real32,0._real32, kind=real32),&
       cmplx(2._real32,0._real32, kind=real32),&
       cmplx(3._real32,0._real32, kind=real32)]),&
       'union, complex real32')
  n = n + 1

  tests(n) = assert(all(union([cmplx(1._real64,0._real64, kind=real64),cmplx(2._real64,0._real64, kind=real64)],&
       [cmplx(2._real64,0._real64, kind=real64),cmplx(3._real64,0._real64, kind=real64)])&
       == [cmplx(1._real64,0._real64, kind=real64),&
       cmplx(2._real64,0._real64, kind=real64),&
       cmplx(3._real64,0._real64, kind=real64)]),&
       'union, complex real64')
  n = n + 1

  tests(n) = assert(all(union([cmplx(1._real128,0._real128, kind=real128),cmplx(2._real128,0._real128, kind=real128)],&
       [cmplx(2._real128,0._real128, kind=real128),cmplx(3._real128,0._real128, kind=real128)])&
       == [cmplx(1._real128,0._real128, kind=real128),&
       cmplx(2._real128,0._real128, kind=real128),&
       cmplx(3._real128,0._real128, kind=real128)]),&
       'union, complex real128')
  n = n + 1

  tests(n) = assert(all(union(arange(1,0),arange(1,0)) == arange(1,0)),&
       'union of empty arrays is an empty array')
  n = n + 1

  tests(n) = assert(all(union([1,2,2,3],arange(1,0)) == set([1,2,2,3])),&
       'union(x,[]) == set(x)')
  n = n + 1

  tests(n) = assert(all(union(arange(1,0),[1,2,2,3]) == set([1,2,2,3])),&
       'union([],x) == set(x)')
  n = n + 1

  tests(n) = assert(all(union([1,2],[3,4]) == ([1,2].union.[3,4])),&
       'union operator, a.union.b')
  n = n + 1

  test_failed = .false.
  call report_tests(tests,test_failed)
  isgood = isgood .and. (.not. test_failed)

  write(*,*) ''
  if (isgood) then
     write(*,*) 'mo_functional ', color('o.k.', c_green)
  else
     write(*,*) 'mo_functional ', color('failed!', c_red)
  endif

end program main

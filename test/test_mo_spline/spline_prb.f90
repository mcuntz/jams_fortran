program main

  !*****************************************************************************80
  !
  !! MAIN is the main program for SPLINE_PRB.
  !
  !  Discussion:
  !
  !    SPLINE_PRB calls the SPLINE tests.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    06 June 2013
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline

  implicit none

  logical :: isok
  character(len=100) :: line1, line2, line11, line22

  !  call timestamp ( )

  open(unit=30, file="../FORTRAN_chs_lib/test/test_mo_spline/spline_prb_output.txt",  status="replace", recl=100)

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'SPLINE_PRB'
  write(30, '(a)' ) '  FORTRAN90 version:'
  write(30, '(a)' ) '  Test the SPLINE library.'

  call test001 ( )
  call test002 ( )
  call test003 ( )
  call test004 ( )
  call test005 ( )
  call test006 ( )

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test115 ( )
  call test116 ( )
  call test12 ( )
  call test125 ( )
  call test126 ( )
  call test127 ( )
  !call test13 ( )
  call test131 ( )
  call test14 ( )
  call test143 ( )
  call test144 ( )
  call test145 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )

  call test20 ( )
  call test205 ( )
  call test21 ( )
  call test215 ( )
  call test22 ( )
  call test225 ( )
  call test23 ( )
  call test235 ( )
  call test24 ( )
  !
  !  Terminate.
  !
  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'SPLINE_PRB'
  write(30, '(a)' ) '  Normal end of execution.'

  write(30, '(a)' ) ' '
  !  call timestamp ( )

  close(30)

  ! Check against standard output
  open(unit=30, file="../FORTRAN_chs_lib/test/test_mo_spline/spline_prb_std_output.txt",  action="read", status="old")
  open(unit=31, file="../FORTRAN_chs_lib/test/test_mo_spline/spline_prb_output.txt",      action="read", status="old")

  ! If you do diff, you see small differences in the last digits with some compiler,
  ! i.e. different between gfortran and nag. -> Compare only the first six charachters.
  isok = .true.
  do
     read(30,*,end=99) line1, line11
     read(31,*,end=99) line2, line22
     if ((trim(line1) /= trim(line2)) .or. (trim(line11(1:6)) /= trim(line22(1:6))))then
        write(*,*) 'Line 1: ', trim(line1), ' ', trim(line11)
        write(*,*) 'Line 2: ', trim(line2), ' ', trim(line22)
        isok = .false.
     endif
  end do

99 continue
  close(30)
  close(31)

  if (isok) then
     write(*,*) 'mo_spline o.k.'
  else
     write(*,*) 'mo_spline failed!'
  endif

end program main

subroutine test001 ( )

  !*****************************************************************************80
  !
  !! TEST001 tests PARABOLA_VAL2.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: dim_num = 1
  integer(i4), parameter :: ndata = 5

  integer(i4) i
  integer(i4) left
  real(dp) xdata(ndata)
  real(dp) xval
  real(dp) ydata(dim_num,ndata)
  real(dp) yval(dim_num)
  real(dp) zdata(ndata)
  real(dp) zval(dim_num)

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST001'
  write(30, '(a)' ) '  PARABOLA_VAL2 evaluates parabolas through'
  write(30, '(a)' ) '    3 points in a table'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  Our data tables will actually be parabolas:'
  write(30, '(a)' ) '    Y: 2*x**2 + 3 * x + 1.'
  write(30, '(a)' ) '    Z: 4*x**2 - 2 * x + 5.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '         I        X         Y           Z'
  write(30, '(a)' ) ' '

  do i = 1, ndata

     xval = 2.0D+00 * real ( i, kind=dp )
     xdata(i) = xval
     ydata(1,i) = 2.0D+00 * xval * xval + 3.0 * xval + 1.0D+00
     zdata(i) = 4.0D+00 * xval * xval - 2.0D+00 * xval + 5.0D+00
     write(30, '(2x,i8,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
          i, xdata(i), ydata(1,i), zdata(i)

  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  Interpolated data:'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '      LEFT        X         Y           Z'
  write(30, '(a)' ) ' '

  do i = 1, 5

     xval = real ( 2 * i - 1, kind=dp )
     left = i

     if ( ndata - 2 < left ) then
        left = ndata - 2
     end if

     if ( left < 1 ) then
        left = 1
     end if

     call parabola_val2 ( dim_num, ndata, xdata, ydata, left, xval, yval )

     call parabola_val2 ( dim_num, ndata, xdata, zdata, left, xval, zval )

     write(30, '(2x,i8,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
          left, xval, yval(1), zval(1)

  end do

  return
end subroutine test001
subroutine test002 ( )

  !*****************************************************************************80
  !
  !! TEST002 tests R8VEC_BRACKET.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 10
  integer(i4), parameter :: test_num = 6

  integer(i4) left
  integer(i4) right
  integer(i4) test
  real(dp) x(n)
  real (  kind=dp ), dimension ( test_num ) :: xtest = (/ &
       -10.0D+00, 1.0D+00, 4.5D+00, 5.0D+00, 10.0D+00, 12.0D+00 /)
  real(dp) xval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST002'
  write(30, '(a)' ) '  R8VEC_BRACKET finds a pair of entries in a'
  write(30, '(a)' ) '    sorted real array which bracket a value.'

  call r8vec_indicator ( n, x )
  x(6) = x(5)

  !call r8vec_print ( n, x, '  Sorted array:' )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    LEFT             RIGHT'
  write(30, '(a)' ) '  X(LEFT)   XVAL   X(RIGHT)'
  write(30, '(a)' ) ' '

  do test = 1, test_num

     xval = xtest(test)

     call r8vec_bracket ( n, x, xval, left, right )

     write(30, '(2x,i14,14x,i14)' ) left, right

     if ( 1 <= left .and. 1 <= right ) then
        write(30, '(2x,3g14.6)' ) x(left), xval, x(right)
     else if ( left < 1 .and. 1 <= right ) then
        write(30, '(2x,14x,2g14.6)' )          xval, x(right)
     else if ( 1 <= left .and. right < 1 ) then
        write(30, '(2x,2g14.6)' ) x(left), xval
     else if ( left < 1 .and. right < 1 ) then
        write(30, '(2x,14x,g14.6)' )          xval
     end if

  end do

  return
end subroutine test002
subroutine test003 ( )

  !*****************************************************************************80
  !
  !! TEST003 tests R8VEC_BRACKET3.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 10
  integer(i4), parameter :: test_num = 6

  integer(i4) left
  integer(i4) test
  real(dp) x(n)
  real (  kind=dp ), dimension ( test_num ) :: xtest = (/ &
       -10.0D+00, 1.0D+00, 4.5D+00, 5.0D+00, 10.0D+00, 12.0D+00 /)
  real(dp) xval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST003'
  write(30, '(a)' ) '  R8VEC_BRACKET3 finds a pair of entries in a'
  write(30, '(a)' ) '    sorted real array which bracket a value.'

  call r8vec_indicator ( n, x )
  x(6) = x(5)

  !call r8vec_print ( n, x, '  Sorted array:' )

  left = ( n + 1 ) / 2

  do test = 1, test_num

     xval = xtest(test)

     write(30, '(a)' ) ' '
     write(30, '(a,g14.6)' ) '  Search for XVAL = ', xval

     write(30, '(a,i8)' ) '  Starting guess for interval is = ', left

     call r8vec_bracket3 ( n, x, xval, left )

     write(30, '(a)' ) '  Nearest interval:'
     write(30, '(2x,a,i8,a,g14.6)' ) '    X[', left,' ]= ', x(left)
     write(30, '(2x,a,i8,a,g14.6)' ) '    X[', left+1, ' ]= ', x(left+1)

  end do

  return
end subroutine test003
subroutine test004 ( )

  !*****************************************************************************80
  !
  !! TEST004 tests R8VEC_ORDER_TYPE.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 4
  integer(i4), parameter :: test_num = 6

  integer(i4) j
  integer(i4) order
  integer(i4) test
  real(dp), dimension(n,test_num) :: x

  x(1,1) = 1.0D+00
  x(2,1) = 3.0D+00
  x(3,1) = 2.0D+00
  x(4,1) = 4.0D+00

  x(1,2) = 2.0D+00
  x(2,2) = 2.0D+00
  x(3,2) = 2.0D+00
  x(4,2) = 2.0D+00

  x(1,3) = 1.0D+00
  x(2,3) = 2.0D+00
  x(3,3) = 2.0D+00
  x(4,3) = 4.0D+00

  x(1,4) = 1.0D+00
  x(2,4) = 2.0D+00
  x(3,4) = 3.0D+00
  x(4,4) = 4.0D+00

  x(1,5) = 4.0D+00
  x(2,5) = 4.0D+00
  x(3,5) = 3.0D+00
  x(4,5) = 1.0D+00

  x(1,6) = 9.0D+00
  x(2,6) = 7.0D+00
  x(3,6) = 3.0D+00
  x(4,6) = 0.0D+00

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST004'
  write(30, '(a)' ) '  R8VEC_ORDER_TYPE classifies a real vector as'
  write(30, '(a)' ) '  -1: no order'
  write(30, '(a)' ) '   0: all equal;'
  write(30, '(a)' ) '   1: ascending;'
  write(30, '(a)' ) '   2: strictly ascending;'
  write(30, '(a)' ) '   3: descending;'
  write(30, '(a)' ) '   4: strictly descending.'
  write(30, '(a)' ) ' '

  do test = 1, test_num

     call r8vec_order_type ( n, x(1,test), order )

     write(30, '(a)' ) ' '
     write(30, '(a,i8)' ) '  The following vector has order type ', order
     write(30, '(a)' ) ' '
     do j = 1, n
        write(30, '(2x,i8,g14.6)' ) j, x(j,test)
     end do

  end do

  return
end subroutine test004
subroutine test005 ( )

  !*****************************************************************************80
  !
  !! TEST005 tests R83_NP_FS.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 10

  real(dp) a(3,n)
  real(dp) b(n)
  integer(i4) :: seed = 123456789
  real(dp) x(n)

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST005'
  write(30, '(a)' ) '  R83_NP_FS factors and solves a tridiagonal'
  write(30, '(a)' ) '    linear system.'
  write(30, '(a)' ) ' '
  write(30, '(a,i8)' ) '  Matrix order N = ', n
  !
  !  Set the matrix elements.
  !
  call r83_uniform ( n, seed, a )
  !
  !  Set the desired solution.
  !
  call r8vec_indicator ( n, x )
  !
  !  Compute b = A * x.
  !
  call r83_mxv ( n, a, x, b )
  !
  !  Wipe out the solution.
  !
  x(1:n) = 0.0D+00
  !
  !  Solve the system.
  !
  call r83_np_fs ( n, a, b, x )

  !call r8vec_print ( n, x, '  Solution:' )

  return
end subroutine test005
subroutine test006 ( )

  !*****************************************************************************80
  !
  !! TEST006 tests DATA_TO_DIF and DIF_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: maxtab = 8

  real(dp) diftab(maxtab)
  real(dp) err
  real(dp) exact
  integer(i4) j
  integer(i4) ntab
  real(dp) xtab(maxtab)
  real(dp) xval
  real(dp) ytab(maxtab)
  real(dp) yval

  xval = 2.5D+00
  exact = exp ( xval )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST006'
  write(30, '(a,i8)' ) '  Approximate Y = EXP(X) using orders 1 to ', maxtab
  write(30, '(a,g14.6)' ) '  Evaluate at X = ', xval
  write(30, '(a,g14.6)' ) '  where EXP(X)=   ', exact

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    Order  Approximate Y     Error'
  write(30, '(a)' ) ' '

  do ntab = 1, maxtab

     do j = 1, ntab
        xtab(j) = real ( j - 1, kind=dp )
        ytab(j) = exp ( xtab(j) )
     end do

     call data_to_dif ( ntab, xtab, ytab, diftab )

     call dif_val ( ntab, xtab, diftab, xval, yval )

     err = yval - exact
     write(30, ' ( 2x, i8, 2g14.6 )' ) ntab, yval, err

  end do

  return
end subroutine test006
subroutine test01 ( )

  !*****************************************************************************80
  !
  !! TEST01 tests BASIS_FUNCTION_B_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: ndata = 5

  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  character mark
  integer(i4), parameter :: nsample = 4
  real(dp), dimension ( ndata ) :: tdata = (/ &
       0.0D+00, 1.0D+00, 4.0D+00, 6.0D+00, 10.0D+00 /)
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST01'
  write(30, '(a)' ) '  BASIS_FUNCTION_B_VAL evaluates the '
  write(30, '(a)' ) '    B spline basis function.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    T            B(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_function_b_val ( tdata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test01
subroutine test02 ( )

  !*****************************************************************************80
  !
  !! TEST02 tests BASIS_FUNCTION_BETA_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: ndata = 5

  real(dp) beta1
  real(dp) beta2
  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  character mark
  integer(i4), parameter :: nsample = 4
  real(dp), dimension ( ndata ) :: tdata = (/ &
       0.0D+00, 1.0D+00, 4.0D+00, 6.0D+00, 10.0D+00 /)
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST02'
  write(30, '(a)' ) '  BASIS_FUNCTION_BETA_VAL evaluates the '
  write(30, '(a)' ) '    Beta spline basis function.'

  beta1 = 1.0D+00
  beta2 = 0.0D+00

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  BETA1 = ', beta1
  write(30, '(a,g14.6)' ) '  BETA2 = ', beta2
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       T            B(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_function_beta_val ( beta1, beta2, tdata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  beta1 = 1.0D+00
  beta2 = 100.0D+00

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  BETA1 = ', beta1
  write(30, '(a,g14.6)' ) '  BETA2 = ', beta2
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       T            B(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_function_beta_val ( beta1, beta2, tdata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  beta1 = 100.0D+00
  beta2 = 0.0D+00

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  BETA1 = ', beta1
  write(30, '(a,g14.6)' ) '  BETA2 = ', beta2
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       T            B(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_function_beta_val ( beta1, beta2, tdata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test02
subroutine test03 ( )

  !*****************************************************************************80
  !
  !! TEST03 tests BASIS_MATRIX_B_UNI and BASIS_MATRIX_TMP.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 4
  integer(i4), parameter :: ndata = 4

  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  integer(i4) left
  character mark
  real(dp) mbasis(n,n)
  integer(i4), parameter :: nsample = 4
  real(dp), dimension (ndata) ::  tdata = (/ &
       -1.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /)
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp), dimension (ndata) ::  ydata = (/ &
       4.0D+00, 7.0D+00, 12.0D+00, 19.0D+00 /)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST03'
  write(30, '(a)' ) '  BASIS_MATRIX_B_UNI sets up the basis matrix'
  write(30, '(a)' ) '    for the uniform B spline.'

  call basis_matrix_b_uni ( mbasis )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       TDATA         YDATA'
  write(30, '(a)' ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '        T            Spline(T)'
  write(30, '(a)' ) ' '

  left = 2

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test03
subroutine test04 ( )

  !*****************************************************************************80
  !
  !! TEST04 tests BASIS_MATRIX_BETA_UNI and BASIS_MATRIX_TMP.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 4
  integer(i4), parameter :: ndata = 4

  real(dp) beta1
  real(dp) beta2
  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  integer(i4) left
  character mark
  real(dp) mbasis(n,n)
  integer(i4), parameter :: nsample = 4
  real(dp), dimension ( ndata ) :: tdata = (/ &
       -1.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /)
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp), dimension ( ndata ) :: ydata = (/ &
       4.0D+00, 7.0D+00, 12.0D+00, 19.0D+00 /)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST04'
  write(30, '(a)' ) '  BASIS_MATRIX_BETA_UNI sets up the basis matrix'
  write(30, '(a)' ) '    for the uniform beta spline.'
  !
  !  First test
  !
  beta1 = 1.0D+00
  beta2 = 0.0D+00

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  BETA1 = ', beta1
  write(30, '(a,g14.6)' ) '  BETA2 = ', beta2

  call basis_matrix_beta_uni ( beta1, beta2, mbasis )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    TDATA, YDATA'
  write(30, '(a)' ) ' '
  do i = 1, ndata
     write(30, '(2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    T, Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do
  !
  !  Second test
  !
  beta1 = 1.0D+00
  beta2 = 100.0D+00

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  BETA1 = ', beta1
  write(30, '(a,g14.6)' ) '  BETA2 = ', beta2

  call basis_matrix_beta_uni ( beta1, beta2, mbasis )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    TDATA, YDATA'
  write(30, '(a)' ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    T, Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do
  !
  !  Third test
  !
  beta1 = 100.0D+00
  beta2 = 0.0D+00

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  BETA1 = ', beta1
  write(30, '(a,g14.6)' ) '  BETA2 = ', beta2

  call basis_matrix_beta_uni ( beta1, beta2, mbasis )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '     TDATA        YDATA'
  write(30, '(a)' ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '        T           Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp  ) * thi ) &
             / real ( nsample,     kind=dp  )

        call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test04
subroutine test05 ( )

  !*****************************************************************************80
  !
  !! TEST05 tests BASIS_MATRIX_BEZIER and BASIS_MATRIX_TMP.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 4
  integer(i4), parameter :: ndata = 4

  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  integer(i4) left
  character mark
  real(dp) mbasis(n,n)
  integer(i4), parameter :: nsample = 4
  real(dp), dimension ( ndata ) :: tdata = (/ &
       0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00 /)
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp), dimension ( ndata ) :: ydata = (/ &
       7.0D+00,  8.3333333D+00,   10.0D+00, 12.0D+00 /)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST05'
  write(30, '(a)' ) '  BASIS_MATRIX_BEZIER sets up the basis'
  write(30, '(a)' ) '    matrix for the uniform Bezier spline.'

  call basis_matrix_bezier ( mbasis )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       TDATA          YDATA'
  write(30, '(a)' ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '        T            Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test05
subroutine test06 ( )

  !*****************************************************************************80
  !
  !! TEST06 tests BASIS_MATRIX_HERMITE and BASIS_MATRIX_TMP.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 4
  integer(i4), parameter :: ndata = 4

  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  integer(i4) left
  character mark
  real(dp) mbasis(n,n)
  integer(i4), parameter :: nsample = 4
  real(dp), dimension ( ndata ) :: tdata = (/ &
       0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00 /)
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp), dimension ( ndata ) :: ydata = (/ &
       7.0D+00, 12.0D+00, 4.0D+00, 6.0D+00 /)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST06'
  write(30, '(a)' ) '  BASIS_MATRIX_HERMITE sets up the basis matrix'
  write(30, '(a)' ) '    for the Hermite spline.'

  call basis_matrix_hermite ( mbasis )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       TDATA        YDATA'
  write(30, '(a)' ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '        T           Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test06
subroutine test07 ( )

  !*****************************************************************************80
  !
  !! TEST07 tests BASIS_MATRIX_OVERHAUSER_UNI and BASIS_MATRIX_TMP.
  !
  !  Discussion:
  !
  !   YDATA(1:NDATA) = ( TDATA(1:NDATA) + 2 )**2 + 3
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 4
  integer(i4), parameter :: ndata = 4

  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  integer(i4) left
  character mark
  real(dp) mbasis(n,n)
  integer(i4), parameter :: nsample = 4
  real(dp), dimension ( ndata ) :: tdata = (/ &
       -1.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /)
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp), dimension ( ndata ) :: ydata = (/ &
       4.0D+00, 7.0D+00, 12.0D+00, 19.0D+00 /)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST07'
  write(30, '(a)' ) '  BASIS_MATRIX_OVERHAUSER_UNI sets up the basis'
  write(30, '(a)' ) '    matrix for the uniform Overhauser spline.'

  call basis_matrix_overhauser_uni ( mbasis )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       TDATA         YDATA'
  write(30, '(a)' ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '        T            Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test07
subroutine test08 ( )

  !*****************************************************************************80
  !
  !! TEST08 tests BASIS_MATRIX_OVERHAUSER_NONUNI and BASIS_MATRIX_TMP.
  !
  !  Discussion:
  !
  !    YDATA(1:NDATA) = ( TDATA(1:NDATA) - 2 )**2 + 3
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 4
  integer(i4), parameter :: ndata = 4

  real(dp) alpha
  real(dp) beta
  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  integer(i4) left
  character mark
  real(dp) mbasis(n,n)
  integer(i4), parameter :: nsample = 4
  real(dp), dimension ( ndata ) :: tdata
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp), dimension ( ndata ) :: ydata
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST08'
  write(30, '(a)' ) '  BASIS_MATRIX_OVERHAUSER_NONUNI sets up the'
  write(30, '(a)' ) '    basis matrix for the nonuniform Overhauser'
  write(30, '(a)' ) '    spline.'

  tdata = (/ 0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00 /)

  alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
  beta = ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  ALPHA = ', alpha
  write(30, '(a,g14.6)' ) '  BETA = ', beta

  call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

  do i = 1, ndata
     ydata(i) = ( tdata(i) - 2.0D+00 )**2 + 3.0D+00
  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       TDATA         YDATA'
  write(30, '(a)' ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       T           Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  tdata(1:4) = (/ 0.0D+00, 1.0D+00, 2.0D+00, 5.0D+00 /)

  alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
  beta = ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  ALPHA = ', alpha
  write(30, '(a,g14.6)' ) '  BETA = ', beta

  call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

  do i = 1, ndata
     ydata(i) = ( tdata(i) - 2.0D+00 )**2 + 3.0D+00
  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    TDATA, YDATA'
  write(30, '(a)' ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    T, Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  tdata(1:4) = (/ 0.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /)

  alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
  beta = ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  ALPHA = ', alpha
  write(30, '(a,g14.6)' ) '  BETA = ', beta

  call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

  do i = 1, ndata
     ydata(i) = ( tdata(i) - 2.0D+00 )**2 + 3.0D+00
  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    TDATA, YDATA'
  write(30, '(a)' ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    T, Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test08
subroutine test09 ( )

  !*****************************************************************************80
  !
  !! TEST09 tests BASIS_MATRIX_OVERHAUSER_NONUNI and BASIS_MATRIX_TMP.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 4
  integer(i4), parameter :: ndata = 4

  real(dp) alpha
  real(dp) beta
  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  integer(i4) left
  character mark
  real(dp) mbasis(n,n)
  integer(i4), parameter :: nsample = 4
  real(dp) tdata(ndata)
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp) ydata(ndata)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST09'
  write(30, '(a)' ) '  BASIS_MATRIX_OVERHAUSER_NONUNI sets up the'
  write(30, '(a)' ) '    basis matrix for the nonuniform Overhauser '
  write(30, '(a)' ) '    spline.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  First test that the nonuniform code can match'
  write(30, '(a)' ) '  the uniform code.  Compare these results with'
  write(30, '(a)' ) '  the uniform output.'
  write(30, '(a)' ) ' '

  tdata(1:4) = (/ -1.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /)

  alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
  beta =  ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  ALPHA = ', alpha
  write(30, '(a,g14.6)' ) '  BETA = ', beta

  call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

  do i = 1, ndata
     ydata(i) = ( tdata(i) + 2.0D+00 )**2 + 3.0D+00
  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    TDATA, YDATA'
  write(30, '(a)' ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    T, Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  Now test that the nonuniform code on a'
  write(30, '(a)' ) '  nonuniform grid.'
  write(30, '(a)' ) ' '

  tdata(1:4) = (/ -4.0D+00, -3.0D+00, -1.0D+00, 2.0D+00 /)

  alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
  beta =  ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  ALPHA = ', alpha
  write(30, '(a,g14.6)' ) '  BETA =  ', beta

  call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

  ydata(1:ndata) = ( tdata(1:ndata) + 2.0D+00 )**2 + 3.0D+00

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       TDATA         YDATA'
  write(30, '(a)' ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '        T            Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test09
subroutine test10 ( )

  !*****************************************************************************80
  !
  !! TEST10 tests BC_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 2

  integer(i4) i
  integer(i4), parameter :: nsample = 101
  real(dp) t
  real(dp), dimension (0:n) :: xcon = (/ 0.0D+00, 0.75D+00, 1.0D+00 /)
  real(dp) xval
  real(dp), dimension (0:n) :: ycon = (/ 1.0D+00, 0.0D+00,  1.0D+00 /)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST10'
  write(30, '(a)' ) '  BC_VAL evaluates a general Bezier function.'
  !
  !  One point on the curve should be about (0.75, 0.536).
  !
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       T            X(T)          Y(T)'
  write(30, '(a)' ) ' '

  do i = 1, nsample
     t = real (       i - 1, kind=dp ) &
          / real ( nsample - 1, kind=dp )
     call bc_val ( n, t, xcon, ycon, xval, yval )
     write(30, '(2x,3g14.6)' ) t, xval, yval
  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  The point ( 0.75, 0.536 ) should be on the curve.'

  return
end subroutine test10
subroutine test11 ( )

  !*****************************************************************************80
  !
  !! TEST11 tests BEZ_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 2

  real(dp) :: a = 0.0D+00
  real(dp) :: b = 1.0D+00
  !real(dp) bez_val
  integer(i4) i
  integer(i4) :: nsample = 21
  real(dp) x
  real(dp), dimension ( 0 : n ) :: y = (/ 1.0D+00, 0.0D+00, 1.0D+00 /)

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST11'
  write(30, '(a)' ) '  BEZ_VAL evaluates a Bezier function.'
  !
  !  One point on the curve should be (0.75, 20/32).
  !
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '         T    X(T)          Y(T)'
  write(30, '(a)' ) ' '

  do i = 1, nsample

     x = ( real ( nsample - i,     kind=dp ) * a   &
          + real (           i - 1, kind=dp ) * b ) &
          / real ( nsample     - 1, kind=dp )

     write(30, '(2x,i8,2g14.6)' ) i, x, bez_val ( n, x, a, b, y )

  end do

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  When X = ', 0.75D+00
  write(30, '(a,g14.6)' ) '  BEZ_VAL(X) should be ', 0.625D+00

  return
end subroutine test11
subroutine test115 ( )

  !*****************************************************************************80
  !
  !! TEST115 tests BP01.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n_max = 3

  real(dp) :: a = 0.0D+00
  real(dp) :: b = 1.0D+00
  real(dp) bern(0:n_max)
  integer(i4) i
  integer(i4) n
  integer(i4) :: nsample = 11
  real(dp) x

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST115'
  write(30, '(a)' ) '  BP01 evaluates the Bernstein basis polynomials'
  write(30, '(a)' ) '  for the interval [0,1].'

  do n = 0, n_max

     write(30, '(a)' ) ' '
     write(30, '(a,i8)' ) '  Degree N = ', n
     write(30, '(a)' ) ' '
     write(30, '(a)' ) &
          '   X         BERN(N,0,X)  BERN(N,1,X)  BERN(N,2,X)  BERN(N,3,X)'
     write(30, '(a)' ) ' '

     do i = 1, nsample

        x = ( real ( nsample - i,     kind=dp ) * a   &
             + real (           i - 1, kind=dp ) * b ) &
             / real ( nsample     - 1, kind=dp )

        call bp01 ( n, x, bern )

        write(30, '(2x,f8.4,4x,5g14.6)' ) x, bern(0:n)

     end do

  end do

  return
end subroutine test115
subroutine test116 ( )

  !*****************************************************************************80
  !
  !! TEST116 tests BPAB.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n_max = 3

  real(dp), parameter :: a = 1.0D+00
  real(dp), parameter :: b = 3.0D+00
  real(dp) bern(0:n_max)
  integer(i4) i
  integer(i4) n
  real(dp) x
  integer(i4), parameter :: nsample = 11

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST116'
  write(30, '(a)' ) '  BPAB evaluates the Bernstein basis polynomials'
  write(30, '(a)' ) '  for the interval [A,B].'
  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  A = ', a
  write(30, '(a,g14.6)' ) '  B = ', b

  do n = 0, n_max

     write(30, '(a)' ) ' '
     write(30, '(a,i8)' ) '  Degree N = ', n
     write(30, '(a)' ) ' '
     write(30, '(a)' ) &
          '   X         BERN(N,0,X)  BERN(N,1,X)  BERN(N,2,X)  BERN(N,3,X)'
     write(30, '(a)' ) ' '

     do i = 1, nsample

        x = ( real ( nsample - i,     kind=dp ) * a   &
             + real (           i - 1, kind=dp ) * b ) &
             / real ( nsample     - 1, kind=dp )

        call bpab ( n, a, b, x, bern )

        write(30, '(2x,f8.4,4x,5g14.6)' ) x, bern(0:n)

     end do

  end do

  return
end subroutine test116
subroutine test12 ( )

  !*****************************************************************************80
  !
  !! TEST12 tests BPAB_APPROX.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: maxdata = 10

  real(dp) a
  real(dp) b
  integer(i4) i
  integer(i4) ndata
  integer(i4) nsample
  real(dp) xdata(0:maxdata)
  real(dp) xval
  real(dp) ydata(0:maxdata)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST12'
  write(30, '(a)' ) '  BPAB_APPROX evaluates the Bernstein polynomial'
  write(30, '(a)' ) '  approximant to a function F(X).'

  a = 1.0D+00
  b = 3.0D+00

  do ndata = 0, 9, 3

     do i = 0, ndata

        if ( ndata == 0 ) then
           xdata(i) = 0.5D+00 * ( a + b )
        else
           xdata(i) = ( real ( ndata - i, kind=dp ) * a   &
                + real (         i, kind=dp ) * b ) &
                / real ( ndata,     kind=dp )
        end if

        ydata(i) = sin ( xdata(i) )

     end do

     write(30, '(a)' ) ' '
     write(30, '(a)' ) '       XDATA        YDATA'
     write(30, '(a)' ) ' '
     do i = 0, ndata
        write(30, '(2x,2g14.6)' ) xdata(i), ydata(i)
     end do

     write(30, '(a)' ) ' '
     write(30, '(a,i8)' ) '  Bernstein approximant of degree N = ', ndata
     write(30, '(a)' ) ' '
     write(30, '(a)' ) '       X            F(X)          BERN(X)        ERROR'
     write(30, '(a)' ) ' '

     nsample = 2 * ndata + 1

     do i = 1, nsample

        if ( nsample == 1 ) then
           xval = 0.5D+00 * ( a + b )
        else
           xval = ( real ( nsample - i,     kind=dp ) * a   &
                + real (           i - 1, kind=dp ) * b ) &
                / real ( nsample     - 1, kind=dp )
        end if

        call bpab_approx ( ndata, a, b, ydata, xval, yval )

        write(30, '(2x,4g14.6)' ) xval, sin(xval), yval, yval - sin(xval)

     end do

  end do

  return
end subroutine test12
subroutine test125 ( )

  !*****************************************************************************80
  !
  !! TEST125 tests LEAST_SET_OLD and LEAST_VAL_OLD.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: maxdeg = 6
  integer(i4), parameter :: ntab = 21

  real(dp) b(1:maxdeg)
  real(dp) c(0:maxdeg)
  real(dp) d(2:maxdeg)
  real(dp) eps
  real(dp) error
  integer(i4) i
  integer(i4) ierror
  integer(i4) j
  integer(i4) jhi
  integer(i4) ndeg
  real(dp) ptab(ntab)
  real(dp) xtab(ntab)
  real(dp) xval
  real(dp) ytab(ntab)
  real(dp) ytrue
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST125'
  write(30, '(a)' ) '  LEAST_SET_OLD sets a least squares polynomial,'
  write(30, '(a)' ) '  LEAST_VAL_OLD evaluates it.'

  do i = 1, ntab
     xtab(i) = ( real ( ntab - i,     kind=dp ) * ( -1.0D+00 )   &
          + real (        i - 1, kind=dp ) * ( +1.0D+00 ) ) &
          / real ( ntab     - 1, kind=dp )
     ytab(i) = real ( int ( exp ( xtab(i) ) * 100.0D+00 + 0.5D+00 ) ) / 100.0D+00
  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', ntab
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       X             Y'
  write(30, '(a)'    ) ' '

  do i = 1, ntab
     write(30, '(2x,2g14.6)' ) xtab(i), ytab(i)
  end do

  do ndeg = 1, maxdeg

     write(30, '(a)' ) ' '
     write(30, '(a,i8)' ) '  Using a polynomial of degree: ', ndeg
     write(30, '(a)' ) ' '

     call least_set_old ( ntab, xtab, ytab, ndeg, ptab, b, c, d, eps, ierror )

     write(30, '(a)' ) ' '
     write(30, '(a,g14.6)' ) '  Total approximation error = ', eps
     write(30, '(a)' ) ' '
     write(30, '(a)' ) '      X            F(X)          P(X)         Error'
     write(30, '(a)' ) ' '

     do i = 1, ntab

        if ( i < ntab ) then
           jhi = 2
        else
           jhi = 0
        end if

        do j = 0, jhi

           if ( i < ntab ) then

              xval = ( real ( 3 - j, kind=dp ) * xtab(i)     &
                   + real (     j, kind=dp ) * xtab(i+1) ) &
                   / real ( 3,     kind=dp )

           else

              xval = xtab(i)

           end if

           call least_val_old ( xval, ndeg, b, c, d, yval )

           ytrue = real ( int ( exp ( xval ) * 100.0D+00 + 0.5D+00 ), kind=dp ) &
                / 100.0D+00

           error = yval - ytrue
           write(30, '(2x,5g14.6)' ) xval, yval, ytrue, error
        end do

     end do

  end do

  return
end subroutine test125
subroutine test126 ( )

  !*****************************************************************************80
  !
  !! TEST126 tests LEAST_SET and LEAST_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 December 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: point_num = 21
  integer(i4), parameter :: nterms = 4

  real(dp) b(nterms)
  real(dp) c(nterms)
  real(dp) d(nterms)
  real(dp) f(point_num)
  integer(i4) i
  integer(i4) nterms2
  real(dp) px
  real(dp) w(point_num)
  real(dp) x(point_num)

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST126'
  write(30, '(a)' ) '  LEAST_SET sets a least squares polynomial,'
  write(30, '(a)' ) '  LEAST_VAL evaluates it.'

  w(1:point_num) = 1.0D+00

  do i = 1, point_num
     x(i) = - 1.0D+00 + real ( i - 1, kind=dp ) / 10.0D+00
     f(i) = real ( int ( exp ( x(i) ) * 100.0D+00 + 0.5D+00 ), kind=dp ) &
          / 100.0D+00
  end do

  call least_set ( point_num, x, f, w, nterms, b, c, d )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  X, F(X), P(X), Error'
  write(30, '(a)' ) ' '

  do nterms2 = 1, nterms
     write(30, '(a)' ) ' '
     write(30, '(a,i8)' ) '  Using polynomial order = ', nterms2
     write(30, '(a)' ) ' '
     do i = 1, point_num
        call least_val ( nterms2, b, c, d, x(i), px )
        write(30, '(5g14.6)' ) x(i), f(i), px, px - f(i)
     end do
  end do

  return
end subroutine test126
subroutine test127 ( )

  !*****************************************************************************80
  !
  !! TEST127 tests LEAST_SET and LEAST_VAL2.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    09 December 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: point_num = 21
  integer(i4), parameter :: nterms = 4

  real(dp) b(nterms)
  real(dp) c(nterms)
  real(dp) d(nterms)
  real(dp) f(point_num)
  real(dp) fp(point_num)
  integer(i4) i
  integer(i4) nterms2
  real(dp) px
  real(dp) pxp
  real(dp) w(point_num)
  real(dp) x(point_num)

  w(1:point_num) = 1.0D+00

  do i = 1, point_num
     x(i) = -1.0D+00 + real ( i - 1, kind=dp ) / 10.0D+00
     f(i) = x(i)**2 - x(i) - 6.0D+00
     fp(i) = 2.0D+00 * x(i) - 1.0D+00
  end do

  call least_set ( point_num, x, f, w, nterms, b, c, d )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST127'
  write(30, '(a)' ) '  LEAST_SET sets a least squares polynomial,'
  write(30, '(a)' ) '  LEAST_VAL2 evaluates it.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  X, F(X), P(X), FP(X), PP(X)'
  write(30, '(a)' ) ' '

  do nterms2 = 1, nterms
     write(30, '(a)' ) ' '
     write(30, '(a,i8)' ) '  Using polynomial order = ', nterms2
     write(30, '(a)' ) ' '
     do i = 1, point_num
        call least_val2 ( nterms2, b, c, d, x(i), px, pxp )
        write(30, '(5g14.6)' ) x(i), f(i), px, fp(i), pxp
     end do
  end do

  return
end subroutine test127
subroutine test13 ( )

  !*****************************************************************************80
  !
  !! TEST13 tests SPLINE_B_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: ndata = 11

  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  character mark
  integer(i4), parameter :: nsample = 4
  real(dp), parameter :: pi = 3.141592653589793D+00
  real(dp) tdata(ndata)
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp) ydata(ndata)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST13'
  write(30, '(a)' ) '  SPLINE_B_VAL evaluates the B spline.'

  do i = 1, ndata
     tdata(i) = real ( i - 1, kind=dp )
     ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / real ( ndata - 1, kind=dp ) )
  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', ndata
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y'
  write(30, '(a)'    ) ' '

  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       T           Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call spline_b_val ( ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test13
subroutine test131 ( )

  !*****************************************************************************80
  !
  !! TEST131 tests SPLINE_B.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: ndata = 11

  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  character mark
  integer(i4), parameter :: nsample = 4
  real(dp), parameter :: pi = 3.141592653589793D+00
  real(dp) tdata(ndata)
  real(dp) thi
  real(dp) tlo
  real(dp) tval(nsample+1)
  real(dp) ydata(ndata)
  real(dp) yval(nsample+1)

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST13'
  write(30, '(a)' ) '  SPLINE_B_VAL evaluates the B spline.'

  do i = 1, ndata
     tdata(i) = real ( i - 1, kind=dp )
     ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / real ( ndata - 1, kind=dp ) )
  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', ndata
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y'
  write(30, '(a)'    ) ' '

  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '       T           Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi
        tval(j+1) = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )
     end do
     yval(1:jhi+1) = spline_b(tval(1:jhi+1), tdata, ydata)

     do j = 0, jhi

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval(j+1), yval(j+1)

     end do

  end do

  return
end subroutine test131
subroutine test14 ( )

  !*****************************************************************************80
  !
  !! TEST14 tests SPLINE_BETA_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: ndata = 11

  real(dp) beta1
  real(dp) beta2
  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  character mark
  integer(i4), parameter :: nsample = 4
  real(dp), parameter :: pi = 3.141592653589793D+00
  real(dp) tdata(ndata)
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp) ydata(ndata)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST14'
  write(30, '(a)' ) '  SPLINE_BETA_VAL evaluates the BETA spline.'

  do i = 1, ndata
     tdata(i) = real ( i - 1, kind=dp )
     ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / real ( ndata - 1, kind=dp ) )
  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', ndata
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y'
  write(30, '(a)'    ) ' '

  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  beta1 = 1.0D+00
  beta2 = 0.0D+00

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  BETA1 = ', beta1
  write(30, '(a,g14.6)' ) '  BETA2 = ', beta2
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    T, Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call spline_beta_val ( beta1, beta2, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  beta1 = 1.0D+00
  beta2 = 100.0D+00

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  BETA1 = ', beta1
  write(30, '(a,g14.6)' ) '  BETA2 = ', beta2
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    T, Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call spline_beta_val ( beta1, beta2, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  beta1 = 100.0D+00
  beta2 = 0.0D+00

  write(30, '(a)' ) ' '
  write(30, '(a,g14.6)' ) '  BETA1 = ', beta1
  write(30, '(a,g14.6)' ) '  BETA2 = ', beta2
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    T, Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call spline_beta_val ( beta1, beta2, ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test14
subroutine test143 ( )

  !*****************************************************************************80
  !
  !! TEST143 tests SPLINE_BEZIER_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: dim_num = 1
  integer(i4), parameter :: interval_num = 3

  real(dp) data_val(dim_num,3*interval_num+1)
  real(dp) dxdt
  integer(i4) interval
  integer(i4) j
  character mark
  integer(i4) point
  integer(i4) point_num
  real(dp), allocatable, dimension ( : ) :: point_t
  real(dp), allocatable, dimension ( :, : ) :: point_val
  real(dp) t
  real(dp) t_max
  real(dp) t_min
  real(dp) x
  real(dp), parameter :: x_max = 2.0D+00 * 3.141592653589793D+00
  real(dp), parameter :: x_min = 0.0D+00
  real(dp) xdata(0:interval_num)
  real(dp) ydata(0:interval_num)
  real(dp) ypdata(0:interval_num)

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST143'
  write(30, '(a)' ) '  SPLINE_BEZIER_VAL evaluates a cubic Bezier spline.'
  !
  !  Construct the data.
  !
  do interval = 0, interval_num

     x = ( real ( interval_num - interval, kind=dp ) * x_min   &
          + real (                interval, kind=dp ) * x_max ) &
          / real ( interval_num,            kind=dp )

     xdata(interval) = x
     ydata(interval) =  sin ( x )
     ypdata(interval) = cos ( x )

  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of intervals = ', interval_num
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       X             Y            dYdX'
  write(30, '(a)'    ) ' '

  do interval = 0, interval_num
     write(30, '(2x,3g14.6)' ) &
          xdata(interval), ydata(interval), ypdata(interval)
  end do
  !
  !  Construct the Bezier control data.
  !
  dxdt = ( x_max - x_min ) / real ( interval_num, kind=dp )
  j = 0

  do interval = 1, interval_num

     if ( interval == 1 ) then
        j = j + 1
        data_val(1,j) = ydata(interval-1)
     end if

     j = j + 1
     data_val(1,j) = ydata(interval-1) + ypdata(interval-1) * dxdt / 3.0D+00

     j = j + 1
     data_val(1,j) = ydata(interval)   - ypdata(interval) * dxdt  / 3.0D+00

     j = j + 1
     data_val(1,j) = ydata(interval)

  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  The control points'
  write(30, '(a)' ) '  Interpolation points are marked with a "*".'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '                     T            P(T)'
  write(30, '(a)' ) ' '

  do j = 1, 3 * interval_num + 1

     t = real ( j - 1, kind=dp ) / 3.0

     if ( abs ( nint ( t ) - t ) < 0.00001D+00 ) then
        mark = '*'
     else
        mark = ' '
     end if

     write(30, '(2x,a,2x,i8,2g14.6)' ) mark, j, t, data_val(1,j)
  end do

  point_num = 6 * interval_num + 1
  allocate ( point_t(1:point_num) )
  allocate ( point_val(1:dim_num,1:point_num) )

  t_min = 0.0D+00
  t_max = real ( interval_num, kind=dp )

  do point = 1, point_num

     t = ( real ( point_num - point,     kind=dp ) * t_min   &
          + real (             point - 1, kind=dp ) * t_max ) &
          / real ( point_num         - 1, kind=dp )

     point_t(point) = t

  end do

  call spline_bezier_val ( dim_num, interval_num, data_val, point_num, &
       point_t, point_val )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  The Bezier spline, sampled at various points.'
  write(30, '(a)' ) '  Interpolation points are marked with a "*".'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '                   T        Spline(T)      F(X(T))'
  write(30, '(a)' ) ' '

  do point = 1, point_num

     if ( abs ( nint ( point_t(point) ) - point_t(point) ) < 0.00001D+00 ) then
        mark = '*'
     else
        mark = ' '
     end if

     x = ( real ( point_num - point,     kind=dp ) * x_min   &
          + real (             point - 1, kind=dp ) * x_max ) &
          / real ( point_num         - 1, kind=dp )

     write(30, '(2x,a,2x,i8,2x,f10.6,2g14.6)' ) &
          mark, point, point_t(point), point_val(1,point), sin ( x )

  end do

  deallocate ( point_t )
  deallocate ( point_val )

  return
end subroutine test143
subroutine test144 ( )

  !*****************************************************************************80
  !
  !! TEST144 tests SPLINE_BEZIER_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: dim_num = 1
  integer(i4), parameter :: interval_num = 3

  real(dp) data_val(dim_num,3*interval_num+1)
  integer(i4) interval
  integer(i4) j
  character mark
  integer(i4) point
  integer(i4) point_num
  real(dp), allocatable, dimension ( : ) :: point_t
  real(dp), allocatable, dimension ( :, : ) :: point_val
  real(dp) t
  real(dp) t_max
  real(dp) t_min
  real(dp) x
  real(dp), parameter :: x_max = 2.0D+00 * 3.141592653589793D+00
  real(dp), parameter :: x_min = 0.0D+00
  real(dp) xdata(0:3*interval_num)
  real(dp) ydata(0:3*interval_num)
  real(dp) y0
  real(dp) y1
  real(dp) y2
  real(dp) y3
  !  real(dp) ypdata(0:interval_num)

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST144'
  write(30, '(a)' ) '  SPLINE_BEZIER_VAL evaluates a cubic Bezier spline.'
  write(30, '(a)' ) '  Normally, the "interior" points of a Bezier spline'
  write(30, '(a)' ) '  are not interpolating points.  Instead, the'
  write(30, '(a)' ) '  derivatives at the interval endpoints are used.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  This example shows, however, that it is possible'
  write(30, '(a)' ) '  to start with only data values, and to "massage"'
  write(30, '(a)' ) '  the data so that we can define a cubic Bezier spline'
  write(30, '(a)' ) '  which interpolates ALL the data.'
  !
  !  Construct the data.
  !
  do interval = 0, 3 * interval_num

     x = ( real ( 3 * interval_num - interval, kind=dp ) * x_min   &
          + real (                    interval, kind=dp ) * x_max ) &
          / real ( 3 * interval_num,            kind=dp )

     xdata(interval) = x
     ydata(interval) = sin ( x )

  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of intervals = ', interval_num
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       X             Y'
  write(30, '(a)'    ) ' '

  do interval = 0, 3 * interval_num
     write(30, '(2x,3g14.6)' ) xdata(interval), ydata(interval)
  end do
  !
  !  Construct the Bezier control data.
  !  The middle points must be "massaged".
  !
  j = 0

  do interval = 1, interval_num

     y0 = ydata(3*(interval-1))
     y1 = ydata(3*(interval-1)+1)
     y2 = ydata(3*(interval-1)+2)
     y3 = ydata(3*(interval-1)+3)

     if ( interval == 1 ) then
        j = j + 1
        data_val(1,j) = y0
     end if

     j = j + 1
     data_val(1,j) = &
          ( -5.0D+00 * y0 + 18.0D+00 * y1 - 9.0D+00 * y2 + 2.0D+00 * y3 ) / 6.0D+00

     j = j + 1
     data_val(1,j) = &
          ( 2.0D+00 * y0 - 9.0D+00 * y1 + 18.0D+00 * y2 - 5.0D+00 * y3 ) / 6.0D+00

     j = j + 1
     data_val(1,j) = y3

  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  The control points'
  write(30, '(a)' ) '  ALL control points will be interpolation points!'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '                     T            P(T)'
  write(30, '(a)' ) ' '

  do j = 1, 3 * interval_num + 1

     t = real ( j - 1, kind=dp ) / 3.0

     mark = '*'

     write(30, '(2x,a,2x,i8,2g14.6)' ) mark, j, t, data_val(1,j)
  end do

  point_num = 6 * interval_num + 1
  allocate ( point_t(1:point_num) )
  allocate ( point_val(1:dim_num,1:point_num) )

  t_min = 0.0D+00
  t_max = real ( interval_num, kind=dp )

  do point = 1, point_num

     t = ( real ( point_num - point,     kind=dp ) * t_min   &
          + real (             point - 1, kind=dp ) * t_max ) &
          / real ( point_num         - 1, kind=dp )

     point_t(point) = t

  end do

  call spline_bezier_val ( dim_num, interval_num, data_val, point_num, &
       point_t, point_val )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  The Bezier spline, sampled at various points.'
  write(30, '(a)' ) '  Interpolation points are marked with a "*".'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '                   T        Spline(T)      F(X(T))'
  write(30, '(a)' ) ' '

  do point = 1, point_num

     if ( abs ( nint ( 3 * point_t(point) ) - 3 * point_t(point) ) &
          < 0.00001D+00 ) then
        mark = '*'
     else
        mark = ' '
     end if

     x = ( real ( point_num - point,     kind=dp ) * x_min   &
          + real (             point - 1, kind=dp ) * x_max ) &
          / real ( point_num         - 1, kind=dp )

     write(30, '(2x,a,2x,i8,2x,f10.6,2g14.6)' ) &
          mark, point, point_t(point), point_val(1,point), sin ( x )

  end do

  deallocate ( point_t )
  deallocate ( point_val )

  return
end subroutine test144
subroutine test145 ( )

  !*****************************************************************************80
  !
  !! TEST145 tests SPLINE_CONSTANT_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: ndata = 12
  integer(i4), parameter :: n_test = 20

  real(dp) ahi
  real(dp) alo
  real(dp) frunge
  real(dp) fval
  integer(i4) i
  integer(i4) j
  integer(i4) :: seed = 123456789
  real(dp) tdata(ndata-1)
  real(dp) thi
  real(dp) tlo
  real(dp) t_test(n_test)
  real(dp) tval
  real(dp) ydata(ndata)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST145'
  write(30, '(a)' ) '  SPLINE_CONSTANT_VAL evaluates a piecewise '
  write(30, '(a)' ) '  constant spline.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  Runge''s function, evenly spaced knots.'
  !
  !  Set the data.
  !
  tlo = -1.0D+00
  thi = +1.0D+00
  call r8vec_even ( ndata-1, tlo, thi, tdata )

  do i = 1, ndata

     if ( i == 1 ) then
        tval = tdata(1)
     else if ( i < ndata ) then
        tval = 0.5D+00 * ( tdata(i-1) + tdata(i) )
     else if ( i == ndata ) then
        tval = tdata(i-1)
     end if

     ydata(i) = frunge ( tval )

  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', ndata
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y'
  write(30, '(a)'    ) ' '

  do i = 1, ndata
     write(30, '(2x,a1,14x,g14.6)' ) '*', ydata(i)
     if ( i < ndata ) then
        write(30, '(2x,a1, g14.6)' ) '*', tdata(i)
     end if
  end do
  !
  !  Sample the spline.
  !
  write(30, * ) 'DEBUG: TLO = ', tlo
  write(30, * ) 'DEBUG: THI = ', thi

  alo = tlo - 1.0D+00
  ahi = thi + 1.0D+00

  call r8vec_uniform_01 ( n_test, seed, t_test )

  t_test(1:n_test) = ( 1.0D+00 - t_test(1:n_test) ) * alo &
       +             t_test(1:n_test)   * ahi

  call r8vec_sort_bubble_a ( n_test, t_test )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '     T     Y(interp)    Y(exact)'
  write(30, '(a)' ) ' '

  j = 0
  write(30, '(2x,a1,14x,g14.6)' ) '*', ydata(j+1)
  j = j + 1

  do i = 1, n_test

     tval = t_test(i)

     call spline_constant_val ( ndata, tdata, ydata, tval, yval )

     if ( j <= ndata - 1 ) then
        do while ( tdata(j) <= tval )
           fval = frunge ( tdata(j) )
           write(30, '(2x,a1,g14.6,14x,g14.6)' ) '*', tdata(j), fval
           write(30, '(2x,a1,14x,g14.6)' ) '*', ydata(j+1)
           j = j + 1
           if ( ndata <= j ) then
              exit
           end if
        end do
     end if

     fval = frunge ( tval )

     write(30, '(2x,a1,3g14.6)' ) ' ', tval, yval, fval

  end do

  return
end subroutine test145
subroutine test15 ( )

  !*****************************************************************************80
  !
  !! TEST15 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 June 2013
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 11

  real(dp) frunge
  real(dp) fprunge
  real(dp) fpprunge
  integer(i4) i
  integer(i4) ibcbeg
  integer(i4) ibcend
  integer(i4) j
  integer(i4) jhi
  integer(i4) k
  real(dp) t(n)
  real(dp) tval
  real(dp) y(n)
  real(dp) ybcbeg
  real(dp) ybcend
  real(dp) ypp(n)
  real(dp) yppval
  real(dp) ypval
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST15'
  write(30, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
  write(30, '(a)' ) '  SPLINE_CUBIC_VAL evaluates it.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  Runge''s function, evenly spaced knots.'

  do i = 1, n

     t(i) = ( real ( n - i,     kind=dp ) * (-1.0D+00)   &
          + real (     i - 1, kind=dp ) * (+1.0D+00) ) &
          / real ( n     - 1, kind=dp )

     y(i) =  frunge ( t(i) )

  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', n
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y'
  write(30, '(a)'    ) ' '

  do i = 1, n
     write(30, '(2x,g14.6,2x,g14.6)' ) t(i), y(i)
  end do
  !
  !  Try various boundary conditions.
  !
  do k = 0, 4

     if ( k == 0 ) then

        ibcbeg = 0
        ybcbeg = 0.0D+00

        ibcend = 0
        ybcend = 0.0D+00

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 0 at both ends:'
        write(30, '(a)' ) '  Spline is quadratic in boundary intervals.'

     else if ( k == 1 ) then

        ibcbeg = 1
        ybcbeg = fprunge ( t(1) )

        ibcend = 1
        ybcend = fprunge ( t(n) )

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 1 at both ends:'
        write(30, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg
        write(30, '(a,g14.6)' ) '  Y''(right) = ', ybcend

     else if ( k == 2 ) then

        ibcbeg = 2
        ybcbeg = fpprunge ( t(1) )

        ibcend = 2
        ybcend = fpprunge ( t(n) )

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 2 at both ends:'
        write(30, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
        write(30, '(a,g14.6)' ) '  YP"(right) = ', ybcend

     else if ( k == 3 ) then

        ibcbeg = 2
        ybcbeg = 0.0D+00

        ibcend = 2
        ybcend = 0.0D+00

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  "Natural" spline:'
        write(30, '(a)' ) '  Boundary condition 2 at both ends:'
        write(30, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
        write(30, '(a,g14.6)' ) '  YP"(right) = ', ybcend

     else if ( k == 4 ) then

        ibcbeg = 3
        ibcend = 3

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  "Not-a-knot" spline:'

     end if

     call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

     if ( k == 3 ) then
        write(30, '(a)' ) ' '
        write(30, '(a)' ) '       I      Y(I)          YPP(I)'
        write(30, '(a)' ) ' '
        do i = 1, n
           write(30, '(2x,i8,2g14.6)' ) i, y(i), ypp(i)
        end do
     end if

     write(30, '(a)' ) ' '
     write(30, '(a)' ) '    SPLINE"(T)      F"(T):'
     write(30, '(a)' ) ' '
     do i = 1, n
        write(30, '(2x,2g14.6)' ) ypp(i), fpprunge ( t(i) )
     end do

     write(30, '(a)' ) ' '
     write(30, '(a)' ) '      T            SPLINE(T)      F(T)'
     write(30, '(a)' ) ' '

     do i = 0, n

        if ( i == 0 ) then
           jhi = 1
        else if ( i < n ) then
           jhi = 2
        else
           jhi = 2
        end if

        do j = 1, jhi

           if ( i == 0 ) then
              tval = t(1) - 1.0D+00
           else if ( i < n ) then
              tval = ( real ( jhi - j + 1, kind=dp ) * t(i)     &
                   + real (       j - 1, kind=dp ) * t(i+1) ) &
                   / real ( jhi,         kind=dp )
           else
              if ( j == 1 ) then
                 tval = t(n)
              else
                 tval = t(n) + 1.0D+00
              end if
           end if

           call spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )

           write(30, '(2x,3g14.6)' ) tval, yval, frunge ( tval )

        end do
     end do

  end do

  return
end subroutine test15
subroutine test16 ( )

  !*****************************************************************************80
  !
  !! TEST16 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL2.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 11

  real(dp) frunge
  real(dp) fprunge
  real(dp) fpprunge
  integer(i4) i
  integer(i4) ibcbeg
  integer(i4) ibcend
  integer(i4) j
  integer(i4) jhi
  integer(i4) k
  integer(i4) left
  integer(i4) left_in
  real(dp) t(n)
  real(dp) tval
  real(dp) y(n)
  real(dp) ybcbeg
  real(dp) ybcend
  real(dp) ypp(n)
  real(dp) yppval
  real(dp) ypval
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST16'
  write(30, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
  write(30, '(a)' ) '  SPLINE_CUBIC_VAL2 evaluates it.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  Runge''s function, evenly spaced knots.'
  do i = 1, n

     t(i) = ( real ( n - i,     kind=dp ) * (-1.0D+00)   &
          + real (     i - 1, kind=dp ) * (+1.0D+00) ) &
          / real ( n     - 1, kind=dp )

     y(i) =  frunge ( t(i) )

  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', n
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y'
  write(30, '(a)'    ) ' '

  do i = 1, n
     write(30, '(2x,2g14.6)' ) t(i), y(i)
  end do
  !
  !  Try boundary condition types 0, 1 and 2.
  !
  do k = 0, 2

     if ( k == 0 ) then

        ibcbeg = 0
        ybcbeg = 0.0D+00

        ibcend = 0
        ybcend = 0.0D+00

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 0 at both ends:'
        write(30, '(a)' ) '  Spline is quadratic in boundary intervals.'

     else if ( k == 1 ) then

        ibcbeg = 1
        ybcbeg = fprunge ( t(1) )

        ibcend = 1
        ybcend = fprunge ( t(n) )

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 1 at both ends:'
        write(30, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg
        write(30, '(a,g14.6)' ) '  Y''(right) = ', ybcend

     else if ( k == 2 ) then

        ibcbeg = 2
        ybcbeg = fpprunge ( t(1) )

        ibcend = 2
        ybcend = fpprunge ( t(n) )

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 2 at both ends:'
        write(30, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
        write(30, '(a,g14.6)' ) '  YP"(right) = ', ybcend

     end if

     call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

     write(30, '(a)' ) ' '
     write(30, '(a)' ) '  SPLINE"(T)      F"(T)'
     write(30, '(a)' ) ' '
     do i = 1, n
        write(30, '(2x,2g14.6)' ) ypp(i), fpprunge(t(i))
     end do

     left = 0

     write(30, '(a)' ) ' '
     write(30, '(a)' ) &
          '      T             SPLINE(T)       F(T)       LEFT_IN  LEFT_OUT'
     write(30, '(a)' ) ' '

     do i = 0, n

        if ( i == 0 ) then
           jhi = 1
        else if ( i < n ) then
           jhi = 2
        else
           jhi = 2
        end if

        do j = 1, jhi

           if ( i == 0 ) then
              tval = t(1) - 1.0D+00
           else if ( i < n ) then
              tval = ( real ( jhi - j + 1, kind=dp ) * t(i)     &
                   + real (       j - 1, kind=dp ) * t(i+1) ) &
                   / real ( jhi,         kind=dp )
           else
              if ( j == 1 ) then
                 tval = t(n)
              else
                 tval = t(n) + 1.0D+00
              end if
           end if

           left_in = left

           call spline_cubic_val2 ( n, t, y, ypp, left, tval, yval, ypval, yppval )

           write(30, '(2x,3g14.6,2i8)' ) tval, yval, frunge ( tval ), &
                left_in, left

        end do
     end do

  end do

  return
end subroutine test16
subroutine test17 ( )

  !*****************************************************************************80
  !
  !! TEST17 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
  !
  !  Discussion:
  !
  !    For boundary condition 0, the spline should come very close within
  !    the interpolation interval.
  !
  !    For conditions 1 and 2, the spline should be essentially exactly equal
  !    to the data, inside and outside the interpolation interval.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 11

  real(dp) fcube
  real(dp) fpcube
  real(dp) fppcube
  integer(i4) i
  integer(i4) ibcbeg
  integer(i4) ibcend
  integer(i4) j
  integer(i4) jhi
  integer(i4) k
  real(dp) t(n)
  real(dp) tval
  real(dp) y(n)
  real(dp) ybcbeg
  real(dp) ybcend
  real(dp) ypp(n)
  real(dp) yppval
  real(dp) ypval
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST17'
  write(30, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
  write(30, '(a)' ) '  SPLINE_CUBIC_VAL evaluates it.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  Cubic data, unevenly spaced knots.'

  do i = 1, n
     t(i) = ( real ( i - 1, kind=dp ) &
          / real ( n - 1, kind=dp ) )**2
  end do

  do i = 1, n
     y(i) = fcube ( t(i) )
  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', n
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y'
  write(30, '(a)'    ) ' '

  do i = 1, n
     write(30, '(2x,g14.6,2x,g14.6)' ) t(i), y(i)
  end do
  !
  !  Try boundary condition types 0, 1 and 2.
  !
  do k = 0, 2

     if ( k == 0 ) then

        ibcbeg = 0
        ybcbeg = 0.0D+00

        ibcend = 0
        ybcend = 0.0D+00

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 0 at both ends:'
        write(30, '(a)' ) '  Spline is quadratic in boundary intervals.'

     else if ( k == 1 ) then

        ibcbeg = 1
        ybcbeg = fpcube ( t(1) )

        ibcend = 1
        ybcend = fpcube ( t(n) )

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 1 at both ends:'
        write(30, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg
        write(30, '(a,g14.6)' ) '  Y''(right) = ', ybcend

     else if ( k == 2 ) then

        ibcbeg = 2
        ybcbeg = fppcube ( t(1) )

        ibcend = 2
        ybcend = fppcube ( t(n) )

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 2 at both ends:'
        write(30, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
        write(30, '(a,g14.6)' ) '  YP"(right) = ', ybcend

     end if

     call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

     write(30, '(a)' ) ' '
     write(30, '(a)' ) '       SPLINE"(T)    F"(T):'
     write(30, '(a)' ) ' '
     do i = 1, n
        write(30, '(2x,2g14.6)' ) ypp(i), fppcube(t(i))
     end do

     write(30, '(a)' ) ' '
     write(30, '(a)' ) '       T           SPLINE(T)      F(T)'
     write(30, '(a)' ) ' '

     do i = 0, n

        if ( i == 0 ) then
           jhi = 1
        else if ( i < n ) then
           jhi = 2
        else
           jhi = 2
        end if

        do j = 1, jhi

           if ( i == 0 ) then
              tval = t(1) - 1.0D+00
           else if ( i < n ) then
              tval = ( real ( jhi - j + 1, kind=dp ) * t(i)     &
                   + real (       j - 1, kind=dp ) * t(i+1) ) &
                   / real ( jhi,         kind=dp )
           else
              if ( j == 1 ) then
                 tval = t(n)
              else
                 tval = t(n) + 1.0D+00
              end if
           end if

           call spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )

           write(30, '(2x,3g14.6)' ) tval, yval, fcube ( tval )

        end do
     end do

  end do

  return
end subroutine test17
subroutine test18 ( )

  !*****************************************************************************80
  !
  !! TEST18 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
  !
  !  Discussion:
  !
  !    For boundary condition 0, the spline should come very close within
  !    the interpolation interval.
  !
  !    For conditions 1 and 2, the spline should be essentially exactly equal
  !    to the data, inside and outside the interpolation interval.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 11

  real(dp) fcube
  real(dp) fpcube
  real(dp) fppcube
  integer(i4) i
  integer(i4) ibcbeg
  integer(i4) ibcend
  integer(i4) j
  integer(i4) jhi
  integer(i4) k
  real(dp) t(n)
  real(dp) tval
  real(dp) y(n)
  real(dp) ybcbeg
  real(dp) ybcend
  real(dp) ypp(n)
  real(dp) yppval
  real(dp) ypval
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST18'
  write(30, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
  write(30, '(a)' ) '  SPLINE_CUBIC_VAL evaluates it.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  Cubic data, evenly spaced knots.'
  do i = 1, n
     t(i) = real ( i - 1, kind=dp ) / real ( n - 1, kind=dp )
     y(i) =  fcube ( t(i) )
  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', n
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y'
  write(30, '(a)'    ) ' '

  do i = 1, n
     write(30, '(2x,g14.6,2x,g14.6)' ) t(i), y(i)
  end do
  !
  !  Try boundary condition types 0, 1 and 2.
  !
  do k = 0, 2

     if ( k == 0 ) then

        ibcbeg = 0
        ybcbeg = 0.0D+00

        ibcend = 0
        ybcend = 0.0D+00

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 0 at both ends:'
        write(30, '(a)' ) '  Spline is quadratic in boundary intervals.'

     else if ( k == 1 ) then

        ibcbeg = 1
        ybcbeg = fpcube ( t(1) )

        ibcend = 1
        ybcend = fpcube ( t(n) )

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 1 at both ends:'
        write(30, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg
        write(30, '(a,g14.6)' ) '  Y''(right) = ', ybcend

     else if ( k == 2 ) then

        ibcbeg = 2
        ybcbeg = fppcube ( t(1) )

        ibcend = 2
        ybcend = fppcube ( t(n) )

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 2 at both ends:'
        write(30, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
        write(30, '(a,g14.6)' ) '  YP"(right) = ', ybcend

     end if

     call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

     write(30, '(a)' ) ' '
     write(30, '(a)' ) '   SPLINE"(T)      F"(T):'
     write(30, '(a)' ) ' '
     do i = 1, n
        write(30, '(2x,2g14.6)' ) ypp(i), fppcube(t(i))
     end do

     write(30, '(a)' ) ' '
     write(30, '(a)' ) '        T      SPLINE(T)    F(T)'
     write(30, '(a)' ) ' '

     do i = 0, n

        if ( i == 0 ) then
           jhi = 1
        else if ( i < n ) then
           jhi = 2
        else
           jhi = 2
        end if

        do j = 1, jhi

           if ( i == 0 ) then
              tval = t(1) - 1.0D+00
           else if ( i < n ) then
              tval = ( real ( jhi - j + 1, kind=dp ) * t(i)     &
                   + real (       j - 1, kind=dp ) * t(i+1) ) &
                   / real ( jhi,         kind=dp )
           else
              if ( j == 1 ) then
                 tval = t(n)
              else
                 tval = t(n) + 1.0D+00
              end if
           end if

           call spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )

           write(30, '(2x,f10.4)' ) tval
           write(30, '(2x,10x,2f10.4)' )   yval, fcube ( tval )
           write(30, '(2x,10x,2f10.4)' )   ypval, fpcube ( tval )
           write(30, '(2x,10x,2f10.4)' )   yppval, fppcube ( tval )

        end do
     end do

  end do

  return
end subroutine test18
subroutine test19 ( )

  !*****************************************************************************80
  !
  !! TEST19 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 2

  real(dp) fcube
  real(dp) fpcube
  real(dp) fppcube
  integer(i4) i
  integer(i4) ibcbeg
  integer(i4) ibcend
  integer(i4) j
  integer(i4) jhi
  integer(i4) k1
  integer(i4) k2
  real(dp) t(n)
  real(dp) tval
  real(dp) y(n)
  real(dp) ybcbeg
  real(dp) ybcend
  real(dp) ypp(n)
  real(dp) yppval
  real(dp) ypval
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST19'
  write(30, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
  write(30, '(a)' ) '  SPLINE_CUBIC_VAL evaluates it.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  Cubic data, evenly spaced knots.'
  write(30, '(a)' ) '  ONLY TWO KNOTS!'

  do i = 1, n
     t(i) = real ( i - 1, kind=dp ) / real ( n - 1, kind=dp )
     y(i) = fcube ( t(i) )
  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', n
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y'
  write(30, '(a)'    ) ' '

  do i = 1, n
     write(30, '(2x,g14.6,2x,g14.6)' ) t(i), y(i)
  end do
  !
  !  Try all 9 pairs of boundary condition types 0, 1 and 2.
  !
  do k1 = 0, 2

     if ( k1 == 0 ) then

        ibcbeg = 0
        ybcbeg = 0.0D+00

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 0 at left end.'

     else if ( k1 == 1 ) then

        ibcbeg = 1
        ybcbeg = fpcube ( t(1) )

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 1 at left end.'
        write(30, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg

     else if ( k1 == 2 ) then

        ibcbeg = 2
        ybcbeg = fppcube ( t(1) )

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  Boundary condition 2 at left ends:'
        write(30, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg

     end if

     do k2 = 0, 2

        if ( k2 == 0 ) then

           ibcend = 0
           ybcend = 0.0D+00

           write(30, '(a)' ) ' '
           write(30, '(a)' ) '  Boundary condition 0 at right end.'

        else if ( k2 == 1 ) then

           ibcend = 1
           ybcend = fpcube ( t(n) )

           write(30, '(a)' ) ' '
           write(30, '(a)' ) '  Boundary condition 1 at right end.'
           write(30, '(a,g14.6)' ) '  Y''(right) = ', ybcend

        else if ( k2 == 2 ) then

           ibcend = 2
           ybcend = fppcube ( t(n) )

           write(30, '(a)' ) ' '
           write(30, '(a)' ) '  Boundary condition 2 at right end.'
           write(30, '(a,g14.6)' ) '  YP"(right) = ', ybcend

        end if

        call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  SPLINE"(T), F"(T):'
        write(30, '(a)' ) ' '
        do i = 1, n
           write(30, '(2x,2g14.6)' ) ypp(i), fppcube(t(i))
        end do

        write(30, '(a)' ) ' '
        write(30, '(a)' ) '  T, SPLINE(T), F(T)'
        write(30, '(a)' ) ' '

        do i = 0, n

           if ( i == 0 ) then
              jhi = 1
           else if ( i < n ) then
              jhi = 2
           else
              jhi = 2
           end if

           do j = 1, jhi

              if ( i == 0 ) then
                 tval = t(1) - 1.0D+00
              else if ( i < n ) then
                 tval = ( real ( jhi - j + 1, kind=dp ) * t(i)     &
                      + real (       j - 1, kind=dp ) * t(i+1) ) &
                      / real ( jhi,         kind=dp )
              else
                 if ( j == 1 ) then
                    tval = t(n)
                 else
                    tval = t(n) + 1.0D+00
                 end if
              end if

              call spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )

              write(30, '(2x,f10.4)' ) tval
              write(30, '(2x,10x,2f10.4)' )   yval, fcube ( tval )
              write(30, '(2x,10x,2f10.4)' )   ypval, fpcube ( tval )
              write(30, '(2x,10x,2f10.4)' )   yppval, fppcube ( tval )

           end do
        end do

     end do

  end do

  return
end subroutine test19
subroutine test20 ( )

  !*****************************************************************************80
  !
  !! TEST20 tests SPLINE_HERMITE_SET and SPLINE_HERMITE_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: ndata = 4

  real(dp) c(4,ndata)
  real(dp) fpval
  real(dp) fval
  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  character mark
  real(dp), parameter :: pi = 3.141592653589793D+00
  real(dp) tdata(ndata)
  real(dp) tval
  real(dp) ydata(ndata)
  real(dp) ypdata(ndata)
  real(dp) ypval
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST20'
  write(30, '(a)' ) '  SPLINE_HERMITE_SET sets up a Hermite spline;'
  write(30, '(a)' ) '  SPLINE_HERMITE_VAL evaluates it.'
  !
  !  Set the data.
  !
  do i = 1, ndata
     tdata(i) = ( real ( ndata - i,     kind=dp ) *   0.0D+00          &
          + real (         i - 1, kind=dp ) * ( 0.5D+00 * pi ) ) &
          / real ( ndata     - 1, kind=dp )
     ydata(i) = sin ( tdata(i) )
     ypdata(i) = cos ( tdata(i) )
  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', ndata
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y             Y'''
  write(30, '(a)'    ) ' '
  do i = 1, ndata
     write(30, '(2x,3g14.6)' ) tdata(i), ydata(i), ypdata(i)
  end do
  !
  !  Set up the spline.
  !
  call spline_hermite_set ( ndata, tdata, ydata, ypdata, c )
  !
  !  Now evaluate the spline all over the place.
  !
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T        Y(hermite)     ' // &
       'Y(exact)      Y''(hermite)   Y''(exact)'
  write(30, '(a)'    ) ' '

  do i = 1, ndata

     if ( i == ndata ) then
        jhi = 0
     else
        jhi = 2
     end if

     do j = 0, jhi

        tval = real ( 3 * ( i - 1 ) + j, kind=dp ) * ( 0.5D+00 * pi ) &
             / real ( 3 * ( ndata - 1 ), kind=dp )

        fval = sin ( tval )
        fpval = cos ( tval )

        call spline_hermite_val ( ndata, tdata, c, tval, yval, ypval )

        if (j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,5g14.6)' ) mark, tval, yval, fval, ypval, fpval

     end do

  end do

  return
end subroutine test20
subroutine test205 ( )

  !*****************************************************************************80
  !
  !! TEST205 tests SPLINE_LINEAR_INT and SPLINE_LINEAR_INTSET.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 4

  real(dp) a
  real(dp) b
  real(dp) data_x(n)
  real(dp) data_y(n)
  integer(i4) i
  real(dp) int_x(n+1)
  real(dp) int_v(n)
  real(dp) value

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST205'
  write(30, '(a)' ) '  SPLINE_LINEAR_INTSET is given some interval endpoints,'
  write(30, '(a)' ) '  and a value associated with each interval.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  It determines a linear spline, with breakpoints'
  write(30, '(a)' ) '  at the centers of each interval, whose integral'
  write(30, '(a)' ) '  over each interval is equal to the given value.'

  int_x(1:n+1) = (/ 0.0D+00, 1.0D+00, 4.0D+00, 5.0D+00, 10.0D+00 /)
  int_v(1:n) = (/ 10.0D+00, 2.0D+00, 8.0D+00, 27.5D+00 /)

  !  call r8vec_print ( n+1, int_x, '  The interval end points:' )

  !  call r8vec_print ( n, int_v, '  The desired interval integral values:' )

  call spline_linear_intset ( n, int_x, int_v, data_x, data_y )

  !  call r8vec_print ( n, data_x, '  The spline break points:' )
  !  call r8vec_print ( n, data_y, '  The spline data values: ' )

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  As a check, call SPLINE_LINEAR_INT to compute'
  write(30, '(a)' ) '  the integral of the spline over each interval,'
  write(30, '(a)' ) '  and compare to the desired value.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '     A       B    Desired      Computed'
  write(30, '(a)' ) ' '

  do i = 1, n
     a = int_x(i)
     b = int_x(i+1)
     call spline_linear_int ( n, data_x, data_y, a, b, value )
     write(30, '(2x,2f8.2,2g14.6)' ) a, b, int_v(i), value
  end do

  return
end subroutine test205
subroutine test21 ( )

  !*****************************************************************************80
  !
  !! TEST21 tests SPLINE_LINEAR_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 11

  real(dp) frunge
  real(dp) fval
  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  real(dp) t(n)
  real(dp) tval
  real(dp) y(n)
  real(dp) ypval
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST21'
  write(30, '(a)' ) '  SPLINE_LINEAR_VAL evaluates a linear spline.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  Runge''s function, evenly spaced knots.'

  do i = 1, n

     t(i) = ( real ( n - i,     kind=dp ) * (-1.0D+00)   &
          + real (     i - 1, kind=dp ) * (+1.0D+00) ) &
          / real ( n     - 1, kind=dp )

     y(i) =  frunge ( t(i) )

  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', n
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y'
  write(30, '(a)'    ) ' '
  do i = 1, n
     write(30, '(2x,2g14.6)' ) t(i), y(i)
  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  Interpolation:'
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y            Yexact'
  write(30, '(a)'    ) ' '

  do i = 0, n

     if ( i == 0 ) then
        jhi = 1
     else if ( i < n ) then
        jhi = 2
     else
        jhi = 2
     end if

     do j = 1, jhi

        if ( i == 0 ) then
           tval = t(1) - 1.0D+00
        else if ( i < n ) then
           tval = ( real ( jhi - j + 1, kind=dp ) * t(i)     &
                + real (       j - 1, kind=dp ) * t(i+1) ) &
                / real ( jhi,         kind=dp )
        else
           if ( j == 1 ) then
              tval = t(n)
           else
              tval = t(n) + 1.0D+00
           end if
        end if

        call spline_linear_val ( n, t, y, tval, yval, ypval )

        fval = frunge ( tval )

        write(30, '(2x,3g14.6)' ) tval, yval, fval

     end do

  end do

  return
end subroutine test21
subroutine test215 ( )

  !*****************************************************************************80
  !
  !! TEST215 tests SPLINE_LINEAR_INT.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 3

  real(dp) a
  real(dp) b
  integer(i4) i
  real(dp) int_val
  real(dp), dimension ( n ) :: t = (/ 2.0D+00, 4.5D+00, 7.5D+00 /)
  real(dp), dimension ( n ) :: y = (/ 3.0D+00, 3.75D+00, 5.5D+00 /)

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST215'
  write(30, '(a)' ) '  SPLINE_LINEAR_INT computes the integral '
  write(30, '(a)' ) '  of a linear spline.'
  write(30, '(a)' ) ' '

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', n
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y'
  write(30, '(a)'    ) ' '
  do i = 1, n
     write(30, '(2x,2g14.6)' ) t(i), y(i)
  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '    A             B           Integral'
  write(30, '(a)' ) ' '

  do i = 1, 5

     if ( i == 1 ) then
        a = 0.0D+00
        b = 4.0D+00
     else if ( i == 2 ) then
        a = 4.0D+00
        b = 5.0D+00
     else if ( i == 3 ) then
        a = 5.0D+00
        b = 10.0D+00
     else if ( i == 4 ) then
        a = 0.0D+00
        b = 10.0D+00
     else
        a = 10.0D+00
        b = 0.0D+00
     end if

     call spline_linear_int ( n, t, y, a, b, int_val )

     write(30, '(2x,3g14.6)' ) a, b, int_val

  end do

  return
end subroutine test215
subroutine test22 ( )

  !*****************************************************************************80
  !
  !! TEST22 tests SPLINE_OVERHAUSER_UNI_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: ndata = 11

  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  character mark
  integer(i4), parameter :: nsample = 4
  real(dp), parameter :: pi = 3.141592653589793D+00
  real(dp) tdata(ndata)
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp) ydata(ndata)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST22'
  write(30, '(a)' ) '  SPLINE_OVERHAUSER_UNI_VAL evaluates the'
  write(30, '(a)' ) '    uniform Overhauser spline.'

  do i = 1, ndata
     tdata(i) = real ( i - 1, kind=dp )
     ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / real ( ndata - 1, kind=dp ) )
  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', ndata
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '         T             Y'
  write(30, '(a)'    ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '      T, Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call spline_overhauser_uni_val ( ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test22
subroutine test225 ( )

  !*****************************************************************************80
  !
  !! TEST225 tests SPLINE_OVERHAUSER_NONUNI_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: ndata = 11

  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  character mark
  integer(i4), parameter :: nsample = 4
  real(dp), parameter :: pi = 3.141592653589793D+00
  real(dp) tdata(ndata)
  real(dp) thi
  real(dp) tlo
  real(dp) tval
  real(dp) ydata(ndata)
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST225'
  write(30, '(a)' ) '  SPLINE_OVERHAUSER_NONUNI_VAL evaluates the'
  write(30, '(a)' ) '    nonuniform Overhauser spline.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  In this initial draft of a test, we simply'
  write(30, '(a)' ) '  use uniform nodes.'

  do i = 1, ndata
     tdata(i) = real ( i - 1, kind=dp )
     ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / real ( ndata - 1, kind=dp ) )
  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', ndata
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '         T             Y'
  write(30, '(a)'    ) ' '
  do i = 1, ndata
     write(30, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  write(30, '(a)' ) ' '
  write(30, '(a)' ) '      T, Spline(T)'
  write(30, '(a)' ) ' '

  do i = 0, ndata

     if ( i == 0 ) then
        tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
        thi = tdata(1)
     else if ( i < ndata ) then
        tlo = tdata(i)
        thi = tdata(i+1)
     else if ( ndata <= i ) then
        tlo = tdata(ndata)
        thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
     end if

     if ( i < ndata ) then
        jhi = nsample - 1
     else
        jhi = nsample
     end if

     do j = 0, jhi

        tval = ( real ( nsample - j, kind=dp ) * tlo   &
             + real (           j, kind=dp ) * thi ) &
             / real ( nsample,     kind=dp )

        call spline_overhauser_nonuni_val ( ndata, tdata, ydata, tval, yval )

        if ( 0 < i .and. j == 0 ) then
           mark = '*'
        else
           mark = ' '
        end if

        write(30, '(2x,a1,2g14.6)' ) mark, tval, yval

     end do

  end do

  return
end subroutine test225
subroutine test23 ( )

  !*****************************************************************************80
  !
  !! TEST23 tests SPLINE_OVERHAUSER_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: ndata = 4
  integer(i4), parameter :: dim_num = 2

  integer(i4) i
  real(dp) tdata(ndata)
  real(dp) tval
  real(dp) ydata(dim_num,ndata)
  real(dp) yval(dim_num)

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST23'
  write(30, '(a)' ) '  SPLINE_OVERHAUSER_VAL evaluates the'
  write(30, '(a)' ) '    Overhauser spline.'
  !
  !  Set the data.
  !
  tdata(1) = 1.0D+00
  ydata(1,1) =   0.0D+00
  ydata(2,1) =   0.0D+00

  tdata(2) = 2.0D+00
  ydata(1,2) =   1.0D+00
  ydata(2,2) =   1.0D+00

  tdata(3) = 3.0D+00
  ydata(1,3) =   2.0D+00
  ydata(2,3) = - 1.0D+00

  tdata(4) = 4.0D+00
  ydata(1,4) =   3.0D+00
  ydata(2,4) =   0.0D+00

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', ndata
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y'
  write(30, '(a)'    ) ' '

  do i = 1, ndata
     write(30, '(2x,3g14.6)' ) tdata(i), ydata(1:dim_num,i)
  end do
  !
  !  Now evaluate the spline all over the place.
  !
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  T, Spline value'
  write(30, '(a)' ) ' '

  do i = 0, 6 * ndata + 3

     tval = real ( i, kind=dp ) / 6.0D+00
     call spline_overhauser_val ( dim_num, ndata, tdata, ydata, tval, yval )
     write(30, '(2x,3g14.6)' ) tval, yval(1:dim_num)

  end do

  return
end subroutine test23
subroutine test235 ( )

  !*****************************************************************************80
  !
  !! TEST235 tests SPLINE_PCHIP_SET and SPLINE_PCHIP_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 21
  integer(i4), parameter :: ne = 101

  real(dp) d(n)
  real(dp) diff
  real(dp) f(n)
  !  real(dp) fd(ne)
  real(dp) fe(ne)
  integer(i4) i
  real(dp), external :: frunge
  real(dp) x(n)
  real(dp) xe(ne)

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST235'
  write(30, '(a)' ) '  SPLINE_PCHIP_SET carries out piecewise cubic '
  write(30, '(a)' ) '    Hermite interpolation.'
  write(30, '(a)' ) '  SPLINE_PCHIP_VAL evaluates the interpolant.'
  write(30, '(a)' ) ' '
  !
  !  Compute Runge's function at N points in [-1,1].
  !
  do i = 1, n
     x(i) = -1.0D+00 + real ( i - 1, kind=dp ) / 10.0D+00
     f(i) = frunge ( x(i) )
  end do
  !
  !  SPLINE_PCHIP_SET takes the data in X and F, and constructs a table in D
  !  that defines the interpolant.
  !
  call spline_pchip_set ( n, x, f, d )
  !
  !  Evaluate the interpolant and derivative at NE points from -1 to 0.
  !
  do i = 1, ne
     xe(i) = -1.0D+00 + real ( i - 1, kind=dp ) / real ( ne - 1, kind=dp )
  end do

  call spline_pchip_val ( n, x, f, d, ne, xe, fe )
  !
  !  Print the table of X, F(exact) and F(interpolated)
  !
  do i = 1, ne
     diff = fe(i) - frunge ( xe(i) )
     write(30, '(2x,f8.4,2x,f10.6,2x,f10.6,2x,g14.6)' ) &
          xe(i), frunge ( xe(i) ), fe(i), diff
  end do

  return
end subroutine test235
subroutine test24 ( )

  !*****************************************************************************80
  !
  !! TEST24 tests SPLINE_QUADRATIC_VAL.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    17 June 2006
  !
  !  Author:
  !
  !    John Burkardt
  !
  use mo_kind
  use mo_spline
  implicit none
  integer(i4), parameter :: n = 11

  real(dp) frunge
  real(dp) fval
  integer(i4) i
  integer(i4) j
  integer(i4) jhi
  real(dp) t(n)
  real(dp) tval
  real(dp) y(n)
  real(dp) ypval
  real(dp) yval

  write(30, '(a)' ) ' '
  write(30, '(a)' ) 'TEST24'
  write(30, '(a)' ) '  SPLINE_QUADRATIC_VAL evaluates a '
  write(30, '(a)' ) '    quadratic spline.'
  write(30, '(a)' ) ' '
  write(30, '(a)' ) '  Runge''s function, evenly spaced knots.'

  do i = 1, n

     t(i) = ( real ( n - i,     kind=dp ) * (-1.0D+00)   &
          + real (     i - 1, kind=dp ) * (+1.0D+00) ) &
          / real ( n     - 1, kind=dp )

     y(i) =  frunge ( t(i) )

  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  The data to be interpolated:'
  write(30, '(a)'    ) ' '
  write(30, '(a,i8)' ) '  Number of data values = ', n
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '         T             Y'
  write(30, '(a)'    ) ' '
  do i = 1, n
     write(30, '(2x,2g14.6)' ) t(i), y(i)
  end do

  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '  Interpolated values'
  write(30, '(a)'    ) ' '
  write(30, '(a)'    ) '       T             Y           Y(exact)'
  write(30, '(a)'    ) ' '

  do i = 0, n

     if ( i == 0 ) then
        jhi = 1
     else if ( i < n ) then
        jhi = 2
     else
        jhi = 2
     end if

     do j = 1, jhi

        if ( i == 0 ) then
           tval = t(1) - 1.0D+00
        else if ( i < n ) then
           tval = ( real ( jhi - j + 1, kind=dp ) * t(i)     &
                + real (       j - 1, kind=dp ) * t(i+1) ) &
                / real ( jhi,         kind=dp )
        else
           if ( j == 1 ) then
              tval = t(n)
           else
              tval = t(n) + 1.0D+00
           end if
        end if

        call spline_quadratic_val ( n, t, y, tval, yval, ypval )

        fval = frunge ( tval )

        write(30, '(2x,3g14.6)' ) tval, yval, fval

     end do

  end do

  return
end subroutine test24
function fcube ( x )

  !*****************************************************************************80
  !
  !! FCUBE evaluates a cubic function.
  !
  !  Discussion:
  !
  !    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 February 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real X, the argument.
  !
  !    Output, real FCUBE, the value of the function.
  !
  use mo_kind
  use mo_spline
  implicit none
  real(dp) fcube
  real(dp) x

  fcube = ( ( (       1.0D+00 ) &
       * x + 2.0D+00 ) &
       * x + 3.0D+00 ) &
       * x + 4.0D+00

  return
end function fcube
function fpcube ( x )

  !*****************************************************************************80
  !
  !! FPCUBE evaluates the derivative of a cubic function.
  !
  !  Discussion:
  !
  !    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
  !
  !  Modified:
  !
  !    10 February 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real X, the argument.
  !
  !    Output, real FPCUBE, the value of the derivative of the cubic function.
  !
  use mo_kind
  use mo_spline
  implicit none
  real(dp) fpcube
  real(dp) x

  fpcube = ( 3.0D+00 * x + 4.0D+00 ) * x + 3.0D+00

  return
end function fpcube
function fppcube ( x )

  !*****************************************************************************80
  !
  !! FPPCUBE evaluates the second derivative of a cubic function.
  !
  !  Discussion:
  !
  !    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    10 February 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real X, the argument.
  !
  !    Output, real FPPCUBE, the second derivative of the cubic function.
  !
  use mo_kind
  use mo_spline
  implicit none
  real(dp) fppcube
  real(dp) x

  fppcube = 6.0D+00 * x + 4.0D+00

  return
end function fppcube
function frunge ( x )

  !*****************************************************************************80
  !
  !! FRUNGE evaluates the Runge function.
  !
  !  Discussion:
  !
  !    Interpolation of the Runge function at evenly spaced points in [-1,1]
  !    is a common test.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    24 January 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real X, the argument.
  !
  !    Output, real FRUNGE, the value of the function.
  !
  use mo_kind
  use mo_spline
  implicit none
  real(dp) frunge
  real(dp) x

  frunge = 1.0D+00 / ( 1.0D+00 + 25.0D+00 * x * x )

  return
end function frunge
function fprunge ( x )

  !*****************************************************************************80
  !
  !! FPRUNGE evaluates the derivative of the Runge function.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    24 January 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real X, the argument.
  !
  !    Output, real FPRUNGE, the value of the derivative of the Runge function.
  !
  use mo_kind
  use mo_spline
  implicit none
  real(dp) fprunge
  real(dp) x

  fprunge = - 50.0D+00 * x / ( 1.0D+00 + 25.0D+00 * x * x )**2

  return
end function fprunge
function fpprunge ( x )

  !*****************************************************************************80
  !
  !! FPPRUNGE evaluates the second derivative of the Runge function.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    24 January 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real X, the argument.
  !
  !    Output, real FPPRUNGE, the value of the second derivative of
  !    the Runge function.
  !
  use mo_kind
  use mo_spline
  implicit none
  real(dp) fpprunge
  real(dp) x

  fpprunge = ( - 50.0D+00 + 3750.0D+00 * x * x ) &
       / ( 1.0D+00 + 25.0D+00 * x * x )**3

  return
end function fpprunge

FUNCTION expint_sp(n,x)
    USE mo_kind,            ONLY : i4, sp
    USE mo_nrutil,          ONLY : arth,assert,nrerror
    USE mo_constants,       ONLY : EULER

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: n
    REAL(SP), INTENT(IN) :: x
    REAL(SP) :: expint_sp
    INTEGER(i4), PARAMETER :: MAXIT=100_i4
    REAL(SP), PARAMETER :: EPS=epsilon(x),BIG=huge(x)*EPS
    INTEGER(i4) :: i,nm1
    REAL(SP) :: a,b,c,d,del,fact,h

    call assert(n >= 0_i4, x >= 0.0_sp, (x > 0.0_sp .or. n > 1_i4), &
    'expint_sp args')

    if (n == 0_i4) then
        expint_sp=exp(-x)/x
        RETURN
    end if

    nm1=n-1_i4

    if (x == 0.0_sp) then
        expint_sp=1.0_sp/nm1

    else if (x > 1.0_sp) then
        b=x+n
        c=BIG
        d=1.0_sp/b
        h=d

        do i=1_i4,MAXIT
            a=-i*(nm1+i)
            b=b+2.0_sp
            d=1.0_sp/(a*d+b)
            c=b+a/c
            del=c*d
            h=h*del
            if (abs(del-1.0_sp) <= EPS) exit
        end do

        if (i > MAXIT) call nrerror('expint_sp: continued fraction failed')
        expint_sp=h*exp(-x)

    else
        if (nm1 /= 0_i4) then
        expint_sp=1.0_sp/nm1
        else
            expint_sp=-log(x)-EULER
        end if

        fact=1.0_sp

        do i=1_i4,MAXIT
            fact=-fact*x/i
            if (i /= nm1) then
                del=-fact/(i-nm1)
            else
                del=fact*(-log(x)-EULER+sum(1.0_sp/arth(1_i4,1_i4,nm1)))
            end if

            expint_sp=expint_sp+del

            if (abs(del) < abs(expint_sp)*EPS) exit
        end do

        if (i > MAXIT) call nrerror('expint_sp: series failed')
    end if
END FUNCTION expint_sp

FUNCTION expint_dp(n,x)
    USE mo_kind,            ONLY : i4, dp
    USE mo_nrutil,          ONLY : arth,assert,nrerror
    USE mo_constants,       ONLY : EULER_D

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: n
    REAL(dp), INTENT(IN) :: x
    REAL(dp) :: expint_dp
    INTEGER(i4), PARAMETER :: MAXIT=100_i4
    REAL(dp), PARAMETER :: EPS=epsilon(x),BIG=huge(x)*EPS
    INTEGER(i4) :: i,nm1
    REAL(dp) :: a,b,c,d,del,fact,h

    call assert(n >= 0_i4, x >= 0.0_dp, (x > 0.0_dp .or. n > 1_i4), &
    'expint_dp args')

    if (n == 0_i4) then
        expint_dp=exp(-x)/x
        RETURN
    end if

    nm1=n-1_i4

    if (x == 0.0_dp) then
        expint_dp=1.0_dp/nm1

    else if (x > 1.0_dp) then
        b=x+n
        c=BIG
        d=1.0_dp/b
        h=d

        do i=1_i4,MAXIT
            a=-i*(nm1+i)
            b=b+2.0_dp
            d=1.0_dp/(a*d+b)
            c=b+a/c
            del=c*d
            h=h*del
            if (abs(del-1.0_dp) <= EPS) exit
        end do

        if (i > MAXIT) call nrerror('expint_dp: continued fraction failed')
        expint_dp=h*exp(-x)

    else
        if (nm1 /= 0_i4) then
            expint_dp=1.0_dp/nm1
        else
            expint_dp=-log(x)-EULER_D
        end if

        fact=1.0_dp

        do i=1_i4,MAXIT
            fact=-fact*x/i
            if (i /= nm1) then
                del=-fact/(i-nm1)
            else
                del=fact*(-log(x)-EULER_D+sum(1.0_dp/arth(1_i4,1_i4,nm1)))
            end if

            expint_dp=expint_dp+del

            if (abs(del) < abs(expint_dp)*EPS) exit
        end do

        if (i > MAXIT) call nrerror('expint_dp: series failed')
    end if
END FUNCTION expint_dp

FUNCTION ei_sp(x)

    USE mo_kind,            ONLY : i4, sp
    USE mo_nrutil,          ONLY : assert,nrerror
    USE mo_constants,       ONLY : EULER

    IMPLICIT NONE

    REAL(SP), INTENT(IN)    :: x
    REAL(SP)                :: ei_sp
    INTEGER(i4), PARAMETER  :: MAXIT=100_i4
    REAL(SP), PARAMETER     :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
    INTEGER(i4)             :: k
    REAL(SP)                :: fact,prev,sm,term

    call assert(x > 0.0_sp, 'ei_sp arg')

    if (x < FPMIN) then
        ei_sp=log(x)+EULER

    else if (x <= -log(EPS)) then
        sm=0.0_sp
        fact=1.0_sp

        do k=1_i4,MAXIT
            fact=fact*x/k
            term=fact/k
            sm=sm+term
            if (term < EPS*sm) exit
        end do

        if (k > MAXIT) call nrerror('series failed in ei_sp')
        ei_sp=sm+log(x)+EULER

    else
        sm=0.0_sp
        term=1.0_sp

        do k=1_i4,MAXIT
            prev=term
            term=term*k/x
            if (term < EPS) exit
            if (term < prev) then
                sm=sm+term
            else
                sm=sm-prev
                exit
            end if
        end do

        if (k > MAXIT) call nrerror('asymptotic failed in ei_sp')
        ei_sp=exp(x)*(1.0_sp+sm)/x

    end if

END FUNCTION ei_sp

FUNCTION ei_dp(x)

    USE mo_kind,            ONLY : i4, dp
    USE mo_nrutil,          ONLY : assert,nrerror
    USE mo_constants,       ONLY : EULER

    IMPLICIT NONE

    REAL(dp), INTENT(IN)    :: x
    REAL(dp)                :: ei_dp
    INTEGER(i4), PARAMETER  :: MAXIT=100_i4
    REAL(dp), PARAMETER     :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
    INTEGER(i4)             :: k
    REAL(dp)                :: fact,prev,sm,term

    call assert(x > 0.0_dp, 'ei_dp arg')

    if (x < FPMIN) then
        ei_dp=log(x)+EULER

    else if (x <= -log(EPS)) then
        sm=0.0_dp
        fact=1.0_dp

        do k=1_i4,MAXIT
            fact=fact*x/k
            term=fact/k
            sm=sm+term
            if (term < EPS*sm) exit
        end do

        if (k > MAXIT) call nrerror('series failed in ei_dp')
        ei_dp=sm+log(x)+EULER

    else
        sm=0.0_dp
        term=1.0_dp

        do k=1_i4,MAXIT
            prev=term
            term=term*k/x
            if (term < EPS) exit
            if (term < prev) then
                sm=sm+term
            else
                sm=sm-prev
                exit
            end if
        end do

        if (k > MAXIT) call nrerror('asymptotic failed in ei_dp')
        ei_dp=exp(x)*(1.0_dp+sm)/x

    end if

END FUNCTION ei_dp

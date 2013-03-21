SUBROUTINE mnbrak_sp(ax,bx,cx,fa,fb,fc,func)
  use mo_kind
  use mo_nrutil, ONLY : swap
  IMPLICIT NONE
  REAL(SP), INTENT(INOUT)   :: ax,bx
  REAL(SP), INTENT(OUT)     :: cx,fa,fb,fc
  INTERFACE
     FUNCTION func(x)
       use mo_kind
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: x
       REAL(SP)             :: func
     END FUNCTION func
  END INTERFACE
  REAL(SP), PARAMETER       :: GOLD=1.618034_sp,GLIMIT=100.0_sp,TINY=1.0e-20_sp
  REAL(SP)                  :: fu,q,r,u,ulim
  fa=func(ax)
  fb=func(bx)
  if (fb > fa) then
     call swap(ax,bx)
     call swap(fa,fb)
  end if
  cx=bx+GOLD*(bx-ax)
  fc=func(cx)
  do
     if (fb < fc) RETURN
     r=(bx-ax)*(fb-fc)
     q=(bx-cx)*(fb-fa)
     u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_sp*sign(max(abs(q-r),TINY),q-r))
     ulim=bx+GLIMIT*(cx-bx)
     if ((bx-u)*(u-cx) > 0.0) then
        fu=func(u)
        if (fu < fc) then
           ax=bx
           fa=fb
           bx=u
           fb=fu
           RETURN
        else if (fu > fb) then
           cx=u
           fc=fu
           RETURN
        end if
        u=cx+GOLD*(cx-bx)
        fu=func(u)
     else if ((cx-u)*(u-ulim) > 0.0) then
        fu=func(u)
        if (fu < fc) then
           bx=cx
           cx=u
           u=cx+GOLD*(cx-bx)
           call shft(fb,fc,fu,func(u))
        end if
     else if ((u-ulim)*(ulim-cx) >= 0.0) then
        u=ulim
        fu=func(u)
     else
        u=cx+GOLD*(cx-bx)
        fu=func(u)
     end if
     call shft(ax,bx,cx,u)
     call shft(fa,fb,fc,fu)
  end do
CONTAINS
  !BL
  SUBROUTINE shft(a,b,c,d)
    REAL(SP), INTENT(OUT)   :: a
    REAL(SP), INTENT(INOUT) :: b,c
    REAL(SP), INTENT(IN)    :: d
    a=b
    b=c
    c=d
  END SUBROUTINE shft
END SUBROUTINE mnbrak_sp

SUBROUTINE mnbrak_dp(ax,bx,cx,fa,fb,fc,func)
  use mo_kind
  use mo_nrutil, ONLY : swap
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT)   :: ax,bx
  REAL(DP), INTENT(OUT)     :: cx,fa,fb,fc
  INTERFACE
     FUNCTION func(x)
       use mo_kind
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP)             :: func
     END FUNCTION func
  END INTERFACE
  REAL(DP), PARAMETER       :: GOLD=1.618034_dp,GLIMIT=100.0_dp,TINY=1.0e-20_dp
  REAL(DP)                  :: fu,q,r,u,ulim
  fa=func(ax)
  fb=func(bx)
  if (fb > fa) then
     call swap(ax,bx)
     call swap(fa,fb)
  end if
  cx=bx+GOLD*(bx-ax)
  fc=func(cx)
  do
     if (fb < fc) RETURN
     r=(bx-ax)*(fb-fc)
     q=(bx-cx)*(fb-fa)
     u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_dp*sign(max(abs(q-r),TINY),q-r))
     ulim=bx+GLIMIT*(cx-bx)
     if ((bx-u)*(u-cx) > 0.0) then
        fu=func(u)
        if (fu < fc) then
           ax=bx
           fa=fb
           bx=u
           fb=fu
           RETURN
        else if (fu > fb) then
           cx=u
           fc=fu
           RETURN
        end if
        u=cx+GOLD*(cx-bx)
        fu=func(u)
     else if ((cx-u)*(u-ulim) > 0.0) then
        fu=func(u)
        if (fu < fc) then
           bx=cx
           cx=u
           u=cx+GOLD*(cx-bx)
           call shft(fb,fc,fu,func(u))
        end if
     else if ((u-ulim)*(ulim-cx) >= 0.0) then
        u=ulim
        fu=func(u)
     else
        u=cx+GOLD*(cx-bx)
        fu=func(u)
     end if
     call shft(ax,bx,cx,u)
     call shft(fa,fb,fc,fu)
  end do
CONTAINS
  !BL
  SUBROUTINE shft(a,b,c,d)
    REAL(DP), INTENT(OUT)   :: a
    REAL(DP), INTENT(INOUT) :: b,c
    REAL(DP), INTENT(IN)    :: d
    a=b
    b=c
    c=d
  END SUBROUTINE shft
END SUBROUTINE mnbrak_dp

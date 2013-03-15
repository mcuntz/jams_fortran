subroutine laguer_sp(a,x,its)
  use mo_kind,   only : i4, sp, spc
  use mo_nrutil, only : nrerror, poly, poly_term

  implicit none

  integer(i4),                intent(out)   :: its
  complex(spc),               intent(inout) :: x
  complex(spc), dimension(:), intent(in)    :: a
  real(sp),     parameter                   :: eps=epsilon(1.0_sp)
  integer(i4),  parameter                   :: mr=8,mt=10,maxit=mt*mr
  integer(i4)                               :: iter,m
  real(sp)                                  :: abx,abp,abm,err
  complex(spc)                              :: dx,x1,f,g,h,sq,gp,gm,g2
  complex(spc), dimension(size(a))          :: b,d
  real(sp),     dimension(mr)               :: frac = &
       (/ 0.5_sp,0.25_sp,0.75_sp,0.13_sp,0.38_sp,0.62_sp,0.88_sp,1.0_sp /)

  m=size(a)-1
  do iter=1,maxit
     its=iter
     abx=abs(x)
     b(m+1:1:-1)=poly_term(a(m+1:1:-1),x)
     d(m:1:-1)=poly_term(b(m+1:2:-1),x)
     f=poly(x,d(2:m))
     err=eps*poly(abx,abs(b(1:m+1)))
     if (abs(b(1)) <= err) return
     g=d(1)/b(1)
     g2=g*g
     h=g2-2.0_sp*f/b(1)
     sq=sqrt((m-1)*(m*h-g2))
     gp=g+sq
     gm=g-sq
     abp=abs(gp)
     abm=abs(gm)
     if (abp < abm) gp=gm
     if (max(abp,abm) > 0.0_sp) then
        dx=m/gp
     else
        dx=exp(cmplx(log(1.0_sp+abx),iter,kind=spc))
     end if
     x1=x-dx
     if (x == x1) return
     if (mod(iter,mt) /= 0) then
        x=x1
     else
        x=x-dx*frac(iter/mt) 
     end if
  end do
  call nrerror('laguer: too many iterations')
end subroutine laguer_sp

subroutine laguer_dp(a,x,its)
  use mo_kind,   only : i4, dp, dpc
  use mo_nrutil, only : nrerror, poly, poly_term

  implicit none

  integer(i4),                intent(out)   :: its
  complex(dpc),               intent(inout) :: x
  complex(dpc), dimension(:), intent(in)    :: a
  real(dp),     parameter                   :: eps=epsilon(1.0_dp)
  integer(i4),  parameter                   :: mr=8,mt=10,maxit=mt*mr
  integer(i4)                               :: iter,m
  real(dp)                                  :: abx,abp,abm,err
  complex(dpc)                              :: dx,x1,f,g,h,sq,gp,gm,g2
  complex(dpc), dimension(size(a))          :: b,d
  real(dp),     dimension(mr)               :: frac = &
       (/ 0.5_dp,0.25_dp,0.75_dp,0.13_dp,0.38_dp,0.62_dp,0.88_dp,1.0_dp /)

  m=size(a)-1
  do iter=1,maxit
     its=iter
     abx=abs(x)
     b(m+1:1:-1) = poly_term(a(m+1:1:-1),x)
     d(m:1:-1)   = poly_term(b(m+1:2:-1),x)
     f           = poly(x,d(2:m))
     err=eps*poly(abx,abs(b(1:m+1)))
     if (abs(b(1)) <= err) return
     g=d(1)/b(1)
     g2=g*g
     h=g2-2.0_dp*f/b(1)
     sq=sqrt((m-1)*(m*h-g2))
     gp=g+sq
     gm=g-sq
     abp=abs(gp)
     abm=abs(gm)
     if (abp < abm) gp=gm
     if (max(abp,abm) > 0.0_dp) then
        dx=m/gp
     else
        dx=exp(cmplx(log(1.0_dp+abx),iter,kind=dpc))
     end if
     x1=x-dx
     if (x == x1) return
     if (mod(iter,mt) /= 0) then
        x=x1
     else
        x=x-dx*frac(iter/mt)
     end if
  end do
  call nrerror('laguer: too many iterations')
end subroutine laguer_dp

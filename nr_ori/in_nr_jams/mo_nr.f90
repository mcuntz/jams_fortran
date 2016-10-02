MODULE mo_nr
	INTERFACE
		SUBROUTINE airy(x,ai,bi,aip,bip)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: ai,bi,aip,bip
		END SUBROUTINE airy
	END INTERFACE
	INTERFACE
		SUBROUTINE amebsa(p,y,pb,yb,ftol,func,iter,temptr)
		use mo_nrtype
		INTEGER(I4B), INTENT(INOUT) :: iter
		REAL(SP), INTENT(INOUT) :: yb
		REAL(SP), INTENT(IN) :: ftol,temptr
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y,pb
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE amebsa
	END INTERFACE
	INTERFACE
		SUBROUTINE amoeba(p,y,ftol,func,iter)
		use mo_nrtype
		INTEGER(I4B), INTENT(OUT) :: iter
		REAL(SP), INTENT(IN) :: ftol
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE amoeba
	END INTERFACE
	INTERFACE
		SUBROUTINE anneal(x,y,iorder)
		use mo_nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: iorder
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		END SUBROUTINE anneal
	END INTERFACE
	INTERFACE
		SUBROUTINE asolve(b,x,itrnsp)
		use mo_nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: b
		REAL(DP), DIMENSION(:), INTENT(OUT) :: x
		INTEGER(I4B), INTENT(IN) :: itrnsp
		END SUBROUTINE asolve
	END INTERFACE
	INTERFACE
		SUBROUTINE atimes(x,r,itrnsp)
		use mo_nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(OUT) :: r
		INTEGER(I4B), INTENT(IN) :: itrnsp
		END SUBROUTINE atimes
	END INTERFACE
	INTERFACE
		SUBROUTINE avevar(data,ave,var)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data
		REAL(SP), INTENT(OUT) :: ave,var
		END SUBROUTINE avevar
	END INTERFACE
	INTERFACE
		SUBROUTINE balanc(a)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		END SUBROUTINE balanc
	END INTERFACE
	INTERFACE
		SUBROUTINE banbks(a,m1,m2,al,indx,b)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: m1,m2
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,al
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
		END SUBROUTINE banbks
	END INTERFACE
	INTERFACE
		SUBROUTINE bandec(a,m1,m2,al,indx,d)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: m1,m2
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
		REAL(SP), INTENT(OUT) :: d
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: al
		END SUBROUTINE bandec
	END INTERFACE
	INTERFACE
		SUBROUTINE banmul(a,m1,m2,x,b)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: m1,m2
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(OUT) :: b
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
		END SUBROUTINE banmul
	END INTERFACE
	INTERFACE
		SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: d1,d2
		REAL(SP), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
		REAL(SP), DIMENSION(4,4), INTENT(OUT) :: c
		END SUBROUTINE bcucof
	END INTERFACE
	INTERFACE
		SUBROUTINE bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,&
			ansy1,ansy2)
		use mo_nrtype
		REAL(SP), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
		REAL(SP), INTENT(IN) :: x1l,x1u,x2l,x2u,x1,x2
		REAL(SP), INTENT(OUT) :: ansy,ansy1,ansy2
		END SUBROUTINE bcuint
	END INTERFACE
	INTERFACE beschb
		SUBROUTINE beschb_s(x,gam1,gam2,gampl,gammi)
		use mo_nrtype
		REAL(DP), INTENT(IN) :: x
		REAL(DP), INTENT(OUT) :: gam1,gam2,gampl,gammi
		END SUBROUTINE beschb_s
!BL
		SUBROUTINE beschb_v(x,gam1,gam2,gampl,gammi)
		use mo_nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(OUT) :: gam1,gam2,gampl,gammi
		END SUBROUTINE beschb_v
	END INTERFACE
	INTERFACE bessi
		FUNCTION bessi_s(n,x)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessi_s
		END FUNCTION bessi_s
!BL
		FUNCTION bessi_v(n,x)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessi_v
		END FUNCTION bessi_v
	END INTERFACE
	INTERFACE bessi0
		FUNCTION bessi0_s(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessi0_s
		END FUNCTION bessi0_s
!BL
		FUNCTION bessi0_v(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessi0_v
		END FUNCTION bessi0_v
	END INTERFACE
	INTERFACE bessi1
		FUNCTION bessi1_s(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessi1_s
		END FUNCTION bessi1_s
!BL
		FUNCTION bessi1_v(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessi1_v
		END FUNCTION bessi1_v
	END INTERFACE
	INTERFACE
		SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,xnu
		REAL(SP), INTENT(OUT) :: ri,rk,rip,rkp
		END SUBROUTINE bessik
	END INTERFACE
	INTERFACE bessj
		FUNCTION bessj_s(n,x)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessj_s
		END FUNCTION bessj_s
!BL
		FUNCTION bessj_v(n,x)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessj_v
		END FUNCTION bessj_v
	END INTERFACE
	INTERFACE bessj0
		FUNCTION bessj0_s(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessj0_s
		END FUNCTION bessj0_s
!BL
		FUNCTION bessj0_v(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessj0_v
		END FUNCTION bessj0_v
	END INTERFACE
	INTERFACE bessj1
		FUNCTION bessj1_s(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessj1_s
		END FUNCTION bessj1_s
!BL
		FUNCTION bessj1_v(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessj1_v
		END FUNCTION bessj1_v
	END INTERFACE
	INTERFACE bessjy
		SUBROUTINE bessjy_s(x,xnu,rj,ry,rjp,ryp)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,xnu
		REAL(SP), INTENT(OUT) :: rj,ry,rjp,ryp
		END SUBROUTINE bessjy_s
!BL
		SUBROUTINE bessjy_v(x,xnu,rj,ry,rjp,ryp)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: xnu
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
		END SUBROUTINE bessjy_v
	END INTERFACE
	INTERFACE bessk
		FUNCTION bessk_s(n,x)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessk_s
		END FUNCTION bessk_s
!BL
		FUNCTION bessk_v(n,x)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessk_v
		END FUNCTION bessk_v
	END INTERFACE
	INTERFACE bessk0
		FUNCTION bessk0_s(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessk0_s
		END FUNCTION bessk0_s
!BL
		FUNCTION bessk0_v(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessk0_v
		END FUNCTION bessk0_v
	END INTERFACE
	INTERFACE bessk1
		FUNCTION bessk1_s(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessk1_s
		END FUNCTION bessk1_s
!BL
		FUNCTION bessk1_v(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessk1_v
		END FUNCTION bessk1_v
	END INTERFACE
	INTERFACE bessy
		FUNCTION bessy_s(n,x)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessy_s
		END FUNCTION bessy_s
!BL
		FUNCTION bessy_v(n,x)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessy_v
		END FUNCTION bessy_v
	END INTERFACE
	INTERFACE bessy0
		FUNCTION bessy0_s(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessy0_s
		END FUNCTION bessy0_s
!BL
		FUNCTION bessy0_v(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessy0_v
		END FUNCTION bessy0_v
	END INTERFACE
	INTERFACE bessy1
		FUNCTION bessy1_s(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: bessy1_s
		END FUNCTION bessy1_s
!BL
		FUNCTION bessy1_v(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: bessy1_v
		END FUNCTION bessy1_v
	END INTERFACE
	INTERFACE beta
		FUNCTION beta_s(z,w)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: z,w
		REAL(SP) :: beta_s
		END FUNCTION beta_s
!BL
		FUNCTION beta_v(z,w)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: z,w
		REAL(SP), DIMENSION(size(z)) :: beta_v
		END FUNCTION beta_v
	END INTERFACE
	INTERFACE betacf
		FUNCTION betacf_s(a,b,x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b,x
		REAL(SP) :: betacf_s
		END FUNCTION betacf_s
!BL
		FUNCTION betacf_v(a,b,x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,x
		REAL(SP), DIMENSION(size(x)) :: betacf_v
		END FUNCTION betacf_v
	END INTERFACE
	INTERFACE betai
		FUNCTION betai_s(a,b,x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b,x
		REAL(SP) :: betai_s
		END FUNCTION betai_s
!BL
		FUNCTION betai_v(a,b,x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,x
		REAL(SP), DIMENSION(size(a)) :: betai_v
		END FUNCTION betai_v
	END INTERFACE
	INTERFACE bico
		FUNCTION bico_s(n,k)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n,k
		REAL(SP) :: bico_s
		END FUNCTION bico_s
!BL
		FUNCTION bico_v(n,k)
		use mo_nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n,k
		REAL(SP), DIMENSION(size(n)) :: bico_v
		END FUNCTION bico_v
	END INTERFACE
	INTERFACE
		FUNCTION bnldev(pp,n)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: pp
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP) :: bnldev
		END FUNCTION bnldev
	END INTERFACE
	INTERFACE
		FUNCTION brent(ax,bx,cx,func,tol,xmin)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: ax,bx,cx,tol
		REAL(SP), INTENT(OUT) :: xmin
		REAL(SP) :: brent
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION brent
	END INTERFACE
	INTERFACE
		SUBROUTINE broydn(x,check)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
		LOGICAL(LGT), INTENT(OUT) :: check
		END SUBROUTINE broydn
	END INTERFACE
	INTERFACE
		SUBROUTINE bsstep(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
		REAL(SP), INTENT(INOUT) :: x
		REAL(SP), INTENT(IN) :: htry,eps
		REAL(SP), INTENT(OUT) :: hdid,hnext
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE bsstep
	END INTERFACE
	INTERFACE
		SUBROUTINE caldat(julian,mm,id,iyyy)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: julian
		INTEGER(I4B), INTENT(OUT) :: mm,id,iyyy
		END SUBROUTINE caldat
	END INTERFACE
	INTERFACE
		FUNCTION chder(a,b,c)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(:), INTENT(IN) :: c
		REAL(SP), DIMENSION(size(c)) :: chder
		END FUNCTION chder
	END INTERFACE
	INTERFACE chebev
		FUNCTION chebev_s(a,b,c,x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b,x
		REAL(SP), DIMENSION(:), INTENT(IN) :: c
		REAL(SP) :: chebev_s
		END FUNCTION chebev_s
!BL
		FUNCTION chebev_v(a,b,c,x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(:), INTENT(IN) :: c,x
		REAL(SP), DIMENSION(size(x)) :: chebev_v
		END FUNCTION chebev_v
	END INTERFACE
	INTERFACE
		FUNCTION chebft(a,b,n,func)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(n) :: chebft
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION chebft
	END INTERFACE
	INTERFACE
		FUNCTION chebpc(c)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: c
		REAL(SP), DIMENSION(size(c)) :: chebpc
		END FUNCTION chebpc
	END INTERFACE
	INTERFACE
		FUNCTION chint(a,b,c)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(:), INTENT(IN) :: c
		REAL(SP), DIMENSION(size(c)) :: chint
		END FUNCTION chint
	END INTERFACE
	INTERFACE
		SUBROUTINE choldc(a,p)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: p
		END SUBROUTINE choldc
	END INTERFACE
	INTERFACE
		SUBROUTINE cholsl(a,p,b,x)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
		REAL(SP), DIMENSION(:), INTENT(IN) :: p,b
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
		END SUBROUTINE cholsl
	END INTERFACE
	INTERFACE
		SUBROUTINE chsone(bins,ebins,knstrn,df,chsq,prob)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: knstrn
		REAL(SP), INTENT(OUT) :: df,chsq,prob
		REAL(SP), DIMENSION(:), INTENT(IN) :: bins,ebins
		END SUBROUTINE chsone
	END INTERFACE
	INTERFACE
		SUBROUTINE chstwo(bins1,bins2,knstrn,df,chsq,prob)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: knstrn
		REAL(SP), INTENT(OUT) :: df,chsq,prob
		REAL(SP), DIMENSION(:), INTENT(IN) :: bins1,bins2
		END SUBROUTINE chstwo
	END INTERFACE
	INTERFACE
		SUBROUTINE cisi(x,ci,si)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: ci,si
		END SUBROUTINE cisi
	END INTERFACE
	INTERFACE
		SUBROUTINE cntab1(nn,chisq,df,prob,cramrv,ccc)
		use mo_nrtype
		INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: nn
		REAL(SP), INTENT(OUT) :: chisq,df,prob,cramrv,ccc
		END SUBROUTINE cntab1
	END INTERFACE
	INTERFACE
		SUBROUTINE cntab2(nn,h,hx,hy,hygx,hxgy,uygx,uxgy,uxy)
		use mo_nrtype
		INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: nn
		REAL(SP), INTENT(OUT) :: h,hx,hy,hygx,hxgy,uygx,uxgy,uxy
		END SUBROUTINE cntab2
	END INTERFACE
	INTERFACE
		FUNCTION convlv(data,respns,isign)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data
		REAL(SP), DIMENSION(:), INTENT(IN) :: respns
		INTEGER(I4B), INTENT(IN) :: isign
		REAL(SP), DIMENSION(size(data)) :: convlv
		END FUNCTION convlv
	END INTERFACE
	INTERFACE
		FUNCTION correl(data1,data2)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		REAL(SP), DIMENSION(size(data1)) :: correl
		END FUNCTION correl
	END INTERFACE
	INTERFACE
		SUBROUTINE cosft1(y)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		END SUBROUTINE cosft1
	END INTERFACE
	INTERFACE
		SUBROUTINE cosft2(y,isign)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE cosft2
	END INTERFACE
	INTERFACE
		SUBROUTINE covsrt(covar,maska)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
		END SUBROUTINE covsrt
	END INTERFACE
	INTERFACE
		SUBROUTINE cyclic(a,b,c,alpha,beta,r,x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN):: a,b,c,r
		REAL(SP), INTENT(IN) :: alpha,beta
		REAL(SP), DIMENSION(:), INTENT(OUT):: x
		END SUBROUTINE cyclic
	END INTERFACE
	INTERFACE
		SUBROUTINE daub4(a,isign)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE daub4
	END INTERFACE
	INTERFACE dawson
		FUNCTION dawson_s(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: dawson_s
		END FUNCTION dawson_s
!BL
		FUNCTION dawson_v(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: dawson_v
		END FUNCTION dawson_v
	END INTERFACE
	INTERFACE
		FUNCTION dbrent(ax,bx,cx,func,dbrent_dfunc,tol,xmin)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: ax,bx,cx,tol
		REAL(SP), INTENT(OUT) :: xmin
		REAL(SP) :: dbrent
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
!BL
			FUNCTION dbrent_dfunc(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: dbrent_dfunc
			END FUNCTION dbrent_dfunc
		END INTERFACE
		END FUNCTION dbrent
	END INTERFACE
	INTERFACE
		SUBROUTINE ddpoly(c,x,pd)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: c
		REAL(SP), DIMENSION(:), INTENT(OUT) :: pd
		END SUBROUTINE ddpoly
	END INTERFACE
	INTERFACE
		FUNCTION decchk(string,ch)
		use mo_nrtype
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: string
		CHARACTER(1), INTENT(OUT) :: ch
		LOGICAL(LGT) :: decchk
		END FUNCTION decchk
	END INTERFACE
	INTERFACE
		SUBROUTINE dfpmin(p,gtol,iter,fret,func,dfunc)
		use mo_nrtype
		INTEGER(I4B), INTENT(OUT) :: iter
		REAL(SP), INTENT(IN) :: gtol
		REAL(SP), INTENT(OUT) :: fret
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
		INTERFACE
			FUNCTION func(p)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: p
			REAL(SP) :: func
			END FUNCTION func
!BL
			FUNCTION dfunc(p)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: p
			REAL(SP), DIMENSION(size(p)) :: dfunc
			END FUNCTION dfunc
		END INTERFACE
		END SUBROUTINE dfpmin
	END INTERFACE
	INTERFACE
		FUNCTION dfridr(func,x,h,err)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,h
		REAL(SP), INTENT(OUT) :: err
		REAL(SP) :: dfridr
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION dfridr
	END INTERFACE
	INTERFACE
		SUBROUTINE dftcor(w,delta,a,b,endpts,corre,corim,corfac)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: w,delta,a,b
		REAL(SP), INTENT(OUT) :: corre,corim,corfac
		REAL(SP), DIMENSION(:), INTENT(IN) :: endpts
		END SUBROUTINE dftcor
	END INTERFACE
	INTERFACE
		SUBROUTINE dftint(func,a,b,w,cosint,sinint)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b,w
		REAL(SP), INTENT(OUT) :: cosint,sinint
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE dftint
	END INTERFACE
	INTERFACE
		SUBROUTINE difeq(k,k1,k2,jsf,is1,isf,indexv,s,y)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: is1,isf,jsf,k,k1,k2
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indexv
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: s
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: y
		END SUBROUTINE difeq
	END INTERFACE
	INTERFACE
		FUNCTION eclass(lista,listb,n)
		use mo_nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: lista,listb
		INTEGER(I4B), INTENT(IN) :: n
		INTEGER(I4B), DIMENSION(n) :: eclass
		END FUNCTION eclass
	END INTERFACE
	INTERFACE
		FUNCTION eclazz(equiv,n)
		use mo_nrtype
		INTERFACE
			FUNCTION equiv(i,j)
			use mo_nrtype
			LOGICAL(LGT) :: equiv
			INTEGER(I4B), INTENT(IN) :: i,j
			END FUNCTION equiv
		END INTERFACE
		INTEGER(I4B), INTENT(IN) :: n
		INTEGER(I4B), DIMENSION(n) :: eclazz
		END FUNCTION eclazz
	END INTERFACE
	INTERFACE
		FUNCTION ei(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: ei
		END FUNCTION ei
	END INTERFACE
	INTERFACE
		SUBROUTINE eigsrt(d,v)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: v
		END SUBROUTINE eigsrt
	END INTERFACE
	INTERFACE elle
		FUNCTION elle_s(phi,ak)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: phi,ak
		REAL(SP) :: elle_s
		END FUNCTION elle_s
!BL
		FUNCTION elle_v(phi,ak)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: phi,ak
		REAL(SP), DIMENSION(size(phi)) :: elle_v
		END FUNCTION elle_v
	END INTERFACE
	INTERFACE ellf
		FUNCTION ellf_s(phi,ak)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: phi,ak
		REAL(SP) :: ellf_s
		END FUNCTION ellf_s
!BL
		FUNCTION ellf_v(phi,ak)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: phi,ak
		REAL(SP), DIMENSION(size(phi)) :: ellf_v
		END FUNCTION ellf_v
	END INTERFACE
	INTERFACE ellpi
		FUNCTION ellpi_s(phi,en,ak)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: phi,en,ak
		REAL(SP) :: ellpi_s
		END FUNCTION ellpi_s
!BL
		FUNCTION ellpi_v(phi,en,ak)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: phi,en,ak
		REAL(SP), DIMENSION(size(phi)) :: ellpi_v
		END FUNCTION ellpi_v
	END INTERFACE
	INTERFACE
		SUBROUTINE elmhes(a)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		END SUBROUTINE elmhes
	END INTERFACE
	INTERFACE erf
		FUNCTION erf_s(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: erf_s
		END FUNCTION erf_s
!BL
		FUNCTION erf_v(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: erf_v
		END FUNCTION erf_v
	END INTERFACE
	INTERFACE erfc
		FUNCTION erfc_s(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: erfc_s
		END FUNCTION erfc_s
!BL
		FUNCTION erfc_v(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: erfc_v
		END FUNCTION erfc_v
	END INTERFACE
	INTERFACE erfcc
		FUNCTION erfcc_s(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: erfcc_s
		END FUNCTION erfcc_s
!BL
		FUNCTION erfcc_v(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: erfcc_v
		END FUNCTION erfcc_v
	END INTERFACE
	INTERFACE
		SUBROUTINE eulsum(sum,term,jterm)
		use mo_nrtype
		REAL(SP), INTENT(INOUT) :: sum
		REAL(SP), INTENT(IN) :: term
		INTEGER(I4B), INTENT(IN) :: jterm
		END SUBROUTINE eulsum
	END INTERFACE
	INTERFACE
		FUNCTION evlmem(fdt,d,xms)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: fdt,xms
		REAL(SP), DIMENSION(:), INTENT(IN) :: d
		REAL(SP) :: evlmem
		END FUNCTION evlmem
	END INTERFACE
	INTERFACE expdev
		SUBROUTINE expdev_s(harvest)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE expdev_s
!BL
		SUBROUTINE expdev_v(harvest)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE expdev_v
	END INTERFACE
	INTERFACE
		FUNCTION expint(n,x)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: expint
		END FUNCTION expint
	END INTERFACE
	INTERFACE factln
		FUNCTION factln_s(n)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP) :: factln_s
		END FUNCTION factln_s
!BL
		FUNCTION factln_v(n)
		use mo_nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
		REAL(SP), DIMENSION(size(n)) :: factln_v
		END FUNCTION factln_v
	END INTERFACE
	INTERFACE factrl
		FUNCTION factrl_s(n)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP) :: factrl_s
		END FUNCTION factrl_s
!BL
		FUNCTION factrl_v(n)
		use mo_nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
		REAL(SP), DIMENSION(size(n)) :: factrl_v
		END FUNCTION factrl_v
	END INTERFACE
	INTERFACE
		SUBROUTINE fasper(x,y,ofac,hifac,px,py,jmax,prob)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), INTENT(IN) :: ofac,hifac
		INTEGER(I4B), INTENT(OUT) :: jmax
		REAL(SP), INTENT(OUT) :: prob
		REAL(SP), DIMENSION(:), POINTER :: px,py
		END SUBROUTINE fasper
	END INTERFACE
	INTERFACE
		SUBROUTINE fdjac(x,fvec,df)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: fvec
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: df
		END SUBROUTINE fdjac
	END INTERFACE
	INTERFACE
		SUBROUTINE fgauss(x,a,y,dyda)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: y
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
		END SUBROUTINE fgauss
	END INTERFACE
	INTERFACE
		SUBROUTINE fit(x,y,a,b,siga,sigb,chi2,q,sig)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
		REAL(SP), DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
		END SUBROUTINE fit
	END INTERFACE
	INTERFACE
		SUBROUTINE fitexy(x,y,sigx,sigy,a,b,siga,sigb,chi2,q)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sigx,sigy
		REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
		END SUBROUTINE fitexy
	END INTERFACE
	INTERFACE
		SUBROUTINE fixrts(d)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
		END SUBROUTINE fixrts
	END INTERFACE
	INTERFACE
		FUNCTION fleg(x,n)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(n) :: fleg
		END FUNCTION fleg
	END INTERFACE
	INTERFACE
		SUBROUTINE flmoon(n,nph,jd,frac)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n,nph
		INTEGER(I4B), INTENT(OUT) :: jd
		REAL(SP), INTENT(OUT) :: frac
		END SUBROUTINE flmoon
	END INTERFACE
	INTERFACE four1
		SUBROUTINE four1_dp(data,isign)
		use mo_nrtype
		COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE four1_dp
!BL
		SUBROUTINE four1_sp(data,isign)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE four1_sp
	END INTERFACE
	INTERFACE
		SUBROUTINE four1_alt(data,isign)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE four1_alt
	END INTERFACE
	INTERFACE
		SUBROUTINE four1_gather(data,isign)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE four1_gather
	END INTERFACE
	INTERFACE
		SUBROUTINE four2(data,isign)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
		INTEGER(I4B),INTENT(IN) :: isign
		END SUBROUTINE four2
	END INTERFACE
	INTERFACE
		SUBROUTINE four2_alt(data,isign)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE four2_alt
	END INTERFACE
	INTERFACE
		SUBROUTINE four3(data,isign)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
		INTEGER(I4B),INTENT(IN) :: isign
		END SUBROUTINE four3
	END INTERFACE
	INTERFACE
		SUBROUTINE four3_alt(data,isign)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE four3_alt
	END INTERFACE
	INTERFACE
		SUBROUTINE fourcol(data,isign)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE fourcol
	END INTERFACE
	INTERFACE
		SUBROUTINE fourcol_3d(data,isign)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE fourcol_3d
	END INTERFACE
	INTERFACE
		SUBROUTINE fourn_gather(data,nn,isign)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nn
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE fourn_gather
	END INTERFACE
	INTERFACE fourrow
		SUBROUTINE fourrow_dp(data,isign)
		use mo_nrtype
		COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE fourrow_dp
!BL
		SUBROUTINE fourrow_sp(data,isign)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE fourrow_sp
	END INTERFACE
	INTERFACE
		SUBROUTINE fourrow_3d(data,isign)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE fourrow_3d
	END INTERFACE
	INTERFACE
		FUNCTION fpoly(x,n)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(n) :: fpoly
		END FUNCTION fpoly
	END INTERFACE
	INTERFACE
		SUBROUTINE fred2(a,b,t,f,w,g,ak)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(:), INTENT(OUT) :: t,f,w
		INTERFACE
			FUNCTION g(t)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: t
			REAL(SP), DIMENSION(size(t)) :: g
			END FUNCTION g
!BL
			FUNCTION ak(t,s)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
			REAL(SP), DIMENSION(size(t),size(s)) :: ak
			END FUNCTION ak
		END INTERFACE
		END SUBROUTINE fred2
	END INTERFACE
	INTERFACE
		FUNCTION fredin(x,a,b,t,f,w,g,ak)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,t,f,w
		REAL(SP), DIMENSION(size(x)) :: fredin
		INTERFACE
			FUNCTION g(t)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: t
			REAL(SP), DIMENSION(size(t)) :: g
			END FUNCTION g
!BL
			FUNCTION ak(t,s)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
			REAL(SP), DIMENSION(size(t),size(s)) :: ak
			END FUNCTION ak
		END INTERFACE
		END FUNCTION fredin
	END INTERFACE
	INTERFACE
		SUBROUTINE frenel(x,s,c)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: s,c
		END SUBROUTINE frenel
	END INTERFACE
	INTERFACE
		SUBROUTINE frprmn(p,ftol,iter,fret)
		use mo_nrtype
		INTEGER(I4B), INTENT(OUT) :: iter
		REAL(SP), INTENT(IN) :: ftol
		REAL(SP), INTENT(OUT) :: fret
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
		END SUBROUTINE frprmn
	END INTERFACE
	INTERFACE
		SUBROUTINE ftest(data1,data2,f,prob)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: f,prob
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		END SUBROUTINE ftest
	END INTERFACE
	INTERFACE
		FUNCTION gamdev(ia)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: ia
		REAL(SP) :: gamdev
		END FUNCTION gamdev
	END INTERFACE
	INTERFACE gammln
		FUNCTION gammln_s(xx)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: xx
		REAL(SP) :: gammln_s
		END FUNCTION gammln_s
!BL
		FUNCTION gammln_v(xx)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xx
		REAL(SP), DIMENSION(size(xx)) :: gammln_v
		END FUNCTION gammln_v
	END INTERFACE
	INTERFACE gammp
		FUNCTION gammp_s(a,x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,x
		REAL(SP) :: gammp_s
		END FUNCTION gammp_s
!BL
		FUNCTION gammp_v(a,x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
		REAL(SP), DIMENSION(size(a)) :: gammp_v
		END FUNCTION gammp_v
	END INTERFACE
	INTERFACE gammq
		FUNCTION gammq_s(a,x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,x
		REAL(SP) :: gammq_s
		END FUNCTION gammq_s
!BL
		FUNCTION gammq_v(a,x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
		REAL(SP), DIMENSION(size(a)) :: gammq_v
		END FUNCTION gammq_v
	END INTERFACE
	INTERFACE gasdev
		SUBROUTINE gasdev_s(harvest)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE gasdev_s
!BL
		SUBROUTINE gasdev_v(harvest)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE gasdev_v
	END INTERFACE
	INTERFACE
		SUBROUTINE gaucof(a,b,amu0,x,w)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: amu0
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
		END SUBROUTINE gaucof
	END INTERFACE
	INTERFACE
		SUBROUTINE gauher(x,w)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
		END SUBROUTINE gauher
	END INTERFACE
	INTERFACE
		SUBROUTINE gaujac(x,w,alf,bet)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: alf,bet
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
		END SUBROUTINE gaujac
	END INTERFACE
	INTERFACE
		SUBROUTINE gaulag(x,w,alf)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: alf
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
		END SUBROUTINE gaulag
	END INTERFACE
	INTERFACE
		SUBROUTINE gauleg(x1,x2,x,w)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x1,x2
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
		END SUBROUTINE gauleg
	END INTERFACE
	INTERFACE
		SUBROUTINE gaussj(a,b)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
		END SUBROUTINE gaussj
	END INTERFACE
	INTERFACE gcf
		FUNCTION gcf_s(a,x,gln)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,x
		REAL(SP), OPTIONAL, INTENT(OUT) :: gln
		REAL(SP) :: gcf_s
		END FUNCTION gcf_s
!BL
		FUNCTION gcf_v(a,x,gln)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
		REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
		REAL(SP), DIMENSION(size(a)) :: gcf_v
		END FUNCTION gcf_v
	END INTERFACE
	INTERFACE
		FUNCTION golden(ax,bx,cx,func,tol,xmin)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: ax,bx,cx,tol
		REAL(SP), INTENT(OUT) :: xmin
		REAL(SP) :: golden
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION golden
	END INTERFACE
	INTERFACE gser
		FUNCTION gser_s(a,x,gln)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,x
		REAL(SP), OPTIONAL, INTENT(OUT) :: gln
		REAL(SP) :: gser_s
		END FUNCTION gser_s
!BL
		FUNCTION gser_v(a,x,gln)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
		REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
		REAL(SP), DIMENSION(size(a)) :: gser_v
		END FUNCTION gser_v
	END INTERFACE
	INTERFACE
		SUBROUTINE hqr(a,wr,wi)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: wr,wi
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		END SUBROUTINE hqr
	END INTERFACE
	INTERFACE
		SUBROUTINE hunt(xx,x,jlo)
		use mo_nrtype
		INTEGER(I4B), INTENT(INOUT) :: jlo
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: xx
		END SUBROUTINE hunt
	END INTERFACE
	INTERFACE
		SUBROUTINE hypdrv(s,ry,rdyds)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: s
		REAL(SP), DIMENSION(:), INTENT(IN) :: ry
		REAL(SP), DIMENSION(:), INTENT(OUT) :: rdyds
		END SUBROUTINE hypdrv
	END INTERFACE
	INTERFACE
		FUNCTION hypgeo(a,b,c,z)
		use mo_nrtype
		COMPLEX(SPC), INTENT(IN) :: a,b,c,z
		COMPLEX(SPC) :: hypgeo
		END FUNCTION hypgeo
	END INTERFACE
	INTERFACE
		SUBROUTINE hypser(a,b,c,z,series,deriv)
		use mo_nrtype
		COMPLEX(SPC), INTENT(IN) :: a,b,c,z
		COMPLEX(SPC), INTENT(OUT) :: series,deriv
		END SUBROUTINE hypser
	END INTERFACE
	INTERFACE
		FUNCTION icrc(crc,buf,jinit,jrev)
		use mo_nrtype
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: buf
		INTEGER(I2B), INTENT(IN) :: crc,jinit
		INTEGER(I4B), INTENT(IN) :: jrev
		INTEGER(I2B) :: icrc
		END FUNCTION icrc
	END INTERFACE
	INTERFACE
		FUNCTION igray(n,is)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n,is
		INTEGER(I4B) :: igray
		END FUNCTION igray
	END INTERFACE
	INTERFACE
		RECURSIVE SUBROUTINE index_bypack(arr,index,partial)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: index
		INTEGER, OPTIONAL, INTENT(IN) :: partial
		END SUBROUTINE index_bypack
	END INTERFACE
	INTERFACE indexx
		SUBROUTINE indexx_sp(arr,index)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
		END SUBROUTINE indexx_sp
		SUBROUTINE indexx_i4b(iarr,index)
		use mo_nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
		END SUBROUTINE indexx_i4b
	END INTERFACE
	INTERFACE
		FUNCTION interp(uc)
		use mo_nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: uc
		REAL(DP), DIMENSION(2*size(uc,1)-1,2*size(uc,1)-1) :: interp
		END FUNCTION interp
	END INTERFACE
	INTERFACE
		FUNCTION rank(indx)
		use mo_nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
		INTEGER(I4B), DIMENSION(size(indx)) :: rank
		END FUNCTION rank
	END INTERFACE
	INTERFACE
		FUNCTION irbit1(iseed)
		use mo_nrtype
		INTEGER(I4B), INTENT(INOUT) :: iseed
		INTEGER(I4B) :: irbit1
		END FUNCTION irbit1
	END INTERFACE
	INTERFACE
		FUNCTION irbit2(iseed)
		use mo_nrtype
		INTEGER(I4B), INTENT(INOUT) :: iseed
		INTEGER(I4B) :: irbit2
		END FUNCTION irbit2
	END INTERFACE
	INTERFACE
		SUBROUTINE jacobi(a,d,v,nrot)
		use mo_nrtype
		INTEGER(I4B), INTENT(OUT) :: nrot
		REAL(SP), DIMENSION(:), INTENT(OUT) :: d
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
		END SUBROUTINE jacobi
	END INTERFACE
	INTERFACE
		SUBROUTINE jacobn(x,y,dfdx,dfdy)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dfdx
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dfdy
		END SUBROUTINE jacobn
	END INTERFACE
	INTERFACE
		FUNCTION julday(mm,id,iyyy)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: mm,id,iyyy
		INTEGER(I4B) :: julday
		END FUNCTION julday
	END INTERFACE
	INTERFACE
		SUBROUTINE kendl1(data1,data2,tau,z,prob)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: tau,z,prob
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		END SUBROUTINE kendl1
	END INTERFACE
	INTERFACE
		SUBROUTINE kendl2(tab,tau,z,prob)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: tab
		REAL(SP), INTENT(OUT) :: tau,z,prob
		END SUBROUTINE kendl2
	END INTERFACE
	INTERFACE
		FUNCTION kermom(y,m)
		use mo_nrtype
		REAL(DP), INTENT(IN) :: y
		INTEGER(I4B), INTENT(IN) :: m
		REAL(DP), DIMENSION(m) :: kermom
		END FUNCTION kermom
	END INTERFACE
	INTERFACE
		SUBROUTINE ks2d1s(x1,y1,quadvl,d1,prob)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x1,y1
		REAL(SP), INTENT(OUT) :: d1,prob
		INTERFACE
			SUBROUTINE quadvl(x,y,fa,fb,fc,fd)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x,y
			REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
			END SUBROUTINE quadvl
		END INTERFACE
		END SUBROUTINE ks2d1s
	END INTERFACE
	INTERFACE
		SUBROUTINE ks2d2s(x1,y1,x2,y2,d,prob)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x1,y1,x2,y2
		REAL(SP), INTENT(OUT) :: d,prob
		END SUBROUTINE ks2d2s
	END INTERFACE
	INTERFACE
		SUBROUTINE ksone(data,func,d,prob)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: d,prob
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE ksone
	END INTERFACE
	INTERFACE
		SUBROUTINE kstwo(data1,data2,d,prob)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: d,prob
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		END SUBROUTINE kstwo
	END INTERFACE
	INTERFACE
		SUBROUTINE laguer(a,x,its)
		use mo_nrtype
		INTEGER(I4B), INTENT(OUT) :: its
		COMPLEX(SPC), INTENT(INOUT) :: x
		COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
		END SUBROUTINE laguer
	END INTERFACE
	INTERFACE
		SUBROUTINE lfit(x,y,sig,a,maska,covar,chisq,funcs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
		REAL(SP), INTENT(OUT) :: chisq
		INTERFACE
			SUBROUTINE funcs(x,arr)
			use mo_nrtype
			REAL(SP),INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(OUT) :: arr
			END SUBROUTINE funcs
		END INTERFACE
		END SUBROUTINE lfit
	END INTERFACE
	INTERFACE
		SUBROUTINE linbcg(b,x,itol,tol,itmax,iter,err)
		use mo_nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: b
		REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
		INTEGER(I4B), INTENT(IN) :: itol,itmax
		REAL(DP), INTENT(IN) :: tol
		INTEGER(I4B), INTENT(OUT) :: iter
		REAL(DP), INTENT(OUT) :: err
		END SUBROUTINE linbcg
	END INTERFACE
	INTERFACE
		SUBROUTINE linmin(p,xi,fret)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: fret
		REAL(SP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
		END SUBROUTINE linmin
	END INTERFACE
	INTERFACE
		SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xold,g
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
		REAL(SP), INTENT(IN) :: fold,stpmax
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x
		REAL(SP), INTENT(OUT) :: f
		LOGICAL(LGT), INTENT(OUT) :: check
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP) :: func
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE lnsrch
	END INTERFACE
	INTERFACE
		FUNCTION locate(xx,x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xx
		REAL(SP), INTENT(IN) :: x
		INTEGER(I4B) :: locate
		END FUNCTION locate
	END INTERFACE
	INTERFACE
		FUNCTION lop(u)
		use mo_nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: u
		REAL(DP), DIMENSION(size(u,1),size(u,1)) :: lop
		END FUNCTION lop
	END INTERFACE
	INTERFACE
		SUBROUTINE lubksb(a,indx,b)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
		END SUBROUTINE lubksb
	END INTERFACE
	INTERFACE
		SUBROUTINE ludcmp(a,indx,d)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
		REAL(SP), INTENT(OUT) :: d
		END SUBROUTINE ludcmp
	END INTERFACE
	INTERFACE
		SUBROUTINE machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,&
			maxexp,eps,epsneg,xmin,xmax)
		use mo_nrtype
		INTEGER(I4B), INTENT(OUT) :: ibeta,iexp,irnd,it,machep,maxexp,&
			minexp,negep,ngrd
		REAL(SP), INTENT(OUT) :: eps,epsneg,xmax,xmin
		END SUBROUTINE machar
	END INTERFACE
	INTERFACE
		SUBROUTINE medfit(x,y,a,b,abdev)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), INTENT(OUT) :: a,b,abdev
		END SUBROUTINE medfit
	END INTERFACE
	INTERFACE
		SUBROUTINE memcof(data,xms,d)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: xms
		REAL(SP), DIMENSION(:), INTENT(IN) :: data
		REAL(SP), DIMENSION(:), INTENT(OUT) :: d
		END SUBROUTINE memcof
	END INTERFACE
	INTERFACE
		SUBROUTINE mgfas(u,maxcyc)
		use mo_nrtype
		REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
		INTEGER(I4B), INTENT(IN) :: maxcyc
		END SUBROUTINE mgfas
	END INTERFACE
	INTERFACE
		SUBROUTINE mglin(u,ncycle)
		use mo_nrtype
		REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
		INTEGER(I4B), INTENT(IN) :: ncycle
		END SUBROUTINE mglin
	END INTERFACE
	INTERFACE
		SUBROUTINE midexp(funk,aa,bb,s,n)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: aa,bb
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION funk(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: funk
			END FUNCTION funk
		END INTERFACE
		END SUBROUTINE midexp
	END INTERFACE
	INTERFACE
		SUBROUTINE midinf(funk,aa,bb,s,n)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: aa,bb
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION funk(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: funk
			END FUNCTION funk
		END INTERFACE
		END SUBROUTINE midinf
	END INTERFACE
	INTERFACE
		SUBROUTINE midpnt(func,a,b,s,n)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE midpnt
	END INTERFACE
	INTERFACE
		SUBROUTINE midsql(funk,aa,bb,s,n)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: aa,bb
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION funk(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: funk
			END FUNCTION funk
		END INTERFACE
		END SUBROUTINE midsql
	END INTERFACE
	INTERFACE
		SUBROUTINE midsqu(funk,aa,bb,s,n)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: aa,bb
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION funk(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: funk
			END FUNCTION funk
		END INTERFACE
		END SUBROUTINE midsqu
	END INTERFACE
	INTERFACE
		RECURSIVE SUBROUTINE miser(func,regn,ndim,npts,dith,ave,var)
		use mo_nrtype
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP) :: func
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			END FUNCTION func
		END INTERFACE
		REAL(SP), DIMENSION(:), INTENT(IN) :: regn
		INTEGER(I4B), INTENT(IN) :: ndim,npts
		REAL(SP), INTENT(IN) :: dith
		REAL(SP), INTENT(OUT) :: ave,var
		END SUBROUTINE miser
	END INTERFACE
	INTERFACE
		SUBROUTINE mmid(y,dydx,xs,htot,nstep,yout,derivs)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: nstep
		REAL(SP), INTENT(IN) :: xs,htot
		REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE mmid
	END INTERFACE
	INTERFACE
		SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
		use mo_nrtype
		REAL(SP), INTENT(INOUT) :: ax,bx
		REAL(SP), INTENT(OUT) :: cx,fa,fb,fc
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE mnbrak
	END INTERFACE
	INTERFACE
		SUBROUTINE mnewt(ntrial,x,tolx,tolf,usrfun)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: ntrial
		REAL(SP), INTENT(IN) :: tolx,tolf
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
		INTERFACE
			SUBROUTINE usrfun(x,fvec,fjac)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(OUT) :: fvec
			REAL(SP), DIMENSION(:,:), INTENT(OUT) :: fjac
			END SUBROUTINE usrfun
		END INTERFACE
		END SUBROUTINE mnewt
	END INTERFACE
	INTERFACE
		SUBROUTINE moment(data,ave,adev,sdev,var,skew,curt)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: ave,adev,sdev,var,skew,curt
		REAL(SP), DIMENSION(:), INTENT(IN) :: data
		END SUBROUTINE moment
	END INTERFACE
	INTERFACE
		SUBROUTINE mp2dfr(a,s,n,m)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		INTEGER(I4B), INTENT(OUT) :: m
		CHARACTER(1), DIMENSION(:), INTENT(INOUT) :: a
		CHARACTER(1), DIMENSION(:), INTENT(OUT) :: s
		END SUBROUTINE mp2dfr
	END INTERFACE
	INTERFACE
		SUBROUTINE mpdiv(q,r,u,v,n,m)
		use mo_nrtype
		CHARACTER(1), DIMENSION(:), INTENT(OUT) :: q,r
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
		INTEGER(I4B), INTENT(IN) :: n,m
		END SUBROUTINE mpdiv
	END INTERFACE
	INTERFACE
		SUBROUTINE mpinv(u,v,n,m)
		use mo_nrtype
		CHARACTER(1), DIMENSION(:), INTENT(OUT) :: u
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
		INTEGER(I4B), INTENT(IN) :: n,m
		END SUBROUTINE mpinv
	END INTERFACE
	INTERFACE
		SUBROUTINE mpmul(w,u,v,n,m)
		use mo_nrtype
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
		CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
		INTEGER(I4B), INTENT(IN) :: n,m
		END SUBROUTINE mpmul
	END INTERFACE
	INTERFACE
		SUBROUTINE mppi(n)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		END SUBROUTINE mppi
	END INTERFACE
	INTERFACE
		SUBROUTINE mprove(a,alud,indx,b,x)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,alud
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
		REAL(SP), DIMENSION(:), INTENT(IN) :: b
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
		END SUBROUTINE mprove
	END INTERFACE
	INTERFACE
		SUBROUTINE mpsqrt(w,u,v,n,m)
		use mo_nrtype
		CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w,u
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
		INTEGER(I4B), INTENT(IN) :: n,m
		END SUBROUTINE mpsqrt
	END INTERFACE
	INTERFACE
		SUBROUTINE mrqcof(x,y,sig,a,maska,alpha,beta,chisq,funcs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,a,sig
		REAL(SP), DIMENSION(:), INTENT(OUT) :: beta
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: alpha
		REAL(SP), INTENT(OUT) :: chisq
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
		INTERFACE
			SUBROUTINE funcs(x,a,yfit,dyda)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
			REAL(SP), DIMENSION(:), INTENT(OUT) :: yfit
			REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
			END SUBROUTINE funcs
		END INTERFACE
		END SUBROUTINE mrqcof
	END INTERFACE
	INTERFACE
		SUBROUTINE mrqmin(x,y,sig,a,maska,covar,alpha,chisq,funcs,alamda)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: covar,alpha
		REAL(SP), INTENT(OUT) :: chisq
		REAL(SP), INTENT(INOUT) :: alamda
		LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
		INTERFACE
			SUBROUTINE funcs(x,a,yfit,dyda)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
			REAL(SP), DIMENSION(:), INTENT(OUT) :: yfit
			REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
			END SUBROUTINE funcs
		END INTERFACE
		END SUBROUTINE mrqmin
	END INTERFACE
	INTERFACE
		SUBROUTINE newt(x,check)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
		LOGICAL(LGT), INTENT(OUT) :: check
		END SUBROUTINE newt
	END INTERFACE
	INTERFACE
		SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: ystart
		REAL(SP), INTENT(IN) :: x1,x2,eps,h1,hmin
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
!BL
			SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
			REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
			REAL(SP), INTENT(INOUT) :: x
			REAL(SP), INTENT(IN) :: htry,eps
			REAL(SP), INTENT(OUT) :: hdid,hnext
				INTERFACE
				SUBROUTINE derivs(x,y,dydx)
					use mo_nrtype
					REAL(SP), INTENT(IN) :: x
					REAL(SP), DIMENSION(:), INTENT(IN) :: y
					REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
					END SUBROUTINE derivs
				END INTERFACE
			END SUBROUTINE rkqs
		END INTERFACE
		END SUBROUTINE odeint
	END INTERFACE
	INTERFACE
		SUBROUTINE orthog(anu,alpha,beta,a,b)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: anu,alpha,beta
		REAL(SP), DIMENSION(:), INTENT(OUT) :: a,b
		END SUBROUTINE orthog
	END INTERFACE
	INTERFACE
		SUBROUTINE pade(cof,resid)
		use mo_nrtype
		REAL(DP), DIMENSION(:), INTENT(INOUT) :: cof
		REAL(SP), INTENT(OUT) :: resid
		END SUBROUTINE pade
	END INTERFACE
	INTERFACE
		FUNCTION pccheb(d)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: d
		REAL(SP), DIMENSION(size(d)) :: pccheb
		END FUNCTION pccheb
	END INTERFACE
	INTERFACE
		SUBROUTINE pcshft(a,b,d)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
		END SUBROUTINE pcshft
	END INTERFACE
	INTERFACE
		SUBROUTINE pearsn(x,y,r,prob,z)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: r,prob,z
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		END SUBROUTINE pearsn
	END INTERFACE
	INTERFACE
		SUBROUTINE period(x,y,ofac,hifac,px,py,jmax,prob)
		use mo_nrtype
		INTEGER(I4B), INTENT(OUT) :: jmax
		REAL(SP), INTENT(IN) :: ofac,hifac
		REAL(SP), INTENT(OUT) :: prob
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), DIMENSION(:), POINTER :: px,py
		END SUBROUTINE period
	END INTERFACE
	INTERFACE plgndr
		FUNCTION plgndr_s(l,m,x)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: l,m
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: plgndr_s
		END FUNCTION plgndr_s
!BL
		FUNCTION plgndr_v(l,m,x)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: l,m
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: plgndr_v
		END FUNCTION plgndr_v
	END INTERFACE
	INTERFACE
		FUNCTION poidev(xm)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: xm
		REAL(SP) :: poidev
		END FUNCTION poidev
	END INTERFACE
	INTERFACE
		FUNCTION polcoe(x,y)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), DIMENSION(size(x)) :: polcoe
		END FUNCTION polcoe
	END INTERFACE
	INTERFACE
		FUNCTION polcof(xa,ya)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
		REAL(SP), DIMENSION(size(xa)) :: polcof
		END FUNCTION polcof
	END INTERFACE
	INTERFACE
		SUBROUTINE poldiv(u,v,q,r)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: u,v
		REAL(SP), DIMENSION(:), INTENT(OUT) :: q,r
		END SUBROUTINE poldiv
	END INTERFACE
	INTERFACE
		SUBROUTINE polin2(x1a,x2a,ya,x1,x2,y,dy)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
		REAL(SP), INTENT(IN) :: x1,x2
		REAL(SP), INTENT(OUT) :: y,dy
		END SUBROUTINE polin2
	END INTERFACE
	INTERFACE
		SUBROUTINE polint(xa,ya,x,y,dy)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: y,dy
		END SUBROUTINE polint
	END INTERFACE
	INTERFACE
		SUBROUTINE powell(p,xi,ftol,iter,fret)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: xi
		INTEGER(I4B), INTENT(OUT) :: iter
		REAL(SP), INTENT(IN) :: ftol
		REAL(SP), INTENT(OUT) :: fret
		END SUBROUTINE powell
	END INTERFACE
	INTERFACE
		FUNCTION predic(data,d,nfut)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data,d
		INTEGER(I4B), INTENT(IN) :: nfut
		REAL(SP), DIMENSION(nfut) :: predic
		END FUNCTION predic
	END INTERFACE
	INTERFACE
		FUNCTION probks(alam)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: alam
		REAL(SP) :: probks
		END FUNCTION probks
	END INTERFACE
	INTERFACE psdes
		SUBROUTINE psdes_s(lword,rword)
		use mo_nrtype
		INTEGER(I4B), INTENT(INOUT) :: lword,rword
		END SUBROUTINE psdes_s
!BL
		SUBROUTINE psdes_v(lword,rword)
		use mo_nrtype
		INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: lword,rword
		END SUBROUTINE psdes_v
	END INTERFACE
	INTERFACE
		SUBROUTINE pwt(a,isign)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE pwt
	END INTERFACE
	INTERFACE
		SUBROUTINE pwtset(n)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		END SUBROUTINE pwtset
	END INTERFACE
	INTERFACE pythag
		FUNCTION pythag_dp(a,b)
		use mo_nrtype
		REAL(DP), INTENT(IN) :: a,b
		REAL(DP) :: pythag_dp
		END FUNCTION pythag_dp
!BL
		FUNCTION pythag_sp(a,b)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP) :: pythag_sp
		END FUNCTION pythag_sp
	END INTERFACE
	INTERFACE
		SUBROUTINE pzextr(iest,xest,yest,yz,dy)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: iest
		REAL(SP), INTENT(IN) :: xest
		REAL(SP), DIMENSION(:), INTENT(IN) :: yest
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
		END SUBROUTINE pzextr
	END INTERFACE
	INTERFACE
		SUBROUTINE qrdcmp(a,c,d,sing)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: c,d
		LOGICAL(LGT), INTENT(OUT) :: sing
		END SUBROUTINE qrdcmp
	END INTERFACE
	INTERFACE
		FUNCTION qromb(func,a,b)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP) :: qromb
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION qromb
	END INTERFACE
	INTERFACE
		FUNCTION qromo(func,a,b,choose)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP) :: qromo
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		INTERFACE
			SUBROUTINE choose(funk,aa,bb,s,n)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: aa,bb
			REAL(SP), INTENT(INOUT) :: s
			INTEGER(I4B), INTENT(IN) :: n
			INTERFACE
				FUNCTION funk(x)
				use mo_nrtype
				REAL(SP), DIMENSION(:), INTENT(IN) :: x
				REAL(SP), DIMENSION(size(x)) :: funk
				END FUNCTION funk
			END INTERFACE
			END SUBROUTINE choose
		END INTERFACE
		END FUNCTION qromo
	END INTERFACE
	INTERFACE
		SUBROUTINE qroot(p,b,c,eps)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP), INTENT(INOUT) :: b,c
		REAL(SP), INTENT(IN) :: eps
		END SUBROUTINE qroot
	END INTERFACE
	INTERFACE
		SUBROUTINE qrsolv(a,c,d,b)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
		REAL(SP), DIMENSION(:), INTENT(IN) :: c,d
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
		END SUBROUTINE qrsolv
	END INTERFACE
	INTERFACE
		SUBROUTINE qrupdt(r,qt,u,v)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: r,qt
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: u
		REAL(SP), DIMENSION(:), INTENT(IN) :: v
		END SUBROUTINE qrupdt
	END INTERFACE
	INTERFACE
		FUNCTION qsimp(func,a,b)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP) :: qsimp
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION qsimp
	END INTERFACE
	INTERFACE
		FUNCTION qtrap(func,a,b)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP) :: qtrap
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION qtrap
	END INTERFACE
	INTERFACE
		SUBROUTINE quadct(x,y,xx,yy,fa,fb,fc,fd)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP), DIMENSION(:), INTENT(IN) :: xx,yy
		REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
		END SUBROUTINE quadct
	END INTERFACE
	INTERFACE
		SUBROUTINE quadmx(a)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: a
		END SUBROUTINE quadmx
	END INTERFACE
	INTERFACE
		SUBROUTINE quadvl(x,y,fa,fb,fc,fd)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
		END SUBROUTINE quadvl
	END INTERFACE
	INTERFACE
		FUNCTION ran(idum)
		INTEGER(selected_int_kind(9)), INTENT(INOUT) :: idum
		REAL :: ran
		END FUNCTION ran
	END INTERFACE
	INTERFACE ran0
		SUBROUTINE ran0_s(harvest)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE ran0_s
!BL
		SUBROUTINE ran0_v(harvest)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE ran0_v
	END INTERFACE
	INTERFACE ran1
		SUBROUTINE ran1_s(harvest)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE ran1_s
!BL
		SUBROUTINE ran1_v(harvest)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE ran1_v
	END INTERFACE
	INTERFACE ran2
		SUBROUTINE ran2_s(harvest)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE ran2_s
!BL
		SUBROUTINE ran2_v(harvest)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE ran2_v
	END INTERFACE
	INTERFACE ran3
		SUBROUTINE ran3_s(harvest)
		use mo_nrtype
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE ran3_s
!BL
		SUBROUTINE ran3_v(harvest)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE ran3_v
	END INTERFACE
	INTERFACE
		SUBROUTINE ratint(xa,ya,x,y,dy)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: y,dy
		END SUBROUTINE ratint
	END INTERFACE
	INTERFACE
		SUBROUTINE ratlsq(func,a,b,mm,kk,cof,dev)
		use mo_nrtype
		REAL(DP), INTENT(IN) :: a,b
		INTEGER(I4B), INTENT(IN) :: mm,kk
		REAL(DP), DIMENSION(:), INTENT(OUT) :: cof
		REAL(DP), INTENT(OUT) :: dev
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(DP), DIMENSION(:), INTENT(IN) :: x
			REAL(DP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE ratlsq
	END INTERFACE
	INTERFACE ratval
		FUNCTION ratval_s(x,cof,mm,kk)
		use mo_nrtype
		REAL(DP), INTENT(IN) :: x
		INTEGER(I4B), INTENT(IN) :: mm,kk
		REAL(DP), DIMENSION(mm+kk+1), INTENT(IN) :: cof
		REAL(DP) :: ratval_s
		END FUNCTION ratval_s
!BL
		FUNCTION ratval_v(x,cof,mm,kk)
		use mo_nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		INTEGER(I4B), INTENT(IN) :: mm,kk
		REAL(DP), DIMENSION(mm+kk+1), INTENT(IN) :: cof
		REAL(DP), DIMENSION(size(x)) :: ratval_v
		END FUNCTION ratval_v
	END INTERFACE
	INTERFACE rc
		FUNCTION rc_s(x,y)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP) :: rc_s
		END FUNCTION rc_s
!BL
		FUNCTION rc_v(x,y)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), DIMENSION(size(x)) :: rc_v
		END FUNCTION rc_v
	END INTERFACE
	INTERFACE rd
		FUNCTION rd_s(x,y,z)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,y,z
		REAL(SP) :: rd_s
		END FUNCTION rd_s
!BL
		FUNCTION rd_v(x,y,z)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z
		REAL(SP), DIMENSION(size(x)) :: rd_v
		END FUNCTION rd_v
	END INTERFACE
	INTERFACE realft
		SUBROUTINE realft_dp(data,isign,zdata)
		use mo_nrtype
		REAL(DP), DIMENSION(:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		COMPLEX(DPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
		END SUBROUTINE realft_dp
!BL
		SUBROUTINE realft_sp(data,isign,zdata)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
		INTEGER(I4B), INTENT(IN) :: isign
		COMPLEX(SPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
		END SUBROUTINE realft_sp
	END INTERFACE
	INTERFACE
		RECURSIVE FUNCTION recur1(a,b) RESULT(u)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
		REAL(SP), DIMENSION(size(a)) :: u
		END FUNCTION recur1
	END INTERFACE
	INTERFACE
		FUNCTION recur2(a,b,c)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c
		REAL(SP), DIMENSION(size(a)) :: recur2
		END FUNCTION recur2
	END INTERFACE
	INTERFACE
		SUBROUTINE relax(u,rhs)
		use mo_nrtype
		REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
		END SUBROUTINE relax
	END INTERFACE
	INTERFACE
		SUBROUTINE relax2(u,rhs)
		use mo_nrtype
		REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
		END SUBROUTINE relax2
	END INTERFACE
	INTERFACE
	FUNCTION resid(u,rhs)
		use mo_nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,rhs
		REAL(DP), DIMENSION(size(u,1),size(u,1)) :: resid
		END FUNCTION resid
	END INTERFACE
	INTERFACE rf
		FUNCTION rf_s(x,y,z)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,y,z
		REAL(SP) :: rf_s
		END FUNCTION rf_s
!BL
		FUNCTION rf_v(x,y,z)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z
		REAL(SP), DIMENSION(size(x)) :: rf_v
		END FUNCTION rf_v
	END INTERFACE
	INTERFACE rj
		FUNCTION rj_s(x,y,z,p)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,y,z,p
		REAL(SP) :: rj_s
		END FUNCTION rj_s
!BL
		FUNCTION rj_v(x,y,z,p)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z,p
		REAL(SP), DIMENSION(size(x)) :: rj_v
		END FUNCTION rj_v
	END INTERFACE
	INTERFACE
		SUBROUTINE rk4(y,dydx,x,h,yout,derivs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
		REAL(SP), INTENT(IN) :: x,h
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE rk4
	END INTERFACE
	INTERFACE
		SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
		REAL(SP), INTENT(IN) :: x,h
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yout,yerr
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE rkck
	END INTERFACE
	INTERFACE
		SUBROUTINE rkdumb(vstart,x1,x2,nstep,derivs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: vstart
		REAL(SP), INTENT(IN) :: x1,x2
		INTEGER(I4B), INTENT(IN) :: nstep
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE rkdumb
	END INTERFACE
	INTERFACE
		SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
		REAL(SP), INTENT(INOUT) :: x
		REAL(SP), INTENT(IN) :: htry,eps
		REAL(SP), INTENT(OUT) :: hdid,hnext
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE rkqs
	END INTERFACE
	INTERFACE
		SUBROUTINE rlft2(data,spec,speq,isign)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: data
		COMPLEX(SPC), DIMENSION(:,:), INTENT(OUT) :: spec
		COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: speq
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE rlft2
	END INTERFACE
	INTERFACE
		SUBROUTINE rlft3(data,spec,speq,isign)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:,:), INTENT(INOUT) :: data
		COMPLEX(SPC), DIMENSION(:,:,:), INTENT(OUT) :: spec
		COMPLEX(SPC), DIMENSION(:,:), INTENT(OUT) :: speq
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE rlft3
	END INTERFACE
	INTERFACE
		SUBROUTINE rotate(r,qt,i,a,b)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), TARGET, INTENT(INOUT) :: r,qt
		INTEGER(I4B), INTENT(IN) :: i
		REAL(SP), INTENT(IN) :: a,b
		END SUBROUTINE rotate
	END INTERFACE
	INTERFACE
		SUBROUTINE rsolv(a,d,b)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
		REAL(SP), DIMENSION(:), INTENT(IN) :: d
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
		END SUBROUTINE rsolv
	END INTERFACE
	INTERFACE
		FUNCTION rstrct(uf)
		use mo_nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: uf
		REAL(DP), DIMENSION((size(uf,1)+1)/2,(size(uf,1)+1)/2) :: rstrct
		END FUNCTION rstrct
	END INTERFACE
	INTERFACE
		FUNCTION rtbis(func,x1,x2,xacc)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x1,x2,xacc
		REAL(SP) :: rtbis
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION rtbis
	END INTERFACE
	INTERFACE
		FUNCTION rtflsp(func,x1,x2,xacc)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x1,x2,xacc
		REAL(SP) :: rtflsp
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION rtflsp
	END INTERFACE
	INTERFACE
		FUNCTION rtnewt(funcd,x1,x2,xacc)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x1,x2,xacc
		REAL(SP) :: rtnewt
		INTERFACE
			SUBROUTINE funcd(x,fval,fderiv)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), INTENT(OUT) :: fval,fderiv
			END SUBROUTINE funcd
		END INTERFACE
		END FUNCTION rtnewt
	END INTERFACE
	INTERFACE
		FUNCTION rtsafe(funcd,x1,x2,xacc)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x1,x2,xacc
		REAL(SP) :: rtsafe
		INTERFACE
			SUBROUTINE funcd(x,fval,fderiv)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), INTENT(OUT) :: fval,fderiv
			END SUBROUTINE funcd
		END INTERFACE
		END FUNCTION rtsafe
	END INTERFACE
	INTERFACE
		FUNCTION rtsec(func,x1,x2,xacc)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x1,x2,xacc
		REAL(SP) :: rtsec
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION rtsec
	END INTERFACE
	INTERFACE
		SUBROUTINE rzextr(iest,xest,yest,yz,dy)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: iest
		REAL(SP), INTENT(IN) :: xest
		REAL(SP), DIMENSION(:), INTENT(IN) :: yest
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
		END SUBROUTINE rzextr
	END INTERFACE
	INTERFACE
		FUNCTION savgol(nl,nrr,ld,m)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: nl,nrr,ld,m
		REAL(SP), DIMENSION(nl+nrr+1) :: savgol
		END FUNCTION savgol
	END INTERFACE
	INTERFACE
		SUBROUTINE scrsho(func)
		use mo_nrtype
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE scrsho
	END INTERFACE
	INTERFACE
		FUNCTION select(k,arr)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: k
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		REAL(SP) :: select
		END FUNCTION select
	END INTERFACE
	INTERFACE
		FUNCTION select_bypack(k,arr)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: k
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		REAL(SP) :: select_bypack
		END FUNCTION select_bypack
	END INTERFACE
	INTERFACE
		SUBROUTINE select_heap(arr,heap)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		REAL(SP), DIMENSION(:), INTENT(OUT) :: heap
		END SUBROUTINE select_heap
	END INTERFACE
	INTERFACE
		FUNCTION select_inplace(k,arr)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: k
		REAL(SP), DIMENSION(:), INTENT(IN) :: arr
		REAL(SP) :: select_inplace
		END FUNCTION select_inplace
	END INTERFACE
	INTERFACE
		SUBROUTINE simplx(a,m1,m2,m3,icase,izrov,iposv)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		INTEGER(I4B), INTENT(IN) :: m1,m2,m3
		INTEGER(I4B), INTENT(OUT) :: icase
		INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: izrov,iposv
		END SUBROUTINE simplx
	END INTERFACE
	INTERFACE
		SUBROUTINE simpr(y,dydx,dfdx,dfdy,xs,htot,nstep,yout,derivs)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: xs,htot
		REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx,dfdx
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: dfdy
		INTEGER(I4B), INTENT(IN) :: nstep
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE simpr
	END INTERFACE
	INTERFACE
		SUBROUTINE sinft(y)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		END SUBROUTINE sinft
	END INTERFACE
	INTERFACE
		SUBROUTINE slvsm2(u,rhs)
		use mo_nrtype
		REAL(DP), DIMENSION(3,3), INTENT(OUT) :: u
		REAL(DP), DIMENSION(3,3), INTENT(IN) :: rhs
		END SUBROUTINE slvsm2
	END INTERFACE
	INTERFACE
		SUBROUTINE slvsml(u,rhs)
		use mo_nrtype
		REAL(DP), DIMENSION(3,3), INTENT(OUT) :: u
		REAL(DP), DIMENSION(3,3), INTENT(IN) :: rhs
		END SUBROUTINE slvsml
	END INTERFACE
	INTERFACE
		SUBROUTINE sncndn(uu,emmc,sn,cn,dn)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: uu,emmc
		REAL(SP), INTENT(OUT) :: sn,cn,dn
		END SUBROUTINE sncndn
	END INTERFACE
	INTERFACE
		FUNCTION snrm(sx,itol)
		use mo_nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: sx
		INTEGER(I4B), INTENT(IN) :: itol
		REAL(DP) :: snrm
		END FUNCTION snrm
	END INTERFACE
	INTERFACE
		SUBROUTINE sobseq(x,init)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x
		INTEGER(I4B), OPTIONAL, INTENT(IN) :: init
		END SUBROUTINE sobseq
	END INTERFACE
	INTERFACE
		SUBROUTINE solvde(itmax,conv,slowc,scalv,indexv,nb,y)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: itmax,nb
		REAL(SP), INTENT(IN) :: conv,slowc
		REAL(SP), DIMENSION(:), INTENT(IN) :: scalv
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indexv
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: y
		END SUBROUTINE solvde
	END INTERFACE
	INTERFACE
		SUBROUTINE sor(a,b,c,d,e,f,u,rjac)
		use mo_nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: a,b,c,d,e,f
		REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
		REAL(DP), INTENT(IN) :: rjac
		END SUBROUTINE sor
	END INTERFACE
	INTERFACE
		SUBROUTINE sort(arr)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort
	END INTERFACE
	INTERFACE
		SUBROUTINE sort2(arr,slave)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave
		END SUBROUTINE sort2
	END INTERFACE
	INTERFACE
		SUBROUTINE sort3(arr,slave1,slave2)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave1,slave2
		END SUBROUTINE sort3
	END INTERFACE
	INTERFACE
		SUBROUTINE sort_bypack(arr)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort_bypack
	END INTERFACE
	INTERFACE
		SUBROUTINE sort_byreshape(arr)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort_byreshape
	END INTERFACE
	INTERFACE
		SUBROUTINE sort_heap(arr)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort_heap
	END INTERFACE
	INTERFACE
		SUBROUTINE sort_pick(arr)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort_pick
	END INTERFACE
	INTERFACE
		SUBROUTINE sort_radix(arr)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort_radix
	END INTERFACE
	INTERFACE
		SUBROUTINE sort_shell(arr)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
		END SUBROUTINE sort_shell
	END INTERFACE
	INTERFACE
		SUBROUTINE spctrm(p,k,ovrlap,unit,n_window)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(OUT) :: p
		INTEGER(I4B), INTENT(IN) :: k
		LOGICAL(LGT), INTENT(IN) :: ovrlap
		INTEGER(I4B), OPTIONAL, INTENT(IN) :: n_window,unit
		END SUBROUTINE spctrm
	END INTERFACE
	INTERFACE
		SUBROUTINE spear(data1,data2,d,zd,probd,rs,probrs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		REAL(SP), INTENT(OUT) :: d,zd,probd,rs,probrs
		END SUBROUTINE spear
	END INTERFACE
	INTERFACE sphbes
		SUBROUTINE sphbes_s(n,x,sj,sy,sjp,syp)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: sj,sy,sjp,syp
		END SUBROUTINE sphbes_s
!BL
		SUBROUTINE sphbes_v(n,x,sj,sy,sjp,syp)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(OUT) :: sj,sy,sjp,syp
		END SUBROUTINE sphbes_v
	END INTERFACE
	INTERFACE
		SUBROUTINE splie2(x1a,x2a,ya,y2a)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: y2a
		END SUBROUTINE splie2
	END INTERFACE
	INTERFACE
		FUNCTION splin2(x1a,x2a,ya,y2a,x1,x2)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya,y2a
		REAL(SP), INTENT(IN) :: x1,x2
		REAL(SP) :: splin2
		END FUNCTION splin2
	END INTERFACE
	INTERFACE
		SUBROUTINE spline(x,y,yp1,ypn,y2)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
		REAL(SP), INTENT(IN) :: yp1,ypn
		REAL(SP), DIMENSION(:), INTENT(OUT) :: y2
		END SUBROUTINE spline
	END INTERFACE
	INTERFACE
		FUNCTION splint(xa,ya,y2a,x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: splint
		END FUNCTION splint
	END INTERFACE
	INTERFACE sprsax
		SUBROUTINE sprsax_dp(sa,x,b)
		use mo_nrtype
		TYPE(sprs2_dp), INTENT(IN) :: sa
		REAL(DP), DIMENSION (:), INTENT(IN) :: x
		REAL(DP), DIMENSION (:), INTENT(OUT) :: b
		END SUBROUTINE sprsax_dp
!BL
		SUBROUTINE sprsax_sp(sa,x,b)
		use mo_nrtype
		TYPE(sprs2_sp), INTENT(IN) :: sa
		REAL(SP), DIMENSION (:), INTENT(IN) :: x
		REAL(SP), DIMENSION (:), INTENT(OUT) :: b
		END SUBROUTINE sprsax_sp
	END INTERFACE
	INTERFACE sprsdiag
		SUBROUTINE sprsdiag_dp(sa,b)
		use mo_nrtype
		TYPE(sprs2_dp), INTENT(IN) :: sa
		REAL(DP), DIMENSION(:), INTENT(OUT) :: b
		END SUBROUTINE sprsdiag_dp
!BL
		SUBROUTINE sprsdiag_sp(sa,b)
		use mo_nrtype
		TYPE(sprs2_sp), INTENT(IN) :: sa
		REAL(SP), DIMENSION(:), INTENT(OUT) :: b
		END SUBROUTINE sprsdiag_sp
	END INTERFACE
	INTERFACE sprsin
		SUBROUTINE sprsin_sp(a,thresh,sa)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
		REAL(SP), INTENT(IN) :: thresh
		TYPE(sprs2_sp), INTENT(OUT) :: sa
		END SUBROUTINE sprsin_sp
!BL
		SUBROUTINE sprsin_dp(a,thresh,sa)
		use mo_nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
		REAL(DP), INTENT(IN) :: thresh
		TYPE(sprs2_dp), INTENT(OUT) :: sa
		END SUBROUTINE sprsin_dp
	END INTERFACE
	INTERFACE
		SUBROUTINE sprstp(sa)
		use mo_nrtype
		TYPE(sprs2_sp), INTENT(INOUT) :: sa
		END SUBROUTINE sprstp
	END INTERFACE
	INTERFACE sprstx
		SUBROUTINE sprstx_dp(sa,x,b)
		use mo_nrtype
		TYPE(sprs2_dp), INTENT(IN) :: sa
		REAL(DP), DIMENSION (:), INTENT(IN) :: x
		REAL(DP), DIMENSION (:), INTENT(OUT) :: b
		END SUBROUTINE sprstx_dp
!BL
		SUBROUTINE sprstx_sp(sa,x,b)
		use mo_nrtype
		TYPE(sprs2_sp), INTENT(IN) :: sa
		REAL(SP), DIMENSION (:), INTENT(IN) :: x
		REAL(SP), DIMENSION (:), INTENT(OUT) :: b
		END SUBROUTINE sprstx_sp
	END INTERFACE
	INTERFACE
		SUBROUTINE stifbs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
		REAL(SP), INTENT(IN) :: htry,eps
		REAL(SP), INTENT(INOUT) :: x
		REAL(SP), INTENT(OUT) :: hdid,hnext
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE stifbs
	END INTERFACE
	INTERFACE
		SUBROUTINE stiff(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
		REAL(SP), INTENT(INOUT) :: x
		REAL(SP), INTENT(IN) :: htry,eps
		REAL(SP), INTENT(OUT) :: hdid,hnext
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE stiff
	END INTERFACE
	INTERFACE
		SUBROUTINE stoerm(y,d2y,xs,htot,nstep,yout,derivs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: y,d2y
		REAL(SP), INTENT(IN) :: xs,htot
		INTEGER(I4B), INTENT(IN) :: nstep
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE stoerm
	END INTERFACE
	INTERFACE svbksb
		SUBROUTINE svbksb_dp(u,w,v,b,x)
		use mo_nrtype
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,v
		REAL(DP), DIMENSION(:), INTENT(IN) :: w,b
		REAL(DP), DIMENSION(:), INTENT(OUT) :: x
		END SUBROUTINE svbksb_dp
!BL
		SUBROUTINE svbksb_sp(u,w,v,b,x)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: u,v
		REAL(SP), DIMENSION(:), INTENT(IN) :: w,b
		REAL(SP), DIMENSION(:), INTENT(OUT) :: x
		END SUBROUTINE svbksb_sp
	END INTERFACE
	INTERFACE svdcmp
		SUBROUTINE svdcmp_dp(a,w,v)
		use mo_nrtype
		REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(DP), DIMENSION(:), INTENT(OUT) :: w
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: v
		END SUBROUTINE svdcmp_dp
!BL
		SUBROUTINE svdcmp_sp(a,w,v)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: w
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
		END SUBROUTINE svdcmp_sp
	END INTERFACE
	INTERFACE
		SUBROUTINE svdfit(x,y,sig,a,v,w,chisq,funcs)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
		REAL(SP), DIMENSION(:), INTENT(OUT) :: a,w
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
		REAL(SP), INTENT(OUT) :: chisq
		INTERFACE
			FUNCTION funcs(x,n)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			INTEGER(I4B), INTENT(IN) :: n
			REAL(SP), DIMENSION(n) :: funcs
			END FUNCTION funcs
		END INTERFACE
		END SUBROUTINE svdfit
	END INTERFACE
	INTERFACE
		SUBROUTINE svdvar(v,w,cvm)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(IN) :: v
		REAL(SP), DIMENSION(:), INTENT(IN) :: w
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: cvm
		END SUBROUTINE svdvar
	END INTERFACE
	INTERFACE
		FUNCTION toeplz(r,y)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: r,y
		REAL(SP), DIMENSION(size(y)) :: toeplz
		END FUNCTION toeplz
	END INTERFACE
	INTERFACE
		SUBROUTINE tptest(data1,data2,t,prob)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		REAL(SP), INTENT(OUT) :: t,prob
		END SUBROUTINE tptest
	END INTERFACE
	INTERFACE
		SUBROUTINE tqli(d,e,z)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: d,e
		REAL(SP), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
		END SUBROUTINE tqli
	END INTERFACE
	INTERFACE
		SUBROUTINE trapzd(func,a,b,s,n)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: a,b
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE trapzd
	END INTERFACE
	INTERFACE
		SUBROUTINE tred2(a,d,e,novectors)
		use mo_nrtype
		REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: d,e
		LOGICAL(LGT), OPTIONAL, INTENT(IN) :: novectors
		END SUBROUTINE tred2
	END INTERFACE
!	On a purely serial machine, for greater efficiency, remove
!	the generic name tridag from the following interface,
!	and put it on the next one after that.
	INTERFACE tridag
		RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
		REAL(SP), DIMENSION(:), INTENT(OUT) :: u
		END SUBROUTINE tridag_par
	END INTERFACE
	INTERFACE
		SUBROUTINE tridag_ser(a,b,c,r,u)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
		REAL(SP), DIMENSION(:), INTENT(OUT) :: u
		END SUBROUTINE tridag_ser
	END INTERFACE
	INTERFACE
		SUBROUTINE ttest(data1,data2,t,prob)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		REAL(SP), INTENT(OUT) :: t,prob
		END SUBROUTINE ttest
	END INTERFACE
	INTERFACE
		SUBROUTINE tutest(data1,data2,t,prob)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		REAL(SP), INTENT(OUT) :: t,prob
		END SUBROUTINE tutest
	END INTERFACE
	INTERFACE
		SUBROUTINE twofft(data1,data2,fft1,fft2)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
		COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: fft1,fft2
		END SUBROUTINE twofft
	END INTERFACE
	INTERFACE
		FUNCTION vander(x,q)
		use mo_nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x,q
		REAL(DP), DIMENSION(size(x)) :: vander
		END FUNCTION vander
	END INTERFACE
	INTERFACE
		SUBROUTINE vegas(region,func,init,ncall,itmx,nprn,tgral,sd,chi2a)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: region
		INTEGER(I4B), INTENT(IN) :: init,ncall,itmx,nprn
		REAL(SP), INTENT(OUT) :: tgral,sd,chi2a
		INTERFACE
			FUNCTION func(pt,wgt)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(IN) :: pt
			REAL(SP), INTENT(IN) :: wgt
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE vegas
	END INTERFACE
	INTERFACE
		SUBROUTINE voltra(t0,h,t,f,g,ak)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: t0,h
		REAL(SP), DIMENSION(:), INTENT(OUT) :: t
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: f
		INTERFACE
			FUNCTION g(t)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: t
			REAL(SP), DIMENSION(:), POINTER :: g
			END FUNCTION g
!BL
			FUNCTION ak(t,s)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: t,s
			REAL(SP), DIMENSION(:,:), POINTER :: ak
			END FUNCTION ak
		END INTERFACE
		END SUBROUTINE voltra
	END INTERFACE
	INTERFACE
		SUBROUTINE wt1(a,isign,wtstep)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		INTEGER(I4B), INTENT(IN) :: isign
		INTERFACE
			SUBROUTINE wtstep(a,isign)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
			INTEGER(I4B), INTENT(IN) :: isign
			END SUBROUTINE wtstep
		END INTERFACE
		END SUBROUTINE wt1
	END INTERFACE
	INTERFACE
		SUBROUTINE wtn(a,nn,isign,wtstep)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nn
		INTEGER(I4B), INTENT(IN) :: isign
		INTERFACE
			SUBROUTINE wtstep(a,isign)
			use mo_nrtype
			REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
			INTEGER(I4B), INTENT(IN) :: isign
			END SUBROUTINE wtstep
		END INTERFACE
		END SUBROUTINE wtn
	END INTERFACE
	INTERFACE
		FUNCTION wwghts(n,h,kermom)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), INTENT(IN) :: h
		REAL(SP), DIMENSION(n) :: wwghts
		INTERFACE
			FUNCTION kermom(y,m)
			use mo_nrtype
			REAL(DP), INTENT(IN) :: y
			INTEGER(I4B), INTENT(IN) :: m
			REAL(DP), DIMENSION(m) :: kermom
			END FUNCTION kermom
		END INTERFACE
		END FUNCTION wwghts
	END INTERFACE
	INTERFACE
		SUBROUTINE zbrac(func,x1,x2,succes)
		use mo_nrtype
		REAL(SP), INTENT(INOUT) :: x1,x2
		LOGICAL(LGT), INTENT(OUT) :: succes
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE zbrac
	END INTERFACE
	INTERFACE
		SUBROUTINE zbrak(func,x1,x2,n,xb1,xb2,nb)
		use mo_nrtype
		INTEGER(I4B), INTENT(IN) :: n
		INTEGER(I4B), INTENT(OUT) :: nb
		REAL(SP), INTENT(IN) :: x1,x2
		REAL(SP), DIMENSION(:), POINTER :: xb1,xb2
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END SUBROUTINE zbrak
	END INTERFACE
	INTERFACE
		FUNCTION zbrent(func,x1,x2,tol)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x1,x2,tol
		REAL(SP) :: zbrent
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION zbrent
	END INTERFACE
	INTERFACE
		SUBROUTINE zrhqr(a,rtr,rti)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: rtr,rti
		END SUBROUTINE zrhqr
	END INTERFACE
	INTERFACE
		FUNCTION zriddr(func,x1,x2,xacc)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x1,x2,xacc
		REAL(SP) :: zriddr
		INTERFACE
			FUNCTION func(x)
			use mo_nrtype
			REAL(SP), INTENT(IN) :: x
			REAL(SP) :: func
			END FUNCTION func
		END INTERFACE
		END FUNCTION zriddr
	END INTERFACE
	INTERFACE
		SUBROUTINE zroots(a,roots,polish)
		use mo_nrtype
		COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
		COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: roots
		LOGICAL(LGT), INTENT(IN) :: polish
		END SUBROUTINE zroots
	END INTERFACE

CONTAINS

! -----------------------------------------------------------

	SUBROUTINE airy(x,ai,bi,aip,bip)
	use mo_nrtype
	use mo_nr, ONLY : bessik,bessjy
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), INTENT(OUT) :: ai,bi,aip,bip
	REAL(SP) :: absx,ri,rip,rj,rjp,rk,rkp,rootx,ry,ryp,z
	REAL(SP), PARAMETER :: THIRD=1.0_sp/3.0_sp,TWOTHR=2.0_sp/3.0_sp, &
		ONOVRT=0.5773502691896258_sp
	absx=abs(x)
	rootx=sqrt(absx)
	z=TWOTHR*absx*rootx
	if (x > 0.0) then
		call bessik(z,THIRD,ri,rk,rip,rkp)
		ai=rootx*ONOVRT*rk/PI
		bi=rootx*(rk/PI+2.0_sp*ONOVRT*ri)
		call bessik(z,TWOTHR,ri,rk,rip,rkp)
		aip=-x*ONOVRT*rk/PI
		bip=x*(rk/PI+2.0_sp*ONOVRT*ri)
	else if (x < 0.0) then
		call bessjy(z,THIRD,rj,ry,rjp,ryp)
		ai=0.5_sp*rootx*(rj-ONOVRT*ry)
		bi=-0.5_sp*rootx*(ry+ONOVRT*rj)
		call bessjy(z,TWOTHR,rj,ry,rjp,ryp)
		aip=0.5_sp*absx*(ONOVRT*ry+rj)
		bip=0.5_sp*absx*(ONOVRT*rj-ry)
	else
		ai=0.3550280538878172_sp
		bi=ai/ONOVRT
		aip=-0.2588194037928068_sp
		bip=-aip/ONOVRT
	end if
	END SUBROUTINE airy

! -----------------------------------------------------------

	SUBROUTINE amebsa(p,y,pb,yb,ftol,func,iter,temptr)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,imaxloc,iminloc,swap
	use mo_nr, ONLY : ran1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: iter
	REAL(SP), INTENT(INOUT) :: yb
	REAL(SP), INTENT(IN) :: ftol,temptr
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y,pb
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: NMAX=200
	INTEGER(I4B) :: ihi,ndim
	REAL(SP) :: yhi
	REAL(SP), DIMENSION(size(p,2)) :: psum
	call amebsa_private
	CONTAINS
!BL
	SUBROUTINE amebsa_private
	INTEGER(I4B) :: i,ilo,inhi
	REAL(SP) :: rtol,ylo,ynhi,ysave,ytry
	REAL(SP), DIMENSION(size(y)) :: yt,harvest
	ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,size(pb),'amebsa')
	psum(:)=sum(p(:,:),dim=1)
	do
		call ran1(harvest)
		yt(:)=y(:)-temptr*log(harvest)
		ilo=iminloc(yt(:))
		ylo=yt(ilo)
		ihi=imaxloc(yt(:))
		yhi=yt(ihi)
		yt(ihi)=ylo
		inhi=imaxloc(yt(:))
		ynhi=yt(inhi)
		rtol=2.0_sp*abs(yhi-ylo)/(abs(yhi)+abs(ylo))
		if (rtol < ftol .or. iter < 0) then
			call swap(y(1),y(ilo))
			call swap(p(1,:),p(ilo,:))
			RETURN
		end if
		ytry=amotsa(-1.0_sp)
		iter=iter-1
		if (ytry <= ylo) then
			ytry=amotsa(2.0_sp)
			iter=iter-1
		else if (ytry >= ynhi) then
			ysave=yhi
			ytry=amotsa(0.5_sp)
			iter=iter-1
			if (ytry >= ysave) then
				p(:,:)=0.5_sp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
				do i=1,ndim+1
					if (i /= ilo) y(i)=func(p(i,:))
				end do
				iter=iter-ndim
				psum(:)=sum(p(:,:),dim=1)
			end if
		end if
	end do
	END SUBROUTINE amebsa_private
!BL
	FUNCTION amotsa(fac)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: fac
	REAL(SP) :: amotsa
	REAL(SP) :: fac1,fac2,yflu,ytry,harv
	REAL(SP), DIMENSION(size(p,2)) :: ptry
	fac1=(1.0_sp-fac)/ndim
	fac2=fac1-fac
	ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
	ytry=func(ptry)
	if (ytry <= yb) then
		pb(:)=ptry(:)
		yb=ytry
	end if
	call ran1(harv)
	yflu=ytry+temptr*log(harv)
	if (yflu < yhi) then
		y(ihi)=ytry
		yhi=yflu
		psum(:)=psum(:)-p(ihi,:)+ptry(:)
		p(ihi,:)=ptry(:)
	end if
	amotsa=yflu
	END FUNCTION amotsa
	END SUBROUTINE amebsa

! -----------------------------------------------------------

	SUBROUTINE amoeba(p,y,ftol,func,iter)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(SP), INTENT(IN) :: ftol
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=5000
	INTEGER(I4B) :: ihi,ndim
	REAL(SP), DIMENSION(size(p,2)) :: psum
	call amoeba_private
	CONTAINS
!BL
	SUBROUTINE amoeba_private
	IMPLICIT NONE
	INTEGER(I4B) :: i,ilo,inhi
	REAL(SP) :: rtol,ysave,ytry,ytmp
	ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
	iter=0
	psum(:)=sum(p(:,:),dim=1)
	do
		ilo=iminloc(y(:))
		ihi=imaxloc(y(:))
		ytmp=y(ihi)
		y(ihi)=y(ilo)
		inhi=imaxloc(y(:))
		y(ihi)=ytmp
		rtol=2.0_sp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
		if (rtol < ftol) then
			call swap(y(1),y(ilo))
			call swap(p(1,:),p(ilo,:))
			RETURN
		end if
		if (iter >= ITMAX) call nrerror('ITMAX exceeded in amoeba')
		ytry=amotry(-1.0_sp)
		iter=iter+1
		if (ytry <= y(ilo)) then
			ytry=amotry(2.0_sp)
			iter=iter+1
		else if (ytry >= y(inhi)) then
			ysave=y(ihi)
			ytry=amotry(0.5_sp)
			iter=iter+1
			if (ytry >= ysave) then
				p(:,:)=0.5_sp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
				do i=1,ndim+1
					if (i /= ilo) y(i)=func(p(i,:))
				end do
				iter=iter+ndim
				psum(:)=sum(p(:,:),dim=1)
			end if
		end if
	end do
	END SUBROUTINE amoeba_private
!BL
	FUNCTION amotry(fac)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: fac
	REAL(SP) :: amotry
	REAL(SP) :: fac1,fac2,ytry
	REAL(SP), DIMENSION(size(p,2)) :: ptry
	fac1=(1.0_sp-fac)/ndim
	fac2=fac1-fac
	ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
	ytry=func(ptry)
	if (ytry < y(ihi)) then
		y(ihi)=ytry
		psum(:)=psum(:)-p(ihi,:)+ptry(:)
		p(ihi,:)=ptry(:)
	end if
	amotry=ytry
	END FUNCTION amotry
	END SUBROUTINE amoeba

! -----------------------------------------------------------

	SUBROUTINE anneal(x,y,iorder)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq,swap
	use mo_nr, ONLY : ran1
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: iorder
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	INTEGER(I4B), DIMENSION(6) :: n
	INTEGER(I4B) :: i1,i2,j,k,nlimit,ncity,nn,nover,nsucc
	REAL(SP) :: de,harvest,path,t,tfactr
	LOGICAL(LGT) :: ans
	ncity=assert_eq(size(x),size(y),size(iorder),'anneal')
	nover=100*ncity
	nlimit=10*ncity
	tfactr=0.9_sp
	t=0.5_sp
	path=sum(alen_v(x(iorder(1:ncity-1)),x(iorder(2:ncity)),&
		y(iorder(1:ncity-1)),y(iorder(2:ncity))))
	i1=iorder(ncity)
	i2=iorder(1)
	path=path+alen(x(i1),x(i2),y(i1),y(i2))
	do j=1,100
		nsucc=0
		do k=1,nover
			do
				call ran1(harvest)
				n(1)=1+int(ncity*harvest)
				call ran1(harvest)
				n(2)=1+int((ncity-1)*harvest)
				if (n(2) >= n(1)) n(2)=n(2)+1
				nn=1+mod((n(1)-n(2)+ncity-1),ncity)
				if (nn >= 3) exit
			end do
			call ran1(harvest)
			if (harvest < 0.5_sp) then
				call ran1(harvest)
				n(3)=n(2)+int(abs(nn-2)*harvest)+1
				n(3)=1+mod(n(3)-1,ncity)
				call trncst(x,y,iorder,n,de)
				call metrop(de,t,ans)
				if (ans) then
					nsucc=nsucc+1
					path=path+de
					call trnspt(iorder,n)
				end if
			else
				call revcst(x,y,iorder,n,de)
				call metrop(de,t,ans)
				if (ans) then
					nsucc=nsucc+1
					path=path+de
					call revers(iorder,n)
				end if
			end if
			if (nsucc >= nlimit) exit
		end do
		write(*,*)
		write(*,*) 'T =',t,' Path Length =',path
		write(*,*) 'Successful Moves: ',nsucc
		t=t*tfactr
		if (nsucc == 0) RETURN
	end do
	CONTAINS
!BL
	FUNCTION alen(x1,x2,y1,y2)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,y1,y2
	REAL(SP) :: alen
	alen=sqrt((x2-x1)**2+(y2-y1)**2)
	END FUNCTION alen
!BL
	FUNCTION alen_v(x1,x2,y1,y2)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x1,x2,y1,y2
	REAL(SP), DIMENSION(size(x1)) :: alen_v
	alen_v=sqrt((x2-x1)**2+(y2-y1)**2)
	END FUNCTION alen_v
!BL
	SUBROUTINE metrop(de,t,ans)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: de,t
	LOGICAL(LGT), INTENT(OUT) :: ans
	call ran1(harvest)
	ans=(de < 0.0) .or. (harvest < exp(-de/t))
	END SUBROUTINE metrop
!BL
	SUBROUTINE revcst(x,y,iorder,n,de)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iorder
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: n
	REAL(SP), INTENT(OUT) :: de
	INTEGER(I4B) :: ncity
	REAL(SP), DIMENSION(4) :: xx,yy
	ncity=size(x)
	n(3)=1+mod((n(1)+ncity-2),ncity)
	n(4)=1+mod(n(2),ncity)
	xx(1:4)=x(iorder(n(1:4)))
	yy(1:4)=y(iorder(n(1:4)))
	de=-alen(xx(1),xx(3),yy(1),yy(3))&
		-alen(xx(2),xx(4),yy(2),yy(4))&
		+alen(xx(1),xx(4),yy(1),yy(4))&
		+alen(xx(2),xx(3),yy(2),yy(3))
	END SUBROUTINE revcst
!BL
	SUBROUTINE revers(iorder,n)
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: iorder
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
	INTEGER(I4B) :: j,k,l,nn,ncity
	ncity=size(iorder)
	nn=(1+mod(n(2)-n(1)+ncity,ncity))/2
	do j=1,nn
		k=1+mod((n(1)+j-2),ncity)
		l=1+mod((n(2)-j+ncity),ncity)
		call swap(iorder(k),iorder(l))
	end do
	END SUBROUTINE revers
!BL
	SUBROUTINE trncst(x,y,iorder,n,de)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iorder
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: n
	REAL(SP), INTENT(OUT) :: de
	INTEGER(I4B) :: ncity
	REAL(SP), DIMENSION(6) :: xx,yy
	ncity=size(x)
	n(4)=1+mod(n(3),ncity)
	n(5)=1+mod((n(1)+ncity-2),ncity)
	n(6)=1+mod(n(2),ncity)
	xx(1:6)=x(iorder(n(1:6)))
	yy(1:6)=y(iorder(n(1:6)))
	de=-alen(xx(2),xx(6),yy(2),yy(6))&
		-alen(xx(1),xx(5),yy(1),yy(5))&
		-alen(xx(3),xx(4),yy(3),yy(4))&
		+alen(xx(1),xx(3),yy(1),yy(3))&
		+alen(xx(2),xx(4),yy(2),yy(4))&
		+alen(xx(5),xx(6),yy(5),yy(6))
	END SUBROUTINE trncst
!BL
	SUBROUTINE trnspt(iorder,n)
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: iorder
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
	INTEGER(I4B) :: m1,m2,m3,nn,ncity
	INTEGER(I4B), DIMENSION(size(iorder)) :: jorder
	ncity=size(iorder)
	m1=1+mod((n(2)-n(1)+ncity),ncity)
	m2=1+mod((n(5)-n(4)+ncity),ncity)
	m3=1+mod((n(3)-n(6)+ncity),ncity)
	jorder(1:m1)=iorder(1+mod((arth(1,1,m1)+n(1)-2),ncity))
	nn=m1
	jorder(nn+1:nn+m2)=iorder(1+mod((arth(1,1,m2)+n(4)-2),ncity))
	nn=nn+m2
	jorder(nn+1:nn+m3)=iorder(1+mod((arth(1,1,m3)+n(6)-2),ncity))
	iorder(1:ncity)=jorder(1:ncity)
	END SUBROUTINE trnspt
	END SUBROUTINE anneal

! -----------------------------------------------------------

MODULE arcode_info
	use mo_nrtype
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NWK=20
	TYPE arithcode
		INTEGER(I4B), DIMENSION(:), POINTER :: ilob,iupb,ncumfq
		INTEGER(I4B) :: jdif,nc,minint,nch,ncum,nrad
	END TYPE arithcode
CONTAINS
	SUBROUTINE arcode_allocate(acode,mc)
	use mo_nrtype
	IMPLICIT NONE
	TYPE(arithcode) :: acode
	INTEGER(I4B) :: mc
	allocate(acode%ilob(NWK),acode%iupb(NWK),acode%ncumfq(mc+2))
	END SUBROUTINE arcode_allocate
!BL
	SUBROUTINE arcode_deallocate(acode)
	use mo_nrtype
	IMPLICIT NONE
	TYPE(arithcode) :: acode
	deallocate(acode%ncumfq,acode%iupb,acode%ilob)
	nullify(acode%ilob)
	nullify(acode%iupb)
	nullify(acode%ncumfq)
	END SUBROUTINE arcode_deallocate
END MODULE arcode_info

	SUBROUTINE arcmak(nfreq,nradd,acode)
	use mo_nrtype
          use mo_nrutil, ONLY : cumsum,nrerror
	USE arcode_info
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: nradd
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nfreq
	TYPE(arithcode) :: acode
	INTEGER(I4B), PARAMETER :: MAXINT=huge(nradd)
	if (nradd > 256) call nrerror('output radix may not exceed 256 in arcmak')
	acode%minint=MAXINT/nradd
	acode%nch=size(nfreq)
	acode%nrad=nradd
	call arcode_allocate(acode,acode%nch)
	acode%ncumfq(1)=0
	acode%ncumfq(2:acode%nch+1)=cumsum(max(nfreq(1:acode%nch),1))
	acode%ncumfq(acode%nch+2)=acode%ncumfq(acode%nch+1)+1
	acode%ncum=acode%ncumfq(acode%nch+2)
	END SUBROUTINE arcmak

! -----------------------------------------------------------

	SUBROUTINE arcode(ich,codep,lcd,isign,acode)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror,reallocate
	USE arcode_info
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: ich,lcd
	INTEGER(I4B), INTENT(IN) :: isign
	CHARACTER(1), DIMENSION(:), POINTER :: codep
	TYPE(arithcode) :: acode
	INTEGER(I4B) :: ihi,j,ja,jh,jl,m
	if (isign == 0) then
		acode%jdif=acode%nrad-1
		acode%ilob(:)=0
		acode%iupb(:)=acode%nrad-1
		do j=NWK,1,-1
			acode%nc=j
			if (acode%jdif > acode%minint) RETURN
			acode%jdif=(acode%jdif+1)*acode%nrad-1
		end do
		call nrerror('NWK too small in arcode')
	else
		if (isign > 0) then
			if (ich > acode%nch .or. ich < 0) call nrerror('bad ich in arcode')
		else
			ja=ichar(codep(lcd))-acode%ilob(acode%nc)
			do j=acode%nc+1,NWK
				ja=ja*acode%nrad+(ichar(codep(j+lcd-acode%nc))-acode%ilob(j))
			end do
			ich=0
			ihi=acode%nch+1
			do
				if (ihi-ich <= 1) exit
				m=(ich+ihi)/2
				if (ja >= jtry(acode%jdif,acode%ncumfq(m+1),acode%ncum)) then
					ich=m
				else
					ihi=m
				end if
			end do
			if (ich == acode%nch) RETURN
		end if
		jh=jtry(acode%jdif,acode%ncumfq(ich+2),acode%ncum)
		jl=jtry(acode%jdif,acode%ncumfq(ich+1),acode%ncum)
		acode%jdif=jh-jl
		call arcsum(acode%ilob,acode%iupb,jh,NWK,acode%nrad,acode%nc)
		call arcsum(acode%ilob,acode%ilob,jl,NWK,acode%nrad,acode%nc)
		do j=acode%nc,NWK
			if (ich /= acode%nch .and. acode%iupb(j) /= acode%ilob(j)) exit
			if (acode%nc > size(codep)) codep=>reallocate(codep,2*size(codep))
			if (isign > 0) codep(lcd)=char(acode%ilob(j))
			lcd=lcd+1
		end do
		if (j > NWK) RETURN
		acode%nc=j
		j=0
		do
			if (acode%jdif >= acode%minint) exit
			j=j+1
			acode%jdif=acode%jdif*acode%nrad
		end do
		if (acode%nc-j < 1) call nrerror('NWK too small in arcode')
		if (j /= 0) then
			acode%iupb((acode%nc-j):(NWK-j))=acode%iupb(acode%nc:NWK)
			acode%ilob((acode%nc-j):(NWK-j))=acode%ilob(acode%nc:NWK)
		end if
		acode%nc=acode%nc-j
		acode%iupb((NWK-j+1):NWK)=0
		acode%ilob((NWK-j+1):NWK)=0
	end if
	CONTAINS
!BL
	FUNCTION jtry(m,n,k)
	use mo_nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: m,n,k
	INTEGER(I4B) :: jtry
	jtry=int((real(m,dp)*real(n,dp))/real(k,dp))
	END FUNCTION jtry
!BL
	SUBROUTINE arcsum(iin,iout,ja,nwk,nrad,nc)
	use mo_nrtype
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iin
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: iout
	INTEGER(I4B), INTENT(IN) :: nwk,nrad,nc
	INTEGER(I4B), INTENT(INOUT) :: ja
	INTEGER(I4B) :: j,jtmp,karry
	karry=0
	do j=nwk,nc+1,-1
		jtmp=ja
		ja=ja/nrad
		iout(j)=iin(j)+(jtmp-ja*nrad)+karry
		if (iout(j) >= nrad) then
			iout(j)=iout(j)-nrad
			karry=1
		else
			karry=0
		end if
	end do
	iout(nc)=iin(nc)+ja+karry
	END SUBROUTINE arcsum
	END SUBROUTINE arcode

! -----------------------------------------------------------

	SUBROUTINE asolve(b,x,itrnsp)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : sprsdiag
	USE xlinbcg_data
	REAL(DP), DIMENSION(:), INTENT(IN) :: b
	REAL(DP), DIMENSION(:), INTENT(OUT) :: x
	INTEGER(I4B), INTENT(IN) :: itrnsp
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(b),size(x),'asolve')
	call sprsdiag(sa,x)
	if (any(x == 0.0)) call nrerror('asolve: singular diagonal matrix')
	x=b/x
	END SUBROUTINE asolve

! -----------------------------------------------------------

	SUBROUTINE atimes(x,r,itrnsp)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : sprsax,sprstx
	USE xlinbcg_data
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(:), INTENT(OUT) :: r
	INTEGER(I4B), INTENT(IN) :: itrnsp
	INTEGER(I4B) :: n
	n=assert_eq(size(x),size(r),'atimes')
	if (itrnsp == 0) then
		call sprsax(sa,x,r)
	else
		call sprstx(sa,x,r)
	end if
	END SUBROUTINE atimes

! -----------------------------------------------------------

	SUBROUTINE avevar(data,ave,var)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data
	REAL(SP), INTENT(OUT) :: ave,var
	INTEGER(I4B) :: n
	REAL(SP), DIMENSION(size(data)) :: s
	n=size(data)
	ave=sum(data(:))/n
	s(:)=data(:)-ave
	var=dot_product(s,s)
	var=(var-sum(s)**2/n)/(n-1)
	END SUBROUTINE avevar

! -----------------------------------------------------------

	PROGRAM badluk
	use mo_nrtype
	use mo_nr, ONLY : flmoon,julday
	IMPLICIT NONE
	INTEGER(I4B) :: ic,icon,idwk,ifrac,im,iyyy,jd,jday,n
	INTEGER(I4B) :: iybeg=1900,iyend=2000
	REAL(SP) :: frac
	REAL(SP), PARAMETER :: TIMZON=-5.0_sp/24.0_sp
	write (*,'(1x,a,i5,a,i5)') 'Full moons on Friday the 13th from',&
		iybeg,' to',iyend
	do iyyy=iybeg,iyend
		do im=1,12
			jday=julday(im,13,iyyy)
			idwk=mod(jday+1,7)
			if (idwk == 5) then
				n=12.37_sp*(iyyy-1900+(im-0.5_sp)/12.0_sp)
				icon=0
				do
					call flmoon(n,2,jd,frac)
					ifrac=nint(24.0_sp*(frac+TIMZON))
					if (ifrac < 0) then
						jd=jd-1
						ifrac=ifrac+24
					end if
					if (ifrac > 12) then
						jd=jd+1
						ifrac=ifrac-12
					else
						ifrac=ifrac+12
					end if
					if (jd == jday) then
						write (*,'(/1x,i2,a,i2,a,i4)') im,'/',13,'/',iyyy
						write (*,'(1x,a,i2,a)') 'Full moon ',ifrac,&
							' hrs after midnight (EST).'
						exit
					else
						ic=isign(1,jday-jd)
						if (ic == -icon) exit
						icon=ic
						n=n+ic
					end if
				end do
			end if
		end do
	end do
	END PROGRAM badluk

! -----------------------------------------------------------

	SUBROUTINE balanc(a)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), PARAMETER :: RADX=radix(a),SQRADX=RADX**2
	INTEGER(I4B) :: i,last,ndum
	REAL(SP) :: c,f,g,r,s
	ndum=assert_eq(size(a,1),size(a,2),'balanc')
	do
		last=1
		do i=1,size(a,1)
			c=sum(abs(a(:,i)))-a(i,i)
			r=sum(abs(a(i,:)))-a(i,i)
			if (c /= 0.0 .and. r /= 0.0) then
				g=r/RADX
				f=1.0
				s=c+r
				do
					if (c >= g) exit
					f=f*RADX
					c=c*SQRADX
				end do
				g=r*RADX
				do
					if (c <= g) exit
					f=f/RADX
					c=c/SQRADX
				end do
				if ((c+r)/f < 0.95_sp*s) then
					last=0
					g=1.0_sp/f
					a(i,:)=a(i,:)*g
					a(:,i)=a(:,i)*f
				end if
			end if
		end do
		if (last /= 0) exit
	end do
	END SUBROUTINE balanc

! -----------------------------------------------------------

	SUBROUTINE banbks(a,m1,m2,al,indx,b)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,al
	INTEGER(I4B), INTENT(IN) :: m1,m2
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: i,k,l,mdum,mm,n
	n=assert_eq(size(a,1),size(al,1),size(b),size(indx),'banbks: n')
	mm=assert_eq(size(a,2),m1+m2+1,'banbks: mm')
	mdum=assert_eq(size(al,2),m1,'banbks: mdum')
	do k=1,n
		l=min(n,m1+k)
		i=indx(k)
		if (i /= k) call swap(b(i),b(k))
		b(k+1:l)=b(k+1:l)-al(k,1:l-k)*b(k)
	end do
	do i=n,1,-1
		l=min(mm,n-i+1)
		b(i)=(b(i)-dot_product(a(i,2:l),b(1+i:i+l-1)))/a(i,1)
	end do
	END SUBROUTINE banbks

! -----------------------------------------------------------

	SUBROUTINE bandec(a,m1,m2,al,indx,d)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,imaxloc,swap,arth
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B), INTENT(IN) :: m1,m2
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: al
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
	REAL(SP), INTENT(OUT) :: d
	REAL(SP), PARAMETER :: TINY=1.0e-20_sp
	INTEGER(I4B) :: i,k,l,mdum,mm,n
	REAL(SP) :: dum
	n=assert_eq(size(a,1),size(al,1),size(indx),'bandec: n')
	mm=assert_eq(size(a,2),m1+m2+1,'bandec: mm')
	mdum=assert_eq(size(al,2),m1,'bandec: mdum')
	a(1:m1,:)=eoshift(a(1:m1,:),dim=2,shift=arth(m1,-1,m1))
	d=1.0
	do k=1,n
		l=min(m1+k,n)
		i=imaxloc(abs(a(k:l,1)))+k-1
		dum=a(i,1)
		if (dum == 0.0) a(k,1)=TINY
		indx(k)=i
		if (i /= k) then
			d=-d
			call swap(a(k,1:mm),a(i,1:mm))
		end if
		do i=k+1,l
			dum=a(i,1)/a(k,1)
			al(k,i-k)=dum
			a(i,1:mm-1)=a(i,2:mm)-dum*a(k,2:mm)
			a(i,mm)=0.0
		end do
	end do
	END SUBROUTINE bandec

! -----------------------------------------------------------

	SUBROUTINE banmul(a,m1,m2,x,b)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,arth
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	INTEGER(I4B), INTENT(IN) :: m1,m2
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(OUT) :: b
	INTEGER(I4B) :: m,n
	n=assert_eq(size(a,1),size(b),size(x),'banmul: n')
	m=assert_eq(size(a,2),m1+m2+1,'banmul: m')
	b=sum(a*eoshift(spread(x,dim=2,ncopies=m), &
		dim=1,shift=arth(-m1,1,m)),dim=2)
	END SUBROUTINE banmul

! -----------------------------------------------------------

	SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: d1,d2
	REAL(SP), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
	REAL(SP), DIMENSION(4,4), INTENT(OUT) :: c
	REAL(SP), DIMENSION(16) :: x
	REAL(SP), DIMENSION(16,16) :: wt
	DATA wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
		8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,&
		2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,&
		2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0,&
		-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,&
		-1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,&
		-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
	x(1:4)=y
	x(5:8)=y1*d1
	x(9:12)=y2*d2
	x(13:16)=y12*d1*d2
	x=matmul(wt,x)
	c=reshape(x,(/4,4/),order=(/2,1/))
	END SUBROUTINE bcucof

! -----------------------------------------------------------

	SUBROUTINE bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : bcucof
	IMPLICIT NONE
	REAL(SP), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
	REAL(SP), INTENT(IN) :: x1l,x1u,x2l,x2u,x1,x2
	REAL(SP), INTENT(OUT) :: ansy,ansy1,ansy2
	INTEGER(I4B) :: i
	REAL(SP) :: t,u
	REAL(SP), DIMENSION(4,4) :: c
	call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
	if (x1u == x1l .or. x2u == x2l) call &
		nrerror('bcuint: problem with input values - boundary pair equal?')
	t=(x1-x1l)/(x1u-x1l)
	u=(x2-x2l)/(x2u-x2l)
	ansy=0.0
	ansy2=0.0
	ansy1=0.0
	do i=4,1,-1
		ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
		ansy2=t*ansy2+(3.0_sp*c(i,4)*u+2.0_sp*c(i,3))*u+c(i,2)
		ansy1=u*ansy1+(3.0_sp*c(4,i)*t+2.0_sp*c(3,i))*t+c(2,i)
	end do
	ansy1=ansy1/(x1u-x1l)
	ansy2=ansy2/(x2u-x2l)
	END SUBROUTINE bcuint

! -----------------------------------------------------------

	SUBROUTINE beschb_s(x,gam1,gam2,gampl,gammi)
	use mo_nrtype
	use mo_nr, ONLY : chebev
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	REAL(DP), INTENT(OUT) :: gam1,gam2,gampl,gammi
	INTEGER(I4B), PARAMETER :: NUSE1=5,NUSE2=5
	REAL(SP) :: xx
	REAL(SP), DIMENSION(7) :: c1=(/-1.142022680371168_sp,&
		6.5165112670737e-3_sp,3.087090173086e-4_sp,-3.4706269649e-6_sp,&
		6.9437664e-9_sp,3.67795e-11_sp,-1.356e-13_sp/)
	REAL(SP), DIMENSION(8) :: c2=(/1.843740587300905_sp,&
		-7.68528408447867e-2_sp,1.2719271366546e-3_sp,&
		-4.9717367042e-6_sp, -3.31261198e-8_sp,2.423096e-10_sp,&
		-1.702e-13_sp,-1.49e-15_sp/)
	xx=8.0_dp*x*x-1.0_dp
	gam1=chebev(-1.0_sp,1.0_sp,c1(1:NUSE1),xx)
	gam2=chebev(-1.0_sp,1.0_sp,c2(1:NUSE2),xx)
	gampl=gam2-x*gam1
	gammi=gam2+x*gam1
	END SUBROUTINE beschb_s


	SUBROUTINE beschb_v(x,gam1,gam2,gampl,gammi)
	use mo_nrtype
	use mo_nr, ONLY : chebev
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(:), INTENT(OUT) :: gam1,gam2,gampl,gammi
	INTEGER(I4B), PARAMETER :: NUSE1=5,NUSE2=5
	REAL(SP), DIMENSION(size(x)) :: xx
	REAL(SP), DIMENSION(7) :: c1=(/-1.142022680371168_sp,&
		6.5165112670737e-3_sp,3.087090173086e-4_sp,-3.4706269649e-6_sp,&
		6.9437664e-9_sp,3.67795e-11_sp,-1.356e-13_sp/)
	REAL(SP), DIMENSION(8) :: c2=(/1.843740587300905_sp,&
		-7.68528408447867e-2_sp,1.2719271366546e-3_sp,&
		-4.9717367042e-6_sp, -3.31261198e-8_sp,2.423096e-10_sp,&
		-1.702e-13_sp,-1.49e-15_sp/)
	xx=8.0_dp*x*x-1.0_dp
	gam1=chebev(-1.0_sp,1.0_sp,c1(1:NUSE1),xx)
	gam2=chebev(-1.0_sp,1.0_sp,c2(1:NUSE2),xx)
	gampl=gam2-x*gam1
	gammi=gam2+x*gam1
	END SUBROUTINE beschb_v

! -----------------------------------------------------------

	FUNCTION bessi_s(n,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessi0
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessi_s
	INTEGER(I4B), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
	INTEGER(I4B) :: j,m
	REAL(SP) :: bi,bim,bip,tox
	call assert(n >= 2, 'bessi_s args')
	bessi_s=0.0
	if (x*x <= 8.0_sp*tiny(x)) RETURN
	tox=2.0_sp/abs(x)
	bip=0.0
	bi=1.0
	m=2*((n+int(sqrt(real(IACC*n,sp)))))
	do j=m,1,-1
		bim=bip+j*tox*bi
		bip=bi
		bi=bim
		if (exponent(bi) > IEXP) then
			bessi_s=scale(bessi_s,-IEXP)
			bi=scale(bi,-IEXP)
			bip=scale(bip,-IEXP)
		end if
		if (j == n) bessi_s=bip
	end do
	bessi_s=bessi_s*bessi0(x)/bi
	if (x < 0.0 .and. mod(n,2) == 1) bessi_s=-bessi_s
	END FUNCTION bessi_s


	FUNCTION bessi_v(n,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessi0
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessi_v
	INTEGER(I4B), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
	INTEGER(I4B) :: j,m
	REAL(SP), DIMENSION(size(x)) :: bi,bim,bip,tox
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	call assert(n >= 2, 'bessi_v args')
	bessi_v=0.0
	mask = (x <= 8.0_sp*tiny(x))
	tox=2.0_sp/merge(2.0_sp,abs(x),mask)
	bip=0.0
	bi=1.0_sp
	m=2*((n+int(sqrt(real(IACC*n,sp)))))
	do j=m,1,-1
		bim=bip+j*tox*bi
		bip=bi
		bi=bim
		where (exponent(bi) > IEXP)
			bessi_v=scale(bessi_v,-IEXP)
			bi=scale(bi,-IEXP)
			bip=scale(bip,-IEXP)
		end where
		if (j == n) bessi_v=bip
	end do
	bessi_v=bessi_v*bessi0(x)/bi
	where (mask) bessi_v=0.0_sp
	where (x < 0.0 .and. mod(n,2) == 1) bessi_v=-bessi_v
	END FUNCTION bessi_v

! -----------------------------------------------------------

	FUNCTION bessi0_s(x)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessi0_s
	REAL(SP) :: ax
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,3.5156229_dp,&
		3.0899424_dp,1.2067492_dp,0.2659732_dp,0.360768e-1_dp,&
		0.45813e-2_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,0.1328592e-1_dp,&
		0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,&
		-0.2057706e-1_dp,0.2635537e-1_dp,-0.1647633e-1_dp,&
		0.392377e-2_dp/)
	ax=abs(x)
	if (ax < 3.75) then
		bessi0_s=poly(real((x/3.75_sp)**2,dp),p)
	else
		bessi0_s=(exp(ax)/sqrt(ax))*poly(real(3.75_sp/ax,dp),q)
	end if
	END FUNCTION bessi0_s


	FUNCTION bessi0_v(x)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessi0_v
	REAL(SP), DIMENSION(size(x)) :: ax
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,3.5156229_dp,&
		3.0899424_dp,1.2067492_dp,0.2659732_dp,0.360768e-1_dp,&
		0.45813e-2_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,0.1328592e-1_dp,&
		0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,&
		-0.2057706e-1_dp,0.2635537e-1_dp,-0.1647633e-1_dp,&
		0.392377e-2_dp/)
	ax=abs(x)
	mask = (ax < 3.75)
	where (mask)
		bessi0_v=poly(real((x/3.75_sp)**2,dp),p,mask)
	elsewhere
		y=3.75_sp/ax
		bessi0_v=(exp(ax)/sqrt(ax))*poly(real(y,dp),q,.not. mask)
	end where
	END FUNCTION bessi0_v

! -----------------------------------------------------------

	FUNCTION bessi1_s(x)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessi1_s
	REAL(SP) :: ax
	REAL(DP), DIMENSION(7) :: p = (/0.5_dp,0.87890594_dp,&
		0.51498869_dp,0.15084934_dp,0.2658733e-1_dp,&
		0.301532e-2_dp,0.32411e-3_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,-0.3988024e-1_dp,&
		-0.362018e-2_dp,0.163801e-2_dp,-0.1031555e-1_dp,&
		0.2282967e-1_dp,-0.2895312e-1_dp,0.1787654e-1_dp,&
		-0.420059e-2_dp/)
	ax=abs(x)
	if (ax < 3.75) then
		bessi1_s=ax*poly(real((x/3.75_sp)**2,dp),p)
	else
		bessi1_s=(exp(ax)/sqrt(ax))*poly(real(3.75_sp/ax,dp),q)
	end if
	if (x < 0.0) bessi1_s=-bessi1_s
	END FUNCTION bessi1_s


	FUNCTION bessi1_v(x)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessi1_v
	REAL(SP), DIMENSION(size(x)) :: ax
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(7) :: p = (/0.5_dp,0.87890594_dp,&
		0.51498869_dp,0.15084934_dp,0.2658733e-1_dp,&
		0.301532e-2_dp,0.32411e-3_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,-0.3988024e-1_dp,&
		-0.362018e-2_dp,0.163801e-2_dp,-0.1031555e-1_dp,&
		0.2282967e-1_dp,-0.2895312e-1_dp,0.1787654e-1_dp,&
		-0.420059e-2_dp/)
	ax=abs(x)
	mask = (ax < 3.75)
	where (mask)
		bessi1_v=ax*poly(real((x/3.75_sp)**2,dp),p,mask)
	elsewhere
		y=3.75_sp/ax
		bessi1_v=(exp(ax)/sqrt(ax))*poly(real(y,dp),q,.not. mask)
	end where
	where (x < 0.0) bessi1_v=-bessi1_v
	END FUNCTION bessi1_v

! -----------------------------------------------------------

	SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,nrerror
	use mo_nr, ONLY : beschb
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,xnu
	REAL(SP), INTENT(OUT) :: ri,rk,rip,rkp
	INTEGER(I4B), PARAMETER :: MAXIT=10000
	REAL(SP), PARAMETER :: XMIN=2.0
	REAL(DP), PARAMETER :: EPS=1.0e-10_dp,FPMIN=1.0e-30_dp
	INTEGER(I4B) :: i,l,nl
	REAL(DP) :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,&
		gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,&
		ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup,rktemp,&
		s,sum,sum1,x2,xi,xi2,xmu,xmu2
	call assert(x > 0.0, xnu >= 0.0, 'bessik args')
	nl=int(xnu+0.5_dp)
	xmu=xnu-nl
	xmu2=xmu*xmu
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	h=xnu*xi
	if (h < FPMIN) h=FPMIN
	b=xi2*xnu
	d=0.0
	c=h
	do i=1,MAXIT
		b=b+xi2
		d=1.0_dp/(b+d)
		c=b+1.0_dp/c
		del=c*d
		h=del*h
		if (abs(del-1.0_dp) < EPS) exit
	end do
	if (i > MAXIT) call nrerror('x too large in bessik; try asymptotic expansion')
	ril=FPMIN
	ripl=h*ril
	ril1=ril
	rip1=ripl
	fact=xnu*xi
	do l=nl,1,-1
		ritemp=fact*ril+ripl
		fact=fact-xi
		ripl=fact*ritemp+ril
		ril=ritemp
	end do
	f=ripl/ril
	if (x < XMIN) then
		x2=0.5_dp*x
		pimu=PI_D*xmu
		if (abs(pimu) < EPS) then
			fact=1.0
		else
			fact=pimu/sin(pimu)
		end if
		d=-log(x2)
		e=xmu*d
		if (abs(e) < EPS) then
			fact2=1.0
		else
			fact2=sinh(e)/e
		end if
		call beschb(xmu,gam1,gam2,gampl,gammi)
		ff=fact*(gam1*cosh(e)+gam2*fact2*d)
		sum=ff
		e=exp(e)
		p=0.5_dp*e/gampl
		q=0.5_dp/(e*gammi)
		c=1.0
		d=x2*x2
		sum1=p
		do i=1,MAXIT
			ff=(i*ff+p+q)/(i*i-xmu2)
			c=c*d/i
			p=p/(i-xmu)
			q=q/(i+xmu)
			del=c*ff
			sum=sum+del
			del1=c*(p-i*ff)
			sum1=sum1+del1
			if (abs(del) < abs(sum)*EPS) exit
		end do
		if (i > MAXIT) call nrerror('bessk series failed to converge')
		rkmu=sum
		rk1=sum1*xi2
	else
		b=2.0_dp*(1.0_dp+x)
		d=1.0_dp/b
		delh=d
		h=delh
		q1=0.0
		q2=1.0
		a1=0.25_dp-xmu2
		c=a1
		q=c
		a=-a1
		s=1.0_dp+q*delh
		do i=2,MAXIT
			a=a-2*(i-1)
			c=-a*c/i
			qnew=(q1-b*q2)/a
			q1=q2
			q2=qnew
			q=q+c*qnew
			b=b+2.0_dp
			d=1.0_dp/(b+a*d)
			delh=(b*d-1.0_dp)*delh
			h=h+delh
			dels=q*delh
			s=s+dels
			if (abs(dels/s) < EPS) exit
		end do
		if (i > MAXIT) call nrerror('bessik: failure to converge in cf2')
		h=a1*h
		rkmu=sqrt(PI_D/(2.0_dp*x))*exp(-x)/s
		rk1=rkmu*(xmu+x+0.5_dp-h)*xi
	end if
	rkmup=xmu*xi*rkmu-rk1
	rimu=xi/(f*rkmu-rkmup)
	ri=(rimu*ril1)/ril
	rip=(rimu*rip1)/ril
	do i=1,nl
		rktemp=(xmu+i)*xi2*rk1+rkmu
		rkmu=rk1
		rk1=rktemp
	end do
	rk=rkmu
	rkp=xnu*xi*rkmu-rk1
	END SUBROUTINE bessik

! -----------------------------------------------------------

	FUNCTION bessj_s(n,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessj0,bessj1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessj_s
	INTEGER(I4B), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
	INTEGER(I4B) :: j,jsum,m
	REAL(SP) :: ax,bj,bjm,bjp,summ,tox
	call assert(n >= 2, 'bessj_s args')
	ax=abs(x)
	if (ax*ax <= 8.0_sp*tiny(x)) then
		bessj_s=0.0
	else if (ax > real(n,sp)) then
		tox=2.0_sp/ax
		bjm=bessj0(ax)
		bj=bessj1(ax)
		do j=1,n-1
			bjp=j*tox*bj-bjm
			bjm=bj
			bj=bjp
		end do
		bessj_s=bj
	else
		tox=2.0_sp/ax
		m=2*((n+int(sqrt(real(IACC*n,sp))))/2)
		bessj_s=0.0
		jsum=0
		summ=0.0
		bjp=0.0
		bj=1.0
		do j=m,1,-1
			bjm=j*tox*bj-bjp
			bjp=bj
			bj=bjm
			if (exponent(bj) > IEXP) then
				bj=scale(bj,-IEXP)
				bjp=scale(bjp,-IEXP)
				bessj_s=scale(bessj_s,-IEXP)
				summ=scale(summ,-IEXP)
			end if
			if (jsum /= 0) summ=summ+bj
			jsum=1-jsum
			if (j == n) bessj_s=bjp
		end do
		summ=2.0_sp*summ-bj
		bessj_s=bessj_s/summ
	end if
	if (x < 0.0 .and. mod(n,2) == 1) bessj_s=-bessj_s
	END FUNCTION bessj_s

	FUNCTION bessj_v(n,xx)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessj0,bessj1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	REAL(SP), DIMENSION(size(xx)) :: bessj_v
	INTEGER(I4B), PARAMETER :: IACC=40,IEXP=maxexponent(xx)/2
	REAL(SP), DIMENSION(size(xx)) :: ax
	LOGICAL(LGT), DIMENSION(size(xx)) :: mask,mask0
	REAL(SP), DIMENSION(:), ALLOCATABLE :: x,bj,bjm,bjp,summ,tox,bessjle
	LOGICAL(LGT), DIMENSION(:), ALLOCATABLE :: renorm
	INTEGER(I4B) :: j,jsum,m,npak
	call assert(n >= 2, 'bessj_v args')
	ax=abs(xx)
	mask = (ax <= real(n,sp))
	mask0 = (ax*ax <= 8.0_sp*tiny(xx))
	bessj_v=bessjle_v(n,ax,logical(mask .and. .not.mask0, kind=lgt))
	bessj_v=merge(bessjgt_v(n,ax,.not. mask),bessj_v,.not. mask)
	where (mask0) bessj_v=0.0
	where (xx < 0.0 .and. mod(n,2) == 1) bessj_v=-bessj_v
	CONTAINS
!BL
!BL
	FUNCTION bessjgt_v(n,xx,mask)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	LOGICAL(LGT), DIMENSION(size(xx)), INTENT(IN) :: mask
	REAL(SP), DIMENSION(size(xx)) :: bessjgt_v
	npak=count(mask)
	if (npak == 0) RETURN
	allocate(x(npak),bj(npak),bjm(npak),bjp(npak),tox(npak))
	x=pack(xx,mask)
	tox=2.0_sp/x
	bjm=bessj0(x)
	bj=bessj1(x)
	do j=1,n-1
		bjp=j*tox*bj-bjm
		bjm=bj
		bj=bjp
	end do
	bessjgt_v=unpack(bj,mask,0.0_sp)
	deallocate(x,bj,bjm,bjp,tox)
	END FUNCTION bessjgt_v
!BL
	FUNCTION bessjle_v(n,xx,mask)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	LOGICAL(LGT), DIMENSION(size(xx)), INTENT(IN) :: mask
	REAL(SP), DIMENSION(size(xx)) :: bessjle_v
	npak=count(mask)
	if (npak == 0) RETURN
	allocate(x(npak),bj(npak),bjm(npak),bjp(npak),summ(npak), &
		bessjle(npak),tox(npak),renorm(npak))
	x=pack(xx,mask)
	tox=2.0_sp/x
	m=2*((n+int(sqrt(real(IACC*n,sp))))/2)
	bessjle=0.0
	jsum=0
	summ=0.0
	bjp=0.0
	bj=1.0
	do j=m,1,-1
		bjm=j*tox*bj-bjp
		bjp=bj
		bj=bjm
		renorm = (exponent(bj)>IEXP)
		bj=merge(scale(bj,-IEXP),bj,renorm)
		bjp=merge(scale(bjp,-IEXP),bjp,renorm)
		bessjle=merge(scale(bessjle,-IEXP),bessjle,renorm)
		summ=merge(scale(summ,-IEXP),summ,renorm)
		if (jsum /= 0) summ=summ+bj
		jsum=1-jsum
		if (j == n) bessjle=bjp
	end do
	summ=2.0_sp*summ-bj
	bessjle=bessjle/summ
	bessjle_v=unpack(bessjle,mask,0.0_sp)
	deallocate(x,bj,bjm,bjp,summ,bessjle,tox,renorm)
	END FUNCTION bessjle_v
	END FUNCTION bessj_v

! -----------------------------------------------------------

	FUNCTION bessj0_s(x)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessj0_s
	REAL(SP) :: ax,xx,z
	REAL(DP) :: y
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
		0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/-0.1562499995e-1_dp,&
		0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
		-0.934945152e-7_dp/)
	REAL(DP), DIMENSION(6) :: r = (/57568490574.0_dp,-13362590354.0_dp,&
		651619640.7_dp,-11214424.18_dp,77392.33017_dp,&
		-184.9052456_dp/)
	REAL(DP), DIMENSION(6) :: s = (/57568490411.0_dp,1029532985.0_dp,&
		9494680.718_dp,59272.64853_dp,267.8532712_dp,1.0_dp/)
	if (abs(x) < 8.0) then
		y=x**2
		bessj0_s=poly(y,r)/poly(y,s)
	else
		ax=abs(x)
		z=8.0_sp/ax
		y=z**2
		xx=ax-0.785398164_sp
		bessj0_s=sqrt(0.636619772_sp/ax)*(cos(xx)*&
			poly(y,p)-z*sin(xx)*poly(y,q))
	end if
	END FUNCTION bessj0_s



	FUNCTION bessj0_v(x)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessj0_v
	REAL(SP), DIMENSION(size(x)) :: ax,xx,z
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
		0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/-0.1562499995e-1_dp,&
		0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
		-0.934945152e-7_dp/)
	REAL(DP), DIMENSION(6) :: r = (/57568490574.0_dp,-13362590354.0_dp,&
		651619640.7_dp,-11214424.18_dp,77392.33017_dp,&
		-184.9052456_dp/)
	REAL(DP), DIMENSION(6) :: s = (/57568490411.0_dp,1029532985.0_dp,&
		9494680.718_dp,59272.64853_dp,267.8532712_dp,1.0_dp/)
	mask = (abs(x) < 8.0)
	where (mask)
		y=x**2
		bessj0_v=poly(y,r,mask)/poly(y,s,mask)
	elsewhere
		ax=abs(x)
		z=8.0_sp/ax
		y=z**2
		xx=ax-0.785398164_sp
		bessj0_v=sqrt(0.636619772_sp/ax)*(cos(xx)*&
			poly(y,p,.not. mask)-z*sin(xx)*poly(y,q,.not. mask))
	end where
	END FUNCTION bessj0_v

! -----------------------------------------------------------

	FUNCTION bessj1_s(x)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessj1_s
	REAL(SP) :: ax,xx,z
	REAL(DP) :: y
	REAL(DP), DIMENSION(6) :: r = (/72362614232.0_dp,&
		-7895059235.0_dp,242396853.1_dp,-2972611.439_dp,&
		15704.48260_dp,-30.16036606_dp/)
	REAL(DP), DIMENSION(6) :: s = (/144725228442.0_dp,2300535178.0_dp,&
		18583304.74_dp,99447.43394_dp,376.9991397_dp,1.0_dp/)
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,0.183105e-2_dp,&
		-0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/0.04687499995_dp,&
		-0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
		0.105787412e-6_dp/)
	if (abs(x) < 8.0) then
		y=x**2
		bessj1_s=x*(poly(y,r)/poly(y,s))
	else
		ax=abs(x)
		z=8.0_sp/ax
		y=z**2
		xx=ax-2.356194491_sp
		bessj1_s=sqrt(0.636619772_sp/ax)*(cos(xx)*&
			poly(y,p)-z*sin(xx)*poly(y,q))*sign(1.0_sp,x)
	end if
	END FUNCTION bessj1_s


	FUNCTION bessj1_v(x)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessj1_v
	REAL(SP), DIMENSION(size(x)) :: ax,xx,z
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(6) :: r = (/72362614232.0_dp,&
		-7895059235.0_dp,242396853.1_dp,-2972611.439_dp,&
		15704.48260_dp,-30.16036606_dp/)
	REAL(DP), DIMENSION(6) :: s = (/144725228442.0_dp,2300535178.0_dp,&
		18583304.74_dp,99447.43394_dp,376.9991397_dp,1.0_dp/)
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,0.183105e-2_dp,&
		-0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/0.04687499995_dp,&
		-0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
		0.105787412e-6_dp/)
	mask = (abs(x) < 8.0)
	where (mask)
		y=x**2
		bessj1_v=x*(poly(y,r,mask)/poly(y,s,mask))
	elsewhere
		ax=abs(x)
		z=8.0_sp/ax
		y=z**2
		xx=ax-2.356194491_sp
		bessj1_v=sqrt(0.636619772_sp/ax)*(cos(xx)*&
			poly(y,p,.not. mask)-z*sin(xx)*poly(y,q,.not. mask))*&
			sign(1.0_sp,x)
	end where
	END FUNCTION bessj1_v

! -----------------------------------------------------------

	SUBROUTINE bessjy_s(x,xnu,rj,ry,rjp,ryp)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,nrerror
	use mo_nr, ONLY : beschb
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,xnu
	REAL(SP), INTENT(OUT) :: rj,ry,rjp,ryp
	INTEGER(I4B), PARAMETER :: MAXIT=10000
	REAL(DP), PARAMETER :: XMIN=2.0_dp,EPS=1.0e-10_dp,FPMIN=1.0e-30_dp
	INTEGER(I4B) :: i,isign,l,nl
	REAL(DP) :: a,b,c,d,del,del1,e,f,fact,fact2,fact3,ff,gam,gam1,gam2,&
		gammi,gampl,h,p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,&
		ry1,rymu,rymup,rytemp,sum,sum1,w,x2,xi,xi2,xmu,xmu2
	COMPLEX(DPC) :: aa,bb,cc,dd,dl,pq
	call assert(x > 0.0, xnu >= 0.0, 'bessjy args')
	nl=merge(int(xnu+0.5_dp), max(0,int(xnu-x+1.5_dp)), x < XMIN)
	xmu=xnu-nl
	xmu2=xmu*xmu
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	w=xi2/PI_D
	isign=1
	h=xnu*xi
	if (h < FPMIN) h=FPMIN
	b=xi2*xnu
	d=0.0
	c=h
	do i=1,MAXIT
		b=b+xi2
		d=b-d
		if (abs(d) < FPMIN) d=FPMIN
		c=b-1.0_dp/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_dp/d
		del=c*d
		h=del*h
		if (d < 0.0) isign=-isign
		if (abs(del-1.0_dp) < EPS) exit
	end do
	if (i > MAXIT) call nrerror('x too large in bessjy; try asymptotic expansion')
	rjl=isign*FPMIN
	rjpl=h*rjl
	rjl1=rjl
	rjp1=rjpl
	fact=xnu*xi
	do l=nl,1,-1
		rjtemp=fact*rjl+rjpl
		fact=fact-xi
		rjpl=fact*rjtemp-rjl
		rjl=rjtemp
	end do
	if (rjl == 0.0) rjl=EPS
	f=rjpl/rjl
	if (x < XMIN) then
		x2=0.5_dp*x
		pimu=PI_D*xmu
		if (abs(pimu) < EPS) then
			fact=1.0
		else
			fact=pimu/sin(pimu)
		end if
		d=-log(x2)
		e=xmu*d
		if (abs(e) < EPS) then
			fact2=1.0
		else
			fact2=sinh(e)/e
		end if
		call beschb(xmu,gam1,gam2,gampl,gammi)
		ff=2.0_dp/PI_D*fact*(gam1*cosh(e)+gam2*fact2*d)
		e=exp(e)
		p=e/(gampl*PI_D)
		q=1.0_dp/(e*PI_D*gammi)
		pimu2=0.5_dp*pimu
		if (abs(pimu2) < EPS) then
			fact3=1.0
		else
			fact3=sin(pimu2)/pimu2
		end if
		r=PI_D*pimu2*fact3*fact3
		c=1.0
		d=-x2*x2
		sum=ff+r*q
		sum1=p
		do i=1,MAXIT
			ff=(i*ff+p+q)/(i*i-xmu2)
			c=c*d/i
			p=p/(i-xmu)
			q=q/(i+xmu)
			del=c*(ff+r*q)
			sum=sum+del
			del1=c*p-i*del
			sum1=sum1+del1
			if (abs(del) < (1.0_dp+abs(sum))*EPS) exit
		end do
		if (i > MAXIT) call nrerror('bessy series failed to converge')
		rymu=-sum
		ry1=-sum1*xi2
		rymup=xmu*xi*rymu-ry1
		rjmu=w/(rymup-f*rymu)
	else
		a=0.25_dp-xmu2
		pq=cmplx(-0.5_dp*xi,1.0_dp,kind=dpc)
		aa=cmplx(0.0_dp,xi*a,kind=dpc)
		bb=cmplx(2.0_dp*x,2.0_dp,kind=dpc)
		cc=bb+aa/pq
		dd=1.0_dp/bb
		pq=cc*dd*pq
		do i=2,MAXIT
			a=a+2*(i-1)
			bb=bb+cmplx(0.0_dp,2.0_dp,kind=dpc)
			dd=a*dd+bb
			if (absc(dd) < FPMIN) dd=FPMIN
			cc=bb+a/cc
			if (absc(cc) < FPMIN) cc=FPMIN
			dd=1.0_dp/dd
			dl=cc*dd
			pq=pq*dl
			if (absc(dl-1.0_dp) < EPS) exit
		end do
		if (i > MAXIT) call nrerror('cf2 failed in bessjy')
		p=real(pq)
		q=aimag(pq)
		gam=(p-f)/q
		rjmu=sqrt(w/((p-f)*gam+q))
		rjmu=sign(rjmu,rjl)
		rymu=rjmu*gam
		rymup=rymu*(p+q/gam)
		ry1=xmu*xi*rymu-rymup
	end if
	fact=rjmu/rjl
	rj=rjl1*fact
	rjp=rjp1*fact
	do i=1,nl
		rytemp=(xmu+i)*xi2*ry1-rymu
		rymu=ry1
		ry1=rytemp
	end do
	ry=rymu
	ryp=xnu*xi*rymu-ry1
	CONTAINS
!BL
	FUNCTION absc(z)
	IMPLICIT NONE
	COMPLEX(DPC), INTENT(IN) :: z
	REAL(DP) :: absc
	absc=abs(real(z))+abs(aimag(z))
	END FUNCTION absc
	END SUBROUTINE bessjy_s

	SUBROUTINE bessjy_v(x,xnu,rj,ry,rjp,ryp)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: xnu
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
	INTEGER(I4B), PARAMETER :: MAXIT=10000
	REAL(DP), PARAMETER :: XMIN=2.0_dp,EPS=1.0e-10_dp,FPMIN=1.0e-30_dp
	INTEGER(I4B) :: i,n
	INTEGER(I4B), DIMENSION(size(x)) :: isign
	REAL(DP), DIMENSION(size(x)) :: b,c,d,del,&
		h,xi,xi2,&
		rj_lt,ry_lt,rjp_lt,ryp_lt,rj_ge,ry_ge,rjp_ge,ryp_ge
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	call assert(all(x > 0.0), xnu >= 0.0, 'bessjy args')
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	isign=1
	h=xnu*xi
	where (h < FPMIN) h=FPMIN
	b=xi2*xnu
	d=0.0
	c=h
	converged=.false.
	do i=1,MAXIT
		where (.not. converged)
			b=b+xi2
			d=b-d
			d=merge(FPMIN,d, abs(d) < FPMIN )
			c=b-1.0_dp/c
			c=merge(FPMIN,c, abs(c) < FPMIN )
			d=1.0_dp/d
			del=c*d
			h=del*h
			isign=merge(-isign,isign, d < 0.0 )
			converged=(abs(del-1.0_dp) < EPS)
		end where
		if (all(converged)) exit
	end do
	if (i > MAXIT) call nrerror('bessjy: x too large; try asymptotic expansion')
	converged=(x < XMIN)
	n=count(converged)
	call bessjy_xltxmin(pack(x,converged),xnu,pack(isign*FPMIN,converged),pack(h,converged),&
		rj_lt(1:n),ry_lt(1:n),rjp_lt(1:n),ryp_lt(1:n))
	n=size(x)-n
	call bessjy_xgexmin(pack(x,.not. converged),xnu,pack(isign,.not. converged),pack(h,.not. converged),&
		rj_ge(1:n),ry_ge(1:n),rjp_ge(1:n),ryp_ge(1:n))
	rj=unpacked(converged,rj_lt,rj_ge)
	rjp=unpacked(converged,rjp_lt,rjp_ge)
	ry=unpacked(converged,ry_lt,ry_ge)
	ryp=unpacked(converged,ryp_lt,ryp_ge)
	CONTAINS
!BL
	SUBROUTINE bessjy_xltxmin(x,xnu,rjli,hi,rj,ry,rjp,ryp)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : beschb
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), INTENT(IN) :: xnu
	REAL(DP), DIMENSION(:), INTENT(IN) :: hi,rjli
	REAL(DP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
	INTEGER(I4B) :: i,nl
	REAL(DP) :: xmu,xmu2,gam1,gam2,gammi,gampl
	REAL(DP), DIMENSION(size(x)) :: &
		c,d,del,del1,e,f,fact,fact2,fact3,ff,h,&
		p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjpl,rjp1,rjtemp,&
		ry1,rymu,rymup,rytemp,sum,sum1,w,x2,xi,xi2
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	nl=int(xnu+0.5_dp)
	xmu=xnu-nl
	xmu2=xmu*xmu
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	w=xi2/PI_D
	h=hi
	where (h < FPMIN) h=FPMIN
	rjl=rjli
	rjpl=h*rjl
	rjl1=rjl
	rjp1=rjpl
	fact=xnu*xi
	do i=int(xnu+0.5_dp),1,-1
		rjtemp=fact*rjl+rjpl
		fact=fact-xi
		rjpl=fact*rjtemp-rjl
		rjl=rjtemp
	end do
	where (rjl == 0.0) rjl=EPS
	f=rjpl/rjl
	x2=0.5_dp*x
	pimu=PI_D*xmu
	where (abs(pimu) < EPS)
		fact=1.0
	elsewhere
		fact=pimu/sin(pimu)
	end where
	d=-log(x2)
	e=xmu*d
	where (abs(e) < EPS)
		fact2=1.0
	elsewhere
		fact2=sinh(e)/e
	end where
	call beschb(xmu,gam1,gam2,gampl,gammi)
	ff=2.0_dp/PI_D*fact*(gam1*cosh(e)+gam2*fact2*d)
	e=exp(e)
	p=e/(gampl*PI_D)
	q=1.0_dp/(e*PI_D*gammi)
	pimu2=0.5_dp*pimu
	where (abs(pimu2) < EPS)
		fact3=1.0
	elsewhere
		fact3=sin(pimu2)/pimu2
	end where
	r=PI_D*pimu2*fact3*fact3
	c=1.0
	d=-x2*x2
	sum=ff+r*q
	sum1=p
	converged=.false.
	do i=1,MAXIT
		where (.not. converged)
			ff=(i*ff+p+q)/(i*i-xmu2)
			c=c*d/i
			p=p/(i-xmu)
			q=q/(i+xmu)
			del=c*(ff+r*q)
			sum=sum+del
			del1=c*p-i*del
			sum1=sum1+del1
			converged=(abs(del) < (1.0_dp+abs(sum))*EPS)
		end where
		if (all(converged)) exit
	end do
	if (i > MAXIT) call nrerror('bessy series failed to converge')
	rymu=-sum
	ry1=-sum1*xi2
	rymup=xmu*xi*rymu-ry1
	rjmu=w/(rymup-f*rymu)
	do i=1,nl
		rytemp=(xmu+i)*xi2*ry1-rymu
		rymu=ry1
		ry1=rytemp
	end do
	fact=rjmu/rjl
	rj=rjl1*fact
	rjp=rjp1*fact
	ry=rymu
	ryp=xnu*xi*rymu-ry1
	END SUBROUTINE bessjy_xltxmin
!BL
	SUBROUTINE bessjy_xgexmin(x,xnu,isign,h,rj,ry,rjp,ryp)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), INTENT(IN) :: xnu
	INTEGER(I4B),DIMENSION(:), INTENT(IN) :: isign
	REAL(DP), DIMENSION(:), INTENT(IN) :: h
	REAL(DP), DIMENSION(:), INTENT(OUT) :: rj,rjp,ry,ryp
	INTEGER(I4B) :: i,nlmax
	INTEGER(I4B),DIMENSION(size(x)) :: nl
	REAL(DP), DIMENSION(size(x)) :: &
		a,f,fact,gam,p,q,rjl,rjl1,rjmu,rjpl,rjp1,rjtemp,ry1,rymu,rymup,rytemp,w,xi,xi2,xmu,xmu2
	COMPLEX(DPC), DIMENSION(size(x)) :: aa,bb,cc,dd,dl,pq
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	nl=max(0,int(xnu-x+1.5_dp))
	nlmax=maxval(nl)
	xmu=xnu-nl
	xi=1.0_dp/x
	xi2=2.0_dp*xi
	w=xi2/PI_D
	xmu2=xmu*xmu
	rjl=isign*FPMIN
	rjpl=h*rjl
	rjl1=rjl
	rjp1=rjpl
	fact=xnu*xi
	do i=nlmax,1,-1
		converged=(i > nl)
		if (all(converged)) exit
		where (.not. converged)
			rjtemp=fact*rjl+rjpl
			fact=fact-xi
			rjpl=fact*rjtemp-rjl
			rjl=rjtemp
		end where
	end do
	where (rjl == 0.0) rjl=EPS
	f=rjpl/rjl
	a=0.25_dp-xmu2
	pq=cmplx(-0.5_dp*xi,1.0_dp,kind=dpc)
	aa=cmplx(0.0_dp,xi*a,kind=dpc)
	bb=cmplx(2.0_dp*x,2.0_dp,kind=dpc)
	cc=bb+aa/pq
	dd=1.0_dp/bb
	pq=cc*dd*pq
	converged=.false.
	do i=2,MAXIT
		where (.not. converged)
			a=a+2*(i-1)
			bb=bb+cmplx(0.0_dp,2.0_dp,kind=dpc)
			dd=a*dd+bb
			dd=merge(cmplx(FPMIN,kind=dpc),dd, absc(dd) < FPMIN )
			cc=bb+a/cc
			cc=merge(cmplx(FPMIN,kind=dpc),cc, absc(cc) < FPMIN )
			dd=1.0_dp/dd
			dl=cc*dd
			pq=pq*dl
			converged=(absc(dl-1.0_dp) < EPS)
		end where
		if (all(converged)) exit
	end do
	if (i > MAXIT) call nrerror('bessjy: cf2 section failed')
	p=real(pq,dp)
	q=aimag(pq)
	gam=(p-f)/q
	rjmu=sqrt(w/((p-f)*gam+q))
	rjmu=sign(rjmu,rjl)
	rymu=rjmu*gam
	rymup=rymu*(p+q/gam)
	ry1=xmu*xi*rymu-rymup
	do i=1,nlmax
		converged=(i > nl)
		if (all(converged)) exit
		where (.not. converged)
			rytemp=(xmu+i)*xi2*ry1-rymu
			rymu=ry1
			ry1=rytemp
		end where
	end do
	fact=rjmu/rjl
	rj=rjl1*fact
	rjp=rjp1*fact
	ry=rymu
	ryp=xnu*xi*rymu-ry1
	END SUBROUTINE bessjy_xgexmin
!BL
	FUNCTION absc(z)
	IMPLICIT NONE
	COMPLEX(DPC), DIMENSION(:), INTENT(IN) :: z
	REAL(DP), DIMENSION(size(z)) :: absc
	absc=abs(real(z))+abs(aimag(z))
	END FUNCTION absc
!BL
	FUNCTION unpacked(mask,vtrue,vfalse)
	IMPLICIT NONE
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
	REAL(DP), DIMENSION(:), INTENT(IN) :: vtrue, vfalse
	REAL(DP), DIMENSION(size(mask)) :: unpacked
	unpacked=unpack(vtrue,converged,0.0_dp)
	unpacked=unpack(vfalse,.not. converged,unpacked)
	END FUNCTION unpacked
	END SUBROUTINE bessjy_v

! -----------------------------------------------------------

	FUNCTION bessk_s(n,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessk0,bessk1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessk_s
	INTEGER(I4B) :: j
	REAL(SP) :: bk,bkm,bkp,tox
	call assert(n >= 2, x > 0.0, 'bessk_s args')
	tox=2.0_sp/x
	bkm=bessk0(x)
	bk=bessk1(x)
	do j=1,n-1
		bkp=bkm+j*tox*bk
		bkm=bk
		bk=bkp
	end do
	bessk_s=bk
	END FUNCTION bessk_s


	FUNCTION bessk_v(n,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessk0,bessk1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessk_v
	INTEGER(I4B) :: j
	REAL(SP), DIMENSION(size(x)) :: bk,bkm,bkp,tox
	call assert(n >= 2, all(x > 0.0), 'bessk_v args')
	tox=2.0_sp/x
	bkm=bessk0(x)
	bk=bessk1(x)
	do j=1,n-1
		bkp=bkm+j*tox*bk
		bkm=bk
		bk=bkp
	end do
	bessk_v=bk
	END FUNCTION bessk_v

! -----------------------------------------------------------

	FUNCTION bessk0_s(x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,poly
	use mo_nr, ONLY : bessi0
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessk0_s
	REAL(DP) :: y
	REAL(DP), DIMENSION(7) :: p = (/-0.57721566_dp,0.42278420_dp,&
		0.23069756_dp,0.3488590e-1_dp,0.262698e-2_dp,0.10750e-3_dp,&
		0.74e-5_dp/)
	REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,-0.7832358e-1_dp,&
		0.2189568e-1_dp,-0.1062446e-1_dp,0.587872e-2_dp,&
		-0.251540e-2_dp,0.53208e-3_dp/)
	call assert(x > 0.0, 'bessk0_s arg')
	if (x <= 2.0) then
		y=x*x/4.0_sp
		bessk0_s=(-log(x/2.0_sp)*bessi0(x))+poly(y,p)
	else
		y=(2.0_sp/x)
		bessk0_s=(exp(-x)/sqrt(x))*poly(y,q)
	end if
	END FUNCTION bessk0_s


	FUNCTION bessk0_v(x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,poly
	use mo_nr, ONLY : bessi0
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessk0_v
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(7) :: p = (/-0.57721566_dp,0.42278420_dp,&
		0.23069756_dp,0.3488590e-1_dp,0.262698e-2_dp,0.10750e-3_dp,&
		0.74e-5_dp/)
	REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,-0.7832358e-1_dp,&
		0.2189568e-1_dp,-0.1062446e-1_dp,0.587872e-2_dp,&
		-0.251540e-2_dp,0.53208e-3_dp/)
	call assert(all(x > 0.0), 'bessk0_v arg')
	mask = (x <= 2.0)
	where (mask)
		y=x*x/4.0_sp
		bessk0_v=(-log(x/2.0_sp)*bessi0(x))+poly(y,p,mask)
	elsewhere
		y=(2.0_sp/x)
		bessk0_v=(exp(-x)/sqrt(x))*poly(y,q,.not. mask)
	end where
	END FUNCTION bessk0_v

! -----------------------------------------------------------

	FUNCTION bessk1_s(x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,poly
	use mo_nr, ONLY : bessi1
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessk1_s
	REAL(DP) :: y
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,0.15443144_dp,&
		-0.67278579_dp,-0.18156897_dp,-0.1919402e-1_dp,&
		-0.110404e-2_dp,-0.4686e-4_dp/)
	REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,0.23498619_dp,&
		-0.3655620e-1_dp,0.1504268e-1_dp,-0.780353e-2_dp,&
		0.325614e-2_dp,-0.68245e-3_dp/)
	call assert(x > 0.0, 'bessk1_s arg')
	if (x <= 2.0) then
		y=x*x/4.0_sp
		bessk1_s=(log(x/2.0_sp)*bessi1(x))+(1.0_sp/x)*poly(y,p)
	else
		y=2.0_sp/x
		bessk1_s=(exp(-x)/sqrt(x))*poly(y,q)
	end if
	END FUNCTION bessk1_s


	FUNCTION bessk1_v(x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,poly
	use mo_nr, ONLY : bessi1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessk1_v
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,0.15443144_dp,&
		-0.67278579_dp,-0.18156897_dp,-0.1919402e-1_dp,&
		-0.110404e-2_dp,-0.4686e-4_dp/)
	REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,0.23498619_dp,&
		-0.3655620e-1_dp,0.1504268e-1_dp,-0.780353e-2_dp,&
		0.325614e-2_dp,-0.68245e-3_dp/)
	call assert(all(x > 0.0), 'bessk1_v arg')
	mask = (x <= 2.0)
	where (mask)
		y=x*x/4.0_sp
		bessk1_v=(log(x/2.0_sp)*bessi1(x))+(1.0_sp/x)*poly(y,p,mask)
	elsewhere
		y=2.0_sp/x
		bessk1_v=(exp(-x)/sqrt(x))*poly(y,q,.not. mask)
	end where
	END FUNCTION bessk1_v

! -----------------------------------------------------------

	FUNCTION bessy_s(n,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessy0,bessy1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessy_s
	INTEGER(I4B) :: j
	REAL(SP) :: by,bym,byp,tox
	call assert(n >= 2, x > 0.0, 'bessy_s args')
	tox=2.0_sp/x
	by=bessy1(x)
	bym=bessy0(x)
	do j=1,n-1
		byp=j*tox*by-bym
		bym=by
		by=byp
	end do
	bessy_s=by
	END FUNCTION bessy_s


	FUNCTION bessy_v(n,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessy0,bessy1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessy_v
	INTEGER(I4B) :: j
	REAL(SP), DIMENSION(size(x)) :: by,bym,byp,tox
	call assert(n >= 2, all(x > 0.0), 'bessy_v args')
	tox=2.0_sp/x
	by=bessy1(x)
	bym=bessy0(x)
	do j=1,n-1
		byp=j*tox*by-bym
		bym=by
		by=byp
	end do
	bessy_v=by
	END FUNCTION bessy_v

! -----------------------------------------------------------

	FUNCTION bessy0_s(x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,poly
	use mo_nr, ONLY : bessj0
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessy0_s
	REAL(SP) :: xx,z
	REAL(DP) :: y
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
		0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/-0.1562499995e-1_dp,&
		0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
		-0.934945152e-7_dp/)
	REAL(DP), DIMENSION(6) :: r = (/-2957821389.0_dp,7062834065.0_dp,&
		-512359803.6_dp,10879881.29_dp,-86327.92757_dp,&
		228.4622733_dp/)
	REAL(DP), DIMENSION(6) :: s = (/40076544269.0_dp,745249964.8_dp,&
		7189466.438_dp,47447.26470_dp,226.1030244_dp,1.0_dp/)
	call assert(x > 0.0, 'bessy0_s arg')
	if (abs(x) < 8.0) then
		y=x**2
		bessy0_s=(poly(y,r)/poly(y,s))+&
			0.636619772_sp*bessj0(x)*log(x)
	else
		z=8.0_sp/x
		y=z**2
		xx=x-0.785398164_sp
		bessy0_s=sqrt(0.636619772_sp/x)*(sin(xx)*&
			poly(y,p)+z*cos(xx)*poly(y,q))
	end if
	END FUNCTION bessy0_s


	FUNCTION bessy0_v(x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,poly
	use mo_nr, ONLY : bessj0
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessy0_v
	REAL(SP), DIMENSION(size(x)) :: xx,z
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
		0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/-0.1562499995e-1_dp,&
		0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
		-0.934945152e-7_dp/)
	REAL(DP), DIMENSION(6) :: r = (/-2957821389.0_dp,7062834065.0_dp,&
		-512359803.6_dp,10879881.29_dp,-86327.92757_dp,&
		228.4622733_dp/)
	REAL(DP), DIMENSION(6) :: s = (/40076544269.0_dp,745249964.8_dp,&
		7189466.438_dp,47447.26470_dp,226.1030244_dp,1.0_dp/)
	call assert(all(x > 0.0), 'bessy0_v arg')
	mask = (abs(x) < 8.0)
	where (mask)
		y=x**2
		bessy0_v=(poly(y,r,mask)/poly(y,s,mask))+&
			0.636619772_sp*bessj0(x)*log(x)
	elsewhere
		z=8.0_sp/x
		y=z**2
		xx=x-0.785398164_sp
		bessy0_v=sqrt(0.636619772_sp/x)*(sin(xx)*&
			poly(y,p,.not. mask)+z*cos(xx)*poly(y,q,.not. mask))
	end where
	END FUNCTION bessy0_v

! -----------------------------------------------------------

	FUNCTION bessy1_s(x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,poly
	use mo_nr, ONLY : bessj1
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessy1_s
	REAL(SP) :: xx,z
	REAL(DP) :: y
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,0.183105e-2_dp,&
		-0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/0.04687499995_dp,&
		-0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
		0.105787412e-6_dp/)
	REAL(DP), DIMENSION(6) :: r = (/-0.4900604943e13_dp,&
		0.1275274390e13_dp,-0.5153438139e11_dp,0.7349264551e9_dp,&
		-0.4237922726e7_dp,0.8511937935e4_dp/)
	REAL(DP), DIMENSION(7) :: s = (/0.2499580570e14_dp,&
		0.4244419664e12_dp,0.3733650367e10_dp,0.2245904002e8_dp,&
		0.1020426050e6_dp,0.3549632885e3_dp,1.0_dp/)
	call assert(x > 0.0, 'bessy1_s arg')
	if (abs(x) < 8.0) then
		y=x**2
		bessy1_s=x*(poly(y,r)/poly(y,s))+&
			0.636619772_sp*(bessj1(x)*log(x)-1.0_sp/x)
	else
		z=8.0_sp/x
		y=z**2
		xx=x-2.356194491_sp
		bessy1_s=sqrt(0.636619772_sp/x)*(sin(xx)*&
			poly(y,p)+z*cos(xx)*poly(y,q))
	end if
	END FUNCTION bessy1_s


	FUNCTION bessy1_v(x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,poly
	use mo_nr, ONLY : bessj1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: bessy1_v
	REAL(SP), DIMENSION(size(x)) :: xx,z
	REAL(DP), DIMENSION(size(x)) :: y
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,0.183105e-2_dp,&
		-0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
	REAL(DP), DIMENSION(5) :: q = (/0.04687499995_dp,&
		-0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
		0.105787412e-6_dp/)
	REAL(DP), DIMENSION(6) :: r = (/-0.4900604943e13_dp,&
		0.1275274390e13_dp,-0.5153438139e11_dp,0.7349264551e9_dp,&
		-0.4237922726e7_dp,0.8511937935e4_dp/)
	REAL(DP), DIMENSION(7) :: s = (/0.2499580570e14_dp,&
		0.4244419664e12_dp,0.3733650367e10_dp,0.2245904002e8_dp,&
		0.1020426050e6_dp,0.3549632885e3_dp,1.0_dp/)
	call assert(all(x > 0.0), 'bessy1_v arg')
	mask = (abs(x) < 8.0)
	where (mask)
		y=x**2
		bessy1_v=x*(poly(y,r,mask)/poly(y,s,mask))+&
			0.636619772_sp*(bessj1(x)*log(x)-1.0_sp/x)
	elsewhere
		z=8.0_sp/x
		y=z**2
		xx=x-2.356194491_sp
		bessy1_v=sqrt(0.636619772_sp/x)*(sin(xx)*&
			poly(y,p,.not. mask)+z*cos(xx)*poly(y,q,.not. mask))
	end where
	END FUNCTION bessy1_v

! -----------------------------------------------------------

	FUNCTION beta_s(z,w)
	use mo_nrtype
	use mo_nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: z,w
	REAL(SP) :: beta_s
	beta_s=exp(gammln(z)+gammln(w)-gammln(z+w))
	END FUNCTION beta_s


	FUNCTION beta_v(z,w)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: z,w
	REAL(SP), DIMENSION(size(z)) :: beta_v
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(z),size(w),'beta_v')
	beta_v=exp(gammln(z)+gammln(w)-gammln(z+w))
	END FUNCTION beta_v

! -----------------------------------------------------------

	FUNCTION betacf_s(a,b,x)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b,x
	REAL(SP) :: betacf_s
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(SP), PARAMETER :: EPS=epsilon(x), FPMIN=tiny(x)/EPS
	REAL(SP) :: aa,c,d,del,h,qab,qam,qap
	INTEGER(I4B) :: m,m2
	qab=a+b
	qap=a+1.0_sp
	qam=a-1.0_sp
	c=1.0
	d=1.0_sp-qab*x/qap
	if (abs(d) < FPMIN) d=FPMIN
	d=1.0_sp/d
	h=d
	do m=1,MAXIT
		m2=2*m
		aa=m*(b-m)*x/((qam+m2)*(a+m2))
		d=1.0_sp+aa*d
		if (abs(d) < FPMIN) d=FPMIN
		c=1.0_sp+aa/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_sp/d
		h=h*d*c
		aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
		d=1.0_sp+aa*d
		if (abs(d) < FPMIN) d=FPMIN
		c=1.0_sp+aa/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_sp/d
		del=d*c
		h=h*del
		if (abs(del-1.0_sp) <= EPS) exit
	end do
	if (m > MAXIT)&
		call nrerror('a or b too big, or MAXIT too small in betacf_s')
	betacf_s=h
	END FUNCTION betacf_s


	FUNCTION betacf_v(a,b,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,x
	REAL(SP), DIMENSION(size(x)) :: betacf_v
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(SP), PARAMETER :: EPS=epsilon(x), FPMIN=tiny(x)/EPS
	REAL(SP), DIMENSION(size(x)) :: aa,c,d,del,h,qab,qam,qap
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	INTEGER(I4B) :: m
	INTEGER(I4B), DIMENSION(size(x)) :: m2
	m=assert_eq(size(a),size(b),size(x),'betacf_v')
	qab=a+b
	qap=a+1.0_sp
	qam=a-1.0_sp
	c=1.0
	d=1.0_sp-qab*x/qap
	where (abs(d) < FPMIN) d=FPMIN
	d=1.0_sp/d
	h=d
	converged=.false.
	do m=1,MAXIT
		where (.not. converged)
			m2=2*m
			aa=m*(b-m)*x/((qam+m2)*(a+m2))
			d=1.0_sp+aa*d
			d=merge(FPMIN,d, abs(d)<FPMIN )
			c=1.0_sp+aa/c
			c=merge(FPMIN,c, abs(c)<FPMIN )
			d=1.0_sp/d
			h=h*d*c
			aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
			d=1.0_sp+aa*d
			d=merge(FPMIN,d, abs(d)<FPMIN )
			c=1.0_sp+aa/c
			c=merge(FPMIN,c, abs(c)<FPMIN )
			d=1.0_sp/d
			del=d*c
			h=h*del
			converged = (abs(del-1.0_sp) <= EPS)
		end where
		if (all(converged)) exit
	end do
	if (m > MAXIT)&
		call nrerror('a or b too big, or MAXIT too small in betacf_v')
	betacf_v=h
	END FUNCTION betacf_v

! -----------------------------------------------------------

	FUNCTION betai_s(a,b,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : betacf,gammln
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b,x
	REAL(SP) :: betai_s
	REAL(SP) :: bt
	call assert(x >= 0.0, x <= 1.0, 'betai_s arg')
	if (x == 0.0 .or. x == 1.0) then
		bt=0.0
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)&
			+a*log(x)+b*log(1.0_sp-x))
	end if
	if (x < (a+1.0_sp)/(a+b+2.0_sp)) then
		betai_s=bt*betacf(a,b,x)/a
	else
		betai_s=1.0_sp-bt*betacf(b,a,1.0_sp-x)/b
	end if
	END FUNCTION betai_s


	FUNCTION betai_v(a,b,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq
	use mo_nr, ONLY : betacf,gammln
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,x
	REAL(SP), DIMENSION(size(a)) :: betai_v
	REAL(SP), DIMENSION(size(a)) :: bt
	LOGICAL(LGT), DIMENSION(size(a)) :: mask
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(a),size(b),size(x),'betai_v')
	call assert(all(x >= 0.0), all(x <= 1.0), 'betai_v arg')
	where (x == 0.0 .or. x == 1.0)
		bt=0.0
	elsewhere
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)&
			+a*log(x)+b*log(1.0_sp-x))
	end where
	mask=(x < (a+1.0_sp)/(a+b+2.0_sp))
	betai_v=bt*betacf(merge(a,b,mask),merge(b,a,mask),&
		merge(x,1.0_sp-x,mask))/merge(a,b,mask)
	where (.not. mask) betai_v=1.0_sp-betai_v
	END FUNCTION betai_v

! -----------------------------------------------------------

	FUNCTION bico_s(n,k)
	use mo_nrtype
	use mo_nr, ONLY : factln
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n,k
	REAL(SP) :: bico_s
	bico_s=nint(exp(factln(n)-factln(k)-factln(n-k)))
	END FUNCTION bico_s


	FUNCTION bico_v(n,k)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : factln
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n,k
	REAL(SP), DIMENSION(size(n)) :: bico_v
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(n),size(k),'bico_v')
	bico_v=nint(exp(factln(n)-factln(k)-factln(n-k)))
	END FUNCTION bico_v

! -----------------------------------------------------------

	FUNCTION bnldev(pp,n)
	use mo_nrtype
	use mo_nr, ONLY : gammln,ran1
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: pp
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP) :: bnldev
	INTEGER(I4B) :: j
	INTEGER(I4B), SAVE :: nold=-1
	REAL(SP) :: am,em,g,h,p,sq,t,y,arr(24)
	REAL(SP), SAVE :: pc,plog,pclog,en,oldg,pold=-1.0
	p=merge(pp,1.0_sp-pp, pp <= 0.5_sp )
	am=n*p
	if (n < 25) then
		call ran1(arr(1:n))
		bnldev=count(arr(1:n)<p)
	else if (am < 1.0) then
		g=exp(-am)
		t=1.0
		do j=0,n
			call ran1(h)
			t=t*h
			if (t < g) exit
		end do
		bnldev=merge(j,n, j <= n)
	else
		if (n /= nold) then
			en=n
			oldg=gammln(en+1.0_sp)
			nold=n
		end if
		if (p /= pold) then
			pc=1.0_sp-p
			plog=log(p)
			pclog=log(pc)
			pold=p
		end if
		sq=sqrt(2.0_sp*am*pc)
		do
			call ran1(h)
			y=tan(PI*h)
			em=sq*y+am
			if (em < 0.0 .or. em >= en+1.0_sp) cycle
			em=int(em)
			t=1.2_sp*sq*(1.0_sp+y**2)*exp(oldg-gammln(em+1.0_sp)-&
				gammln(en-em+1.0_sp)+em*plog+(en-em)*pclog)
			call ran1(h)
			if (h <= t) exit
		end do
		bnldev=em
	end if
	if (p /= pp) bnldev=n-bnldev
	END FUNCTION bnldev

! -----------------------------------------------------------

	FUNCTION brent(ax,bx,cx,func,tol,xmin)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: ax,bx,cx,tol
	REAL(SP), INTENT(OUT) :: xmin
	REAL(SP) :: brent
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: CGOLD=0.3819660_sp,ZEPS=1.0e-3_sp*epsilon(ax)
	INTEGER(I4B) :: iter
	REAL(SP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
	a=min(ax,cx)
	b=max(ax,cx)
	v=bx
	w=v
	x=v
	e=0.0
	fx=func(x)
	fv=fx
	fw=fx
	do iter=1,ITMAX
		xm=0.5_sp*(a+b)
		tol1=tol*abs(x)+ZEPS
		tol2=2.0_sp*tol1
		if (abs(x-xm) <= (tol2-0.5_sp*(b-a))) then
			xmin=x
			brent=fx
			RETURN
		end if
		if (abs(e) > tol1) then
			r=(x-w)*(fx-fv)
			q=(x-v)*(fx-fw)
			p=(x-v)*q-(x-w)*r
			q=2.0_sp*(q-r)
			if (q > 0.0) p=-p
			q=abs(q)
			etemp=e
			e=d
			if (abs(p) >= abs(0.5_sp*q*etemp) .or. &
				p <= q*(a-x) .or. p >= q*(b-x)) then
				e=merge(a-x,b-x, x >= xm )
				d=CGOLD*e
			else
				d=p/q
				u=x+d
				if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
			end if
		else
			e=merge(a-x,b-x, x >= xm )
			d=CGOLD*e
		end if
		u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
		fu=func(u)
		if (fu <= fx) then
			if (u >= x) then
				a=x
			else
				b=x
			end if
			call shft(v,w,x,u)
			call shft(fv,fw,fx,fu)
		else
			if (u < x) then
				a=u
			else
				b=u
			end if
			if (fu <= fw .or. w == x) then
				v=w
				fv=fw
				w=u
				fw=fu
			else if (fu <= fv .or. v == x .or. v == w) then
				v=u
				fv=fu
			end if
		end if
	end do
	call nrerror('brent: exceed maximum iterations')
	CONTAINS
!BL
	SUBROUTINE shft(a,b,c,d)
	REAL(SP), INTENT(OUT) :: a
	REAL(SP), INTENT(INOUT) :: b,c
	REAL(SP), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft
	END FUNCTION brent

! -----------------------------------------------------------

	SUBROUTINE broydn(x,check)
	use mo_nrtype
          use mo_nrutil, ONLY : get_diag,lower_triangle,nrerror,&
		outerprod,put_diag,unit_matrix,vabs
	use mo_nr, ONLY : fdjac,lnsrch,qrdcmp,qrupdt,rsolv
	USE fminln
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
	LOGICAL(LGT), INTENT(OUT) :: check
	INTEGER(I4B), PARAMETER :: MAXITS=200
	REAL(SP), PARAMETER :: EPS=epsilon(x),TOLF=1.0e-4_sp,TOLMIN=1.0e-6_sp,&
		TOLX=EPS,STPMX=100.0
	INTEGER(I4B) :: i,its,k,n
	REAL(SP) :: f,fold,stpmax
	REAL(SP), DIMENSION(size(x)), TARGET :: fvec
	REAL(SP), DIMENSION(size(x)) :: c,d,fvcold,g,p,s,t,w,xold
	REAL(SP), DIMENSION(size(x),size(x)) :: qt,r
	LOGICAL :: restrt,sing
	fmin_fvecp=>fvec
	n=size(x)
	f=fmin(x)
	if (maxval(abs(fvec(:))) < 0.01_sp*TOLF) then
		check=.false.
		RETURN
	end if
	stpmax=STPMX*max(vabs(x(:)),real(n,sp))
	restrt=.true.
	do its=1,MAXITS
		if (restrt) then
			call fdjac(x,fvec,r)
			call qrdcmp(r,c,d,sing)
			if (sing) call nrerror('singular Jacobian in broydn')
			call unit_matrix(qt)
			do k=1,n-1
				if (c(k) /= 0.0) then
					qt(k:n,:)=qt(k:n,:)-outerprod(r(k:n,k),&
						matmul(r(k:n,k),qt(k:n,:)))/c(k)
				end if
			end do
			where (lower_triangle(n,n)) r(:,:)=0.0
			call put_diag(d(:),r(:,:))
		else
			s(:)=x(:)-xold(:)
			do i=1,n
				t(i)=dot_product(r(i,i:n),s(i:n))
			end do
			w(:)=fvec(:)-fvcold(:)-matmul(t(:),qt(:,:))
			where (abs(w(:)) < EPS*(abs(fvec(:))+abs(fvcold(:)))) &
				w(:)=0.0
			if (any(w(:) /= 0.0)) then
				t(:)=matmul(qt(:,:),w(:))
				s(:)=s(:)/dot_product(s,s)
				call qrupdt(r,qt,t,s)
				d(:)=get_diag(r(:,:))
				if (any(d(:) == 0.0)) &
					call nrerror('r singular in broydn')
			end if
		end if
		p(:)=-matmul(qt(:,:),fvec(:))
		do i=1,n
			g(i)=-dot_product(r(1:i,i),p(1:i))
		end do
		xold(:)=x(:)
		fvcold(:)=fvec(:)
		fold=f
		call rsolv(r,d,p)
		call lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin)
		if (maxval(abs(fvec(:))) < TOLF) then
			check=.false.
			RETURN
		end if
		if (check) then
			if (restrt .or. maxval(abs(g(:))*max(abs(x(:)), &
				1.0_sp)/max(f,0.5_sp*n)) < TOLMIN) RETURN
			restrt=.true.
		else
			restrt=.false.
			if (maxval((abs(x(:)-xold(:)))/max(abs(x(:)), &
				1.0_sp)) < TOLX) RETURN
		end if
	end do
	call nrerror('MAXITS exceeded in broydn')
	END SUBROUTINE broydn

! -----------------------------------------------------------

	SUBROUTINE bsstep(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq,cumsum,iminloc,nrerror,&
		outerdiff,outerprod,upper_triangle
	use mo_nr, ONLY : mmid,pzextr
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
	REAL(SP), INTENT(INOUT) :: x
	REAL(SP), INTENT(IN) :: htry,eps
	REAL(SP), INTENT(OUT) :: hdid,hnext
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B), PARAMETER :: IMAX=9, KMAXX=IMAX-1
	REAL(SP), PARAMETER :: SAFE1=0.25_sp,SAFE2=0.7_sp,REDMAX=1.0e-5_sp,&
		REDMIN=0.7_sp,TINY=1.0e-30_sp,SCALMX=0.1_sp
	INTEGER(I4B) :: k,km,ndum
	INTEGER(I4B), DIMENSION(IMAX) :: nseq = (/ 2,4,6,8,10,12,14,16,18 /)
	INTEGER(I4B), SAVE :: kopt,kmax
	REAL(SP), DIMENSION(KMAXX,KMAXX), SAVE :: alf
	REAL(SP), DIMENSION(KMAXX) :: err
	REAL(SP), DIMENSION(IMAX), SAVE :: a
	REAL(SP), SAVE :: epsold = -1.0_sp,xnew
	REAL(SP) :: eps1,errmax,fact,h,red,scale,wrkmin,xest
	REAL(SP), DIMENSION(size(y)) :: yerr,ysav,yseq
	LOGICAL(LGT) :: reduct
	LOGICAL(LGT), SAVE :: first=.true.
	ndum=assert_eq(size(y),size(dydx),size(yscal),'bsstep')
	if (eps /= epsold) then
		hnext=-1.0e29_sp
		xnew=-1.0e29_sp
		eps1=SAFE1*eps
		a(:)=cumsum(nseq,1)
		where (upper_triangle(KMAXX,KMAXX)) alf=eps1** &
			(outerdiff(a(2:),a(2:))/outerprod(arth( &
			3.0_sp,2.0_sp,KMAXX),(a(2:)-a(1)+1.0_sp)))
		epsold=eps
		do kopt=2,KMAXX-1
			if (a(kopt+1) > a(kopt)*alf(kopt-1,kopt)) exit
		end do
		kmax=kopt
	end if
	h=htry
	ysav(:)=y(:)
	if (h /= hnext .or. x /= xnew) then
		first=.true.
		kopt=kmax
	end if
	reduct=.false.
	main_loop: do
		do k=1,kmax
			xnew=x+h
			if (xnew == x) call nrerror('step size underflow in bsstep')
			call mmid(ysav,dydx,x,h,nseq(k),yseq,derivs)
			xest=(h/nseq(k))**2
			call pzextr(k,xest,yseq,y,yerr)
			if (k /= 1) then
				errmax=maxval(abs(yerr(:)/yscal(:)))
				errmax=max(TINY,errmax)/eps
				km=k-1
				err(km)=(errmax/SAFE1)**(1.0_sp/(2*km+1))
			end if
			if (k /= 1 .and. (k >= kopt-1 .or. first)) then
				if (errmax < 1.0) exit main_loop
				if (k == kmax .or. k == kopt+1) then
					red=SAFE2/err(km)
					exit
				else if (k == kopt) then
					if (alf(kopt-1,kopt) < err(km)) then
						red=1.0_sp/err(km)
						exit
					end if
				else if (kopt == kmax) then
					if (alf(km,kmax-1) < err(km)) then
						red=alf(km,kmax-1)*SAFE2/err(km)
						exit
					end if
				else if (alf(km,kopt) < err(km)) then
					red=alf(km,kopt-1)/err(km)
					exit
				end if
			end if
		end do
		red=max(min(red,REDMIN),REDMAX)
		h=h*red
		reduct=.true.
	end do main_loop
	x=xnew
	hdid=h
	first=.false.
	kopt=1+iminloc(a(2:km+1)*max(err(1:km),SCALMX))
	scale=max(err(kopt-1),SCALMX)
	wrkmin=scale*a(kopt)
	hnext=h/scale
	if (kopt >= k .and. kopt /= kmax .and. .not. reduct) then
		fact=max(scale/alf(kopt-1,kopt),SCALMX)
		if (a(kopt+1)*fact <= wrkmin) then
			hnext=h/fact
			kopt=kopt+1
		end if
	end if
	END SUBROUTINE bsstep

! -----------------------------------------------------------

	SUBROUTINE caldat(julian,mm,id,iyyy)
	use mo_nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: julian
	INTEGER(I4B), INTENT(OUT) :: mm,id,iyyy
	INTEGER(I4B) :: ja,jalpha,jb,jc,jd,je
	INTEGER(I4B), PARAMETER :: IGREG=2299161
	if (julian >= IGREG) then
		jalpha=int(((julian-1867216)-0.25_sp)/36524.25_sp)
		ja=julian+1+jalpha-int(0.25_sp*jalpha)
	else
		ja=julian
	end if
	jb=ja+1524
	jc=int(6680.0_sp+((jb-2439870)-122.1_sp)/365.25_sp)
	jd=365*jc+int(0.25_sp*jc)
	je=int((jb-jd)/30.6001_sp)
	id=jb-jd-int(30.6001_sp*je)
	mm=je-1
	if (mm > 12) mm=mm-12
	iyyy=jc-4715
	if (mm > 2) iyyy=iyyy-1
	if (iyyy <= 0) iyyy=iyyy-1
	END SUBROUTINE caldat

! -----------------------------------------------------------

	FUNCTION chder(a,b,c)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,cumsum
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(:), INTENT(IN) :: c
	REAL(SP), DIMENSION(size(c)) :: chder
	INTEGER(I4B) :: n
	REAL(SP) :: con
	REAL(SP), DIMENSION(size(c)) :: temp
	n=size(c)
	temp(1)=0.0
	temp(2:n)=2.0_sp*arth(n-1,-1,n-1)*c(n:2:-1)
	chder(n:1:-2)=cumsum(temp(1:n:2))
	chder(n-1:1:-2)=cumsum(temp(2:n:2))
	con=2.0_sp/(b-a)
	chder=chder*con
	END FUNCTION chder

! -----------------------------------------------------------

	FUNCTION chebev_s(a,b,c,x)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b,x
	REAL(SP), DIMENSION(:), INTENT(IN) :: c
	REAL(SP) :: chebev_s
	INTEGER(I4B) :: j,m
	REAL(SP) :: d,dd,sv,y,y2
	if ((x-a)*(x-b) > 0.0) call nrerror('x not in range in chebev_s')
	m=size(c)
	d=0.0
	dd=0.0
	y=(2.0_sp*x-a-b)/(b-a)
	y2=2.0_sp*y
	do j=m,2,-1
		sv=d
		d=y2*d-dd+c(j)
		dd=sv
	end do
	chebev_s=y*d-dd+0.5_sp*c(1)
	END FUNCTION chebev_s


	FUNCTION chebev_v(a,b,c,x)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(:), INTENT(IN) :: c,x
	REAL(SP), DIMENSION(size(x)) :: chebev_v
	INTEGER(I4B) :: j,m
	REAL(SP), DIMENSION(size(x)) :: d,dd,sv,y,y2
	if (any((x-a)*(x-b) > 0.0)) call nrerror('x not in range in chebev_v')
	m=size(c)
	d=0.0
	dd=0.0
	y=(2.0_sp*x-a-b)/(b-a)
	y2=2.0_sp*y
	do j=m,2,-1
		sv=d
		d=y2*d-dd+c(j)
		dd=sv
	end do
	chebev_v=y*d-dd+0.5_sp*c(1)
	END FUNCTION chebev_v

! -----------------------------------------------------------

	FUNCTION chebft(a,b,n,func)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,outerprod
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(n) :: chebft
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(DP) :: bma,bpa
	REAL(DP), DIMENSION(n) :: theta
	bma=0.5_dp*(b-a)
	bpa=0.5_dp*(b+a)
	theta(:)=PI_D*arth(0.5_dp,1.0_dp,n)/n
	chebft(:)=matmul(cos(outerprod(arth(0.0_dp,1.0_dp,n),theta)), &
		func(real(cos(theta)*bma+bpa,sp)))*2.0_dp/n
	END FUNCTION chebft

! -----------------------------------------------------------

	FUNCTION chebpc(c)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: c
	REAL(SP), DIMENSION(size(c)) :: chebpc
	INTEGER(I4B) :: j,n
	REAL(SP), DIMENSION(size(c)) :: dd,sv
	n=size(c)
	chebpc=0.0
	dd=0.0
	chebpc(1)=c(n)
	do j=n-1,2,-1
		sv(2:n-j+1)=chebpc(2:n-j+1)
		chebpc(2:n-j+1)=2.0_sp*chebpc(1:n-j)-dd(2:n-j+1)
		dd(2:n-j+1)=sv(2:n-j+1)
		sv(1)=chebpc(1)
		chebpc(1)=-dd(1)+c(j)
		dd(1)=sv(1)
	end do
	chebpc(2:n)=chebpc(1:n-1)-dd(2:n)
	chebpc(1)=-dd(1)+0.5_sp*c(1)
	END FUNCTION chebpc

! -----------------------------------------------------------

	FUNCTION chint(a,b,c)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(:), INTENT(IN) :: c
	REAL(SP), DIMENSION(size(c)) :: chint
	INTEGER(I4B) :: n
	REAL(SP) :: con
	n=size(c)
	con=0.25_sp*(b-a)
	chint(2:n-1)=con*(c(1:n-2)-c(3:n))/arth(1,1,n-2)
	chint(n)=con*c(n-1)/(n-1)
	chint(1)=2.0_sp*(sum(chint(2:n:2))-sum(chint(3:n:2)))
	END FUNCTION chint

! -----------------------------------------------------------

MODULE chixyfit
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	REAL(SP), DIMENSION(:), POINTER :: xxp,yyp,sxp,syp,wwp
	REAL(SP) :: aa,offs
CONTAINS
!BL
	FUNCTION chixy(bang)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: bang
	REAL(SP) :: chixy
	REAL(SP), PARAMETER :: BIG=1.0e30_sp
	REAL(SP) :: avex,avey,sumw,b
	if (.not. associated(wwp)) call nrerror("chixy: bad pointers")
	b=tan(bang)
	wwp(:)=(b*sxp(:))**2+syp(:)**2
	where (wwp(:) < 1.0/BIG)
		wwp(:)=BIG
	elsewhere
		wwp(:)=1.0_sp/wwp(:)
	end where
	sumw=sum(wwp)
	avex=dot_product(wwp,xxp)/sumw
	avey=dot_product(wwp,yyp)/sumw
	aa=avey-b*avex
	chixy=sum(wwp(:)*(yyp(:)-aa-b*xxp(:))**2)-offs
	END FUNCTION chixy
END MODULE chixyfit

! -----------------------------------------------------------

	SUBROUTINE choldc(a,p)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: p
	INTEGER(I4B) :: i,n
	REAL(SP) :: summ
	n=assert_eq(size(a,1),size(a,2),size(p),'choldc')
	do i=1,n
		summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
		if (summ <= 0.0) call nrerror('choldc failed')
		p(i)=sqrt(summ)
		a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
	end do
	END SUBROUTINE choldc

! -----------------------------------------------------------

	SUBROUTINE cholsl(a,p,b,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	REAL(SP), DIMENSION(:), INTENT(IN) :: p,b
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
	INTEGER(I4B) :: i,n
	n=assert_eq((/size(a,1),size(a,2),size(p),size(b),size(x)/),'cholsl')
	do i=1,n
		x(i)=(b(i)-dot_product(a(i,1:i-1),x(1:i-1)))/p(i)
	end do
	do i=n,1,-1
		x(i)=(x(i)-dot_product(a(i+1:n,i),x(i+1:n)))/p(i)
	end do
	END SUBROUTINE cholsl

! -----------------------------------------------------------

	SUBROUTINE chsone(bins,ebins,knstrn,df,chsq,prob)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : gammq
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: knstrn
	REAL(SP), INTENT(OUT) :: df,chsq,prob
	REAL(SP), DIMENSION(:), INTENT(IN) :: bins,ebins
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(bins),size(ebins),'chsone')
	if (any(ebins(:) <= 0.0)) call nrerror('bad expected number in chsone')
	df=size(bins)-knstrn
	chsq=sum((bins(:)-ebins(:))**2/ebins(:))
	prob=gammq(0.5_sp*df,0.5_sp*chsq)
	END SUBROUTINE chsone

! -----------------------------------------------------------

	SUBROUTINE chstwo(bins1,bins2,knstrn,df,chsq,prob)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : gammq
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: knstrn
	REAL(SP), INTENT(OUT) :: df,chsq,prob
	REAL(SP), DIMENSION(:), INTENT(IN) :: bins1,bins2
	INTEGER(I4B) :: ndum
	LOGICAL(LGT), DIMENSION(size(bins1)) :: nzeromask
	ndum=assert_eq(size(bins1),size(bins2),'chstwo')
	nzeromask = bins1(:) /= 0.0 .or. bins2(:) /= 0.0
	chsq=sum((bins1(:)-bins2(:))**2/(bins1(:)+bins2(:)),mask=nzeromask)
	df=count(nzeromask)-knstrn
	prob=gammq(0.5_sp*df,0.5_sp*chsq)
	END SUBROUTINE chstwo

! -----------------------------------------------------------

	SUBROUTINE cisi(x,ci,si)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), INTENT(OUT) :: ci,si
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=4.0_sp*tiny(x),&
		BIG=huge(x)*EPS,TMIN=2.0
	INTEGER(I4B) :: i,k
	REAL(SP) :: a,err,fact,sign,sum,sumc,sums,t,term
	COMPLEX(SPC) :: h,b,c,d,del
	LOGICAL(LGT) :: odd
	t=abs(x)
	if (t == 0.0) then
		si=0.0
		ci=-BIG
		RETURN
	end if
	if (t > TMIN) then
		b=cmplx(1.0_sp,t,kind=spc)
		c=BIG
		d=1.0_sp/b
		h=d
		do i=2,MAXIT
			a=-(i-1)**2
			b=b+2.0_sp
			d=1.0_sp/(a*d+b)
			c=b+a/c
			del=c*d
			h=h*del
			if (absc(del-1.0_sp) <= EPS) exit
		end do
		if (i > MAXIT) call nrerror('continued fraction failed in cisi')
		h=cmplx(cos(t),-sin(t),kind=spc)*h
		ci=-real(h)
		si=PIO2+aimag(h)
	else
		if (t < sqrt(FPMIN)) then
			sumc=0.0
			sums=t
		else
			sum=0.0
			sums=0.0
			sumc=0.0
			sign=1.0
			fact=1.0
			odd=.true.
			do k=1,MAXIT
				fact=fact*t/k
				term=fact/k
				sum=sum+sign*term
				err=term/abs(sum)
				if (odd) then
					sign=-sign
					sums=sum
					sum=sumc
				else
					sumc=sum
					sum=sums
				end if
				if (err < EPS) exit
				odd=.not. odd
			end do
			if (k > MAXIT) call nrerror('MAXIT exceeded in cisi')
		end if
		si=sums
		ci=sumc+log(t)+EULER
	end if
	if (x < 0.0) si=-si
	CONTAINS
!BL
	FUNCTION absc(z)
	IMPLICIT NONE
	COMPLEX(SPC), INTENT(IN) :: z
	REAL(SP) :: absc
	absc=abs(real(z))+abs(aimag(z))
	END FUNCTION absc
	END SUBROUTINE cisi

! -----------------------------------------------------------

	SUBROUTINE cntab1(nn,chisq,df,prob,cramrv,ccc)
	use mo_nrtype
          use mo_nrutil, ONLY : outerprod
	use mo_nr, ONLY : gammq
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: nn
	REAL(SP), INTENT(OUT) :: chisq,df,prob,cramrv,ccc
	REAL(SP), PARAMETER :: TINY=1.0e-30_sp
	INTEGER(I4B) :: nni,nnj
	REAL(SP) :: sumn
	REAL(SP), DIMENSION(size(nn,1)) :: sumi
	REAL(SP), DIMENSION(size(nn,2)) :: sumj
	REAL(SP), DIMENSION(size(nn,1),size(nn,2)) :: expctd
	sumi(:)=sum(nn(:,:),dim=2)
	sumj(:)=sum(nn(:,:),dim=1)
	sumn=sum(sumi(:))
	nni=size(sumi)-count(sumi(:) == 0.0)
	nnj=size(sumj)-count(sumj(:) == 0.0)
	df=nni*nnj-nni-nnj+1
	expctd(:,:)=outerprod(sumi(:),sumj(:))/sumn
	chisq=sum((nn(:,:)-expctd(:,:))**2/(expctd(:,:)+TINY))
	prob=gammq(0.5_sp*df,0.5_sp*chisq)
	cramrv=sqrt(chisq/(sumn*min(nni-1,nnj-1)))
	ccc=sqrt(chisq/(chisq+sumn))
	END SUBROUTINE cntab1

! -----------------------------------------------------------

	SUBROUTINE cntab2(nn,h,hx,hy,hygx,hxgy,uygx,uxgy,uxy)
	use mo_nrtype
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: nn
	REAL(SP), INTENT(OUT) :: h,hx,hy,hygx,hxgy,uygx,uxgy,uxy
	REAL(SP), PARAMETER :: TINY=1.0e-30_sp
	REAL(SP) :: sumn
	REAL(SP), DIMENSION(size(nn,1)) :: sumi
	REAL(SP), DIMENSION(size(nn,2)) :: sumj
	sumi(:)=sum(nn(:,:),dim=2)
	sumj(:)=sum(nn(:,:),dim=1)
	sumn=sum(sumi(:))
	hx=-sum(sumi(:)*log(sumi(:)/sumn), mask=(sumi(:) /= 0.0) )/sumn
	hy=-sum(sumj(:)*log(sumj(:)/sumn), mask=(sumj(:) /= 0.0) )/sumn
	h=-sum(nn(:,:)*log(nn(:,:)/sumn), mask=(nn(:,:) /= 0) )/sumn
	hygx=h-hx
	hxgy=h-hy
	uygx=(hy-hygx)/(hy+TINY)
	uxgy=(hx-hxgy)/(hx+TINY)
	uxy=2.0_sp*(hx+hy-h)/(hx+hy+TINY)
	END SUBROUTINE cntab2

! -----------------------------------------------------------

	FUNCTION convlv(data,respns,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,nrerror
	use mo_nr, ONLY : realft
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
	REAL(SP), DIMENSION(:), INTENT(IN) :: respns
	INTEGER(I4B), INTENT(IN) :: isign
	REAL(SP), DIMENSION(size(data)) :: convlv
	INTEGER(I4B) :: no2,n,m
	COMPLEX(SPC), DIMENSION(size(data)/2) :: tmpd,tmpr
	n=size(data)
	m=size(respns)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in convlv')
	call assert(mod(m,2)==1, 'm must be odd in convlv')
	convlv(1:m)=respns(:)
	convlv(n-(m-3)/2:n)=convlv((m+3)/2:m)
	convlv((m+3)/2:n-(m-1)/2)=0.0
	no2=n/2
	call realft(data,1,tmpd)
	call realft(convlv,1,tmpr)
	if (isign == 1) then
		tmpr(1)=cmplx(real(tmpd(1))*real(tmpr(1))/no2, &
			aimag(tmpd(1))*aimag(tmpr(1))/no2, kind=spc)
		tmpr(2:)=tmpd(2:)*tmpr(2:)/no2
	else if (isign == -1) then
		if (any(abs(tmpr(2:)) == 0.0) .or. real(tmpr(1)) == 0.0 &
			.or. aimag(tmpr(1)) == 0.0) call nrerror &
			('deconvolving at response zero in convlv')
		tmpr(1)=cmplx(real(tmpd(1))/real(tmpr(1))/no2, &
			aimag(tmpd(1))/aimag(tmpr(1))/no2, kind=spc)
		tmpr(2:)=tmpd(2:)/tmpr(2:)/no2
	else
		call nrerror('no meaning for isign in convlv')
	end if
	call realft(convlv,-1,tmpr)
	END FUNCTION convlv

! -----------------------------------------------------------

	FUNCTION correl(data1,data2)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq
	use mo_nr, ONLY : realft
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: data1,data2
	REAL(SP), DIMENSION(size(data1)) :: correl
	COMPLEX(SPC), DIMENSION(size(data1)/2) :: cdat1,cdat2
	INTEGER(I4B) :: no2,n
	n=assert_eq(size(data1),size(data2),'correl')
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in correl')
	no2=n/2
	call realft(data1,1,cdat1)
	call realft(data2,1,cdat2)
	cdat1(1)=cmplx(real(cdat1(1))*real(cdat2(1))/no2, &
		aimag(cdat1(1))*aimag(cdat2(1))/no2, kind=spc)
	cdat1(2:)=cdat1(2:)*conjg(cdat2(2:))/no2
	call realft(correl,-1,cdat1)
	END FUNCTION correl

! -----------------------------------------------------------

	SUBROUTINE cosft1(y)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,cumsum,zroots_unity
	use mo_nr, ONLY : realft
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	COMPLEX(SPC), DIMENSION((size(y)-1)/2) :: w
	REAL(SP), DIMENSION((size(y)-1)/2-1) :: y1,y2
	REAL(SP) :: summ
	INTEGER(I4B) :: n,nh
	n=size(y)-1
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in cosft1')
	nh=n/2
	w=zroots_unity(n+n,nh)
	summ=0.5_sp*(y(1)-y(n+1))
	y(1)=0.5_sp*(y(1)+y(n+1))
	y1=0.5_sp*(y(2:nh)+y(n:nh+2:-1))
	y2=y(2:nh)-y(n:nh+2:-1)
	summ=summ+sum(real(w(2:nh))*y2)
	y2=y2*aimag(w(2:nh))
	y(2:nh)=y1-y2
	y(n:nh+2:-1)=y1+y2
	call realft(y(1:n),1)
	y(n+1)=y(2)
	y(2)=summ
	y(2:n:2)=cumsum(y(2:n:2))
	END SUBROUTINE cosft1

! -----------------------------------------------------------

	SUBROUTINE cosft2(y,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,cumsum,zroots_unity
	use mo_nr, ONLY : realft
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(size(y)) :: w
	REAL(SP), DIMENSION(size(y)/2) :: y1,y2
	REAL(SP) :: ytemp
	INTEGER(I4B) :: n,nh
	n=size(y)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in cosft2')
	nh=n/2
	w=zroots_unity(4*n,n)
	if (isign == 1) then
		y1=0.5_sp*(y(1:nh)+y(n:nh+1:-1))
		y2=aimag(w(2:n:2))*(y(1:nh)-y(n:nh+1:-1))
		y(1:nh)=y1+y2
		y(n:nh+1:-1)=y1-y2
		call realft(y,1)
		y1(1:nh-1)=y(3:n-1:2)*real(w(3:n-1:2)) &
			-y(4:n:2)*aimag(w(3:n-1:2))
		y2(1:nh-1)=y(4:n:2)*real(w(3:n-1:2)) &
			+y(3:n-1:2)*aimag(w(3:n-1:2))
		y(3:n-1:2)=y1(1:nh-1)
		y(4:n:2)=y2(1:nh-1)
		ytemp=0.5_sp*y(2)
		y(n-2:2:-2)=cumsum(y(n:4:-2),ytemp)
		y(n)=ytemp
	else if (isign == -1) then
		ytemp=y(n)
		y(4:n:2)=y(2:n-2:2)-y(4:n:2)
		y(2)=2.0_sp*ytemp
		y1(1:nh-1)=y(3:n-1:2)*real(w(3:n-1:2)) &
			+y(4:n:2)*aimag(w(3:n-1:2))
		y2(1:nh-1)=y(4:n:2)*real(w(3:n-1:2)) &
			-y(3:n-1:2)*aimag(w(3:n-1:2))
		y(3:n-1:2)=y1(1:nh-1)
		y(4:n:2)=y2(1:nh-1)
		call realft(y,-1)
		y1=y(1:nh)+y(n:nh+1:-1)
		y2=(0.5_sp/aimag(w(2:n:2)))*(y(1:nh)-y(n:nh+1:-1))
		y(1:nh)=0.5_sp*(y1+y2)
		y(n:nh+1:-1)=0.5_sp*(y1-y2)
	end if
	END SUBROUTINE cosft2

! -----------------------------------------------------------

	SUBROUTINE covsrt(covar,maska)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
	INTEGER(I4B) :: ma,mfit,j,k
	ma=assert_eq(size(covar,1),size(covar,2),size(maska),'covsrt')
	mfit=count(maska)
	covar(mfit+1:ma,1:ma)=0.0
	covar(1:ma,mfit+1:ma)=0.0
	k=mfit
	do j=ma,1,-1
		if (maska(j)) then
			call swap(covar(1:ma,k),covar(1:ma,j))
			call swap(covar(k,1:ma),covar(j,1:ma))
			k=k-1
		end if
	end do
	END SUBROUTINE covsrt

! -----------------------------------------------------------

	SUBROUTINE cyclic(a,b,c,alpha,beta,r,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq
	use mo_nr, ONLY : tridag
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN):: a,b,c,r
	REAL(SP), INTENT(IN) :: alpha,beta
	REAL(SP), DIMENSION(:), INTENT(OUT):: x
	INTEGER(I4B) :: n
	REAL(SP) :: fact,gamma
	REAL(SP), DIMENSION(size(x)) :: bb,u,z
	n=assert_eq((/size(a),size(b),size(c),size(r),size(x)/),'cyclic')
	call assert(n > 2, 'cyclic arg')
	gamma=-b(1)
	bb(1)=b(1)-gamma
	bb(n)=b(n)-alpha*beta/gamma
	bb(2:n-1)=b(2:n-1)
	call tridag(a(2:n),bb,c(1:n-1),r,x)
	u(1)=gamma
	u(n)=alpha
	u(2:n-1)=0.0
	call tridag(a(2:n),bb,c(1:n-1),u,z)
	fact=(x(1)+beta*x(n)/gamma)/(1.0_sp+z(1)+beta*z(n)/gamma)
	x=x-fact*z
	END SUBROUTINE cyclic

! -----------------------------------------------------------

	SUBROUTINE daub4(a,isign)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
	INTEGER(I4B), INTENT(IN) :: isign
	REAL(SP), DIMENSION(size(a)) :: wksp
	REAL(SP), PARAMETER :: C0=0.4829629131445341_sp,&
		C1=0.8365163037378079_sp,C2=0.2241438680420134_sp,&
		C3=-0.1294095225512604_sp
	INTEGER(I4B) :: n,nh,nhp,nhm
	n=size(a)
	if (n < 4) RETURN
	nh=n/2
	nhp=nh+1
	nhm=nh-1
	if (isign >= 0) then
		wksp(1:nhm) = C0*a(1:n-3:2)+C1*a(2:n-2:2) &
			+C2*a(3:n-1:2)+C3*a(4:n:2)
		wksp(nh)=C0*a(n-1)+C1*a(n)+C2*a(1)+C3*a(2)
		wksp(nhp:n-1) = C3*a(1:n-3:2)-C2*a(2:n-2:2) &
			+C1*a(3:n-1:2)-C0*a(4:n:2)
		wksp(n)=C3*a(n-1)-C2*a(n)+C1*a(1)-C0*a(2)
	else
		wksp(1)=C2*a(nh)+C1*a(n)+C0*a(1)+C3*a(nhp)
		wksp(2)=C3*a(nh)-C0*a(n)+C1*a(1)-C2*a(nhp)
		wksp(3:n-1:2) = C2*a(1:nhm)+C1*a(nhp:n-1) &
			+C0*a(2:nh)+C3*a(nh+2:n)
		wksp(4:n:2) = C3*a(1:nhm)-C0*a(nhp:n-1) &
			+C1*a(2:nh)-C2*a(nh+2:n)
	end if
	a(1:n)=wksp(1:n)
	END SUBROUTINE daub4

! -----------------------------------------------------------

	FUNCTION dawson_s(x)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,geop
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: dawson_s
	INTEGER(I4B), PARAMETER :: NMAX=6
	REAL(SP), PARAMETER :: H=0.4_sp,A1=2.0_sp/3.0_sp,A2=0.4_sp,&
		A3=2.0_sp/7.0_sp
	INTEGER(I4B) :: i,n0
	REAL(SP) :: ec,x2,xp,xx
	REAL(SP), DIMENSION(NMAX) :: d1,d2,e1
	REAL(SP), DIMENSION(NMAX), SAVE :: c=(/ (0.0_sp,i=1,NMAX) /)
	if (c(1) == 0.0) c(1:NMAX)=exp(-(arth(1,2,NMAX)*H)**2)
	if (abs(x) < 0.2_sp) then
		x2=x**2
		dawson_s=x*(1.0_sp-A1*x2*(1.0_sp-A2*x2*(1.0_sp-A3*x2)))
	else
		xx=abs(x)
		n0=2*nint(0.5_sp*xx/H)
		xp=xx-real(n0,sp)*H
		ec=exp(2.0_sp*xp*H)
		d1=arth(n0+1,2,NMAX)
		d2=arth(n0-1,-2,NMAX)
		e1=geop(ec,ec**2,NMAX)
		dawson_s=0.5641895835477563_sp*sign(exp(-xp**2),x)*&
			sum(c*(e1/d1+1.0_sp/(d2*e1)))
	end if
	END FUNCTION dawson_s

	FUNCTION dawson_v(x)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: dawson_v
	INTEGER(I4B), PARAMETER :: NMAX=6
	REAL(SP), PARAMETER :: H=0.4_sp,A1=2.0_sp/3.0_sp,A2=0.4_sp,&
		A3=2.0_sp/7.0_sp
	INTEGER(I4B) :: i,n
	REAL(SP), DIMENSION(size(x)) :: x2
	REAL(SP), DIMENSION(NMAX), SAVE :: c=(/ (0.0_sp,i=1,NMAX) /)
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	if (c(1) == 0.0) c(1:NMAX)=exp(-(arth(1,2,NMAX)*H)**2)
	mask = (abs(x) >= 0.2_sp)
	dawson_v=dawsonseries_v(x,mask)
	where (.not. mask)
		x2=x**2
		dawson_v=x*(1.0_sp-A1*x2*(1.0_sp-A2*x2*(1.0_sp-A3*x2)))
	end where
	CONTAINS
!BL
	FUNCTION dawsonseries_v(xin,mask)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xin
	LOGICAL(LGT), DIMENSION(size(xin)), INTENT(IN) :: mask
	REAL(SP), DIMENSION(size(xin)) :: dawsonseries_v
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: n0
	REAL(SP), DIMENSION(:), ALLOCATABLE :: d1,d2,e1,e2,sm,xp,xx,x
	n=count(mask)
	if (n == 0) RETURN
	allocate(n0(n),d1(n),d2(n),e1(n),e2(n),sm(n),xp(n),xx(n),x(n))
	x=pack(xin,mask)
	xx=abs(x)
	n0=2*nint(0.5_sp*xx/H)
	xp=xx-real(n0,sp)*H
	e1=exp(2.0_sp*xp*H)
	e2=e1**2
	d1=n0+1.0_sp
	d2=d1-2.0_sp
	sm=0.0
	do i=1,NMAX
		sm=sm+c(i)*(e1/d1+1.0_sp/(d2*e1))
		d1=d1+2.0_sp
		d2=d2-2.0_sp
		e1=e2*e1
	end do
	sm=0.5641895835477563_sp*sign(exp(-xp**2),x)*sm
	dawsonseries_v=unpack(sm,mask,0.0_sp)
	deallocate(n0,d1,d2,e1,e2,sm,xp,xx)
	END FUNCTION dawsonseries_v
	END FUNCTION dawson_v

! -----------------------------------------------------------

	FUNCTION dbrent(ax,bx,cx,func,dfunc,tol,xmin)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: ax,bx,cx,tol
	REAL(SP), INTENT(OUT) :: xmin
	REAL(SP) :: dbrent
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
!BL
		FUNCTION dfunc(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: ZEPS=1.0e-3_sp*epsilon(ax)
	INTEGER(I4B) :: iter
	REAL(SP) :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,&
		u,u1,u2,v,w,x,xm
	LOGICAL :: ok1,ok2
	a=min(ax,cx)
	b=max(ax,cx)
	v=bx
	w=v
	x=v
	e=0.0
	fx=func(x)
	fv=fx
	fw=fx
	dx=dfunc(x)
	dv=dx
	dw=dx
	do iter=1,ITMAX
		xm=0.5_sp*(a+b)
		tol1=tol*abs(x)+ZEPS
		tol2=2.0_sp*tol1
		if (abs(x-xm) <= (tol2-0.5_sp*(b-a))) exit
		if (abs(e) > tol1) then
			d1=2.0_sp*(b-a)
			d2=d1
			if (dw /= dx) d1=(w-x)*dx/(dx-dw)
			if (dv /= dx) d2=(v-x)*dx/(dx-dv)
			u1=x+d1
			u2=x+d2
			ok1=((a-u1)*(u1-b) > 0.0) .and. (dx*d1 <= 0.0)
			ok2=((a-u2)*(u2-b) > 0.0) .and. (dx*d2 <= 0.0)
			olde=e
			e=d
			if (ok1 .or. ok2) then
				if (ok1 .and. ok2) then
					d=merge(d1,d2, abs(d1) < abs(d2))
				else
					d=merge(d1,d2,ok1)
				end if
				if (abs(d) <= abs(0.5_sp*olde)) then
					u=x+d
					if (u-a < tol2 .or. b-u < tol2) &
						d=sign(tol1,xm-x)
				else
					e=merge(a,b, dx >= 0.0)-x
					d=0.5_sp*e
				end if
			else
				e=merge(a,b, dx >= 0.0)-x
				d=0.5_sp*e
			end if
		else
			e=merge(a,b, dx >= 0.0)-x
			d=0.5_sp*e
		end if
		if (abs(d) >= tol1) then
			u=x+d
			fu=func(u)
		else
			u=x+sign(tol1,d)
			fu=func(u)
			if (fu > fx) exit
		end if
		du=dfunc(u)
		if (fu <= fx) then
			if (u >= x) then
				a=x
			else
				b=x
			end if
			call mov3(v,fv,dv,w,fw,dw)
			call mov3(w,fw,dw,x,fx,dx)
			call mov3(x,fx,dx,u,fu,du)
		else
			if (u < x) then
				a=u
			else
				b=u
			end if
			if (fu <= fw .or. w == x) then
				call mov3(v,fv,dv,w,fw,dw)
				call mov3(w,fw,dw,u,fu,du)
			else if (fu <= fv .or. v == x .or. v == w) then
				call mov3(v,fv,dv,u,fu,du)
			end if
		end if
	end do
	if (iter > ITMAX) call nrerror('dbrent: exceeded maximum iterations')
	xmin=x
	dbrent=fx
	CONTAINS
!BL
	SUBROUTINE mov3(a,b,c,d,e,f)
	REAL(SP), INTENT(IN) :: d,e,f
	REAL(SP), INTENT(OUT) :: a,b,c
	a=d
	b=e
	c=f
	END SUBROUTINE mov3
	END FUNCTION dbrent

! -----------------------------------------------------------

	SUBROUTINE ddpoly(c,x,pd)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,cumprod,poly_term
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: c
	REAL(SP), DIMENSION(:), INTENT(OUT) :: pd
	INTEGER(I4B) :: i,nc,nd
	REAL(SP), DIMENSION(size(pd)) :: fac
	REAL(SP), DIMENSION(size(c)) :: d
	nc=size(c)
	nd=size(pd)
	d(nc:1:-1)=poly_term(c(nc:1:-1),x)
	do i=2,min(nd,nc)
		d(nc:i:-1)=poly_term(d(nc:i:-1),x)
	end do
	pd=d(1:nd)
	fac=cumprod(arth(1.0_sp,1.0_sp,nd))
	pd(3:nd)=fac(2:nd-1)*pd(3:nd)
	END SUBROUTINE ddpoly

! -----------------------------------------------------------

	FUNCTION decchk(string,ch)
	use mo_nrtype
          use mo_nrutil, ONLY : ifirstloc
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: string
	CHARACTER(1), INTENT(OUT) :: ch
	LOGICAL(LGT) :: decchk
	INTEGER(I4B) :: i,j,k,m
	INTEGER(I4B) :: ip(0:9,0:7) = reshape((/ &
		0,1,2,3,4,5,6,7,8,9,1,5,7,6,2,8,3,0,9,4,&
		5,8,0,3,7,9,6,1,4,2,8,9,1,6,0,4,3,5,2,7,9,4,5,3,1,2,6,8,7,0,&
		4,2,8,6,5,7,3,9,0,1,2,7,9,3,8,0,6,4,1,5,7,0,4,6,9,1,3,2,5,8 /),&
		(/ 10,8 /) )
	INTEGER(I4B) :: ij(0:9,0:9) = reshape((/ &
		0,1,2,3,4,5,6,7,8,9,1,2,3,4,0,9,5,6,7,8,2,3,4,0,1,8,9,5,6,&
		7,3,4,0,1,2,7,8,9,5,6,4,0,1,2,3,6,7,8,9,5,5,6,7,8,9,0,1,2,3,&
		4,6,7,8,9,5,4,0,1,2,3,7,8,9,5,6,3,4,0,1,2,8,9,5,6,7,2,3,4,0,&
		1,9,5,6,7,8,1,2,3,4,0 /),(/ 10,10 /))
	k=0
	m=0
	do j=1,size(string)
		i=ichar(string(j))
		if (i >= 48 .and. i <= 57) then
			k=ij(k,ip(mod(i+2,10),mod(m,8)))
			m=m+1
		end if
	end do
	decchk=logical(k == 0,kind=lgt)
	i=mod(m,8)
	i=ifirstloc(ij(k,ip(0:9,i)) == 0)-1
	ch=char(i+48)
	END FUNCTION decchk

! -----------------------------------------------------------

	SUBROUTINE dfpmin(p,gtol,iter,fret,func,dfunc)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror,outerprod,unit_matrix,vabs
	use mo_nr, ONLY : lnsrch
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(SP), INTENT(IN) :: gtol
	REAL(SP), INTENT(OUT) :: fret
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
	INTERFACE
		FUNCTION func(p)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP) :: func
		END FUNCTION func
!BL
		FUNCTION dfunc(p)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP), DIMENSION(size(p)) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=200
	REAL(SP), PARAMETER :: STPMX=100.0_sp,EPS=epsilon(p),TOLX=4.0_sp*EPS
	INTEGER(I4B) :: its
	LOGICAL :: check
	REAL(SP) :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
	REAL(SP), DIMENSION(size(p)) :: dg,g,hdg,pnew,xi
	REAL(SP), DIMENSION(size(p),size(p)) :: hessin
	fp=func(p)
	g=dfunc(p)
	call unit_matrix(hessin)
	xi=-g
	stpmax=STPMX*max(vabs(p),real(size(p),sp))
	do its=1,ITMAX
		iter=its
		call lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,func)
		fp=fret
		xi=pnew-p
		p=pnew
		if (maxval(abs(xi)/max(abs(p),1.0_sp)) < TOLX) RETURN
		dg=g
		g=dfunc(p)
		den=max(fret,1.0_sp)
		if (maxval(abs(g)*max(abs(p),1.0_sp)/den) < gtol) RETURN
		dg=g-dg
		hdg=matmul(hessin,dg)
		fac=dot_product(dg,xi)
		fae=dot_product(dg,hdg)
		sumdg=dot_product(dg,dg)
		sumxi=dot_product(xi,xi)
		if (fac**2 > EPS*sumdg*sumxi) then
			fac=1.0_sp/fac
			fad=1.0_sp/fae
			dg=fac*xi-fad*hdg
			hessin=hessin+fac*outerprod(xi,xi)-&
				fad*outerprod(hdg,hdg)+fae*outerprod(dg,dg)
		end if
		xi=-matmul(hessin,g)
	end do
	call nrerror('dfpmin: too many iterations')
	END SUBROUTINE dfpmin

! -----------------------------------------------------------

	FUNCTION dfridr(func,x,h,err)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,geop,iminloc
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,h
	REAL(SP), INTENT(OUT) :: err
	REAL(SP) :: dfridr
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B),PARAMETER :: NTAB=10
	REAL(SP), PARAMETER :: CON=1.4_sp,CON2=CON*CON,BIG=huge(x),SAFE=2.0
	INTEGER(I4B) :: ierrmin,i,j
	REAL(SP) :: hh
	REAL(SP), DIMENSION(NTAB-1) :: errt,fac
	REAL(SP), DIMENSION(NTAB,NTAB) :: a
	call assert(h /= 0.0, 'dfridr arg')
	hh=h
	a(1,1)=(func(x+hh)-func(x-hh))/(2.0_sp*hh)
	err=BIG
	fac(1:NTAB-1)=geop(CON2,CON2,NTAB-1)
	do i=2,NTAB
		hh=hh/CON
		a(1,i)=(func(x+hh)-func(x-hh))/(2.0_sp*hh)
		do j=2,i
			a(j,i)=(a(j-1,i)*fac(j-1)-a(j-1,i-1))/(fac(j-1)-1.0_sp)
		end do
		errt(1:i-1)=max(abs(a(2:i,i)-a(1:i-1,i)),abs(a(2:i,i)-a(1:i-1,i-1)))
		ierrmin=iminloc(errt(1:i-1))
		if (errt(ierrmin) <= err) then
			err=errt(ierrmin)
			dfridr=a(1+ierrmin,i)
		end if
		if (abs(a(i,i)-a(i-1,i-1)) >= SAFE*err) RETURN
	end do
	END FUNCTION dfridr

! -----------------------------------------------------------

	SUBROUTINE dftcor(w,delta,a,b,endpts,corre,corim,corfac)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: w,delta,a,b
	REAL(SP), INTENT(OUT) :: corre,corim,corfac
	REAL(SP), DIMENSION(:), INTENT(IN) :: endpts
	REAL(SP) :: a0i,a0r,a1i,a1r,a2i,a2r,a3i,a3r,arg,c,cl,cr,s,sl,sr,t,&
		t2,t4,t6
	REAL(DP) :: cth,ctth,spth2,sth,sth4i,stth,th,th2,th4,&
		tmth2,tth4i
	th=w*delta
	call assert(a < b, th >= 0.0, th <= PI_D, 'dftcor args')
	if (abs(th) < 5.0e-2_dp) then
		t=th
		t2=t*t
		t4=t2*t2
		t6=t4*t2
		corfac=1.0_sp-(11.0_sp/720.0_sp)*t4+(23.0_sp/15120.0_sp)*t6
		a0r=(-2.0_sp/3.0_sp)+t2/45.0_sp+(103.0_sp/15120.0_sp)*t4-&
			(169.0_sp/226800.0_sp)*t6
		a1r=(7.0_sp/24.0_sp)-(7.0_sp/180.0_sp)*t2+(5.0_sp/3456.0_sp)*t4&
			-(7.0_sp/259200.0_sp)*t6
		a2r=(-1.0_sp/6.0_sp)+t2/45.0_sp-(5.0_sp/6048.0_sp)*t4+t6/64800.0_sp
		a3r=(1.0_sp/24.0_sp)-t2/180.0_sp+(5.0_sp/24192.0_sp)*t4-t6/259200.0_sp
		a0i=t*(2.0_sp/45.0_sp+(2.0_sp/105.0_sp)*t2-&
			(8.0_sp/2835.0_sp)*t4+(86.0_sp/467775.0_sp)*t6)
		a1i=t*(7.0_sp/72.0_sp-t2/168.0_sp+(11.0_sp/72576.0_sp)*t4-&
			(13.0_sp/5987520.0_sp)*t6)
		a2i=t*(-7.0_sp/90.0_sp+t2/210.0_sp-(11.0_sp/90720.0_sp)*t4+&
			(13.0_sp/7484400.0_sp)*t6)
		a3i=t*(7.0_sp/360.0_sp-t2/840.0_sp+(11.0_sp/362880.0_sp)*t4-&
			(13.0_sp/29937600.0_sp)*t6)
	else
		cth=cos(th)
		sth=sin(th)
		ctth=cth**2-sth**2
		stth=2.0_dp*sth*cth
		th2=th*th
		th4=th2*th2
		tmth2=3.0_dp-th2
		spth2=6.0_dp+th2
		sth4i=1.0_sp/(6.0_dp*th4)
		tth4i=2.0_dp*sth4i
		corfac=tth4i*spth2*(3.0_sp-4.0_dp*cth+ctth)
		a0r=sth4i*(-42.0_dp+5.0_dp*th2+spth2*(8.0_dp*cth-ctth))
		a0i=sth4i*(th*(-12.0_dp+6.0_dp*th2)+spth2*stth)
		a1r=sth4i*(14.0_dp*tmth2-7.0_dp*spth2*cth)
		a1i=sth4i*(30.0_dp*th-5.0_dp*spth2*sth)
		a2r=tth4i*(-4.0_dp*tmth2+2.0_dp*spth2*cth)
		a2i=tth4i*(-12.0_dp*th+2.0_dp*spth2*sth)
		a3r=sth4i*(2.0_dp*tmth2-spth2*cth)
		a3i=sth4i*(6.0_dp*th-spth2*sth)
	end if
	cl=a0r*endpts(1)+a1r*endpts(2)+a2r*endpts(3)+a3r*endpts(4)
	sl=a0i*endpts(1)+a1i*endpts(2)+a2i*endpts(3)+a3i*endpts(4)
	cr=a0r*endpts(8)+a1r*endpts(7)+a2r*endpts(6)+a3r*endpts(5)
	sr=-a0i*endpts(8)-a1i*endpts(7)-a2i*endpts(6)-a3i*endpts(5)
	arg=w*(b-a)
	c=cos(arg)
	s=sin(arg)
	corre=cl+c*cr-s*sr
	corim=sl+s*cr+c*sr
	END SUBROUTINE dftcor

! -----------------------------------------------------------

	SUBROUTINE dftint(func,a,b,w,cosint,sinint)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	use mo_nr, ONLY : dftcor,polint,realft
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b,w
	REAL(SP), INTENT(OUT) :: cosint,sinint
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: M=64,NDFT=1024,MPOL=6
	INTEGER(I4B) :: nn
	INTEGER(I4B), SAVE :: init=0
	INTEGER(I4B), DIMENSION(MPOL) :: nnmpol
	REAL(SP) :: c,cdft,cerr,corfac,corim,corre,en,s,sdft,serr
	REAL(SP), SAVE :: delta
	REAL(SP), DIMENSION(MPOL) :: cpol,spol,xpol
	REAL(SP), DIMENSION(NDFT), SAVE :: data
	REAL(SP), DIMENSION(8), SAVE :: endpts
	REAL(SP), SAVE :: aold=-1.0e30_sp,bold=-1.0e30_sp
	if (init /= 1 .or. a /= aold .or. b /= bold) then
		init=1
		aold=a
		bold=b
		delta=(b-a)/M
		data(1:M+1)=func(a+arth(0,1,M+1)*delta)
		data(M+2:NDFT)=0.0
		endpts(1:4)=data(1:4)
		endpts(5:8)=data(M-2:M+1)
		call realft(data(1:NDFT),1)
		data(2)=0.0
	end if
	en=w*delta*NDFT/TWOPI+1.0_sp
	nn=min(max(int(en-0.5_sp*MPOL+1.0_sp),1),NDFT/2-MPOL+1)
	nnmpol=arth(nn,1,MPOL)
	cpol(1:MPOL)=data(2*nnmpol(:)-1)
	spol(1:MPOL)=data(2*nnmpol(:))
	xpol(1:MPOL)=nnmpol(:)
	call polint(xpol,cpol,en,cdft,cerr)
	call polint(xpol,spol,en,sdft,serr)
	call dftcor(w,delta,a,b,endpts,corre,corim,corfac)
	cdft=cdft*corfac+corre
	sdft=sdft*corfac+corim
	c=delta*cos(w*a)
	s=delta*sin(w*a)
	cosint=c*cdft-s*sdft
	sinint=s*cdft+c*sdft
	END SUBROUTINE dftint

! -----------------------------------------------------------

	SUBROUTINE difeq(k,k1,k2,jsf,is1,isf,indexv,s,y)
	use mo_nrtype
	USE sfroid_data
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: is1,isf,jsf,k,k1,k2
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indexv
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: s
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: y
	REAL(SP) :: temp,temp2
	INTEGER(I4B), DIMENSION(3) :: indexv3
	indexv3(1:3)=3+indexv(1:3)
	if (k == k1) then
		if (mod(n+mm,2) == 1) then
			s(3,indexv3(1:3))= (/ 1.0_sp, 0.0_sp, 0.0_sp /)
			s(3,jsf)=y(1,1)
		else
			s(3,indexv3(1:3))= (/ 0.0_sp, 1.0_sp, 0.0_sp /)
			s(3,jsf)=y(2,1)
		end if
	else if (k > k2) then
		s(1,indexv3(1:3))= (/ -(y(3,M)-c2)/(2.0_sp*(mm+1.0_sp)),&
			1.0_sp, -y(1,M)/(2.0_sp*(mm+1.0_sp)) /)
		s(1,jsf)=y(2,M)-(y(3,M)-c2)*y(1,M)/(2.0_sp*(mm+1.0_sp))
		s(2,indexv3(1:3))=(/ 1.0_sp, 0.0_sp, 0.0_sp /)
		s(2,jsf)=y(1,M)-anorm
	else
		s(1,indexv(1:3))=(/ -1.0_sp, -0.5_sp*h, 0.0_sp /)
		s(1,indexv3(1:3))=(/ 1.0_sp, -0.5_sp*h, 0.0_sp /)
		temp=h/(1.0_sp-(x(k)+x(k-1))**2*0.25_sp)
		temp2=0.5_sp*(y(3,k)+y(3,k-1))-c2*0.25_sp*(x(k)+x(k-1))**2
		s(2,indexv(1:3))=(/ temp*temp2*0.5_sp,&
			-1.0_sp-0.5_sp*temp*(mm+1.0_sp)*(x(k)+x(k-1)),&
			0.25_sp*temp*(y(1,k)+y(1,k-1)) /)
		s(2,indexv3(1:3))=s(2,indexv(1:3))
		s(2,indexv3(2))=s(2,indexv3(2))+2.0_sp
		s(3,indexv(1:3))=(/ 0.0_sp, 0.0_sp, -1.0_sp /)
		s(3,indexv3(1:3))=(/ 0.0_sp, 0.0_sp, 1.0_sp /)
		s(1,jsf)=y(1,k)-y(1,k-1)-0.5_sp*h*(y(2,k)+y(2,k-1))
		s(2,jsf)=y(2,k)-y(2,k-1)-temp*((x(k)+x(k-1))*&
			0.5_sp*(mm+1.0_sp)*(y(2,k)+y(2,k-1))-temp2*&
			0.5_sp*(y(1,k)+y(1,k-1)))
		s(3,jsf)=y(3,k)-y(3,k-1)
	end if
	END SUBROUTINE difeq

! -----------------------------------------------------------

MODULE df1dim_mod
	use mo_nrtype
	INTEGER(I4B) :: ncom
	REAL(SP), DIMENSION(:), POINTER :: pcom,xicom
CONTAINS
!BL
	FUNCTION f1dim(x)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: f1dim
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), DIMENSION(:), ALLOCATABLE :: xt
	allocate(xt(ncom))
	xt(:)=pcom(:)+x*xicom(:)
	f1dim=func(xt)
	deallocate(xt)
	END FUNCTION f1dim
!BL
	FUNCTION df1dim(x)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: df1dim
	INTERFACE
		FUNCTION dfunc(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	REAL(SP), DIMENSION(:), ALLOCATABLE :: xt,df
	allocate(xt(ncom),df(ncom))
	xt(:)=pcom(:)+x*xicom(:)
	df(:)=dfunc(xt)
	df1dim=dot_product(df,xicom)
	deallocate(xt,df)
	END FUNCTION df1dim
END MODULE df1dim_mod

	SUBROUTINE dlinmin(p,xi,fret)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : mnbrak,dbrent
	USE df1dim_mod
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: fret
	REAL(SP), DIMENSION(:), TARGET :: p,xi
	REAL(SP), PARAMETER :: TOL=1.0e-4_sp
	REAL(SP) :: ax,bx,fa,fb,fx,xmin,xx
	ncom=assert_eq(size(p),size(xi),'dlinmin')
	pcom=>p
	xicom=>xi
	ax=0.0
	xx=1.0
	call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
	fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin)
	xi=xmin*xi
	p=p+xi
	END SUBROUTINE dlinmin

! -----------------------------------------------------------

	FUNCTION eclass(lista,listb,n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: lista,listb
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B), DIMENSION(n) :: eclass
	INTEGER :: j,k,l,m
	m=assert_eq(size(lista),size(listb),'eclass')
	eclass(1:n)=arth(1,1,n)
	do l=1,m
		j=lista(l)
		do
			if (eclass(j) == j) exit
			j=eclass(j)
		end do
		k=listb(l)
		do
			if (eclass(k) == k) exit
			k=eclass(k)
		end do
		if (j /= k) eclass(j)=k
	end do
	do j=1,n
		do
			if (eclass(j) == eclass(eclass(j))) exit
			eclass(j)=eclass(eclass(j))
		end do
	end do
	END FUNCTION eclass

! -----------------------------------------------------------

	FUNCTION eclazz(equiv,n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	INTERFACE
		FUNCTION equiv(i,j)
		use mo_nrtype
		IMPLICIT NONE
		LOGICAL(LGT) :: equiv
		INTEGER(I4B), INTENT(IN) :: i,j
		END FUNCTION equiv
	END INTERFACE
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B), DIMENSION(n) :: eclazz
	INTEGER :: i,j
	eclazz(1:n)=arth(1,1,n)
	do i=2,n
		do j=1,i-1
			eclazz(j)=eclazz(eclazz(j))
			if (equiv(i,j)) eclazz(eclazz(eclazz(j)))=i
		end do
	end do
	do i=1,n
		eclazz(i)=eclazz(eclazz(i))
	end do
	END FUNCTION eclazz

! -----------------------------------------------------------

	FUNCTION ei(x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: ei
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
	INTEGER(I4B) :: k
	REAL(SP) :: fact,prev,sm,term
	call assert(x > 0.0, 'ei arg')
	if (x < FPMIN) then
		ei=log(x)+EULER
	else if (x <= -log(EPS)) then
		sm=0.0
		fact=1.0
		do k=1,MAXIT
			fact=fact*x/k
			term=fact/k
			sm=sm+term
			if (term < EPS*sm) exit
		end do
		if (k > MAXIT) call nrerror('series failed in ei')
		ei=sm+log(x)+EULER
	else
		sm=0.0
		term=1.0
		do k=1,MAXIT
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
		if (k > MAXIT) call nrerror('asymptotic failed in ei')
		ei=exp(x)*(1.0_sp+sm)/x
	end if
	END FUNCTION ei

! -----------------------------------------------------------

	SUBROUTINE eigsrt(d,v)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,imaxloc,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: v
	INTEGER(I4B) :: i,j,n
	n=assert_eq(size(d),size(v,1),size(v,2),'eigsrt')
	do i=1,n-1
		j=imaxloc(d(i:n))+i-1
		if (j /= i) then
			call swap(d(i),d(j))
			call swap(v(:,i),v(:,j))
		end if
	end do
	END SUBROUTINE eigsrt

! -----------------------------------------------------------

	FUNCTION elle_s(phi,ak)
	use mo_nrtype
	use mo_nr, ONLY : rd,rf
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: phi,ak
	REAL(SP) :: elle_s
	REAL(SP) :: cc,q,s
	s=sin(phi)
	cc=cos(phi)**2
	q=(1.0_sp-s*ak)*(1.0_sp+s*ak)
	elle_s=s*(rf(cc,q,1.0_sp)-((s*ak)**2)*rd(cc,q,1.0_sp)/3.0_sp)
	END FUNCTION elle_s


	FUNCTION elle_v(phi,ak)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : rd,rf
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: phi,ak
	REAL(SP), DIMENSION(size(phi)) :: elle_v
	REAL(SP), DIMENSION(size(phi)) :: cc,q,s
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(phi),size(ak),'elle_v')
	s=sin(phi)
	cc=cos(phi)**2
	q=(1.0_sp-s*ak)*(1.0_sp+s*ak)
	elle_v=s*(rf(cc,q,spread(1.0_sp,1,size(phi)))-((s*ak)**2)*&
		rd(cc,q,spread(1.0_sp,1,size(phi)))/3.0_sp)
	END FUNCTION elle_v

! -----------------------------------------------------------

	FUNCTION ellf_s(phi,ak)
	use mo_nrtype
	use mo_nr, ONLY : rf
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: phi,ak
	REAL(SP) :: ellf_s
	REAL(SP) :: s
	s=sin(phi)
	ellf_s=s*rf(cos(phi)**2,(1.0_sp-s*ak)*(1.0_sp+s*ak),1.0_sp)
	END FUNCTION ellf_s


	FUNCTION ellf_v(phi,ak)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : rf
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: phi,ak
	REAL(SP), DIMENSION(size(phi)) :: ellf_v
	REAL(SP), DIMENSION(size(phi)) :: s
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(phi),size(ak),'ellf_v')
	s=sin(phi)
	ellf_v=s*rf(cos(phi)**2,(1.0_sp-s*ak)*(1.0_sp+s*ak),&
		spread(1.0_sp,1,size(phi)))
	END FUNCTION ellf_v

! -----------------------------------------------------------

	FUNCTION ellpi_s(phi,en,ak)
	use mo_nrtype
	use mo_nr, ONLY : rf,rj
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: phi,en,ak
	REAL(SP) :: ellpi_s
	REAL(SP) :: cc,enss,q,s
	s=sin(phi)
	enss=en*s*s
	cc=cos(phi)**2
	q=(1.0_sp-s*ak)*(1.0_sp+s*ak)
	ellpi_s=s*(rf(cc,q,1.0_sp)-enss*rj(cc,q,1.0_sp,1.0_sp+enss)/3.0_sp)
	END FUNCTION ellpi_s


	FUNCTION ellpi_v(phi,en,ak)
	use mo_nrtype
	use mo_nr, ONLY : rf,rj
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: phi,en,ak
	REAL(SP), DIMENSION(size(phi)) :: ellpi_v
	REAL(SP), DIMENSION(size(phi)) :: cc,enss,q,s
	s=sin(phi)
	enss=en*s*s
	cc=cos(phi)**2
	q=(1.0_sp-s*ak)*(1.0_sp+s*ak)
	ellpi_v=s*(rf(cc,q,spread(1.0_sp,1,size(phi)))-enss*&
		rj(cc,q,spread(1.0_sp,1,size(phi)),1.0_sp+enss)/3.0_sp)
	END FUNCTION ellpi_v

! -----------------------------------------------------------

	SUBROUTINE elmhes(a)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,imaxloc,outerprod,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B) :: i,m,n
	REAL(SP) :: x
	REAL(SP), DIMENSION(size(a,1)) :: y
	n=assert_eq(size(a,1),size(a,2),'elmhes')
	do m=2,n-1
		i=imaxloc(abs(a(m:n,m-1)))+m-1
		x=a(i,m-1)
		if (i /= m) then
			call swap(a(i,m-1:n),a(m,m-1:n))
			call swap(a(:,i),a(:,m))
		end if
		if (x /= 0.0) then
			y(m+1:n)=a(m+1:n,m-1)/x
			a(m+1:n,m-1)=y(m+1:n)
			a(m+1:n,m:n)=a(m+1:n,m:n)-outerprod(y(m+1:n),a(m,m:n))
			a(:,m)=a(:,m)+matmul(a(:,m+1:n),y(m+1:n))
		end if
	end do
	END SUBROUTINE elmhes

! -----------------------------------------------------------

	FUNCTION erf_s(x)
	use mo_nrtype
	use mo_nr, ONLY : gammp
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: erf_s
	erf_s=gammp(0.5_sp,x**2)
	if (x < 0.0) erf_s=-erf_s
	END FUNCTION erf_s


	FUNCTION erf_v(x)
	use mo_nrtype
	use mo_nr, ONLY : gammp
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: erf_v
	erf_v=gammp(spread(0.5_sp,1,size(x)),x**2)
	where (x < 0.0) erf_v=-erf_v
	END FUNCTION erf_v

! -----------------------------------------------------------

	FUNCTION erfc_s(x)
	use mo_nrtype
	use mo_nr, ONLY : gammp,gammq
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: erfc_s
	erfc_s=merge(1.0_sp+gammp(0.5_sp,x**2),gammq(0.5_sp,x**2), x < 0.0)
	END FUNCTION erfc_s


	FUNCTION erfc_v(x)
	use mo_nrtype
	use mo_nr, ONLY : gammp,gammq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: erfc_v
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	mask = (x < 0.0)
	erfc_v=merge(1.0_sp+gammp(spread(0.5_sp,1,size(x)), &
		merge(x,0.0_sp,mask)**2),gammq(spread(0.5_sp,1,size(x)), &
		merge(x,0.0_sp,.not. mask)**2),mask)
	END FUNCTION erfc_v

! -----------------------------------------------------------

	FUNCTION erfcc_s(x)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: erfcc_s
	REAL(SP) :: t,z
	REAL(SP), DIMENSION(10) :: coef = (/-1.26551223_sp,1.00002368_sp,&
		0.37409196_sp,0.09678418_sp,-0.18628806_sp,0.27886807_sp,&
		-1.13520398_sp,1.48851587_sp,-0.82215223_sp,0.17087277_sp/)
	z=abs(x)
	t=1.0_sp/(1.0_sp+0.5_sp*z)
	erfcc_s=t*exp(-z*z+poly(t,coef))
	if (x < 0.0) erfcc_s=2.0_sp-erfcc_s
	END FUNCTION erfcc_s


	FUNCTION erfcc_v(x)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: erfcc_v,t,z
	REAL(SP), DIMENSION(10) :: coef = (/-1.26551223_sp,1.00002368_sp,&
		0.37409196_sp,0.09678418_sp,-0.18628806_sp,0.27886807_sp,&
		-1.13520398_sp,1.48851587_sp,-0.82215223_sp,0.17087277_sp/)
	z=abs(x)
	t=1.0_sp/(1.0_sp+0.5_sp*z)
	erfcc_v=t*exp(-z*z+poly(t,coef))
	where (x < 0.0) erfcc_v=2.0_sp-erfcc_v
	END FUNCTION erfcc_v

! -----------------------------------------------------------

	SUBROUTINE eulsum(sum,term,jterm)
	use mo_nrtype
          use mo_nrutil, ONLY : poly_term,reallocate
	IMPLICIT NONE
	REAL(SP), INTENT(INOUT) :: sum
	REAL(SP), INTENT(IN) :: term
	INTEGER(I4B), INTENT(IN) :: jterm
	REAL(SP), DIMENSION(:), POINTER, SAVE :: wksp
	INTEGER(I4B), SAVE :: nterm
	LOGICAL(LGT), SAVE :: init=.true.
	if (init) then
		init=.false.
		nullify(wksp)
	end if
	if (jterm == 1) then
		nterm=1
		wksp=>reallocate(wksp,100)
		wksp(1)=term
		sum=0.5_sp*term
	else
		if (nterm+1 > size(wksp)) wksp=>reallocate(wksp,2*size(wksp))
		wksp(2:nterm+1)=0.5_sp*wksp(1:nterm)
		wksp(1)=term
		wksp(1:nterm+1)=poly_term(wksp(1:nterm+1),0.5_sp)
		if (abs(wksp(nterm+1)) <= abs(wksp(nterm))) then
			sum=sum+0.5_sp*wksp(nterm+1)
			nterm=nterm+1
		else
			sum=sum+wksp(nterm+1)
		end if
	end if
	END SUBROUTINE eulsum

! -----------------------------------------------------------

	FUNCTION evlmem(fdt,d,xms)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: fdt,xms
	REAL(SP), DIMENSION(:), INTENT(IN) :: d
	REAL(SP) :: evlmem
	COMPLEX(SPC) :: z,zz
	REAL(DP) :: theta
	theta=TWOPI_D*fdt
	z=cmplx(cos(theta),sin(theta),kind=spc)
	zz=1.0_sp-z*poly(z,d)
	evlmem=xms/abs(zz)**2
	END FUNCTION evlmem

! -----------------------------------------------------------

	SUBROUTINE expdev_s(harvest)
	use mo_nrtype
	use mo_nr, ONLY : ran1
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: harvest
	REAL(SP) :: dum
	call ran1(dum)
	harvest=-log(dum)
	END SUBROUTINE expdev_s

	SUBROUTINE expdev_v(harvest)
	use mo_nrtype
	use mo_nr, ONLY : ran1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
	REAL(SP), DIMENSION(size(harvest)) :: dum
	call ran1(dum)
	harvest=-log(dum)
	END SUBROUTINE expdev_v

! -----------------------------------------------------------

	FUNCTION expint(n,x)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert,nrerror
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: expint
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(SP), PARAMETER :: EPS=epsilon(x),BIG=huge(x)*EPS
	INTEGER(I4B) :: i,nm1
	REAL(SP) :: a,b,c,d,del,fact,h
	call assert(n >= 0, x >= 0.0, (x > 0.0 .or. n > 1), &
		'expint args')
	if (n == 0) then
		expint=exp(-x)/x
		RETURN
	end if
	nm1=n-1
	if (x == 0.0) then
		expint=1.0_sp/nm1
	else if (x > 1.0) then
		b=x+n
		c=BIG
		d=1.0_sp/b
		h=d
		do i=1,MAXIT
			a=-i*(nm1+i)
			b=b+2.0_sp
			d=1.0_sp/(a*d+b)
			c=b+a/c
			del=c*d
			h=h*del
			if (abs(del-1.0_sp) <= EPS) exit
		end do
		if (i > MAXIT) call nrerror('expint: continued fraction failed')
		expint=h*exp(-x)
	else
		if (nm1 /= 0) then
			expint=1.0_sp/nm1
		else
			expint=-log(x)-EULER
		end if
		fact=1.0
		do i=1,MAXIT
			fact=-fact*x/i
			if (i /= nm1) then
				del=-fact/(i-nm1)
			else
				del=fact*(-log(x)-EULER+sum(1.0_sp/arth(1,1,nm1)))
			end if
			expint=expint+del
			if (abs(del) < abs(expint)*EPS) exit
		end do
		if (i > MAXIT) call nrerror('expint: series failed')
	end if
	END FUNCTION expint

! -----------------------------------------------------------

	FUNCTION factln_s(n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert
	use mo_nr, ONLY : gammln
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP) :: factln_s
	INTEGER(I4B), PARAMETER :: TMAX=100
	REAL(SP), DIMENSION(TMAX), SAVE :: a
	LOGICAL(LGT), SAVE :: init=.true.
	if (init) then
		a(1:TMAX)=gammln(arth(1.0_sp,1.0_sp,TMAX))
		init=.false.
	end if
	call assert(n >= 0, 'factln_s arg')
	if (n < TMAX) then
		factln_s=a(n+1)
	else
		factln_s=gammln(n+1.0_sp)
	end if
	END FUNCTION factln_s


	FUNCTION factln_v(n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert
	use mo_nr, ONLY : gammln
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
	REAL(SP), DIMENSION(size(n)) :: factln_v
	LOGICAL(LGT), DIMENSION(size(n)) :: mask
	INTEGER(I4B), PARAMETER :: TMAX=100
	REAL(SP), DIMENSION(TMAX), SAVE :: a
	LOGICAL(LGT), SAVE :: init=.true.
	if (init) then
		a(1:TMAX)=gammln(arth(1.0_sp,1.0_sp,TMAX))
		init=.false.
	end if
	call assert(all(n >= 0), 'factln_v arg')
	mask = (n >= TMAX)
	factln_v=unpack(gammln(pack(n,mask)+1.0_sp),mask,0.0_sp)
	where (.not. mask) factln_v=a(n+1)
	END FUNCTION factln_v

! -----------------------------------------------------------

	FUNCTION factrl_s(n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert,cumprod
	use mo_nr, ONLY : gammln
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP) :: factrl_s
	INTEGER(I4B), SAVE :: ntop=0
	INTEGER(I4B), PARAMETER :: NMAX=32
	REAL(SP), DIMENSION(NMAX), SAVE :: a
	call assert(n >= 0, 'factrl_s arg')
	if (n < ntop) then
		factrl_s=a(n+1)
	else if (n < NMAX) then
		ntop=NMAX
		a(1)=1.0
		a(2:NMAX)=cumprod(arth(1.0_sp,1.0_sp,NMAX-1))
		factrl_s=a(n+1)
	else
		factrl_s=exp(gammln(n+1.0_sp))
	end if
	END FUNCTION factrl_s

	FUNCTION factrl_v(n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert,cumprod
	use mo_nr, ONLY : gammln
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
	REAL(SP), DIMENSION(size(n)) :: factrl_v
	LOGICAL(LGT), DIMENSION(size(n)) :: mask
	INTEGER(I4B), SAVE :: ntop=0
	INTEGER(I4B), PARAMETER :: NMAX=32
	REAL(SP), DIMENSION(NMAX), SAVE :: a
	call assert(all(n >= 0), 'factrl_v arg')
	if (ntop == 0) then
		ntop=NMAX
		a(1)=1.0
		a(2:NMAX)=cumprod(arth(1.0_sp,1.0_sp,NMAX-1))
	end if
	mask = (n >= NMAX)
	factrl_v=unpack(exp(gammln(pack(n,mask)+1.0_sp)),mask,0.0_sp)
	where (.not. mask) factrl_v=a(n+1)
	END FUNCTION factrl_v

! -----------------------------------------------------------

	SUBROUTINE fasper(x,y,ofac,hifac,px,py,jmax,prob)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq,imaxloc,nrerror
	use mo_nr, ONLY : avevar,realft
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), INTENT(IN) :: ofac,hifac
	INTEGER(I4B), INTENT(OUT) :: jmax
	REAL(SP), INTENT(OUT) :: prob
	REAL(SP), DIMENSION(:), POINTER :: px,py
	INTEGER(I4B), PARAMETER :: MACC=4
	INTEGER(I4B) :: j,k,n,ndim,nfreq,nfreqt,nout
	REAL(SP) :: ave,ck,ckk,cterm,cwt,den,df,effm,expy,fac,fndim,hc2wt,&
		hs2wt,hypo,sterm,swt,var,xdif,xmax,xmin
	REAL(SP), DIMENSION(:), ALLOCATABLE :: wk1,wk2
	LOGICAL(LGT), SAVE :: init=.true.
	n=assert_eq(size(x),size(y),'fasper')
	if (init) then
		init=.false.
		nullify(px,py)
	else
		if (associated(px)) deallocate(px)
		if (associated(py)) deallocate(py)
	end if
	nfreqt=ofac*hifac*n*MACC
	nfreq=64
	do
		if (nfreq >= nfreqt) exit
		nfreq=nfreq*2
	end do
	ndim=2*nfreq
	allocate(wk1(ndim),wk2(ndim))
	call avevar(y(1:n),ave,var)
	xmax=maxval(x(:))
	xmin=minval(x(:))
	xdif=xmax-xmin
	wk1(1:ndim)=0.0
	wk2(1:ndim)=0.0
	fac=ndim/(xdif*ofac)
	fndim=ndim
	do j=1,n
		ck=1.0_sp+mod((x(j)-xmin)*fac,fndim)
		ckk=1.0_sp+mod(2.0_sp*(ck-1.0_sp),fndim)
		call spreadval(y(j)-ave,wk1,ck,MACC)
		call spreadval(1.0_sp,wk2,ckk,MACC)
	end do
	call realft(wk1(1:ndim),1)
	call realft(wk2(1:ndim),1)
	df=1.0_sp/(xdif*ofac)
	nout=0.5_sp*ofac*hifac*n
	allocate(px(nout),py(nout))
	k=3
	do j=1,nout
		hypo=sqrt(wk2(k)**2+wk2(k+1)**2)
		hc2wt=0.5_sp*wk2(k)/hypo
		hs2wt=0.5_sp*wk2(k+1)/hypo
		cwt=sqrt(0.5_sp+hc2wt)
		swt=sign(sqrt(0.5_sp-hc2wt),hs2wt)
		den=0.5_sp*n+hc2wt*wk2(k)+hs2wt*wk2(k+1)
		cterm=(cwt*wk1(k)+swt*wk1(k+1))**2/den
		sterm=(cwt*wk1(k+1)-swt*wk1(k))**2/(n-den)
		px(j)=j*df
		py(j)=(cterm+sterm)/(2.0_sp*var)
		k=k+2
	end do
	deallocate(wk1,wk2)
	jmax=imaxloc(py(1:nout))
	expy=exp(-py(jmax))
	effm=2.0_sp*nout/ofac
	prob=effm*expy
	if (prob > 0.01_sp) prob=1.0_sp-(1.0_sp-expy)**effm
	CONTAINS
!BL
	SUBROUTINE spreadval(y,yy,x,m)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: y,x
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: yy
	INTEGER(I4B), INTENT(IN) :: m
	INTEGER(I4B) :: ihi,ilo,ix,j,nden,n
	REAL(SP) :: fac
	INTEGER(I4B), DIMENSION(10) :: nfac = (/ &
		1,1,2,6,24,120,720,5040,40320,362880 /)
	if (m > 10) call nrerror('factorial table too small in spreadval')
	n=size(yy)
	ix=x
	if (x == real(ix,sp)) then
		yy(ix)=yy(ix)+y
	else
		ilo=min(max(int(x-0.5_sp*m+1.0_sp),1),n-m+1)
		ihi=ilo+m-1
		nden=nfac(m)
		fac=product(x-arth(ilo,1,m))
		yy(ihi)=yy(ihi)+y*fac/(nden*(x-ihi))
		do j=ihi-1,ilo,-1
			nden=(nden/(j+1-ilo))*(j-ihi)
			yy(j)=yy(j)+y*fac/(nden*(x-j))
		end do
	end if
	END SUBROUTINE spreadval
	END SUBROUTINE fasper

! -----------------------------------------------------------

	SUBROUTINE fdjac(x,fvec,df)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: fvec
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: df
	INTERFACE
		FUNCTION funcv(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: funcv
		END FUNCTION funcv
	END INTERFACE
	REAL(SP), PARAMETER :: EPS=1.0e-4_sp
	INTEGER(I4B) :: j,n
	REAL(SP), DIMENSION(size(x)) :: xsav,xph,h
	n=assert_eq(size(x),size(fvec),size(df,1),size(df,2),'fdjac')
	xsav=x
	h=EPS*abs(xsav)
	where (h == 0.0) h=EPS
	xph=xsav+h
	h=xph-xsav
	do j=1,n
		x(j)=xph(j)
		df(:,j)=(funcv(x)-fvec(:))/h(j)
		x(j)=xsav(j)
	end do
	END SUBROUTINE fdjac

! -----------------------------------------------------------

	SUBROUTINE fgauss(x,a,y,dyda)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: y
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
	INTEGER(I4B) :: i,na,nx
	REAL(SP), DIMENSION(size(x)) :: arg,ex,fac
	nx=assert_eq(size(x),size(y),size(dyda,1),'fgauss: nx')
	na=assert_eq(size(a),size(dyda,2),'fgauss: na')
	y(:)=0.0
	do i=1,na-1,3
		arg(:)=(x(:)-a(i+1))/a(i+2)
		ex(:)=exp(-arg(:)**2)
		fac(:)=a(i)*ex(:)*2.0_sp*arg(:)
		y(:)=y(:)+a(i)*ex(:)
		dyda(:,i)=ex(:)
		dyda(:,i+1)=fac(:)/a(i+2)
		dyda(:,i+2)=fac(:)*arg(:)/a(i+2)
	end do
	END SUBROUTINE fgauss

! -----------------------------------------------------------

	SUBROUTINE fit(x,y,a,b,siga,sigb,chi2,q,sig)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : gammq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
	REAL(SP), DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
	INTEGER(I4B) :: ndum
	REAL(SP) :: sigdat,ss,sx,sxoss,sy,st2
	REAL(SP), DIMENSION(size(x)), TARGET :: t
	REAL(SP), DIMENSION(:), POINTER :: wt
	if (present(sig)) then
		ndum=assert_eq(size(x),size(y),size(sig),'fit')
		wt=>t
		wt(:)=1.0_sp/(sig(:)**2)
		ss=sum(wt(:))
		sx=dot_product(wt,x)
		sy=dot_product(wt,y)
	else
		ndum=assert_eq(size(x),size(y),'fit')
		ss=real(size(x),sp)
		sx=sum(x)
		sy=sum(y)
	end if
	sxoss=sx/ss
	t(:)=x(:)-sxoss
	if (present(sig)) then
		t(:)=t(:)/sig(:)
		b=dot_product(t/sig,y)
	else
		b=dot_product(t,y)
	end if
	st2=dot_product(t,t)
	b=b/st2
	a=(sy-sx*b)/ss
	siga=sqrt((1.0_sp+sx*sx/(ss*st2))/ss)
	sigb=sqrt(1.0_sp/st2)
	t(:)=y(:)-a-b*x(:)
	if (present(sig)) then
		t(:)=t(:)/sig(:)
		chi2=dot_product(t,t)
		q=gammq(0.5_sp*(size(x)-2),0.5_sp*chi2)
	else
		chi2=dot_product(t,t)
		q=1.0
		sigdat=sqrt(chi2/(size(x)-2))
		siga=siga*sigdat
		sigb=sigb*sigdat
	end if
	END SUBROUTINE fit

! -----------------------------------------------------------

	SUBROUTINE fitexy(x,y,sigx,sigy,a,b,siga,sigb,chi2,q)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,swap
	use mo_nr, ONLY : avevar,brent,fit,gammq,mnbrak,zbrent
	USE chixyfit
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sigx,sigy
	REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
	REAL(SP), PARAMETER :: POTN=1.571000_sp,BIG=1.0e30_sp,ACC=1.0e-3_sp
	INTEGER(I4B) :: j,n
	REAL(SP), DIMENSION(size(x)), TARGET :: xx,yy,sx,sy,ww
	REAL(SP), DIMENSION(6) :: ang,ch
	REAL(SP) :: amx,amn,varx,vary,scale,bmn,bmx,d1,d2,r2,&
		dum1,dum2,dum3,dum4,dum5
	n=assert_eq(size(x),size(y),size(sigx),size(sigy),'fitexy')
	xxp=>xx
	yyp=>yy
	sxp=>sx
	syp=>sy
	wwp=>ww
	call avevar(x,dum1,varx)
	call avevar(y,dum1,vary)
	scale=sqrt(varx/vary)
	xx(:)=x(:)
	yy(:)=y(:)*scale
	sx(:)=sigx(:)
	sy(:)=sigy(:)*scale
	ww(:)=sqrt(sx(:)**2+sy(:)**2)
	call fit(xx,yy,dum1,b,dum2,dum3,dum4,dum5,ww)
	offs=0.0
	ang(1)=0.0
	ang(2)=atan(b)
	ang(4)=0.0
	ang(5)=ang(2)
	ang(6)=POTN
	do j=4,6
		ch(j)=chixy(ang(j))
	end do
	call mnbrak(ang(1),ang(2),ang(3),ch(1),ch(2),ch(3),chixy)
	chi2=brent(ang(1),ang(2),ang(3),chixy,ACC,b)
	chi2=chixy(b)
	a=aa
	q=gammq(0.5_sp*(n-2),0.5_sp*chi2)
	r2=1.0_sp/sum(ww(:))
	bmx=BIG
	bmn=BIG
	offs=chi2+1.0_sp
	do j=1,6
		if (ch(j) > offs) then
			d1=mod(abs(ang(j)-b),PI)
			d2=PI-d1
			if (ang(j) < b) call swap(d1,d2)
			if (d1 < bmx) bmx=d1
			if (d2 < bmn) bmn=d2
		end if
	end do
	if (bmx <  BIG) then
		bmx=zbrent(chixy,b,b+bmx,ACC)-b
		amx=aa-a
		bmn=zbrent(chixy,b,b-bmn,ACC)-b
		amn=aa-a
		sigb=sqrt(0.5_sp*(bmx**2+bmn**2))/(scale*cos(b)**2)
		siga=sqrt(0.5_sp*(amx**2+amn**2)+r2)/scale
	else
		sigb=BIG
		siga=BIG
	end if
	a=a/scale
	b=tan(b)/scale
	END SUBROUTINE fitexy

! -----------------------------------------------------------

	SUBROUTINE fixrts(d)
	use mo_nrtype
	use mo_nr, ONLY : zroots
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
	INTEGER(I4B) :: i,m
	LOGICAL(LGT) :: polish
	COMPLEX(SPC), DIMENSION(size(d)+1) :: a
	COMPLEX(SPC), DIMENSION(size(d)) :: roots
	m=size(d)
	a(m+1)=cmplx(1.0_sp,kind=spc)
	a(m:1:-1)=cmplx(-d(1:m),kind=spc)
	polish=.true.
	call zroots(a(1:m+1),roots,polish)
	where (abs(roots) > 1.0) roots=1.0_sp/conjg(roots)
	a(1)=-roots(1)
	a(2:m+1)=cmplx(1.0_sp,kind=spc)
	do i=2,m
		a(2:i)=a(1:i-1)-roots(i)*a(2:i)
		a(1)=-roots(i)*a(1)
	end do
	d(m:1:-1)=-real(a(1:m))
	END SUBROUTINE fixrts

! -----------------------------------------------------------

	FUNCTION fleg(x,nl)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	INTEGER(I4B), INTENT(IN) :: nl
	REAL(SP), DIMENSION(nl) :: fleg
	INTEGER(I4B) :: j
	REAL(SP) :: d,f1,f2,twox
	fleg(1)=1.0
	fleg(2)=x
	if (nl > 2) then
		twox=2.0_sp*x
		f2=x
		d=1.0
		do j=3,nl
			f1=d
			f2=f2+twox
			d=d+1.0_sp
			fleg(j)=(f2*fleg(j-1)-f1*fleg(j-2))/d
		end do
	end if
	END FUNCTION fleg

! -----------------------------------------------------------

	SUBROUTINE flmoon(n,nph,jd,frac)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n,nph
	INTEGER(I4B), INTENT(OUT) :: jd
	REAL(SP), INTENT(OUT) :: frac
	REAL(SP), PARAMETER :: RAD=PI/180.0_sp
	INTEGER(I4B) :: i
	REAL(SP) :: am,as,c,t,t2,xtra
	c=n+nph/4.0_sp
	t=c/1236.85_sp
	t2=t**2
	as=359.2242_sp+29.105356_sp*c
	am=306.0253_sp+385.816918_sp*c+0.010730_sp*t2
	jd=2415020+28*n+7*nph
	xtra=0.75933_sp+1.53058868_sp*c+(1.178e-4_sp-1.55e-7_sp*t)*t2
	select case(nph)
		case(0,2)
			xtra=xtra+(0.1734_sp-3.93e-4_sp*t)*sin(RAD*as)-0.4068_sp*sin(RAD*am)
		case(1,3)
			xtra=xtra+(0.1721_sp-4.0e-4_sp*t)*sin(RAD*as)-0.6280_sp*sin(RAD*am)
		case default
			call nrerror('flmoon: nph is unknown')
	end select
	i=int(merge(xtra,xtra-1.0_sp, xtra >= 0.0))
	jd=jd+i
	frac=xtra-i
	END SUBROUTINE flmoon

! -----------------------------------------------------------

MODULE fminln
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	REAL(SP), DIMENSION(:), POINTER :: fmin_fvecp
CONTAINS
!BL
	FUNCTION fmin(x)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP) :: fmin
	INTERFACE
		FUNCTION funcv(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: funcv
		END FUNCTION funcv
	END INTERFACE
	if (.not. associated(fmin_fvecp)) call &
		nrerror('fmin: problem with pointer for returned values')
	fmin_fvecp=funcv(x)
	fmin=0.5_sp*dot_product(fmin_fvecp,fmin_fvecp)
	END FUNCTION fmin
END MODULE fminln

! -----------------------------------------------------------

	SUBROUTINE four1_sp(data,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert
	use mo_nr, ONLY : fourrow
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
	COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
	REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
	INTEGER(I4B) :: n,m1,m2,j
	n=size(data)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_sp')
	m1=2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
	m2=n/m1
	allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
	dat=reshape(data,shape(dat))
	call fourrow(dat,isign)
	theta=arth(0,isign,m1)*TWOPI_D/n
	wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
	w=cmplx(1.0_dp,0.0_dp,kind=dpc)
	do j=2,m2
		w=w*wp+w
		dat(:,j)=dat(:,j)*w
	end do
	temp=transpose(dat)
	call fourrow(temp,isign)
	data=reshape(temp,shape(data))
	deallocate(dat,w,wp,theta,temp)
	END SUBROUTINE four1_sp

	SUBROUTINE four1_dp(data,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert
	use mo_nr, ONLY : fourrow
	IMPLICIT NONE
	COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
	COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
	REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
	INTEGER(I4B) :: n,m1,m2,j
	n=size(data)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_dp')
	m1=2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
	m2=n/m1
	allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
	dat=reshape(data,shape(dat))
	call fourrow(dat,isign)
	theta=arth(0,isign,m1)*TWOPI_D/n
	wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
	w=cmplx(1.0_dp,0.0_dp,kind=dpc)
	do j=2,m2
		w=w*wp+w
		dat(:,j)=dat(:,j)*w
	end do
	temp=transpose(dat)
	call fourrow(temp,isign)
	data=reshape(temp,shape(data))
	deallocate(dat,w,wp,theta,temp)
	END SUBROUTINE four1_dp

! -----------------------------------------------------------

	SUBROUTINE four1_alt(data,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert
	use mo_nr, ONLY : fourcol
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
	COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
	REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
	INTEGER(I4B) :: n,m1,m2,j
	n=size(data)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_alt')
	m1=2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
	m2=n/m1
	allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
	dat=reshape(data,shape(dat))
	temp=transpose(dat)
	call fourcol(temp,isign)
	theta=arth(0,isign,m1)*TWOPI_D/n
	wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
	w=cmplx(1.0_dp,0.0_dp,kind=dpc)
	do j=2,m2
		w=w*wp+w
		temp(j,:)=temp(j,:)*w
	end do
	dat=transpose(temp)
	call fourcol(dat,isign)
	temp=transpose(dat)
	data=reshape(temp,shape(data))
	deallocate(dat,w,wp,theta,temp)
	END SUBROUTINE four1_alt

! -----------------------------------------------------------

	SUBROUTINE four1_gather(data,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER(I4B) :: n,n2,m,mm
	REAL(DP) :: theta
	COMPLEX(SPC) :: wp
	INTEGER(I4B), DIMENSION(size(data)) :: jarr
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: jrev
	COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: wtab,dtemp
	n=size(data)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_gather')
	if (n <= 1) RETURN
	allocate(jrev(n))
	jarr=arth(0,1,n)
	jrev=0
	n2=n/2
	m=n2
	do
		where (iand(jarr,1) /= 0) jrev=jrev+m
		jarr=jarr/2
		m=m/2
		if (m == 0) exit
	end do
	data=data(jrev+1)
	deallocate(jrev)
	allocate(dtemp(n),wtab(n2))
	jarr=arth(0,1,n)
	m=1
	mm=n2
	wtab(1)=(1.0_sp,0.0_sp)
	do
		where (iand(jarr,m) /= 0)
			dtemp=data*wtab(mm*iand(jarr,m-1)+1)
			data=eoshift(data,-m)-dtemp
		elsewhere
			data=data+eoshift(dtemp,m)
		end where
		m=m*2
		if (m >= n) exit
		mm=mm/2
		theta=PI_D/(isign*m)
		wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2, sin(theta),kind=spc)
		wtab(mm+1:n2:2*mm)=wtab(1:n2-mm:2*mm)*wp+wtab(1:n2-mm:2*mm)
	end do
	deallocate(dtemp,wtab)
	END SUBROUTINE four1_gather

! -----------------------------------------------------------

	SUBROUTINE four2(data,isign)
	use mo_nrtype
	use mo_nr, ONLY : fourrow
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(size(data,2),size(data,1)) :: temp
	call fourrow(data,isign)
	temp=transpose(data)
	call fourrow(temp,isign)
	data=transpose(temp)
	END SUBROUTINE four2

! -----------------------------------------------------------

	SUBROUTINE four2_alt(data,isign)
	use mo_nrtype
	use mo_nr, ONLY : fourcol
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(size(data,2),size(data,1)) :: temp
	temp=transpose(data)
	call fourcol(temp,isign)
	data=transpose(temp)
	call fourcol(data,isign)
	END SUBROUTINE four2_alt

! -----------------------------------------------------------

	SUBROUTINE four3(data,isign)
	use mo_nrtype
	use mo_nr, ONLY : fourrow_3d
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(:,:,:), ALLOCATABLE :: dat2,dat3
	call fourrow_3d(data,isign)
	allocate(dat2(size(data,2),size(data,3),size(data,1)))
	dat2=reshape(data,shape=shape(dat2),order=(/3,1,2/))
	call fourrow_3d(dat2,isign)
	allocate(dat3(size(data,3),size(data,1),size(data,2)))
	dat3=reshape(dat2,shape=shape(dat3),order=(/3,1,2/))
	deallocate(dat2)
	call fourrow_3d(dat3,isign)
	data=reshape(dat3,shape=shape(data),order=(/3,1,2/))
	deallocate(dat3)
	END SUBROUTINE four3

! -----------------------------------------------------------

	SUBROUTINE four3_alt(data,isign)
	use mo_nrtype
	use mo_nr, ONLY : fourcol_3d
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(:,:,:), ALLOCATABLE :: dat2,dat3
	call fourcol_3d(data,isign)
	allocate(dat2(size(data,2),size(data,3),size(data,1)))
	dat2=reshape(data,shape=shape(dat2),order=(/3,1,2/))
	call fourcol_3d(dat2,isign)
	allocate(dat3(size(data,3),size(data,1),size(data,2)))
	dat3=reshape(dat2,shape=shape(dat3),order=(/3,1,2/))
	deallocate(dat2)
	call fourcol_3d(dat3,isign)
	data=reshape(dat3,shape=shape(data),order=(/3,1,2/))
	deallocate(dat3)
	END SUBROUTINE four3_alt

! -----------------------------------------------------------

	SUBROUTINE fourcol(data,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,swap
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
	REAL(DP) :: theta
	COMPLEX(SPC), DIMENSION(size(data,2)) :: temp
	COMPLEX(DPC) :: w,wp
	COMPLEX(SPC) :: ws
	n=size(data,1)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourcol')
	n2=n/2
	j=n2
	do i=1,n-2
		if (j > i) call swap(data(j+1,:),data(i+1,:))
		m=n2
		do
			if (m < 2 .or. j < m) exit
			j=j-m
			m=m/2
		end do
		j=j+m
	end do
	mmax=1
	do
		if (n <= mmax) exit
		istep=2*mmax
		theta=PI_D/(isign*mmax)
		wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
		w=cmplx(1.0_dp,0.0_dp,kind=dpc)
		do m=1,mmax
			ws=w
			do i=m,n,istep
				j=i+mmax
				temp=ws*data(j,:)
				data(j,:)=data(i,:)-temp
				data(i,:)=data(i,:)+temp
			end do
			w=w*wp+w
		end do
		mmax=istep
	end do
	END SUBROUTINE fourcol

! -----------------------------------------------------------

	SUBROUTINE fourcol_3d(data,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,swap
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
	REAL(DP) :: theta
	COMPLEX(SPC), DIMENSION(size(data,2),size(data,3)) :: temp
	COMPLEX(DPC) :: w,wp
	COMPLEX(SPC) :: ws
	n=size(data,1)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourcol_3d')
	n2=n/2
	j=n2
	do i=1,n-2
		if (j > i) call swap(data(j+1,:,:),data(i+1,:,:))
		m=n2
		do
			if (m < 2 .or. j < m) exit
			j=j-m
			m=m/2
		end do
		j=j+m
	end do
	mmax=1
	do
		if (n <= mmax) exit
		istep=2*mmax
		theta=PI_D/(isign*mmax)
		wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
		w=cmplx(1.0_dp,0.0_dp,kind=dpc)
		do m=1,mmax
			ws=w
			do i=m,n,istep
				j=i+mmax
				temp=ws*data(j,:,:)
				data(j,:,:)=data(i,:,:)-temp
				data(i,:,:)=data(i,:,:)+temp
			end do
			w=w*wp+w
		end do
		mmax=istep
	end do
	END SUBROUTINE fourcol_3d

! -----------------------------------------------------------

	SUBROUTINE fourn_gather(data,nn,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
	INTEGER(I4B), DIMENSION(:) :: nn
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: jarr
	INTEGER(I4B) :: ndim,idim,ntot,nprev,n,n2,msk0,msk1,msk2,m,mm,mn
	REAL(DP) :: theta
	COMPLEX(SPC) :: wp
	COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: wtab,dtemp
	call assert(iand(nn,nn-1)==0, &
		'each dimension must be a power of 2 in fourn_gather')
	ndim=size(nn)
	ntot=product(nn)
	nprev=1
	allocate(jarr(ntot))
	do idim=1,ndim
		jarr=arth(0,1,ntot)
		n=nn(idim)
		n2=n/2
		msk0=nprev
		msk1=nprev*n2
		msk2=msk0+msk1
		do
			if (msk1 <= msk0) exit
			where (iand(jarr,msk0) == 0 .neqv. iand(jarr,msk1) == 0) &
				jarr=ieor(jarr,msk2)
			msk0=msk0*2
			msk1=msk1/2
			msk2=msk0+msk1
		end do
		data=data(jarr+1)
		allocate(dtemp(ntot),wtab(n2))
		jarr=iand(n-1,arth(0,1,ntot)/nprev)
		m=1
		mm=n2
		mn=m*nprev
		wtab(1)=(1.0_sp,0.0_sp)
		do
			if (mm == 0) exit
			where (iand(jarr,m) /= 0)
				dtemp=data*wtab(mm*iand(jarr,m-1)+1)
				data=eoshift(data,-mn)-dtemp
			elsewhere
				data=data+eoshift(dtemp,mn)
			end where
			m=m*2
			if (m >= n) exit
			mn=m*nprev
			mm=mm/2
			theta=PI_D/(isign*m)
			wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
			wtab(mm+1:n2:2*mm)=wtab(1:n2-mm:2*mm)*wp &
				+wtab(1:n2-mm:2*mm)
		end do
		deallocate(dtemp,wtab)
		nprev=n*nprev
	end do
	deallocate(jarr)
	END SUBROUTINE fourn_gather

! -----------------------------------------------------------

	SUBROUTINE fourrow_sp(data,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,swap
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
	REAL(DP) :: theta
	COMPLEX(SPC), DIMENSION(size(data,1)) :: temp
	COMPLEX(DPC) :: w,wp
	COMPLEX(SPC) :: ws
	n=size(data,2)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_sp')
	n2=n/2
	j=n2
	do i=1,n-2
		if (j > i) call swap(data(:,j+1),data(:,i+1))
		m=n2
		do
			if (m < 2 .or. j < m) exit
			j=j-m
			m=m/2
		end do
		j=j+m
	end do
	mmax=1
	do
		if (n <= mmax) exit
		istep=2*mmax
		theta=PI_D/(isign*mmax)
		wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
		w=cmplx(1.0_dp,0.0_dp,kind=dpc)
		do m=1,mmax
			ws=w
			do i=m,n,istep
				j=i+mmax
				temp=ws*data(:,j)
				data(:,j)=data(:,i)-temp
				data(:,i)=data(:,i)+temp
			end do
			w=w*wp+w
		end do
		mmax=istep
	end do
	END SUBROUTINE fourrow_sp

	SUBROUTINE fourrow_dp(data,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,swap
	IMPLICIT NONE
	COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
	REAL(DP) :: theta
	COMPLEX(DPC), DIMENSION(size(data,1)) :: temp
	COMPLEX(DPC) :: w,wp
	COMPLEX(DPC) :: ws
	n=size(data,2)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_dp')
	n2=n/2
	j=n2
	do i=1,n-2
		if (j > i) call swap(data(:,j+1),data(:,i+1))
		m=n2
		do
			if (m < 2 .or. j < m) exit
			j=j-m
			m=m/2
		end do
		j=j+m
	end do
	mmax=1
	do
		if (n <= mmax) exit
		istep=2*mmax
		theta=PI_D/(isign*mmax)
		wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
		w=cmplx(1.0_dp,0.0_dp,kind=dpc)
		do m=1,mmax
			ws=w
			do i=m,n,istep
				j=i+mmax
				temp=ws*data(:,j)
				data(:,j)=data(:,i)-temp
				data(:,i)=data(:,i)+temp
			end do
			w=w*wp+w
		end do
		mmax=istep
	end do
	END SUBROUTINE fourrow_dp

! -----------------------------------------------------------

	SUBROUTINE fourrow_3d(data,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,swap
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
	REAL(DP) :: theta
	COMPLEX(SPC), DIMENSION(size(data,1),size(data,2)) :: temp
	COMPLEX(DPC) :: w,wp
	COMPLEX(SPC) :: ws
	n=size(data,3)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_3d')
	n2=n/2
	j=n2
	do i=1,n-2
		if (j > i) call swap(data(:,:,j+1),data(:,:,i+1))
		m=n2
		do
			if (m < 2 .or. j < m) exit
			j=j-m
			m=m/2
		end do
		j=j+m
	end do
	mmax=1
	do
		if (n <= mmax) exit
		istep=2*mmax
		theta=PI_D/(isign*mmax)
		wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
		w=cmplx(1.0_dp,0.0_dp,kind=dpc)
		do m=1,mmax
			ws=w
			do i=m,n,istep
				j=i+mmax
				temp=ws*data(:,:,j)
				data(:,:,j)=data(:,:,i)-temp
				data(:,:,i)=data(:,:,i)+temp
			end do
			w=w*wp+w
		end do
		mmax=istep
	end do
	END SUBROUTINE fourrow_3d

! -----------------------------------------------------------

	FUNCTION fpoly(x,n)
	use mo_nrtype
          use mo_nrutil, ONLY : geop
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(n) :: fpoly
	fpoly=geop(1.0_sp,x,n)
	END FUNCTION fpoly

! -----------------------------------------------------------

	SUBROUTINE fred2(a,b,t,f,w,g,ak)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,unit_matrix
	use mo_nr, ONLY : gauleg,lubksb,ludcmp
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(:), INTENT(OUT) :: t,f,w
	INTERFACE
		FUNCTION g(t)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: t
		REAL(SP), DIMENSION(size(t)) :: g
		END FUNCTION g
!BL
		FUNCTION ak(t,s)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
		REAL(SP), DIMENSION(size(t),size(s)) :: ak
		END FUNCTION ak
	END INTERFACE
	INTEGER(I4B) :: n
	INTEGER(I4B), DIMENSION(size(f)) :: indx
	REAL(SP) :: d
	REAL(SP), DIMENSION(size(f),size(f)) :: omk
	n=assert_eq(size(f),size(t),size(w),'fred2')
	call gauleg(a,b,t,w)
	call unit_matrix(omk)
	omk=omk-ak(t,t)*spread(w,dim=1,ncopies=n)
	f=g(t)
	call ludcmp(omk,indx,d)
	call lubksb(omk,indx,f)
	END SUBROUTINE fred2

! -----------------------------------------------------------

	PROGRAM fredex
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	use mo_nr, ONLY : quadmx,ludcmp,lubksb
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=40
	INTEGER(I4B) :: j
	INTEGER(I4B), DIMENSION(N) :: indx
	REAL(SP) :: d
	REAL(SP), DIMENSION(N) :: g,x
	REAL(SP), DIMENSION(N,N) :: a
	call quadmx(a)
	call ludcmp(a,indx,d)
	x(:)=arth(0,1,n)*PI/(n-1)
	g(:)=sin(x(:))
	call lubksb(a,indx,g)
	do j=1,n
		write (*,*) j,x(j),g(j)
	end do
	write (*,*) 'normal completion'
	END PROGRAM fredex

! -----------------------------------------------------------

	FUNCTION fredin(x,a,b,t,f,w,g,ak)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,t,f,w
	REAL(SP), DIMENSION(size(x)) :: fredin
	INTERFACE
		FUNCTION g(t)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: t
		REAL(SP), DIMENSION(size(t)) :: g
		END FUNCTION g
!BL
		FUNCTION ak(t,s)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
		REAL(SP), DIMENSION(size(t),size(s)) :: ak
		END FUNCTION ak
	END INTERFACE
	INTEGER(I4B) :: n
	n=assert_eq(size(f),size(t),size(w),'fredin')
	fredin=g(x)+matmul(ak(x,t),w*f)
	END FUNCTION fredin

! -----------------------------------------------------------

	SUBROUTINE frenel(x,s,c)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), INTENT(OUT) :: s,c
	INTEGER(I4B), PARAMETER :: MAXIT=100
	REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x),BIG=huge(x)*EPS,&
		XMIN=1.5
	INTEGER(I4B) :: k,n
	REAL(SP) :: a,ax,fact,pix2,sign,sum,sumc,sums,term,test
	COMPLEX(SPC) :: b,cc,d,h,del,cs
	LOGICAL(LGT) :: odd
	ax=abs(x)
	if (ax < sqrt(FPMIN)) then
		s=0.0
		c=ax
	else if (ax <= XMIN) then
		sum=0.0
		sums=0.0
		sumc=ax
		sign=1.0
		fact=PIO2*ax*ax
		odd=.true.
		term=ax
		n=3
		do k=1,MAXIT
			term=term*fact/k
			sum=sum+sign*term/n
			test=abs(sum)*EPS
			if (odd) then
				sign=-sign
				sums=sum
				sum=sumc
			else
				sumc=sum
				sum=sums
			end if
			if (term < test) exit
			odd=.not. odd
			n=n+2
		end do
		if (k > MAXIT) call nrerror('frenel: series failed')
		s=sums
		c=sumc
	else
		pix2=PI*ax*ax
		b=cmplx(1.0_sp,-pix2,kind=spc)
		cc=BIG
		d=1.0_sp/b
		h=d
		n=-1
		do k=2,MAXIT
			n=n+2
			a=-n*(n+1)
			b=b+4.0_sp
			d=1.0_sp/(a*d+b)
			cc=b+a/cc
			del=cc*d
			h=h*del
			if (absc(del-1.0_sp) <= EPS) exit
		end do
		if (k > MAXIT) call nrerror('cf failed in frenel')
		h=h*cmplx(ax,-ax,kind=spc)
		cs=cmplx(0.5_sp,0.5_sp,kind=spc)*(1.0_sp-&
			cmplx(cos(0.5_sp*pix2),sin(0.5_sp*pix2),kind=spc)*h)
		c=real(cs)
		s=aimag(cs)
	end if
	if (x < 0.0) then
		c=-c
		s=-s
	end if
	CONTAINS
!BL
	FUNCTION absc(z)
	IMPLICIT NONE
	COMPLEX(SPC), INTENT(IN) :: z
	REAL(SP) :: absc
	absc=abs(real(z))+abs(aimag(z))
	END FUNCTION absc
	END SUBROUTINE frenel

! -----------------------------------------------------------

	SUBROUTINE frprmn(p,ftol,iter,fret)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : linmin
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(SP), INTENT(IN) :: ftol
	REAL(SP), INTENT(OUT) :: fret
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
	INTERFACE
		FUNCTION func(p)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP) :: func
		END FUNCTION func
!BL
		FUNCTION dfunc(p)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP), DIMENSION(size(p)) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=200
	REAL(SP), PARAMETER :: EPS=1.0e-10_sp
	INTEGER(I4B) :: its
	REAL(SP) :: dgg,fp,gam,gg
	REAL(SP), DIMENSION(size(p)) :: g,h,xi
	fp=func(p)
	xi=dfunc(p)
	g=-xi
	h=g
	xi=h
	do its=1,ITMAX
		iter=its
		call linmin(p,xi,fret)
		if (2.0_sp*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) RETURN
		fp=func(p)
		xi=dfunc(p)
		gg=dot_product(g,g)
!		dgg=dot_product(xi,xi)
		dgg=dot_product(xi+g,xi)
		if (gg == 0.0) RETURN
		gam=dgg/gg
		g=-xi
		h=g+gam*h
		xi=h
	end do
	call nrerror('frprmn: maximum iterations exceeded')
	END SUBROUTINE frprmn

! -----------------------------------------------------------

	SUBROUTINE ftest(data1,data2,f,prob)
	use mo_nrtype
	use mo_nr, ONLY : avevar,betai
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: f,prob
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	INTEGER(I4B) :: n1,n2
	REAL(SP) :: ave1,ave2,df1,df2,var1,var2
	n1=size(data1)
	n2=size(data2)
	call avevar(data1,ave1,var1)
	call avevar(data2,ave2,var2)
	if (var1 > var2) then
		f=var1/var2
		df1=n1-1
		df2=n2-1
	else
		f=var2/var1
		df1=n2-1
		df2=n1-1
	end if
	prob=2.0_sp*betai(0.5_sp*df2,0.5_sp*df1,df2/(df2+df1*f))
	if (prob > 1.0) prob=2.0_sp-prob
	END SUBROUTINE ftest

! -----------------------------------------------------------

	FUNCTION gamdev(ia)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : ran1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ia
	REAL(SP) :: gamdev
	REAL(SP) :: am,e,h,s,x,y,v(2),arr(5)
	call assert(ia >= 1, 'gamdev arg')
	if (ia < 6) then
		call ran1(arr(1:ia))
		x=-log(product(arr(1:ia)))
	else
		do
			call ran1(v)
			v(2)=2.0_sp*v(2)-1.0_sp
			if (dot_product(v,v) > 1.0) cycle
			y=v(2)/v(1)
			am=ia-1
			s=sqrt(2.0_sp*am+1.0_sp)
			x=s*y+am
			if (x <= 0.0) cycle
			e=(1.0_sp+y**2)*exp(am*log(x/am)-s*y)
			call ran1(h)
			if (h <= e) exit
		end do
	end if
	gamdev=x
	END FUNCTION gamdev

! -----------------------------------------------------------

	FUNCTION gammln_s(xx)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: xx
	REAL(SP) :: gammln_s
	REAL(DP) :: tmp,x
	REAL(DP) :: stp = 2.5066282746310005_dp
	REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
		-86.50532032941677_dp,24.01409824083091_dp,&
		-1.231739572450155_dp,0.1208650973866179e-2_dp,&
		-0.5395239384953e-5_dp/)
	call assert(xx > 0.0, 'gammln_s arg')
	x=xx
	tmp=x+5.5_dp
	tmp=(x+0.5_dp)*log(tmp)-tmp
	gammln_s=tmp+log(stp*(1.000000000190015_dp+&
		sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
	END FUNCTION gammln_s


	FUNCTION gammln_v(xx)
	use mo_nrtype
          use mo_nrutil, ONLY: assert
	IMPLICIT NONE
	INTEGER(I4B) :: i
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	REAL(SP), DIMENSION(size(xx)) :: gammln_v
	REAL(DP), DIMENSION(size(xx)) :: ser,tmp,x,y
	REAL(DP) :: stp = 2.5066282746310005_dp
	REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
		-86.50532032941677_dp,24.01409824083091_dp,&
		-1.231739572450155_dp,0.1208650973866179e-2_dp,&
		-0.5395239384953e-5_dp/)
	if (size(xx) == 0) RETURN
	call assert(all(xx > 0.0), 'gammln_v arg')
	x=xx
	tmp=x+5.5_dp
	tmp=(x+0.5_dp)*log(tmp)-tmp
	ser=1.000000000190015_dp
	y=x
	do i=1,size(coef)
		y=y+1.0_dp
		ser=ser+coef(i)/y
	end do
	gammln_v=tmp+log(stp*ser/x)
	END FUNCTION gammln_v

! -----------------------------------------------------------

	FUNCTION gammp_s(a,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,x
	REAL(SP) :: gammp_s
	call assert( x >= 0.0,  a > 0.0, 'gammp_s args')
	if (x<a+1.0_sp) then
		gammp_s=gser(a,x)
	else
		gammp_s=1.0_sp-gcf(a,x)
	end if
	END FUNCTION gammp_s


	FUNCTION gammp_v(a,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq
	use mo_nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(SP), DIMENSION(size(x)) :: gammp_v
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(a),size(x),'gammp_v')
	call assert( all(x >= 0.0),  all(a > 0.0), 'gammp_v args')
	mask = (x<a+1.0_sp)
	gammp_v=merge(gser(a,merge(x,0.0_sp,mask)), &
		1.0_sp-gcf(a,merge(x,0.0_sp,.not. mask)),mask)
	END FUNCTION gammp_v

! -----------------------------------------------------------

	FUNCTION gammq_s(a,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,x
	REAL(SP) :: gammq_s
	call assert( x >= 0.0,  a > 0.0, 'gammq_s args')
	if (x<a+1.0_sp) then
		gammq_s=1.0_sp-gser(a,x)
	else
		gammq_s=gcf(a,x)
	end if
	END FUNCTION gammq_s


	FUNCTION gammq_v(a,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq
	use mo_nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(SP), DIMENSION(size(a)) :: gammq_v
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(a),size(x),'gammq_v')
	call assert( all(x >= 0.0),  all(a > 0.0), 'gammq_v args')
	mask = (x<a+1.0_sp)
	gammq_v=merge(1.0_sp-gser(a,merge(x,0.0_sp,mask)), &
		gcf(a,merge(x,0.0_sp,.not. mask)),mask)
	END FUNCTION gammq_v

! -----------------------------------------------------------

	SUBROUTINE gasdev_s(harvest)
	use mo_nrtype
	use mo_nr, ONLY : ran1
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: harvest
	REAL(SP) :: rsq,v1,v2
	REAL(SP), SAVE :: g
	LOGICAL, SAVE :: gaus_stored=.false.
	if (gaus_stored) then
		harvest=g
		gaus_stored=.false.
	else
		do
			call ran1(v1)
			call ran1(v2)
			v1=2.0_sp*v1-1.0_sp
			v2=2.0_sp*v2-1.0_sp
			rsq=v1**2+v2**2
			if (rsq > 0.0 .and. rsq < 1.0) exit
		end do
		rsq=sqrt(-2.0_sp*log(rsq)/rsq)
		harvest=v1*rsq
		g=v2*rsq
		gaus_stored=.true.
	end if
	END SUBROUTINE gasdev_s

	SUBROUTINE gasdev_v(harvest)
	use mo_nrtype
          use mo_nrutil, ONLY : array_copy
	use mo_nr, ONLY : ran1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
	REAL(SP), DIMENSION(size(harvest)) :: rsq,v1,v2
	REAL(SP), ALLOCATABLE, DIMENSION(:), SAVE :: g
	INTEGER(I4B) :: n,ng,nn,m
	INTEGER(I4B), SAVE :: last_allocated=0
	LOGICAL, SAVE :: gaus_stored=.false.
	LOGICAL, DIMENSION(size(harvest)) :: mask
	n=size(harvest)
	if (n /= last_allocated) then
		if (last_allocated /= 0) deallocate(g)
		allocate(g(n))
		last_allocated=n
		gaus_stored=.false.
	end if
	if (gaus_stored) then
		harvest=g
		gaus_stored=.false.
	else
		ng=1
		do
			if (ng > n) exit
			call ran1(v1(ng:n))
			call ran1(v2(ng:n))
			v1(ng:n)=2.0_sp*v1(ng:n)-1.0_sp
			v2(ng:n)=2.0_sp*v2(ng:n)-1.0_sp
			rsq(ng:n)=v1(ng:n)**2+v2(ng:n)**2
			mask(ng:n)=(rsq(ng:n)>0.0 .and. rsq(ng:n)<1.0)
			call array_copy(pack(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
			v2(ng:ng+nn-1)=pack(v2(ng:n),mask(ng:n))
			rsq(ng:ng+nn-1)=pack(rsq(ng:n),mask(ng:n))
			ng=ng+nn
		end do
		rsq=sqrt(-2.0_sp*log(rsq)/rsq)
		harvest=v1*rsq
		g=v2*rsq
		gaus_stored=.true.
	end if
	END SUBROUTINE gasdev_v

! -----------------------------------------------------------

	SUBROUTINE gaucof(a,b,amu0,x,w)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,unit_matrix
	use mo_nr, ONLY : eigsrt,tqli
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
	REAL(SP), INTENT(IN) :: amu0
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(SP), DIMENSION(size(a),size(a)) :: z
	INTEGER(I4B) :: n
	n=assert_eq(size(a),size(b),size(x),size(w),'gaucof')
	b(2:n)=sqrt(b(2:n))
	call unit_matrix(z)
	call tqli(a,b,z)
	call eigsrt(a,z)
	x=a
	w=amu0*z(1,:)**2
	END SUBROUTINE gaucof

! -----------------------------------------------------------

	SUBROUTINE gauher(x,w)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-13_dp,PIM4=0.7511255444649425_dp
	INTEGER(I4B) :: its,j,m,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(SP) :: anu
	REAL(SP), PARAMETER :: C1=9.084064e-01_sp,C2=5.214976e-02_sp,&
		C3=2.579930e-03_sp,C4=3.986126e-03_sp
	REAL(SP), DIMENSION((size(x)+1)/2) :: rhs,r2,r3,theta
	REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION((size(x)+1)/2) :: unfinished
	n=assert_eq(size(x),size(w),'gauher')
	m=(n+1)/2
	anu=2.0_sp*n+1.0_sp
	rhs=arth(3,4,m)*PI/anu
	r3=rhs**(1.0_sp/3.0_sp)
	r2=r3**2
	theta=r3*(C1+r2*(C2+r2*(C3+r2*C4)))
	z=sqrt(anu)*cos(theta)
	unfinished=.true.
	do its=1,MAXIT
		where (unfinished)
			p1=PIM4
			p2=0.0
		end where
		do j=1,n
			where (unfinished)
				p3=p2
				p2=p1
				p1=z*sqrt(2.0_dp/j)*p2-sqrt(real(j-1,dp)/real(j,dp))*p3
			end where
		end do
		where (unfinished)
			pp=sqrt(2.0_dp*n)*p2
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gauher')
	x(1:m)=z
	x(n:n-m+1:-1)=-z
	w(1:m)=2.0_dp/pp**2
	w(n:n-m+1:-1)=w(1:m)
	END SUBROUTINE gauher

! -----------------------------------------------------------

	SUBROUTINE gaujac(x,w,alf,bet)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq,nrerror
	use mo_nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: alf,bet
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-14_dp
	INTEGER(I4B) :: its,j,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(DP) :: alfbet,a,c,temp
	REAL(DP), DIMENSION(size(x)) :: b,p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION(size(x)) :: unfinished
	n=assert_eq(size(x),size(w),'gaujac')
	alfbet=alf+bet
	z=cos(PI*(arth(1,1,n)-0.25_dp+0.5_dp*alf)/(n+0.5_dp*(alfbet+1.0_dp)))
	unfinished=.true.
	do its=1,MAXIT
		temp=2.0_dp+alfbet
		where (unfinished)
			p1=(alf-bet+temp*z)/2.0_dp
			p2=1.0
		end where
		do j=2,n
			a=2*j*(j+alfbet)*temp
			temp=temp+2.0_dp
			c=2.0_dp*(j-1.0_dp+alf)*(j-1.0_dp+bet)*temp
			where (unfinished)
				p3=p2
				p2=p1
				b=(temp-1.0_dp)*(alf*alf-bet*bet+temp*&
					(temp-2.0_dp)*z)
				p1=(b*p2-c*p3)/a
			end where
		end do
		where (unfinished)
			pp=(n*(alf-bet-temp*z)*p1+2.0_dp*(n+alf)*&
				(n+bet)*p2)/(temp*(1.0_dp-z*z))
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gaujac')
	x=z
	w=exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0_sp)-&
		gammln(n+alf+bet+1.0_sp))*temp*2.0_sp**alfbet/(pp*p2)
	END SUBROUTINE gaujac

! -----------------------------------------------------------

	SUBROUTINE gaulag(x,w,alf)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq,nrerror
	use mo_nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: alf
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-13_dp
	INTEGER(I4B) :: its,j,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(SP) :: anu
	REAL(SP), PARAMETER :: C1=9.084064e-01_sp,C2=5.214976e-02_sp,&
		C3=2.579930e-03_sp,C4=3.986126e-03_sp
	REAL(SP), DIMENSION(size(x)) :: rhs,r2,r3,theta
	REAL(DP), DIMENSION(size(x)) :: p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION(size(x)) :: unfinished
	n=assert_eq(size(x),size(w),'gaulag')
	anu=4.0_sp*n+2.0_sp*alf+2.0_sp
	rhs=arth(4*n-1,-4,n)*PI/anu
	r3=rhs**(1.0_sp/3.0_sp)
	r2=r3**2
	theta=r3*(C1+r2*(C2+r2*(C3+r2*C4)))
	z=anu*cos(theta)**2
	unfinished=.true.
	do its=1,MAXIT
		where (unfinished)
			p1=1.0
			p2=0.0
		end where
		do j=1,n
			where (unfinished)
				p3=p2
				p2=p1
				p1=((2.0_dp*j-1.0_dp+alf-z)*p2-(j-1.0_dp+alf)*p3)/j
			end where
		end do
		where (unfinished)
			pp=(n*p1-(n+alf)*p2)/z
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS*z)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gaulag')
	x=z
	w=-exp(gammln(alf+n)-gammln(real(n,sp)))/(pp*n*p2)
	END SUBROUTINE gaulag

! -----------------------------------------------------------

	SUBROUTINE gauleg(x1,x2,x,w)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq,nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-14_dp
	INTEGER(I4B) :: its,j,m,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(DP) :: xl,xm
	REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION((size(x)+1)/2) :: unfinished
	n=assert_eq(size(x),size(w),'gauleg')
	m=(n+1)/2
	xm=0.5_dp*(x2+x1)
	xl=0.5_dp*(x2-x1)
	z=cos(PI_D*(arth(1,1,m)-0.25_dp)/(n+0.5_dp))
	unfinished=.true.
	do its=1,MAXIT
		where (unfinished)
			p1=1.0
			p2=0.0
		end where
		do j=1,n
			where (unfinished)
				p3=p2
				p2=p1
				p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
			end where
		end do
		where (unfinished)
			pp=n*(z*p1-p2)/(z*z-1.0_dp)
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gauleg')
	x(1:m)=xm-xl*z
	x(n:n-m+1:-1)=xm+xl*z
	w(1:m)=2.0_dp*xl/((1.0_dp-z**2)*pp**2)
	w(n:n-m+1:-1)=w(1:m)
	END SUBROUTINE gauleg

! -----------------------------------------------------------

	SUBROUTINE gaussj(a,b)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror,outerand,outerprod,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
	INTEGER(I4B), DIMENSION(size(a,1)) :: ipiv,indxr,indxc
	LOGICAL(LGT), DIMENSION(size(a,1)) :: lpiv
	REAL(SP) :: pivinv
	REAL(SP), DIMENSION(size(a,1)) :: dumc
	INTEGER(I4B), TARGET :: irc(2)
	INTEGER(I4B) :: i,l,n
	INTEGER(I4B), POINTER :: irow,icol
	n=assert_eq(size(a,1),size(a,2),size(b,1),'gaussj')
	irow => irc(1)
	icol => irc(2)
	ipiv=0
	do i=1,n
		lpiv = (ipiv == 0)
		irc=maxloc(abs(a),outerand(lpiv,lpiv))
		ipiv(icol)=ipiv(icol)+1
		if (ipiv(icol) > 1) call nrerror('gaussj: singular matrix (1)')
		if (irow /= icol) then
			call swap(a(irow,:),a(icol,:))
			call swap(b(irow,:),b(icol,:))
		end if
		indxr(i)=irow
		indxc(i)=icol
		if (a(icol,icol) == 0.0) &
			call nrerror('gaussj: singular matrix (2)')
		pivinv=1.0_sp/a(icol,icol)
		a(icol,icol)=1.0
		a(icol,:)=a(icol,:)*pivinv
		b(icol,:)=b(icol,:)*pivinv
		dumc=a(:,icol)
		a(:,icol)=0.0
		a(icol,icol)=pivinv
		a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
		b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
		a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
		b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
	end do
	do l=n,1,-1
		call swap(a(:,indxr(l)),a(:,indxc(l)))
	end do
	END SUBROUTINE gaussj

! -----------------------------------------------------------

	FUNCTION gcf_s(a,x,gln)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,x
	REAL(SP), OPTIONAL, INTENT(OUT) :: gln
	REAL(SP) :: gcf_s
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
	INTEGER(I4B) :: i
	REAL(SP) :: an,b,c,d,del,h
	if (x == 0.0) then
		gcf_s=1.0
		RETURN
	end if
	b=x+1.0_sp-a
	c=1.0_sp/FPMIN
	d=1.0_sp/b
	h=d
	do i=1,ITMAX
		an=-i*(i-a)
		b=b+2.0_sp
		d=an*d+b
		if (abs(d) < FPMIN) d=FPMIN
		c=b+an/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_sp/d
		del=d*c
		h=h*del
		if (abs(del-1.0_sp) <= EPS) exit
	end do
	if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_s')
	if (present(gln)) then
		gln=gammln(a)
		gcf_s=exp(-x+a*log(x)-gln)*h
	else
		gcf_s=exp(-x+a*log(x)-gammln(a))*h
	end if
	END FUNCTION gcf_s


	FUNCTION gcf_v(a,x,gln)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
	REAL(SP), DIMENSION(size(a)) :: gcf_v
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
	INTEGER(I4B) :: i
	REAL(SP), DIMENSION(size(a)) :: an,b,c,d,del,h
	LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
	i=assert_eq(size(a),size(x),'gcf_v')
	zero=(x == 0.0)
	where (zero)
		gcf_v=1.0
	elsewhere
		b=x+1.0_sp-a
		c=1.0_sp/FPMIN
		d=1.0_sp/b
		h=d
	end where
	converged=zero
	do i=1,ITMAX
		where (.not. converged)
			an=-i*(i-a)
			b=b+2.0_sp
			d=an*d+b
			d=merge(FPMIN,d, abs(d)<FPMIN )
			c=b+an/c
			c=merge(FPMIN,c, abs(c)<FPMIN )
			d=1.0_sp/d
			del=d*c
			h=h*del
			converged = (abs(del-1.0_sp)<=EPS)
		end where
		if (all(converged)) exit
	end do
	if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_v')
	if (present(gln)) then
		if (size(gln) < size(a)) call &
			nrerror('gser: Not enough space for gln')
		gln=gammln(a)
		where (.not. zero) gcf_v=exp(-x+a*log(x)-gln)*h
	else
		where (.not. zero) gcf_v=exp(-x+a*log(x)-gammln(a))*h
	end if
	END FUNCTION gcf_v

! -----------------------------------------------------------

	FUNCTION golden(ax,bx,cx,func,tol,xmin)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: ax,bx,cx,tol
	REAL(SP), INTENT(OUT) :: xmin
	REAL(SP) :: golden
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), PARAMETER :: R=0.61803399_sp,C=1.0_sp-R
	REAL(SP) :: f1,f2,x0,x1,x2,x3
	x0=ax
	x3=cx
	if (abs(cx-bx) > abs(bx-ax)) then
		x1=bx
		x2=bx+C*(cx-bx)
	else
		x2=bx
		x1=bx-C*(bx-ax)
	end if
	f1=func(x1)
	f2=func(x2)
	do
		if (abs(x3-x0) <= tol*(abs(x1)+abs(x2))) exit
		if (f2 < f1) then
			call shft3(x0,x1,x2,R*x2+C*x3)
			call shft2(f1,f2,func(x2))
		else
			call shft3(x3,x2,x1,R*x1+C*x0)
			call shft2(f2,f1,func(x1))
		end if
	end do
	if (f1 < f2) then
		golden=f1
		xmin=x1
	else
		golden=f2
		xmin=x2
	end if
	CONTAINS
!BL
	SUBROUTINE shft2(a,b,c)
	REAL(SP), INTENT(OUT) :: a
	REAL(SP), INTENT(INOUT) :: b
	REAL(SP), INTENT(IN) :: c
	a=b
	b=c
	END SUBROUTINE shft2
!BL
	SUBROUTINE shft3(a,b,c,d)
	REAL(SP), INTENT(OUT) :: a
	REAL(SP), INTENT(INOUT) :: b,c
	REAL(SP), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft3
	END FUNCTION golden

! -----------------------------------------------------------

	FUNCTION gser_s(a,x,gln)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,x
	REAL(SP), OPTIONAL, INTENT(OUT) :: gln
	REAL(SP) :: gser_s
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: EPS=epsilon(x)
	INTEGER(I4B) :: n
	REAL(SP) :: ap,del,summ
	if (x == 0.0) then
		gser_s=0.0
		RETURN
	end if
	ap=a
	summ=1.0_sp/a
	del=summ
	do n=1,ITMAX
		ap=ap+1.0_sp
		del=del*x/ap
		summ=summ+del
		if (abs(del) < abs(summ)*EPS) exit
	end do
	if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_s')
	if (present(gln)) then
		gln=gammln(a)
		gser_s=summ*exp(-x+a*log(x)-gln)
	else
		gser_s=summ*exp(-x+a*log(x)-gammln(a))
	end if
	END FUNCTION gser_s


	FUNCTION gser_v(a,x,gln)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : gammln
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
	REAL(SP), DIMENSION(size(a)) :: gser_v
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: EPS=epsilon(x)
	INTEGER(I4B) :: n
	REAL(SP), DIMENSION(size(a)) :: ap,del,summ
	LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
	n=assert_eq(size(a),size(x),'gser_v')
	zero=(x == 0.0)
	where (zero) gser_v=0.0
	ap=a
	summ=1.0_sp/a
	del=summ
	converged=zero
	do n=1,ITMAX
		where (.not. converged)
			ap=ap+1.0_sp
			del=del*x/ap
			summ=summ+del
			converged = (abs(del) < abs(summ)*EPS)
		end where
		if (all(converged)) exit
	end do
	if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_v')
	if (present(gln)) then
		if (size(gln) < size(a)) call &
			nrerror('gser: Not enough space for gln')
		gln=gammln(a)
		where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gln)
	else
		where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gammln(a))
	end if
	END FUNCTION gser_v

! -----------------------------------------------------------

	SUBROUTINE hqr(a,wr,wi)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,diagadd,nrerror,upper_triangle
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: wr,wi
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B) :: i,its,k,l,m,n,nn,mnnk
	REAL(SP) :: anorm,p,q,r,s,t,u,v,w,x,y,z
	REAL(SP), DIMENSION(size(a,1)) :: pp
	n=assert_eq(size(a,1),size(a,2),size(wr),size(wi),'hqr')
	anorm=sum(abs(a),mask=upper_triangle(n,n,extra=2))
	nn=n
	t=0.0
	do
		if (nn < 1) exit
		its=0
		iterate: do
			do l=nn,2,-1
				s=abs(a(l-1,l-1))+abs(a(l,l))
				if (s == 0.0) s=anorm
				if (abs(a(l,l-1))+s == s) exit
			end do
			x=a(nn,nn)
			if (l == nn) then
				wr(nn)=x+t
				wi(nn)=0.0
				nn=nn-1
				exit iterate
			end if
			y=a(nn-1,nn-1)
			w=a(nn,nn-1)*a(nn-1,nn)
			if (l == nn-1) then
				p=0.5_sp*(y-x)
				q=p**2+w
				z=sqrt(abs(q))
				x=x+t
				if (q >= 0.0) then
					z=p+sign(z,p)
					wr(nn)=x+z
					wr(nn-1)=wr(nn)
					if (z /= 0.0) wr(nn)=x-w/z
					wi(nn)=0.0
					wi(nn-1)=0.0
				else
					wr(nn)=x+p
					wr(nn-1)=wr(nn)
					wi(nn)=z
					wi(nn-1)=-z
				end if
				nn=nn-2
				exit iterate
			end if
			if (its == 30) call nrerror('too many iterations in hqr')
			if (its == 10 .or. its == 20) then
				t=t+x
				call diagadd(a(1:nn,1:nn),-x)
				s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
				x=0.75_sp*s
				y=x
				w=-0.4375_sp*s**2
			end if
			its=its+1
			do m=nn-2,l,-1
				z=a(m,m)
				r=x-z
				s=y-z
				p=(r*s-w)/a(m+1,m)+a(m,m+1)
				q=a(m+1,m+1)-z-r-s
				r=a(m+2,m+1)
				s=abs(p)+abs(q)+abs(r)
				p=p/s
				q=q/s
				r=r/s
				if (m == l) exit
				u=abs(a(m,m-1))*(abs(q)+abs(r))
				v=abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1)))
				if (u+v == v) exit
			end do
			do i=m+2,nn
				a(i,i-2)=0.0
				if (i /= m+2) a(i,i-3)=0.0
			end do
			do k=m,nn-1
				if (k /= m) then
					p=a(k,k-1)
					q=a(k+1,k-1)
					r=0.0
					if (k /= nn-1) r=a(k+2,k-1)
					x=abs(p)+abs(q)+abs(r)
					if (x /= 0.0) then
						p=p/x
						q=q/x
						r=r/x
					end if
				end if
				s=sign(sqrt(p**2+q**2+r**2),p)
				if (s /= 0.0) then
					if (k == m) then
						if (l /= m) a(k,k-1)=-a(k,k-1)
					else
						a(k,k-1)=-s*x
					end if
					p=p+s
					x=p/s
					y=q/s
					z=r/s
					q=q/p
					r=r/p
					pp(k:nn)=a(k,k:nn)+q*a(k+1,k:nn)
					if (k /= nn-1) then
						pp(k:nn)=pp(k:nn)+r*a(k+2,k:nn)
						a(k+2,k:nn)=a(k+2,k:nn)-pp(k:nn)*z
					end if
					a(k+1,k:nn)=a(k+1,k:nn)-pp(k:nn)*y
					a(k,k:nn)=a(k,k:nn)-pp(k:nn)*x
					mnnk=min(nn,k+3)
					pp(l:mnnk)=x*a(l:mnnk,k)+y*a(l:mnnk,k+1)
					if (k /= nn-1) then
						pp(l:mnnk)=pp(l:mnnk)+z*a(l:mnnk,k+2)
						a(l:mnnk,k+2)=a(l:mnnk,k+2)-pp(l:mnnk)*r
					end if
					a(l:mnnk,k+1)=a(l:mnnk,k+1)-pp(l:mnnk)*q
					a(l:mnnk,k)=a(l:mnnk,k)-pp(l:mnnk)
				end if
			end do
		end do iterate
	end do
	END SUBROUTINE hqr

! -----------------------------------------------------------

	SUBROUTINE hufdec(ich,code,nb,hcode)
	use mo_nrtype
	USE huf_info
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: ich
	INTEGER(I4B), INTENT(INOUT) :: nb
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: code
	TYPE(huffcode) :: hcode
	INTEGER(I4B) :: l,nc,node
	node=hcode%nodemax
	do
		nc=nb/8+1
		if (nc > size(code)) then
			ich=hcode%nch
			RETURN
		end if
		l=mod(nb,8)
		nb=nb+1
		if (btest(ichar(code(nc)),l)) then
			node=hcode%iright(node)
		else
			node=hcode%left(node)
		end if
		if (node <= hcode%nch) then
			ich=node-1
			RETURN
		end if
	end do
	END SUBROUTINE hufdec

! -----------------------------------------------------------

	SUBROUTINE hufenc(ich,codep,nb,hcode)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror,reallocate
	USE huf_info
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ich
	INTEGER(I4B), INTENT(INOUT) :: nb
	CHARACTER(1), DIMENSION(:), POINTER :: codep
	TYPE(huffcode) :: hcode
	INTEGER(I4B) :: k,l,n,nc,ntmp
	k=ich+1
	if (k > hcode%nch .or. k < 1) call &
		nrerror('hufenc: ich out of range')
	do n=hcode%ncode(k),1,-1
		nc=nb/8+1
		if (nc > size(codep)) codep=>reallocate(codep,2*size(codep))
		l=mod(nb,8)
		if (l == 0) codep(nc)=char(0)
		if (btest(hcode%icode(k),n-1)) then
			ntmp=ibset(ichar(codep(nc)),l)
			codep(nc)=char(ntmp)
		end if
		nb=nb+1
	end do
	END SUBROUTINE hufenc

! -----------------------------------------------------------

MODULE huf_info
	use mo_nrtype
	IMPLICIT NONE
	TYPE huffcode
		INTEGER(I4B) :: nch,nodemax
		INTEGER(I4B), DIMENSION(:), POINTER :: icode,left,iright,ncode
	END TYPE huffcode
CONTAINS
	SUBROUTINE huff_allocate(hcode,mc)
	use mo_nrtype
	IMPLICIT NONE
	TYPE(huffcode) :: hcode
	INTEGER(I4B) :: mc
	INTEGER(I4B) :: mq
	mq=2*mc-1
	allocate(hcode%icode(mq),hcode%ncode(mq),hcode%left(mq),hcode%iright(mq))
	hcode%icode(:)=0
	hcode%ncode(:)=0
	END SUBROUTINE huff_allocate
!BL
	SUBROUTINE huff_deallocate(hcode)
	use mo_nrtype
	IMPLICIT NONE
	TYPE(huffcode) :: hcode
	deallocate(hcode%iright,hcode%left,hcode%ncode,hcode%icode)
	nullify(hcode%icode)
	nullify(hcode%ncode)
	nullify(hcode%left)
	nullify(hcode%iright)
	END SUBROUTINE huff_deallocate
END MODULE huf_info

	SUBROUTINE hufmak(nfreq,ilong,nlong,hcode)
	use mo_nrtype
          use mo_nrutil, ONLY : array_copy,arth,imaxloc,nrerror
	USE huf_info
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: ilong,nlong
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nfreq
	TYPE(huffcode) :: hcode
	INTEGER(I4B) :: ibit,j,k,n,node,nused,nerr
	INTEGER(I4B), DIMENSION(2*size(nfreq)-1) :: indx,iup,nprob
	hcode%nch=size(nfreq)
	call huff_allocate(hcode,size(nfreq))
	nused=0
	nprob(1:hcode%nch)=nfreq(1:hcode%nch)
	call array_copy(pack(arth(1,1,hcode%nch), nfreq(1:hcode%nch) /= 0 ),&
		indx,nused,nerr)
	do j=nused,1,-1
		call hufapp(j)
	end do
	k=hcode%nch
	do
		if (nused <= 1) exit
		node=indx(1)
		indx(1)=indx(nused)
		nused=nused-1
		call hufapp(1)
		k=k+1
		nprob(k)=nprob(indx(1))+nprob(node)
		hcode%left(k)=node
		hcode%iright(k)=indx(1)
		iup(indx(1))=-k
		iup(node)=k
		indx(1)=k
		call hufapp(1)
	end do
	hcode%nodemax=k
	iup(hcode%nodemax)=0
	do j=1,hcode%nch
		if (nprob(j) /= 0) then
			n=0
			ibit=0
			node=iup(j)
			do
				if (node == 0) exit
				if (node < 0) then
					n=ibset(n,ibit)
					node=-node
				end if
				node=iup(node)
				ibit=ibit+1
			end do
			hcode%icode(j)=n
			hcode%ncode(j)=ibit
		end if
	end do
	ilong=imaxloc(hcode%ncode(1:hcode%nch))
	nlong=hcode%ncode(ilong)
	if (nlong > bit_size(1_i4b)) call &
		nrerror('hufmak: Number of possible bits for code exceeded')
	CONTAINS
!BL
	SUBROUTINE hufapp(l)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: l
	INTEGER(I4B) :: i,j,k,n
	n=nused
	i=l
	k=indx(i)
	do
		if (i > n/2) exit
		j=i+i
		if (j < n .and. nprob(indx(j)) > nprob(indx(j+1))) &
			j=j+1
		if (nprob(k) <= nprob(indx(j))) exit
		indx(i)=indx(j)
		i=j
	end do
	indx(i)=k
	END SUBROUTINE hufapp
	END SUBROUTINE hufmak

! -----------------------------------------------------------

	SUBROUTINE hunt(xx,x,jlo)
	use mo_nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: jlo
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	INTEGER(I4B) :: n,inc,jhi,jm
	LOGICAL :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	if (jlo <= 0 .or. jlo > n) then
		jlo=0
		jhi=n+1
	else
		inc=1
		if (x >= xx(jlo) .eqv. ascnd) then
			do
				jhi=jlo+inc
				if (jhi > n) then
					jhi=n+1
					exit
				else
					if (x < xx(jhi) .eqv. ascnd) exit
					jlo=jhi
					inc=inc+inc
				end if
			end do
		else
			jhi=jlo
			do
				jlo=jhi-inc
				if (jlo < 1) then
					jlo=0
					exit
				else
					if (x >= xx(jlo) .eqv. ascnd) exit
					jhi=jlo
					inc=inc+inc
				end if
			end do
		end if
	end if
	do
		if (jhi-jlo <= 1) then
			if (x == xx(n)) jlo=n-1
			if (x == xx(1)) jlo=1
			exit
		else
			jm=(jhi+jlo)/2
			if (x >= xx(jm) .eqv. ascnd) then
				jlo=jm
			else
				jhi=jm
			end if
		end if
	end do
	END SUBROUTINE hunt

! -----------------------------------------------------------

	SUBROUTINE hypdrv(s,ry,rdyds)
	use mo_nrtype
	USE hypgeo_info
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: s
	REAL(SP), DIMENSION(:), INTENT(IN) :: ry
	REAL(SP), DIMENSION(:), INTENT(OUT) :: rdyds
	COMPLEX(SPC), DIMENSION(2) :: y,dyds
	COMPLEX(SPC) :: z
	y=cmplx(ry(1:4:2),ry(2:4:2),kind=spc)
	z=hypgeo_z0+s*hypgeo_dz
	dyds(1)=y(2)*hypgeo_dz
	dyds(2)=((hypgeo_aa*hypgeo_bb)*y(1)-(hypgeo_cc-&
		((hypgeo_aa+hypgeo_bb)+1.0_sp)*z)*y(2))*hypgeo_dz/(z*(1.0_sp-z))
	rdyds(1:4:2)=real(dyds)
	rdyds(2:4:2)=aimag(dyds)
	END SUBROUTINE hypdrv

! -----------------------------------------------------------

MODULE hypgeo_info
	use mo_nrtype
	COMPLEX(SPC) :: hypgeo_aa,hypgeo_bb,hypgeo_cc,hypgeo_dz,hypgeo_z0
END MODULE hypgeo_info


	FUNCTION hypgeo(a,b,c,z)
	use mo_nrtype
	USE hypgeo_info
	use mo_nr, ONLY : bsstep,hypdrv,hypser,odeint
	IMPLICIT NONE
	COMPLEX(SPC), INTENT(IN) :: a,b,c,z
	COMPLEX(SPC) :: hypgeo
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	COMPLEX(SPC), DIMENSION(2) :: y
	REAL(SP), DIMENSION(4) :: ry
	if (real(z)**2+aimag(z)**2 <= 0.25) then
		call hypser(a,b,c,z,hypgeo,y(2))
		RETURN
	else if (real(z) < 0.0) then
		hypgeo_z0=cmplx(-0.5_sp,0.0_sp,kind=spc)
	else if (real(z) <= 1.0) then
		hypgeo_z0=cmplx(0.5_sp,0.0_sp,kind=spc)
	else
		hypgeo_z0=cmplx(0.0_sp,sign(0.5_sp,aimag(z)),kind=spc)
	end if
	hypgeo_aa=a
	hypgeo_bb=b
	hypgeo_cc=c
	hypgeo_dz=z-hypgeo_z0
	call hypser(hypgeo_aa,hypgeo_bb,hypgeo_cc,hypgeo_z0,y(1),y(2))
	ry(1:4:2)=real(y)
	ry(2:4:2)=aimag(y)
	call odeint(ry,0.0_sp,1.0_sp,EPS,0.1_sp,0.0001_sp,hypdrv,bsstep)
	y=cmplx(ry(1:4:2),ry(2:4:2),kind=spc)
	hypgeo=y(1)
	END FUNCTION hypgeo

! -----------------------------------------------------------

	SUBROUTINE hypser(a,b,c,z,series,deriv)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	COMPLEX(SPC), INTENT(IN) :: a,b,c,z
	COMPLEX(SPC), INTENT(OUT) :: series,deriv
	INTEGER(I4B) :: n
	INTEGER(I4B), PARAMETER :: MAXIT=1000
	COMPLEX(SPC) :: aa,bb,cc,fac,temp
	deriv=cmplx(0.0_sp,0.0_sp,kind=spc)
	fac=cmplx(1.0_sp,0.0_sp,kind=spc)
	temp=fac
	aa=a
	bb=b
	cc=c
	do n=1,MAXIT
		fac=((aa*bb)/cc)*fac
		deriv=deriv+fac
		fac=fac*z/n
		series=temp+fac
		if (series == temp) RETURN
		temp=series
		aa=aa+1.0
		bb=bb+1.0
		cc=cc+1.0
	end do
	call nrerror('hypser: convergence failure')
	END SUBROUTINE hypser

! -----------------------------------------------------------

	FUNCTION icrc(crc,buf,jinit,jrev)
	use mo_nrtype
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: buf
	INTEGER(I2B), INTENT(IN) :: crc,jinit
	INTEGER(I4B), INTENT(IN) :: jrev
	INTEGER(I2B) :: icrc
	INTEGER(I4B), SAVE :: init=0
	INTEGER(I2B) :: j,cword,ich
	INTEGER(I2B), DIMENSION(0:255), SAVE :: icrctb,rchr
	INTEGER(I2B), DIMENSION(0:15) :: it = &
		(/ 0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15 /)
	if (init == 0) then
		init=1
		do j=0,255
			icrctb(j)=icrc1(ishft(j,8),char(0))
			rchr(j)=ishft(it(iand(j,15_I2B)),4)+it(ishft(j,-4))
		end do
	end if
	cword=crc
	if (jinit >= 0) then
		cword=ior(jinit,ishft(jinit,8))
	else if (jrev < 0) then
		cword=ior(rchr(hibyte()),ishft(rchr(lobyte()),8))
	end if
	do j=1,size(buf)
		ich=ichar(buf(j))
		if (jrev < 0) ich=rchr(ich)
		cword=ieor(icrctb(ieor(ich,hibyte())),ishft(lobyte(),8))
	end do
	icrc=merge(cword, &
		ior(rchr(hibyte()),ishft(rchr(lobyte()),8)), jrev >= 0)
	CONTAINS
!BL
	FUNCTION hibyte()
	INTEGER(I2B) :: hibyte
	hibyte = ishft(cword,-8)
	END FUNCTION hibyte
!BL
	FUNCTION lobyte()
	INTEGER(I2B) :: lobyte
	lobyte = iand(cword,255_I2B)
	END FUNCTION lobyte
!BL
	FUNCTION icrc1(crc,onech)
	INTEGER(I2B), INTENT(IN) :: crc
	CHARACTER(1), INTENT(IN) :: onech
	INTEGER(I2B) :: icrc1
	INTEGER(I2B) :: i,ich, bit16, ccitt
	DATA bit16,ccitt /Z'8000', Z'1021'/
	ich=ichar(onech)
	icrc1=ieor(crc,ishft(ich,8))
	do i=1,8
		icrc1=merge(ieor(ccitt,ishft(icrc1,1)), &
			ishft(icrc1,1), iand(icrc1,bit16) /= 0)
	end do
	END FUNCTION icrc1
	END FUNCTION icrc

! -----------------------------------------------------------

	FUNCTION igray(n,is)
	use mo_nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n,is
	INTEGER(I4B) :: igray
	INTEGER(I4B) :: idiv,ish
	if (is >= 0) then
		igray=ieor(n,n/2)
	else
		ish=-1
		igray=n
		do
			idiv=ishft(igray,ish)
			igray=ieor(igray,idiv)
			if (idiv <= 1 .or. ish == -16) RETURN
			ish=ish+ish
		end do
	end if
	END FUNCTION igray

! -----------------------------------------------------------

	RECURSIVE SUBROUTINE index_bypack(arr,index,partial)
	use mo_nrtype
          use mo_nrutil, ONLY : array_copy,arth,assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: index
	INTEGER, OPTIONAL, INTENT(IN) :: partial
	REAL(SP) :: a
	INTEGER(I4B) :: n,k,nl,indext,nerr
	INTEGER(I4B), SAVE :: level=0
	LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: mask
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE, SAVE :: temp
	if (present(partial)) then
		n=size(index)
	else
		n=assert_eq(size(index),size(arr),'indexx_bypack')
		index=arth(1,1,n)
	end if
	if (n <= 1) RETURN
	k=(1+n)/2
	call icomp_xchg(index(1),index(k))
	call icomp_xchg(index(k),index(n))
	call icomp_xchg(index(1),index(k))
	if (n <= 3) RETURN
	level=level+1
	if (level == 1) allocate(mask(n),temp(n))
	indext=index(k)
	a=arr(indext)
	mask(1:n) = (arr(index) <= a)
	mask(k) = .false.
	call array_copy(pack(index,mask(1:n)),temp,nl,nerr)
	mask(k) = .true.
	temp(nl+2:n)=pack(index,.not. mask(1:n))
	temp(nl+1)=indext
	index=temp(1:n)
	call index_bypack(arr,index(1:nl),partial=1)
	call index_bypack(arr,index(nl+2:n),partial=1)
	if (level == 1) deallocate(mask,temp)
	level=level-1
	CONTAINS
!BL
	SUBROUTINE icomp_xchg(i,j)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: i,j
	INTEGER(I4B) :: swp
	if (arr(j) < arr(i)) then
		swp=i
		i=j
		j=swp
	end if
	END SUBROUTINE icomp_xchg
	END SUBROUTINE index_bypack

! -----------------------------------------------------------

	SUBROUTINE indexx_sp(arr,index)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq,nrerror,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
	INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
	REAL(SP) :: a
	INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
	INTEGER(I4B), DIMENSION(NSTACK) :: istack
	n=assert_eq(size(index),size(arr),'indexx_sp')
	index=arth(1,1,n)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				indext=index(j)
				a=arr(indext)
				do i=j-1,l,-1
					if (arr(index(i)) <= a) exit
					index(i+1)=index(i)
				end do
				index(i+1)=indext
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(index(k),index(l+1))
			call icomp_xchg(index(l),index(r))
			call icomp_xchg(index(l+1),index(r))
			call icomp_xchg(index(l),index(l+1))
			i=l+1
			j=r
			indext=index(l+1)
			a=arr(indext)
			do
				do
					i=i+1
					if (arr(index(i)) >= a) exit
				end do
				do
					j=j-1
					if (arr(index(j)) <= a) exit
				end do
				if (j < i) exit
				call swap(index(i),index(j))
			end do
			index(l+1)=index(j)
			index(j)=indext
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	CONTAINS
!BL
	SUBROUTINE icomp_xchg(i,j)
	INTEGER(I4B), INTENT(INOUT) :: i,j
	INTEGER(I4B) :: swp
	if (arr(j) < arr(i)) then
		swp=i
		i=j
		j=swp
	end if
	END SUBROUTINE icomp_xchg
	END SUBROUTINE indexx_sp

	SUBROUTINE indexx_i4b(iarr,index)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq,nrerror,swap
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
	INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
	INTEGER(I4B) :: a
	INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
	INTEGER(I4B), DIMENSION(NSTACK) :: istack
	n=assert_eq(size(index),size(iarr),'indexx_sp')
	index=arth(1,1,n)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				indext=index(j)
				a=iarr(indext)
				do i=j-1,1,-1
					if (iarr(index(i)) <= a) exit
					index(i+1)=index(i)
				end do
				index(i+1)=indext
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(index(k),index(l+1))
			call icomp_xchg(index(l),index(r))
			call icomp_xchg(index(l+1),index(r))
			call icomp_xchg(index(l),index(l+1))
			i=l+1
			j=r
			indext=index(l+1)
			a=iarr(indext)
			do
				do
					i=i+1
					if (iarr(index(i)) >= a) exit
				end do
				do
					j=j-1
					if (iarr(index(j)) <= a) exit
				end do
				if (j < i) exit
				call swap(index(i),index(j))
			end do
			index(l+1)=index(j)
			index(j)=indext
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	CONTAINS
!BL
	SUBROUTINE icomp_xchg(i,j)
	INTEGER(I4B), INTENT(INOUT) :: i,j
	INTEGER(I4B) :: swp
	if (iarr(j) < iarr(i)) then
		swp=i
		i=j
		j=swp
	end if
	END SUBROUTINE icomp_xchg
	END SUBROUTINE indexx_i4b

! -----------------------------------------------------------

	FUNCTION interp(uc)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: uc
	REAL(DP), DIMENSION(2*size(uc,1)-1,2*size(uc,1)-1) :: interp
	INTEGER(I4B) :: nc,nf
	nc=assert_eq(size(uc,1),size(uc,2),'interp')
	nf=2*nc-1
	interp(1:nf:2,1:nf:2)=uc(1:nc,1:nc)
	interp(2:nf-1:2,1:nf:2)=0.5_dp*(interp(3:nf:2,1:nf:2)+ &
		interp(1:nf-2:2,1:nf:2))
	interp(1:nf,2:nf-1:2)=0.5_dp*(interp(1:nf,3:nf:2)+interp(1:nf,1:nf-2:2))
	END FUNCTION interp

! -----------------------------------------------------------

	FUNCTION irbit1(iseed)
	use mo_nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: iseed
	INTEGER(I4B) :: irbit1
	if (btest(iseed,17) .neqv. btest(iseed,4) .neqv. btest(iseed,1) &
		.neqv. btest(iseed,0)) then
		iseed=ibset(ishft(iseed,1),0)
		irbit1=1
	else
		iseed=ishft(iseed,1)
		irbit1=0
	end if
	END FUNCTION irbit1

! -----------------------------------------------------------

	FUNCTION irbit2(iseed)
	use mo_nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: iseed
	INTEGER(I4B) :: irbit2
	INTEGER(I4B), PARAMETER :: IB1=1,IB2=2,IB5=16,MASK=IB1+IB2+IB5
	if (btest(iseed,17)) then
		iseed=ibset(ishft(ieor(iseed,MASK),1),0)
		irbit2=1
	else
		iseed=ibclr(ishft(iseed,1),0)
		irbit2=0
	end if
	END FUNCTION irbit2

! -----------------------------------------------------------

	SUBROUTINE jacobi(a,d,v,nrot)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,get_diag,nrerror,unit_matrix,&
		upper_triangle
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: nrot
	REAL(SP), DIMENSION(:), INTENT(OUT) :: d
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
	INTEGER(I4B) :: i,ip,iq,n
	REAL(SP) :: c,g,h,s,sm,t,tau,theta,tresh
	REAL(SP), DIMENSION(size(d)) :: b,z
	n=assert_eq((/size(a,1),size(a,2),size(d),size(v,1),size(v,2)/),'jacobi')
	call unit_matrix(v(:,:))
	b(:)=get_diag(a(:,:))
	d(:)=b(:)
	z(:)=0.0
	nrot=0
	do i=1,50
		sm=sum(abs(a),mask=upper_triangle(n,n))
		if (sm == 0.0) RETURN
		tresh=merge(0.2_sp*sm/n**2,0.0_sp, i < 4 )
		do ip=1,n-1
			do iq=ip+1,n
				g=100.0_sp*abs(a(ip,iq))
				if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
					.and. (abs(d(iq))+g == abs(d(iq)))) then
					a(ip,iq)=0.0
				else if (abs(a(ip,iq)) > tresh) then
					h=d(iq)-d(ip)
					if (abs(h)+g == abs(h)) then
						t=a(ip,iq)/h
					else
						theta=0.5_sp*h/a(ip,iq)
						t=1.0_sp/(abs(theta)+sqrt(1.0_sp+theta**2))
						if (theta < 0.0) t=-t
					end if
					c=1.0_sp/sqrt(1+t**2)
					s=t*c
					tau=s/(1.0_sp+c)
					h=t*a(ip,iq)
					z(ip)=z(ip)-h
					z(iq)=z(iq)+h
					d(ip)=d(ip)-h
					d(iq)=d(iq)+h
					a(ip,iq)=0.0
					call jrotate(a(1:ip-1,ip),a(1:ip-1,iq))
					call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq))
					call jrotate(a(ip,iq+1:n),a(iq,iq+1:n))
					call jrotate(v(:,ip),v(:,iq))
					nrot=nrot+1
				end if
			end do
		end do
		b(:)=b(:)+z(:)
		d(:)=b(:)
		z(:)=0.0
	end do
	call nrerror('too many iterations in jacobi')
	CONTAINS
!BL
	SUBROUTINE jrotate(a1,a2)
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a1,a2
	REAL(SP), DIMENSION(size(a1)) :: wk1
	wk1(:)=a1(:)
	a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
	a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
	END SUBROUTINE jrotate
	END SUBROUTINE jacobi

! -----------------------------------------------------------

	SUBROUTINE jacobn(x,y,dfdx,dfdy)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: dfdx
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dfdy
	dfdx(:)=0.0
	dfdy(1,1)=-0.013_sp-1000.0_sp*y(3)
	dfdy(1,2)=0.0
	dfdy(1,3)=-1000.0_sp*y(1)
	dfdy(2,1)=0.0
	dfdy(2,2)=-2500.0_sp*y(3)
	dfdy(2,3)=-2500.0_sp*y(2)
	dfdy(3,1)=-0.013_sp-1000.0_sp*y(3)
	dfdy(3,2)=-2500.0_sp*y(3)
	dfdy(3,3)=-1000.0_sp*y(1)-2500.0_sp*y(2)
	END SUBROUTINE jacobn

	SUBROUTINE derivs(x,y,dydx)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
	dydx(1)=-0.013_sp*y(1)-1000.0_sp*y(1)*y(3)
	dydx(2)=-2500.0_sp*y(2)*y(3)
	dydx(3)=-0.013_sp*y(1)-1000.0_sp*y(1)*y(3)-2500.0_sp*y(2)*y(3)
	END SUBROUTINE derivs

! -----------------------------------------------------------

	FUNCTION julday(mm,id,iyyy)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: mm,id,iyyy
	INTEGER(I4B) :: julday
	INTEGER(I4B), PARAMETER :: IGREG=15+31*(10+12*1582)
	INTEGER(I4B) :: ja,jm,jy
	jy=iyyy
	if (jy == 0) call nrerror('julday: there is no year zero')
	if (jy < 0) jy=jy+1
	if (mm > 2) then
		jm=mm+1
	else
		jy=jy-1
		jm=mm+13
	end if
	julday=int(365.25_sp*jy)+int(30.6001_sp*jm)+id+1720995
	if (id+31*(mm+12*iyyy) >= IGREG) then
		ja=int(0.01_sp*jy)
		julday=julday+2-ja+int(0.25_sp*ja)
	end if
	END FUNCTION julday

! -----------------------------------------------------------

	SUBROUTINE kendl1(data1,data2,tau,z,prob)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : erfcc
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: tau,z,prob
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	INTEGER(I4B) :: is,j,n,n1,n2
	REAL(SP) :: var
	REAL(SP), DIMENSION(size(data1)) :: a1,a2
	n=assert_eq(size(data1),size(data2),'kendl1')
	n1=0
	n2=0
	is=0
	do j=1,n-1
		a1(j+1:n)=data1(j)-data1(j+1:n)
		a2(j+1:n)=data2(j)-data2(j+1:n)
		n1=n1+count(a1(j+1:n) /= 0.0)
		n2=n2+count(a2(j+1:n) /= 0.0)
		is=is+count((a1(j+1:n) > 0.0 .and. a2(j+1:n) > 0.0) &
			.or. (a1(j+1:n) < 0.0 .and. a2(j+1:n) < 0.0)) - &
			count((a1(j+1:n) > 0.0 .and. a2(j+1:n) < 0.0) &
			.or. (a1(j+1:n) < 0.0 .and. a2(j+1:n) > 0.0))
	end do
	tau=real(is,sp)/sqrt(real(n1,sp)*real(n2,sp))
	var=(4.0_sp*n+10.0_sp)/(9.0_sp*n*(n-1.0_sp))
	z=tau/sqrt(var)
	prob=erfcc(abs(z)/SQRT2)
	END SUBROUTINE kendl1

! -----------------------------------------------------------

	SUBROUTINE kendl2(tab,tau,z,prob)
	use mo_nrtype
          use mo_nrutil, ONLY : cumsum
	use mo_nr, ONLY : erfcc
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: tab
	REAL(SP), INTENT(OUT) :: tau,z,prob
	REAL(SP), DIMENSION(size(tab,1),size(tab,2)) :: cum,cumt
	INTEGER(I4B) :: i,j,ii,jj
	REAL(SP) :: sc,sd,en1,en2,points,var
	ii=size(tab,1)
	jj=size(tab,2)
	do i=1,ii
		cumt(i,jj:1:-1)=cumsum(tab(i,jj:1:-1))
	end do
	en2=sum(tab(1:ii,1:jj-1)*cumt(1:ii,2:jj))
	do j=1,jj
		cum(ii:1:-1,j)=cumsum(cumt(ii:1:-1,j))
	end do
	points=cum(1,1)
	sc=sum(tab(1:ii-1,1:jj-1)*cum(2:ii,2:jj))
	do j=1,jj
		cum(1:ii,j)=cumsum(cumt(1:ii,j))
	end do
	sd=sum(tab(2:ii,1:jj-1)*cum(1:ii-1,2:jj))
	do j=1,jj
		cumt(ii:1:-1,j)=cumsum(tab(ii:1:-1,j))
	end do
	en1=sum(tab(1:ii-1,1:jj)*cumt(2:ii,1:jj))
	tau=(sc-sd)/sqrt((en1+sc+sd)*(en2+sc+sd))
	var=(4.0_sp*points+10.0_sp)/(9.0_sp*points*(points-1.0_sp))
	z=tau/sqrt(var)
	prob=erfcc(abs(z)/SQRT2)
	END SUBROUTINE kendl2

! -----------------------------------------------------------

MODULE kermom_info
	use mo_nrtype
	REAL(DP) :: kermom_x
END MODULE kermom_info

	FUNCTION kermom(y,m)
	use mo_nrtype
	USE kermom_info
	IMPLICIT NONE
	REAL(DP),  INTENT(IN) :: y
	INTEGER(I4B), INTENT(IN) :: m
	REAL(DP), DIMENSION(m) :: kermom
	REAL(DP) :: x,d,df,clog,x2,x3,x4
	x=kermom_x
	if (y >= x) then
		d=y-x
		df=2.0_dp*sqrt(d)*d
		kermom(1:4) = (/ df/3.0_dp, df*(x/3.0_dp+d/5.0_dp),&
			df*((x/3.0_dp + 0.4_dp*d)*x + d**2/7.0_dp),&
			df*(((x/3.0_dp + 0.6_dp*d)*x + 3.0_dp*d**2/7.0_dp)*x&
				+ d**3/9.0_dp) /)
	else
		x2=x**2
		x3=x2*x
		x4=x2*x2
		d=x-y
		clog=log(d)
		kermom(1:4) = (/ d*(clog-1.0_dp),&
			-0.25_dp*(3.0_dp*x+y-2.0_dp*clog*(x+y))*d,&
			(-11.0_dp*x3+y*(6.0_dp*x2+y*(3.0_dp*x+2.0_dp*y))&
				+6.0_dp*clog*(x3-y**3))/18.0_dp,&
			(-25.0_dp*x4+y*(12.0_dp*x3+y*(6.0_dp*x2+y*&
				(4.0_dp*x+3.0_dp*y)))+12.0_dp*clog*(x4-y**4))/48.0_dp /)
	end if
	END FUNCTION kermom

! -----------------------------------------------------------

	SUBROUTINE ks2d1s(x1,y1,quadvl,d1,prob)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : pearsn,probks,quadct
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x1,y1
	REAL(SP), INTENT(OUT) :: d1,prob
	INTERFACE
		SUBROUTINE quadvl(x,y,fa,fb,fc,fd)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
		END SUBROUTINE quadvl
	END INTERFACE
	INTEGER(I4B) :: j,n1
	REAL(SP) :: dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,rr,sqen
	n1=assert_eq(size(x1),size(y1),'ks2d1s')
	d1=0.0
	do j=1,n1
		call quadct(x1(j),y1(j),x1,y1,fa,fb,fc,fd)
		call quadvl(x1(j),y1(j),ga,gb,gc,gd)
		d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
	end do
	call pearsn(x1,y1,r1,dum,dumm)
	sqen=sqrt(real(n1,sp))
	rr=sqrt(1.0_sp-r1**2)
	prob=probks(d1*sqen/(1.0_sp+rr*(0.25_sp-0.75_sp/sqen)))
	END SUBROUTINE ks2d1s

! -----------------------------------------------------------

	SUBROUTINE ks2d2s(x1,y1,x2,y2,d,prob)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : pearsn,probks,quadct
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x1,y1,x2,y2
	REAL(SP), INTENT(OUT) :: d,prob
	INTEGER(I4B) :: j,n1,n2
	REAL(SP) :: d1,d2,dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,r2,rr,sqen
	n1=assert_eq(size(x1),size(y1),'ks2d2s: n1')
	n2=assert_eq(size(x2),size(y2),'ks2d2s: n2')
	d1=0.0
	do j=1,n1
		call quadct(x1(j),y1(j),x1,y1,fa,fb,fc,fd)
		call quadct(x1(j),y1(j),x2,y2,ga,gb,gc,gd)
		d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
	end do
	d2=0.0
	do j=1,n2
		call quadct(x2(j),y2(j),x1,y1,fa,fb,fc,fd)
		call quadct(x2(j),y2(j),x2,y2,ga,gb,gc,gd)
		d2=max(d2,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
	end do
	d=0.5_sp*(d1+d2)
	sqen=sqrt(real(n1,sp)*real(n2,sp)/real(n1+n2,sp))
	call pearsn(x1,y1,r1,dum,dumm)
	call pearsn(x2,y2,r2,dum,dumm)
	rr=sqrt(1.0_sp-0.5_sp*(r1**2+r2**2))
	prob=probks(d*sqen/(1.0_sp+rr*(0.25_sp-0.75_sp/sqen)))
	END SUBROUTINE ks2d2s

! -----------------------------------------------------------

	SUBROUTINE ksone(data,func,d,prob)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	use mo_nr, ONLY : probks,sort
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: d,prob
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B) :: n
	REAL(SP) :: en
	REAL(SP), DIMENSION(size(data)) :: fvals
	REAL(SP), DIMENSION(size(data)+1) :: temp
	call sort(data)
	n=size(data)
	en=n
	fvals(:)=func(data(:))
	temp=arth(0,1,n+1)/en
	d=maxval(max(abs(temp(1:n)-fvals(:)), &
		abs(temp(2:n+1)-fvals(:))))
	en=sqrt(en)
	prob=probks((en+0.12_sp+0.11_sp/en)*d)
	END SUBROUTINE ksone

! -----------------------------------------------------------

	SUBROUTINE kstwo(data1,data2,d,prob)
	use mo_nrtype
          use mo_nrutil, ONLY : cumsum
	use mo_nr, ONLY : probks,sort2
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: d,prob
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	INTEGER(I4B) :: n1,n2
	REAL(SP) :: en1,en2,en
	REAL(SP), DIMENSION(size(data1)+size(data2)) :: dat,org
	n1=size(data1)
	n2=size(data2)
	en1=n1
	en2=n2
	dat(1:n1)=data1
	dat(n1+1:)=data2
	org(1:n1)=0.0
	org(n1+1:)=1.0
	call sort2(dat,org)
	d=maxval(abs(cumsum(org)/en2-cumsum(1.0_sp-org)/en1))
	en=sqrt(en1*en2/(en1+en2))
	prob=probks((en+0.12_sp+0.11_sp/en)*d)
	END SUBROUTINE kstwo

! -----------------------------------------------------------

	SUBROUTINE laguer(a,x,its)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror,poly,poly_term
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: its
	COMPLEX(SPC), INTENT(INOUT) :: x
	COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
	REAL(SP), PARAMETER :: EPS=epsilon(1.0_sp)
	INTEGER(I4B), PARAMETER :: MR=8,MT=10,MAXIT=MT*MR
	INTEGER(I4B) :: iter,m
	REAL(SP) :: abx,abp,abm,err
	COMPLEX(SPC) :: dx,x1,f,g,h,sq,gp,gm,g2
	COMPLEX(SPC), DIMENSION(size(a)) :: b,d
	REAL(SP), DIMENSION(MR) :: frac = &
		(/ 0.5_sp,0.25_sp,0.75_sp,0.13_sp,0.38_sp,0.62_sp,0.88_sp,1.0_sp /)
	m=size(a)-1
	do iter=1,MAXIT
		its=iter
		abx=abs(x)
		b(m+1:1:-1)=poly_term(a(m+1:1:-1),x)
		d(m:1:-1)=poly_term(b(m+1:2:-1),x)
		f=poly(x,d(2:m))
		err=EPS*poly(abx,abs(b(1:m+1)))
		if (abs(b(1)) <= err) RETURN
		g=d(1)/b(1)
		g2=g*g
		h=g2-2.0_sp*f/b(1)
		sq=sqrt((m-1)*(m*h-g2))
		gp=g+sq
		gm=g-sq
		abp=abs(gp)
		abm=abs(gm)
		if (abp < abm) gp=gm
		if (max(abp,abm) > 0.0) then
			dx=m/gp
		else
			dx=exp(cmplx(log(1.0_sp+abx),iter,kind=spc))
		end if
		x1=x-dx
		if (x == x1) RETURN
		if (mod(iter,MT) /= 0) then
			x=x1
		else
			x=x-dx*frac(iter/MT)
		end if
	end do
	call nrerror('laguer: too many iterations')
	END SUBROUTINE laguer

! -----------------------------------------------------------

	SUBROUTINE lfit(x,y,sig,a,maska,covar,chisq,funcs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,diagmult,nrerror
	use mo_nr, ONLY :covsrt,gaussj
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
	REAL(SP), INTENT(OUT) :: chisq
	INTERFACE
		SUBROUTINE funcs(x,arr)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP),INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(OUT) :: arr
		END SUBROUTINE funcs
	END INTERFACE
	INTEGER(I4B) :: i,j,k,l,ma,mfit,n
	REAL(SP) :: sig2i,wt,ym
	REAL(SP), DIMENSION(size(maska)) :: afunc
	REAL(SP), DIMENSION(size(maska),1) :: beta
	n=assert_eq(size(x),size(y),size(sig),'lfit: n')
	ma=assert_eq(size(maska),size(a),size(covar,1),size(covar,2),'lfit: ma')
	mfit=count(maska)
	if (mfit == 0) call nrerror('lfit: no parameters to be fitted')
	covar(1:mfit,1:mfit)=0.0
	beta(1:mfit,1)=0.0
	do i=1,n
		call funcs(x(i),afunc)
		ym=y(i)
		if (mfit < ma) ym=ym-sum(a(1:ma)*afunc(1:ma), mask=.not. maska)
		sig2i=1.0_sp/sig(i)**2
		j=0
		do l=1,ma
			if (maska(l)) then
				j=j+1
				wt=afunc(l)*sig2i
				k=count(maska(1:l))
				covar(j,1:k)=covar(j,1:k)+wt*pack(afunc(1:l),maska(1:l))
				beta(j,1)=beta(j,1)+ym*wt
			end if
		end do
	end do
	call diagmult(covar(1:mfit,1:mfit),0.5_sp)
	covar(1:mfit,1:mfit)= &
		covar(1:mfit,1:mfit)+transpose(covar(1:mfit,1:mfit))
	call gaussj(covar(1:mfit,1:mfit),beta(1:mfit,1:1))
	a(1:ma)=unpack(beta(1:ma,1),maska,a(1:ma))
	chisq=0.0
	do i=1,n
		call funcs(x(i),afunc)
		chisq=chisq+((y(i)-dot_product(a(1:ma),afunc(1:ma)))/sig(i))**2
	end do
	call covsrt(covar,maska)
	END SUBROUTINE lfit

! -----------------------------------------------------------

	SUBROUTINE linbcg(b,x,itol,tol,itmax,iter,err)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : atimes,asolve,snrm
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: b
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
	INTEGER(I4B), INTENT(IN) :: itol,itmax
	REAL(DP), INTENT(IN) :: tol
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(DP), INTENT(OUT) :: err
	REAL(DP), PARAMETER :: EPS=1.0e-14_dp
	INTEGER(I4B) :: n
	REAL(DP) :: ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm
	REAL(DP), DIMENSION(size(b)) :: p,pp,r,rr,z,zz
	n=assert_eq(size(b),size(x),'linbcg')
	iter=0
	call atimes(x,r,0)
	r=b-r
	rr=r
!	call atimes(r,rr,0)
	select case(itol)
		case(1)
			bnrm=snrm(b,itol)
			call asolve(r,z,0)
		case(2)
			call asolve(b,z,0)
			bnrm=snrm(z,itol)
			call asolve(r,z,0)
		case(3:4)
			call asolve(b,z,0)
			bnrm=snrm(z,itol)
			call asolve(r,z,0)
			znrm=snrm(z,itol)
		case default
			call nrerror('illegal itol in linbcg')
	end select
	do
		if (iter > itmax) exit
		iter=iter+1
		call asolve(rr,zz,1)
		bknum=dot_product(z,rr)
		if (iter == 1) then
			p=z
			pp=zz
		else
			bk=bknum/bkden
			p=bk*p+z
			pp=bk*pp+zz
		end if
		bkden=bknum
		call atimes(p,z,0)
		akden=dot_product(z,pp)
		ak=bknum/akden
		call atimes(pp,zz,1)
		x=x+ak*p
		r=r-ak*z
		rr=rr-ak*zz
		call asolve(r,z,0)
		select case(itol)
			case(1)
				err=snrm(r,itol)/bnrm
			case(2)
				err=snrm(z,itol)/bnrm
			case(3:4)
				zm1nrm=znrm
				znrm=snrm(z,itol)
				if (abs(zm1nrm-znrm) > EPS*znrm) then
					dxnrm=abs(ak)*snrm(p,itol)
					err=znrm/abs(zm1nrm-znrm)*dxnrm
				else
					err=znrm/bnrm
					cycle
				end if
				xnrm=snrm(x,itol)
				if (err <= 0.5_dp*xnrm) then
					err=err/xnrm
				else
					err=znrm/bnrm
					cycle
				end if
		end select
		write (*,*) ' iter=',iter,' err=',err
		if (err <= tol) exit
	end do
	END SUBROUTINE linbcg

! -----------------------------------------------------------

MODULE f1dim_mod
	use mo_nrtype
	INTEGER(I4B) :: ncom
	REAL(SP), DIMENSION(:), POINTER :: pcom,xicom
CONTAINS
!BL
	FUNCTION f1dim(x)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: f1dim
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), DIMENSION(:), ALLOCATABLE :: xt
	allocate(xt(ncom))
	xt(:)=pcom(:)+x*xicom(:)
	f1dim=func(xt)
	deallocate(xt)
	END FUNCTION f1dim
END MODULE f1dim_mod

	SUBROUTINE linmin(p,xi,fret)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : mnbrak,brent
	USE f1dim_mod
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: fret
	REAL(SP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
	REAL(SP), PARAMETER :: TOL=1.0e-4_sp
	REAL(SP) :: ax,bx,fa,fb,fx,xmin,xx
	ncom=assert_eq(size(p),size(xi),'linmin')
	pcom=>p
	xicom=>xi
	ax=0.0
	xx=1.0
	call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
	fret=brent(ax,xx,bx,f1dim,TOL,xmin)
	xi=xmin*xi
	p=p+xi
	END SUBROUTINE linmin

! -----------------------------------------------------------

	SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror,vabs
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xold,g
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
	REAL(SP), INTENT(IN) :: fold,stpmax
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x
	REAL(SP), INTENT(OUT) :: f
	LOGICAL(LGT), INTENT(OUT) :: check
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP) :: func
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		END FUNCTION func
	END INTERFACE
	REAL(SP), PARAMETER :: ALF=1.0e-4_sp,TOLX=epsilon(x)
	INTEGER(I4B) :: ndum
	REAL(SP) :: a,alam,alam2,alamin,b,disc,f2,fold2,pabs,rhs1,rhs2,slope,&
		tmplam
	ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
	check=.false.
	pabs=vabs(p(:))
	if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
	slope=dot_product(g,p)
	alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_sp))
	alam=1.0
	do
		x(:)=xold(:)+alam*p(:)
		f=func(x)
		if (alam < alamin) then
			x(:)=xold(:)
			check=.true.
			RETURN
		else if (f <= fold+ALF*alam*slope) then
			RETURN
		else
			if (alam == 1.0) then
				tmplam=-slope/(2.0_sp*(f-fold-slope))
			else
				rhs1=f-fold-alam*slope
				rhs2=f2-fold2-alam2*slope
				a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
				b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
					(alam-alam2)
				if (a == 0.0) then
					tmplam=-slope/(2.0_sp*b)
				else
					disc=b*b-3.0_sp*a*slope
					if (disc < 0.0) call nrerror('roundoff problem in lnsrch')
					tmplam=(-b+sqrt(disc))/(3.0_sp*a)
				end if
				if (tmplam > 0.5_sp*alam) tmplam=0.5_sp*alam
			end if
		end if
		alam2=alam
		f2=f
		fold2=fold
		alam=max(tmplam,0.1_sp*alam)
	end do
	END SUBROUTINE lnsrch

! -----------------------------------------------------------

	FUNCTION locate(xx,x)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	REAL(SP), INTENT(IN) :: x
	INTEGER(I4B) :: locate
	INTEGER(I4B) :: n,jl,jm,ju
	LOGICAL :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	jl=0
	ju=n+1
	do
		if (ju-jl <= 1) exit
		jm=(ju+jl)/2
		if (ascnd .eqv. (x >= xx(jm))) then
			jl=jm
		else
			ju=jm
		end if
	end do
	if (x == xx(1)) then
		locate=1
	else if (x == xx(n)) then
		locate=n-1
	else
		locate=jl
	end if
	END FUNCTION locate

! -----------------------------------------------------------

	FUNCTION lop(u)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: u
	REAL(DP), DIMENSION(size(u,1),size(u,1)) :: lop
	INTEGER(I4B) :: n
	REAL(DP) :: h,h2i
	n=assert_eq(size(u,1),size(u,2),'lop')
	h=1.0_dp/(n-1)
	h2i=1.0_dp/(h*h)
	lop(2:n-1,2:n-1)=h2i*(u(3:n,2:n-1)+u(1:n-2,2:n-1)+u(2:n-1,3:n)+&
		u(2:n-1,1:n-2)-4.0_dp*u(2:n-1,2:n-1))+u(2:n-1,2:n-1)**2
	lop(1:n,1)=0.0
	lop(1:n,n)=0.0
	lop(1,1:n)=0.0
	lop(n,1:n)=0.0
	END FUNCTION lop

! -----------------------------------------------------------

	SUBROUTINE lubksb(a,indx,b)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: i,n,ii,ll
	REAL(SP) :: summ
	n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
	ii=0
	do i=1,n
		ll=indx(i)
		summ=b(ll)
		b(ll)=b(i)
		if (ii /= 0) then
			summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
		else if (summ /= 0.0) then
			ii=i
		end if
		b(i)=summ
	end do
	do i=n,1,-1
		b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
	end do
	END SUBROUTINE lubksb

! -----------------------------------------------------------

	SUBROUTINE ludcmp(a,indx,d)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
	REAL(SP), INTENT(OUT) :: d
	REAL(SP), DIMENSION(size(a,1)) :: vv
	REAL(SP), PARAMETER :: TINY=1.0e-20_sp
	INTEGER(I4B) :: j,n,imax
	n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
	d=1.0
	vv=maxval(abs(a),dim=2)
	if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
	vv=1.0_sp/vv
	do j=1,n
		imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
		if (j /= imax) then
			call swap(a(imax,:),a(j,:))
			d=-d
			vv(imax)=vv(j)
		end if
		indx(j)=imax
		if (a(j,j) == 0.0) a(j,j)=TINY
		a(j+1:n,j)=a(j+1:n,j)/a(j,j)
		a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
	end do
	END SUBROUTINE ludcmp

! -----------------------------------------------------------

	SUBROUTINE machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,&
		maxexp,eps,epsneg,xmin,xmax)
	use mo_nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: ibeta,iexp,irnd,it,machep,maxexp,minexp,negep,ngrd
	REAL(SP), INTENT(OUT) :: eps,epsneg,xmax,xmin
	REAL(SP), PARAMETER :: RX=1.0
	REAL(SP) :: a,beta,betah,one,temp,tempa,two,zero
	ibeta=radix(RX)
	it=digits(RX)
	machep=exponent(nearest(RX,RX)-RX)-1
	negep=exponent(nearest(RX,-RX)-RX)-1
	minexp=minexponent(RX)-1
	maxexp=maxexponent(RX)
	iexp=nint(log(real(maxexp-minexp+2,sp))/log(2.0_sp))
	eps=real(ibeta,sp)**machep
	epsneg=real(ibeta,sp)**negep
	xmax=huge(RX)
	xmin=tiny(RX)
	one=RX
	two=one+one
	zero=one-one
	beta=real(ibeta,sp)
	a=beta**(-negep)
	irnd=0
	betah=beta/two
	temp=a+betah
	if (temp-a /= zero) irnd=1
	tempa=a+beta
	temp=tempa+betah
	if ((irnd == 0) .and. (temp-tempa /= zero)) irnd=2
	ngrd=0
	temp=one+eps
	if ((irnd == 0) .and. (temp*one-one /= zero)) ngrd=1
	temp=xmin/two
	if (temp /= zero) irnd=irnd+3
	END SUBROUTINE machar

! -----------------------------------------------------------

	SUBROUTINE medfit(x,y,a,b,abdev)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : select
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), INTENT(OUT) :: a,b,abdev
	INTEGER(I4B) :: ndata
	REAL(SP) :: aa
	call medfit_private
	CONTAINS
!BL
	SUBROUTINE medfit_private
	IMPLICIT NONE
	REAL(SP) :: b1,b2,bb,chisq,del,f,f1,f2,sigb,sx,sxx,sxy,sy
	REAL(SP), DIMENSION(size(x)) :: tmp
	ndata=assert_eq(size(x),size(y),'medfit')
	sx=sum(x)
	sy=sum(y)
	sxy=dot_product(x,y)
	sxx=dot_product(x,x)
	del=ndata*sxx-sx**2
	aa=(sxx*sy-sx*sxy)/del
	bb=(ndata*sxy-sx*sy)/del
	tmp(:)=y(:)-(aa+bb*x(:))
	chisq=dot_product(tmp,tmp)
	sigb=sqrt(chisq/del)
	b1=bb
	f1=rofunc(b1)
	b2=bb+sign(3.0_sp*sigb,f1)
	f2=rofunc(b2)
	if (b2 == b1) then
		a=aa
		b=bb
		RETURN
	endif
	do
		if (f1*f2 <= 0.0) exit
		bb=b2+1.6_sp*(b2-b1)
		b1=b2
		f1=f2
		b2=bb
		f2=rofunc(b2)
	end do
	sigb=0.01_sp*sigb
	do
		if (abs(b2-b1) <= sigb) exit
		bb=b1+0.5_sp*(b2-b1)
		if (bb == b1 .or. bb == b2) exit
		f=rofunc(bb)
		if (f*f1 >= 0.0) then
			f1=f
			b1=bb
		else
			f2=f
			b2=bb
		end if
	end do
	a=aa
	b=bb
	abdev=abdev/ndata
	END SUBROUTINE medfit_private
!BL
	FUNCTION rofunc(b)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: b
	REAL(SP) :: rofunc
	REAL(SP), PARAMETER :: EPS=epsilon(b)
	INTEGER(I4B) :: j
	REAL(SP), DIMENSION(size(x)) :: arr,d
	arr(:)=y(:)-b*x(:)
	if (mod(ndata,2) == 0) then
		j=ndata/2
		aa=0.5_sp*(select(j,arr)+select(j+1,arr))
	else
		aa=select((ndata+1)/2,arr)
	end if
	d(:)=y(:)-(b*x(:)+aa)
	abdev=sum(abs(d))
	where (y(:) /= 0.0) d(:)=d(:)/abs(y(:))
	rofunc=sum(x(:)*sign(1.0_sp,d(:)), mask=(abs(d(:)) > EPS) )
	END FUNCTION rofunc
	END SUBROUTINE medfit

! -----------------------------------------------------------

	SUBROUTINE memcof(data,xms,d)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: xms
	REAL(SP), DIMENSION(:), INTENT(IN) :: data
	REAL(SP), DIMENSION(:), INTENT(OUT) :: d
	INTEGER(I4B) :: k,m,n
	REAL(SP) :: denom,pneum
	REAL(SP), DIMENSION(size(data)) :: wk1,wk2,wktmp
	REAL(SP), DIMENSION(size(d)) :: wkm
	m=size(d)
	n=size(data)
	xms=dot_product(data,data)/n
	wk1(1:n-1)=data(1:n-1)
	wk2(1:n-1)=data(2:n)
	do k=1,m
		pneum=dot_product(wk1(1:n-k),wk2(1:n-k))
		denom=dot_product(wk1(1:n-k),wk1(1:n-k))+ &
			dot_product(wk2(1:n-k),wk2(1:n-k))
		d(k)=2.0_sp*pneum/denom
		xms=xms*(1.0_sp-d(k)**2)
		d(1:k-1)=wkm(1:k-1)-d(k)*wkm(k-1:1:-1)
		if (k == m) RETURN
		wkm(1:k)=d(1:k)
		wktmp(2:n-k)=wk1(2:n-k)
		wk1(1:n-k-1)=wk1(1:n-k-1)-wkm(k)*wk2(1:n-k-1)
		wk2(1:n-k-1)=wk2(2:n-k)-wkm(k)*wktmp(2:n-k)
	end do
	call nrerror('never get here in memcof')
	END SUBROUTINE memcof

! -----------------------------------------------------------

	SUBROUTINE mgfas(u,maxcyc)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : interp,lop,rstrct,slvsm2
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	INTEGER(I4B), INTENT(IN) :: maxcyc
	INTEGER(I4B) :: j,jcycle,n,ng,ngrid,nn
	REAL(DP) :: res,trerr
	TYPE ptr2d
		REAL(DP), POINTER :: a(:,:)
	END TYPE ptr2d
	TYPE(ptr2d), ALLOCATABLE :: rho(:)
	REAL(DP), DIMENSION(:,:), POINTER :: uj,uj_1
	n=assert_eq(size(u,1),size(u,2),'mgfas')
	ng=nint(log(n-1.0)/log(2.0))
	if (n /= 2**ng+1) call nrerror('n-1 must be a power of 2 in mgfas')
	allocate(rho(ng))
	nn=n
	ngrid=ng
	allocate(rho(ngrid)%a(nn,nn))
	rho(ngrid)%a=u
	do
		if (nn <= 3) exit
		nn=nn/2+1
		ngrid=ngrid-1
		allocate(rho(ngrid)%a(nn,nn))
		rho(ngrid)%a=rstrct(rho(ngrid+1)%a)
	end do
	nn=3
	allocate(uj(nn,nn))
	call slvsm2(uj,rho(1)%a)
	do j=2,ng
		nn=2*nn-1
		uj_1=>uj
		allocate(uj(nn,nn))
		uj=interp(uj_1)
		deallocate(uj_1)
		do jcycle=1,maxcyc
			call mg(j,uj,trerr=trerr)
			res=sqrt(sum((lop(uj)-rho(j)%a)**2))/nn
			if (res < trerr) exit
		end do
	end do
	u=uj
	deallocate(uj)
	do j=1,ng
		deallocate(rho(j)%a)
	end do
	deallocate(rho)
	CONTAINS
!BL
	RECURSIVE SUBROUTINE mg(j,u,rhs,trerr)
	use mo_nrtype
	use mo_nr, ONLY : interp,lop,relax2,rstrct,slvsm2
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: j
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	REAL(DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rhs
	REAL(DP), INTENT(OUT), OPTIONAL :: trerr
	INTEGER(I4B), PARAMETER :: NPRE=1,NPOST=1
	REAL(DP), PARAMETER :: ALPHA=0.33_dp
	INTEGER(I4B) :: jpost,jpre
	REAL(DP), DIMENSION((size(u,1)+1)/2,(size(u,1)+1)/2) :: v,ut,tau
	if (j == 1) then
		call slvsm2(u,rhs+rho(j)%a)
	else
		do jpre=1,NPRE
			if (present(rhs)) then
				call relax2(u,rhs+rho(j)%a)
			else
				call relax2(u,rho(j)%a)
			end if
		end do
		ut=rstrct(u)
		v=ut
		tau=lop(ut)-rstrct(lop(u))
		if (present(rhs)) then
			tau=rstrct(rhs)+tau
		else
			trerr=ALPHA*sqrt(sum(tau**2))/size(tau,1)
		end if
		call mg(j-1,v,tau)
		u=u+interp(v-ut)
		do jpost=1,NPOST
			if (present(rhs)) then
				call relax2(u,rhs+rho(j)%a)
			else
				call relax2(u,rho(j)%a)
			end if
		end do
	end if
	END SUBROUTINE mg
	END SUBROUTINE mgfas

! -----------------------------------------------------------

	SUBROUTINE mglin(u,ncycle)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : interp,rstrct,slvsml
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	INTEGER(I4B), INTENT(IN) :: ncycle
	INTEGER(I4B) :: j,jcycle,n,ng,ngrid,nn
	TYPE ptr2d
		REAL(DP), POINTER :: a(:,:)
	END TYPE ptr2d
	TYPE(ptr2d), ALLOCATABLE :: rho(:)
	REAL(DP), DIMENSION(:,:), POINTER :: uj,uj_1
	n=assert_eq(size(u,1),size(u,2),'mglin')
	ng=nint(log(n-1.0)/log(2.0))
	if (n /= 2**ng+1) call nrerror('n-1 must be a power of 2 in mglin')
	allocate(rho(ng))
	nn=n
	ngrid=ng
	allocate(rho(ngrid)%a(nn,nn))
	rho(ngrid)%a=u
	do
		if (nn <= 3) exit
		nn=nn/2+1
		ngrid=ngrid-1
		allocate(rho(ngrid)%a(nn,nn))
		rho(ngrid)%a=rstrct(rho(ngrid+1)%a)
	end do
	nn=3
	allocate(uj(nn,nn))
	call slvsml(uj,rho(1)%a)
	do j=2,ng
		nn=2*nn-1
		uj_1=>uj
		allocate(uj(nn,nn))
		uj=interp(uj_1)
		deallocate(uj_1)
		do jcycle=1,ncycle
			call mg(j,uj,rho(j)%a)
		end do
	end do
	u=uj
	deallocate(uj)
	do j=1,ng
		deallocate(rho(j)%a)
	end do
	deallocate(rho)
	CONTAINS
!BL
	RECURSIVE SUBROUTINE mg(j,u,rhs)
	use mo_nrtype
	use mo_nr, ONLY : interp,relax,resid,rstrct,slvsml
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: j
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
	INTEGER(I4B), PARAMETER :: NPRE=1,NPOST=1
	INTEGER(I4B) :: jpost,jpre
	REAL(DP), DIMENSION((size(u,1)+1)/2,(size(u,1)+1)/2) :: res,v
	if (j == 1) then
		call slvsml(u,rhs)
	else
		do jpre=1,NPRE
			call relax(u,rhs)
		end do
		res=rstrct(resid(u,rhs))
		v=0.0
		call mg(j-1,v,res)
		u=u+interp(v)
		do jpost=1,NPOST
			call relax(u,rhs)
		end do
	end if
	END SUBROUTINE mg
	END SUBROUTINE mglin

! -----------------------------------------------------------

	SUBROUTINE midexp(funk,aa,bb,s,n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: aa,bb
	REAL(SP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION funk(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: funk
		END FUNCTION funk
	END INTERFACE
	REAL(SP) :: a,b,del
	INTEGER(I4B) :: it
	REAL(SP), DIMENSION(2*3**(n-2)) :: x
	b=exp(-aa)
	a=0.0
	if (n == 1) then
		s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
	else
		it=3**(n-2)
		del=(b-a)/(3.0_sp*it)
		x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
		x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
		s=s/3.0_sp+del*sum(func(x))
	end if
	CONTAINS
!BL
		FUNCTION func(x)
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		func=funk(-log(x))/x
		END FUNCTION func
	END SUBROUTINE midexp

! -----------------------------------------------------------

	SUBROUTINE midinf(funk,aa,bb,s,n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: aa,bb
	REAL(SP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION funk(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: funk
		END FUNCTION funk
	END INTERFACE
	REAL(SP) :: a,b,del
	INTEGER(I4B) :: it
	REAL(SP), DIMENSION(2*3**(n-2)) :: x
	call assert(aa*bb > 0.0, 'midinf args')
	b=1.0_sp/aa
	a=1.0_sp/bb
	if (n == 1) then
		s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
	else
		it=3**(n-2)
		del=(b-a)/(3.0_sp*it)
		x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
		x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
		s=s/3.0_sp+del*sum(func(x))
	end if
	CONTAINS
!BL
		FUNCTION func(x)
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		func=funk(1.0_sp/x)/x**2
		END FUNCTION func
	END SUBROUTINE midinf

! -----------------------------------------------------------

	SUBROUTINE midpnt(func,a,b,s,n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP) :: del
	INTEGER(I4B) :: it
	REAL(SP), DIMENSION(2*3**(n-2)) :: x
	if (n == 1) then
		s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
	else
		it=3**(n-2)
		del=(b-a)/(3.0_sp*it)
		x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
		x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
		s=s/3.0_sp+del*sum(func(x))
	end if
	END SUBROUTINE midpnt

! -----------------------------------------------------------

	SUBROUTINE midsql(funk,aa,bb,s,n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: aa,bb
	REAL(SP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION funk(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: funk
		END FUNCTION funk
	END INTERFACE
	REAL(SP) :: a,b,del
	INTEGER(I4B) :: it
	REAL(SP), DIMENSION(2*3**(n-2)) :: x
	b=sqrt(bb-aa)
	a=0.0
	if (n == 1) then
		s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
	else
		it=3**(n-2)
		del=(b-a)/(3.0_sp*it)
		x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
		x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
		s=s/3.0_sp+del*sum(func(x))
	end if
	CONTAINS
!BL
		FUNCTION func(x)
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		func=2.0_sp*x*funk(aa+x**2)
		END FUNCTION func
	END SUBROUTINE midsql

! -----------------------------------------------------------

	SUBROUTINE midsqu(funk,aa,bb,s,n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: aa,bb
	REAL(SP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION funk(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: funk
		END FUNCTION funk
	END INTERFACE
	REAL(SP) :: a,b,del
	INTEGER(I4B) :: it
	REAL(SP), DIMENSION(2*3**(n-2)) :: x
	b=sqrt(bb-aa)
	a=0.0
	if (n == 1) then
		s=(b-a)*sum(func( (/0.5_sp*(a+b)/) ))
	else
		it=3**(n-2)
		del=(b-a)/(3.0_sp*it)
		x(1:2*it-1:2)=arth(a+0.5_sp*del,3.0_sp*del,it)
		x(2:2*it:2)=x(1:2*it-1:2)+2.0_sp*del
		s=s/3.0_sp+del*sum(func(x))
	end if
	CONTAINS
!BL
		FUNCTION func(x)
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		func=2.0_sp*x*funk(bb-x**2)
		END FUNCTION func
	END SUBROUTINE midsqu

! -----------------------------------------------------------

	RECURSIVE SUBROUTINE miser(func,regn,ndim,npts,dith,ave,var)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP) :: func
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		END FUNCTION func
	END INTERFACE
	REAL(SP), DIMENSION(:), INTENT(IN) :: regn
	INTEGER(I4B), INTENT(IN) :: ndim,npts
	REAL(SP), INTENT(IN) :: dith
	REAL(SP), INTENT(OUT) :: ave,var
	REAL(SP), PARAMETER :: PFAC=0.1_sp,TINY=1.0e-30_sp,BIG=1.0e30_sp
	INTEGER(I4B), PARAMETER :: MNPT=15,MNBS=60
	REAL(SP), DIMENSION(:), ALLOCATABLE :: regn_temp
	INTEGER(I4B) :: j,jb,n,ndum,npre,nptl,nptr
	INTEGER(I4B), SAVE :: iran=0
	REAL(SP) :: avel,varl,fracl,fval,rgl,rgm,rgr,&
		s,sigl,siglb,sigr,sigrb,sm,sm2,sumb,sumr
	REAL(SP), DIMENSION(:), ALLOCATABLE :: fmaxl,fmaxr,fminl,fminr,pt,rmid
	ndum=assert_eq(size(regn),2*ndim,'miser')
	allocate(pt(ndim))
	if (npts < MNBS) then
		sm=0.0
		sm2=0.0
		do n=1,npts
			call ranpt(pt,regn)
			fval=func(pt)
			sm=sm+fval
			sm2=sm2+fval**2
		end do
		ave=sm/npts
		var=max(TINY,(sm2-sm**2/npts)/npts**2)
	else
		npre=max(int(npts*PFAC),MNPT)
		allocate(rmid(ndim),fmaxl(ndim),fmaxr(ndim),fminl(ndim),fminr(ndim))
		fminl(:)=BIG
		fminr(:)=BIG
		fmaxl(:)=-BIG
		fmaxr(:)=-BIG
		do j=1,ndim
			iran=mod(iran*2661+36979,175000)
			s=sign(dith,real(iran-87500,sp))
			rmid(j)=(0.5_sp+s)*regn(j)+(0.5_sp-s)*regn(ndim+j)
		end do
		do n=1,npre
			call ranpt(pt,regn)
			fval=func(pt)
			where (pt <= rmid)
				fminl=min(fminl,fval)
				fmaxl=max(fmaxl,fval)
			elsewhere
				fminr=min(fminr,fval)
				fmaxr=max(fmaxr,fval)
			end where
		end do
		sumb=BIG
		jb=0
		siglb=1.0
		sigrb=1.0
		do j=1,ndim
			if (fmaxl(j) > fminl(j) .and. fmaxr(j) > fminr(j)) then
				sigl=max(TINY,(fmaxl(j)-fminl(j))**(2.0_sp/3.0_sp))
				sigr=max(TINY,(fmaxr(j)-fminr(j))**(2.0_sp/3.0_sp))
				sumr=sigl+sigr
				if (sumr <= sumb) then
					sumb=sumr
					jb=j
					siglb=sigl
					sigrb=sigr
				end if
			end if
		end do
		deallocate(fminr,fminl,fmaxr,fmaxl)
		if (jb == 0) jb=1+(ndim*iran)/175000
		rgl=regn(jb)
		rgm=rmid(jb)
		rgr=regn(ndim+jb)
		fracl=abs((rgm-rgl)/(rgr-rgl))
		nptl=(MNPT+(npts-npre-2*MNPT)*fracl*siglb/ &
			(fracl*siglb+(1.0_sp-fracl)*sigrb))
		nptr=npts-npre-nptl
		allocate(regn_temp(2*ndim))
		regn_temp(:)=regn(:)
		regn_temp(ndim+jb)=rmid(jb)
		call miser(func,regn_temp,ndim,nptl,dith,avel,varl)
		regn_temp(jb)=rmid(jb)
		regn_temp(ndim+jb)=regn(ndim+jb)
		call miser(func,regn_temp,ndim,nptr,dith,ave,var)
		deallocate(regn_temp)
		ave=fracl*avel+(1-fracl)*ave
		var=fracl*fracl*varl+(1-fracl)*(1-fracl)*var
		deallocate(rmid)
	end if
	deallocate(pt)
	CONTAINS
!BL
	SUBROUTINE ranpt(pt,region)
	use mo_nr, ONLY : ran1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: pt
	REAL(SP), DIMENSION(:), INTENT(IN) :: region
	INTEGER(I4B) :: n
	call ran1(pt)
	n=size(pt)
	pt(1:n)=region(1:n)+(region(n+1:2*n)-region(1:n))*pt(1:n)
	END SUBROUTINE ranpt
	END SUBROUTINE miser

! -----------------------------------------------------------

	SUBROUTINE mmid(y,dydx,xs,htot,nstep,yout,derivs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,swap
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: nstep
	REAL(SP), INTENT(IN) :: xs,htot
	REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
	REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: n,ndum
	REAL(SP) :: h,h2,x
	REAL(SP), DIMENSION(size(y)) :: ym,yn
	ndum=assert_eq(size(y),size(dydx),size(yout),'mmid')
	h=htot/nstep
	ym=y
	yn=y+h*dydx
	x=xs+h
	call derivs(x,yn,yout)
	h2=2.0_sp*h
	do n=2,nstep
		call swap(ym,yn)
		yn=yn+h2*yout
		x=x+h
		call derivs(x,yn,yout)
	end do
	yout=0.5_sp*(ym+yn+h*yout)
	END SUBROUTINE mmid

! -----------------------------------------------------------

	SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
	use mo_nrtype
          use mo_nrutil, ONLY : swap
	IMPLICIT NONE
	REAL(SP), INTENT(INOUT) :: ax,bx
	REAL(SP), INTENT(OUT) :: cx,fa,fb,fc
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), PARAMETER :: GOLD=1.618034_sp,GLIMIT=100.0_sp,TINY=1.0e-20_sp
	REAL(SP) :: fu,q,r,u,ulim
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
	REAL(SP), INTENT(OUT) :: a
	REAL(SP), INTENT(INOUT) :: b,c
	REAL(SP), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft
	END SUBROUTINE mnbrak

! -----------------------------------------------------------

	SUBROUTINE mnewt(ntrial,x,tolx,tolf,usrfun)
	use mo_nrtype
	use mo_nr, ONLY : lubksb,ludcmp
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ntrial
	REAL(SP), INTENT(IN) :: tolx,tolf
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
	INTERFACE
		SUBROUTINE usrfun(x,fvec,fjac)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(OUT) :: fvec
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: fjac
		END SUBROUTINE usrfun
	END INTERFACE
	INTEGER(I4B) :: i
	INTEGER(I4B), DIMENSION(size(x)) :: indx
	REAL(SP) :: d
	REAL(SP), DIMENSION(size(x)) :: fvec,p
	REAL(SP), DIMENSION(size(x),size(x)) :: fjac
	do  i=1,ntrial
		call usrfun(x,fvec,fjac)
		if (sum(abs(fvec)) <= tolf) RETURN
		p=-fvec
		call ludcmp(fjac,indx,d)
		call lubksb(fjac,indx,p)
		x=x+p
		if (sum(abs(p)) <= tolx) RETURN
	end do
	END SUBROUTINE mnewt

! -----------------------------------------------------------

	SUBROUTINE moment(data,ave,adev,sdev,var,skew,curt)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: ave,adev,sdev,var,skew,curt
	REAL(SP), DIMENSION(:), INTENT(IN) :: data
	INTEGER(I4B) :: n
	REAL(SP) :: ep
	REAL(SP), DIMENSION(size(data)) :: p,s
	n=size(data)
	if (n <= 1) call nrerror('moment: n must be at least 2')
	ave=sum(data(:))/n
	s(:)=data(:)-ave
	ep=sum(s(:))
	adev=sum(abs(s(:)))/n
	p(:)=s(:)*s(:)
	var=sum(p(:))
	p(:)=p(:)*s(:)
	skew=sum(p(:))
	p(:)=p(:)*s(:)
	curt=sum(p(:))
	var=(var-ep**2/n)/(n-1)
	sdev=sqrt(var)
	if (var /= 0.0) then
		skew=skew/(n*sdev**3)
		curt=curt/(n*var**2)-3.0_sp
	else
		call nrerror('moment: no skew or kurtosis when zero variance')
	end if
	END SUBROUTINE moment

! -----------------------------------------------------------

	SUBROUTINE mp2dfr(a,s,n,m)
	use mo_nrtype
	USE mpops, ONLY : mplsh,mpsmu
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B), INTENT(OUT) :: m
	CHARACTER(1), DIMENSION(:), INTENT(INOUT) :: a
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: s
	INTEGER(I4B), PARAMETER :: IAZ=48
	INTEGER(I4B) :: j
	m=int(2.408_sp*n)
	do j=1,m
		call mpsmu(a,a,n,10)
		s(j)=char(ichar(a(1))+IAZ)
		call mplsh(a,n)
	end do
	END SUBROUTINE mp2dfr

! -----------------------------------------------------------

	SUBROUTINE mpdiv(q,r,u,v,n,m)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : mpinv,mpmul
	USE mpops, ONLY : mpmov,mpsub
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: q,r
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
!	The logical dimensions are: CHARACTER(1) :: q(n-m+1),r(m),u(n),v(m)
	INTEGER(I4B), INTENT(IN) :: n,m
	INTEGER(I4B), PARAMETER :: MACC=3
	INTEGER(I4B) :: is
	CHARACTER(1), DIMENSION(:), ALLOCATABLE :: s
	CHARACTER(1), DIMENSION(:), ALLOCATABLE, TARGET :: rr
	CHARACTER(1), DIMENSION(:), POINTER :: rr2
	allocate(rr(2*(n+MACC)),s(n+MACC))
	rr2=>rr(2:)
	call mpinv(s,v,n-m+MACC,m)
	call mpmul(rr,s,u,n-m+MACC,n)
	call mpmov(q,rr2,n-m+1)
	call mpmul(rr,q,v,n-m+1,m)
	call mpsub(is,rr2,u,rr2,n)
	if (is /= 0) call nrerror('MACC too small in mpdiv')
	call mpmov(r,rr(n-m+2:),m)
	deallocate(rr,s)
	END SUBROUTINE mpdiv

! -----------------------------------------------------------

	SUBROUTINE mpinv(u,v,n,m)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	use mo_nr, ONLY : mpmul
	USE mpops, ONLY : mpmov,mpneg
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: u
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
	INTEGER(I4B), INTENT(IN) :: n,m
	INTEGER(I4B), PARAMETER :: MF=4
	REAL(SP), PARAMETER :: BI=1.0_sp/256.0_sp
	INTEGER(I4B) :: i,j,mm
	REAL(SP) :: fu
	CHARACTER(1), DIMENSION(:), ALLOCATABLE :: rr,s
	allocate(rr(n+m+1),s(n))
	mm=min(MF,m)
	fu=1.0_sp/poly(BI,real(ichar(v(:)),sp))
	do j=1,n
		i=int(fu)
		u(j)=char(i)
		fu=256.0_sp*(fu-i)
	end do
	do
		call mpmul(rr,u,v,n,m)
		call mpmov(s,rr(2:),n)
		call mpneg(s,n)
		s(1)=char(ichar(s(1))-254)
		call mpmul(rr,s,u,n,n)
		call mpmov(u,rr(2:),n)
		if (all(ichar(s(2:n-1)) == 0)) exit
	end do
	deallocate(rr,s)
	END SUBROUTINE mpinv

! -----------------------------------------------------------

	SUBROUTINE mpmul(w,u,v,n,m)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : realft
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n,m
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
!	The logical dimensions are: CHARACTER(1) :: w(n+m),u(n),v(m)
	REAL(DP), PARAMETER :: RX=256.0
	INTEGER(I4B) :: j,mn,nn
	REAL(DP) :: cy,t
	REAL(DP), DIMENSION(:), ALLOCATABLE :: a,b,tb
	mn=max(m,n)
	nn=1
	do
		if (nn >= mn) exit
		nn=nn+nn
	end do
	nn=nn+nn
	allocate(a(nn),b(nn),tb((nn-1)/2))
	a(1:n)=ichar(u(1:n))
	a(n+1:nn)=0.0
	b(1:m)=ichar(v(1:m))
	b(m+1:nn)=0.0
	call realft(a(1:nn),1)
	call realft(b(1:nn),1)
	b(1)=b(1)*a(1)
	b(2)=b(2)*a(2)
	tb=b(3:nn:2)
	b(3:nn:2)=tb*a(3:nn:2)-b(4:nn:2)*a(4:nn:2)
	b(4:nn:2)=tb*a(4:nn:2)+b(4:nn:2)*a(3:nn:2)
	call realft(b(1:nn),-1)
	b(:)=b(:)/(nn/2)
	cy=0.0
	do j=nn,1,-1
		t=b(j)+cy+0.5_dp
		b(j)=mod(t,RX)
		cy=int(t/RX)
	end do
	if (cy >= RX) call nrerror('mpmul: sanity check failed in fftmul')
	w(1)=char(int(cy))
	w(2:(n+m))=char(int(b(1:(n+m-1))))
	deallocate(a,b,tb)
	END SUBROUTINE mpmul

! -----------------------------------------------------------

MODULE mpops
	use mo_nrtype
	INTEGER(I4B), PARAMETER :: NPAR_ICARRY=64
	CONTAINS
!BL
	SUBROUTINE icarry(karry,isum,nbits)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: karry
	INTEGER(I2B), DIMENSION(:), INTENT(INOUT) :: isum
	INTEGER(I4B), INTENT(IN) :: nbits
	INTEGER(I4B) :: n,j
	INTEGER(I2B), DIMENSION(size(isum)) :: ihi
	INTEGER(I2B) :: mb,ihh
	n=size(isum)
	mb=ishft(1,nbits)-1
	karry=0
	if (n < NPAR_ICARRY ) then
		do j=n,2,-1
			ihh=ishft(isum(j),-nbits)
			if (ihh /= 0) then
				isum(j)=iand(isum(j),mb)
				isum(j-1)=isum(j-1)+ihh
			end if
		end do
		ihh=ishft(isum(1),-nbits)
		isum(1)=iand(isum(1),mb)
		karry=karry+ihh
	else
		do
			ihi=ishft(isum,-nbits)
			if (all(ihi == 0)) exit
			where (ihi /= 0) isum=iand(isum,mb)
			where (ihi(2:n) /= 0) isum(1:n-1)=isum(1:n-1)+ihi(2:n)
			karry=karry+ihi(1)
		end do
	end if
	END SUBROUTINE icarry
!BL
	SUBROUTINE mpadd(w,u,v,n)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I2B), DIMENSION(n) :: isum
	INTEGER(I4B) :: karry
	isum=ichar(u(1:n))+ichar(v(1:n))
	call icarry(karry,isum,8_I4B)
	w(2:n+1)=char(isum)
	w(1)=char(karry)
	END SUBROUTINE mpadd
!BL
	SUBROUTINE mpsub(is,w,u,v,n)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: is
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u,v
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B) :: karry
	INTEGER(I2B), DIMENSION(n) :: isum
	isum=255+ichar(u(1:n))-ichar(v(1:n))
	isum(n)=isum(n)+1
	call icarry(karry,isum,8_I4B)
	w(1:n)=char(isum)
	is=karry-1
	END SUBROUTINE mpsub
!BL
	SUBROUTINE mpsad(w,u,n,iv)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u
	INTEGER(I4B), INTENT(IN) :: n,iv
	INTEGER(I4B) :: karry
	INTEGER(I2B), DIMENSION(n) :: isum
	isum=ichar(u(1:n))
	isum(n)=isum(n)+iv
	call icarry(karry,isum,8_I4B)
	w(2:n+1)=char(isum)
	w(1)=char(karry)
	END SUBROUTINE mpsad
!BL
	SUBROUTINE mpsmu(w,u,n,iv)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u
	INTEGER(I4B), INTENT(IN) :: n,iv
	INTEGER(I4B) :: karry
	INTEGER(I2B), DIMENSION(n) :: isum
	isum=ichar(u(1:n))*iv
	call icarry(karry,isum,8_I4B)
	w(2:n+1)=char(isum)
	w(1)=char(karry)
	END SUBROUTINE mpsmu
!BL
	SUBROUTINE mpneg(u,n)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(INOUT) :: u
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B) :: karry
	INTEGER(I2B), DIMENSION(n) :: isum
	isum=255-ichar(u(1:n))
	isum(n)=isum(n)+1
	call icarry(karry,isum,8_I4B)
	u(1:n)=char(isum)
	END SUBROUTINE mpneg
!BL
	SUBROUTINE mplsh(u,n)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(INOUT) :: u
	INTEGER(I4B), INTENT(IN) :: n
	u(1:n)=u(2:n+1)
	END SUBROUTINE mplsh
!BL
	SUBROUTINE mpmov(u,v,n)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: u
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
	INTEGER(I4B), INTENT(IN) :: n
	u(1:n)=v(1:n)
	END SUBROUTINE mpmov
!BL
	SUBROUTINE mpsdv(w,u,n,iv,ir)
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: u
	INTEGER(I4B), INTENT(IN) :: n,iv
	INTEGER(I4B), INTENT(OUT) :: ir
	INTEGER(I4B) :: i,j
	ir=0
	do j=1,n
		i=256*ir+ichar(u(j))
		w(j)=char(i/iv)
		ir=mod(i,iv)
	end do
	END SUBROUTINE mpsdv
END MODULE mpops

! -----------------------------------------------------------

	SUBROUTINE mppi(n)
	use mo_nrtype
	use mo_nr, ONLY : mp2dfr,mpinv,mpmul,mpsqrt
	USE mpops, ONLY : mpadd,mplsh,mpmov,mpsdv
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B), PARAMETER :: IAOFF=48
	INTEGER(I4B) :: ir,j,m
	CHARACTER(1), DIMENSION(n) :: sx,sxi
	CHARACTER(1), DIMENSION(2*n) :: t,y
	CHARACTER(1), DIMENSION(3*n) :: s
	CHARACTER(1), DIMENSION(n+1) :: x,bigpi
	t(1)=char(2)
	t(2:n)=char(0)
	call mpsqrt(x,x,t,n,n)
	call mpadd(bigpi,t,x,n)
	call mplsh(bigpi,n)
	call mpsqrt(sx,sxi,x,n,n)
	call mpmov(y,sx,n)
	do
		call mpadd(x,sx,sxi,n)
		call mpsdv(x,x(2:),n,2,ir)
		call mpsqrt(sx,sxi,x,n,n)
		call mpmul(t,y,sx,n,n)
		call mpadd(t(2:),t(2:),sxi,n)
		x(1)=char(ichar(x(1))+1)
		y(1)=char(ichar(y(1))+1)
		call mpinv(s,y,n,n)
		call mpmul(y,t(3:),s,n,n)
		call mplsh(y,n)
		call mpmul(t,x,s,n,n)
		m=mod(255+ichar(t(2)),256)
		if (abs(ichar(t(n+1))-m) > 1 .or. any(ichar(t(3:n)) /= m)) then
			call mpmul(s,bigpi,t(2:),n,n)
			call mpmov(bigpi,s(2:),n)
			cycle
		end if
		write (*,*) 'pi='
		s(1)=char(ichar(bigpi(1))+IAOFF)
		s(2)='.'
		call mp2dfr(bigpi(2:),s(3:),n-1,m)
		write (*,'(1x,64a1)') (s(j),j=1,m+1)
		RETURN
	end do
	END SUBROUTINE mppi

! -----------------------------------------------------------

	SUBROUTINE mprove(a,alud,indx,b,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : lubksb
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a,alud
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
	REAL(SP), DIMENSION(:), INTENT(IN) :: b
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
	INTEGER(I4B) :: ndum
	REAL(SP), DIMENSION(size(a,1)) :: r
	ndum=assert_eq((/size(a,1),size(a,2),size(alud,1),size(alud,2),size(b),&
		size(x),size(indx)/),'mprove')
	r=matmul(real(a,dp),real(x,dp))-real(b,dp)
	call lubksb(alud,indx,r)
	x=x-r
	END SUBROUTINE mprove

! -----------------------------------------------------------

	SUBROUTINE mpsqrt(w,u,v,n,m)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	use mo_nr, ONLY : mpmul
	USE mpops, ONLY : mplsh,mpmov,mpneg,mpsdv
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(OUT) :: w,u
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: v
	INTEGER(I4B), INTENT(IN) :: n,m
	INTEGER(I4B), PARAMETER :: MF=3
	REAL(SP), PARAMETER :: BI=1.0_sp/256.0_sp
	INTEGER(I4B) :: i,ir,j,mm
	REAL(SP) :: fu
	CHARACTER(1), DIMENSION(:), ALLOCATABLE :: r,s
	allocate(r(2*n),s(2*n))
	mm=min(m,MF)
	fu=1.0_sp/sqrt(poly(BI,real(ichar(v(:)),sp)))
	do j=1,n
		i=int(fu)
		u(j)=char(i)
		fu=256.0_sp*(fu-i)
	end do
	do
		call mpmul(r,u,u,n,n)
		call mplsh(r,n)
		call mpmul(s,r,v,n,m)
		call mplsh(s,n)
		call mpneg(s,n)
		s(1)=char(ichar(s(1))-253)
		call mpsdv(s,s,n,2,ir)
		if (any(ichar(s(2:n-1)) /= 0)) then
			call mpmul(r,s,u,n,n)
			call mpmov(u,r(2:),n)
			cycle
		end if
		call mpmul(r,u,v,n,m)
		call mpmov(w,r(2:),n)
		deallocate(r,s)
		RETURN
	end do
	END SUBROUTINE mpsqrt

! -----------------------------------------------------------

	SUBROUTINE mrqmin(x,y,sig,a,maska,covar,alpha,chisq,funcs,alamda)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,diagmult
	use mo_nr, ONLY : covsrt,gaussj
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: covar,alpha
	REAL(SP), INTENT(OUT) :: chisq
	REAL(SP), INTENT(INOUT) :: alamda
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
	INTERFACE
		SUBROUTINE funcs(x,a,yfit,dyda)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yfit
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
		END SUBROUTINE funcs
	END INTERFACE
	INTEGER(I4B) :: ma,ndata
	INTEGER(I4B), SAVE :: mfit
	call mrqmin_private
	CONTAINS
!BL
	SUBROUTINE mrqmin_private
	REAL(SP), SAVE :: ochisq
	REAL(SP), DIMENSION(:), ALLOCATABLE, SAVE :: atry,beta
	REAL(SP), DIMENSION(:,:), ALLOCATABLE, SAVE :: da
	ndata=assert_eq(size(x),size(y),size(sig),'mrqmin: ndata')
	ma=assert_eq((/size(a),size(maska),size(covar,1),size(covar,2),&
		size(alpha,1),size(alpha,2)/),'mrqmin: ma')
	mfit=count(maska)
	if (alamda < 0.0) then
		allocate(atry(ma),beta(ma),da(ma,1))
		alamda=0.001_sp
		call mrqcof(a,alpha,beta)
		ochisq=chisq
		atry=a
	end if
	covar(1:mfit,1:mfit)=alpha(1:mfit,1:mfit)
	call diagmult(covar(1:mfit,1:mfit),1.0_sp+alamda)
	da(1:mfit,1)=beta(1:mfit)
	call gaussj(covar(1:mfit,1:mfit),da(1:mfit,1:1))
	if (alamda == 0.0) then
		call covsrt(covar,maska)
		deallocate(atry,beta,da)
		RETURN
	end if
	atry=a+unpack(da(1:mfit,1),maska,0.0_sp)
	call mrqcof(atry,covar,da(1:mfit,1))
	if (chisq < ochisq) then
		alamda=0.1_sp*alamda
		ochisq=chisq
		alpha(1:mfit,1:mfit)=covar(1:mfit,1:mfit)
		beta(1:mfit)=da(1:mfit,1)
		a=atry
	else
		alamda=10.0_sp*alamda
		chisq=ochisq
	end if
	END SUBROUTINE mrqmin_private
!BL
	SUBROUTINE mrqcof(a,alpha,beta)
	REAL(SP), DIMENSION(:), INTENT(IN) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: beta
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: alpha
	INTEGER(I4B) :: j,k,l,m
	REAL(SP), DIMENSION(size(x),size(a)) :: dyda
	REAL(SP), DIMENSION(size(x)) :: dy,sig2i,wt,ymod
	call funcs(x,a,ymod,dyda)
	sig2i=1.0_sp/(sig**2)
	dy=y-ymod
	j=0
	do l=1,ma
		if (maska(l)) then
			j=j+1
			wt=dyda(:,l)*sig2i
			k=0
			do m=1,l
				if (maska(m)) then
					k=k+1
					alpha(j,k)=dot_product(wt,dyda(:,m))
					alpha(k,j)=alpha(j,k)
				end if
			end do
			beta(j)=dot_product(dy,wt)
		end if
	end do
	chisq=dot_product(dy**2,sig2i)
	END SUBROUTINE mrqcof
	END SUBROUTINE mrqmin

! -----------------------------------------------------------

	SUBROUTINE newt(x,check)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror,vabs
	use mo_nr, ONLY : fdjac,lnsrch,lubksb,ludcmp
	USE fminln
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
	LOGICAL(LGT), INTENT(OUT) :: check
	INTEGER(I4B), PARAMETER :: MAXITS=200
	REAL(SP), PARAMETER :: TOLF=1.0e-4_sp,TOLMIN=1.0e-6_sp,TOLX=epsilon(x),&
		STPMX=100.0
	INTEGER(I4B) :: its
	INTEGER(I4B), DIMENSION(size(x)) :: indx
	REAL(SP) :: d,f,fold,stpmax
	REAL(SP), DIMENSION(size(x)) :: g,p,xold
	REAL(SP), DIMENSION(size(x)), TARGET :: fvec
	REAL(SP), DIMENSION(size(x),size(x)) :: fjac
	fmin_fvecp=>fvec
	f=fmin(x)
	if (maxval(abs(fvec(:))) < 0.01_sp*TOLF) then
		check=.false.
		RETURN
	end if
	stpmax=STPMX*max(vabs(x(:)),real(size(x),sp))
	do its=1,MAXITS
		call fdjac(x,fvec,fjac)
		g(:)=matmul(fvec(:),fjac(:,:))
		xold(:)=x(:)
		fold=f
		p(:)=-fvec(:)
		call ludcmp(fjac,indx,d)
		call lubksb(fjac,indx,p)
		call lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin)
		if (maxval(abs(fvec(:))) < TOLF) then
			check=.false.
			RETURN
		end if
		if (check) then
			check=(maxval(abs(g(:))*max(abs(x(:)),1.0_sp) / &
				max(f,0.5_sp*size(x))) < TOLMIN)
			RETURN
		end if
		if (maxval(abs(x(:)-xold(:))/max(abs(x(:)),1.0_sp)) < TOLX) &
			RETURN
	end do
	call nrerror('MAXITS exceeded in newt')
	END SUBROUTINE newt

! -----------------------------------------------------------

MODULE ode_path
	use mo_nrtype
	INTEGER(I4B) :: nok,nbad,kount
	LOGICAL(LGT), SAVE :: save_steps=.false.
	REAL(SP) :: dxsav
	REAL(SP), DIMENSION(:), POINTER :: xp
	REAL(SP), DIMENSION(:,:), POINTER :: yp
END MODULE ode_path

	SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror,reallocate
	USE ode_path
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: ystart
	REAL(SP), INTENT(IN) :: x1,x2,eps,h1,hmin
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
!BL
		SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
		REAL(SP), INTENT(INOUT) :: x
		REAL(SP), INTENT(IN) :: htry,eps
		REAL(SP), INTENT(OUT) :: hdid,hnext
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			use mo_nrtype
			IMPLICIT NONE
			REAL(SP), INTENT(IN) :: x
			REAL(SP), DIMENSION(:), INTENT(IN) :: y
			REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
		END SUBROUTINE rkqs
	END INTERFACE
	REAL(SP), PARAMETER :: TINY=1.0e-30_sp
	INTEGER(I4B), PARAMETER :: MAXSTP=10000
	INTEGER(I4B) :: nstp
	REAL(SP) :: h,hdid,hnext,x,xsav
	REAL(SP), DIMENSION(size(ystart)) :: dydx,y,yscal
	x=x1
	h=sign(h1,x2-x1)
	nok=0
	nbad=0
	kount=0
	y(:)=ystart(:)
	if (save_steps) then
		xsav=x-2.0_sp*dxsav
		nullify(xp,yp)
		allocate(xp(256))
		allocate(yp(size(ystart),size(xp)))
	end if
	do nstp=1,MAXSTP
		call derivs(x,y,dydx)
		yscal(:)=abs(y(:))+abs(h*dydx(:))+TINY
		if (save_steps .and. (abs(x-xsav) > abs(dxsav))) &
			call save_a_step
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x
		call rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)
		if (hdid == h) then
			nok=nok+1
		else
			nbad=nbad+1
		end if
		if ((x-x2)*(x2-x1) >= 0.0) then
			ystart(:)=y(:)
			if (save_steps) call save_a_step
			RETURN
		end if
		if (abs(hnext) < hmin)&
			call nrerror('stepsize smaller than minimum in odeint')
		h=hnext
	end do
	call nrerror('too many steps in odeint')
	CONTAINS
!BL
	SUBROUTINE save_a_step
	kount=kount+1
	if (kount > size(xp)) then
		xp=>reallocate(xp,2*size(xp))
		yp=>reallocate(yp,size(yp,1),size(xp))
	end if
	xp(kount)=x
	yp(:,kount)=y(:)
	xsav=x
	END SUBROUTINE save_a_step
	END SUBROUTINE odeint

! -----------------------------------------------------------

	SUBROUTINE orthog(anu,alpha,beta,a,b)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: anu,alpha,beta
	REAL(SP), DIMENSION(:), INTENT(OUT) :: a,b
	INTEGER(I4B) :: k,n,ndum
	REAL(SP), DIMENSION(2*size(a)+1,2*size(a)+1) :: sig
	n=assert_eq(size(a),size(b),'orthog: n')
	ndum=assert_eq(2*n,size(alpha)+1,size(anu),size(beta)+1,'orthog: ndum')
	sig(1,3:2*n)=0.0
	sig(2,2:2*n+1)=anu(1:2*n)
	a(1)=alpha(1)+anu(2)/anu(1)
	b(1)=0.0
	do k=3,n+1
		sig(k,k:2*n-k+3)=sig(k-1,k+1:2*n-k+4)+(alpha(k-1:2*n-k+2) &
			-a(k-2))*sig(k-1,k:2*n-k+3)-b(k-2)*sig(k-2,k:2*n-k+3) &
			+beta(k-1:2*n-k+2)*sig(k-1,k-1:2*n-k+2)
		a(k-1)=alpha(k-1)+sig(k,k+1)/sig(k,k)-sig(k-1,k)/sig(k-1,k-1)
		b(k-1)=sig(k,k)/sig(k-1,k-1)
	end do
	END SUBROUTINE orthog

! -----------------------------------------------------------

	SUBROUTINE pade(cof,resid)
	use mo_nrtype
	use mo_nr, ONLY : lubksb,ludcmp,mprove
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: cof
	REAL(SP), INTENT(OUT) :: resid
	INTEGER(I4B) :: k,n
	INTEGER(I4B), DIMENSION((size(cof)-1)/2) :: indx
	REAL(SP), PARAMETER :: BIG=1.0e30_sp
	REAL(SP) :: d,rr,rrold
	REAL(SP), DIMENSION((size(cof)-1)/2) :: x,y,z
	REAL(SP), DIMENSION((size(cof)-1)/2,(size(cof)-1)/2) :: q,qlu
	n=(size(cof)-1)/2
	x=cof(n+2:2*n+1)
	y=x
	do k=1,n
		q(:,k)=cof(n+2-k:2*n+1-k)
	end do
	qlu=q
	call ludcmp(qlu,indx,d)
	call lubksb(qlu,indx,x)
	rr=BIG
	do
		rrold=rr
		z=x
		call mprove(q,qlu,indx,y,x)
		rr=sum((z-x)**2)
		if (rr >= rrold) exit
	end do
	resid=sqrt(rrold)
	do k=1,n
		y(k)=cof(k+1)-dot_product(z(1:k),cof(k:1:-1))
	end do
	cof(2:n+1)=y
	cof(n+2:2*n+1)=-z
	END SUBROUTINE pade

! -----------------------------------------------------------

	FUNCTION pccheb(d)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,cumprod,geop
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: d
	REAL(SP), DIMENSION(size(d)) :: pccheb
	INTEGER(I4B) :: k,n
	REAL(SP), DIMENSION(size(d)) :: denom,numer,pow
	n=size(d)
	pccheb(1)=2.0_sp*d(1)
	pow=geop(1.0_sp,2.0_sp,n)
	numer(1)=1.0
	denom(1)=1.0
	denom(2:(n+3)/2)=cumprod(arth(1.0_sp,1.0_sp,(n+1)/2))
	pccheb(2:n)=0.0
	do k=2,n
		numer(2:(k+3)/2)=cumprod(arth(k-1.0_sp,-1.0_sp,(k+1)/2))
		pccheb(k:1:-2)=pccheb(k:1:-2)+&
			d(k)/pow(k-1)*numer(1:(k+1)/2)/denom(1:(k+1)/2)
	end do
	END FUNCTION pccheb

! -----------------------------------------------------------

	SUBROUTINE pcshft(a,b,d)
	use mo_nrtype
          use mo_nrutil, ONLY : geop
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
	INTEGER(I4B) :: j,n
	REAL(SP), DIMENSION(size(d)) :: dd
	REAL(SP) :: x
	n=size(d)
	dd=d*geop(1.0_sp,2.0_sp/(b-a),n)
	x=-0.5_sp*(a+b)
	d(1)=dd(n)
	d(2:n)=0.0
	do j=n-1,1,-1
		d(2:n+1-j)=d(2:n+1-j)*x+d(1:n-j)
		d(1)=d(1)*x+dd(j)
	end do
	END SUBROUTINE pcshft

! -----------------------------------------------------------

	SUBROUTINE pearsn(x,y,r,prob,z)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : betai
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: r,prob,z
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), PARAMETER :: TINY=1.0e-20_sp
	REAL(SP), DIMENSION(size(x)) :: xt,yt
	REAL(SP) :: ax,ay,df,sxx,sxy,syy,t
	INTEGER(I4B) :: n
	n=assert_eq(size(x),size(y),'pearsn')
	ax=sum(x)/n
	ay=sum(y)/n
	xt(:)=x(:)-ax
	yt(:)=y(:)-ay
	sxx=dot_product(xt,xt)
	syy=dot_product(yt,yt)
	sxy=dot_product(xt,yt)
	r=sxy/(sqrt(sxx*syy)+TINY)
	z=0.5_sp*log(((1.0_sp+r)+TINY)/((1.0_sp-r)+TINY))
	df=n-2
	t=r*sqrt(df/(((1.0_sp-r)+TINY)*((1.0_sp+r)+TINY)))
	prob=betai(0.5_sp*df,0.5_sp,df/(df+t**2))
!	prob=erfcc(abs(z*sqrt(n-1.0_sp))/SQRT2)
	END SUBROUTINE pearsn

! -----------------------------------------------------------

	SUBROUTINE period(x,y,ofac,hifac,px,py,jmax,prob)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,imaxloc
	use mo_nr, ONLY : avevar
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: jmax
	REAL(SP), INTENT(IN) :: ofac,hifac
	REAL(SP), INTENT(OUT) :: prob
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), DIMENSION(:), POINTER :: px,py
	INTEGER(I4B) :: i,n,nout
	REAL(SP) :: ave,cwtau,effm,expy,pnow,sumc,sumcy,&
		sums,sumsh,sumsy,swtau,var,wtau,xave,xdif,xmax,xmin
	REAL(DP), DIMENSION(size(x)) :: tmp1,tmp2,wi,wpi,wpr,wr
	LOGICAL(LGT), SAVE :: init=.true.
	n=assert_eq(size(x),size(y),'period')
	if (init) then
		init=.false.
		nullify(px,py)
	else
		if (associated(px)) deallocate(px)
		if (associated(py)) deallocate(py)
	end if
	nout=0.5_sp*ofac*hifac*n
	allocate(px(nout),py(nout))
	call avevar(y(:),ave,var)
	xmax=maxval(x(:))
	xmin=minval(x(:))
	xdif=xmax-xmin
	xave=0.5_sp*(xmax+xmin)
	pnow=1.0_sp/(xdif*ofac)
	tmp1(:)=TWOPI_D*((x(:)-xave)*pnow)
	wpr(:)=-2.0_dp*sin(0.5_dp*tmp1)**2
	wpi(:)=sin(tmp1(:))
	wr(:)=cos(tmp1(:))
	wi(:)=wpi(:)
	do i=1,nout
		px(i)=pnow
		sumsh=dot_product(wi,wr)
		sumc=dot_product(wr(:)-wi(:),wr(:)+wi(:))
		wtau=0.5_sp*atan2(2.0_sp*sumsh,sumc)
		swtau=sin(wtau)
		cwtau=cos(wtau)
		tmp1(:)=wi(:)*cwtau-wr(:)*swtau
		tmp2(:)=wr(:)*cwtau+wi(:)*swtau
		sums=dot_product(tmp1,tmp1)
		sumc=dot_product(tmp2,tmp2)
		sumsy=dot_product(y(:)-ave,tmp1)
		sumcy=dot_product(y(:)-ave,tmp2)
		tmp1(:)=wr(:)
		wr(:)=(wr(:)*wpr(:)-wi(:)*wpi(:))+wr(:)
		wi(:)=(wi(:)*wpr(:)+tmp1(:)*wpi(:))+wi(:)
		py(i)=0.5_sp*(sumcy**2/sumc+sumsy**2/sums)/var
		pnow=pnow+1.0_sp/(ofac*xdif)
	end do
	jmax=imaxloc(py(1:nout))
	expy=exp(-py(jmax))
	effm=2.0_sp*nout/ofac
	prob=effm*expy
	if (prob > 0.01_sp) prob=1.0_sp-(1.0_sp-expy)**effm
	END SUBROUTINE period

! -----------------------------------------------------------

	FUNCTION plgndr_s(l,m,x)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: l,m
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: plgndr_s
	INTEGER(I4B) :: ll
	REAL(SP) :: pll,pmm,pmmp1,somx2
	call assert(m >= 0, m <= l, abs(x) <= 1.0, 'plgndr_s args')
	pmm=1.0
	if (m > 0) then
		somx2=sqrt((1.0_sp-x)*(1.0_sp+x))
		pmm=product(arth(1.0_sp,2.0_sp,m))*somx2**m
		if (mod(m,2) == 1) pmm=-pmm
	end if
	if (l == m) then
		plgndr_s=pmm
	else
		pmmp1=x*(2*m+1)*pmm
		if (l == m+1) then
			plgndr_s=pmmp1
		else
			do ll=m+2,l
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
				pmm=pmmp1
				pmmp1=pll
			end do
			plgndr_s=pll
		end if
	end if
	END FUNCTION plgndr_s


	FUNCTION plgndr_v(l,m,x)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: l,m
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: plgndr_v
	INTEGER(I4B) :: ll
	REAL(SP), DIMENSION(size(x)) :: pll,pmm,pmmp1,somx2
	call assert(m >= 0, m <= l, all(abs(x) <= 1.0), 'plgndr_v args')
	pmm=1.0
	if (m > 0) then
		somx2=sqrt((1.0_sp-x)*(1.0_sp+x))
		pmm=product(arth(1.0_sp,2.0_sp,m))*somx2**m
		if (mod(m,2) == 1) pmm=-pmm
	end if
	if (l == m) then
		plgndr_v=pmm
	else
		pmmp1=x*(2*m+1)*pmm
		if (l == m+1) then
			plgndr_v=pmmp1
		else
			do ll=m+2,l
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
				pmm=pmmp1
				pmmp1=pll
			end do
			plgndr_v=pll
		end if
	end if
	END FUNCTION plgndr_v

! -----------------------------------------------------------

	FUNCTION poidev(xm)
	use mo_nrtype
	use mo_nr, ONLY : gammln,ran1
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: xm
	REAL(SP) :: poidev
	REAL(SP) :: em,harvest,t,y
	REAL(SP), SAVE :: alxm,g,oldm=-1.0_sp,sq
	if (xm < 12.0) then
		if (xm /= oldm) then
			oldm=xm
			g=exp(-xm)
		end if
		em=-1
		t=1.0
		do
			em=em+1.0_sp
			call ran1(harvest)
			t=t*harvest
			if (t <= g) exit
		end do
	else
		if (xm /= oldm) then
			oldm=xm
			sq=sqrt(2.0_sp*xm)
			alxm=log(xm)
			g=xm*alxm-gammln(xm+1.0_sp)
		end if
		do
			do
				call ran1(harvest)
				y=tan(PI*harvest)
				em=sq*y+xm
				if (em >= 0.0) exit
			end do
			em=int(em)
			t=0.9_sp*(1.0_sp+y**2)*exp(em*alxm-gammln(em+1.0_sp)-g)
			call ran1(harvest)
			if (harvest <= t) exit
		end do
	end if
	poidev=em
	END FUNCTION poidev

! -----------------------------------------------------------

	FUNCTION polcoe(x,y)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,outerdiff
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), DIMENSION(size(x)) :: polcoe
	INTEGER(I4B) :: i,k,n
	REAL(SP), DIMENSION(size(x)) :: s
	REAL(SP), DIMENSION(size(x),size(x)) :: a
	n=assert_eq(size(x),size(y),'polcoe')
	s=0.0
	s(n)=-x(1)
	do i=2,n
		s(n+1-i:n-1)=s(n+1-i:n-1)-x(i)*s(n+2-i:n)
		s(n)=s(n)-x(i)
	end do
	a=outerdiff(x,x)
	polcoe=product(a,dim=2,mask=a /= 0.0)
	a(:,1)=-s(1)/x(:)
	do k=2,n
		a(:,k)=-(s(k)-a(:,k-1))/x(:)
	end do
	s=y/polcoe
	polcoe=matmul(s,a)
	END FUNCTION polcoe

! -----------------------------------------------------------

	FUNCTION polcof(xa,ya)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,iminloc
	use mo_nr, ONLY : polint
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
	REAL(SP), DIMENSION(size(xa)) :: polcof
	INTEGER(I4B) :: j,k,m,n
	REAL(SP) :: dy
	REAL(SP), DIMENSION(size(xa)) :: x,y
	n=assert_eq(size(xa),size(ya),'polcof')
	x=xa
	y=ya
	do j=1,n
		m=n+1-j
		call polint(x(1:m),y(1:m),0.0_sp,polcof(j),dy)
		k=iminloc(abs(x(1:m)))
		where (x(1:m) /= 0.0) y(1:m)=(y(1:m)-polcof(j))/x(1:m)
		y(k:m-1)=y(k+1:m)
		x(k:m-1)=x(k+1:m)
	end do
	END FUNCTION polcof

! -----------------------------------------------------------

	SUBROUTINE poldiv(u,v,q,r)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: u,v
	REAL(SP), DIMENSION(:), INTENT(OUT) :: q,r
	INTEGER(I4B) :: i,n,nv
	n=assert_eq(size(u),size(q),size(r),'poldiv')
	nv=size(v)
	r(:)=u(:)
	q(:)=0.0
	do i=n-nv,0,-1
		q(i+1)=r(nv+i)/v(nv)
		r(i+1:nv+i-1)=r(i+1:nv+i-1)-q(i+1)*v(1:nv-1)
	end do
	r(nv:n)=0.0
	END SUBROUTINE poldiv

! -----------------------------------------------------------

	SUBROUTINE polin2(x1a,x2a,ya,x1,x2,y,dy)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : polint
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
	REAL(SP), INTENT(IN) :: x1,x2
	REAL(SP), INTENT(OUT) :: y,dy
	INTEGER(I4B) :: j,m,ndum
	REAL(SP), DIMENSION(size(x1a)) :: ymtmp
	REAL(SP), DIMENSION(size(x2a)) :: yntmp
	m=assert_eq(size(x1a),size(ya,1),'polin2: m')
	ndum=assert_eq(size(x2a),size(ya,2),'polin2: ndum')
	do j=1,m
		yntmp=ya(j,:)
		call polint(x2a,yntmp,x2,ymtmp(j),dy)
	end do
	call polint(x1a,ymtmp,x1,y,dy)
	END SUBROUTINE polin2

! -----------------------------------------------------------

	SUBROUTINE polint(xa,ya,x,y,dy)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,iminloc,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
	REAL(SP), INTENT(IN) :: x
	REAL(SP), INTENT(OUT) :: y,dy
	INTEGER(I4B) :: m,n,ns
	REAL(SP), DIMENSION(size(xa)) :: c,d,den,ho
	n=assert_eq(size(xa),size(ya),'polint')
	c=ya
	d=ya
	ho=xa-x
	ns=iminloc(abs(x-xa))
	y=ya(ns)
	ns=ns-1
	do m=1,n-1
		den(1:n-m)=ho(1:n-m)-ho(1+m:n)
		if (any(den(1:n-m) == 0.0)) &
			call nrerror('polint: calculation failure')
		den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
		d(1:n-m)=ho(1+m:n)*den(1:n-m)
		c(1:n-m)=ho(1:n-m)*den(1:n-m)
		if (2*ns < n-m) then
			dy=c(ns+1)
		else
			dy=d(ns)
			ns=ns-1
		end if
		y=y+dy
	end do
	END SUBROUTINE polint

! -----------------------------------------------------------

	SUBROUTINE powell(p,xi,ftol,iter,fret)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : linmin
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: xi
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(SP), INTENT(IN) :: ftol
	REAL(SP), INTENT(OUT) :: fret
	INTERFACE
		FUNCTION func(p)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=200
	INTEGER(I4B) :: i,ibig,n
	REAL(SP) :: del,fp,fptt,t
	REAL(SP), DIMENSION(size(p)) :: pt,ptt,xit
	n=assert_eq(size(p),size(xi,1),size(xi,2),'powell')
	fret=func(p)
	pt(:)=p(:)
	iter=0
	do
		iter=iter+1
		fp=fret
		ibig=0
		del=0.0
		do i=1,n
			xit(:)=xi(:,i)
			fptt=fret
			call linmin(p,xit,fret)
			if (abs(fptt-fret) > del) then
				del=abs(fptt-fret)
				ibig=i
			end if
		end do
		if (2.0_sp*abs(fp-fret) <= ftol*(abs(fp)+abs(fret))) RETURN
		if (iter == ITMAX) call &
			nrerror('powell exceeding maximum iterations')
		ptt(:)=2.0_sp*p(:)-pt(:)
		xit(:)=p(:)-pt(:)
		pt(:)=p(:)
		fptt=func(ptt)
		if (fptt >= fp) cycle
		t=2.0_sp*(fp-2.0_sp*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
		if (t >= 0.0) cycle
		call linmin(p,xit,fret)
		xi(:,ibig)=xi(:,n)
		xi(:,n)=xit(:)
	end do
	END SUBROUTINE powell

! -----------------------------------------------------------

	FUNCTION predic(data,d,nfut)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data,d
	INTEGER(I4B), INTENT(IN) :: nfut
	REAL(SP), DIMENSION(nfut) :: predic
	INTEGER(I4B) :: j,ndata,m
	REAL(SP) :: discrp,sm
	REAL(SP), DIMENSION(size(d)) :: reg
	m=size(d)
	ndata=size(data)
	reg(1:m)=data(ndata:ndata+1-m:-1)
	do j=1,nfut
		discrp=0.0
		sm=discrp+dot_product(d,reg)
		reg=eoshift(reg,-1,sm)
		predic(j)=sm
	end do
	END FUNCTION predic

! -----------------------------------------------------------

	FUNCTION probks(alam)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: alam
	REAL(SP) :: probks
	REAL(SP), PARAMETER :: EPS1=0.001_sp,EPS2=1.0e-8_sp
	INTEGER(I4B), PARAMETER :: NITER=100
	INTEGER(I4B) :: j
	REAL(SP) :: a2,fac,term,termbf
	a2=-2.0_sp*alam**2
	fac=2.0
	probks=0.0
	termbf=0.0
	do j=1,NITER
		term=fac*exp(a2*j**2)
		probks=probks+term
		if (abs(term) <= EPS1*termbf .or. abs(term) <= EPS2*probks) RETURN
		fac=-fac
		termbf=abs(term)
	end do
	probks=1.0
	END FUNCTION probks

! -----------------------------------------------------------

	SUBROUTINE psdes_s(lword,rword)
	use mo_nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: lword,rword
	INTEGER(I4B), PARAMETER :: NITER=4
	INTEGER(I4B), DIMENSION(4), SAVE :: C1,C2
	DATA C1 /Z'BAA96887',Z'1E17D32C',Z'03BCDC3C',Z'0F33D1B2'/
	DATA C2 /Z'4B0F3B58',Z'E874F0C3',Z'6955C5A6',Z'55A7CA46'/
	INTEGER(I4B) :: i,ia,ib,iswap,itmph,itmpl
	do i=1,NITER
		iswap=rword
		ia=ieor(rword,C1(i))
		itmpl=iand(ia,65535)
		itmph=iand(ishft(ia,-16),65535)
		ib=itmpl**2+not(itmph**2)
		ia=ior(ishft(ib,16),iand(ishft(ib,-16),65535))
		rword=ieor(lword,ieor(C2(i),ia)+itmpl*itmph)
		lword=iswap
	end do
	END SUBROUTINE psdes_s

	SUBROUTINE psdes_v(lword,rword)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: lword,rword
	INTEGER(I4B), PARAMETER :: NITER=4
	INTEGER(I4B), DIMENSION(4), SAVE :: C1,C2
	DATA C1 /Z'BAA96887',Z'1E17D32C',Z'03BCDC3C',Z'0F33D1B2'/
	DATA C2 /Z'4B0F3B58',Z'E874F0C3',Z'6955C5A6',Z'55A7CA46'/
	INTEGER(I4B), DIMENSION(size(lword)) :: ia,ib,iswap,itmph,itmpl
	INTEGER(I4B) :: i
	i=assert_eq(size(lword),size(rword),'psdes_v')
	do i=1,NITER
		iswap=rword
		ia=ieor(rword,C1(i))
		itmpl=iand(ia,65535)
		itmph=iand(ishft(ia,-16),65535)
		ib=itmpl**2+not(itmph**2)
		ia=ior(ishft(ib,16),iand(ishft(ib,-16),65535))
		rword=ieor(lword,ieor(C2(i),ia)+itmpl*itmph)
		lword=iswap
	end do
	END SUBROUTINE psdes_v

	SUBROUTINE psdes_safe_s(lword,rword)
	use mo_nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: lword,rword
	INTEGER(I4B), PARAMETER :: NITER=4
	INTEGER(I4B), DIMENSION(4), SAVE :: C1,C2
	DATA C1 /Z'BAA96887',Z'1E17D32C',Z'03BCDC3C',Z'0F33D1B2'/
	DATA C2 /Z'4B0F3B58',Z'E874F0C3',Z'6955C5A6',Z'55A7CA46'/
	INTEGER(I4B) :: i,ia,ib,iswap
	REAL(DP) :: alo,ahi
	do i=1,NITER
		iswap=rword
		ia=ieor(rword,C1(i))
		alo=real(iand(ia,65535),dp)
		ahi=real(iand(ishft(ia,-16),65535),dp)
		ib=modint(alo*alo+real(not(modint(ahi*ahi)),dp))
		ia=ior(ishft(ib,16),iand(ishft(ib,-16),65535))
		rword=ieor(lword,modint(real(ieor(C2(i),ia),dp)+alo*ahi))
		lword=iswap
	end do
	CONTAINS
!BL
	FUNCTION modint(x)
	REAL(DP), INTENT(IN) :: x
	INTEGER(I4B) :: modint
	REAL(DP) :: a
	REAL(DP), PARAMETER :: big=huge(modint), base=big+big+2.0_dp
	a=modulo(x,base)
	if (a > big) a=a-base
	modint=nint(a,kind=i4b)
	END FUNCTION modint
	END SUBROUTINE psdes_safe_s

	SUBROUTINE psdes_safe_v(lword,rword)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: lword,rword
	INTEGER(I4B), PARAMETER :: NITER=4
	INTEGER(I4B), SAVE :: C1(4),C2(4)
	DATA C1 /Z'BAA96887',Z'1E17D32C',Z'03BCDC3C',Z'0F33D1B2'/
	DATA C2 /Z'4B0F3B58',Z'E874F0C3',Z'6955C5A6',Z'55A7CA46'/
	INTEGER(I4B), DIMENSION(size(lword)) :: ia,ib,iswap
	REAL(DP), DIMENSION(size(lword)) :: alo,ahi
	INTEGER(I4B) :: i
	i=assert_eq(size(lword),size(rword),'psdes_safe_v')
	do i=1,NITER
		iswap=rword
		ia=ieor(rword,C1(i))
		alo=real(iand(ia,65535),dp)
		ahi=real(iand(ishft(ia,-16),65535),dp)
		ib=modint(alo*alo+real(not(modint(ahi*ahi)),dp))
		ia=ior(ishft(ib,16),iand(ishft(ib,-16),65535))
		rword=ieor(lword,modint(real(ieor(C2(i),ia),dp)+alo*ahi))
		lword=iswap
	end do
	CONTAINS
!BL
	FUNCTION modint(x)
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	INTEGER(I4B), DIMENSION(size(x)) :: modint
	REAL(DP), DIMENSION(size(x)) :: a
	REAL(DP), PARAMETER :: big=huge(modint), base=big+big+2.0_dp
	a=modulo(x,base)
	where (a > big) a=a-base
	modint=nint(a,kind=i4b)
	END FUNCTION modint
	END SUBROUTINE psdes_safe_v

! -----------------------------------------------------------

	SUBROUTINE pwt(a,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,nrerror
	USE pwtcom
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
	INTEGER(I4B), INTENT(IN) :: isign
	REAL(SP), DIMENSION(size(a)) :: wksp
	INTEGER(I4B), DIMENSION(size(a)/2) :: jf,jr
	INTEGER(I4B) :: k,n,nh,nmod
	n=size(a)
	if (n < 4) RETURN
	if (ncof == 0) call nrerror('pwt: must call pwtset before pwt')
	nmod=ncof*n
	nh=n/2
	wksp(:)=0.0
	jf=iand(n-1,arth(2+nmod+ioff,2,nh))
	jr=iand(n-1,arth(2+nmod+joff,2,nh))
	do k=1,ncof
		if (isign >= 0) then
			wksp(1:nh)=wksp(1:nh)+cc(k)*a(jf+1)
			wksp(nh+1:n)=wksp(nh+1:n)+cr(k)*a(jr+1)
		else
			wksp(jf+1)=wksp(jf+1)+cc(k)*a(1:nh)
			wksp(jr+1)=wksp(jr+1)+cr(k)*a(nh+1:n)
		end if
		if (k == ncof) exit
		jf=iand(n-1,jf+1)
		jr=iand(n-1,jr+1)
	end do
	a(:)=wksp(:)
	END SUBROUTINE pwt

! -----------------------------------------------------------

MODULE pwtcom
	use mo_nrtype
	INTEGER(I4B), SAVE :: ncof=0,ioff,joff
	REAL(SP), DIMENSION(:), ALLOCATABLE, SAVE :: cc,cr
END MODULE pwtcom

	SUBROUTINE pwtset(n)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	USE pwtcom
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP) :: sig
	REAL(SP), PARAMETER :: &
		c4(4)=(/&
		0.4829629131445341_sp, 0.8365163037378079_sp, &
		0.2241438680420134_sp,-0.1294095225512604_sp /), &
		c12(12)=(/&
		0.111540743350_sp, 0.494623890398_sp, 0.751133908021_sp, &
		0.315250351709_sp,-0.226264693965_sp,-0.129766867567_sp, &
		0.097501605587_sp, 0.027522865530_sp,-0.031582039318_sp, &
		0.000553842201_sp, 0.004777257511_sp,-0.001077301085_sp /), &
		c20(20)=(/&
		0.026670057901_sp, 0.188176800078_sp, 0.527201188932_sp, &
		0.688459039454_sp, 0.281172343661_sp,-0.249846424327_sp, &
		-0.195946274377_sp, 0.127369340336_sp, 0.093057364604_sp, &
		-0.071394147166_sp,-0.029457536822_sp, 0.033212674059_sp, &
		0.003606553567_sp,-0.010733175483_sp, 0.001395351747_sp, &
		0.001992405295_sp,-0.000685856695_sp,-0.000116466855_sp, &
		0.000093588670_sp,-0.000013264203_sp /)
	if (allocated(cc)) deallocate(cc)
	if (allocated(cr)) deallocate(cr)
	allocate(cc(n),cr(n))
	ncof=n
	ioff=-n/2
	joff=-n/2
	sig=-1.0
	select case(n)
		case(4)
			cc=c4
		case(12)
			cc=c12
		case(20)
			cc=c20
		case default
			call nrerror('unimplemented value n in pwtset')
	end select
	cr(n:1:-1) = cc
	cr(n:1:-2) = -cr(n:1:-2)
	END SUBROUTINE pwtset

! -----------------------------------------------------------

	FUNCTION pythag_sp(a,b)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: pythag_sp
	REAL(SP) :: absa,absb
	absa=abs(a)
	absb=abs(b)
	if (absa > absb) then
		pythag_sp=absa*sqrt(1.0_sp+(absb/absa)**2)
	else
		if (absb == 0.0) then
			pythag_sp=0.0
		else
			pythag_sp=absb*sqrt(1.0_sp+(absa/absb)**2)
		end if
	end if
	END FUNCTION pythag_sp

	FUNCTION pythag_dp(a,b)
	use mo_nrtype
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,b
	REAL(DP) :: pythag_dp
	REAL(DP) :: absa,absb
	absa=abs(a)
	absb=abs(b)
	if (absa > absb) then
		pythag_dp=absa*sqrt(1.0_dp+(absb/absa)**2)
	else
		if (absb == 0.0) then
			pythag_dp=0.0
		else
			pythag_dp=absb*sqrt(1.0_dp+(absa/absb)**2)
		end if
	end if
	END FUNCTION pythag_dp

! -----------------------------------------------------------

	SUBROUTINE pzextr(iest,xest,yest,yz,dy)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: iest
	REAL(SP), INTENT(IN) :: xest
	REAL(SP), DIMENSION(:), INTENT(IN) :: yest
	REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
	INTEGER(I4B), PARAMETER :: IEST_MAX=16
	INTEGER(I4B) :: j,nv
	INTEGER(I4B), SAVE :: nvold=-1
	REAL(SP) :: delta,f1,f2
	REAL(SP), DIMENSION(size(yz)) :: d,tmp,q
	REAL(SP), DIMENSION(IEST_MAX), SAVE :: x
	REAL(SP), DIMENSION(:,:), ALLOCATABLE, SAVE :: qcol
	nv=assert_eq(size(yz),size(yest),size(dy),'pzextr')
	if (iest > IEST_MAX) call &
		nrerror('pzextr: probable misuse, too much extrapolation')
	if (nv /= nvold) then
		if (allocated(qcol)) deallocate(qcol)
		allocate(qcol(nv,IEST_MAX))
		nvold=nv
	end if
	x(iest)=xest
	dy(:)=yest(:)
	yz(:)=yest(:)
	if (iest == 1) then
		qcol(:,1)=yest(:)
	else
		d(:)=yest(:)
		do j=1,iest-1
			delta=1.0_sp/(x(iest-j)-xest)
			f1=xest*delta
			f2=x(iest-j)*delta
			q(:)=qcol(:,j)
			qcol(:,j)=dy(:)
			tmp(:)=d(:)-q(:)
			dy(:)=f1*tmp(:)
			d(:)=f2*tmp(:)
			yz(:)=yz(:)+dy(:)
		end do
		qcol(:,iest)=dy(:)
	end if
	END SUBROUTINE pzextr

! -----------------------------------------------------------

	SUBROUTINE qrdcmp(a,c,d,sing)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,outerprod,vabs
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: c,d
	LOGICAL(LGT), INTENT(OUT) :: sing
	INTEGER(I4B) :: k,n
	REAL(SP) :: scale,sigma
	n=assert_eq(size(a,1),size(a,2),size(c),size(d),'qrdcmp')
	sing=.false.
	do k=1,n-1
		scale=maxval(abs(a(k:n,k)))
		if (scale == 0.0) then
			sing=.true.
			c(k)=0.0
			d(k)=0.0
		else
			a(k:n,k)=a(k:n,k)/scale
			sigma=sign(vabs(a(k:n,k)),a(k,k))
			a(k,k)=a(k,k)+sigma
			c(k)=sigma*a(k,k)
			d(k)=-scale*sigma
			a(k:n,k+1:n)=a(k:n,k+1:n)-outerprod(a(k:n,k),&
				matmul(a(k:n,k),a(k:n,k+1:n)))/c(k)
		end if
	end do
	d(n)=a(n,n)
	if (d(n) == 0.0) sing=.true.
	END SUBROUTINE qrdcmp

! -----------------------------------------------------------

	FUNCTION qromb(func,a,b)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : polint,trapzd
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qromb
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	REAL(SP), DIMENSION(JMAXP) :: h,s
	REAL(SP) :: dqromb
	INTEGER(I4B) :: j
	h(1)=1.0
	do j=1,JMAX
		call trapzd(func,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromb,dqromb)
			if (abs(dqromb) <= EPS*abs(qromb)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=0.25_sp*h(j)
	end do
	call nrerror('qromb: too many steps')
	END FUNCTION qromb

! -----------------------------------------------------------

	FUNCTION qromo(func,a,b,choose)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : polint
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qromo
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
!BL
		SUBROUTINE choose(funk,aa,bb,s,n)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: aa,bb
		REAL(SP), INTENT(INOUT) :: s
		INTEGER(I4B), INTENT(IN) :: n
		INTERFACE
			FUNCTION funk(x)
			use mo_nrtype
			IMPLICIT NONE
			REAL(SP), DIMENSION(:), INTENT(IN) :: x
			REAL(SP), DIMENSION(size(x)) :: funk
			END FUNCTION funk
		END INTERFACE
		END SUBROUTINE choose
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=14,JMAXP=JMAX+1,K=5,KM=K-1
	REAL(SP), PARAMETER :: EPS=1.0e-6
	REAL(SP), DIMENSION(JMAXP) :: h,s
	REAL(SP) :: dqromo
	INTEGER(I4B) :: j
	h(1)=1.0
	do j=1,JMAX
		call choose(func,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromo,dqromo)
			if (abs(dqromo) <= EPS*abs(qromo)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=h(j)/9.0_sp
	end do
	call nrerror('qromo: too many steps')
	END FUNCTION qromo

! -----------------------------------------------------------

	SUBROUTINE qroot(p,b,c,eps)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : poldiv
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: p
	REAL(SP), INTENT(INOUT) :: b,c
	REAL(SP), INTENT(IN) :: eps
	INTEGER(I4B), PARAMETER :: ITMAX=20
	REAL(SP), PARAMETER :: TINY=1.0e-6_sp
	INTEGER(I4B) :: iter,n
	REAL(SP) :: delb,delc,div,r,rb,rc,s,sb,sc
	REAL(SP), DIMENSION(3) :: d
	REAL(SP), DIMENSION(size(p)) :: q,qq,rem
	n=size(p)
	d(3)=1.0
	do iter=1,ITMAX
		d(2)=b
		d(1)=c
		call poldiv(p,d,q,rem)
		s=rem(1)
		r=rem(2)
		call poldiv(q(1:n-1),d(:),qq(1:n-1),rem(1:n-1))
		sc=-rem(1)
		rc=-rem(2)
		sb=-c*rc
		rb=sc-b*rc
		div=1.0_sp/(sb*rc-sc*rb)
		delb=(r*sc-s*rc)*div
		delc=(-r*sb+s*rb)*div
		b=b+delb
		c=c+delc
		if ((abs(delb) <= eps*abs(b) .or. abs(b) < TINY) .and. &
			(abs(delc) <= eps*abs(c) .or. abs(c) < TINY)) RETURN
	end do
	call nrerror('qroot: too many iterations')
	END SUBROUTINE qroot

! -----------------------------------------------------------

	SUBROUTINE qrsolv(a,c,d,b)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : rsolv
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	REAL(SP), DIMENSION(:), INTENT(IN) :: c,d
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: j,n
	REAL(SP) :: tau
	n=assert_eq((/size(a,1),size(a,2),size(b),size(c),size(d)/),'qrsolv')
	do j=1,n-1
		tau=dot_product(a(j:n,j),b(j:n))/c(j)
		b(j:n)=b(j:n)-tau*a(j:n,j)
	end do
	call rsolv(a,d,b)
	END SUBROUTINE qrsolv

! -----------------------------------------------------------

	SUBROUTINE qrupdt(r,qt,u,v)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,ifirstloc
	use mo_nr, ONLY : rotate,pythag
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: r,qt
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: u
	REAL(SP), DIMENSION(:), INTENT(IN) :: v
	INTEGER(I4B) :: i,k,n
	n=assert_eq((/size(r,1),size(r,2),size(qt,1),size(qt,2),size(u),&
		size(v)/),'qrupdt')
	k=n+1-ifirstloc(u(n:1:-1) /= 0.0)
	if (k < 1) k=1
	do i=k-1,1,-1
		call rotate(r,qt,i,u(i),-u(i+1))
		u(i)=pythag(u(i),u(i+1))
	end do
	r(1,:)=r(1,:)+u(1)*v
	do i=1,k-1
		call rotate(r,qt,i,r(i,i),-r(i+1,i))
	end do
	END SUBROUTINE qrupdt

! -----------------------------------------------------------

	FUNCTION qsimp(func,a,b)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : trapzd
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qsimp
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp, UNLIKELY=-1.0e30_sp
	INTEGER(I4B) :: j
	REAL(SP) :: os,ost,st
	ost=UNLIKELY
	os= UNLIKELY
	do j=1,JMAX
		call trapzd(func,a,b,st,j)
		qsimp=(4.0_sp*st-ost)/3.0_sp
		if (abs(qsimp-os) < EPS*abs(os)) RETURN
		if (qsimp == 0.0 .and. os == 0.0 .and. j > 6) RETURN
		os=qsimp
		ost=st
	end do
	call nrerror('qsimp: too many steps')
	END FUNCTION qsimp

! -----------------------------------------------------------

	FUNCTION qtrap(func,a,b)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : trapzd
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qtrap
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp, UNLIKELY=-1.0e30_sp
	REAL(SP) :: olds
	INTEGER(I4B) :: j
	olds=UNLIKELY
	do j=1,JMAX
		call trapzd(func,a,b,qtrap,j)
		if (abs(qtrap-olds) < EPS*abs(olds)) RETURN
		if (qtrap == 0.0 .and. olds == 0.0 .and. j > 6) RETURN
		olds=qtrap
	end do
	call nrerror('qtrap: too many steps')
	END FUNCTION qtrap

! -----------------------------------------------------------

	MODULE quad3d_qgaus_mod
	use mo_nrtype
	PRIVATE
	PUBLIC quad3d_qgaus
	REAL(SP) :: xsav,ysav
	INTERFACE
		FUNCTION func(x,y,z)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP), DIMENSION(:), INTENT(IN) :: z
		REAL(SP), DIMENSION(size(z)) :: func
		END FUNCTION func
!BL
		FUNCTION y1(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: y1
		END FUNCTION y1
!BL
		FUNCTION y2(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: y2
		END FUNCTION y2
!BL
		FUNCTION z1(x,y)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP) :: z1
		END FUNCTION z1
!BL
		FUNCTION z2(x,y)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP) :: z2
		END FUNCTION z2
	END INTERFACE
	CONTAINS
!BL
	FUNCTION h(x)
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: h
	INTEGER(I4B) :: i
	do i=1,size(x)
		xsav=x(i)
		h(i)=qgaus(g,y1(xsav),y2(xsav))
	end do
	END FUNCTION h
!BL
	FUNCTION g(y)
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(size(y)) :: g
	INTEGER(I4B) :: j
	do j=1,size(y)
		ysav=y(j)
		g(j)=qgaus(f,z1(xsav,ysav),z2(xsav,ysav))
	end do
	END FUNCTION g
!BL
	FUNCTION f(z)
	REAL(SP), DIMENSION(:), INTENT(IN) :: z
	REAL(SP), DIMENSION(size(z)) :: f
	f=func(xsav,ysav,z)
	END FUNCTION f
!BL
	RECURSIVE FUNCTION qgaus(func,a,b)
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qgaus
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP) :: xm,xr
	REAL(SP), DIMENSION(5) :: dx, w = (/ 0.2955242247_sp,0.2692667193_sp,&
		0.2190863625_sp,0.1494513491_sp,0.0666713443_sp /),&
		x = (/ 0.1488743389_sp,0.4333953941_sp,0.6794095682_sp,&
		0.8650633666_sp,0.9739065285_sp /)
	xm=0.5_sp*(b+a)
	xr=0.5_sp*(b-a)
	dx(:)=xr*x(:)
	qgaus=xr*sum(w(:)*(func(xm+dx)+func(xm-dx)))
	END FUNCTION qgaus
!BL
	SUBROUTINE quad3d_qgaus(x1,x2,ss)
	REAL(SP), INTENT(IN) :: x1,x2
	REAL(SP), INTENT(OUT) :: ss
	ss=qgaus(h,x1,x2)
	END SUBROUTINE quad3d_qgaus
	END MODULE quad3d_qgaus_mod

	MODULE quad3d_qromb_mod
	use mo_nrtype
	PRIVATE
	PUBLIC quad3d_qromb
	REAL(SP) :: xsav,ysav
	INTERFACE
		FUNCTION func(x,y,z)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP), DIMENSION(:), INTENT(IN) :: z
		REAL(SP), DIMENSION(size(z)) :: func
		END FUNCTION func
!BL
		FUNCTION y1(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: y1
		END FUNCTION y1
!BL
		FUNCTION y2(x)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: y2
		END FUNCTION y2
!BL
		FUNCTION z1(x,y)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP) :: z1
		END FUNCTION z1
!BL
		FUNCTION z2(x,y)
		use mo_nrtype
		REAL(SP), INTENT(IN) :: x,y
		REAL(SP) :: z2
		END FUNCTION z2
	END INTERFACE
	CONTAINS
!BL
	FUNCTION h(x)
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: h
	INTEGER(I4B) :: i
	do i=1,size(x)
		xsav=x(i)
		h(i)=qromb(g,y1(xsav),y2(xsav))
	end do
	END FUNCTION h
!BL
	FUNCTION g(y)
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(size(y)) :: g
	INTEGER(I4B) :: j
	do j=1,size(y)
		ysav=y(j)
		g(j)=qromb(f,z1(xsav,ysav),z2(xsav,ysav))
	end do
	END FUNCTION g
!BL
	FUNCTION f(z)
	REAL(SP), DIMENSION(:), INTENT(IN) :: z
	REAL(SP), DIMENSION(size(z)) :: f
	f=func(xsav,ysav,z)
	END FUNCTION f
!BL
	RECURSIVE FUNCTION qromb(func,a,b)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : polint
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP) :: qromb
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
	REAL(SP), PARAMETER :: EPS=3.0e-6_sp
	REAL(SP), DIMENSION(JMAXP) :: h,s
	REAL(SP) :: dqromb
	INTEGER(I4B) :: j
	h(1)=1.0
	do j=1,JMAX
		call trapzd(func,a,b,s(j),j)
		if (j >= K) then
			call polint(h(j-KM:j),s(j-KM:j),0.0_sp,qromb,dqromb)
			if (abs(dqromb) <= EPS*abs(qromb)) RETURN
		end if
		s(j+1)=s(j)
		h(j+1)=0.25_sp*h(j)
	end do
	call nrerror('qromb: too many steps')
	END FUNCTION qromb
!BL
	RECURSIVE SUBROUTINE trapzd(func,a,b,s,n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP) :: del,fsum
	INTEGER(I4B) :: it
	if (n == 1) then
		s=0.5_sp*(b-a)*sum(func( (/ a,b /) ))
	else
		it=2**(n-2)
		del=(b-a)/it
		fsum=sum(func(arth(a+0.5_sp*del,del,it)))
		s=0.5_sp*(s+del*fsum)
	end if
	END SUBROUTINE trapzd
!BL
	SUBROUTINE quad3d_qromb(x1,x2,ss)
	REAL(SP), INTENT(IN) :: x1,x2
	REAL(SP), INTENT(OUT) :: ss
	ss=qromb(h,x1,x2)
	END SUBROUTINE quad3d_qromb
	END MODULE quad3d_qromb_mod

! -----------------------------------------------------------

	SUBROUTINE quadct(x,y,xx,yy,fa,fb,fc,fd)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx,yy
	REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
	INTEGER(I4B) :: na,nb,nc,nd,nn
	REAL(SP) :: ff
	nn=assert_eq(size(xx),size(yy),'quadct')
	na=count(yy(:) > y .and. xx(:) > x)
	nb=count(yy(:) > y .and. xx(:) <= x)
	nc=count(yy(:) <= y .and. xx(:) <= x)
	nd=nn-na-nb-nc
	ff=1.0_sp/nn
	fa=ff*na
	fb=ff*nb
	fc=ff*nc
	fd=ff*nd
	END SUBROUTINE quadct

! -----------------------------------------------------------

	SUBROUTINE quadmx(a)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq,diagadd,outerprod
	use mo_nr, ONLY : wwghts,kermom
	USE kermom_info
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: a
	INTEGER(I4B) :: j,n
	REAL(SP) :: h,x
	REAL(SP), DIMENSION(size(a,1)) :: wt
	n=assert_eq(size(a,1),size(a,2),'quadmx')
	h=PI/(n-1)
	do j=1,n
		x=(j-1)*h
		kermom_x=x
		wt(:)=wwghts(n,h,kermom)
		a(j,:)=wt(:)
	end do
	wt(:)=cos(arth(0,1,n)*h)
	a(:,:)=a(:,:)*outerprod(wt(:),wt(:))
	call diagadd(a,1.0_sp)
	END SUBROUTINE quadmx

! -----------------------------------------------------------

	SUBROUTINE quadvl(x,y,fa,fb,fc,fd)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y
	REAL(SP), INTENT(OUT) :: fa,fb,fc,fd
	REAL(SP) :: qa,qb,qc,qd
	qa=min(2.0_sp,max(0.0_sp,1.0_sp-x))
	qb=min(2.0_sp,max(0.0_sp,1.0_sp-y))
	qc=min(2.0_sp,max(0.0_sp,x+1.0_sp))
	qd=min(2.0_sp,max(0.0_sp,y+1.0_sp))
	fa=0.25_sp*qa*qb
	fb=0.25_sp*qb*qc
	fc=0.25_sp*qc*qd
	fd=0.25_sp*qd*qa
	END SUBROUTINE quadvl

! -----------------------------------------------------------

	FUNCTION ran(idum)
	IMPLICIT NONE
	INTEGER, PARAMETER :: K4B=selected_int_kind(9)
	INTEGER(K4B), INTENT(INOUT) :: idum
	REAL :: ran
	INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
	REAL, SAVE :: am
	INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
	if (idum <= 0 .or. iy < 0) then
		am=nearest(1.0,-1.0)/IM
		iy=ior(ieor(888889999,abs(idum)),1)
		ix=ieor(777755555,abs(idum))
		idum=abs(idum)+1
	end if
	ix=ieor(ix,ishft(ix,13))
	ix=ieor(ix,ishft(ix,-17))
	ix=ieor(ix,ishft(ix,5))
	k=iy/IQ
	iy=IA*(iy-k*IQ)-IR*k
	if (iy < 0) iy=iy+IM
	ran=am*ior(iand(IM,ieor(ix,iy)),1)
	END FUNCTION ran

! -----------------------------------------------------------

	SUBROUTINE ran0_s(harvest)
	use mo_nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init,iran0,jran0,kran0,nran0,rans
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: harvest
	if (lenran < 1) call ran_init(1)
	rans=iran0-kran0
	if (rans < 0) rans=rans+2147483579_k4b
	iran0=jran0
	jran0=kran0
	kran0=rans
	nran0=ieor(nran0,ishft(nran0,13))
	nran0=ieor(nran0,ishft(nran0,-17))
	nran0=ieor(nran0,ishft(nran0,5))
	rans=ieor(nran0,rans)
	harvest=amm*merge(rans,not(rans), rans<0 )
	END SUBROUTINE ran0_s

	SUBROUTINE ran0_v(harvest)
	use mo_nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init,iran,jran,kran,nran,ranv
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
	INTEGER(K4B) :: n
	n=size(harvest)
	if (lenran < n+1) call ran_init(n+1)
	ranv(1:n)=iran(1:n)-kran(1:n)
	where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
	iran(1:n)=jran(1:n)
	jran(1:n)=kran(1:n)
	kran(1:n)=ranv(1:n)
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
	ranv(1:n)=ieor(nran(1:n),ranv(1:n))
	harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
	END SUBROUTINE ran0_v

! -----------------------------------------------------------

	SUBROUTINE ran1_s(harvest)
	use mo_nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
		iran0,jran0,kran0,nran0,mran0,rans
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: harvest
	if (lenran < 1) call ran_init(1)
	rans=iran0-kran0
	if (rans < 0) rans=rans+2147483579_k4b
	iran0=jran0
	jran0=kran0
	kran0=rans
	nran0=ieor(nran0,ishft(nran0,13))
	nran0=ieor(nran0,ishft(nran0,-17))
	nran0=ieor(nran0,ishft(nran0,5))
	if (nran0 == 1) nran0=270369_k4b
	mran0=ieor(mran0,ishft(mran0,5))
	mran0=ieor(mran0,ishft(mran0,-13))
	mran0=ieor(mran0,ishft(mran0,6))
	rans=ieor(nran0,rans)+mran0
	harvest=amm*merge(rans,not(rans), rans<0 )
	END SUBROUTINE ran1_s

	SUBROUTINE ran1_v(harvest)
	use mo_nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
		iran,jran,kran,nran,mran,ranv
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
	INTEGER(K4B) :: n
	n=size(harvest)
	if (lenran < n+1) call ran_init(n+1)
	ranv(1:n)=iran(1:n)-kran(1:n)
	where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
	iran(1:n)=jran(1:n)
	jran(1:n)=kran(1:n)
	kran(1:n)=ranv(1:n)
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
	where (nran(1:n) == 1) nran(1:n)=270369_k4b
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
	ranv(1:n)=ieor(nran(1:n),ranv(1:n))+mran(1:n)
	harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
	END SUBROUTINE ran1_v

! -----------------------------------------------------------

	SUBROUTINE ran2_s(harvest)
	use mo_nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
		iran0,jran0,kran0,nran0,mran0,rans
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: harvest
	if (lenran < 1) call ran_init(1)
	rans=iran0-kran0
	if (rans < 0) rans=rans+2147483579_k4b
	iran0=jran0
	jran0=kran0
	kran0=rans
	nran0=ieor(nran0,ishft(nran0,13))
	nran0=ieor(nran0,ishft(nran0,-17))
	nran0=ieor(nran0,ishft(nran0,5))
	rans=iand(mran0,65535)
	mran0=ishft(3533*ishft(mran0,-16)+rans,16)+ &
		3533*rans+820265819_k4b
	rans=ieor(nran0,kran0)+mran0
	harvest=amm*merge(rans,not(rans), rans<0 )
	END SUBROUTINE ran2_s

	SUBROUTINE ran2_v(harvest)
	use mo_nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
		iran,jran,kran,nran,mran,ranv
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
	INTEGER(K4B) :: n
	n=size(harvest)
	if (lenran < n+1) call ran_init(n+1)
	ranv(1:n)=iran(1:n)-kran(1:n)
	where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
	iran(1:n)=jran(1:n)
	jran(1:n)=kran(1:n)
	kran(1:n)=ranv(1:n)
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
	ranv(1:n)=iand(mran(1:n),65535)
	mran(1:n)=ishft(3533*ishft(mran(1:n),-16)+ranv(1:n),16)+ &
		3533*ranv(1:n)+820265819_k4b
	ranv(1:n)=ieor(nran(1:n),kran(1:n))+mran(1:n)
	harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
	END SUBROUTINE ran2_v

! -----------------------------------------------------------

	SUBROUTINE ran3_s(harvest)
	use mo_nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init,ran_hash,mran0,nran0,rans
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: harvest
	INTEGER(K4B) :: temp
	if (lenran < 1) call ran_init(1)
	nran0=ieor(nran0,ishft(nran0,13))
	nran0=ieor(nran0,ishft(nran0,-17))
	nran0=ieor(nran0,ishft(nran0,5))
	if (nran0 == 1) nran0=270369_k4b
	rans=nran0
	mran0=ieor(mran0,ishft(mran0,5))
	mran0=ieor(mran0,ishft(mran0,-13))
	mran0=ieor(mran0,ishft(mran0,6))
	temp=mran0
	call ran_hash(temp,rans)
	harvest=amm*merge(rans,not(rans), rans<0 )
	END SUBROUTINE ran3_s

	SUBROUTINE ran3_v(harvest)
	use mo_nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init,ran_hash,mran,nran,ranv
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
	INTEGER(K4B), DIMENSION(size(harvest)) :: temp
	INTEGER(K4B) :: n
	n=size(harvest)
	if (lenran < n+1) call ran_init(n+1)
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
	where (nran(1:n) == 1) nran(1:n)=270369_k4b
	ranv(1:n)=nran(1:n)
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
	temp=mran(1:n)
	call ran_hash(temp,ranv(1:n))
	harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
	END SUBROUTINE ran3_v

! -----------------------------------------------------------

MODULE ran_state
	use mo_nrtype
	IMPLICIT NONE
	INTEGER, PARAMETER :: K4B=selected_int_kind(9)
	INTEGER(K4B), PARAMETER :: hg=huge(1_K4B), hgm=-hg, hgng=hgm-1
	INTEGER(K4B), SAVE :: lenran=0, seq=0
	INTEGER(K4B), SAVE :: iran0,jran0,kran0,nran0,mran0,rans
	INTEGER(K4B), DIMENSION(:,:), POINTER, SAVE :: ranseeds
	INTEGER(K4B), DIMENSION(:), POINTER, SAVE :: iran,jran,kran, &
		nran,mran,ranv
	REAL(SP), SAVE :: amm
	INTERFACE ran_hash
		MODULE PROCEDURE ran_hash_s, ran_hash_v
	END INTERFACE
CONTAINS
!BL
	SUBROUTINE ran_init(length)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,nrerror,reallocate
	IMPLICIT NONE
	INTEGER(K4B), INTENT(IN) :: length
	INTEGER(K4B) :: new,j,hgt
	if (length < lenran) RETURN
	hgt=hg
	if (hg /= 2147483647) call nrerror('ran_init: arith assump 1 fails')
	if (hgng >= 0)        call nrerror('ran_init: arith assump 2 fails')
	if (hgt+1 /= hgng)    call nrerror('ran_init: arith assump 3 fails')
	if (not(hg) >= 0)     call nrerror('ran_init: arith assump 4 fails')
	if (not(hgng) < 0)    call nrerror('ran_init: arith assump 5 fails')
	if (hg+hgng >= 0)     call nrerror('ran_init: arith assump 6 fails')
	if (not(-1_k4b) < 0)  call nrerror('ran_init: arith assump 7 fails')
	if (not(0_k4b) >= 0)  call nrerror('ran_init: arith assump 8 fails')
	if (not(1_k4b) >= 0)  call nrerror('ran_init: arith assump 9 fails')
	if (lenran > 0) then
		ranseeds=>reallocate(ranseeds,length,5)
		ranv=>reallocate(ranv,length-1)
		new=lenran+1
	else
		allocate(ranseeds(length,5))
		allocate(ranv(length-1))
		new=1
		amm=nearest(1.0_sp,-1.0_sp)/hgng
		if (amm*hgng >= 1.0 .or. amm*hgng <= 0.0) &
			call nrerror('ran_init: arth assump 10 fails')
	end if
	ranseeds(new:,1)=seq
	ranseeds(new:,2:5)=spread(arth(new,1,size(ranseeds(new:,1))),2,4)
	do j=1,4
		call ran_hash(ranseeds(new:,j),ranseeds(new:,j+1))
	end do
	where (ranseeds(new:,1:3) < 0) &
		ranseeds(new:,1:3)=not(ranseeds(new:,1:3))
	where (ranseeds(new:,4:5) == 0) ranseeds(new:,4:5)=1
	if (new == 1) then
		iran0=ranseeds(1,1)
		jran0=ranseeds(1,2)
		kran0=ranseeds(1,3)
		mran0=ranseeds(1,4)
		nran0=ranseeds(1,5)
		rans=nran0
	end if
	if (length > 1) then
		iran => ranseeds(2:,1)
		jran => ranseeds(2:,2)
		kran => ranseeds(2:,3)
		mran => ranseeds(2:,4)
		nran => ranseeds(2:,5)
		ranv = nran
	end if
	lenran=length
	END SUBROUTINE ran_init
!BL
	SUBROUTINE ran_deallocate
	if (lenran > 0) then
		deallocate(ranseeds,ranv)
		nullify(ranseeds,ranv,iran,jran,kran,mran,nran)
		lenran = 0
	end if
	END SUBROUTINE ran_deallocate
!BL
	SUBROUTINE ran_seed(sequence,size,put,get)
	IMPLICIT NONE
	INTEGER, OPTIONAL, INTENT(IN) :: sequence
	INTEGER, OPTIONAL, INTENT(OUT) :: size
	INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: put
	INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: get
	if (present(size)) then
		size=5*lenran
	else if (present(put)) then
		if (lenran == 0) RETURN
		ranseeds=reshape(put,shape(ranseeds))
		where (ranseeds(:,1:3) < 0) ranseeds(:,1:3)=not(ranseeds(:,1:3))
		where (ranseeds(:,4:5) == 0) ranseeds(:,4:5)=1
		iran0=ranseeds(1,1)
		jran0=ranseeds(1,2)
		kran0=ranseeds(1,3)
		mran0=ranseeds(1,4)
		nran0=ranseeds(1,5)
	else if (present(get)) then
		if (lenran == 0) RETURN
		ranseeds(1,1:5)=(/ iran0,jran0,kran0,mran0,nran0 /)
		get=reshape(ranseeds,shape(get))
	else if (present(sequence)) then
		call ran_deallocate
		seq=sequence
	end if
	END SUBROUTINE ran_seed
!BL
	SUBROUTINE ran_hash_s(il,ir)
	IMPLICIT NONE
	INTEGER(K4B), INTENT(INOUT) :: il,ir
	INTEGER(K4B) :: is,j
	do j=1,4
		is=ir
		ir=ieor(ir,ishft(ir,5))+1422217823
		ir=ieor(ir,ishft(ir,-16))+1842055030
		ir=ieor(ir,ishft(ir,9))+80567781
		ir=ieor(il,ir)
		il=is
	end do
	END SUBROUTINE ran_hash_s
!BL
	SUBROUTINE ran_hash_v(il,ir)
	IMPLICIT NONE
	INTEGER(K4B), DIMENSION(:), INTENT(INOUT) :: il,ir
	INTEGER(K4B), DIMENSION(size(il)) :: is
	INTEGER(K4B) :: j
	do j=1,4
		is=ir
		ir=ieor(ir,ishft(ir,5))+1422217823
		ir=ieor(ir,ishft(ir,-16))+1842055030
		ir=ieor(ir,ishft(ir,9))+80567781
		ir=ieor(il,ir)
		il=is
	end do
	END SUBROUTINE ran_hash_v
END MODULE ran_state

! -----------------------------------------------------------

	FUNCTION rank(index)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: index
	INTEGER(I4B), DIMENSION(size(index)) :: rank
	rank(index(:))=arth(1,1,size(index))
	END FUNCTION rank

! -----------------------------------------------------------

	SUBROUTINE ratint(xa,ya,x,y,dy)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,iminloc,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
	REAL(SP), INTENT(IN) :: x
	REAL(SP), INTENT(OUT) :: y,dy
	INTEGER(I4B) :: m,n,ns
	REAL(SP), DIMENSION(size(xa)) :: c,d,dd,h,t
	REAL(SP), PARAMETER :: TINY=1.0e-25_sp
	n=assert_eq(size(xa),size(ya),'ratint')
	h=xa-x
	ns=iminloc(abs(h))
	y=ya(ns)
	if (x == xa(ns)) then
		dy=0.0
		RETURN
	end if
	c=ya
	d=ya+TINY
	ns=ns-1
	do m=1,n-1
		t(1:n-m)=(xa(1:n-m)-x)*d(1:n-m)/h(1+m:n)
		dd(1:n-m)=t(1:n-m)-c(2:n-m+1)
		if (any(dd(1:n-m) == 0.0)) &
			call nrerror('failure in ratint')
		dd(1:n-m)=(c(2:n-m+1)-d(1:n-m))/dd(1:n-m)
		d(1:n-m)=c(2:n-m+1)*dd(1:n-m)
		c(1:n-m)=t(1:n-m)*dd(1:n-m)
		if (2*ns < n-m) then
			dy=c(ns+1)
		else
			dy=d(ns)
			ns=ns-1
		end if
		y=y+dy
	end do
	END SUBROUTINE ratint

! -----------------------------------------------------------

	SUBROUTINE ratlsq(func,a,b,mm,kk,cof,dev)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,geop
	use mo_nr, ONLY : ratval,svbksb,svdcmp
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,b
	INTEGER(I4B), INTENT(IN) :: mm,kk
	REAL(DP), DIMENSION(:), INTENT(OUT) :: cof
	REAL(DP), INTENT(OUT) :: dev
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: NPFAC=8,MAXIT=5
	REAL(DP), PARAMETER :: BIG=1.0e30_dp
	INTEGER(I4B) :: it,ncof,npt,npth
	REAL(DP) :: devmax,e,theta
	REAL(DP), DIMENSION((mm+kk+1)*NPFAC) :: bb,ee,fs,wt,xs
	REAL(DP), DIMENSION(mm+kk+1) :: coff,w
	REAL(DP), DIMENSION(mm+kk+1,mm+kk+1) :: v
	REAL(DP), DIMENSION((mm+kk+1)*NPFAC,mm+kk+1) :: u,temp
	ncof=mm+kk+1
	npt=NPFAC*ncof
	npth=npt/2
	dev=BIG
	theta=PIO2_D/(npt-1)
	xs(1:npth-1)=a+(b-a)*sin(theta*arth(0,1,npth-1))**2
	xs(npth:npt)=b-(b-a)*sin(theta*arth(npt-npth,-1,npt-npth+1))**2
	fs=func(xs)
	wt=1.0
	ee=1.0
	e=0.0
	do it=1,MAXIT
		bb=wt*(fs+sign(e,ee))
		temp=geop(spread(1.0_dp,1,npt),xs,ncof)
		u(:,1:mm+1)=temp(:,1:mm+1)*spread(wt,2,mm+1)
		u(:,mm+2:ncof)=-temp(:,2:ncof-mm)*spread(bb,2,ncof-mm-1)
		call svdcmp(u,w,v)
		call svbksb(u,w,v,bb,coff)
		ee=ratval(xs,coff,mm,kk)-fs
		wt=abs(ee)
		devmax=maxval(wt)
		e=sum(wt)/npt
		if (devmax <= dev) then
			cof=coff
			dev=devmax
		end if
		write(*,10) it,devmax
	end do
10	format (' ratlsq iteration=',i2,' max error=',1p,e10.3)
	END SUBROUTINE ratlsq

! -----------------------------------------------------------

	FUNCTION ratval_s(x,cof,mm,kk)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	INTEGER(I4B), INTENT(IN) :: mm,kk
	REAL(DP), DIMENSION(mm+kk+1), INTENT(IN) :: cof
	REAL(DP) :: ratval_s
	ratval_s=poly(x,cof(1:mm+1))/(1.0_dp+x*poly(x,cof(mm+2:mm+kk+1)))
	END FUNCTION ratval_s

	FUNCTION ratval_v(x,cof,mm,kk)
	use mo_nrtype
          use mo_nrutil, ONLY : poly
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	INTEGER(I4B), INTENT(IN) :: mm,kk
	REAL(DP), DIMENSION(mm+kk+1), INTENT(IN) :: cof
	REAL(DP), DIMENSION(size(x)) :: ratval_v
	ratval_v=poly(x,cof(1:mm+1))/(1.0_dp+x*poly(x,cof(mm+2:mm+kk+1)))
	END FUNCTION ratval_v

! -----------------------------------------------------------

	FUNCTION rc_s(x,y)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y
	REAL(SP) :: rc_s
	REAL(SP), PARAMETER :: ERRTOL=0.04_sp,TINY=1.69e-38_sp,&
		SQRTNY=1.3e-19_sp,BIG=3.0e37_sp,TNBG=TINY*BIG,&
		COMP1=2.236_sp/SQRTNY,COMP2=TNBG*TNBG/25.0_sp,&
		THIRD=1.0_sp/3.0_sp,&
		C1=0.3_sp,C2=1.0_sp/7.0_sp,C3=0.375_sp,C4=9.0_sp/22.0_sp
	REAL(SP) :: alamb,ave,s,w,xt,yt
	call assert( (/x >= 0.0,y /= 0.0,x+abs(y) >= TINY,x+abs(y) <= BIG, &
		y >= -COMP1 .or. x <= 0.0 .or. x >= COMP2/),'rc_s')
	if (y > 0.0) then
		xt=x
		yt=y
		w=1.0
	else
		xt=x-y
		yt=-y
		w=sqrt(x)/sqrt(xt)
	end if
	do
		alamb=2.0_sp*sqrt(xt)*sqrt(yt)+yt
		xt=0.25_sp*(xt+alamb)
		yt=0.25_sp*(yt+alamb)
		ave=THIRD*(xt+yt+yt)
		s=(yt-ave)/ave
		if (abs(s) <= ERRTOL) exit
	end do
	rc_s=w*(1.0_sp+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
	END FUNCTION rc_s


	FUNCTION rc_v(x,y)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), DIMENSION(size(x)) :: rc_v
	REAL(SP), PARAMETER :: ERRTOL=0.04_sp,TINY=1.69e-38_sp,&
		SQRTNY=1.3e-19_sp,BIG=3.0e37_sp,TNBG=TINY*BIG,&
		COMP1=2.236_sp/SQRTNY,COMP2=TNBG*TNBG/25.0_sp,&
		THIRD=1.0_sp/3.0_sp,&
		C1=0.3_sp,C2=1.0_sp/7.0_sp,C3=0.375_sp,C4=9.0_sp/22.0_sp
	REAL(SP), DIMENSION(size(x)) :: alamb,ave,s,w,xt,yt
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(x),size(y),'rc_v')
	call assert( (/all(x >= 0.0),all(y /= 0.0),all(x+abs(y) >= TINY), &
		all(x+abs(y) <= BIG),all(y >= -COMP1 .or. x <= 0.0  &
		.or. x >= COMP2) /),'rc_v')
	where (y > 0.0)
		xt=x
		yt=y
		w=1.0
	elsewhere
		xt=x-y
		yt=-y
		w=sqrt(x)/sqrt(xt)
	end where
	converged=.false.
	do
		where (.not. converged)
			alamb=2.0_sp*sqrt(xt)*sqrt(yt)+yt
			xt=0.25_sp*(xt+alamb)
			yt=0.25_sp*(yt+alamb)
			ave=THIRD*(xt+yt+yt)
			s=(yt-ave)/ave
			converged = (abs(s) <= ERRTOL)
		end where
		if (all(converged)) exit
	end do
	rc_v=w*(1.0_sp+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
	END FUNCTION rc_v

! -----------------------------------------------------------

	FUNCTION rd_s(x,y,z)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y,z
	REAL(SP) :: rd_s
	REAL(SP), PARAMETER :: ERRTOL=0.05_sp,TINY=1.0e-25_sp,BIG=4.5e21_sp,&
		C1=3.0_sp/14.0_sp,C2=1.0_sp/6.0_sp,C3=9.0_sp/22.0_sp,&
		C4=3.0_sp/26.0_sp,C5=0.25_sp*C3,C6=1.5_sp*C4
	REAL(SP) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
		ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
	call assert(min(x,y) >= 0.0, min(x+y,z) >= TINY, max(x,y,z) <= BIG, &
		'rd_s args')
	xt=x
	yt=y
	zt=z
	sum=0.0
	fac=1.0
	do
		sqrtx=sqrt(xt)
		sqrty=sqrt(yt)
		sqrtz=sqrt(zt)
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
		sum=sum+fac/(sqrtz*(zt+alamb))
		fac=0.25_sp*fac
		xt=0.25_sp*(xt+alamb)
		yt=0.25_sp*(yt+alamb)
		zt=0.25_sp*(zt+alamb)
		ave=0.2_sp*(xt+yt+3.0_sp*zt)
		delx=(ave-xt)/ave
		dely=(ave-yt)/ave
		delz=(ave-zt)/ave
		if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
	end do
	ea=delx*dely
	eb=delz*delz
	ec=ea-eb
	ed=ea-6.0_sp*eb
	ee=ed+ec+ec
	rd_s=3.0_sp*sum+fac*(1.0_sp+ed*(-C1+C5*ed-C6*delz*ee)&
		+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
	END FUNCTION rd_s


	FUNCTION rd_v(x,y,z)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z
	REAL(SP), DIMENSION(size(x)) :: rd_v
	REAL(SP), PARAMETER :: ERRTOL=0.05_sp,TINY=1.0e-25_sp,BIG=4.5e21_sp,&
		C1=3.0_sp/14.0_sp,C2=1.0_sp/6.0_sp,C3=9.0_sp/22.0_sp,&
		C4=3.0_sp/26.0_sp,C5=0.25_sp*C3,C6=1.5_sp*C4
	REAL(SP), DIMENSION(size(x)) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
		ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(x),size(y),size(z),'rd_v')
	call assert(all(min(x,y) >= 0.0), all(min(x+y,z) >= TINY), &
		all(max(x,y,z) <= BIG), 'rd_v args')
	xt=x
	yt=y
	zt=z
	sum=0.0
	fac=1.0
	converged=.false.
	do
		where (.not. converged)
			sqrtx=sqrt(xt)
			sqrty=sqrt(yt)
			sqrtz=sqrt(zt)
			alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
			sum=sum+fac/(sqrtz*(zt+alamb))
			fac=0.25_sp*fac
			xt=0.25_sp*(xt+alamb)
			yt=0.25_sp*(yt+alamb)
			zt=0.25_sp*(zt+alamb)
			ave=0.2_sp*(xt+yt+3.0_sp*zt)
			delx=(ave-xt)/ave
			dely=(ave-yt)/ave
			delz=(ave-zt)/ave
			converged = (all(max(abs(delx),abs(dely),abs(delz)) <= ERRTOL))
		end where
		if (all(converged)) exit
	end do
	ea=delx*dely
	eb=delz*delz
	ec=ea-eb
	ed=ea-6.0_sp*eb
	ee=ed+ec+ec
	rd_v=3.0_sp*sum+fac*(1.0_sp+ed*(-C1+C5*ed-C6*delz*ee)&
		+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
	END FUNCTION rd_v

! -----------------------------------------------------------

	SUBROUTINE realft_sp(data,isign,zdata)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq,zroots_unity
	use mo_nr, ONLY : four1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
	INTEGER(I4B) :: n,ndum,nh,nq
	COMPLEX(SPC), DIMENSION(size(data)/4) :: w
	COMPLEX(SPC), DIMENSION(size(data)/4-1) :: h1,h2
	COMPLEX(SPC), DIMENSION(:), POINTER :: cdata
	COMPLEX(SPC) :: z
	REAL(SP) :: c1=0.5_sp,c2
	n=size(data)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_sp')
	nh=n/2
	nq=n/4
	if (present(zdata)) then
		ndum=assert_eq(n/2,size(zdata),'realft_sp')
		cdata=>zdata
		if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
	else
		allocate(cdata(n/2))
		cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
	end if
	if (isign == 1) then
		c2=-0.5_sp
		call four1(cdata,+1)
	else
		c2=0.5_sp
	end if
	w=zroots_unity(sign(n,isign),n/4)
	w=cmplx(-aimag(w),real(w),kind=spc)
	h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
	h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
	cdata(2:nq)=h1+w(2:nq)*h2
	cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
	z=cdata(1)
	if (isign == 1) then
		cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=spc)
	else
		cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=spc)
		call four1(cdata,-1)
	end if
	if (present(zdata)) then
		if (isign /= 1) then
			data(1:n-1:2)=real(cdata)
			data(2:n:2)=aimag(cdata)
		end if
	else
		data(1:n-1:2)=real(cdata)
		data(2:n:2)=aimag(cdata)
		deallocate(cdata)
	end if
	END SUBROUTINE realft_sp


	SUBROUTINE realft_dp(data,isign,zdata)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq,zroots_unity
	use mo_nr, ONLY : four1
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(DPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
	INTEGER(I4B) :: n,ndum,nh,nq
	COMPLEX(DPC), DIMENSION(size(data)/4) :: w
	COMPLEX(DPC), DIMENSION(size(data)/4-1) :: h1,h2
	COMPLEX(DPC), DIMENSION(:), POINTER :: cdata
	COMPLEX(DPC) :: z
	REAL(DP) :: c1=0.5_dp,c2
	n=size(data)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_dp')
	nh=n/2
	nq=n/4
	if (present(zdata)) then
		ndum=assert_eq(n/2,size(zdata),'realft_dp')
		cdata=>zdata
		if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
	else
		allocate(cdata(n/2))
		cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
	end if
	if (isign == 1) then
		c2=-0.5_dp
		call four1(cdata,+1)
	else
		c2=0.5_dp
	end if
	w=zroots_unity(sign(n,isign),n/4)
	w=cmplx(-aimag(w),real(w),kind=dpc)
	h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
	h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
	cdata(2:nq)=h1+w(2:nq)*h2
	cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
	z=cdata(1)
	if (isign == 1) then
		cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=dpc)
	else
		cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=dpc)
		call four1(cdata,-1)
	end if
	if (present(zdata)) then
		if (isign /= 1) then
			data(1:n-1:2)=real(cdata)
			data(2:n:2)=aimag(cdata)
		end if
	else
		data(1:n-1:2)=real(cdata)
		data(2:n:2)=aimag(cdata)
		deallocate(cdata)
	end if
	END SUBROUTINE realft_dp

! -----------------------------------------------------------

	RECURSIVE FUNCTION recur1(a,b) RESULT(u)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(size(a)) :: u
	INTEGER(I4B), PARAMETER :: NPAR_RECUR1=8
	INTEGER(I4B) :: n,j
	n=assert_eq(size(a),size(b)+1,'recur1')
	u(1)=a(1)
	if (n < NPAR_RECUR1) then
		do j=2,n
			u(j)=a(j)+b(j-1)*u(j-1)
		end do
	else
		u(2:n:2)=recur1(a(2:n:2)+a(1:n-1:2)*b(1:n-1:2), &
				b(3:n-1:2)*b(2:n-2:2))
		u(3:n:2)=a(3:n:2)+b(2:n-1:2)*u(2:n-1:2)
	end if
	END FUNCTION recur1

! -----------------------------------------------------------

	FUNCTION recur2(a,b,c)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c
	REAL(SP), DIMENSION(size(a)) :: recur2
	INTEGER(I4B) :: n
	REAL(SP), DIMENSION(size(a)-1) :: a1,a2,u1,u2
	REAL(SP), DIMENSION(size(a)-2) :: b11,b12,b21,b22
	n=assert_eq(size(a),size(b)+2,size(c)+2,'recur2')
	a1(1)=a(1)
	a2(1)=a(2)
	a1(2:n-1)=0.0
	a2(2:n-1)=a(3:n)
	b11(1:n-2)=0.0
	b12(1:n-2)=1.0
	b21(1:n-2)=c(1:n-2)
	b22(1:n-2)=b(1:n-2)
	call recur1_v(a1,a2,b11,b12,b21,b22,u1,u2)
	recur2(1:n-1)=u1(1:n-1)
	recur2(n)=u2(n-1)
	CONTAINS
!BL
	RECURSIVE SUBROUTINE recur1_v(a1,a2,b11,b12,b21,b22,u1,u2)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a1,a2,b11,b12,b21,b22
	REAL(SP), DIMENSION(:), INTENT(OUT) :: u1,u2
	INTEGER(I4B), PARAMETER :: NPAR_RECUR2=8
	INTEGER(I4B) :: n,j,nn,nn1
	REAL(SP), DIMENSION(size(a1)/2) :: aa1,aa2
	REAL(SP), DIMENSION(size(a1)/2-1) :: bb11,bb12,bb21,bb22
	n=assert_eq((/size(a1),size(a2),size(b11)+1,size(b12)+1,size(b21)+1,&
		size(b22)+1,size(u1),size(u2)/),'recur1_v')
	u1(1)=a1(1)
	u2(1)=a2(1)
	if (n < NPAR_RECUR2) then
		do j=2,n
			u1(j)=a1(j)+b11(j-1)*u1(j-1)+b12(j-1)*u2(j-1)
			u2(j)=a2(j)+b21(j-1)*u1(j-1)+b22(j-1)*u2(j-1)
		end do
	else
		nn=n/2
		nn1=nn-1
		aa1(1:nn)=a1(2:n:2)+b11(1:n-1:2)*a1(1:n-1:2)+&
			b12(1:n-1:2)*a2(1:n-1:2)
		aa2(1:nn)=a2(2:n:2)+b21(1:n-1:2)*a1(1:n-1:2)+&
				b22(1:n-1:2)*a2(1:n-1:2)
		bb11(1:nn1)=b11(3:n-1:2)*b11(2:n-2:2)+&
				b12(3:n-1:2)*b21(2:n-2:2)
		bb12(1:nn1)=b11(3:n-1:2)*b12(2:n-2:2)+&
				b12(3:n-1:2)*b22(2:n-2:2)
		bb21(1:nn1)=b21(3:n-1:2)*b11(2:n-2:2)+&
				b22(3:n-1:2)*b21(2:n-2:2)
		bb22(1:nn1)=b21(3:n-1:2)*b12(2:n-2:2)+&
				b22(3:n-1:2)*b22(2:n-2:2)
		call recur1_v(aa1,aa2,bb11,bb12,bb21,bb22,u1(2:n:2),u2(2:n:2))
		u1(3:n:2)=a1(3:n:2)+b11(2:n-1:2)*u1(2:n-1:2)+&
				b12(2:n-1:2)*u2(2:n-1:2)
		u2(3:n:2)=a2(3:n:2)+b21(2:n-1:2)*u1(2:n-1:2)+&
				b22(2:n-1:2)*u2(2:n-1:2)
	end if
	END SUBROUTINE recur1_v
	END FUNCTION recur2

! -----------------------------------------------------------

	SUBROUTINE relax(u,rhs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
	INTEGER(I4B) :: n
	REAL(DP) :: h,h2
	n=assert_eq(size(u,1),size(u,2),size(rhs,1),size(rhs,2),'relax')
	h=1.0_dp/(n-1)
	h2=h*h
	u(2:n-1:2,2:n-1:2)=0.25_dp*(u(3:n:2,2:n-1:2)+u(1:n-2:2,2:n-1:2)+&
		u(2:n-1:2,3:n:2)+u(2:n-1:2,1:n-2:2)-h2*rhs(2:n-1:2,2:n-1:2))
	u(3:n-2:2,3:n-2:2)=0.25_dp*(u(4:n-1:2,3:n-2:2)+u(2:n-3:2,3:n-2:2)+&
		u(3:n-2:2,4:n-1:2)+u(3:n-2:2,2:n-3:2)-h2*rhs(3:n-2:2,3:n-2:2))
	u(3:n-2:2,2:n-1:2)=0.25_dp*(u(4:n-1:2,2:n-1:2)+u(2:n-3:2,2:n-1:2)+&
		u(3:n-2:2,3:n:2)+u(3:n-2:2,1:n-2:2)-h2*rhs(3:n-2:2,2:n-1:2))
	u(2:n-1:2,3:n-2:2)=0.25_dp*(u(3:n:2,3:n-2:2)+u(1:n-2:2,3:n-2:2)+&
		u(2:n-1:2,4:n-1:2)+u(2:n-1:2,2:n-3:2)-h2*rhs(2:n-1:2,3:n-2:2))
	END SUBROUTINE relax

! -----------------------------------------------------------

	SUBROUTINE relax2(u,rhs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: rhs
	INTEGER(I4B) :: n
	REAL(DP) :: foh2,h,h2i
	REAL(DP) :: res(size(u,1),size(u,1))
	n=assert_eq(size(u,1),size(u,2),size(rhs,1),size(rhs,2),'relax2')
	h=1.0_dp/(n-1)
	h2i=1.0_dp/(h*h)
	foh2=-4.0_dp*h2i
	res(2:n-1:2,2:n-1:2)=h2i*(u(3:n:2,2:n-1:2)+u(1:n-2:2,2:n-1:2)+&
		u(2:n-1:2,3:n:2)+u(2:n-1:2,1:n-2:2)-4.0_dp*u(2:n-1:2,2:n-1:2))&
		+u(2:n-1:2,2:n-1:2)**2-rhs(2:n-1:2,2:n-1:2)
	u(2:n-1:2,2:n-1:2)=u(2:n-1:2,2:n-1:2)-res(2:n-1:2,2:n-1:2)/&
		(foh2+2.0_dp*u(2:n-1:2,2:n-1:2))
	res(3:n-2:2,3:n-2:2)=h2i*(u(4:n-1:2,3:n-2:2)+u(2:n-3:2,3:n-2:2)+&
		u(3:n-2:2,4:n-1:2)+u(3:n-2:2,2:n-3:2)-4.0_dp*u(3:n-2:2,3:n-2:2))&
		+u(3:n-2:2,3:n-2:2)**2-rhs(3:n-2:2,3:n-2:2)
	u(3:n-2:2,3:n-2:2)=u(3:n-2:2,3:n-2:2)-res(3:n-2:2,3:n-2:2)/&
		(foh2+2.0_dp*u(3:n-2:2,3:n-2:2))
	res(3:n-2:2,2:n-1:2)=h2i*(u(4:n-1:2,2:n-1:2)+u(2:n-3:2,2:n-1:2)+&
		u(3:n-2:2,3:n:2)+u(3:n-2:2,1:n-2:2)-4.0_dp*u(3:n-2:2,2:n-1:2))&
		+u(3:n-2:2,2:n-1:2)**2-rhs(3:n-2:2,2:n-1:2)
	u(3:n-2:2,2:n-1:2)=u(3:n-2:2,2:n-1:2)-res(3:n-2:2,2:n-1:2)/&
		(foh2+2.0_dp*u(3:n-2:2,2:n-1:2))
	res(2:n-1:2,3:n-2:2)=h2i*(u(3:n:2,3:n-2:2)+u(1:n-2:2,3:n-2:2)+&
		u(2:n-1:2,4:n-1:2)+u(2:n-1:2,2:n-3:2)-4.0_dp*u(2:n-1:2,3:n-2:2))&
		+u(2:n-1:2,3:n-2:2)**2-rhs(2:n-1:2,3:n-2:2)
	u(2:n-1:2,3:n-2:2)=u(2:n-1:2,3:n-2:2)-res(2:n-1:2,3:n-2:2)/&
		(foh2+2.0_dp*u(2:n-1:2,3:n-2:2))
	END SUBROUTINE relax2

! -----------------------------------------------------------

	FUNCTION resid(u,rhs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,rhs
	REAL(DP), DIMENSION(size(u,1),size(u,1)) :: resid
	INTEGER(I4B) :: n
	REAL(DP) :: h,h2i
	n=assert_eq((/size(u,1),size(u,2),size(rhs,1),size(rhs,2)/),'resid')
	n=size(u,1)
	h=1.0_dp/(n-1)
	h2i=1.0_dp/(h*h)
	resid(2:n-1,2:n-1)=-h2i*(u(3:n,2:n-1)+u(1:n-2,2:n-1)+u(2:n-1,3:n)+&
		u(2:n-1,1:n-2)-4.0_dp*u(2:n-1,2:n-1))+rhs(2:n-1,2:n-1)
	resid(1:n,1)=0.0
	resid(1:n,n)=0.0
	resid(1,1:n)=0.0
	resid(n,1:n)=0.0
	END FUNCTION resid

! -----------------------------------------------------------

	FUNCTION rf_s(x,y,z)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y,z
	REAL(SP) :: rf_s
	REAL(SP), PARAMETER :: ERRTOL=0.08_sp,TINY=1.5e-38_sp,BIG=3.0e37_sp,&
		THIRD=1.0_sp/3.0_sp,&
		C1=1.0_sp/24.0_sp,C2=0.1_sp,C3=3.0_sp/44.0_sp,C4=1.0_sp/14.0_sp
	REAL(SP) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
	call assert(min(x,y,z) >= 0.0, min(x+y,x+z,y+z) >= TINY, &
		max(x,y,z) <= BIG, 'rf_s args')
	xt=x
	yt=y
	zt=z
	do
		sqrtx=sqrt(xt)
		sqrty=sqrt(yt)
		sqrtz=sqrt(zt)
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
		xt=0.25_sp*(xt+alamb)
		yt=0.25_sp*(yt+alamb)
		zt=0.25_sp*(zt+alamb)
		ave=THIRD*(xt+yt+zt)
		delx=(ave-xt)/ave
		dely=(ave-yt)/ave
		delz=(ave-zt)/ave
		if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
	end do
	e2=delx*dely-delz**2
	e3=delx*dely*delz
	rf_s=(1.0_sp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
	END FUNCTION rf_s


	FUNCTION rf_v(x,y,z)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z
	REAL(SP), DIMENSION(size(x)) :: rf_v
	REAL(SP), PARAMETER :: ERRTOL=0.08_sp,TINY=1.5e-38_sp,BIG=3.0e37_sp,&
		THIRD=1.0_sp/3.0_sp,&
		C1=1.0_sp/24.0_sp,C2=0.1_sp,C3=3.0_sp/44.0_sp,C4=1.0_sp/14.0_sp
	REAL(SP), DIMENSION(size(x)) :: alamb,ave,delx,dely,delz,e2,e3,&
		sqrtx,sqrty,sqrtz,xt,yt,zt
	LOGICAL(LGT), DIMENSION(size(x)) :: converged
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(x),size(y),size(z),'rf_v')
	call assert(all(min(x,y,z) >= 0.0), all(min(x+y,x+z,y+z) >= TINY), &
		all(max(x,y,z) <= BIG), 'rf_v args')
	xt=x
	yt=y
	zt=z
	converged=.false.
	do
		where (.not. converged)
			sqrtx=sqrt(xt)
			sqrty=sqrt(yt)
			sqrtz=sqrt(zt)
			alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
			xt=0.25_sp*(xt+alamb)
			yt=0.25_sp*(yt+alamb)
			zt=0.25_sp*(zt+alamb)
			ave=THIRD*(xt+yt+zt)
			delx=(ave-xt)/ave
			dely=(ave-yt)/ave
			delz=(ave-zt)/ave
			converged = (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL)
		end where
		if (all(converged)) exit
	end do
	e2=delx*dely-delz**2
	e3=delx*dely*delz
	rf_v=(1.0_sp+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
	END FUNCTION rf_v

! -----------------------------------------------------------

	FUNCTION rj_s(x,y,z,p)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : rc,rf
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y,z,p
	REAL(SP) :: rj_s
	REAL(SP), PARAMETER :: ERRTOL=0.05_sp,TINY=2.5e-13_sp,BIG=9.0e11_sp,&
		C1=3.0_sp/14.0_sp,C2=1.0_sp/3.0_sp,C3=3.0_sp/22.0_sp,&
		C4=3.0_sp/26.0_sp,C5=0.75_sp*C3,C6=1.5_sp*C4,C7=0.5_sp*C2,&
		C8=C3+C3
	REAL(SP) :: a,alamb,alpha,ave,b,bet,delp,delx,&
		dely,delz,ea,eb,ec,ed,ee,fac,pt,rho,sqrtx,sqrty,sqrtz,&
		sm,tau,xt,yt,zt
	call assert(min(x,y,z) >= 0.0, min(x+y,x+z,y+z,abs(p)) >= TINY, &
		max(x,y,z,abs(p)) <= BIG, 'rj_s args')
	sm=0.0
	fac=1.0
	if (p > 0.0) then
		xt=x
		yt=y
		zt=z
		pt=p
	else
		xt=min(x,y,z)
		zt=max(x,y,z)
		yt=x+y+z-xt-zt
		a=1.0_sp/(yt-p)
		b=a*(zt-yt)*(yt-xt)
		pt=yt+b
		rho=xt*zt/yt
		tau=p*pt/yt
	end if
	do
		sqrtx=sqrt(xt)
		sqrty=sqrt(yt)
		sqrtz=sqrt(zt)
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
		alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
		bet=pt*(pt+alamb)**2
		sm=sm+fac*rc(alpha,bet)
		fac=0.25_sp*fac
		xt=0.25_sp*(xt+alamb)
		yt=0.25_sp*(yt+alamb)
		zt=0.25_sp*(zt+alamb)
		pt=0.25_sp*(pt+alamb)
		ave=0.2_sp*(xt+yt+zt+pt+pt)
		delx=(ave-xt)/ave
		dely=(ave-yt)/ave
		delz=(ave-zt)/ave
		delp=(ave-pt)/ave
		if (max(abs(delx),abs(dely),abs(delz),abs(delp)) <= ERRTOL) exit
	end do
	ea=delx*(dely+delz)+dely*delz
	eb=delx*dely*delz
	ec=delp**2
	ed=ea-3.0_sp*ec
	ee=eb+2.0_sp*delp*(ea-ec)
	rj_s=3.0_sp*sm+fac*(1.0_sp+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8&
		+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
	if (p <= 0.0) rj_s=a*(b*rj_s+3.0_sp*(rc(rho,tau)-rf(xt,yt,zt)))
	END FUNCTION rj_s


	FUNCTION rj_v(x,y,z,p)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq
	use mo_nr, ONLY : rc,rf
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,z,p
	REAL(SP), DIMENSION(size(x)) :: rj_v
	REAL(SP), PARAMETER :: ERRTOL=0.05_sp,TINY=2.5e-13_sp,BIG=9.0e11_sp,&
		C1=3.0_sp/14.0_sp,C2=1.0_sp/3.0_sp,C3=3.0_sp/22.0_sp,&
		C4=3.0_sp/26.0_sp,C5=0.75_sp*C3,C6=1.5_sp*C4,C7=0.5_sp*C2,&
		C8=C3+C3
	REAL(SP), DIMENSION(size(x)) :: a,alamb,alpha,ave,b,bet,delp,delx,&
		dely,delz,ea,eb,ec,ed,ee,fac,pt,rho,sqrtx,sqrty,sqrtz,&
		sm,tau,xt,yt,zt
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(x),size(y),size(z),size(p),'rj_v')
	call assert(all(min(x,y,z) >= 0.0), all(min(x+y,x+z,y+z,abs(p)) >= TINY), &
		all(max(x,y,z,abs(p)) <= BIG), 'rj_v args')
	sm=0.0
	fac=1.0
	where (p > 0.0)
		xt=x
		yt=y
		zt=z
		pt=p
	elsewhere
		xt=min(x,y,z)
		zt=max(x,y,z)
		yt=x+y+z-xt-zt
		a=1.0_sp/(yt-p)
		b=a*(zt-yt)*(yt-xt)
		pt=yt+b
		rho=xt*zt/yt
		tau=p*pt/yt
	end where
	mask=.false.
	do
		where (.not. mask)
			sqrtx=sqrt(xt)
			sqrty=sqrt(yt)
			sqrtz=sqrt(zt)
			alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
			alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
			bet=pt*(pt+alamb)**2
			sm=sm+fac*rc(alpha,bet)
			fac=0.25_sp*fac
			xt=0.25_sp*(xt+alamb)
			yt=0.25_sp*(yt+alamb)
			zt=0.25_sp*(zt+alamb)
			pt=0.25_sp*(pt+alamb)
			ave=0.2_sp*(xt+yt+zt+pt+pt)
			delx=(ave-xt)/ave
			dely=(ave-yt)/ave
			delz=(ave-zt)/ave
			delp=(ave-pt)/ave
			mask = (max(abs(delx),abs(dely),abs(delz),abs(delp)) <= ERRTOL)
		end where
		if (all(mask)) exit
	end do
	ea=delx*(dely+delz)+dely*delz
	eb=delx*dely*delz
	ec=delp**2
	ed=ea-3.0_sp*ec
	ee=eb+2.0_sp*delp*(ea-ec)
	rj_v=3.0_sp*sm+fac*(1.0_sp+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8&
		+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
	mask = (p <= 0.0)
	where (mask) rj_v=a*(b*rj_v+&
		unpack(3.0_sp*(rc(pack(rho,mask),pack(tau,mask))-&
		rf(pack(xt,mask),pack(yt,mask),pack(zt,mask))),mask,0.0_sp))
	END FUNCTION rj_v

! -----------------------------------------------------------

	SUBROUTINE rk4(y,dydx,x,h,yout,derivs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
	REAL(SP), INTENT(IN) :: x,h
	REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: ndum
	REAL(SP) :: h6,hh,xh
	REAL(SP), DIMENSION(size(y)) :: dym,dyt,yt
	ndum=assert_eq(size(y),size(dydx),size(yout),'rk4')
	hh=h*0.5_sp
	h6=h/6.0_sp
	xh=x+hh
	yt=y+hh*dydx
	call derivs(xh,yt,dyt)
	yt=y+hh*dyt
	call derivs(xh,yt,dym)
	yt=y+h*dym
	dym=dyt+dym
	call derivs(x+h,yt,dyt)
	yout=y+h6*(dydx+dyt+2.0_sp*dym)
	END SUBROUTINE rk4

! -----------------------------------------------------------

	SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
	REAL(SP), INTENT(IN) :: x,h
	REAL(SP), DIMENSION(:), INTENT(OUT) :: yout,yerr
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: ndum
	REAL(SP), DIMENSION(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
	REAL(SP), PARAMETER :: A2=0.2_sp,A3=0.3_sp,A4=0.6_sp,A5=1.0_sp,&
		A6=0.875_sp,B21=0.2_sp,B31=3.0_sp/40.0_sp,B32=9.0_sp/40.0_sp,&
		B41=0.3_sp,B42=-0.9_sp,B43=1.2_sp,B51=-11.0_sp/54.0_sp,&
		B52=2.5_sp,B53=-70.0_sp/27.0_sp,B54=35.0_sp/27.0_sp,&
		B61=1631.0_sp/55296.0_sp,B62=175.0_sp/512.0_sp,&
		B63=575.0_sp/13824.0_sp,B64=44275.0_sp/110592.0_sp,&
		B65=253.0_sp/4096.0_sp,C1=37.0_sp/378.0_sp,&
		C3=250.0_sp/621.0_sp,C4=125.0_sp/594.0_sp,&
		C6=512.0_sp/1771.0_sp,DC1=C1-2825.0_sp/27648.0_sp,&
		DC3=C3-18575.0_sp/48384.0_sp,DC4=C4-13525.0_sp/55296.0_sp,&
		DC5=-277.0_sp/14336.0_sp,DC6=C6-0.25_sp
	ndum=assert_eq(size(y),size(dydx),size(yout),size(yerr),'rkck')
	ytemp=y+B21*h*dydx
	call derivs(x+A2*h,ytemp,ak2)
	ytemp=y+h*(B31*dydx+B32*ak2)
	call derivs(x+A3*h,ytemp,ak3)
	ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
	call derivs(x+A4*h,ytemp,ak4)
	ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
	call derivs(x+A5*h,ytemp,ak5)
	ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
	call derivs(x+A6*h,ytemp,ak6)
	yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
	yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
	END SUBROUTINE rkck

! -----------------------------------------------------------

MODULE rkdumb_path
	use mo_nrtype
	REAL(SP), DIMENSION(:), ALLOCATABLE:: xx
	REAL(SP), DIMENSION(:,:), ALLOCATABLE :: y
END MODULE rkdumb_path

	SUBROUTINE rkdumb(vstart,x1,x2,nstep,derivs)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	use mo_nr, ONLY : rk4
	USE rkdumb_path
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: vstart
	REAL(SP), INTENT(IN) :: x1,x2
	INTEGER(I4B), INTENT(IN) :: nstep
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: k
	REAL(SP) :: h,x
	REAL(SP), DIMENSION(size(vstart)) :: dv,v
	v(:)=vstart(:)
	if (allocated(xx)) deallocate(xx)
	if (allocated(y)) deallocate(y)
	allocate(xx(nstep+1))
	allocate(y(size(vstart),nstep+1))
	y(:,1)=v(:)
	xx(1)=x1
	x=x1
	h=(x2-x1)/nstep
	do k=1,nstep
		call derivs(x,v,dv)
		call rk4(v,dv,x,h,v,derivs)
		if (x+h == x) call nrerror('stepsize not significant in rkdumb')
		x=x+h
		xx(k+1)=x
		y(:,k+1)=v(:)
	end do
	END SUBROUTINE rkdumb

! -----------------------------------------------------------

	SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : rkck
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
	REAL(SP), INTENT(INOUT) :: x
	REAL(SP), INTENT(IN) :: htry,eps
	REAL(SP), INTENT(OUT) :: hdid,hnext
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: ndum
	REAL(SP) :: errmax,h,htemp,xnew
	REAL(SP), DIMENSION(size(y)) :: yerr,ytemp
	REAL(SP), PARAMETER :: SAFETY=0.9_sp,PGROW=-0.2_sp,PSHRNK=-0.25_sp,&
		ERRCON=1.89e-4
	ndum=assert_eq(size(y),size(dydx),size(yscal),'rkqs')
	h=htry
	do
		call rkck(y,dydx,x,h,ytemp,yerr,derivs)
		errmax=maxval(abs(yerr(:)/yscal(:)))/eps
		if (errmax <= 1.0) exit
		htemp=SAFETY*h*(errmax**PSHRNK)
		h=sign(max(abs(htemp),0.1_sp*abs(h)),h)
		xnew=x+h
		if (xnew == x) call nrerror('stepsize underflow in rkqs')
	end do
	if (errmax > ERRCON) then
		hnext=SAFETY*h*(errmax**PGROW)
	else
		hnext=5.0_sp*h
	end if
	hdid=h
	x=x+h
	y(:)=ytemp(:)
	END SUBROUTINE rkqs

! -----------------------------------------------------------

	SUBROUTINE rlft2(data,spec,speq,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq
	use mo_nr, ONLY : four2
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: data
	COMPLEX(SPC), DIMENSION(:,:), INTENT(OUT) :: spec
	COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: speq
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER :: i1,j1,nn1,nn2
	REAL(DP) :: theta
	COMPLEX(SPC) :: c1=(0.5_sp,0.0_sp),c2,h1,h2,w
	COMPLEX(SPC), DIMENSION(size(data,2)-1) :: h1a,h2a
	COMPLEX(DPC) :: ww,wp
	nn1=assert_eq(size(data,1),2*size(spec,1),'rlft2: nn1')
	nn2=assert_eq(size(data,2),size(spec,2),size(speq),'rlft2: nn2')
	call assert(iand((/nn1,nn2/),(/nn1,nn2/)-1)==0, &
		'dimensions must be powers of 2 in rlft2')
	c2=cmplx(0.0_sp,-0.5_sp*isign,kind=spc)
	theta=TWOPI_D/(isign*nn1)
	wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=spc)
	if (isign == 1) then
		spec(:,:)=cmplx(data(1:nn1:2,:),data(2:nn1:2,:),kind=spc)
		call four2(spec,isign)
		speq=spec(1,:)
	end if
	h1=c1*(spec(1,1)+conjg(speq(1)))
	h1a=c1*(spec(1,2:nn2)+conjg(speq(nn2:2:-1)))
	h2=c2*(spec(1,1)-conjg(speq(1)))
	h2a=c2*(spec(1,2:nn2)-conjg(speq(nn2:2:-1)))
	spec(1,1)=h1+h2
	spec(1,2:nn2)=h1a+h2a
	speq(1)=conjg(h1-h2)
	speq(nn2:2:-1)=conjg(h1a-h2a)
	ww=cmplx(1.0_dp,0.0_dp,kind=dpc)
	do i1=2,nn1/4+1
		j1=nn1/2-i1+2
		ww=ww*wp+ww
		w=ww
		h1=c1*(spec(i1,1)+conjg(spec(j1,1)))
		h1a=c1*(spec(i1,2:nn2)+conjg(spec(j1,nn2:2:-1)))
		h2=c2*(spec(i1,1)-conjg(spec(j1,1)))
		h2a=c2*(spec(i1,2:nn2)-conjg(spec(j1,nn2:2:-1)))
		spec(i1,1)=h1+w*h2
		spec(i1,2:nn2)=h1a+w*h2a
		spec(j1,1)=conjg(h1-w*h2)
		spec(j1,nn2:2:-1)=conjg(h1a-w*h2a)
	end do
	if (isign == -1) then
		call four2(spec,isign)
		data(1:nn1:2,:)=real(spec)
		data(2:nn1:2,:)=aimag(spec)
	end if
	END SUBROUTINE rlft2

! -----------------------------------------------------------

	SUBROUTINE rlft3(data,spec,speq,isign)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq
	use mo_nr, ONLY : four3
	REAL(SP), DIMENSION(:,:,:), INTENT(INOUT) :: data
	COMPLEX(SPC), DIMENSION(:,:,:), INTENT(OUT) :: spec
	COMPLEX(SPC), DIMENSION(:,:), INTENT(OUT) :: speq
	INTEGER(I4B), INTENT(IN) :: isign
	INTEGER :: i1,i3,j1,j3,nn1,nn2,nn3
	REAL(DP) :: theta
	COMPLEX(SPC) :: c1=(0.5_sp,0.0_sp),c2,h1,h2,w
	COMPLEX(SPC), DIMENSION(size(data,2)-1) :: h1a,h2a
	COMPLEX(DPC) :: ww,wp
	c2=cmplx(0.0_sp,-0.5_sp*isign,kind=spc)
	nn1=assert_eq(size(data,1),2*size(spec,1),'rlft2: nn1')
	nn2=assert_eq(size(data,2),size(spec,2),size(speq,1),'rlft2: nn2')
	nn3=assert_eq(size(data,3),size(spec,3),size(speq,2),'rlft2: nn3')
	call assert(iand((/nn1,nn2,nn3/),(/nn1,nn2,nn3/)-1)==0, &
		'dimensions must be powers of 2 in rlft3')
	theta=TWOPI_D/(isign*nn1)
	wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
	if (isign == 1) then
		spec(:,:,:)=cmplx(data(1:nn1:2,:,:),data(2:nn1:2,:,:),kind=spc)
		call four3(spec,isign)
		speq=spec(1,:,:)
	end if
	do i3=1,nn3
		j3=1
		if (i3 /= 1) j3=nn3-i3+2
		h1=c1*(spec(1,1,i3)+conjg(speq(1,j3)))
		h1a=c1*(spec(1,2:nn2,i3)+conjg(speq(nn2:2:-1,j3)))
		h2=c2*(spec(1,1,i3)-conjg(speq(1,j3)))
		h2a=c2*(spec(1,2:nn2,i3)-conjg(speq(nn2:2:-1,j3)))
		spec(1,1,i3)=h1+h2
		spec(1,2:nn2,i3)=h1a+h2a
		speq(1,j3)=conjg(h1-h2)
		speq(nn2:2:-1,j3)=conjg(h1a-h2a)
		ww=cmplx(1.0_dp,0.0_dp,kind=dpc)
		do i1=2,nn1/4+1
			j1=nn1/2-i1+2
			ww=ww*wp+ww
			w=ww
			h1=c1*(spec(i1,1,i3)+conjg(spec(j1,1,j3)))
			h1a=c1*(spec(i1,2:nn2,i3)+conjg(spec(j1,nn2:2:-1,j3)))
			h2=c2*(spec(i1,1,i3)-conjg(spec(j1,1,j3)))
			h2a=c2*(spec(i1,2:nn2,i3)-conjg(spec(j1,nn2:2:-1,j3)))
			spec(i1,1,i3)=h1+w*h2
			spec(i1,2:nn2,i3)=h1a+w*h2a
			spec(j1,1,j3)=conjg(h1-w*h2)
			spec(j1,nn2:2:-1,j3)=conjg(h1a-w*h2a)
		end do
	end do
	if (isign == -1) then
		call four3(spec,isign)
		data(1:nn1:2,:,:)=real(spec)
		data(2:nn1:2,:,:)=aimag(spec)
	end if
	END SUBROUTINE rlft3

! -----------------------------------------------------------

	SUBROUTINE rotate(r,qt,i,a,b)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), TARGET, INTENT(INOUT) :: r,qt
	INTEGER(I4B), INTENT(IN) :: i
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(size(r,1)) :: temp
	INTEGER(I4B) :: n
	REAL(SP) :: c,fact,s
	n=assert_eq(size(r,1),size(r,2),size(qt,1),size(qt,2),'rotate')
	if (a == 0.0) then
		c=0.0
		s=sign(1.0_sp,b)
	else if (abs(a) > abs(b)) then
		fact=b/a
		c=sign(1.0_sp/sqrt(1.0_sp+fact**2),a)
		s=fact*c
	else
		fact=a/b
		s=sign(1.0_sp/sqrt(1.0_sp+fact**2),b)
		c=fact*s
	end if
	temp(i:n)=r(i,i:n)
	r(i,i:n)=c*temp(i:n)-s*r(i+1,i:n)
	r(i+1,i:n)=s*temp(i:n)+c*r(i+1,i:n)
	temp=qt(i,:)
	qt(i,:)=c*temp-s*qt(i+1,:)
	qt(i+1,:)=s*temp+c*qt(i+1,:)
	END SUBROUTINE rotate

! -----------------------------------------------------------

	SUBROUTINE rsolv(a,d,b)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	REAL(SP), DIMENSION(:), INTENT(IN) :: d
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: i,n
	n=assert_eq(size(a,1),size(a,2),size(b),size(d),'rsolv')
	b(n)=b(n)/d(n)
	do i=n-1,1,-1
		b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/d(i)
	end do
	END SUBROUTINE rsolv

! -----------------------------------------------------------

	FUNCTION rstrct(uf)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: uf
	REAL(DP), DIMENSION((size(uf,1)+1)/2,(size(uf,1)+1)/2) :: rstrct
	INTEGER(I4B) :: nc,nf
	nf=assert_eq(size(uf,1),size(uf,2),'rstrct')
	nc=(nf+1)/2
	rstrct(2:nc-1,2:nc-1)=0.5_dp*uf(3:nf-2:2,3:nf-2:2)+0.125_dp*(&
		uf(4:nf-1:2,3:nf-2:2)+uf(2:nf-3:2,3:nf-2:2)+&
		uf(3:nf-2:2,4:nf-1:2)+uf(3:nf-2:2,2:nf-3:2))
	rstrct(1:nc,1)=uf(1:nf:2,1)
	rstrct(1:nc,nc)=uf(1:nf:2,nf)
	rstrct(1,1:nc)=uf(1,1:nf:2)
	rstrct(nc,1:nc)=uf(nf,1:nf:2)
	END FUNCTION rstrct

! -----------------------------------------------------------

	FUNCTION rtbis(func,x1,x2,xacc)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,xacc
	REAL(SP) :: rtbis
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXIT=40
	INTEGER(I4B) :: j
	REAL(SP) :: dx,f,fmid,xmid
	fmid=func(x2)
	f=func(x1)
	if (f*fmid >= 0.0) call nrerror('rtbis: root must be bracketed')
	if (f < 0.0) then
		rtbis=x1
		dx=x2-x1
	else
		rtbis=x2
		dx=x1-x2
	end if
	do j=1,MAXIT
		dx=dx*0.5_sp
		xmid=rtbis+dx
		fmid=func(xmid)
		if (fmid <= 0.0) rtbis=xmid
		if (abs(dx) < xacc .or. fmid == 0.0) RETURN
	end do
	call nrerror('rtbis: too many bisections')
	END FUNCTION rtbis

! -----------------------------------------------------------

	FUNCTION rtflsp(func,x1,x2,xacc)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror,swap
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,xacc
	REAL(SP) :: rtflsp
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXIT=30
	INTEGER(I4B) :: j
	REAL(SP) :: del,dx,f,fh,fl,xh,xl
	fl=func(x1)
	fh=func(x2)
	if ((fl > 0.0 .and. fh > 0.0) .or. &
		(fl < 0.0 .and. fh < 0.0)) call &
		nrerror('rtflsp: root must be bracketed between arguments')
	if (fl < 0.0) then
		xl=x1
		xh=x2
	else
		xl=x2
		xh=x1
		call swap(fl,fh)
	end if
	dx=xh-xl
	do j=1,MAXIT
		rtflsp=xl+dx*fl/(fl-fh)
		f=func(rtflsp)
		if (f < 0.0) then
			del=xl-rtflsp
			xl=rtflsp
			fl=f
		else
			del=xh-rtflsp
			xh=rtflsp
			fh=f
		end if
		dx=xh-xl
		if (abs(del) < xacc .or. f == 0.0) RETURN
	end do
	call nrerror('rtflsp exceed maximum iterations')
	END FUNCTION rtflsp

! -----------------------------------------------------------

	FUNCTION rtnewt(funcd,x1,x2,xacc)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,xacc
	REAL(SP) :: rtnewt
	INTERFACE
		SUBROUTINE funcd(x,fval,fderiv)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: fval,fderiv
		END SUBROUTINE funcd
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXIT=20
	INTEGER(I4B) :: j
	REAL(SP) :: df,dx,f
	rtnewt=0.5_sp*(x1+x2)
	do j=1,MAXIT
		call funcd(rtnewt,f,df)
		dx=f/df
		rtnewt=rtnewt-dx
		if ((x1-rtnewt)*(rtnewt-x2) < 0.0)&
			call nrerror('rtnewt: values jumped out of brackets')
		if (abs(dx) < xacc) RETURN
	end do
	call nrerror('rtnewt exceeded maximum iterations')
	END FUNCTION rtnewt

! -----------------------------------------------------------

	FUNCTION rtsafe(funcd,x1,x2,xacc)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,xacc
	REAL(SP) :: rtsafe
	INTERFACE
		SUBROUTINE funcd(x,fval,fderiv)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), INTENT(OUT) :: fval,fderiv
		END SUBROUTINE funcd
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXIT=100
	INTEGER(I4B) :: j
	REAL(SP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
	call funcd(x1,fl,df)
	call funcd(x2,fh,df)
	if ((fl > 0.0 .and. fh > 0.0) .or. &
		(fl < 0.0 .and. fh < 0.0)) &
		call nrerror('root must be bracketed in rtsafe')
	if (fl == 0.0) then
		rtsafe=x1
		RETURN
	else if (fh == 0.0) then
		rtsafe=x2
		RETURN
	else if (fl < 0.0) then
		xl=x1
		xh=x2
	else
		xh=x1
		xl=x2
	end if
	rtsafe=0.5_sp*(x1+x2)
	dxold=abs(x2-x1)
	dx=dxold
	call funcd(rtsafe,f,df)
	do j=1,MAXIT
		if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) >= 0.0 .or. &
			abs(2.0_sp*f) > abs(dxold*df) ) then
			dxold=dx
			dx=0.5_sp*(xh-xl)
			rtsafe=xl+dx
			if (xl == rtsafe) RETURN
		else
			dxold=dx
			dx=f/df
			temp=rtsafe
			rtsafe=rtsafe-dx
			if (temp == rtsafe) RETURN
		end if
		if (abs(dx) < xacc) RETURN
		call funcd(rtsafe,f,df)
		if (f < 0.0) then
			xl=rtsafe
		else
			xh=rtsafe
		end if
	end do
	call nrerror('rtsafe: exceeded maximum iterations')
	END FUNCTION rtsafe

! -----------------------------------------------------------

	FUNCTION rtsec(func,x1,x2,xacc)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror,swap
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,xacc
	REAL(SP) :: rtsec
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXIT=30
	INTEGER(I4B) :: j
	REAL(SP) :: dx,f,fl,xl
	fl=func(x1)
	f=func(x2)
	if (abs(fl) < abs(f)) then
		rtsec=x1
		xl=x2
		call swap(fl,f)
	else
		xl=x1
		rtsec=x2
	end if
	do j=1,MAXIT
		dx=(xl-rtsec)*f/(f-fl)
		xl=rtsec
		fl=f
		rtsec=rtsec+dx
		f=func(rtsec)
		if (abs(dx) < xacc .or. f == 0.0) RETURN
	end do
	call nrerror('rtsec: exceed maximum iterations')
	END FUNCTION rtsec

! -----------------------------------------------------------

	SUBROUTINE rzextr(iest,xest,yest,yz,dy)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: iest
	REAL(SP), INTENT(IN) :: xest
	REAL(SP), DIMENSION(:), INTENT(IN) :: yest
	REAL(SP), DIMENSION(:), INTENT(OUT) :: yz,dy
	INTEGER(I4B), PARAMETER :: IEST_MAX=16
	INTEGER(I4B) :: k,nv
	INTEGER(I4B), SAVE :: nvold=-1
	REAL(SP), DIMENSION(size(yz)) :: yy,v,c,b,b1,ddy
	REAL(SP), DIMENSION(:,:), ALLOCATABLE, SAVE :: d
	REAL(SP), DIMENSION(IEST_MAX), SAVE :: fx,x
	nv=assert_eq(size(yz),size(dy),size(yest),'rzextr')
	if (iest > IEST_MAX) call &
		nrerror('rzextr: probable misuse, too much extrapolation')
	if (nv /= nvold) then
		if (allocated(d)) deallocate(d)
		allocate(d(nv,IEST_MAX))
		nvold=nv
	end if
	x(iest)=xest
	if (iest == 1) then
		yz=yest
		d(:,1)=yest
		dy=yest
	else
		fx(2:iest)=x(iest-1:1:-1)/xest
		yy=yest
		v=d(1:nv,1)
		c=yy
		d(1:nv,1)=yy
		do k=2,iest
			b1=fx(k)*v
			b=b1-c
			where (b /= 0.0)
				b=(c-v)/b
				ddy=c*b
				c=b1*b
			elsewhere
				ddy=v
			end where
			if (k /= iest) v=d(1:nv,k)
			d(1:nv,k)=ddy
			yy=yy+ddy
		end do
		dy=ddy
		yz=yy
	end if
	END SUBROUTINE rzextr

! -----------------------------------------------------------

	FUNCTION savgol(nl,nrr,ld,m)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert,poly
	use mo_nr, ONLY : lubksb,ludcmp
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: nl,nrr,ld,m
	REAL(SP), DIMENSION(nl+nrr+1) :: savgol
	INTEGER(I4B) :: imj,ipj,mm,np
	INTEGER(I4B), DIMENSION(m+1) :: indx
	REAL(SP) :: d,sm
	REAL(SP), DIMENSION(m+1) :: b
	REAL(SP), DIMENSION(m+1,m+1) :: a
	INTEGER(I4B) :: irng(nl+nrr+1)
	call assert(nl >= 0, nrr >= 0, ld <= m, nl+nrr >= m, 'savgol args')
	do ipj=0,2*m
		sm=sum(arth(1.0_sp,1.0_sp,nrr)**ipj)+&
			sum(arth(-1.0_sp,-1.0_sp,nl)**ipj)
		if (ipj == 0) sm=sm+1.0_sp
		mm=min(ipj,2*m-ipj)
		do imj=-mm,mm,2
			a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sm
		end do
	end do
	call ludcmp(a(:,:),indx(:),d)
	b(:)=0.0
	b(ld+1)=1.0
	call lubksb(a(:,:),indx(:),b(:))
	savgol(:)=0.0
	irng(:)=arth(-nl,1,nrr+nl+1)
	np=nl+nrr+1
	savgol(mod(np-irng(:),np)+1)=poly(real(irng(:),sp),b(:))
	END FUNCTION savgol

! -----------------------------------------------------------

	SUBROUTINE scrsho(func)
	use mo_nrtype
	IMPLICIT NONE
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ISCR=60,JSCR=21
	INTEGER(I4B) :: i,j,jz
	REAL(SP) :: dx,dyj,x,x1,x2,ybig,ysml
	REAL(SP), DIMENSION(ISCR) :: y
	CHARACTER(1), DIMENSION(ISCR,JSCR) :: scr
	CHARACTER(1) :: blank=' ',zero='-',yy='l',xx='-',ff='x'
	do
		write (*,*) ' Enter x1,x2 (= to stop)'
		read (*,*) x1,x2
		if (x1 == x2) RETURN
		scr(1,1:JSCR)=yy
		scr(ISCR,1:JSCR)=yy
		scr(2:ISCR-1,1)=xx
		scr(2:ISCR-1,JSCR)=xx
		scr(2:ISCR-1,2:JSCR-1)=blank
		dx=(x2-x1)/(ISCR-1)
		x=x1
		do i=1,ISCR
			y(i)=func(x)
			x=x+dx
		end do
		ysml=min(minval(y(:)),0.0_sp)
		ybig=max(maxval(y(:)),0.0_sp)
		if (ybig == ysml) ybig=ysml+1.0
		dyj=(JSCR-1)/(ybig-ysml)
		jz=1-ysml*dyj
		scr(1:ISCR,jz)=zero
		do i=1,ISCR
			j=1+(y(i)-ysml)*dyj
			scr(i,j)=ff
		end do
		write (*,'(1x,1p,e10.3,1x,80a1)') ybig,(scr(i,JSCR),i=1,ISCR)
		do j=JSCR-1,2,-1
			write (*,'(12x,80a1)') (scr(i,j),i=1,ISCR)
		end do
		write (*,'(1x,1p,e10.3,1x,80a1)') ysml,(scr(i,1),i=1,ISCR)
		write (*,'(12x,1p,e10.3,40x,e10.3)') x1,x2
	end do
	END SUBROUTINE scrsho

! -----------------------------------------------------------

	FUNCTION select(k,arr)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,swap
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: k
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	REAL(SP) :: select
	INTEGER(I4B) :: i,r,j,l,n
	REAL(SP) :: a
	n=size(arr)
	call assert(k >= 1, k <= n, 'select args')
	l=1
	r=n
	do
		if (r-l <= 1) then
			if (r-l == 1) call swap(arr(l),arr(r),arr(l)>arr(r))
			select=arr(k)
			RETURN
		else
			i=(l+r)/2
			call swap(arr(i),arr(l+1))
			call swap(arr(l),arr(r),arr(l)>arr(r))
			call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
			call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
			i=l+1
			j=r
			a=arr(l+1)
			do
				do
					i=i+1
					if (arr(i) >= a) exit
				end do
				do
					j=j-1
					if (arr(j) <= a) exit
				end do
				if (j < i) exit
				call swap(arr(i),arr(j))
			end do
			arr(l+1)=arr(j)
			arr(j)=a
			if (j >= k) r=j-1
			if (j <= k) l=i
		end if
	end do
	END FUNCTION select

! -----------------------------------------------------------

	FUNCTION select_bypack(k,arr)
	use mo_nrtype
          use mo_nrutil, ONLY : array_copy,assert,swap
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: k
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	REAL(SP) :: select_bypack
	LOGICAL, DIMENSION(size(arr)) :: mask
	REAL(SP), DIMENSION(size(arr)) :: temp
	INTEGER(I4B) :: i,r,j,l,n,nl,nerr
	REAL(SP) :: a
	n=size(arr)
	call assert(k >= 1, k <= n, 'select_bypack args')
	l=1
	r=n
	do
		if (r-l <= 1) exit
		i=(l+r)/2
		call swap(arr(l),arr(i),arr(l)>arr(i))
		call swap(arr(i),arr(r),arr(i)>arr(r))
		call swap(arr(l),arr(i),arr(l)>arr(i))
		a=arr(i)
		mask(l:r) = (arr(l:r) <= a)
		mask(i) = .false.
		call array_copy(pack(arr(l:r),mask(l:r)),temp(l:),nl,nerr)
		j=l+nl
		mask(i) = .true.
		temp(j+1:r)=pack(arr(l:r),.not. mask(l:r))
		temp(j)=a
		arr(l:r)=temp(l:r)
		if (k > j) then
			l=j+1
		else if (k < j) then
			r=j-1
		else
			l=j
			r=j
		end if
	end do
	if (r-l == 1) call swap(arr(l),arr(r),arr(l)>arr(r))
	select_bypack=arr(k)
	END FUNCTION select_bypack

! -----------------------------------------------------------

	SUBROUTINE select_heap(arr,heap)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror,swap
	use mo_nr, ONLY : sort
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	REAL(SP), DIMENSION(:), INTENT(OUT) :: heap
	INTEGER(I4B) :: i,j,k,m,n
	m=size(heap)
	n=size(arr)
	if (m > n/2 .or. m < 1) call nrerror('probable misuse of select_heap')
	heap=arr(1:m)
	call sort(heap)
	do i=m+1,n
		if (arr(i) > heap(1)) then
			heap(1)=arr(i)
			j=1
			do
				k=2*j
				if (k > m) exit
				if (k /= m) then
					if (heap(k) > heap(k+1)) k=k+1
				end if
				if (heap(j) <= heap(k)) exit
				call swap(heap(k),heap(j))
				j=k
			end do
		end if
	end do
	END SUBROUTINE select_heap

! -----------------------------------------------------------

	FUNCTION select_inplace(k,arr)
	use mo_nrtype
	use mo_nr, ONLY : select
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: k
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	REAL(SP) :: select_inplace
	REAL(SP), DIMENSION(size(arr)) :: tarr
	tarr=arr
	select_inplace=select(k,tarr)
	END FUNCTION select_inplace

! -----------------------------------------------------------

MODULE sfroid_data
	use mo_nrtype
	INTEGER(I4B), PARAMETER :: M=41
	INTEGER(I4B) :: mm,n
	REAL(SP) :: anorm,c2,h
	REAL(SP), DIMENSION(M) :: x
END MODULE sfroid_data

	PROGRAM sfroid
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	use mo_nr, ONLY : plgndr,solvde
	USE sfroid_data
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NE=3,NB=1
	INTEGER(I4B) :: itmax
	INTEGER(I4B), DIMENSION(NE) :: indexv
	REAL(SP) :: conv,slowc
	REAL(SP), DIMENSION(M) :: deriv,fac1,fac2
	REAL(SP), DIMENSION(NE) :: scalv
	REAL(SP), DIMENSION(NE,M) :: y
	itmax=100
	conv=5.0e-6_sp
	slowc=1.0
	h=1.0_sp/(M-1)
	c2=0.0
	write(*,*) 'ENTER M,N'
	read(*,*) mm,n
	indexv(1:3)=merge( (/ 1, 2, 3 /), (/ 2, 1, 3 /), (mod(n+mm,2) == 1) )
	anorm=1.0
	if (mm /= 0) then
		anorm=(-0.5_sp)**mm*product(&
			arth(n+1,1,mm)*arth(real(n,sp),-1.0_sp,mm)/arth(1,1,mm))
	end if
	x(1:M-1)=arth(0,1,M-1)*h
	fac1(1:M-1)=1.0_sp-x(1:M-1)**2
	fac2(1:M-1)=fac1(1:M-1)**(-mm/2.0_sp)
	y(1,1:M-1)=plgndr(n,mm,x(1:M-1))*fac2(1:M-1)
	deriv(1:M-1)=-((n-mm+1)*plgndr(n+1,mm,x(1:M-1))-(n+1)*&
		x(1:M-1)*plgndr(n,mm,x(1:M-1)))/fac1(1:M-1)
	y(2,1:M-1)=mm*x(1:M-1)*y(1,1:M-1)/fac1(1:M-1)+deriv(1:M-1)*fac2(1:M-1)
	y(3,1:M-1)=n*(n+1)-mm*(mm+1)
	x(M)=1.0
	y(1,M)=anorm
	y(3,M)=n*(n+1)-mm*(mm+1)
	y(2,M)=(y(3,M)-c2)*y(1,M)/(2.0_sp*(mm+1.0_sp))
	scalv(1:3)=(/ abs(anorm), max(abs(anorm),y(2,M)), max(1.0_sp,y(3,M)) /)
	do
		write (*,*) 'ENTER C**2 OR 999 TO END'
		read (*,*) c2
		if (c2 == 999.0) exit
		call solvde(itmax,conv,slowc,scalv,indexv,NB,y)
		write (*,*) ' M = ',mm,'  N = ',n,&
			'  C**2 = ',c2,'  LAMBDA = ',y(3,1)+mm*(mm+1)
	end do
	END PROGRAM sfroid

! -----------------------------------------------------------

!	FUNCTION shoot(v) is named "funcv" for use with "newt"
	FUNCTION funcv(v)
	use mo_nrtype
	use mo_nr, ONLY : odeint,rkqs
	USE sphoot_caller, ONLY: nvar,x1,x2
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: v
	REAL(SP), DIMENSION(size(v)) :: funcv
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	REAL(SP) :: h1,hmin
	REAL(SP), DIMENSION(nvar) :: y
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
!BL
		SUBROUTINE load(x1,v,y)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x1
		REAL(SP), DIMENSION(:), INTENT(IN) :: v
		REAL(SP), DIMENSION(:), INTENT(OUT) :: y
		END SUBROUTINE load
!BL
		SUBROUTINE score(x2,y,f)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x2
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: f
		END SUBROUTINE score
	END INTERFACE
	h1=(x2-x1)/100.0_sp
	hmin=0.0
	call load(x1,v,y)
	call odeint(y,x1,x2,EPS,h1,hmin,derivs,rkqs)
	call score(x2,y,funcv)
	END FUNCTION funcv

! -----------------------------------------------------------

!	FUNCTION shootf(v) is named "funcv" for use with "newt"
	FUNCTION funcv(v)
	use mo_nrtype
	use mo_nr, ONLY : odeint,rkqs
	USE sphfpt_caller, ONLY : x1,x2,xf,nn2
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: v
	REAL(SP), DIMENSION(size(v)) :: funcv
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	REAL(SP) :: h1,hmin
	REAL(SP), DIMENSION(size(v)) :: f1,f2,y
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
!BL
		SUBROUTINE load1(x1,v1,y)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x1
		REAL(SP), DIMENSION(:), INTENT(IN) :: v1
		REAL(SP), DIMENSION(:), INTENT(OUT) :: y
		END SUBROUTINE load1
!BL
		SUBROUTINE load2(x2,v2,y)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x2
		REAL(SP), DIMENSION(:), INTENT(IN) :: v2
		REAL(SP), DIMENSION(:), INTENT(OUT) :: y
		END SUBROUTINE load2
!BL
		SUBROUTINE score(x2,y,f)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x2
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: f
		END SUBROUTINE score
	END INTERFACE
	h1=(x2-x1)/100.0_sp
	hmin=0.0
	call load1(x1,v,y)
	call odeint(y,x1,xf,EPS,h1,hmin,derivs,rkqs)
	call score(xf,y,f1)
	call load2(x2,v(nn2+1:),y)
	call odeint(y,x2,xf,EPS,h1,hmin,derivs,rkqs)
	call score(xf,y,f2)
	funcv(:)=f1(:)-f2(:)
	END FUNCTION funcv

! -----------------------------------------------------------

	SUBROUTINE simplx(a,m1,m2,m3,icase,izrov,iposv)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq,ifirstloc,imaxloc,&
		nrerror,outerprod,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B), INTENT(IN) :: m1,m2,m3
	INTEGER(I4B), INTENT(OUT) :: icase
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: izrov,iposv
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	INTEGER(I4B) :: ip,k,kh,kp,m12,nl1,nl2,m,n
	INTEGER(I4B), DIMENSION(size(a,2)) :: l1
	INTEGER(I4B), DIMENSION(size(a,1)) :: l2,l3
	REAL(SP) :: bmax
	LOGICAL(LGT) :: init,proceed
	m=assert_eq(size(a,1)-2,size(iposv),'simplx: m')
	n=assert_eq(size(a,2)-1,size(izrov),'simplx: n')
	if (m /= m1+m2+m3) call nrerror('simplx: bad input constraint counts')
	if (any(a(2:m+1,1) < 0.0)) call nrerror('bad input tableau in simplx')
	nl1=n
	l1(1:n)=arth(1,1,n)
	izrov(:)=l1(1:n)
	nl2=m
	l2(1:m)=arth(1,1,m)
	iposv(:)=n+l2(1:m)
	l3(1:m2)=1
	init=.true.
	main_procedure: do
	phase1: do
		if (init) then
			init=.false.
			if (m2+m3 == 0) then
				proceed=.true.
				exit phase1
			else
				proceed=.false.
			end if
			a(m+2,1:n+1)=-sum(a(m1+2:m+1,1:n+1),dim=1)
		end if
		kp=l1(imaxloc(a(m+2,l1(1:nl1)+1)))
		bmax=a(m+2,kp+1)
		phase1a: do
			if (bmax <= EPS .and. a(m+2,1) < -EPS) then
				icase=-1
				RETURN
			else if (bmax <= EPS .and. a(m+2,1) <= EPS) then
				m12=m1+m2+1
				do ip=m12,m
					if (iposv(ip) == ip+n) then
						kp=l1(imaxloc(abs(a(ip+1,l1(1:nl1)+1))))
						bmax=a(ip+1,kp+1)
						if (bmax > 0.0) exit phase1a
					end if
				end do
				proceed=.true.
				m12=m12-1
				where (spread(l3(1:m12-m1),2,n+1) == 1) &
					a(m1+2:m12+1,1:n+1)=-a(m1+2:m12+1,1:n+1)
				exit phase1
			end if
			call simp1
			if (ip == 0) then
				icase=-1
				RETURN
			end if
			exit phase1a
		end do phase1a
		phase1b: do
			call simp2(m+1,n)
			if (iposv(ip) >= n+m1+m2+1) then
				k=ifirstloc(l1(1:nl1) == kp)
				nl1=nl1-1
				l1(k:nl1)=l1(k+1:nl1+1)
			else
				if (iposv(ip) < n+m1+1) exit phase1b
				kh=iposv(ip)-m1-n
				if (l3(kh) == 0) exit phase1b
				l3(kh)=0
			end if
			a(m+2,kp+1)=a(m+2,kp+1)+1.0_sp
			a(1:m+2,kp+1)=-a(1:m+2,kp+1)
			exit phase1b
		end do phase1b
		call swap(izrov(kp),iposv(ip))
		if (proceed) exit phase1
	end do phase1
	phase2: do
		kp=l1(imaxloc(a(1,l1(1:nl1)+1)))
		bmax=a(1,kp+1)
		if (bmax <= 0.0) then
			icase=0
			RETURN
		end if
		call simp1
		if (ip == 0) then
			icase=1
			RETURN
		end if
		call simp2(m,n)
		call swap(izrov(kp),iposv(ip))
		if (.not. proceed) cycle main_procedure
	end do phase2
	exit main_procedure
	end do main_procedure
	CONTAINS
!BL
	SUBROUTINE simp1
	IMPLICIT NONE
	INTEGER(I4B) :: i,ii,k
	REAL(SP) :: q,q0,q1,qp
	ip=0
	i=ifirstloc(a(l2(1:nl2)+1,kp+1) < -EPS)
	if (i > nl2) RETURN
	q1=-a(l2(i)+1,1)/a(l2(i)+1,kp+1)
	ip=l2(i)
	do i=i+1,nl2
		ii=l2(i)
		if (a(ii+1,kp+1) < -EPS) then
			q=-a(ii+1,1)/a(ii+1,kp+1)
			if (q < q1) then
				ip=ii
				q1=q
			else if (q == q1) then
				do k=1,n
					qp=-a(ip+1,k+1)/a(ip+1,kp+1)
					q0=-a(ii+1,k+1)/a(ii+1,kp+1)
					if (q0 /= qp) exit
				end do
				if (q0 < qp) ip=ii
			end if
		end if
	end do
	END SUBROUTINE simp1
!BL
	SUBROUTINE simp2(i1,k1)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: i1,k1
	INTEGER(I4B) :: ip1,kp1
	REAL(SP) :: piv
	INTEGER(I4B), DIMENSION(k1) :: icol
	INTEGER(I4B), DIMENSION(i1) :: irow
	INTEGER(I4B), DIMENSION(max(i1,k1)+1) :: itmp
	ip1=ip+1
	kp1=kp+1
	piv=1.0_sp/a(ip1,kp1)
	itmp(1:k1+1)=arth(1,1,k1+1)
	icol=pack(itmp(1:k1+1),itmp(1:k1+1) /= kp1)
	itmp(1:i1+1)=arth(1,1,i1+1)
	irow=pack(itmp(1:i1+1),itmp(1:i1+1) /= ip1)
	a(irow,kp1)=a(irow,kp1)*piv
	a(irow,icol)=a(irow,icol)-outerprod(a(irow,kp1),a(ip1,icol))
	a(ip1,icol)=-a(ip1,icol)*piv
	a(ip1,kp1)=piv
	END SUBROUTINE simp2
	END SUBROUTINE simplx

! -----------------------------------------------------------

	SUBROUTINE simpr(y,dydx,dfdx,dfdy,xs,htot,nstep,yout,derivs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,diagadd
	use mo_nr, ONLY : lubksb,ludcmp
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: xs,htot
	REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx,dfdx
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: dfdy
	INTEGER(I4B), INTENT(IN) :: nstep
	REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: ndum,nn
	INTEGER(I4B), DIMENSION(size(y)) :: indx
	REAL(SP) :: d,h,x
	REAL(SP), DIMENSION(size(y)) :: del,ytemp
	REAL(SP), DIMENSION(size(y),size(y)) :: a
	ndum=assert_eq((/size(y),size(dydx),size(dfdx),size(dfdy,1),&
		size(dfdy,2),size(yout)/),'simpr')
	h=htot/nstep
	a(:,:)=-h*dfdy(:,:)
	call diagadd(a,1.0_sp)
	call ludcmp(a,indx,d)
	yout=h*(dydx+h*dfdx)
	call lubksb(a,indx,yout)
	del=yout
	ytemp=y+del
	x=xs+h
	call derivs(x,ytemp,yout)
	do nn=2,nstep
		yout=h*yout-del
		call lubksb(a,indx,yout)
		del=del+2.0_sp*yout
		ytemp=ytemp+del
		x=x+h
		call derivs(x,ytemp,yout)
	end do
	yout=h*yout-del
	call lubksb(a,indx,yout)
	yout=ytemp+yout
	END SUBROUTINE simpr

! -----------------------------------------------------------

	SUBROUTINE sinft(y)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,cumsum,zroots_unity
	use mo_nr, ONLY : realft
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	REAL(SP), DIMENSION(size(y)/2+1) :: wi
	REAL(SP), DIMENSION(size(y)/2) :: y1,y2
	INTEGER(I4B) :: n,nh
	n=size(y)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in sinft')
	nh=n/2
	wi=aimag(zroots_unity(n+n,nh+1))
	y(1)=0.0
	y1=wi(2:nh+1)*(y(2:nh+1)+y(n:nh+1:-1))
	y2=0.5_sp*(y(2:nh+1)-y(n:nh+1:-1))
	y(2:nh+1)=y1+y2
	y(n:nh+1:-1)=y1-y2
	call realft(y,+1)
	y(1)=0.5_sp*y(1)
	y(2)=0.0
	y1=cumsum(y(1:n-1:2))
	y(1:n-1:2)=y(2:n:2)
	y(2:n:2)=y1
	END SUBROUTINE sinft

! -----------------------------------------------------------

	SUBROUTINE slvsm2(u,rhs)
	use mo_nrtype
	IMPLICIT NONE
	REAL(DP), DIMENSION(3,3), INTENT(OUT) :: u
	REAL(DP), DIMENSION(3,3), INTENT(IN) :: rhs
	REAL(DP) :: disc,fact,h
	u=0.0
	h=0.5_dp
	fact=2.0_dp/h**2
	disc=sqrt(fact**2+rhs(2,2))
	u(2,2)=-rhs(2,2)/(fact+disc)
	END SUBROUTINE slvsm2

! -----------------------------------------------------------

	SUBROUTINE slvsml(u,rhs)
	use mo_nrtype
	IMPLICIT NONE
	REAL(DP), DIMENSION(3,3), INTENT(OUT) :: u
	REAL(DP), DIMENSION(3,3), INTENT(IN) :: rhs
	REAL(DP) :: h
	u=0.0
	h=0.5_dp
	u(2,2)=-h*h*rhs(2,2)/4.0_dp
	END SUBROUTINE slvsml

! -----------------------------------------------------------

	SUBROUTINE sncndn(uu,emmc,sn,cn,dn)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: uu,emmc
	REAL(SP), INTENT(OUT) :: sn,cn,dn
	REAL(SP), PARAMETER :: CA=0.0003_sp
	INTEGER(I4B), PARAMETER :: MAXIT=13
	INTEGER(I4B) :: i,ii,l
	REAL(SP) :: a,b,c,d,emc,u
	REAL(SP), DIMENSION(MAXIT) :: em,en
	LOGICAL(LGT) :: bo
	emc=emmc
	u=uu
	if (emc /= 0.0) then
		bo=(emc < 0.0)
		if (bo) then
			d=1.0_sp-emc
			emc=-emc/d
			d=sqrt(d)
			u=d*u
		end if
		a=1.0
		dn=1.0
		do i=1,MAXIT
			l=i
			em(i)=a
			emc=sqrt(emc)
			en(i)=emc
			c=0.5_sp*(a+emc)
			if (abs(a-emc) <= CA*a) exit
			emc=a*emc
			a=c
		end do
		if (i > MAXIT) call nrerror('sncndn: convergence failed')
		u=c*u
		sn=sin(u)
		cn=cos(u)
		if (sn /= 0.0) then
			a=cn/sn
			c=a*c
			do ii=l,1,-1
				b=em(ii)
				a=c*a
				c=dn*c
				dn=(en(ii)+a)/(b+a)
				a=c/b
			end do
			a=1.0_sp/sqrt(c**2+1.0_sp)
			sn=sign(a,sn)
			cn=c*sn
		end if
		if (bo) then
			a=dn
			dn=cn
			cn=a
			sn=sn/d
		end if
	else
		cn=1.0_sp/cosh(u)
		dn=cn
		sn=tanh(u)
	end if
	END SUBROUTINE sncndn

! -----------------------------------------------------------

	FUNCTION snrm(sx,itol)
	use mo_nrtype
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: sx
	INTEGER(I4B), INTENT(IN) :: itol
	REAL(DP) :: snrm
	if (itol <= 3) then
		snrm=sqrt(dot_product(sx,sx))
	else
		snrm=maxval(abs(sx))
	end if
	END FUNCTION snrm

! -----------------------------------------------------------

	SUBROUTINE sobseq(x,init)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x
	INTEGER(I4B), OPTIONAL, INTENT(IN) :: init
	INTEGER(I4B), PARAMETER :: MAXBIT=30,MAXDIM=6
	REAL(SP), SAVE :: fac
	INTEGER(I4B) :: i,im,ipp,j,k,l
	INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE:: iu
	INTEGER(I4B), SAVE :: in
	INTEGER(I4B), DIMENSION(MAXDIM), SAVE :: ip,ix,mdeg
	INTEGER(I4B), DIMENSION(MAXDIM*MAXBIT), SAVE :: iv
	DATA ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
	DATA iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
	if (present(init)) then
		ix=0
		in=0
		if (iv(1) /= 1) RETURN
		fac=1.0_sp/2.0_sp**MAXBIT
		allocate(iu(MAXDIM,MAXBIT))
		iu=reshape(iv,shape(iu))
		do k=1,MAXDIM
			do j=1,mdeg(k)
				iu(k,j)=iu(k,j)*2**(MAXBIT-j)
			end do
			do j=mdeg(k)+1,MAXBIT
				ipp=ip(k)
				i=iu(k,j-mdeg(k))
				i=ieor(i,i/2**mdeg(k))
				do l=mdeg(k)-1,1,-1
					if (btest(ipp,0)) i=ieor(i,iu(k,j-l))
					ipp=ipp/2
				end do
				iu(k,j)=i
			end do
		end do
		iv=reshape(iu,shape(iv))
		deallocate(iu)
	else
		im=in
		do j=1,MAXBIT
			if (.not. btest(im,0)) exit
			im=im/2
		end do
		if (j > MAXBIT) call nrerror('MAXBIT too small in sobseq')
		im=(j-1)*MAXDIM
		j=min(size(x),MAXDIM)
		ix(1:j)=ieor(ix(1:j),iv(1+im:j+im))
		x(1:j)=ix(1:j)*fac
		in=in+1
	end if
	END SUBROUTINE sobseq

! -----------------------------------------------------------

	SUBROUTINE solvde(itmax,conv,slowc,scalv,indexv,nb,y)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,imaxloc,nrerror
	use mo_nr, ONLY : difeq
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: itmax,nb
	REAL(SP), INTENT(IN) :: conv,slowc
	REAL(SP), DIMENSION(:), INTENT(IN) :: scalv
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indexv
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: y
	INTEGER(I4B) :: ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8,&
		j9,jc1,jcf,jv,k,k1,k2,km,kp,m,ne,nvars
	INTEGER(I4B), DIMENSION(size(scalv)) :: kmax
	REAL(SP) :: err,fac
	REAL(SP), DIMENSION(size(scalv)) :: ermax
	REAL(SP), DIMENSION(size(scalv),2*size(scalv)+1) :: s
	REAL(SP), DIMENSION(size(scalv),size(scalv)-nb+1,size(y,2)+1) :: c
	ne=assert_eq(size(scalv),size(indexv),size(y,1),'solvde: ne')
	m=size(y,2)
	k1=1
	k2=m
	nvars=ne*m
	j1=1
	j2=nb
	j3=nb+1
	j4=ne
	j5=j4+j1
	j6=j4+j2
	j7=j4+j3
	j8=j4+j4
	j9=j8+j1
	ic1=1
	ic2=ne-nb
	ic3=ic2+1
	ic4=ne
	jc1=1
	jcf=ic3
	do it=1,itmax
		k=k1
		call difeq(k,k1,k2,j9,ic3,ic4,indexv,s,y)
		call pinvs(ic3,ic4,j5,j9,jc1,k1,c,s)
		do k=k1+1,k2
			kp=k-1
			call difeq(k,k1,k2,j9,ic1,ic4,indexv,s,y)
			call red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,s)
			call pinvs(ic1,ic4,j3,j9,jc1,k,c,s)
		end do
		k=k2+1
		call difeq(k,k1,k2,j9,ic1,ic2,indexv,s,y)
		call red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,c,s)
		call pinvs(ic1,ic2,j7,j9,jcf,k2+1,c,s)
		call bksub(ne,nb,jcf,k1,k2,c)
		do j=1,ne
			jv=indexv(j)
			km=imaxloc(abs(c(jv,1,k1:k2)))+k1-1
			ermax(j)=c(jv,1,km)
			kmax(j)=km
		end do
		ermax(:)=ermax(:)/scalv(:)
		err=sum(sum(abs(c(indexv(:),1,k1:k2)),dim=2)/scalv(:))/nvars
		fac=slowc/max(slowc,err)
		y(:,k1:k2)=y(:,k1:k2)-fac*c(indexv(:),1,k1:k2)
		write(*,'(1x,i4,2f12.6)') it,err,fac
		if (err < conv) RETURN
	end do
	call nrerror('itmax exceeded in solvde')
	CONTAINS
!BL
	SUBROUTINE bksub(ne,nb,jf,k1,k2,c)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ne,nb,jf,k1,k2
	REAL(SP), DIMENSION(:,:,:), INTENT(INOUT) :: c
	INTEGER(I4B) :: im,k,nbf
	nbf=ne-nb
	im=1
	do k=k2,k1,-1
		if (k == k1) im=nbf+1
		c(im:ne,jf,k)=c(im:ne,jf,k)-matmul(c(im:ne,1:nbf,k),c(1:nbf,jf,k+1))
	end do
	c(1:nb,1,k1:k2)=c(1+nbf:nb+nbf,jf,k1:k2)
	c(1+nb:nbf+nb,1,k1:k2)=c(1:nbf,jf,k1+1:k2+1)
	END SUBROUTINE bksub
!BL
	SUBROUTINE pinvs(ie1,ie2,je1,jsf,jc1,k,c,s)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ie1,ie2,je1,jsf,jc1,k
	REAL(SP), DIMENSION(:,:,:), INTENT(OUT) :: c
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: s
	INTEGER(I4B) :: i,icoff,id,ipiv,jcoff,je2,jp,jpiv,js1
	INTEGER(I4B), DIMENSION(ie2) :: indxr
	REAL(SP) :: big,piv,pivinv
	REAL(SP), DIMENSION(ie2) :: pscl
	je2=je1+ie2-ie1
	js1=je2+1
	pscl(ie1:ie2)=maxval(abs(s(ie1:ie2,je1:je2)),dim=2)
	if (any(pscl(ie1:ie2) == 0.0)) &
		call nrerror('singular matrix, row all 0 in pinvs')
	pscl(ie1:ie2)=1.0_sp/pscl(ie1:ie2)
	indxr(ie1:ie2)=0
	do id=ie1,ie2
		piv=0.0
		do i=ie1,ie2
			if (indxr(i) == 0) then
				jp=imaxloc(abs(s(i,je1:je2)))+je1-1
				big=abs(s(i,jp))
				if (big*pscl(i) > piv) then
					ipiv=i
					jpiv=jp
					piv=big*pscl(i)
				end if
			end if
		end do
		if (s(ipiv,jpiv) == 0.0) call nrerror('singular matrix in pinvs')
		indxr(ipiv)=jpiv
		pivinv=1.0_sp/s(ipiv,jpiv)
		s(ipiv,je1:jsf)=s(ipiv,je1:jsf)*pivinv
		s(ipiv,jpiv)=1.0
		do i=ie1,ie2
			if (indxr(i) /= jpiv .and. s(i,jpiv) /= 0.0) then
				s(i,je1:jsf)=s(i,je1:jsf)-s(i,jpiv)*s(ipiv,je1:jsf)
				s(i,jpiv)=0.0
			end if
		end do
	end do
	jcoff=jc1-js1
	icoff=ie1-je1
	c(indxr(ie1:ie2)+icoff,js1+jcoff:jsf+jcoff,k)=s(ie1:ie2,js1:jsf)
	END SUBROUTINE pinvs
!BL
	SUBROUTINE red(iz1,iz2,jz1,jz2,jm1,jm2,jmf,ic1,jc1,jcf,kc,c,s)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: iz1,iz2,jz1,jz2,jm1,jm2,jmf,ic1,jc1,jcf,kc
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: s
	REAL(SP), DIMENSION(:,:,:), INTENT(IN) :: c
	INTEGER(I4B) :: ic,l,loff
	loff=jc1-jm1
	ic=ic1
	do j=jz1,jz2
		do l=jm1,jm2
			s(iz1:iz2,l)=s(iz1:iz2,l)-s(iz1:iz2,j)*c(ic,l+loff,kc)
		end do
		s(iz1:iz2,jmf)=s(iz1:iz2,jmf)-s(iz1:iz2,j)*c(ic,jcf,kc)
		ic=ic+1
	end do
	END SUBROUTINE red
	END SUBROUTINE solvde

! -----------------------------------------------------------

	SUBROUTINE sor(a,b,c,d,e,f,u,rjac)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: a,b,c,d,e,f
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: u
	REAL(DP), INTENT(IN) :: rjac
	INTEGER(I4B), PARAMETER :: MAXITS=1000
	REAL(DP), PARAMETER :: EPS=1.0e-5_dp
	REAL(DP), DIMENSION(size(a,1),size(a,1)) :: resid
	INTEGER(I4B) :: jmax,jm1,jm2,jm3,n
	REAL(DP) :: anorm,anormf,omega
	jmax=assert_eq((/size(a,1),size(a,2),size(b,1),size(b,2), &
		size(c,1),size(c,2),size(d,1),size(d,2),size(e,1), &
		size(e,2),size(f,1),size(f,2),size(u,1),size(u,2)/),'sor')
	jm1=jmax-1
	jm2=jmax-2
	jm3=jmax-3
	anormf=sum(abs(f(2:jm1,2:jm1)))
	omega=1.0
	do n=1,MAXITS
		resid(2:jm1:2,2:jm1:2)=a(2:jm1:2,2:jm1:2)*u(3:jmax:2,2:jm1:2)+&
			b(2:jm1:2,2:jm1:2)*u(1:jm2:2,2:jm1:2)+&
			c(2:jm1:2,2:jm1:2)*u(2:jm1:2,3:jmax:2)+&
			d(2:jm1:2,2:jm1:2)*u(2:jm1:2,1:jm2:2)+&
			e(2:jm1:2,2:jm1:2)*u(2:jm1:2,2:jm1:2)-f(2:jm1:2,2:jm1:2)
		u(2:jm1:2,2:jm1:2)=u(2:jm1:2,2:jm1:2)-omega*&
			resid(2:jm1:2,2:jm1:2)/e(2:jm1:2,2:jm1:2)
		resid(3:jm2:2,3:jm2:2)=a(3:jm2:2,3:jm2:2)*u(4:jm1:2,3:jm2:2)+&
			b(3:jm2:2,3:jm2:2)*u(2:jm3:2,3:jm2:2)+&
			c(3:jm2:2,3:jm2:2)*u(3:jm2:2,4:jm1:2)+&
			d(3:jm2:2,3:jm2:2)*u(3:jm2:2,2:jm3:2)+&
			e(3:jm2:2,3:jm2:2)*u(3:jm2:2,3:jm2:2)-f(3:jm2:2,3:jm2:2)
		u(3:jm2:2,3:jm2:2)=u(3:jm2:2,3:jm2:2)-omega*&
			resid(3:jm2:2,3:jm2:2)/e(3:jm2:2,3:jm2:2)
		omega=merge(1.0_dp/(1.0_dp-0.5_dp*rjac**2), &
			1.0_dp/(1.0_dp-0.25_dp*rjac**2*omega), n == 1)
		resid(3:jm2:2,2:jm1:2)=a(3:jm2:2,2:jm1:2)*u(4:jm1:2,2:jm1:2)+&
			b(3:jm2:2,2:jm1:2)*u(2:jm3:2,2:jm1:2)+&
			c(3:jm2:2,2:jm1:2)*u(3:jm2:2,3:jmax:2)+&
			d(3:jm2:2,2:jm1:2)*u(3:jm2:2,1:jm2:2)+&
			e(3:jm2:2,2:jm1:2)*u(3:jm2:2,2:jm1:2)-f(3:jm2:2,2:jm1:2)
		u(3:jm2:2,2:jm1:2)=u(3:jm2:2,2:jm1:2)-omega*&
			resid(3:jm2:2,2:jm1:2)/e(3:jm2:2,2:jm1:2)
		resid(2:jm1:2,3:jm2:2)=a(2:jm1:2,3:jm2:2)*u(3:jmax:2,3:jm2:2)+&
			b(2:jm1:2,3:jm2:2)*u(1:jm2:2,3:jm2:2)+&
			c(2:jm1:2,3:jm2:2)*u(2:jm1:2,4:jm1:2)+&
			d(2:jm1:2,3:jm2:2)*u(2:jm1:2,2:jm3:2)+&
			e(2:jm1:2,3:jm2:2)*u(2:jm1:2,3:jm2:2)-f(2:jm1:2,3:jm2:2)
		u(2:jm1:2,3:jm2:2)=u(2:jm1:2,3:jm2:2)-omega*&
			resid(2:jm1:2,3:jm2:2)/e(2:jm1:2,3:jm2:2)
		omega=1.0_dp/(1.0_dp-0.25_dp*rjac**2*omega)
		anorm=sum(abs(resid(2:jm1,2:jm1)))
		if (anorm < EPS*anormf) exit
	end do
	if (n > MAXITS) call nrerror('MAXITS exceeded in sor')
	END SUBROUTINE sor

! -----------------------------------------------------------

	SUBROUTINE sort(arr)
	use mo_nrtype
          use mo_nrutil, ONLY : swap,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
	REAL(SP) :: a
	INTEGER(I4B) :: n,k,i,j,jstack,l,r
	INTEGER(I4B), DIMENSION(NSTACK) :: istack
	n=size(arr)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				a=arr(j)
				do i=j-1,l,-1
					if (arr(i) <= a) exit
					arr(i+1)=arr(i)
				end do
				arr(i+1)=a
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(arr(k),arr(l+1))
			call swap(arr(l),arr(r),arr(l)>arr(r))
			call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
			call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
			i=l+1
			j=r
			a=arr(l+1)
			do
				do
					i=i+1
					if (arr(i) >= a) exit
				end do
				do
					j=j-1
					if (arr(j) <= a) exit
				end do
				if (j < i) exit
				call swap(arr(i),arr(j))
			end do
			arr(l+1)=arr(j)
			arr(j)=a
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('sort: NSTACK too small')
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	END SUBROUTINE sort

! -----------------------------------------------------------

	SUBROUTINE sort2(arr,slave)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : indexx
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave
	INTEGER(I4B) :: ndum
	INTEGER(I4B), DIMENSION(size(arr)) :: index
	ndum=assert_eq(size(arr),size(slave),'sort2')
	call indexx(arr,index)
	arr=arr(index)
	slave=slave(index)
	END SUBROUTINE sort2

! -----------------------------------------------------------

	SUBROUTINE sort3(arr,slave1,slave2)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : indexx
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave1,slave2
	INTEGER(I4B) :: ndum
	INTEGER(I4B), DIMENSION(size(arr)) :: index
	ndum=assert_eq(size(arr),size(slave1),size(slave2),'sort3')
	call indexx(arr,index)
	arr=arr(index)
	slave1=slave1(index)
	slave2=slave2(index)
	END SUBROUTINE sort3

! -----------------------------------------------------------

	RECURSIVE SUBROUTINE sort_bypack(arr)
	use mo_nrtype
          use mo_nrutil, ONLY : array_copy,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	REAL(SP) :: a
	INTEGER(I4B) :: n,k,nl,nerr
	INTEGER(I4B), SAVE :: level=0
	LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: mask
	REAL(SP), DIMENSION(:), ALLOCATABLE, SAVE :: temp
	n=size(arr)
	if (n <= 1) RETURN
	k=(1+n)/2
	call swap(arr(1),arr(k),arr(1)>arr(k))
	call swap(arr(k),arr(n),arr(k)>arr(n))
	call swap(arr(1),arr(k),arr(1)>arr(k))
	if (n <= 3) RETURN
	level=level+1
	if (level == 1) allocate(mask(n),temp(n))
	a=arr(k)
	mask(1:n) = (arr <= a)
	mask(k) = .false.
	call array_copy(pack(arr,mask(1:n)),temp,nl,nerr)
	mask(k) = .true.
	temp(nl+2:n)=pack(arr,.not. mask(1:n))
	temp(nl+1)=a
	arr=temp(1:n)
	call sort_bypack(arr(1:nl))
	call sort_bypack(arr(nl+2:n))
	if (level == 1) deallocate(mask,temp)
	level=level-1
	END SUBROUTINE sort_bypack

! -----------------------------------------------------------

	SUBROUTINE sort_byreshape(arr)
	use mo_nrtype
          use mo_nrutil, ONLY : swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	REAL(SP), DIMENSION(:,:), ALLOCATABLE :: tab
	REAL(SP), PARAMETER :: big=huge(arr)
	INTEGER(I4B) :: inc,n,m
	n=size(arr)
	inc=1
	do
		inc=2*inc+1
		if (inc > n) exit
	end do
	do
		inc=inc/2
		m=(n+inc-1)/inc
		allocate(tab(inc,m))
		tab=reshape(arr, (/inc,m/) , (/big/) )
		do
			call swap(tab(:,1:m-1:2),tab(:,2:m:2), &
				tab(:,1:m-1:2)>tab(:,2:m:2))
			call swap(tab(:,2:m-1:2),tab(:,3:m:2), &
				tab(:,2:m-1:2)>tab(:,3:m:2))
			if (all(tab(:,1:m-1) <= tab(:,2:m))) exit
		end do
		arr=reshape(tab,shape(arr))
		deallocate(tab)
		if (inc <= 1) exit
	end do
	END SUBROUTINE sort_byreshape

! -----------------------------------------------------------

	SUBROUTINE sort_heap(arr)
	use mo_nrtype
	use mo_nrutil, ONLY : swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	INTEGER(I4B) :: i,n
	n=size(arr)
	do i=n/2,1,-1
		call sift_down(i,n)
	end do
	do i=n,2,-1
		call swap(arr(1),arr(i))
		call sift_down(1,i-1)
	end do
	CONTAINS
!BL
	SUBROUTINE sift_down(l,r)
	INTEGER(I4B), INTENT(IN) :: l,r
	INTEGER(I4B) :: j,jold
	REAL(SP) :: a
	a=arr(l)
	jold=l
	j=l+l
	do
		if (j > r) exit
		if (j < r) then
			if (arr(j) < arr(j+1)) j=j+1
		end if
		if (a >= arr(j)) exit
		arr(jold)=arr(j)
		jold=j
		j=j+j
	end do
	arr(jold)=a
	END SUBROUTINE sift_down
	END SUBROUTINE sort_heap

! -----------------------------------------------------------

	SUBROUTINE sort_pick(arr)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	INTEGER(I4B) :: i,j,n
	REAL(SP) :: a
	n=size(arr)
	do j=2,n
		a=arr(j)
		do i=j-1,1,-1
			if (arr(i) <= a) exit
			arr(i+1)=arr(i)
		end do
		arr(i+1)=a
	end do
	END SUBROUTINE sort_pick

! -----------------------------------------------------------

	SUBROUTINE sort_radix(arr)
	use mo_nrtype
          use mo_nrutil, ONLY : array_copy,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	INTEGER(I4B), DIMENSION(size(arr)) :: narr,temp
	LOGICAL, DIMENSION(size(arr)) :: msk
	INTEGER(I4B) :: k,negm,ib,ia,n,nl,nerr
	ib=bit_size(narr)
	ia=ceiling(log(real(maxexponent(arr)-minexponent(arr),sp))/log(2.0_sp)) &
		+ digits(arr)
	if (ib /= ia) call nrerror('sort_radix: bit sizes not compatible')
	negm=not(ishftc(1,-1))
	n=size(arr)
	narr=transfer(arr,narr,n)
	where (btest(narr,ib-1)) narr=ieor(narr,negm)
	do k=0,ib-2
		msk=btest(narr,k)
		call array_copy(pack(narr,.not. msk),temp,nl,nerr)
		temp(nl+1:n)=pack(narr,msk)
		narr=temp
	end do
	msk=btest(narr,ib-1)
	call array_copy(pack(narr,msk),temp,nl,nerr)
	temp(nl+1:n)=pack(narr,.not. msk)
	narr=temp
	where (btest(narr,ib-1)) narr=ieor(narr,negm)
	arr=transfer(narr,arr,n)
	END SUBROUTINE sort_radix

! -----------------------------------------------------------

	SUBROUTINE sort_shell(arr)
	use mo_nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	INTEGER(I4B) :: i,j,inc,n
	REAL(SP) :: v
	n=size(arr)
	inc=1
	do
		inc=3*inc+1
		if (inc > n) exit
	end do
	do
		inc=inc/3
		do i=inc+1,n
			v=arr(i)
			j=i
			do
				if (arr(j-inc) <= v) exit
				arr(j)=arr(j-inc)
				j=j-inc
				if (j <= inc) exit
			end do
			arr(j)=v
		end do
		if (inc <= 1) exit
	end do
	END SUBROUTINE sort_shell

! -----------------------------------------------------------

	SUBROUTINE spctrm(p,k,ovrlap,unit,n_window)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,nrerror
	use mo_nr, ONLY : four1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: p
	INTEGER(I4B), INTENT(IN) :: k
	LOGICAL(LGT), INTENT(IN) :: ovrlap
	INTEGER(I4B), OPTIONAL, INTENT(IN) :: n_window,unit
	INTEGER(I4B) :: j,joff,joffn,kk,m,m4,m43,m44,mm,iunit,nn_window
	REAL(SP) :: den,facm,facp,sumw
	REAL(SP), DIMENSION(2*size(p)) :: w
	REAL(SP), DIMENSION(4*size(p)) :: w1
	REAL(SP), DIMENSION(size(p)) :: w2
	COMPLEX(SPC), DIMENSION(2*size(p)) :: cw1
	m=size(p)
	if (present(n_window)) then
		nn_window=n_window
	else
		nn_window=1
	end if
	if (present(unit)) then
		iunit=unit
	else
		iunit=9
	end if
	mm=m+m
	m4=mm+mm
	m44=m4+4
	m43=m4+3
	den=0.0
	facm=m
	facp=1.0_sp/m
	w1(1:mm)=window(arth(1,1,mm),facm,facp,nn_window)
	sumw=dot_product(w1(1:mm),w1(1:mm))
	p(:)=0.0
	if (ovrlap) read (iunit,*) (w2(j),j=1,m)
	do kk=1,k
		do joff=-1,0,1
			if (ovrlap) then
				w1(joff+2:joff+mm:2)=w2(1:m)
				read (iunit,*) (w2(j),j=1,m)
				joffn=joff+mm
				w1(joffn+2:joffn+mm:2)=w2(1:m)
			else
				read (iunit,*) (w1(j),j=joff+2,m4,2)
			end if
		end do
		w=window(arth(1,1,mm),facm,facp,nn_window)
		w1(2:m4:2)=w1(2:m4:2)*w
		w1(1:m4:2)=w1(1:m4:2)*w
		cw1(1:mm)=cmplx(w1(1:m4:2),w1(2:m4:2),kind=spc)
		call four1(cw1(1:mm),1)
		w1(1:m4:2)=real(cw1(1:mm))
		w1(2:m4:2)=aimag(cw1(1:mm))
		p(1)=p(1)+w1(1)**2+w1(2)**2
		p(2:m)=p(2:m)+w1(4:2*m:2)**2+w1(3:2*m-1:2)**2+&
			w1(m44-4:m44-2*m:-2)**2+w1(m43-4:m43-2*m:-2)**2
		den=den+sumw
	end do
	p(:)=p(:)/(m4*den)
	CONTAINS
!BL
	FUNCTION window(j,facm,facp,nn_window)
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: j
	INTEGER(I4B), INTENT(IN) :: nn_window
	REAL(SP), INTENT(IN) :: facm,facp
	REAL(SP), DIMENSION(size(j)) :: window
	select case(nn_window)
		case(1)
			window(j)=(1.0_sp-abs(((j-1)-facm)*facp))
		case(2)
			window(j)=1.0
		case(3)
			window(j)=(1.0_sp-(((j-1)-facm)*facp)**2)
		case default
			call nrerror('unimplemented window function in spctrm')
	end select
	END FUNCTION window
	END SUBROUTINE spctrm

! -----------------------------------------------------------

	SUBROUTINE spear(data1,data2,d,zd,probd,rs,probrs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : betai,erfcc,sort2
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	REAL(SP), INTENT(OUT) :: d,zd,probd,rs,probrs
	INTEGER(I4B) :: n
	REAL(SP) :: aved,df,en,en3n,fac,sf,sg,t,vard
	REAL(SP), DIMENSION(size(data1)) :: wksp1,wksp2
	n=assert_eq(size(data1),size(data2),'spear')
	wksp1(:)=data1(:)
	wksp2(:)=data2(:)
	call sort2(wksp1,wksp2)
	call crank(wksp1,sf)
	call sort2(wksp2,wksp1)
	call crank(wksp2,sg)
	wksp1(:)=wksp1(:)-wksp2(:)
	d=dot_product(wksp1,wksp1)
	en=n
	en3n=en**3-en
	aved=en3n/6.0_sp-(sf+sg)/12.0_sp
	fac=(1.0_sp-sf/en3n)*(1.0_sp-sg/en3n)
	vard=((en-1.0_sp)*en**2*(en+1.0_sp)**2/36.0_sp)*fac
	zd=(d-aved)/sqrt(vard)
	probd=erfcc(abs(zd)/SQRT2)
	rs=(1.0_sp-(6.0_sp/en3n)*(d+(sf+sg)/12.0_sp))/sqrt(fac)
	fac=(1.0_sp+rs)*(1.0_sp-rs)
	if (fac > 0.0) then
		t=rs*sqrt((en-2.0_sp)/fac)
		df=en-2.0_sp
		probrs=betai(0.5_sp*df,0.5_sp,df/(df+t**2))
	else
		probrs=0.0
	end if
	CONTAINS
!BL
	SUBROUTINE crank(w,s)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,array_copy
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: s
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: w
	INTEGER(I4B) :: i,n,ndum,nties
	INTEGER(I4B), DIMENSION(size(w)) :: tstart,tend,tie,idx
	n=size(w)
	idx(:)=arth(1,1,n)
	tie(:)=merge(1,0,w == eoshift(w,-1))
	tie(1)=0
	w(:)=idx(:)
	if (all(tie == 0)) then
		s=0.0
		RETURN
	end if
	call array_copy(pack(idx(:),tie(:)<eoshift(tie(:),1)),tstart,nties,ndum)
	tend(1:nties)=pack(idx(:),tie(:)>eoshift(tie(:),1))
	do i=1,nties
		w(tstart(i):tend(i))=(tstart(i)+tend(i))/2.0_sp
	end do
	tend(1:nties)=tend(1:nties)-tstart(1:nties)+1
	s=sum(tend(1:nties)**3-tend(1:nties))
	END SUBROUTINE crank
	END SUBROUTINE spear

! -----------------------------------------------------------

	SUBROUTINE sphbes_s(n,x,sj,sy,sjp,syp)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessjy
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: x
	REAL(SP), INTENT(OUT) :: sj,sy,sjp,syp
	REAL(SP), PARAMETER :: RTPIO2=1.253314137315500_sp
	REAL(SP) :: factor,order,rj,rjp,ry,ryp
	call assert(n >= 0, x > 0.0, 'sphbes_s args')
	order=n+0.5_sp
	call bessjy(x,order,rj,ry,rjp,ryp)
	factor=RTPIO2/sqrt(x)
	sj=factor*rj
	sy=factor*ry
	sjp=factor*rjp-sj/(2.0_sp*x)
	syp=factor*ryp-sy/(2.0_sp*x)
	END SUBROUTINE sphbes_s


	SUBROUTINE sphbes_v(n,x,sj,sy,sjp,syp)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	use mo_nr, ONLY : bessjy
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(OUT) :: sj,sy,sjp,syp
	REAL(SP), PARAMETER :: RTPIO2=1.253314137315500_sp
	REAL(SP) :: order
	REAL(SP), DIMENSION(size(x)) :: factor,rj,rjp,ry,ryp
	call assert(n >= 0,  all(x > 0.0), 'sphbes_v args')
	order=n+0.5_sp
	call bessjy(x,order,rj,ry,rjp,ryp)
	factor=RTPIO2/sqrt(x)
	sj=factor*rj
	sy=factor*ry
	sjp=factor*rjp-sj/(2.0_sp*x)
	syp=factor*ryp-sy/(2.0_sp*x)
	END SUBROUTINE sphbes_v

! -----------------------------------------------------------

MODULE sphfpt_data
	use mo_nrtype
	INTEGER(I4B) :: m,n
	REAL(SP) :: c2,dx,gamma
END MODULE sphfpt_data

MODULE sphfpt_caller
	use mo_nrtype
	INTEGER(I4B) :: nn2
	REAL(SP) :: x1,x2,xf
END MODULE sphfpt_caller

	PROGRAM sphfpt
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	use mo_nr, ONLY : newt
	USE sphfpt_data
	USE sphfpt_caller
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N1=2,N2=1,NTOT=N1+N2
	REAL(SP), PARAMETER :: DXX=1.0e-4_sp
	REAL(SP), DIMENSION(:), POINTER :: v1,v2
	REAL(SP), DIMENSION(NTOT), TARGET :: v
	LOGICAL(LGT) :: check
	v1=>v(1:N2)
	v2=>v(N2+1:NTOT)
	nn2=N2
	dx=DXX
	do
		write(*,*) 'input m,n,c-squared (999 to end)'
		read(*,*) m,n,c2
		if (c2 == 999.0) exit
		if ((n < m) .or. (m < 0)) cycle
		gamma=(-0.5_sp)**m*product(&
			arth(n+1,1,m)*(arth(real(n,sp),-1.0_sp,m)/arth(1,1,m)))
		v1(1)=n*(n+1)-m*(m+1)+c2/2.0_sp
		v2(2)=v1(1)
		v2(1)=gamma*(1.0_sp-(v2(2)-c2)*dx/(2*(m+1)))
		x1=-1.0_sp+dx
		x2=1.0_sp-dx
		xf=0.0
		call newt(v,check)
		if (check) then
			write(*,*) 'shootf failed; bad initial guess'
			exit
		else
			write(*,'(1x,t6,a)') 'mu(m,n)'
			write(*,'(1x,f12.6)') v1(1)
		end if
	end do
	END PROGRAM sphfpt

	SUBROUTINE load1(x1,v1,y)
	use mo_nrtype
	USE sphfpt_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1
	REAL(SP), DIMENSION(:), INTENT(IN) :: v1
	REAL(SP), DIMENSION(:), INTENT(OUT) :: y
	REAL(SP) :: y1
	y(3)=v1(1)
	y1=merge(gamma,-gamma,mod(n-m,2) == 0)
	y(2)=-(y(3)-c2)*y1/(2*(m+1))
	y(1)=y1+y(2)*dx
	END SUBROUTINE load1

	SUBROUTINE load2(x2,v2,y)
	use mo_nrtype
	USE sphfpt_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x2
	REAL(SP), DIMENSION(:), INTENT(IN) :: v2
	REAL(SP), DIMENSION(:), INTENT(OUT) :: y
	y(3)=v2(2)
	y(1)=v2(1)
	y(2)=(y(3)-c2)*y(1)/(2*(m+1))
	END SUBROUTINE load2

	SUBROUTINE score(xf,y,f)
	use mo_nrtype
	USE sphfpt_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: xf
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: f
	f(1:3)=y(1:3)
	END SUBROUTINE score

! -----------------------------------------------------------

MODULE sphoot_data
	use mo_nrtype
	INTEGER(I4B) :: m,n
	REAL(SP) :: c2,dx,gamma
END MODULE sphoot_data

MODULE sphoot_caller
	use mo_nrtype
	INTEGER(I4B) :: nvar
	REAL(SP) :: x1,x2
END MODULE sphoot_caller

	PROGRAM sphoot
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	use mo_nr, ONLY : newt
	USE sphoot_data
	USE sphoot_caller
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NV=3,N2=1
	REAL(SP), DIMENSION(N2) :: v
	LOGICAL(LGT) :: check
	nvar=NV
	dx=1.0e-4_sp
	do
		write(*,*) 'input m,n,c-squared (999 to end)'
		read(*,*) m,n,c2
		if (c2 == 999.0) exit
		if ((n < m) .or. (m < 0)) cycle
		gamma=(-0.5_sp)**m*product(&
			arth(n+1,1,m)*(arth(real(n,sp),-1.0_sp,m)/arth(1,1,m)))
		v(1)=n*(n+1)-m*(m+1)+c2/2.0_sp
		x1=-1.0_sp+dx
		x2=0.0
		call newt(v,check)
		if (check) then
			write(*,*)'shoot failed; bad initial guess'
			exit
		else
			write(*,'(1x,t6,a)') 'mu(m,n)'
			write(*,'(1x,f12.6)') v(1)
		end if
	end do
	END PROGRAM sphoot

	SUBROUTINE load(x1,v,y)
	use mo_nrtype
	USE sphoot_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1
	REAL(SP), DIMENSION(:), INTENT(IN) :: v
	REAL(SP), DIMENSION(:), INTENT(OUT) :: y
	REAL(SP) :: y1
	y(3)=v(1)
	y1=merge(gamma,-gamma, mod(n-m,2) == 0 )
	y(2)=-(y(3)-c2)*y1/(2*(m+1))
	y(1)=y1+y(2)*dx
	END SUBROUTINE load

	SUBROUTINE score(x2,y,f)
	use mo_nrtype
	USE sphoot_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x2
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: f
	f(1)=merge(y(2),y(1), mod(n-m,2) == 0 )
	END SUBROUTINE score

	SUBROUTINE derivs(x,y,dydx)
	use mo_nrtype
	USE sphoot_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
	dydx(1)=y(2)
	dydx(2)=(2.0_sp*x*(m+1.0_sp)*y(2)-(y(3)-c2*x*x)*y(1))/(1.0_sp-x*x)
	dydx(3)=0.0
	END SUBROUTINE derivs

! -----------------------------------------------------------

	SUBROUTINE splie2(x1a,x2a,ya,y2a)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : spline
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: y2a
	INTEGER(I4B) :: j,m,ndum
	m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splie2: m')
	ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splie2: ndum')
	do j=1,m
		call spline(x2a,ya(j,:),1.0e30_sp,1.0e30_sp,y2a(j,:))
	end do
	END SUBROUTINE splie2

! -----------------------------------------------------------

	FUNCTION splin2(x1a,x2a,ya,y2a,x1,x2)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : spline,splint
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x1a,x2a
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: ya,y2a
	REAL(SP), INTENT(IN) :: x1,x2
	REAL(SP) :: splin2
	INTEGER(I4B) :: j,m,ndum
	REAL(SP), DIMENSION(size(x1a)) :: yytmp,y2tmp2
	m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splin2: m')
	ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splin2: ndum')
	do j=1,m
		yytmp(j)=splint(x2a,ya(j,:),y2a(j,:),x2)
	end do
	call spline(x1a,yytmp,1.0e30_sp,1.0e30_sp,y2tmp2)
	splin2=splint(x1a,yytmp,y2tmp2,x1)
	END FUNCTION splin2

! -----------------------------------------------------------

	SUBROUTINE spline(x,y,yp1,ypn,y2)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : tridag
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	REAL(SP), INTENT(IN) :: yp1,ypn
	REAL(SP), DIMENSION(:), INTENT(OUT) :: y2
	INTEGER(I4B) :: n
	REAL(SP), DIMENSION(size(x)) :: a,b,c,r
	n=assert_eq(size(x),size(y),size(y2),'spline')
	c(1:n-1)=x(2:n)-x(1:n-1)
	r(1:n-1)=6.0_sp*((y(2:n)-y(1:n-1))/c(1:n-1))
	r(2:n-1)=r(2:n-1)-r(1:n-2)
	a(2:n-1)=c(1:n-2)
	b(2:n-1)=2.0_sp*(c(2:n-1)+a(2:n-1))
	b(1)=1.0
	b(n)=1.0
	if (yp1 > 0.99e30_sp) then
		r(1)=0.0
		c(1)=0.0
	else
		r(1)=(3.0_sp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
		c(1)=0.5
	end if
	if (ypn > 0.99e30_sp) then
		r(n)=0.0
		a(n)=0.0
	else
		r(n)=(-3.0_sp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
		a(n)=0.5
	end if
	call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
	END SUBROUTINE spline

! -----------------------------------------------------------

	FUNCTION splint(xa,ya,y2a,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY: locate
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: splint
	INTEGER(I4B) :: khi,klo,n
	REAL(SP) :: a,b,h
	n=assert_eq(size(xa),size(ya),size(y2a),'splint')
	klo=max(min(locate(xa,x),n-1),1)
	khi=klo+1
	h=xa(khi)-xa(klo)
	if (h == 0.0) call nrerror('bad xa input in splint')
	a=(xa(khi)-x)/h
	b=(x-xa(klo))/h
	splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_sp
	END FUNCTION splint

! -----------------------------------------------------------

	SUBROUTINE sprsax_sp(sa,x,b)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,scatter_add
	IMPLICIT NONE
	TYPE(sprs2_sp), INTENT(IN) :: sa
	REAL(SP), DIMENSION (:), INTENT(IN) :: x
	REAL(SP), DIMENSION (:), INTENT(OUT) :: b
	INTEGER(I4B) :: ndum
	ndum=assert_eq(sa%n,size(x),size(b),'sprsax_sp')
	b=0.0_sp
	call scatter_add(b,sa%val*x(sa%jcol),sa%irow)
	END SUBROUTINE sprsax_sp

	SUBROUTINE sprsax_dp(sa,x,b)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,scatter_add
	IMPLICIT NONE
	TYPE(sprs2_dp), INTENT(IN) :: sa
	REAL(DP), DIMENSION (:), INTENT(IN) :: x
	REAL(DP), DIMENSION (:), INTENT(OUT) :: b
	INTEGER(I4B) :: ndum
	ndum=assert_eq(sa%n,size(x),size(b),'sprsax_dp')
	b=0.0_dp
	call scatter_add(b,sa%val*x(sa%jcol),sa%irow)
	END SUBROUTINE sprsax_dp

! -----------------------------------------------------------

	SUBROUTINE sprsdiag_sp(sa,b)
	use mo_nrtype
          use mo_nrutil, ONLY : array_copy,assert_eq
	IMPLICIT NONE
	TYPE(sprs2_sp), INTENT(IN) :: sa
	REAL(SP), DIMENSION(:), INTENT(OUT) :: b
	REAL(SP), DIMENSION(size(b)) :: val
	INTEGER(I4B) :: k,l,ndum,nerr
	INTEGER(I4B), DIMENSION(size(b)) :: i
	LOGICAL(LGT), DIMENSION(:), ALLOCATABLE :: mask
	ndum=assert_eq(sa%n,size(b),'sprsdiag_sp')
	l=sa%len
	allocate(mask(l))
	mask = (sa%irow(1:l) == sa%jcol(1:l))
	call array_copy(pack(sa%val(1:l),mask),val,k,nerr)
	i(1:k)=pack(sa%irow(1:l),mask)
	deallocate(mask)
	b=0.0
	b(i(1:k))=val(1:k)
	END SUBROUTINE sprsdiag_sp

	SUBROUTINE sprsdiag_dp(sa,b)
	use mo_nrtype
          use mo_nrutil, ONLY : array_copy,assert_eq
	IMPLICIT NONE
	TYPE(sprs2_dp), INTENT(IN) :: sa
	REAL(DP), DIMENSION(:), INTENT(OUT) :: b
	REAL(DP), DIMENSION(size(b)) :: val
	INTEGER(I4B) :: k,l,ndum,nerr
	INTEGER(I4B), DIMENSION(size(b)) :: i
	LOGICAL(LGT), DIMENSION(:), ALLOCATABLE :: mask
	ndum=assert_eq(sa%n,size(b),'sprsdiag_dp')
	l=sa%len
	allocate(mask(l))
	mask = (sa%irow(1:l) == sa%jcol(1:l))
	call array_copy(pack(sa%val(1:l),mask),val,k,nerr)
	i(1:k)=pack(sa%irow(1:l),mask)
	deallocate(mask)
	b=0.0
	b(i(1:k))=val(1:k)
	END SUBROUTINE sprsdiag_dp

! -----------------------------------------------------------

	SUBROUTINE sprsin_sp(a,thresh,sa)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	REAL(SP), INTENT(IN) :: thresh
	TYPE(sprs2_sp), INTENT(OUT) :: sa
	INTEGER(I4B) :: n,len
	LOGICAL(LGT), DIMENSION(size(a,1),size(a,2)) :: mask
	n=assert_eq(size(a,1),size(a,2),'sprsin_sp')
	mask=abs(a)>thresh
	len=count(mask)
	allocate(sa%val(len),sa%irow(len),sa%jcol(len))
	sa%n=n
	sa%len=len
	sa%val=pack(a,mask)
	sa%irow=pack(spread(arth(1,1,n),2,n),mask)
	sa%jcol=pack(spread(arth(1,1,n),1,n),mask)
	END SUBROUTINE sprsin_sp

	SUBROUTINE sprsin_dp(a,thresh,sa)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
	REAL(DP), INTENT(IN) :: thresh
	TYPE(sprs2_dp), INTENT(OUT) :: sa
	INTEGER(I4B) :: n,len
	LOGICAL(LGT), DIMENSION(size(a,1),size(a,2)) :: mask
	n=assert_eq(size(a,1),size(a,2),'sprsin_dp')
	mask=abs(a)>thresh
	len=count(mask)
	allocate(sa%val(len),sa%irow(len),sa%jcol(len))
	sa%n=n
	sa%len=len
	sa%val=pack(a,mask)
	sa%irow=pack(spread(arth(1,1,n),2,n),mask)
	sa%jcol=pack(spread(arth(1,1,n),1,n),mask)
	END SUBROUTINE sprsin_dp

! -----------------------------------------------------------

	SUBROUTINE sprstp(sa)
	use mo_nrtype
	IMPLICIT NONE
	TYPE(sprs2_sp), INTENT(INOUT) :: sa
	INTEGER(I4B), DIMENSION(:), POINTER :: temp
	temp=>sa%irow
	sa%irow=>sa%jcol
	sa%jcol=>temp
	END SUBROUTINE sprstp

! -----------------------------------------------------------

	SUBROUTINE sprstx_sp(sa,x,b)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,scatter_add
	IMPLICIT NONE
	TYPE(sprs2_sp), INTENT(IN) :: sa
	REAL(SP), DIMENSION (:), INTENT(IN) :: x
	REAL(SP), DIMENSION (:), INTENT(OUT) :: b
	INTEGER(I4B) :: ndum
	ndum=assert_eq(sa%n,size(x),size(b),'sprstx_sp')
	b=0.0_sp
	call scatter_add(b,sa%val*x(sa%irow),sa%jcol)
	END SUBROUTINE sprstx_sp

	SUBROUTINE sprstx_dp(sa,x,b)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,scatter_add
	IMPLICIT NONE
	TYPE(sprs2_dp), INTENT(IN) :: sa
	REAL(DP), DIMENSION (:), INTENT(IN) :: x
	REAL(DP), DIMENSION (:), INTENT(OUT) :: b
	INTEGER(I4B) :: ndum
	ndum=assert_eq(sa%n,size(x),size(b),'sprstx_dp')
	b=0.0_dp
	call scatter_add(b,sa%val*x(sa%irow),sa%jcol)
	END SUBROUTINE sprstx_dp

! -----------------------------------------------------------

	SUBROUTINE stifbs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert_eq,cumsum,iminloc,nrerror,&
		outerdiff,outerprod,upper_triangle
	use mo_nr, ONLY : simpr,pzextr
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
	REAL(SP), INTENT(IN) :: htry,eps
	REAL(SP), INTENT(INOUT) :: x
	REAL(SP), INTENT(OUT) :: hdid,hnext
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
!BL
		SUBROUTINE jacobn(x,y,dfdx,dfdy)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dfdx
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dfdy
		END SUBROUTINE jacobn
	END INTERFACE
	INTEGER(I4B), PARAMETER :: IMAX=8, KMAXX=IMAX-1
	REAL(SP), PARAMETER :: SAFE1=0.25_sp,SAFE2=0.7_sp,REDMAX=1.0e-5_sp,&
		REDMIN=0.7_sp,TINY=1.0e-30_sp,SCALMX=0.1_sp
	INTEGER(I4B) :: k,km,ndum
	INTEGER(I4B), DIMENSION(IMAX) :: nseq = (/ 2,6,10,14,22,34,50,70 /)
	INTEGER(I4B), SAVE :: kopt,kmax,nvold=-1
	REAL(SP), DIMENSION(KMAXX,KMAXX), SAVE :: alf
	REAL(SP), DIMENSION(KMAXX) :: err
	REAL(SP), DIMENSION(IMAX), SAVE :: a
	REAL(SP), SAVE :: epsold = -1.0
	REAL(SP) :: eps1,errmax,fact,h,red,scale,wrkmin,xest
	REAL(SP), SAVE :: xnew
	REAL(SP), DIMENSION(size(y)) :: dfdx,yerr,ysav,yseq
	REAL(SP), DIMENSION(size(y),size(y)) :: dfdy
	LOGICAL(LGT) :: reduct
	LOGICAL(LGT), SAVE :: first=.true.
	ndum=assert_eq(size(y),size(dydx),size(yscal),'stifbs')
	if (eps /= epsold .or. nvold /= size(y)) then
		hnext=-1.0e29_sp
		xnew=-1.0e29_sp
		eps1=SAFE1*eps
		a(:)=cumsum(nseq,1)
		where (upper_triangle(KMAXX,KMAXX)) alf=eps1** &
			(outerdiff(a(2:),a(2:))/outerprod(arth( &
			3.0_sp,2.0_sp,KMAXX),(a(2:)-a(1)+1.0_sp)))
		epsold=eps
		nvold=size(y)
		a(:)=cumsum(nseq,1+nvold)
		do kopt=2,KMAXX-1
			if (a(kopt+1) > a(kopt)*alf(kopt-1,kopt)) exit
		end do
		kmax=kopt
	end if
	h=htry
	ysav(:)=y(:)
	call jacobn(x,y,dfdx,dfdy)
	if (h /= hnext .or. x /= xnew) then
		first=.true.
		kopt=kmax
	end if
	reduct=.false.
	main_loop: do
		do k=1,kmax
			xnew=x+h
			if (xnew == x) call nrerror('step size underflow in stifbs')
			call simpr(ysav,dydx,dfdx,dfdy,x,h,nseq(k),yseq,derivs)
			xest=(h/nseq(k))**2
			call pzextr(k,xest,yseq,y,yerr)
			if (k /= 1) then
				errmax=maxval(abs(yerr(:)/yscal(:)))
				errmax=max(TINY,errmax)/eps
				km=k-1
				err(km)=(errmax/SAFE1)**(1.0_sp/(2*km+1))
			end if
			if (k /= 1 .and. (k >= kopt-1 .or. first)) then
				if (errmax < 1.0) exit main_loop
				if (k == kmax .or. k == kopt+1) then
					red=SAFE2/err(km)
					exit
				else if (k == kopt) then
					if (alf(kopt-1,kopt) < err(km)) then
						red=1.0_sp/err(km)
						exit
					end if
				else if (kopt == kmax) then
					if (alf(km,kmax-1) < err(km)) then
						red=alf(km,kmax-1)*SAFE2/err(km)
						exit
					end if
				else if (alf(km,kopt) < err(km)) then
					red=alf(km,kopt-1)/err(km)
					exit
				end if
			end if
		end do
		red=max(min(red,REDMIN),REDMAX)
		h=h*red
		reduct=.true.
	end do main_loop
	x=xnew
	hdid=h
	first=.false.
	kopt=1+iminloc(a(2:km+1)*max(err(1:km),SCALMX))
	scale=max(err(kopt-1),SCALMX)
	wrkmin=scale*a(kopt)
	hnext=h/scale
	if (kopt >= k .and. kopt /= kmax .and. .not. reduct) then
		fact=max(scale/alf(kopt-1,kopt),SCALMX)
		if (a(kopt+1)*fact <= wrkmin) then
			hnext=h/fact
			kopt=kopt+1
		end if
	end if
	END SUBROUTINE stifbs

! -----------------------------------------------------------

	SUBROUTINE stiff(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,diagadd,nrerror
	use mo_nr, ONLY : lubksb,ludcmp
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	REAL(SP), DIMENSION(:), INTENT(IN) :: dydx,yscal
	REAL(SP), INTENT(INOUT) :: x
	REAL(SP), INTENT(IN) :: htry,eps
	REAL(SP), INTENT(OUT) :: hdid,hnext
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
!BL
		SUBROUTINE jacobn(x,y,dfdx,dfdy)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dfdx
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dfdy
		END SUBROUTINE jacobn
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXTRY=40
	REAL(SP), PARAMETER :: SAFETY=0.9_sp,GROW=1.5_sp,PGROW=-0.25_sp,&
		SHRNK=0.5_sp,PSHRNK=-1.0_sp/3.0_sp,ERRCON=0.1296_sp,&
		GAM=1.0_sp/2.0_sp,&
		A21=2.0_sp,A31=48.0_sp/25.0_sp,A32=6.0_sp/25.0_sp,C21=-8.0_sp,&
		C31=372.0_sp/25.0_sp,C32=12.0_sp/5.0_sp,&
		C41=-112.0_sp/125.0_sp,C42=-54.0_sp/125.0_sp,&
		C43=-2.0_sp/5.0_sp,B1=19.0_sp/9.0_sp,B2=1.0_sp/2.0_sp,&
		B3=25.0_sp/108.0_sp,B4=125.0_sp/108.0_sp,E1=17.0_sp/54.0_sp,&
		E2=7.0_sp/36.0_sp,E3=0.0_sp,E4=125.0_sp/108.0_sp,&
		C1X=1.0_sp/2.0_sp,C2X=-3.0_sp/2.0_sp,C3X=121.0_sp/50.0_sp,&
		C4X=29.0_sp/250.0_sp,A2X=1.0_sp,A3X=3.0_sp/5.0_sp
	INTEGER(I4B) :: jtry,ndum
	INTEGER(I4B), DIMENSION(size(y)) :: indx
	REAL(SP), DIMENSION(size(y)) :: dfdx,dytmp,err,g1,g2,g3,g4,ysav
	REAL(SP), DIMENSION(size(y),size(y)) :: a,dfdy
	REAL(SP) :: d,errmax,h,xsav
	ndum=assert_eq(size(y),size(dydx),size(yscal),'stiff')
	xsav=x
	ysav(:)=y(:)
	call jacobn(xsav,ysav,dfdx,dfdy)
	h=htry
	do jtry=1,MAXTRY
		a(:,:)=-dfdy(:,:)
		call diagadd(a,1.0_sp/(GAM*h))
		call ludcmp(a,indx,d)
		g1=dydx+h*C1X*dfdx
		call lubksb(a,indx,g1)
		y=ysav+A21*g1
		x=xsav+A2X*h
		call derivs(x,y,dytmp)
		g2=dytmp+h*C2X*dfdx+C21*g1/h
		call lubksb(a,indx,g2)
		y=ysav+A31*g1+A32*g2
		x=xsav+A3X*h
		call derivs(x,y,dytmp)
		g3=dytmp+h*C3X*dfdx+(C31*g1+C32*g2)/h
		call lubksb(a,indx,g3)
		g4=dytmp+h*C4X*dfdx+(C41*g1+C42*g2+C43*g3)/h
		call lubksb(a,indx,g4)
		y=ysav+B1*g1+B2*g2+B3*g3+B4*g4
		err=E1*g1+E2*g2+E3*g3+E4*g4
		x=xsav+h
		if (x == xsav) call &
			nrerror('stepsize not significant in stiff')
		errmax=maxval(abs(err/yscal))/eps
		if (errmax <= 1.0) then
			hdid=h
			hnext=merge(SAFETY*h*errmax**PGROW, GROW*h, &
				errmax > ERRCON)
			RETURN
		else
			hnext=SAFETY*h*errmax**PSHRNK
			h=sign(max(abs(hnext),SHRNK*abs(h)),h)
		end if
	end do
	call nrerror('exceeded MAXTRY in stiff')
	END SUBROUTINE stiff

! -----------------------------------------------------------

	SUBROUTINE stoerm(y,d2y,xs,htot,nstep,yout,derivs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: y,d2y
	REAL(SP), INTENT(IN) :: xs,htot
	INTEGER(I4B), INTENT(IN) :: nstep
	REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: neqn,neqn1,nn,nv
	REAL(SP) :: h,h2,halfh,x
	REAL(SP), DIMENSION(size(y)) :: ytemp
	nv=assert_eq(size(y),size(d2y),size(yout),'stoerm')
	neqn=nv/2
	neqn1=neqn+1
	h=htot/nstep
	halfh=0.5_sp*h
	ytemp(neqn1:nv)=h*(y(neqn1:nv)+halfh*d2y(1:neqn))
	ytemp(1:neqn)=y(1:neqn)+ytemp(neqn1:nv)
	x=xs+h
	call derivs(x,ytemp,yout)
	h2=h*h
	do nn=2,nstep
		ytemp(neqn1:nv)=ytemp(neqn1:nv)+h2*yout(1:neqn)
		ytemp(1:neqn)=ytemp(1:neqn)+ytemp(neqn1:nv)
		x=x+h
		call derivs(x,ytemp,yout)
	end do
	yout(neqn1:nv)=ytemp(neqn1:nv)/h+halfh*yout(1:neqn)
	yout(1:neqn)=ytemp(1:neqn)
	END SUBROUTINE stoerm

! -----------------------------------------------------------

	SUBROUTINE svbksb_sp(u,w,v,b,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: u,v
	REAL(SP), DIMENSION(:), INTENT(IN) :: w,b
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x
	INTEGER(I4B) :: mdum,ndum
	REAL(SP), DIMENSION(size(x)) :: tmp
	mdum=assert_eq(size(u,1),size(b),'svbksb_sp: mdum')
	ndum=assert_eq((/size(u,2),size(v,1),size(v,2),size(w),size(x)/),&
		'svbksb_sp: ndum')
	where (w /= 0.0)
		tmp=matmul(b,u)/w
	elsewhere
		tmp=0.0
	end where
	x=matmul(v,tmp)
	END SUBROUTINE svbksb_sp

	SUBROUTINE svbksb_dp(u,w,v,b,x)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,v
	REAL(DP), DIMENSION(:), INTENT(IN) :: w,b
	REAL(DP), DIMENSION(:), INTENT(OUT) :: x
	INTEGER(I4B) :: mdum,ndum
	REAL(DP), DIMENSION(size(x)) :: tmp
	mdum=assert_eq(size(u,1),size(b),'svbksb_dp: mdum')
	ndum=assert_eq((/size(u,2),size(v,1),size(v,2),size(w),size(x)/),&
		'svbksb_dp: ndum')
	where (w /= 0.0)
		tmp=matmul(b,u)/w
	elsewhere
		tmp=0.0
	end where
	x=matmul(v,tmp)
	END SUBROUTINE svbksb_dp

! -----------------------------------------------------------

	SUBROUTINE svdcmp_sp(a,w,v)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror,outerprod
	use mo_nr, ONLY : pythag
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: w
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
	INTEGER(I4B) :: i,its,j,k,l,m,n,nm
	REAL(SP) :: anorm,c,f,g,h,s,scale,x,y,z
	REAL(SP), DIMENSION(size(a,1)) :: tempm
	REAL(SP), DIMENSION(size(a,2)) :: rv1,tempn
	m=size(a,1)
	n=assert_eq(size(a,2),size(v,1),size(v,2),size(w),'svdcmp_sp')
	g=0.0
	scale=0.0
	do i=1,n
		l=i+1
		rv1(i)=scale*g
		g=0.0
		scale=0.0
		if (i <= m) then
			scale=sum(abs(a(i:m,i)))
			if (scale /= 0.0) then
				a(i:m,i)=a(i:m,i)/scale
				s=dot_product(a(i:m,i),a(i:m,i))
				f=a(i,i)
				g=-sign(sqrt(s),f)
				h=f*g-s
				a(i,i)=f-g
				tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
				a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
				a(i:m,i)=scale*a(i:m,i)
			end if
		end if
		w(i)=scale*g
		g=0.0
		scale=0.0
		if ((i <= m) .and. (i /= n)) then
			scale=sum(abs(a(i,l:n)))
			if (scale /= 0.0) then
				a(i,l:n)=a(i,l:n)/scale
				s=dot_product(a(i,l:n),a(i,l:n))
				f=a(i,l)
				g=-sign(sqrt(s),f)
				h=f*g-s
				a(i,l)=f-g
				rv1(l:n)=a(i,l:n)/h
				tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
				a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
				a(i,l:n)=scale*a(i,l:n)
			end if
		end if
	end do
	anorm=maxval(abs(w)+abs(rv1))
	do i=n,1,-1
		if (i < n) then
			if (g /= 0.0) then
				v(l:n,i)=(a(i,l:n)/a(i,l))/g
				tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
				v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
			end if
			v(i,l:n)=0.0
			v(l:n,i)=0.0
		end if
		v(i,i)=1.0
		g=rv1(i)
		l=i
	end do
	do i=min(m,n),1,-1
		l=i+1
		g=w(i)
		a(i,l:n)=0.0
		if (g /= 0.0) then
			g=1.0_sp/g
			tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
			a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
			a(i:m,i)=a(i:m,i)*g
		else
			a(i:m,i)=0.0
		end if
		a(i,i)=a(i,i)+1.0_sp
	end do
	do k=n,1,-1
		do its=1,30
			do l=k,1,-1
				nm=l-1
				if ((abs(rv1(l))+anorm) == anorm) exit
				if ((abs(w(nm))+anorm) == anorm) then
					c=0.0
					s=1.0
					do i=l,k
						f=s*rv1(i)
						rv1(i)=c*rv1(i)
						if ((abs(f)+anorm) == anorm) exit
						g=w(i)
						h=pythag(f,g)
						w(i)=h
						h=1.0_sp/h
						c= (g*h)
						s=-(f*h)
						tempm(1:m)=a(1:m,nm)
						a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
						a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
					end do
					exit
				end if
			end do
			z=w(k)
			if (l == k) then
				if (z < 0.0) then
					w(k)=-z
					v(1:n,k)=-v(1:n,k)
				end if
				exit
			end if
			if (its == 30) call nrerror('svdcmp_sp: no convergence in svdcmp')
			x=w(l)
			nm=k-1
			y=w(nm)
			g=rv1(nm)
			h=rv1(k)
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_sp*h*y)
			g=pythag(f,1.0_sp)
			f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
			c=1.0
			s=1.0
			do j=l,nm
				i=j+1
				g=rv1(i)
				y=w(i)
				h=s*g
				g=c*g
				z=pythag(f,h)
				rv1(j)=z
				c=f/z
				s=h/z
				f= (x*c)+(g*s)
				g=-(x*s)+(g*c)
				h=y*s
				y=y*c
				tempn(1:n)=v(1:n,j)
				v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
				v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
				z=pythag(f,h)
				w(j)=z
				if (z /= 0.0) then
					z=1.0_sp/z
					c=f*z
					s=h*z
				end if
				f= (c*g)+(s*y)
				x=-(s*g)+(c*y)
				tempm(1:m)=a(1:m,j)
				a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
				a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
			end do
			rv1(l)=0.0
			rv1(k)=f
			w(k)=x
		end do
	end do
	END SUBROUTINE svdcmp_sp

	SUBROUTINE svdcmp_dp(a,w,v)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror,outerprod
	use mo_nr, ONLY : pythag
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(DP), DIMENSION(:), INTENT(OUT) :: w
	REAL(DP), DIMENSION(:,:), INTENT(OUT) :: v
	INTEGER(I4B) :: i,its,j,k,l,m,n,nm
	REAL(DP) :: anorm,c,f,g,h,s,scale,x,y,z
	REAL(DP), DIMENSION(size(a,1)) :: tempm
	REAL(DP), DIMENSION(size(a,2)) :: rv1,tempn
	m=size(a,1)
	n=assert_eq(size(a,2),size(v,1),size(v,2),size(w),'svdcmp_dp')
	g=0.0
	scale=0.0
	do i=1,n
		l=i+1
		rv1(i)=scale*g
		g=0.0
		scale=0.0
		if (i <= m) then
			scale=sum(abs(a(i:m,i)))
			if (scale /= 0.0) then
				a(i:m,i)=a(i:m,i)/scale
				s=dot_product(a(i:m,i),a(i:m,i))
				f=a(i,i)
				g=-sign(sqrt(s),f)
				h=f*g-s
				a(i,i)=f-g
				tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
				a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
				a(i:m,i)=scale*a(i:m,i)
			end if
		end if
		w(i)=scale*g
		g=0.0
		scale=0.0
		if ((i <= m) .and. (i /= n)) then
			scale=sum(abs(a(i,l:n)))
			if (scale /= 0.0) then
				a(i,l:n)=a(i,l:n)/scale
				s=dot_product(a(i,l:n),a(i,l:n))
				f=a(i,l)
				g=-sign(sqrt(s),f)
				h=f*g-s
				a(i,l)=f-g
				rv1(l:n)=a(i,l:n)/h
				tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
				a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
				a(i,l:n)=scale*a(i,l:n)
			end if
		end if
	end do
	anorm=maxval(abs(w)+abs(rv1))
	do i=n,1,-1
		if (i < n) then
			if (g /= 0.0) then
				v(l:n,i)=(a(i,l:n)/a(i,l))/g
				tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
				v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
			end if
			v(i,l:n)=0.0
			v(l:n,i)=0.0
		end if
		v(i,i)=1.0
		g=rv1(i)
		l=i
	end do
	do i=min(m,n),1,-1
		l=i+1
		g=w(i)
		a(i,l:n)=0.0
		if (g /= 0.0) then
			g=1.0_dp/g
			tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
			a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
			a(i:m,i)=a(i:m,i)*g
		else
			a(i:m,i)=0.0
		end if
		a(i,i)=a(i,i)+1.0_dp
	end do
	do k=n,1,-1
		do its=1,30
			do l=k,1,-1
				nm=l-1
				if ((abs(rv1(l))+anorm) == anorm) exit
				if ((abs(w(nm))+anorm) == anorm) then
					c=0.0
					s=1.0
					do i=l,k
						f=s*rv1(i)
						rv1(i)=c*rv1(i)
						if ((abs(f)+anorm) == anorm) exit
						g=w(i)
						h=pythag(f,g)
						w(i)=h
						h=1.0_dp/h
						c= (g*h)
						s=-(f*h)
						tempm(1:m)=a(1:m,nm)
						a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
						a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
					end do
					exit
				end if
			end do
			z=w(k)
			if (l == k) then
				if (z < 0.0) then
					w(k)=-z
					v(1:n,k)=-v(1:n,k)
				end if
				exit
			end if
			if (its == 30) call nrerror('svdcmp_dp: no convergence in svdcmp')
			x=w(l)
			nm=k-1
			y=w(nm)
			g=rv1(nm)
			h=rv1(k)
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_dp*h*y)
			g=pythag(f,1.0_dp)
			f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
			c=1.0
			s=1.0
			do j=l,nm
				i=j+1
				g=rv1(i)
				y=w(i)
				h=s*g
				g=c*g
				z=pythag(f,h)
				rv1(j)=z
				c=f/z
				s=h/z
				f= (x*c)+(g*s)
				g=-(x*s)+(g*c)
				h=y*s
				y=y*c
				tempn(1:n)=v(1:n,j)
				v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
				v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
				z=pythag(f,h)
				w(j)=z
				if (z /= 0.0) then
					z=1.0_dp/z
					c=f*z
					s=h*z
				end if
				f= (c*g)+(s*y)
				x=-(s*g)+(c*y)
				tempm(1:m)=a(1:m,j)
				a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
				a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
			end do
			rv1(l)=0.0
			rv1(k)=f
			w(k)=x
		end do
	end do
	END SUBROUTINE svdcmp_dp

! -----------------------------------------------------------

	SUBROUTINE svdfit(x,y,sig,a,v,w,chisq,funcs)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,vabs
	use mo_nr, ONLY : svbksb,svdcmp
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
	REAL(SP), DIMENSION(:), INTENT(OUT) :: a,w
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
	REAL(SP), INTENT(OUT) :: chisq
	INTERFACE
		FUNCTION funcs(x,n)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		INTEGER(I4B), INTENT(IN) :: n
		REAL(SP), DIMENSION(n) :: funcs
		END FUNCTION funcs
	END INTERFACE
	REAL(SP), PARAMETER :: TOL=1.0e-5_sp
	INTEGER(I4B) :: i,ma,n
	REAL(SP), DIMENSION(size(x)) :: b,sigi
	REAL(SP), DIMENSION(size(x),size(a)) :: u,usav
	n=assert_eq(size(x),size(y),size(sig),'svdfit: n')
	ma=assert_eq(size(a),size(v,1),size(v,2),size(w),'svdfit: ma')
	sigi=1.0_sp/sig
	b=y*sigi
	do i=1,n
		usav(i,:)=funcs(x(i),ma)
	end do
	u=usav*spread(sigi,dim=2,ncopies=ma)
	usav=u
	call svdcmp(u,w,v)
	where (w < TOL*maxval(w)) w=0.0
	call svbksb(u,w,v,b,a)
	chisq=vabs(matmul(usav,a)-b)**2
	END SUBROUTINE svdfit

! -----------------------------------------------------------

	SUBROUTINE svdvar(v,w,cvm)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: v
	REAL(SP), DIMENSION(:), INTENT(IN) :: w
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: cvm
	INTEGER(I4B) :: ma
	REAL(SP), DIMENSION(size(w)) :: wti
	ma=assert_eq((/size(v,1),size(v,2),size(w),size(cvm,1),size(cvm,2)/),&
		'svdvar')
	where (w /= 0.0)
		wti=1.0_sp/(w*w)
	elsewhere
		wti=0.0
	end where
	cvm=v*spread(wti,dim=1,ncopies=ma)
	cvm=matmul(cvm,transpose(v))
	END SUBROUTINE svdvar

! -----------------------------------------------------------

	FUNCTION toeplz(r,y)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: r,y
	REAL(SP), DIMENSION(size(y)) :: toeplz
	INTEGER(I4B) :: m,m1,n,ndum
	REAL(SP) :: sd,sgd,sgn,shn,sxn
	REAL(SP), DIMENSION(size(y)) :: g,h,t
	n=size(y)
	ndum=assert_eq(2*n-1,size(r),'toeplz: ndum')
	if (r(n) == 0.0) call nrerror('toeplz: initial singular minor')
	toeplz(1)=y(1)/r(n)
	if (n == 1) RETURN
	g(1)=r(n-1)/r(n)
	h(1)=r(n+1)/r(n)
	do m=1,n
		m1=m+1
		sxn=-y(m1)+dot_product(r(n+1:n+m),toeplz(m:1:-1))
		sd=-r(n)+dot_product(r(n+1:n+m),g(1:m))
		if (sd == 0.0) exit
		toeplz(m1)=sxn/sd
		toeplz(1:m)=toeplz(1:m)-toeplz(m1)*g(m:1:-1)
		if (m1 == n) RETURN
		sgn=-r(n-m1)+dot_product(r(n-m:n-1),g(1:m))
		shn=-r(n+m1)+dot_product(r(n+m:n+1:-1),h(1:m))
		sgd=-r(n)+dot_product(r(n-m:n-1),h(m:1:-1))
		if (sd == 0.0 .or. sgd == 0.0) exit
		g(m1)=sgn/sgd
		h(m1)=shn/sd
		t(1:m)=g(1:m)
		g(1:m)=g(1:m)-g(m1)*h(m:1:-1)
		h(1:m)=h(1:m)-h(m1)*t(m:1:-1)
	end do
	if (m > n) call nrerror('toeplz: sanity check failed in routine')
	call nrerror('toeplz: singular principal minor')
	END FUNCTION toeplz

! -----------------------------------------------------------

	SUBROUTINE tptest(data1,data2,t,prob)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq
	use mo_nr, ONLY : avevar,betai
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	REAL(SP), INTENT(OUT) :: t,prob
	INTEGER(I4B) :: n
	REAL(SP) :: ave1,ave2,cov,df,sd,var1,var2
	n=assert_eq(size(data1),size(data2),'tptest')
	call avevar(data1,ave1,var1)
	call avevar(data2,ave2,var2)
	cov=dot_product(data1(:)-ave1,data2(:)-ave2)
	df=n-1
	cov=cov/df
	sd=sqrt((var1+var2-2.0_sp*cov)/n)
	t=(ave1-ave2)/sd
	prob=betai(0.5_sp*df,0.5_sp,df/(df+t**2))
	END SUBROUTINE tptest

! -----------------------------------------------------------

	SUBROUTINE tqli(d,e,z)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : pythag
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: d,e
	REAL(SP), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
	INTEGER(I4B) :: i,iter,l,m,n,ndum
	REAL(SP) :: b,c,dd,f,g,p,r,s
	REAL(SP), DIMENSION(size(e)) :: ff
	n=assert_eq(size(d),size(e),'tqli: n')
	if (present(z)) ndum=assert_eq(n,size(z,1),size(z,2),'tqli: ndum')
	e(:)=eoshift(e(:),1)
	do l=1,n
		iter=0
		iterate: do
			do m=l,n-1
				dd=abs(d(m))+abs(d(m+1))
				if (abs(e(m))+dd == dd) exit
			end do
			if (m == l) exit iterate
			if (iter == 30) call nrerror('too many iterations in tqli')
			iter=iter+1
			g=(d(l+1)-d(l))/(2.0_sp*e(l))
			r=pythag(g,1.0_sp)
			g=d(m)-d(l)+e(l)/(g+sign(r,g))
			s=1.0
			c=1.0
			p=0.0
			do i=m-1,l,-1
				f=s*e(i)
				b=c*e(i)
				r=pythag(f,g)
				e(i+1)=r
				if (r == 0.0) then
					d(i+1)=d(i+1)-p
					e(m)=0.0
					cycle iterate
				end if
				s=f/r
				c=g/r
				g=d(i+1)-p
				r=(d(i)-g)*s+2.0_sp*c*b
				p=s*r
				d(i+1)=g+p
				g=c*r-b
				if (present(z)) then
					ff(1:n)=z(1:n,i+1)
					z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
					z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
				end if
			end do
			d(l)=d(l)-p
			e(l)=g
			e(m)=0.0
		end do iterate
	end do
	END SUBROUTINE tqli

! -----------------------------------------------------------

	SUBROUTINE trapzd(func,a,b,s,n)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), INTENT(INOUT) :: s
	INTEGER(I4B), INTENT(IN) :: n
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP) :: del,fsum
	INTEGER(I4B) :: it
	if (n == 1) then
		s=0.5_sp*(b-a)*sum(func( (/ a,b /) ))
	else
		it=2**(n-2)
		del=(b-a)/it
		fsum=sum(func(arth(a+0.5_sp*del,del,it)))
		s=0.5_sp*(s+del*fsum)
	end if
	END SUBROUTINE trapzd

! -----------------------------------------------------------

	SUBROUTINE tred2(a,d,e,novectors)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,outerprod
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: d,e
	LOGICAL(LGT), OPTIONAL, INTENT(IN) :: novectors
	INTEGER(I4B) :: i,j,l,n
	REAL(SP) :: f,g,h,hh,scale
	REAL(SP), DIMENSION(size(a,1)) :: gg
	LOGICAL(LGT), SAVE :: yesvec=.true.
	n=assert_eq(size(a,1),size(a,2),size(d),size(e),'tred2')
	if (present(novectors)) yesvec=.not. novectors
	do i=n,2,-1
		l=i-1
		h=0.0
		if (l > 1) then
			scale=sum(abs(a(i,1:l)))
			if (scale == 0.0) then
				e(i)=a(i,l)
			else
				a(i,1:l)=a(i,1:l)/scale
				h=sum(a(i,1:l)**2)
				f=a(i,l)
				g=-sign(sqrt(h),f)
				e(i)=scale*g
				h=h-f*g
				a(i,l)=f-g
				if (yesvec) a(1:l,i)=a(i,1:l)/h
				do j=1,l
					e(j)=(dot_product(a(j,1:j),a(i,1:j)) &
					+dot_product(a(j+1:l,j),a(i,j+1:l)))/h
				end do
				f=dot_product(e(1:l),a(i,1:l))
				hh=f/(h+h)
				e(1:l)=e(1:l)-hh*a(i,1:l)
				do j=1,l
					a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
				end do
			end if
		else
			e(i)=a(i,l)
		end if
		d(i)=h
	end do
	if (yesvec) d(1)=0.0
	e(1)=0.0
	do i=1,n
		if (yesvec) then
			l=i-1
			if (d(i) /= 0.0) then
				gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
				a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
			end if
			d(i)=a(i,i)
			a(i,i)=1.0
			a(i,1:l)=0.0
			a(1:l,i)=0.0
		else
			d(i)=a(i,i)
		end if
	end do
	END SUBROUTINE tred2

! -----------------------------------------------------------

	SUBROUTINE tridag_ser(a,b,c,r,u)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
	REAL(SP), DIMENSION(:), INTENT(OUT) :: u
	REAL(SP), DIMENSION(size(b)) :: gam
	INTEGER(I4B) :: n,j
	REAL(SP) :: bet
	n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
	bet=b(1)
	if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
	u(1)=r(1)/bet
	do j=2,n
		gam(j)=c(j-1)/bet
		bet=b(j)-a(j-1)*gam(j)
		if (bet == 0.0) &
			call nrerror('tridag_ser: Error at code stage 2')
		u(j)=(r(j)-a(j-1)*u(j-1))/bet
	end do
	do j=n-1,1,-1
		u(j)=u(j)-gam(j+1)*u(j+1)
	end do
	END SUBROUTINE tridag_ser

	RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : tridag_ser
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
	REAL(SP), DIMENSION(:), INTENT(OUT) :: u
	INTEGER(I4B), PARAMETER :: NPAR_TRIDAG=4
	INTEGER(I4B) :: n,n2,nm,nx
	REAL(SP), DIMENSION(size(b)/2) :: y,q,piva
	REAL(SP), DIMENSION(size(b)/2-1) :: x,z
	REAL(SP), DIMENSION(size(a)/2) :: pivc
	n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_par')
	if (n < NPAR_TRIDAG) then
		call tridag_ser(a,b,c,r,u)
	else
		if (maxval(abs(b(1:n))) == 0.0) &
			call nrerror('tridag_par: possible singular matrix')
		n2=size(y)
		nm=size(pivc)
		nx=size(x)
		piva = a(1:n-1:2)/b(1:n-1:2)
		pivc = c(2:n-1:2)/b(3:n:2)
		y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
		q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
		if (nm < n2) then
			y(n2) = b(n)-piva(n2)*c(n-1)
			q(n2) = r(n)-piva(n2)*r(n-1)
		end if
		x = -piva(2:n2)*a(2:n-2:2)
		z = -pivc(1:nx)*c(3:n-1:2)
		call tridag_par(x,y,z,q,u(2:n:2))
		u(1) = (r(1)-c(1)*u(2))/b(1)
		u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
			-c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
		if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
	end if
	END SUBROUTINE tridag_par

! -----------------------------------------------------------

	SUBROUTINE ttest(data1,data2,t,prob)
	use mo_nrtype
	use mo_nr, ONLY : avevar,betai
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	REAL(SP), INTENT(OUT) :: t,prob
	INTEGER(I4B) :: n1,n2
	REAL(SP) :: ave1,ave2,df,var,var1,var2
	n1=size(data1)
	n2=size(data2)
	call avevar(data1,ave1,var1)
	call avevar(data2,ave2,var2)
	df=n1+n2-2
	var=((n1-1)*var1+(n2-1)*var2)/df
	t=(ave1-ave2)/sqrt(var*(1.0_sp/n1+1.0_sp/n2))
	prob=betai(0.5_sp*df,0.5_sp,df/(df+t**2))
	END SUBROUTINE ttest

! -----------------------------------------------------------

	SUBROUTINE tutest(data1,data2,t,prob)
	use mo_nrtype
	use mo_nr, ONLY : avevar,betai
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	REAL(SP), INTENT(OUT) :: t,prob
	INTEGER(I4B) :: n1,n2
	REAL(SP) :: ave1,ave2,df,var1,var2
	n1=size(data1)
	n2=size(data2)
	call avevar(data1,ave1,var1)
	call avevar(data2,ave2,var2)
	t=(ave1-ave2)/sqrt(var1/n1+var2/n2)
	df=(var1/n1+var2/n2)**2/((var1/n1)**2/(n1-1)+(var2/n2)**2/(n2-1))
	prob=betai(0.5_sp*df,0.5_sp,df/(df+t**2))
	END SUBROUTINE tutest

! -----------------------------------------------------------

	SUBROUTINE twofft(data1,data2,fft1,fft2)
	use mo_nrtype
          use mo_nrutil, ONLY : assert,assert_eq
	use mo_nr, ONLY : four1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: fft1,fft2
	INTEGER(I4B) :: n,n2
	COMPLEX(SPC), PARAMETER :: C1=(0.5_sp,0.0_sp), C2=(0.0_sp,-0.5_sp)
	COMPLEX, DIMENSION(size(data1)/2+1) :: h1,h2
	n=assert_eq(size(data1),size(data2),size(fft1),size(fft2),'twofft')
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in twofft')
	fft1=cmplx(data1,data2,kind=spc)
	call four1(fft1,1)
	fft2(1)=cmplx(aimag(fft1(1)),0.0_sp,kind=spc)
	fft1(1)=cmplx(real(fft1(1)),0.0_sp,kind=spc)
	n2=n/2+1
	h1(2:n2)=C1*(fft1(2:n2)+conjg(fft1(n:n2:-1)))
	h2(2:n2)=C2*(fft1(2:n2)-conjg(fft1(n:n2:-1)))
	fft1(2:n2)=h1(2:n2)
	fft1(n:n2:-1)=conjg(h1(2:n2))
	fft2(2:n2)=h2(2:n2)
	fft2(n:n2:-1)=conjg(h2(2:n2))
	END SUBROUTINE twofft

! -----------------------------------------------------------

	FUNCTION vander(x,q)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,outerdiff
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x,q
	REAL(DP), DIMENSION(size(x)) :: vander
	REAL(DP), DIMENSION(size(x)) :: c
	REAL(DP), DIMENSION(size(x),size(x)) :: a
	INTEGER(I4B) :: i,n
	n=assert_eq(size(x),size(q),'vander')
	if (n == 1) then
		vander(1)=q(1)
	else
		c(:)=0.0
		c(n)=-x(1)
		do i=2,n
			c(n+1-i:n-1)=c(n+1-i:n-1)-x(i)*c(n+2-i:n)
			c(n)=c(n)-x(i)
		end do
		a(:,:)=outerdiff(x,x)
		vander(:)=product(a,dim=2,mask=(a /= 0.0))
		a(:,1)=-c(1)/x(:)
		do i=2,n
			a(:,i)=-(c(i)-a(:,i-1))/x(:)
		end do
		vander(:)=matmul(a,q)/vander(:)
	end if
	END FUNCTION vander

! -----------------------------------------------------------

	SUBROUTINE vegas(region,func,init,ncall,itmx,nprn,tgral,sd,chi2a)
	use mo_nrtype
	use mo_nr, ONLY : ran1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: region
	INTEGER(I4B), INTENT(IN) :: init,ncall,itmx,nprn
	REAL(SP), INTENT(OUT) :: tgral,sd,chi2a
	INTERFACE
		FUNCTION func(pt,wgt)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: pt
		REAL(SP), INTENT(IN) :: wgt
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), PARAMETER :: ALPH=1.5_sp,TINY=1.0e-30_sp
	INTEGER(I4B), PARAMETER :: MXDIM=10,NDMX=50
	INTEGER(I4B), SAVE :: i,it,j,k,mds,nd,ndim,ndo,ng,npg
	INTEGER(I4B), DIMENSION(MXDIM), SAVE :: ia,kg
	REAL(SP), SAVE :: calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo,harvest
	REAL(SP), DIMENSION(NDMX,MXDIM), SAVE :: d,di,xi
	REAL(SP), DIMENSION(MXDIM), SAVE :: dt,dx,x
	REAL(SP), DIMENSION(NDMX), SAVE :: r,xin
	REAL(DP), SAVE :: schi,si,swgt
	ndim=size(region)/2
	if (init <= 0) then
		mds=1
		ndo=1
		xi(1,:)=1.0
	end if
	if (init <= 1) then
		si=0.0
		swgt=0.0
		schi=0.0
	end if
	if (init <= 2) then
		nd=NDMX
		ng=1
		if (mds /= 0) then
			ng=(ncall/2.0_sp+0.25_sp)**(1.0_sp/ndim)
			mds=1
			if ((2*ng-NDMX) >= 0) then
				mds=-1
				npg=ng/NDMX+1
				nd=ng/npg
				ng=npg*nd
			end if
		end if
		k=ng**ndim
		npg=max(ncall/k,2)
		calls=real(npg,sp)*real(k,sp)
		dxg=1.0_sp/ng
		dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-1.0_sp)
		xnd=nd
		dxg=dxg*xnd
		dx(1:ndim)=region(1+ndim:2*ndim)-region(1:ndim)
		xjac=1.0_sp/calls*product(dx(1:ndim))
		if (nd /= ndo) then
			r(1:max(nd,ndo))=1.0
			do j=1,ndim
				call rebin(ndo/xnd,nd,r,xin,xi(:,j))
			end do
			ndo=nd
		end if
		if (nprn >= 0) write(*,200) ndim,calls,it,itmx,nprn,&
			ALPH,mds,nd,(j,region(j),j,region(j+ndim),j=1,ndim)
	end if
	do it=1,itmx
		ti=0.0
		tsi=0.0
		kg(:)=1
		d(1:nd,:)=0.0
		di(1:nd,:)=0.0
		iterate: do
			fb=0.0
			f2b=0.0
			do k=1,npg
				wgt=xjac
				do j=1,ndim
					call ran1(harvest)
					xn=(kg(j)-harvest)*dxg+1.0_sp
					ia(j)=max(min(int(xn),NDMX),1)
					if (ia(j) > 1) then
						xo=xi(ia(j),j)-xi(ia(j)-1,j)
						rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
					else
						xo=xi(ia(j),j)
						rc=(xn-ia(j))*xo
					end if
					x(j)=region(j)+rc*dx(j)
					wgt=wgt*xo*xnd
				end do
				f=wgt*func(x(1:ndim),wgt)
				f2=f*f
				fb=fb+f
				f2b=f2b+f2
				do j=1,ndim
					di(ia(j),j)=di(ia(j),j)+f
					if (mds >= 0) d(ia(j),j)=d(ia(j),j)+f2
				end do
			end do
			f2b=sqrt(f2b*npg)
			f2b=(f2b-fb)*(f2b+fb)
			if (f2b <= 0.0) f2b=TINY
			ti=ti+fb
			tsi=tsi+f2b
			if (mds < 0) then
				do j=1,ndim
					d(ia(j),j)=d(ia(j),j)+f2b
				end do
			end if
			do k=ndim,1,-1
				kg(k)=mod(kg(k),ng)+1
				if (kg(k) /= 1) cycle iterate
			end do
			exit iterate
		end do iterate
		tsi=tsi*dv2g
		wgt=1.0_sp/tsi
		si=si+real(wgt,dp)*real(ti,dp)
		schi=schi+real(wgt,dp)*real(ti,dp)**2
		swgt=swgt+real(wgt,dp)
		tgral=si/swgt
		chi2a=max((schi-si*tgral)/(it-0.99_dp),0.0_dp)
		sd=sqrt(1.0_sp/swgt)
		tsi=sqrt(tsi)
		if (nprn >= 0) then
			write(*,201) it,ti,tsi,tgral,sd,chi2a
			if (nprn /= 0) then
				do j=1,ndim
					write(*,202) j,(xi(i,j),di(i,j),&
						i=1+nprn/2,nd,nprn)
				end do
			end if
		end if
		do j=1,ndim
			xo=d(1,j)
			xn=d(2,j)
			d(1,j)=(xo+xn)/2.0_sp
			dt(j)=d(1,j)
			do i=2,nd-1
				rc=xo+xn
				xo=xn
				xn=d(i+1,j)
				d(i,j)=(rc+xn)/3.0_sp
				dt(j)=dt(j)+d(i,j)
			end do
			d(nd,j)=(xo+xn)/2.0_sp
			dt(j)=dt(j)+d(nd,j)
		end do
		where (d(1:nd,:) < TINY) d(1:nd,:)=TINY
		do j=1,ndim
			r(1:nd)=((1.0_sp-d(1:nd,j)/dt(j))/(log(dt(j))-log(d(1:nd,j))))**ALPH
			rc=sum(r(1:nd))
			call rebin(rc/xnd,nd,r,xin,xi(:,j))
		end do
	end do
200	format(/' input parameters for vegas:  ndim=',i3,'  ncall=',f8.0&
		/28x,'  it=',i5,'  itmx=',i5&
		/28x,'  nprn=',i3,'  alph=',f5.2/28x,'  mds=',i3,'   nd=',i4&
		/(30x,'xl(',i2,')= ',g11.4,' xu(',i2,')= ',g11.4))
201	format(/' iteration no.',I3,': ','integral =',g14.7,' +/- ',g9.2,&
		/' all iterations:   integral =',g14.7,' +/- ',g9.2,&
		' chi**2/it''n =',g9.2)
202	format(/' data for axis ',I2/'    X       delta i       ',&
		'   x       delta i       ','    x       delta i       ',&
		/(1x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4))
	CONTAINS
!BL
	SUBROUTINE rebin(rc,nd,r,xin,xi)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: rc
	INTEGER(I4B), INTENT(IN) :: nd
	REAL(SP), DIMENSION(:), INTENT(IN) :: r
	REAL(SP), DIMENSION(:), INTENT(OUT) :: xin
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: xi
	INTEGER(I4B) :: i,k
	REAL(SP) :: dr,xn,xo
	k=0
	xo=0.0
	dr=0.0
	do i=1,nd-1
		do
			if (rc <= dr) exit
			k=k+1
			dr=dr+r(k)
		end do
		if (k > 1) xo=xi(k-1)
		xn=xi(k)
		dr=dr-rc
		xin(i)=xn-(xn-xo)*dr/r(k)
	end do
	xi(1:nd-1)=xin(1:nd-1)
	xi(nd)=1.0
	END SUBROUTINE rebin
	END SUBROUTINE vegas

! -----------------------------------------------------------

	SUBROUTINE voltra(t0,h,t,f,g,ak)
	use mo_nrtype
          use mo_nrutil, ONLY : array_copy,assert_eq,unit_matrix
	use mo_nr, ONLY : lubksb,ludcmp
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: t0,h
	REAL(SP), DIMENSION(:), INTENT(OUT) :: t
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: f
	INTERFACE
		FUNCTION g(t)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: t
		REAL(SP), DIMENSION(:), POINTER :: g
		END FUNCTION g
!BL
		FUNCTION ak(t,s)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: t,s
		REAL(SP), DIMENSION(:,:), POINTER :: ak
		END FUNCTION ak
	END INTERFACE
	INTEGER(I4B) :: i,j,n,ncop,nerr,m
	INTEGER(I4B), DIMENSION(size(f,1)) :: indx
	REAL(SP) :: d
	REAL(SP), DIMENSION(size(f,1)) :: b
	REAL(SP), DIMENSION(size(f,1),size(f,1)) :: a
	n=assert_eq(size(f,2),size(t),'voltra: n')
	t(1)=t0
	call array_copy(g(t(1)),f(:,1),ncop,nerr)
	m=assert_eq(size(f,1),ncop,ncop+nerr,'voltra: m')
	do i=2,n
		t(i)=t(i-1)+h
		b=g(t(i))+0.5_sp*h*matmul(ak(t(i),t(1)),f(:,1))
		do j=2,i-1
			b=b+h*matmul(ak(t(i),t(j)),f(:,j))
		end do
		call unit_matrix(a)
		a=a-0.5_sp*h*ak(t(i),t(i))
		call ludcmp(a,indx,d)
		call lubksb(a,indx,b)
		f(:,i)=b(:)
	end do
	END SUBROUTINE voltra

! -----------------------------------------------------------

	SUBROUTINE wt1(a,isign,wtstep)
	use mo_nrtype
          use mo_nrutil, ONLY : assert
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
	INTEGER(I4B), INTENT(IN) :: isign
	INTERFACE
		SUBROUTINE wtstep(a,isign)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE wtstep
	END INTERFACE
	INTEGER(I4B) :: n,nn
	n=size(a)
	call assert(iand(n,n-1)==0, 'n must be a power of 2 in wt1')
	if (n < 4) RETURN
	if (isign >= 0) then
		nn=n
		do
			if (nn < 4) exit
			call wtstep(a(1:nn),isign)
			nn=nn/2
		end do
	else
		nn=4
		do
			if (nn > n) exit
			call wtstep(a(1:nn),isign)
			nn=nn*2
		end do
	end if
	END SUBROUTINE wt1

! -----------------------------------------------------------

	SUBROUTINE wtn(a,nn,isign,wtstep)
	use mo_nrtype
          use mo_nrutil, ONLY : arth,assert
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nn
	INTEGER(I4B), INTENT(IN) :: isign
	INTERFACE
		SUBROUTINE wtstep(a,isign)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE wtstep
	END INTERFACE
	INTEGER(I4B) :: i1,i2,i3,idim,n,ndim,nnew,nprev,nt,ntot
	REAL(SP), DIMENSION(:), ALLOCATABLE :: wksp
	call assert(iand(nn,nn-1)==0, 'each dimension must be a power of 2 in wtn')
	allocate(wksp(maxval(nn)))
	ndim=size(nn)
	ntot=product(nn(:))
	nprev=1
	do idim=1,ndim
		n=nn(idim)
		nnew=n*nprev
		if (n > 4) then
			do i2=0,ntot-1,nnew
				do i1=1,nprev
					i3=i1+i2
					wksp(1:n)=a(arth(i3,nprev,n))
					i3=i3+n*nprev
					if (isign >= 0) then
						nt=n
						do
							if (nt < 4) exit
							call wtstep(wksp(1:nt),isign)
							nt=nt/2
						end do
					else
						nt=4
						do
							if (nt > n) exit
							call wtstep(wksp(1:nt),isign)
							nt=nt*2
						end do
					end if
					i3=i1+i2
					a(arth(i3,nprev,n))=wksp(1:n)
					i3=i3+n*nprev
				end do
			end do
		end if
		nprev=nnew
	end do
	deallocate(wksp)
	END SUBROUTINE wtn

! -----------------------------------------------------------

	FUNCTION wwghts(n,h,kermom)
	use mo_nrtype
          use mo_nrutil, ONLY : geop
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: h
	REAL(SP), DIMENSION(n) :: wwghts
	INTERFACE
		FUNCTION kermom(y,m)
		use mo_nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: y
		INTEGER(I4B), INTENT(IN) :: m
		REAL(DP), DIMENSION(m) :: kermom
		END FUNCTION kermom
	END INTERFACE
	INTEGER(I4B) :: j
	REAL(DP) :: hh,hi,c,a,b
	REAL(DP), DIMENSION(4) :: wold,wnew,w
	hh=h
	hi=1.0_dp/hh
	wwghts(1:n)=0.0
	wold(1:4)=kermom(0.0_dp,4)
	if (n >= 4) then
		b=0.0
		do j=1,n-3
			c=j-1
			a=b
			b=a+hh
			if (j == n-3) b=(n-1)*hh
			wnew(1:4)=kermom(b,4)
			w(1:4)=(wnew(1:4)-wold(1:4))*geop(1.0_dp,hi,4)
			wwghts(j:j+3)=wwghts(j:j+3)+(/&
				((c+1.0_dp)*(c+2.0_dp)*(c+3.0_dp)*w(1)&
				-(11.0_dp+c*(12.0_dp+c*3.0_dp))*w(2)&
					+3.0_dp*(c+2.0_dp)*w(3)-w(4))/6.0_dp,&
				(-c*(c+2.0_dp)*(c+3.0_dp)*w(1)&
				+(6.0_dp+c*(10.0_dp+c*3.0_dp))*w(2)&
					-(3.0_dp*c+5.0_dp)*w(3)+w(4))*0.50_dp,&
				(c*(c+1.0_dp)*(c+3.0_dp)*w(1)&
				-(3.0_dp+c*(8.0_dp+c*3.0_dp))*w(2)&
					+(3.0_dp*c+4.0_dp)*w(3)-w(4))*0.50_dp,&
				(-c*(c+1.0_dp)*(c+2.0_dp)*w(1)&
				+(2.0_dp+c*(6.0_dp+c*3.0_dp))*w(2)&
				-3.0_dp*(c+1.0_dp)*w(3)+w(4))/6.0_dp /)
			wold(1:4)=wnew(1:4)
		end do
	else if (n == 3) then
		wnew(1:3)=kermom(hh+hh,3)
		w(1:3)= (/ wnew(1)-wold(1), hi*(wnew(2)-wold(2)),&
			hi**2*(wnew(3)-wold(3)) /)
		wwghts(1:3)= (/ w(1)-1.50_dp*w(2)+0.50_dp*w(3),&
			2.0_dp*w(2)-w(3), 0.50_dp*(w(3)-w(2)) /)
	else if (n == 2) then
		wnew(1:2)=kermom(hh,2)
		wwghts(2)=hi*(wnew(2)-wold(2))
		wwghts(1)=wnew(1)-wold(1)-wwghts(2)
	end if
	END FUNCTION wwghts

! -----------------------------------------------------------

	SUBROUTINE zbrac(func,x1,x2,succes)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(INOUT) :: x1,x2
	LOGICAL(LGT), INTENT(OUT) :: succes
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: NTRY=50
	REAL(SP), PARAMETER :: FACTOR=1.6_sp
	INTEGER(I4B) :: j
	REAL(SP) :: f1,f2
	if (x1 == x2) call nrerror('zbrac: you have to guess an initial range')
	f1=func(x1)
	f2=func(x2)
	succes=.true.
	do j=1,NTRY
		if ((f1 > 0.0 .and. f2 < 0.0) .or. &
			(f1 < 0.0 .and. f2 > 0.0)) RETURN
		if (abs(f1) < abs(f2)) then
			x1=x1+FACTOR*(x1-x2)
			f1=func(x1)
		else
			x2=x2+FACTOR*(x2-x1)
			f2=func(x2)
		end if
	end do
	succes=.false.
	END SUBROUTINE zbrac

! -----------------------------------------------------------

	SUBROUTINE zbrak(func,x1,x2,n,xb1,xb2,nb)
	use mo_nrtype
          use mo_nrutil, ONLY : arth
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B), INTENT(OUT) :: nb
	REAL(SP), INTENT(IN) :: x1,x2
	REAL(SP), DIMENSION(:), POINTER :: xb1,xb2
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B) :: i
	REAL(SP) :: dx
	REAL(SP), DIMENSION(0:n) :: f,x
	LOGICAL(LGT), DIMENSION(1:n) :: mask
	LOGICAL(LGT), SAVE :: init=.true.
	if (init) then
		init=.false.
		nullify(xb1,xb2)
	end if
	if (associated(xb1)) deallocate(xb1)
	if (associated(xb2)) deallocate(xb2)
	dx=(x2-x1)/n
	x=x1+dx*arth(0,1,n+1)
	do i=0,n
		f(i)=func(x(i))
	end do
	mask=f(1:n)*f(0:n-1) <= 0.0
	nb=count(mask)
	allocate(xb1(nb),xb2(nb))
	xb1(1:nb)=pack(x(0:n-1),mask)
	xb2(1:nb)=pack(x(1:n),mask)
	END SUBROUTINE zbrak

! -----------------------------------------------------------

	FUNCTION zbrent(func,x1,x2,tol)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,tol
	REAL(SP) :: zbrent
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: EPS=epsilon(x1)
	INTEGER(I4B) :: iter
	REAL(SP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
	a=x1
	b=x2
	fa=func(a)
	fb=func(b)
	if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
		call nrerror('root must be bracketed for zbrent')
	c=b
	fc=fb
	do iter=1,ITMAX
		if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
			c=a
			fc=fa
			d=b-a
			e=d
		end if
		if (abs(fc) < abs(fb)) then
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc=fa
		end if
		tol1=2.0_sp*EPS*abs(b)+0.5_sp*tol
		xm=0.5_sp*(c-b)
		if (abs(xm) <= tol1 .or. fb == 0.0) then
			zbrent=b
			RETURN
		end if
		if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
			s=fb/fa
			if (a == c) then
				p=2.0_sp*xm*s
				q=1.0_sp-s
			else
				q=fa/fc
				r=fb/fc
				p=s*(2.0_sp*xm*q*(q-r)-(b-a)*(r-1.0_sp))
				q=(q-1.0_sp)*(r-1.0_sp)*(s-1.0_sp)
			end if
			if (p > 0.0) q=-q
			p=abs(p)
			if (2.0_sp*p  <  min(3.0_sp*xm*q-abs(tol1*q),abs(e*q))) then
				e=d
				d=p/q
			else
				d=xm
				e=d
			end if
		else
			d=xm
			e=d
		end if
		a=b
		fa=fb
		b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
		fb=func(b)
	end do
	call nrerror('zbrent: exceeded maximum iterations')
	zbrent=b
	END FUNCTION zbrent

! -----------------------------------------------------------

	SUBROUTINE zrhqr(a,rtr,rti)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,nrerror
	use mo_nr, ONLY : balanc,hqr,indexx
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: rtr,rti
	INTEGER(I4B) :: k,m
	INTEGER(I4B), DIMENSION(size(rtr)) :: indx
	REAL(SP), DIMENSION(size(a)-1,size(a)-1) :: hess
	m=assert_eq(size(rtr),size(rti),size(a)-1,'zrhqr')
	if (a(m+1) == 0.0) call &
		nrerror('zrhqr: Last value of array a must not be 0')
	hess(1,:)=-a(m:1:-1)/a(m+1)
	hess(2:m,:)=0.0
	do k=1,m-1
		hess(k+1,k)=1.0
	end do
	call balanc(hess)
	call hqr(hess,rtr,rti)
	call indexx(rtr,indx)
	rtr=rtr(indx)
	rti=rti(indx)
	END SUBROUTINE zrhqr

! -----------------------------------------------------------

	FUNCTION zriddr(func,x1,x2,xacc)
	use mo_nrtype
          use mo_nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,xacc
	REAL(SP) :: zriddr
	INTERFACE
		FUNCTION func(x)
		use mo_nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXIT=60
	REAL(SP), PARAMETER :: UNUSED=-1.11e30_sp
	INTEGER(I4B) :: j
	REAL(SP) :: fh,fl,fm,fnew,s,xh,xl,xm,xnew
	fl=func(x1)
	fh=func(x2)
	if ((fl > 0.0 .and. fh < 0.0) .or. (fl < 0.0 .and. fh > 0.0)) then
		xl=x1
		xh=x2
		zriddr=UNUSED
		do j=1,MAXIT
			xm=0.5_sp*(xl+xh)
			fm=func(xm)
			s=sqrt(fm**2-fl*fh)
			if (s == 0.0) RETURN
			xnew=xm+(xm-xl)*(sign(1.0_sp,fl-fh)*fm/s)
			if (abs(xnew-zriddr) <= xacc) RETURN
			zriddr=xnew
			fnew=func(zriddr)
			if (fnew == 0.0) RETURN
			if (sign(fm,fnew) /= fm) then
				xl=xm
				fl=fm
				xh=zriddr
				fh=fnew
			else if (sign(fl,fnew) /= fl) then
				xh=zriddr
				fh=fnew
			else if (sign(fh,fnew) /= fh) then
				xl=zriddr
				fl=fnew
			else
				call nrerror('zriddr: never get here')
			end if
			if (abs(xh-xl) <= xacc) RETURN
		end do
		call nrerror('zriddr: exceeded maximum iterations')
	else if (fl == 0.0) then
		zriddr=x1
	else if (fh == 0.0) then
		zriddr=x2
	else
		call nrerror('zriddr: root must be bracketed')
	end if
	END FUNCTION zriddr

! -----------------------------------------------------------

	SUBROUTINE zroots(a,roots,polish)
	use mo_nrtype
          use mo_nrutil, ONLY : assert_eq,poly_term
	use mo_nr, ONLY : laguer,indexx
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
	COMPLEX(SPC), DIMENSION(:), INTENT(OUT) :: roots
	LOGICAL(LGT), INTENT(IN) :: polish
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	INTEGER(I4B) :: j,its,m
	INTEGER(I4B), DIMENSION(size(roots)) :: indx
	COMPLEX(SPC) :: x
	COMPLEX(SPC), DIMENSION(size(a)) :: ad
	m=assert_eq(size(roots),size(a)-1,'zroots')
	ad(:)=a(:)
	do j=m,1,-1
		x=cmplx(0.0_sp,kind=spc)
		call laguer(ad(1:j+1),x,its)
		if (abs(aimag(x)) <= 2.0_sp*EPS**2*abs(real(x))) &
			x=cmplx(real(x),kind=spc)
		roots(j)=x
		ad(j:1:-1)=poly_term(ad(j+1:2:-1),x)
	end do
	if (polish) then
		do j=1,m
			call laguer(a(:),roots(j),its)
		end do
	end if
	call indexx(real(roots),indx)
	roots=roots(indx)
	END SUBROUTINE zroots
END MODULE mo_nr

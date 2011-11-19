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

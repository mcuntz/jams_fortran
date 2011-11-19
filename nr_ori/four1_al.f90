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

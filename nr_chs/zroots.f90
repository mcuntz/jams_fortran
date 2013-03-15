subroutine zroots_sp(a,roots,polish)
  use mo_kind,   only : i4, spc, sp 
  use mo_nrutil, only : assert_eq, poly_term
  use mo_nr,     only : laguer
  use mo_sort,   only : sort_index 

  implicit none

  complex(spc), dimension(:), intent(in)  :: a
  complex(spc), dimension(:), intent(out) :: roots
  logical,                    intent(in)  :: polish
  real(sp), parameter                     :: eps=1.0e-6_sp
  integer(i4)                             :: j,its,m
  integer(i4),  dimension(size(roots))    :: indx
  complex(spc)                            :: x
  complex(spc), dimension(size(a))        :: ad

  m=assert_eq(size(roots),size(a)-1,'zroots')
  ad(:)=a(:)
  do j=m,1,-1
     x=cmplx(0.0_sp,kind=spc)
     call laguer(ad(1:j+1),x,its)
     if (abs(aimag(x)) <= 2.0_sp*eps**2*abs(real(x))) &
          x=cmplx(real(x),kind=spc)
     roots(j)=x
     ad(j:1:-1)=poly_term(ad(j+1:2:-1),x)
  end do
  if (polish) then
     do j=1,m
        call laguer(a(:),roots(j),its)
     end do
  end if
  indx = sort_index(real(roots))
  !call indexx(real(roots),indx)
  roots=roots(indx)
end subroutine zroots_sp

subroutine zroots_dp(a,roots,polish)
  use mo_kind,   only : i4, dpc, dp 
  use mo_nrutil, only : assert_eq, poly_term
  use mo_nr,     only : laguer
  use mo_sort,   only : sort_index 

  implicit none

  complex(dpc), dimension(:), intent(in)  :: a
  complex(dpc), dimension(:), intent(out) :: roots
  logical,                    intent(in)  :: polish
  real(dp), parameter                     :: eps=1.0e-6_dp
  integer(i4)                             :: j,its,m
  integer(i4),  dimension(size(roots))    :: indx
  complex(dpc)                            :: x
  complex(dpc), dimension(size(a))        :: ad

  m=assert_eq(size(roots),size(a)-1,'zroots')
  ad(:)=a(:)
  do j=m,1,-1
     x=cmplx(0.0_dp,kind=dpc)
     call laguer(ad(1:j+1),x,its)
     if (abs(aimag(x)) <= 2.0_dp*eps**2*abs(real(x))) &
          x=cmplx(real(x),kind=dpc)
     roots(j)=x
     ad(j:1:-1)=poly_term(ad(j+1:2:-1),x)
  end do
  if (polish) then
     do j=1,m
        call laguer(a(:),roots(j),its)
     end do
  end if
  indx = sort_index(real(roots))
  ! call indexx(real(roots),indx)
  roots=roots(indx)
end subroutine zroots_dp

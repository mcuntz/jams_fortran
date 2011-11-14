!******************************************************************************
!  MODULE
!          mo_kind
!  DESCRIPTION
!          Declares number precision
!          File in parts copied from Numerical Recipes
!******************************************************************************
MODULE mo_kind
    	INTEGER, PARAMETER :: I8 = SELECTED_INT_KIND(18)       ! 8 Byte Integer
	INTEGER, PARAMETER :: I4 = SELECTED_INT_KIND(9)        ! 4 Byte Integer
	INTEGER, PARAMETER :: I2 = SELECTED_INT_KIND(4)        ! 2 Byte Integer
	INTEGER, PARAMETER :: I1 = SELECTED_INT_KIND(2)        ! 1 Byte Integer
	INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,37)    ! Single Precision Real
	INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15,307)  ! Double Precision Real
	INTEGER, PARAMETER :: SPC = KIND((1.0_SP,1.0_SP))      ! Single Precision Complex
	INTEGER, PARAMETER :: DPC = KIND((1.0_DP,1.0_DP))      ! Double Precision Complex
	INTEGER, PARAMETER :: LGT = KIND(.true.)               ! Logical

        ! Type for Single Precision Sparse Arrays
	TYPE sprs2_sp
		INTEGER(I4) :: n,len
		REAL(SP), DIMENSION(:), POINTER :: val
		INTEGER(I4), DIMENSION(:), POINTER :: irow
		INTEGER(I4), DIMENSION(:), POINTER :: jcol
	END TYPE sprs2_sp
        ! Type for Double Precision Sparse Arrays
	TYPE sprs2_dp
		INTEGER(I4) :: n,len
		REAL(DP), DIMENSION(:), POINTER :: val
		INTEGER(I4), DIMENSION(:), POINTER :: irow
		INTEGER(I4), DIMENSION(:), POINTER :: jcol
	END TYPE sprs2_dp
END MODULE mo_kind

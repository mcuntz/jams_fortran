MODULE mo_pumpingtests

  ! This module serves different solutions for the groundwater-flow equation.

  ! Four kinds of solution are implemented:
  !     -Thiem's solution for steady-state and homogeneous aquifer                      :   thiem
  !     -Theis' solution for transient flow and homogenous aquifer                      :   theis
  !     -the extended thiem's solution for steady state and heterogeneous aquifer in 2D :   ext_thiem2d
  !     -the extended thiem's solution for steady state and heterogeneous aquifer and
  !      anisotropical impact in 3D                                                     :   ext_thiem3d
  !     -the extended theis' solution for transient flow and heterogeneous aquifer and
  !      anisotropical impact in 3D                                                     :   ext_theis3d
  !

  ! Written Sebastian Mueller, Apr 2014

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! If you use this routine in your work, you should cite the following reference
  ! Zech A, C L Schneider, and S Attinger (2012)
  !     The Extended Thiem's solution: Including the impact of heterogeneity,
  !     Water Resources Research 48, W10535, doi:10.1029/2012WR011852

  ! The UFZ Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The UFZ Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2014 Sebastian Mueller

  USE mo_kind,      ONLY: i4, sp, dp
  USE mo_constants, ONLY: PI_D, PI
  USE mo_utils,     ONLY: eq, le, ne

  IMPLICIT NONE

  PUBLIC :: thiem       ! thiem's solution for the steady-state and homogeneous aquifer
  PUBLIC :: theis       ! theis' solution for transient-flow and homogeneous aquifer
  PUBLIC :: ext_thiem2d ! the extended thiem's solution for the steady-state and heterogeneous aquifer in 2D
  PUBLIC :: ext_thiem3d ! the extended thiem's solution for the steady-state and heterogeneous aquifer in 3D
  PUBLIC :: ext_theis2d ! the extended theis' solution for the transient-flow and heterogeneous aquifer in 2D
  PUBLIC :: ext_theis3d ! the extended theis' solution for the transient-flow and heterogeneous aquifer in 3D


  ! ------------------------------------------------------------------

  !     NAME
  !         thiem

  !     PURPOSE
  !         Calculates the thiem solution h(r) for the steady-state and homogenous aquifers.

  !     CALLING SEQUENCE
  !         out = thiem(rad, params, inits)

  !     INTENT(IN)
  !         REAL(sp/dp),[DIMENSION(:)]  :: rad                          : input radius                  in m
  !         REAL(sp/dp), DIMENSION(2)   :: params       Parameters
  !                                                         params(1)   : hydraulic conductivity        in m/s
  !                                                         params(2)   : aquifer thickness             in m
  !         REAL(sp/dp), DIMENSION(3)   :: inits        Initial values
  !                                                         inits(1)    : reference point radius        in m
  !                                                         inits(2)    : reference head                in m
  !                                                         inits(3)    : pumping rate at the well      in m^3/s

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURNS
  !>       \return     REAL(sp/dp)      :: thiem        hydraulic head in m

  !     RESTRICTIONS
  !         1) Input values must be floating points.
  !         2) rad       >  0       [radius]
  !         3) inits(1)  >  0       [reference radius]
  !         4) params(1) >  0       [conductivity]
  !         5) params(2) >  0       [aqifer thickness]

  !     EXAMPLE
  !         rad     = 0.5
  !         inits   = (/ 32.0 , 0.0 , -0.001 /)
  !         params  = (/ 0.0001 , 1.0 /)

  !         h       = thiem(rad, inits, params )

  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written, Sebastian Mueller, Apr 2014

  INTERFACE thiem
     MODULE PROCEDURE thiem_sp, thiem_dp, thiem_d1_sp, thiem_d1_dp
  END INTERFACE thiem

  ! ------------------------------------------------------------------

  !     NAME
  !         theis

  !     PURPOSE
  !         Calculates the theis solution h(r) for transient flow and homogenous aquifers.

  !     CALLING SEQUENCE
  !         out = theis(rad, time, params, inits)

  !     INTENT(IN)
  !         REAL(sp/dp),[DIMENSION(:)]  :: rad          Radius          : input radius                  in m
  !         REAL(sp/dp),[DIMENSION(:)]  :: time         Time            : input time                    in s
  !         REAL(sp/dp), DIMENSION(3)   :: params       Parameters
  !                                                         params(1)   : specific storage coefficient  in m^3
  !                                                         params(2)   : hydraulic conductivity        in m/s
  !                                                         params(3)   : aquifer thickness             in m
  !         REAL(sp/dp), DIMENSION(1)   :: inits        Initial values
  !                                                         inits(1)    : pumping rate at the well      in m^3/s

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURNS
  !>       \return     REAL(sp/dp)      :: theis        hydraulic head in m

  !     RESTRICTIONS
  !         1) Input values must be floating points.
  !         2) rad       >  0       [radius]
  !         3) time      >  0       [time]
  !         4) params(1) >  0       [storage coefficient]
  !         5) params(2) >  0       [conductivity]
  !         6) params(3) >  0       [aqifer thickness]

  !     EXAMPLE
  !         rad     = 0.5
  !         time    = 0.5
  !         inits   = (/ -0.001 /)
  !         params  = (/ 0.001 , 0.0001 , 1.0 /)

  !         h       = theis(rad, time, inits, params )

  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written, Sebastian Mueller, Apr 2014

  INTERFACE theis
     MODULE PROCEDURE theis_sp, theis_dp, theis_d1_sp, theis_d1_dp
  END INTERFACE theis

  ! ------------------------------------------------------------------

  !     NAME
  !         ext_thiem2d

  !     PURPOSE
  !         Calculates the extended Thiem solution h(r) steady-state and heterogenous aquifers
  !         in 2 dimensions.

  !     CALLING SEQUENCE
  !         out = ext_thiem2d(rad, params, inits, prop)

  !     INTENT(IN)
  !         REAL(sp/dp),[DIMENSION(:)]  :: rad          Radius          : input radius                  in m
  !         REAL(sp/dp), DIMENSION(3)   :: params       Parameters
  !                                                         params(1)   : geometric mean transmissivity in m^2/s
  !                                                         params(2)   : variance of log-conductivity
  !                                                         params(3)   : correlation length            in m
  !         REAL(sp/dp), DIMENSION(3)   :: inits        Initial values
  !                                                         inits(1)    : reference point radius        in m
  !                                                         inits(2)    : reference head                in m
  !                                                         inits(3)    : pumping rate at the well      in m^3/s

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         REAL(sp/dp)                    :: prop         proportionality factor
  !                                                     default=2.0

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURNS
  !>       \return     REAL(sp/dp)      :: ext_thiem2d      hydraulic head in m

  !     RESTRICTIONS
  !         1) Input values must be floating points.
  !         2) rad       >  0       [radius]
  !         3) inits(1)  >  0       [reference radius]
  !         4) params(1) >  0       [transmissivity]
  !         5) params(2) >  0       [variance]
  !         6) params(3) >  0       [correlation length]

  !     EXAMPLE
  !         rad     = 0.5
  !         inits   = (/ 32.0 , 0.0 , -0.001 /)
  !         params  = (/ 0.0001 , 1.0 , 10.0 /)

  !         h       = ext_thiem2d(rad, inits, params )

  !         -> see also example in test directory

  !     LITERATURE
  !         [1] Zech et al. 2012 The extended Thiem's solution -Including the Impact
  !             of heterogeneity, WRR

  !     HISTORY
  !         Written, Sebastian Mueller, Apr 2014

  INTERFACE ext_thiem2d
     MODULE PROCEDURE ext_thiem2d_sp, ext_thiem2d_dp, ext_thiem2d_d1_sp, ext_thiem2d_d1_dp
  END INTERFACE ext_thiem2d

  ! ------------------------------------------------------------------

  !     NAME
  !         ext_thiem3d

  !     PURPOSE
  !         Calculates the extended Thiem solution h(r) steady-state and heterogenous aquifers
  !         in 3 dimensions with anisotropical impact.

  !     CALLING SEQUENCE
  !         out = ext_thiem3d(rad, params, inits, BC, prop)

  !     INTENT(IN)
  !         REAL(sp/dp),[DIMENSION(:)]  :: rad          Radius          : input radius                  in m
  !         REAL(sp/dp), DIMENSION(5)   :: params       Parameters
  !                                                         params(1)   : geometric mean conductivity   in m^2/s
  !                                                         params(2)   : variance of log-conductivity
  !                                                         params(3)   : correlation length            in m
  !                                                         params(4)   : anisotropy ratio
  !                                                         params(5)   : aquifer thickness             in m
  !         REAL(sp/dp), DIMENSION(3)   :: inits        Initial values
  !                                                         inits(1)    : reference point radius        in m
  !                                                         inits(2)    : reference head                in m
  !                                                         inits(3)    : pumping rate at the well      in m^3/s

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         REAL(sp/dp)                 :: prop         prop-factor for harmonic boundary-condition
  !                                                     the arithmetic case will be automatically corrected
  !                                                     default=1.6
  !         LOGICAL                     :: BC           kind of the boundary-condition
  !                                                     (harmonic: true[default]; arithmetic: false)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURNS
  !>       \return     REAL(sp/dp)      :: ext_thiem3d      hydraulic head in m

  !     RESTRICTIONS
  !         1) Input values must be floating points.
  !         2) rad       >  0       [radius]
  !         3) inits(1)  >  0       [reference radius]
  !         4) params(1) >  0       [conductivity]
  !         5) params(2) >  0       [variance]
  !         6) params(3) >  0       [correlation length]
  !         7) params(4) in [0,1]   [anisotropy ratio]
  !         8) params(5) >  0       [aquifer thickness]

  !     EXAMPLE
  !         rad     = 0.5
  !         inits   = (/ 32.0 , 0.0 , -0.001 /)
  !         params  = (/ 0.0001 , 1.0 , 10.0 , 1.0 , 1.0 /)

  !         h       = ext_thiem3d(rad, inits, params )

  !         -> see also example in test directory

  !     LITERATURE
  !         [1] Zech et al. 2012 The extended Thiem's solution -Including the Impact
  !             of heterogeneity, WRR

  !     HISTORY
  !         Written, Sebastian Mueller, Apr 2014

  INTERFACE ext_thiem3d
     MODULE PROCEDURE ext_thiem3d_sp, ext_thiem3d_dp, ext_thiem3d_d1_sp, ext_thiem3d_d1_dp
  END INTERFACE ext_thiem3d

  ! ------------------------------------------------------------------

  !     NAME
  !         ext_theis3d

  !     PURPOSE
  !         Calculates the extended Theis solution h(r) transient-state and heterogenous aquifers
  !         in 3 dimensions with anisotropical impact.

  !     CALLING SEQUENCE
  !         out = ext_theis3d(rad, time, params, inits, BC, prop)

  !     INTENT(IN)
  !         REAL(sp/dp),[DIMENSION(:)]  :: rad          Radius          : input radius                  in m
  !         REAL(sp/dp),[DIMENSION(:)]  :: time         Time            : input radius                  in m
  !         REAL(sp/dp), DIMENSION(5)   :: params       Parameters
  !                                                         params(1)   : geometric mean conductivity   in m^2/s
  !                                                         params(2)   : variance of log-conductivity
  !                                                         params(3)   : correlation length            in m
  !                                                         params(4)   : anisotropy ratio
  !                                                         params(5)   : aquifer thickness             in m
  !                                                         params(6)   : specific storage coefficient  in m^3
  !         REAL(sp/dp), DIMENSION(3)   :: inits        Initial values
  !                                                         inits(1)    : ref.-point for theis-solution in m
  !                                                         inits(3)    : pumping rate at the well      in m^3/s

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         REAL(sp/dp)                 :: prop         prop-factor for harmonic boundary-condition
  !                                                     the arithmetic case will be automatically corrected
  !                                                     default=1.6
  !         LOGICAL                     :: BC           kind of the boundary-condition
  !                                                     (harmonic: true[default]; arithmetic: false)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURNS
  !>       \return     REAL(sp/dp)      :: ext_theis3d      hydraulic head in m

  !     RESTRICTIONS
  !         1) Input values must be floating points.
  !         2) rad       >  0       [radius]
  !         3) time      >  0       [time]
  !         4) inits(1)  >  0       [reference radius]
  !         5) params(1) >  0       [conductivity]
  !         6) params(2) >  0       [variance]
  !         7) params(3) >  0       [correlation length]
  !         8) params(4) in [0,1]   [anisotropy ratio]
  !         9) params(5) >  0       [aquifer thickness]
  !        10) params(6) >  0       [storage coefficient]

  !     EXAMPLE
  !         rad     = 0.5
  !         time    = 50.0
  !         inits   = (/ 0.01 , -0.001 /)
  !         params  = (/0.0001, 1.0, 10.0, 1.0, 1.0, 0.001/)

  !         h       = ext_theis3d(rad, time, inits, params )

  !         -> see also example in test directory

  !     LITERATURE
  !         [1] ?

  !     HISTORY
  !         Written, Sebastian Mueller, June 2014

  INTERFACE ext_theis3d
     MODULE PROCEDURE ext_theis3d_sp, ext_theis3d_dp, ext_theis3d_sp_d1, ext_theis3d_dp_d1
  END INTERFACE ext_theis3d

  ! ------------------------------------------------------------------

  !     NAME
  !         ext_theis2d

  !     PURPOSE
  !         Calculates the extended Theis solution h(r) transient-state and heterogenous aquifers in 2 dimensions.

  !     CALLING SEQUENCE
  !         out = ext_theis2d(rad, time, params, inits, prop)

  !     INTENT(IN)
  !         REAL(sp/dp),[DIMENSION(:)]  :: rad          Radius          : input radius                  in m
  !         REAL(sp/dp),[DIMENSION(:)]  :: time         Time            : input radius                  in m
  !         REAL(sp/dp), DIMENSION(5)   :: params       Parameters
  !                                                         params(1)   : geometric mean transmissivity in m^2/s
  !                                                         params(2)   : variance of log-conductivity
  !                                                         params(3)   : correlation length            in m
  !                                                         params(6)   : specific storage coefficient  in m^3
  !         REAL(sp/dp), DIMENSION(3)   :: inits        Initial values
  !                                                         inits(1)    : ref.-point for theis-solution in m
  !                                                         inits(3)    : pumping rate at the well      in m^3/s

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         REAL(sp/dp)                 :: prop         proportionality-factor
  !                                                     default=2.0

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURNS
  !>       \return     REAL(sp/dp)      :: ext_theis2d      hydraulic head in m

  !     RESTRICTIONS
  !         1) Input values must be floating points.
  !         2) rad       >  0       [radius]
  !         3) time      >  0       [time]
  !         4) inits(1)  >  0       [reference radius]
  !         5) params(1) >  0       [transmissivity]
  !         6) params(2) >  0       [variance]
  !         7) params(3) >  0       [correlation length]
  !         8) params(4) >  0       [storage coefficient]

  !     EXAMPLE
  !         rad     = 0.5
  !         time    = 50.0
  !         inits   = (/ 0.01 , -0.001 /)
  !         params  = (/0.0001, 1.0, 10.0, 0.001/)

  !         h       = ext_theis2d(rad, time, inits, params )

  !         -> see also example in test directory

  !     LITERATURE
  !         [1] ?

  !     HISTORY
  !         Written, Sebastian Mueller, June 2014

  INTERFACE ext_theis2d
     MODULE PROCEDURE ext_theis2d_sp, ext_theis2d_dp, ext_theis2d_sp_d1, ext_theis2d_dp_d1
  END INTERFACE ext_theis2d


  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  FUNCTION thiem_sp(rad, params, inits)

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: rad
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: params               !params=(K,L)
    REAL(sp), DIMENSION(3),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw)

    REAL(sp)                                        :: thiem_sp

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (le(rad, 0.0_sp))            stop 'Radius must be positiv!'
    if (le(inits(1), 0.0_sp))       stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_sp))      stop 'The conductivity must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    thiem_sp        =   (-inits(3)/(2.0_sp*PI*params(2)*params(1))) * log(rad/inits(1)) + inits(2)

  END FUNCTION thiem_sp

  FUNCTION thiem_dp(rad, params, inits)

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: rad
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: params               !params=(K,L)
    REAL(dp), DIMENSION(3),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw)

    REAL(dp)                                        :: thiem_dp

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (le(rad, 0.0_dp))            stop 'Radius must be positiv!'
    if (le(inits(1), 0.0_dp))       stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_dp))      stop 'The conductivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    thiem_dp        =   (-inits(3)/(2.0_dp*PI_D*params(2)*params(1))) * log(rad/inits(1)) + inits(2)

  END FUNCTION thiem_dp

  FUNCTION thiem_d1_sp(rad, params, inits)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: params               !params=(K,L)
    REAL(sp), DIMENSION(3),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw)

    REAL(sp), DIMENSION(size(rad))                  :: thiem_d1_sp

    INTEGER(i4)                                     :: i

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_sp))      stop 'Radius must be positiv!'
    if (le(inits(1),  0.0_sp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_sp))      stop 'The conductivity must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do i=1_i4, size(rad)
       thiem_d1_sp(i)  =   thiem_sp(rad(i), params, inits)
    end do

  END FUNCTION thiem_d1_sp

  FUNCTION thiem_d1_dp(rad, params, inits)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: params               !params=(K,L)
    REAL(dp), DIMENSION(3),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw)

    REAL(dp), DIMENSION(size(rad))                  :: thiem_d1_dp

    INTEGER(i4)                                     :: i

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_dp))      stop 'Radius must be positiv!'
    if (le(inits(1),  0.0_dp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_dp))      stop 'The conductivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do i=1_i4, size(rad)
       thiem_d1_dp(i)  =   thiem_dp(rad(i), params, inits)
    end do

  END FUNCTION thiem_d1_dp

  FUNCTION theis_sp(rad, time, params, inits)

    USE mo_nr,                      ONLY: expint

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: rad, time
    REAL(sp), DIMENSION(3),             INTENT(IN)  :: params               !params=(S,K,L)
    REAL(sp), DIMENSION(1),             INTENT(IN)  :: inits                !inits=(Qw)

    REAL(sp)                                        :: theis_sp

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (le(rad, 0.0_sp))            stop 'The radius must be positiv.'
    if (time .LT. 0.0_sp)           stop 'The time must be non-negativ.'
    if (le(params(1), 0.0_sp))      stop 'The specific storage coefficient must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The conductivity must be positiv.'
    if (le(params(3), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    if (.not. eq(time, 0.0_sp)) then
       theis_sp    =   (inits(1)/(4.0_sp*PI*params(3)*params(2))) &
            * expint(1_i4,((rad**2_i4)*params(1))/(4.0_sp*params(3)*params(2)*time))
    else
       theis_sp    =   0.0_sp
    end if

  END FUNCTION theis_sp

  FUNCTION theis_dp(rad, time, params, inits)

    USE mo_nr,                      ONLY: expint

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: rad, time
    REAL(dp), DIMENSION(3),             INTENT(IN)  :: params               !params=(S,K,L)
    REAL(dp), DIMENSION(1),             INTENT(IN)  :: inits                !inits=(Qw)

    REAL(dp)                                        :: theis_dp

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (le(rad, 0.0_dp))            stop 'The radius must be positiv.'
    if (time .LT. 0.0_dp)           stop 'The time must be non-negativ.'
    if (le(params(1), 0.0_dp))      stop 'The specific storage coefficient must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The conductivity must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    if (.not. eq(time, 0.0_dp)) then
       theis_dp    =   (inits(1)/(4.0_dp*PI_D*params(3)*params(2))) &
            * expint(1_i4,((rad**2_i4)*params(1))/(4.0_dp*params(3)*params(2)*time))
    else
       theis_dp    =   0.0_dp
    end if

  END FUNCTION theis_dp

  FUNCTION theis_d1_sp(rad, time, params, inits)

    USE mo_nr,                      ONLY: expint

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(sp), DIMENSION(:),             INTENT(IN)  :: time
    REAL(sp), DIMENSION(3),             INTENT(IN)  :: params               !params=(S,K,L)
    REAL(sp), DIMENSION(1),             INTENT(IN)  :: inits                !inits=(Qw)

    REAL(sp), DIMENSION(size(rad),size(time))       :: theis_d1_sp

    INTEGER(i4)                                     :: i, j
    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_sp))      stop 'The radius must be positiv.'
    if (     any(time<0.0_sp))      stop 'The time must be non-negativ.'
    if (le(params(1), 0.0_sp))      stop 'The specific storage coefficient must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The conductivity must be positiv.'
    if (le(params(3), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do i=1_i4, size(rad)
       do j=1_i4, size(time)
          theis_d1_sp(i,j)    =   theis_sp(rad(i), time(j), params, inits)
       end do
    end do

  END FUNCTION theis_d1_sp

  FUNCTION theis_d1_dp(rad, time, params, inits)

    USE mo_nr,                      ONLY: expint

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: time
    REAL(dp), DIMENSION(3),             INTENT(IN)  :: params               !params=(S,K,L)
    REAL(dp), DIMENSION(1),             INTENT(IN)  :: inits                !inits=(Qw)

    REAL(dp), DIMENSION(size(rad),size(time))       :: theis_d1_dp

    INTEGER(i4)                                     :: i, j
    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_dp))      stop 'The radius must be positiv.'
    if (     any(time<0.0_dp))      stop 'The time must be non-negativ.'
    if (le(params(1), 0.0_dp))      stop 'The specific storage coefficient must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The conductivity must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do i=1_i4, size(rad)
       do j=1_i4, size(time)
          theis_d1_dp(i,j)    =   theis_dp(rad(i), time(j), params, inits)
       end do
    end do

  END FUNCTION theis_d1_dp

  FUNCTION ext_thiem2d_sp(rad, params, inits, prop)

    USE mo_nr,                      ONLY: ei, expint

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: rad
    REAL(sp), DIMENSION(3),             INTENT(IN)  :: params, inits        !params=(TG,var,corr); inits=(Ref,href,Qw)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality factor [default=2]

    REAL(sp)                                        :: ext_thiem2d_sp

    REAL(sp)                                        :: sub1, sub2, pr

    ! sub1/2:           are just substitutions to make the formula a bit shorter
    ! pr:               is the intern proportionality factor

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (le(rad, 0.0_sp))            stop 'Radius must be positiv.'
    if (le(inits(1), 0.0_sp))       stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_sp))      stop 'The transmissivity must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_sp))      stop 'The correlation length must be positiv.'

    if (present(prop)) then
       if (le(prop, 0.0_sp))       stop 'The proportionality factor must be positiv.'
       pr          =   prop
    else
       pr          =   2.0_sp
    endif

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    sub1            =   (params(2)/2.0_sp)/( 1.0_sp+(pr*rad/params(3))**2_i4 )
    sub2            =   (params(2)/2.0_sp)/( 1.0_sp+(pr*inits(1)/params(3))**2_i4 )

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    ext_thiem2d_sp  =   -(inits(3)/(4.0_sp*PI*params(1)))*&
         (&
         ( ei( ((pr*rad/params(3))**2_i4)*sub1 ) -&
         ei( ((pr*inits(1)/params(3))**2_i4)*sub2 )&
         ) + &
         exp(params(2)/2.0_sp)*&
         (&
         expint( 1_i4, sub1 ) -&
         expint( 1_i4, sub2 )&
         )&
         )&
         + inits(2)

  END FUNCTION ext_thiem2d_sp

  FUNCTION ext_thiem2d_dp(rad, params, inits, prop)

    USE mo_nr,                      ONLY: ei, expint

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: rad
    REAL(dp), DIMENSION(3),             INTENT(IN)  :: params, inits        !params=(TG,var,corr); inits=(Ref,href,Qw)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality factor [default=2]

    REAL(dp)                                        :: ext_thiem2d_dp

    REAL(dp)                                        :: sub1, sub2, pr

    ! sub1/2:           are just substitutions to make the formula a bit shorter
    ! pr:               is the intern proportionality factor

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (le(rad, 0.0_dp))            stop 'Radius must be positiv.'
    if (le(inits(1), 0.0_dp))       stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_dp))      stop 'The transmissivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The correlation length must be positiv.'

    if (present(prop)) then
       if (le(prop, 0.0_dp))       stop 'The proportionality factor must be positiv.'
       pr          =   prop
    else
       pr          =   2.0_dp
    endif

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    sub1            =   (params(2)/2.0_dp)/( 1.0_dp+(pr*rad/params(3))**2_i4 )
    sub2            =   (params(2)/2.0_dp)/( 1.0_dp+(pr*inits(1)/params(3))**2_i4 )

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    ext_thiem2d_dp  =   -(inits(3)/(4.0_dp*PI_D*params(1)))*&
         (&
         ( ei( ((pr*rad/params(3))**2_i4)*sub1 ) -&
         ei( ((pr*inits(1)/params(3))**2_i4)*sub2 )&
         ) + &
         exp(params(2)/2.0_dp)*&
         (&
         expint( 1_i4, sub1 ) -&
         expint( 1_i4, sub2 )&
         )&
         )&
         + inits(2)

  END FUNCTION ext_thiem2d_dp

  FUNCTION ext_thiem2d_d1_sp(rad, params, inits, prop)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(sp), DIMENSION(3),             INTENT(IN)  :: params, inits        !params=(TG,var,corr); inits=(Ref,href,Qw)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality factor [default=2]

    REAL(sp), DIMENSION(size(rad))                  :: ext_thiem2d_d1_sp

    REAL(sp)                                        :: pr
    INTEGER(i4)                                     :: i
    ! sub1/2:           are just substitutions to make the formula a bit shorter
    ! pr:               is the intern proportionality factor

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_sp))      stop 'Radius must be positiv.'
    if (le(inits(1),  0.0_sp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_sp))      stop 'The transmissivity must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_sp))      stop 'The correlation length must be positiv.'

    if (present(prop)) then
       if (le(prop, 0.0_sp))       stop 'The proportionality factor must be positiv.'
       pr          =   prop
    else
       pr          =   2.0_sp
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do i=1_i4, size(rad)
       ext_thiem2d_d1_sp(i)    =   ext_thiem2d_sp(rad(i), params, inits, pr)
    end do

  END FUNCTION ext_thiem2d_d1_sp

  FUNCTION ext_thiem2d_d1_dp(rad, params, inits, prop)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(dp), DIMENSION(3),             INTENT(IN)  :: params, inits        !params=(TG,var,corr); inits=(Ref,href,Qw)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality factor [default=2]

    REAL(dp), DIMENSION(size(rad))                  :: ext_thiem2d_d1_dp

    REAL(dp)                                        :: pr
    INTEGER(i4)                                     :: i
    ! sub1/2:           are just substitutions to make the formula a bit shorter
    ! pr:               is the intern proportionality factor

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_dp))      stop 'Radius must be positiv.'
    if (le(inits(1),  0.0_dp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_dp))      stop 'The transmissivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The correlation length must be positiv.'

    if (present(prop)) then
       if (le(prop, 0.0_dp))       stop 'The proportionality factor must be positiv.'
       pr          =   prop
    else
       pr          =   2.0_dp
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do i=1_i4, size(rad)
       ext_thiem2d_d1_dp(i)    =   ext_thiem2d_dp(rad(i), params, inits, pr)
    end do

  END FUNCTION ext_thiem2d_d1_dp

  FUNCTION ext_thiem3d_sp(rad, params, inits, BC, prop)

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: rad
    REAL(sp), DIMENSION(5),             INTENT(IN)  :: params               !params=(KG,var,corr,e,L)
    REAL(sp), DIMENSION(3),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]

    REAL(sp)                                        :: ext_thiem3d_sp

    REAL(sp)                                        :: sub11, sub12, sub21, sub22, pr, aniso, Kwell, Kefu, chi
    LOGICAL                                         :: BCkind

    ! sub11/12/21/22:   are just substitutions to make the formula a bit shorter
    ! pr:               is the intern proportionality factor
    ! aniso:            is the anisotropy function
    ! Kwell:            conductivity at the well
    ! Kefu:             far field conductivity for uniform flow
    ! chi:              abbreviation for ln(Kefu/Kwell)
    ! BCkind:           is the intern kind of the boundary condition

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (le(rad      , 0.0_sp))      stop 'Radius must be positiv.'
    if (le(inits(1) , 0.0_sp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_sp))      stop 'The geometric mean conductivity must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_sp))      stop 'The correlation length must be positiv.'
    if (params(4) .LT. 0.0_sp)      stop 'The anisotropy ratio must be in [0,1].'
    if (params(4) .GT. 1.0_sp)      stop 'The anisotropy ratio must be in [0,1].'
    if (le(params(5), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

    if (present(BC)) then
       BCkind = BC
    else
       BCkind = .TRUE.
    endif

    if (present(prop)) then
       if (le(prop, 0.0_sp))       stop 'The proportionality factor must be positiv.'
       pr          =   prop
    else
       pr          =   1.6_sp
    endif

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    if (.NOT. BCkind)   pr=0.5_sp*pr

    if (eq(params(4), 1.0_sp)) then
       aniso       =   1.0_sp/3.0_sp
    else if (eq(params(4), 0.0_sp)) then
       aniso       =   0.0_sp
    else
       aniso       =   (params(4)/2.0_sp) * &
            (&
            ( (1.0_sp-params(4)**2_i4)**(-1.5_sp) )*atan((params(4)**2_i4 - 1.0_sp)**(-0.5_sp)) &
            - params(4)/(1.0_sp-params(4)**2_i4) &
            )
    endif

    Kefu            =   params(1)*exp(params(2)*(0.5_sp - aniso))

    if (BCkind) then
       Kwell       =   params(1)*exp(-params(2)/2.0_sp)
    else
       Kwell       =   params(1)*exp(params(2)/2.0_sp)
    endif

    chi             =   log(Kwell/Kefu)

    if (eq(params(4), 0.0_sp)) then
       sub21       =   log(rad/inits(1))
       sub22       =   log(rad/inits(1))
    else
       sub11       =   (1 + ( (pr/(params(3)*(params(4)**(1.0_sp/3.0_sp))) )**2_i4)*(inits(1)**2_i4)   )**0.5_sp
       sub12       =   (1 + ( (pr/(params(3)*(params(4)**(1.0_sp/3.0_sp))) )**2_i4)*(rad**2_i4)        )**0.5_sp

       sub21       =   log((sub12+1.0_sp)/(sub11+1.0_sp)) - 1.0_sp/sub12 + 1.0_sp/sub11
       sub22       =   log(sub12/sub11) - 1.0_sp/(2.0_sp*(sub12**2_i4)) + 1.0_sp/(2.0_sp*(sub11**2_i4)) &
            - 1.0_sp/(4.0_sp*(sub12**4_i4)) + 1.0_sp/(4.0_sp*(sub11**4_i4))
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    ext_thiem3d_sp  =   -(inits(3)/(2.0_sp*PI*params(5)*Kefu))*&
         (&
         exp(-chi)*log(rad/inits(1)) + sinh(chi)*sub21 + (1.0_sp - cosh(chi))*sub22 &
         )&
         + inits(2)

  END FUNCTION ext_thiem3d_sp

  FUNCTION ext_thiem3d_dp(rad, params, inits, BC, prop)

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: rad
    REAL(dp), DIMENSION(5),             INTENT(IN)  :: params               !params=(KG,var,corr,e,L)
    REAL(dp), DIMENSION(3),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]

    REAL(dp)                                        :: ext_thiem3d_dp

    REAL(dp)                                        :: sub11, sub12, sub21, sub22, pr, aniso, Kwell, Kefu, chi
    LOGICAL                                         :: BCkind

    ! sub11/12/21/22:   are just substitutions to make the formula a bit shorter
    ! pr:               is the intern proportionality factor
    ! aniso:            is the anisotropy function
    ! Kwell:            conductivity at the well
    ! Kefu:             far field conductivity for uniform flow
    ! chi:              abbreviation for ln(Kefu/Kwell)
    ! BCkind:           is the intern kind of the boundary condition

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (le(rad      , 0.0_dp))      stop 'Radius must be positiv.'
    if (le(inits(1) , 0.0_dp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_dp))      stop 'The geometric mean conductivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The correlation length must be positiv.'
    if (params(4) .LT. 0.0_dp)      stop 'The anisotropy ratio must be in [0,1].'
    if (params(4) .GT. 1.0_dp)      stop 'The anisotropy ratio must be in [0,1].'
    if (le(params(5), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

    if (present(BC)) then
       BCkind = BC
    else
       BCkind = .TRUE.
    endif

    if (present(prop)) then
       if (le(prop, 0.0_dp))       stop 'The proportionality factor must be positiv.'
       pr          =   prop
    else
       pr          =   1.6_dp      !Standard-value
    endif

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    if (.NOT. BCkind)   pr=0.5_dp*pr

    if (eq(params(4), 1.0_dp)) then
       aniso       =   1.0_dp/3.0_dp
    else if (eq(params(4), 0.0_dp)) then
       aniso       =   0.0_dp
    else
       aniso       =   (params(4)/2.0_dp) * &
            (&
            ( (1.0_dp-params(4)**2_i4)**(-1.5_dp) )*atan((params(4)**2_i4 - 1.0_dp)**(-0.5_dp)) &
            - params(4)/(1.0_dp-params(4)**2_i4) &
            )
    endif

    Kefu            =   params(1)*exp(params(2)*(0.5_dp - aniso))

    if (BCkind) then
       Kwell       =   params(1)*exp(-params(2)/2.0_dp)
    else
       Kwell       =   params(1)*exp(params(2)/2.0_dp)
    endif

    chi             =   log(Kwell/Kefu)


    if (eq(params(4), 0.0_dp)) then
       sub21       =   log(rad/inits(1))
       sub22       =   log(rad/inits(1))
    else
       sub11       =   (1 + ( (pr/(params(3)*(params(4)**(1.0_dp/3.0_dp))) )**2_i4)*(inits(1)**2_i4)   )**0.5_dp
       sub12       =   (1 + ( (pr/(params(3)*(params(4)**(1.0_dp/3.0_dp))) )**2_i4)*(rad**2_i4)        )**0.5_dp

       sub21       =   log((sub12+1.0_dp)/(sub11+1.0_dp)) - 1.0_dp/sub12 + 1.0_dp/sub11
       sub22       =   log(sub12/sub11) - 1.0_dp/(2.0_dp*(sub12**2_i4)) + 1.0_dp/(2.0_dp*(sub11**2_i4)) &
            - 1.0_dp/(4.0_dp*(sub12**4_i4)) + 1.0_dp/(4.0_dp*(sub11**4_i4))
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    ext_thiem3d_dp  =   -(inits(3)/(2.0_dp*PI_D*params(5)*Kefu))*&
         (&
         exp(-chi)*log(rad/inits(1)) + sinh(chi)*sub21 + (1.0_dp - cosh(chi))*sub22 &
         )&
         + inits(2)

  END FUNCTION ext_thiem3d_dp

  FUNCTION ext_thiem3d_d1_sp(rad, params, inits, BC, prop)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(sp), DIMENSION(5),             INTENT(IN)  :: params               !params=(KG,var,corr,e,L)
    REAL(sp), DIMENSION(3),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]

    REAL(sp), DIMENSION(size(rad))                  :: ext_thiem3d_d1_sp

    REAL(sp)                                        :: pr
    integer(i4)                                     :: i
    LOGICAL                                         :: BCkind

    ! pr:               is the intern proportionality factor
    ! BCkind:           is the intern kind of the boundary condition

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_sp))      stop 'Radius must be positiv.'
    if (le(inits(1) , 0.0_sp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_sp))      stop 'The geometric mean conductivity must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_sp))      stop 'The correlation length must be positiv.'
    if (params(4) .LT. 0.0_sp)      stop 'The anisotropy ratio must be in [0,1].'
    if (params(4) .GT. 1.0_sp)      stop 'The anisotropy ratio must be in [0,1].'
    if (le(params(5), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

    if (present(BC)) then
       BCkind = BC
    else
       BCkind = .TRUE.
    endif

    if (present(prop)) then
       if (le(prop, 0.0_sp))       stop 'The proportionality factor must be positiv.'
       pr          =   prop
    else
       pr          =   1.6_sp      !Standard-value
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do i=1_i4, size(rad)
       ext_thiem3d_d1_sp(i)    =   ext_thiem3d_sp(rad(i), params, inits, BCkind, pr)
    end do

  END FUNCTION ext_thiem3d_d1_sp

  FUNCTION ext_thiem3d_d1_dp(rad, params, inits, BC, prop)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(dp), DIMENSION(5),             INTENT(IN)  :: params               !params=(KG,var,corr,e,L)
    REAL(dp), DIMENSION(3),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]

    REAL(dp), DIMENSION(size(rad))                  :: ext_thiem3d_d1_dp

    REAL(dp)                                        :: pr
    integer(i4)                                     :: i
    LOGICAL                                         :: BCkind

    ! pr:               is the intern proportionality factor
    ! BCkind:           is the intern kind of the boundary condition

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_dp))      stop 'Radius must be positiv.'
    if (le(inits(1) , 0.0_dp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_dp))      stop 'The geometric mean conductivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The correlation length must be positiv.'
    if (params(4) .LT. 0.0_dp)      stop 'The anisotropy ratio must be in [0,1].'
    if (params(4) .GT. 1.0_dp)      stop 'The anisotropy ratio must be in [0,1].'
    if (le(params(5), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

    if (present(BC)) then
       BCkind = BC
    else
       BCkind = .TRUE.
    endif

    if (present(prop)) then
       if (le(prop, 0.0_dp))       stop 'The proportionality factor must be positiv.'
       pr          =   prop
    else
       pr          =   1.6_dp      !Standard-value
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do i=1_i4, size(rad)
       ext_thiem3d_d1_dp(i)    =   ext_thiem3d_dp(rad(i), params, inits, BCkind, pr)
    end do

  END FUNCTION ext_thiem3d_d1_dp

  ! extended theis' solution

  FUNCTION ext_theis3d_sp(rad, time, params, inits, BC, prop, grad, steh_n)

    USE mo_laplace_inversion, ONLY: NLInvSteh

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: rad
    REAL(sp),                           INTENT(IN)  :: time
    REAL(sp), DIMENSION(6),             INTENT(IN)  :: params               !params=(KG,var,corr,e,L,S)
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Ref,Qw)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=30]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]

    REAL(sp)                                        :: ext_theis3d_sp

    REAL(sp)                                        :: pr, aniso, Kwell, Kefu, chi, Kref, C
    REAL(sp), DIMENSION(10)                         :: values
    INTEGER(i4)                                     :: stehlim, limit
    LOGICAL                                         :: BCkind

    ! aniso:            is the anisotropy function
    ! Kwell:            conductivity at the well
    ! Kefu:             far field conductivity for uniform flow
    ! Kref:             conductivity at the reference point (Ref=inits(1))
    ! chi:              abbreviation for ln(Kefu/Kwell)
    ! C:                abbreviation for the factor C=prop/(e^(1/3)*corr)
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! BCkind:           is the intern kind of the boundary condition
    ! pr:               is the intern proportionality factor
    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (le(rad      , 0.0_sp))      stop 'Radius must be positiv.'
    if (time .LT.     0.0_sp)       stop 'The time must be non-negativ.'
    if (le(inits(1) , 0.0_sp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_sp))      stop 'The geometric mean conductivity must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_sp))      stop 'The correlation length must be positiv.'
    if (params(4) .LT. 0.0_sp)      stop 'The anisotropy ratio must be in [0,1].'
    if (params(4) .GT. 1.0_sp)      stop 'The anisotropy ratio must be in [0,1].'
    if (le(params(5), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

    if (present(BC)) then
       BCkind = BC
    else
       BCkind = .TRUE.
    endif

    if (present(prop)) then
       if (le(prop, 0.0_sp))        stop 'The proportionality factor must be positiv.'
       pr           =   prop
    else
       pr           =   1.6_sp
    endif

    if (present(steh_n)) then
       if (steh_n<=0_i4)            stop 'The Stehfest-boundary must be positiv.'
       stehlim      =   steh_n
    else
       stehlim      =   6_i4
    endif

    if (present(grad)) then
       if (grad<=3_i4)              stop 'The limit of the series expansion must be at least 4.'
       limit        =   grad
    else
       limit        =   30_i4
    endif

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    if (.NOT. BCkind)   pr=0.5_sp*pr

    C               =   pr/( (params(4)**(1.0_sp/3.0_sp))*params(3) )

    if (.not. (rad<C**(-1.0_sp)))   stop 'Radius out of convergence-area.'

    if (eq(params(4), 1.0_sp)) then
       aniso       =   1.0_sp/3.0_sp
    else if (eq(params(4), 0.0_sp)) then
       aniso       =   0.0_sp
    else
       aniso       =   (params(4)/2.0_sp) * &
            (&
            ( (1.0_sp-params(4)**2_i4)**(-1.5_sp) )*atan((params(4)**2_i4 - 1.0_sp)**(-0.5_sp)) &
            - params(4)/(1.0_sp-params(4)**2_i4) &
            )
    endif

    Kefu            =   params(1)*exp(params(2)*(0.5_sp - aniso))

    if (BCkind) then
       Kwell       =   params(1)*exp(-params(2)/2.0_sp)
    else
       Kwell       =   params(1)*exp(params(2)/2.0_sp)
    endif

    chi             =   log(Kwell/Kefu)

    Kref            =   Kefu*exp( chi*( sqrt(1.0_sp+(C*inits(1))**2_i4) )**(-3.0_sp) )

    values          =   (/Kefu, chi, C, params(6), inits(2), params(5), Kwell, inits(1), Kref, real(limit,sp)/)

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    if (time>0.0_sp) then
       ext_theis3d_sp  =   min(NLInvSteh(lap_ext_theis3d_sp, values, rad, time, stehlim),0.0_sp)
    else
       ext_theis3d_sp  =   0.0_sp
    end if

  END FUNCTION ext_theis3d_sp

  FUNCTION lap_ext_theis3d_sp(r, s, para)

    use mo_functions,       only: factorial
    use mo_nrutil,          only: poly
    use mo_combinatorics,   only: binomcoeffi, nextpart
    use mo_nr,              only: bessk0

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: r
    REAL(sp),                           INTENT(IN)  :: s
    REAL(sp), DIMENSION(:),             INTENT(IN)  :: para                     !para=(Kefu,chi,C,S,Qw,L,Kwell,Ref,Kref,limit)

    REAL(sp)                                        :: lap_ext_theis3d_sp

    !intern variables
    integer(i4), Parameter                          :: stdsiz=31_i4

    REAL(sp), &
         DIMENSION(max(nint(para(10),i4)+1,stdsiz))     :: g, h, alpha, beta, innerfunc   !expansions of the given functions

    integer(i4), &
         DIMENSION((size(beta)*(size(beta)-1_i4))/2_i4) :: partitions               !All partitionarrays in a row to calculate beta

    REAL(sp)                                        :: prod, sub1, sub2, sub3, coef

    INTEGER(i4)                                     :: n,m

    !Calculate the series-expansion of alpha=1-3chi*(C*r)^2*(1+(C*r)^2)^(5/2)
    alpha           =   0.0_sp
    alpha(1)        =   1.0_sp

    do n=0_i4, min(Floor(para(10)/2.0_sp - 1.0_sp), ((stdsiz-1_i4)/2_i4 - 1_i4))
       alpha(2*(n+1)+1)=   -3*para(2)*binomcoeffi(-2.5_sp,n)*para(3)**(2*(n+1))
    end do

    !Calculate the series-expansion of beta=-r^2*S/KCG(r)
    beta            =   0.0_sp
    partitions      =   0_i4

    !set all partitions to the "biggest" one
    do n=1_i4,size(beta)-1_i4
       partitions((n*(n+1_i4))/2_i4)   =   1_i4
    end do

    !Define the series expansion of the inner function
    innerfunc       =   0.0_sp
    innerfunc(1)    =   -para(2)

    do n=1_i4, min(Floor(para(10)/2.0_sp), ((stdsiz-1_i4)/2_i4))
       innerfunc(2_i4*n+1_i4)  =   -para(2)*binomcoeffi(-1.5_sp,n)*para(3)**(2_i4*n)
    end do

    beta(3)         =   1.0_sp

    !Define the series expansion of beta by the rule of Faa di Bruno
    do n=2_i4,min(2_i4*Floor(para(10)/2.0_sp), stdsiz-1_i4),2_i4
       do
          prod=1.0_sp
          do m=1_i4, n

             if ((partitions(((n-1)*n)/2+m) .ne. 0_i4)) then
                if  (ne(innerfunc(m+1), 0.0_sp)) then
                   prod    =&
                        prod*(innerfunc(m+1)**real(partitions(((n-1)*n)/2+m),sp))/real(factorial(partitions(((n-1)*n)/2+m)),sp)
                else
                   prod    =   0.0_sp
                   exit
                end if
             end if

          end do

          if (.not. (isnan(prod) .or. (abs(prod) > huge(1.0_sp)) .or. (abs(prod) < tiny(1.0_sp))))&
               beta(n+3)=beta(n+3)+prod

          if (partitions(((n-1)*n)/2+1)==n) exit

          partitions(((n-1)*n)/2+1:((n+1)*n)/2)    =   nextpart(partitions(((n-1)*n)/2+1:((n+1)*n)/2))
       end do

    end do

    beta     =   -(para(4)/para(1))*exp(-para(2))*beta

    !Define the series expansion of the first fundamental solution g with the aid of alpha and beta by the frobenius-method
    g               =   0.0_sp
    g(1)            =   1.0_sp

    do n=1_i4,size(g)-1_i4
       do m=0_i4,n-1_i4
          g(n+1)  =   g(n+1)+(real(m,sp)*alpha(n+1-m)+s*beta(n+1-m))*g(m+1)
       end do
       g(n+1)      =   -1.0_sp/(real(n,sp)**2_i4)*g(n+1)
    end do

    !Define the series of the the first part of the second fundamental solution h with the aid of alpha, beta and g by the frobenius-m.
    h               =   0.0_sp
    h(1)            =   1.0_sp

    do n=1_i4,size(h)-1_i4
       do m=0_i4,n-1_i4
          h(n+1)  =   h(n+1) + alpha(n+1-m)*(real(m,sp)*h(m+1)+g(m+1)) + s*beta(n+1-m)*h(m+1)
       end do
       h(n+1)      =   -1.0_sp/(real(n,sp)**2_i4)*(h(n+1)+2.0_sp*real(n,sp)*g(n+1))
    end do

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    coef                =   para(5)/(2.0_sp*PI*para(6)*para(7))

    sub3                =   (para(7)/para(9))*coef/s*bessk0(sqrt(para(8)**2_i4*para(4)*s/(para(6)*para(9))))
    !swap the nearfield conductivity with the reference-conductivity Kwell<->Kref in coef
    sub2                =   -coef/s

    sub1                =   (sub3-sub2*( poly(para(8),h)+log(para(8))*poly(para(8),g) ))/&
         poly(para(8),g)

    lap_ext_theis3d_sp  =   (sub1+sub2*log(r))*poly(r,g) + sub2*poly(r,h)


  END FUNCTION lap_ext_theis3d_sp

  FUNCTION ext_theis3d_dp(rad, time, params, inits, BC, prop, grad, steh_n)

    use mo_laplace_inversion,    ONLY: NLInvSteh

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: rad
    REAL(dp),                           INTENT(IN)  :: time
    REAL(dp), DIMENSION(6),             INTENT(IN)  :: params               !params=(KG,var,corr,e,L,S)
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Ref,Qw)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=30]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]

    REAL(dp)                                        :: ext_theis3d_dp

    REAL(dp)                                        :: pr, aniso, Kwell, Kefu, chi, Kref, C
    REAL(dp), DIMENSION(10)                         :: values
    INTEGER(i4)                                     :: stehlim, limit
    LOGICAL                                         :: BCkind

    ! aniso:            is the anisotropy function
    ! Kwell:            conductivity at the well
    ! Kefu:             far field conductivity for uniform flow
    ! Kref:             conductivity at the reference point (Ref=inits(1))
    ! chi:              abbreviation for ln(Kefu/Kwell)
    ! C:                abbreviation for the factor C=prop/(e^(1/3)*corr)
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! BCkind:           is the intern kind of the boundary condition
    ! pr:               is the intern proportionality factor
    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (le(rad      , 0.0_dp))      stop 'Radius must be positiv.'
    if (time .LT.     0.0_dp)       stop 'The time must be non-negativ.'
    if (le(inits(1) , 0.0_dp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_dp))      stop 'The geometric mean conductivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The correlation length must be positiv.'
    if (params(4) .LT. 0.0_dp)      stop 'The anisotropy ratio must be in [0,1].'
    if (params(4) .GT. 1.0_dp)      stop 'The anisotropy ratio must be in [0,1].'
    if (le(params(5), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

    if (present(BC)) then
       BCkind = BC
    else
       BCkind = .TRUE.
    endif

    if (present(prop)) then
       if (le(prop, 0.0_dp))        stop 'The proportionality factor must be positiv.'
       pr           =   prop
    else
       pr           =   1.6_dp
    endif

    if (present(steh_n)) then
       if (steh_n<=0_i4)            stop 'The Stehfest-boundary must be positiv.'
       stehlim      =   steh_n
    else
       stehlim      =   6_i4
    endif

    if (present(grad)) then
       if (grad<=3_i4)              stop 'The limit of the series expansion must be at least 4.'
       limit        =   grad
    else
       limit        =   30_i4
    endif

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    if (.NOT. BCkind)   pr=0.5_dp*pr

    C               =   pr/( (params(4)**(1.0_dp/3.0_dp))*params(3) )

    if (.not. (rad<C**(-1.0_dp)))   stop 'Radius out of convergence-area.'

    if (eq(params(4), 1.0_dp)) then
       aniso       =   1.0_dp/3.0_dp
    else if (eq(params(4), 0.0_dp)) then
       aniso       =   0.0_dp
    else
       aniso       =   (params(4)/2.0_dp) * &
            (&
            ( (1.0_dp-params(4)**2_i4)**(-1.5_dp) )*atan((params(4)**2_i4 - 1.0_dp)**(-0.5_dp)) &
            - params(4)/(1.0_dp-params(4)**2_i4) &
            )
    endif

    Kefu            =   params(1)*exp(params(2)*(0.5_dp - aniso))

    if (BCkind) then
       Kwell       =   params(1)*exp(-params(2)/2.0_dp)
    else
       Kwell       =   params(1)*exp(params(2)/2.0_dp)
    endif

    chi             =   log(Kwell/Kefu)

    Kref            =   Kefu*exp( chi*( sqrt(1.0_dp+(C*inits(1))**2_i4) )**(-3.0_dp) )

    values          =   (/Kefu, chi, C, params(6), inits(2), params(5), Kwell, inits(1), Kref, real(limit,dp)/)

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    if (time>0.0_dp) then
       ext_theis3d_dp  =   min(NLInvSteh(lap_ext_theis3d_dp, values, rad, time, stehlim),0.0_dp)
    else
       ext_theis3d_dp  =   0.0_dp
    end if

  END FUNCTION ext_theis3d_dp

  FUNCTION lap_ext_theis3d_dp(r, s, para)

    use mo_functions,       only: factorial
    use mo_nrutil,          only: poly
    use mo_combinatorics,   only: binomcoeffi, nextpart
    use mo_nr,              only: bessk0

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: r
    REAL(dp),                           INTENT(IN)  :: s
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: para                     !para=(Kefu,chi,C,S,Qw,L,Kwell,Ref,Kref,limit)

    REAL(dp)                                        :: lap_ext_theis3d_dp

    !intern variables
    integer(i4), Parameter                          :: stdsiz=51_i4

    REAL(dp), &
         DIMENSION(max(nint(para(10),i4)+1,stdsiz))     :: g, h, alpha, beta, innerfunc   !expansions of the given functions

    integer(i4), &
         DIMENSION((size(beta)*(size(beta)-1_i4))/2_i4) :: partitions               !All partitionarrays in a row to calculate beta

    REAL(dp)                                        :: prod, sub1, sub2, sub3, coef

    INTEGER(i4)                                     :: n,m

    !Calculate the series-expansion of alpha=1-3chi*(C*r)^2*(1+(C*r)^2)^(5/2)
    alpha           =   0.0_dp
    alpha(1)        =   1.0_dp

    do n=0_i4, min(Floor(para(10)/2.0_dp - 1.0_dp), ((stdsiz-1_i4)/2_i4 - 1_i4))
       alpha(2*(n+1)+1)=   -3*para(2)*binomcoeffi(-2.5_dp,n)*para(3)**(2*(n+1))
    end do

    !Calculate the series-expansion of beta=-r^2*S/KCG(r)
    beta            =   0.0_dp
    partitions      =   0_i4

    !set all partitions to the "biggest" one
    do n=1_i4,size(beta)-1_i4
       partitions((n*(n+1_i4))/2_i4)   =   1_i4
    end do

    !Define the series expansion of the inner function
    innerfunc       =   0.0_dp
    innerfunc(1)    =   -para(2)

    do n=1_i4, min(Floor(para(10)/2.0_dp), ((stdsiz-1_i4)/2_i4))
       innerfunc(2_i4*n+1_i4)  =   -para(2)*binomcoeffi(-1.5_dp,n)*para(3)**(2_i4*n)
    end do

    beta(3)         =   1.0_dp

    !Define the series expansion of beta by the rule of Faa di Bruno
    do n=2_i4,min(2_i4*Floor(para(10)/2.0_dp), stdsiz-1_i4),2_i4
       do
          prod=1.0_dp
          do m=1_i4, n

             if ((partitions(((n-1)*n)/2+m) .ne. 0_i4)) then
                if  (ne(innerfunc(m+1), 0.0_dp)) then
                   prod    =&
                        prod*(innerfunc(m+1)**real(partitions(((n-1)*n)/2+m),dp))/real(factorial(partitions(((n-1)*n)/2+m)),dp)
                else
                   prod    =   0.0_dp
                   exit
                end if
             end if

          end do

          if (.not. (isnan(prod) .or. (abs(prod) > huge(1.0_dp)) .or. (abs(prod) < tiny(1.0_dp))))&
               beta(n+3)=beta(n+3)+prod

          if (partitions(((n-1)*n)/2+1)==n) exit

          partitions(((n-1)*n)/2+1:((n+1)*n)/2)    =   nextpart(partitions(((n-1)*n)/2+1:((n+1)*n)/2))
       end do

    end do

    beta     =   -(para(4)/para(1))*exp(-para(2))*beta

    !Define the series expansion of the first fundamental solution g with the aid of alpha and beta by the frobenius-method
    g               =   0.0_dp
    g(1)            =   1.0_dp

    do n=1_i4,size(g)-1_i4
       do m=0_i4,n-1_i4
          g(n+1)  =   g(n+1)+(real(m,dp)*alpha(n+1-m)+s*beta(n+1-m))*g(m+1)
       end do
       g(n+1)      =   -1.0_dp/(real(n,dp)**2_i4)*g(n+1)
    end do

    !Define the series of the the first part of the second fundamental solution h with the aid of alpha, beta and g by the frobenius-m.
    h               =   0.0_dp
    h(1)            =   1.0_dp

    do n=1_i4,size(h)-1_i4
       do m=0_i4,n-1_i4
          h(n+1)  =   h(n+1) + alpha(n+1-m)*(real(m,dp)*h(m+1)+g(m+1)) + s*beta(n+1-m)*h(m+1)
       end do
       h(n+1)      =   -1.0_dp/(real(n,dp)**2_i4)*(h(n+1)+2.0_dp*real(n,dp)*g(n+1))
    end do

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    coef                =   para(5)/(2.0_dp*PI_D*para(6)*para(7))

    sub3                =   (para(7)/para(9))*coef/s*bessk0(sqrt(para(8)**2_i4*para(4)*s/(para(6)*para(9))))
    !swap the nearfield conductivity with the reference-conductivity Kwell<->Kref in coef
    sub2                =   -coef/s

    sub1                =   (sub3-sub2*( poly(para(8),h)+log(para(8))*poly(para(8),g) ))/&
         poly(para(8),g)

    lap_ext_theis3d_dp  =   (sub1+sub2*log(r))*poly(r,g) + sub2*poly(r,h)


  END FUNCTION lap_ext_theis3d_dp

  FUNCTION ext_theis3d_sp_d1(rad, time, params, inits, BC, prop, grad, steh_n)

    use mo_laplace_inversion,    ONLY: NLInvSteh

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(sp), DIMENSION(:),             INTENT(IN)  :: time
    REAL(sp), DIMENSION(6),             INTENT(IN)  :: params               !params=(KG,var,corr,e,L,S)
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Ref,Qw)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=30]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]

    REAL(sp), DIMENSION(size(rad),size(time))       :: ext_theis3d_sp_d1

    REAL(sp)                                        :: pr, aniso, Kwell, Kefu, chi, Kref, C
    REAL(sp), DIMENSION(10)                         :: values
    INTEGER(i4)                                     :: stehlim, limit, n
    LOGICAL                                         :: BCkind

    LOGICAL, DIMENSION(size(rad), size(time))       :: mask

    ! aniso:            is the anisotropy function
    ! Kwell:            conductivity at the well
    ! Kefu:             far field conductivity for uniform flow
    ! Kref:             conductivity at the reference point (Ref=inits(1))
    ! chi:              abbreviation for ln(Kefu/Kwell)
    ! C:                abbreviation for the factor C=prop/(e^(1/3)*corr)
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! BCkind:           is the intern kind of the boundary condition
    ! pr:               is the intern proportionality factor
    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! mask:             is the mask for all non-zero time-points

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_sp))      stop 'Radius must be positiv.'
    if (any(time <    0.0_sp))      stop 'The time must be non-negativ.'
    if (le(inits(1) , 0.0_sp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_sp))      stop 'The geometric mean conductivity must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_sp))      stop 'The correlation length must be positiv.'
    if (params(4) .LT. 0.0_sp)      stop 'The anisotropy ratio must be in [0,1].'
    if (params(4) .GT. 1.0_sp)      stop 'The anisotropy ratio must be in [0,1].'
    if (le(params(5), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

    if (present(BC)) then
       BCkind = BC
    else
       BCkind = .TRUE.
    endif

    if (present(prop)) then
       if (le(prop, 0.0_sp))        stop 'The proportionality factor must be positiv.'
       pr           =   prop
    else
       pr           =   1.6_sp
    endif

    if (present(steh_n)) then
       if (steh_n<=0_i4)            stop 'The Stehfest-boundary must be positiv.'
       stehlim      =   steh_n
    else
       stehlim      =   6_i4
    endif

    if (present(grad)) then
       if (grad<=3_i4)              stop 'The limit of the series expansion must be at least 4.'
       limit        =   grad
    else
       limit        =   30_i4
    endif

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    if (.NOT. BCkind)   pr=0.5_sp*pr

    C               =   pr/( (params(4)**(1.0_sp/3.0_sp))*params(3) )

    if (.not. all(rad<C**(-1.0_sp)))   stop 'Radius out of convergence-area.'

    if (eq(params(4), 1.0_sp)) then
       aniso       =   1.0_sp/3.0_sp
    else if (eq(params(4), 0.0_sp)) then
       aniso       =   0.0_sp
    else
       aniso       =   (params(4)/2.0_sp) * &
            (&
            ( (1.0_sp-params(4)**2_i4)**(-1.5_sp) )*atan((params(4)**2_i4 - 1.0_sp)**(-0.5_sp)) &
            - params(4)/(1.0_sp-params(4)**2_i4) &
            )
    endif

    Kefu            =   params(1)*exp(params(2)*(0.5_sp - aniso))

    if (BCkind) then
       Kwell       =   params(1)*exp(-params(2)/2.0_sp)
    else
       Kwell       =   params(1)*exp(params(2)/2.0_sp)
    endif

    chi             =   log(Kwell/Kefu)

    Kref            =   Kefu*exp( chi*( sqrt(1.0_sp+(C*inits(1))**2_i4) )**(-3.0_sp) )

    values          =   (/Kefu, chi, C, params(6), inits(2), params(5), Kwell, inits(1), Kref, real(limit,sp)/)

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do n=1_i4, size(rad)
       mask(n,:)       =   (time>0.0_sp)
    end do

    ext_theis3d_sp_d1   =   0.0_sp

    ext_theis3d_sp_d1   =&
         unpack( pack(NLInvSteh(lap_ext_theis3d_sp_d1, values, rad, pack(time, (time>0.0_sp)), stehlim), .true.),&
         mask,&
         ext_theis3d_sp_d1)

    where( ext_theis3d_sp_d1 .gt. 0.0_sp )  ext_theis3d_sp_d1   =   0.0_sp

  END FUNCTION ext_theis3d_sp_d1

  FUNCTION lap_ext_theis3d_sp_d1(r, s, para)

    use mo_functions,       only: factorial
    use mo_nrutil,          only: poly
    use mo_combinatorics,   only: binomcoeffi, nextpart
    use mo_nr,              only: bessk0

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: r
    REAL(sp), DIMENSION(:),             INTENT(IN)  :: s
    REAL(sp), DIMENSION(:),             INTENT(IN)  :: para                     !para=(Kefu,chi,C,S,Qw,L,Kwell,Ref,Kref,limit)

    REAL(sp), DIMENSION(size(r), size(s))           :: lap_ext_theis3d_sp_d1

    !intern variables
    integer(i4), Parameter                          :: stdsiz=31_i4

    REAL(sp), &
         DIMENSION(max(nint(para(10),i4)+1,stdsiz))     :: g, h, alpha, beta, innerfunc   !expansions of the given functions

    integer(i4), &
         DIMENSION((size(beta)*(size(beta)-1_i4))/2_i4) :: partitions               !All partitionarrays in a row to calculate beta

    REAL(sp)                                        :: prod, sub1, sub2, sub3, coef

    INTEGER(i4)                                     :: n, m, rad_n, laps_n

    !Calculate the series-expansion of alpha=1-3chi*(C*r)^2*(1+(C*r)^2)^(5/2)
    alpha           =   0.0_sp
    alpha(1)        =   1.0_sp

    do n=0_i4, min(Floor(para(10)/2.0_sp - 1.0_sp), ((stdsiz-1_i4)/2_i4 - 1_i4))
       alpha(2*(n+1)+1)=   -3*para(2)*binomcoeffi(-2.5_sp,n)*para(3)**(2*(n+1))
    end do

    !Calculate the series-expansion of beta=-r^2*S/KCG(r)
    beta            =   0.0_sp
    partitions      =   0_i4

    !set all partitions to the "biggest" one
    do n=1_i4,size(beta)-1_i4
       partitions((n*(n+1_i4))/2_i4)   =   1_i4
    end do

    !Define the series expansion of the inner function
    innerfunc       =   0.0_sp
    innerfunc(1)    =   -para(2)

    do n=1_i4, min(Floor(para(10)/2.0_sp), ((stdsiz-1_i4)/2_i4))
       innerfunc(2_i4*n+1_i4)  =   -para(2)*binomcoeffi(-1.5_sp,n)*para(3)**(2_i4*n)
    end do

    beta(3)         =   1.0_sp

    !Define the series expansion of beta by the rule of Faa di Bruno
    do n=2_i4,min(2_i4*Floor(para(10)/2.0_sp), stdsiz-1_i4),2_i4
       do
          prod=1.0_sp
          do m=1_i4, n

             if ((partitions(((n-1)*n)/2+m) .ne. 0_i4)) then
                if  (ne(innerfunc(m+1), 0.0_sp)) then
                   prod    =&
                        prod*(innerfunc(m+1)**real(partitions(((n-1)*n)/2+m),sp))/real(factorial(partitions(((n-1)*n)/2+m)),sp)
                else
                   prod    =   0.0_sp
                   exit
                end if
             end if

          end do

          if (.not. (isnan(prod) .or. (abs(prod) > huge(1.0_sp)) .or. (abs(prod) < tiny(1.0_sp))))&
               beta(n+3)=beta(n+3)+prod

          if (partitions(((n-1)*n)/2+1)==n) exit

          partitions(((n-1)*n)/2+1:((n+1)*n)/2)    =   nextpart(partitions(((n-1)*n)/2+1:((n+1)*n)/2))
       end do

    end do

    beta     =   -(para(4)/para(1))*exp(-para(2))*beta


    do laps_n=1_i4, size(s)

       !Define the series expansion of the first fundamental solution g with the aid of alpha and beta by the frobenius-method
       g               =   0.0_sp
       g(1)            =   1.0_sp

       do n=1_i4,size(g)-1_i4
          do m=0_i4,n-1_i4
             g(n+1)  =   g(n+1)+(real(m,sp)*alpha(n+1-m)+s(laps_n)*beta(n+1-m))*g(m+1)
          end do
          g(n+1)      =   -1.0_sp/(real(n,sp)**2_i4)*g(n+1)
       end do

       !Define the series of the the first part of the second fundamental solution h with the aid of alpha, beta and g by the frob.-m.
       h               =   0.0_sp
       h(1)            =   1.0_sp

       do n=1_i4,size(h)-1_i4
          do m=0_i4,n-1_i4
             h(n+1)  =   h(n+1) + alpha(n+1-m)*(real(m,sp)*h(m+1)+g(m+1)) + s(laps_n)*beta(n+1-m)*h(m+1)
          end do
          h(n+1)      =   -1.0_sp/(real(n,sp)**2_i4)*(h(n+1)+2.0_sp*real(n,sp)*g(n+1))
       end do

       ! ------------------------------------------------------------------
       ! calculate the head
       ! ------------------------------------------------------------------

       coef                =   para(5)/(2.0_sp*PI*para(6)*para(7))

       sub3                =   (para(7)/para(9))*coef/s(laps_n)*bessk0(sqrt(para(8)**2_i4*para(4)*s(laps_n)/(para(6)*para(9))))
       !swap the nearfield conductivity with the reference-conductivity Kwell<->Kref in coef
       sub2                =   -coef/s(laps_n)

       sub1                =   (sub3-sub2*( poly(para(8),h)+log(para(8))*poly(para(8),g) ))/&
            poly(para(8),g)

       do rad_n=1_i4, size(r)
          lap_ext_theis3d_sp_d1(rad_n, laps_n)  =   (sub1+sub2*log(r(rad_n)))*poly(r(rad_n),g) + sub2*poly(r(rad_n),h)
       end do

    end do

  END FUNCTION lap_ext_theis3d_sp_d1

  FUNCTION ext_theis3d_dp_d1(rad, time, params, inits, BC, prop, grad, steh_n)

    use mo_laplace_inversion,    ONLY: NLInvSteh

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: time
    REAL(dp), DIMENSION(6),             INTENT(IN)  :: params               !params=(KG,var,corr,e,L,S)
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Ref,Qw)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=30]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]

    REAL(dp), DIMENSION(size(rad),size(time))       :: ext_theis3d_dp_d1

    REAL(dp)                                        :: pr, aniso, Kwell, Kefu, chi, Kref, C !C=prop/(e^(1/3)*corr)
    REAL(dp), DIMENSION(10)                         :: values
    INTEGER(i4)                                     :: stehlim, limit, n
    LOGICAL                                         :: BCkind

    LOGICAL, DIMENSION(size(rad), size(time))       :: mask

    ! aniso:            is the anisotropy function
    ! Kwell:            conductivity at the well
    ! Kefu:             far field conductivity for uniform flow
    ! Kref:             conductivity at the reference point (Ref=inits(1))
    ! chi:              abbreviation for ln(Kefu/Kwell)
    ! C:                abbreviation for the factor C=prop/(e^(1/3)*corr)
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! BCkind:           is the intern kind of the boundary condition
    ! pr:               is the intern proportionality factor
    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! mask:             is the mask for all non-zero time-points

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_dp))      stop 'Radius must be positiv.'
    if (any(time <    0.0_dp))      stop 'The time must be non-negativ.'
    if (le(inits(1) , 0.0_dp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_dp))      stop 'The geometric mean conductivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The correlation length must be positiv.'
    if (params(4) .LT. 0.0_dp)      stop 'The anisotropy ratio must be in [0,1].'
    if (params(4) .GT. 1.0_dp)      stop 'The anisotropy ratio must be in [0,1].'
    if (le(params(5), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

    if (present(BC)) then
       BCkind = BC
    else
       BCkind = .TRUE.
    endif

    if (present(prop)) then
       if (le(prop, 0.0_dp))        stop 'The proportionality factor must be positiv.'
       pr           =   prop
    else
       pr           =   1.6_dp
    endif

    if (present(steh_n)) then
       if (steh_n<=0_i4)            stop 'The Stehfest-boundary must be positiv.'
       stehlim      =   steh_n
    else
       stehlim      =   6_i4
    endif

    if (present(grad)) then
       if (grad<=3_i4)              stop 'The limit of the series expansion must be at least 4.'
       limit        =   grad
    else
       limit        =   30_i4
    endif

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    if (.NOT. BCkind)   pr=0.5_dp*pr

    C               =   pr/( (params(4)**(1.0_dp/3.0_dp))*params(3) )

    if (.not. all(rad<C**(-1.0_dp)))   stop 'Radius out of convergence-area.'

    if (eq(params(4), 1.0_dp)) then
       aniso       =   1.0_dp/3.0_dp
    else if (eq(params(4), 0.0_dp)) then
       aniso       =   0.0_dp
    else
       aniso       =   (params(4)/2.0_dp) * &
            (&
            ( (1.0_dp-params(4)**2_i4)**(-1.5_dp) )*atan((params(4)**2_i4 - 1.0_dp)**(-0.5_dp)) &
            - params(4)/(1.0_dp-params(4)**2_i4) &
            )
    endif

    Kefu            =   params(1)*exp(params(2)*(0.5_dp - aniso))

    if (BCkind) then
       Kwell       =   params(1)*exp(-params(2)/2.0_dp)
    else
       Kwell       =   params(1)*exp(params(2)/2.0_dp)
    endif

    chi             =   log(Kwell/Kefu)

    Kref            =   Kefu*exp( chi*( sqrt(1.0_dp+(C*inits(1))**2_i4) )**(-3.0_dp) )

    values          =   (/Kefu, chi, C, params(6), inits(2), params(5), Kwell, inits(1), Kref, real(limit,dp)/)

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    !set all heads to zero where time is zero
    do n=1_i4, size(rad)
       mask(n,:)       =   (time>0.0_dp)
    end do

    ext_theis3d_dp_d1   =   0.0_dp

    ext_theis3d_dp_d1   =&
         unpack(     pack(NLInvSteh(lap_ext_theis3d_dp_d1, values, rad, pack(time, (time>0.0_dp)), stehlim), .true.),&
         mask,&
         ext_theis3d_dp_d1)

    where( ext_theis3d_dp_d1 .gt. 0.0_dp )  ext_theis3d_dp_d1   =   0.0_dp

  END FUNCTION ext_theis3d_dp_d1

  FUNCTION lap_ext_theis3d_dp_d1(r, s, para)

    use mo_functions,       only: factorial
    use mo_nrutil,          only: poly
    use mo_combinatorics,   only: binomcoeffi, nextpart
    use mo_nr,              only: bessk0

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: r
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: s
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: para                     !para=(Kefu,chi,C,S,Qw,L,Kwell,Ref,Kref,limit)

    REAL(dp), DIMENSION(size(r), size(s))           :: lap_ext_theis3d_dp_d1

    !intern variables
    integer(i4), Parameter                          :: stdsiz=51_i4

    REAL(dp), &
         DIMENSION(max(nint(para(10),i4)+1,stdsiz))     :: g, h, alpha, beta, innerfunc   !expansions of the given functions

    integer(i4), &
         DIMENSION((size(beta)*(size(beta)-1_i4))/2_i4) :: partitions               !All partitionarrays in a row to calculate beta

    REAL(dp)                                        :: prod, sub1, sub2, sub3, coef

    INTEGER(i4)                                     :: n, m, rad_n, laps_n

    !Calculate the series-expansion of alpha=1-3chi*(C*r)^2*(1+(C*r)^2)^(5/2)
    alpha           =   0.0_dp
    alpha(1)        =   1.0_dp

    do n=0_i4, min(Floor(para(10)/2.0_dp - 1.0_dp), ((stdsiz-1_i4)/2_i4 - 1_i4))
       alpha(2*(n+1)+1)=   -3*para(2)*binomcoeffi(-2.5_dp,n)*para(3)**(2*(n+1))
    end do

    !Calculate the series-expansion of beta=-r^2*S/KCG(r)
    beta            =   0.0_dp
    partitions      =   0_i4

    !set all partitions to the "biggest" one
    do n=1_i4,size(beta)-1_i4
       partitions((n*(n+1_i4))/2_i4)   =   1_i4
    end do

    !Define the series expansion of the inner function
    innerfunc       =   0.0_dp
    innerfunc(1)    =   -para(2)

    do n=1_i4, min(Floor(para(10)/2.0_dp), ((stdsiz-1_i4)/2_i4))
       innerfunc(2_i4*n+1_i4)  =   -para(2)*binomcoeffi(-1.5_dp,n)*para(3)**(2_i4*n)
    end do

    beta(3)         =   1.0_dp

    !Define the series expansion of beta by the rule of Faa di Bruno
    do n=2_i4,min(2_i4*Floor(para(10)/2.0_dp), stdsiz-1_i4),2_i4
       do
          prod=1.0_dp
          do m=1_i4, n

             if ((partitions(((n-1)*n)/2+m) .ne. 0_i4)) then
                if  (ne(innerfunc(m+1), 0.0_dp)) then
                   prod    =&
                        prod*(innerfunc(m+1)**real(partitions(((n-1)*n)/2+m),dp))/real(factorial(partitions(((n-1)*n)/2+m)),dp)
                else
                   prod    =   0.0_dp
                   exit
                end if
             end if

          end do

          if (.not. (isnan(prod) .or. (abs(prod) > huge(1.0_dp)) .or. (abs(prod) < tiny(1.0_dp))))&
               beta(n+3)=beta(n+3)+prod

          if (partitions(((n-1)*n)/2+1)==n) exit

          partitions(((n-1)*n)/2+1:((n+1)*n)/2)    =   nextpart(partitions(((n-1)*n)/2+1:((n+1)*n)/2))
       end do

    end do

    beta     =   -(para(4)/para(1))*exp(-para(2))*beta

    !loop over s(:)
    do laps_n=1_i4, size(s)

       !Define the series expansion of the first fundamental solution g with the aid of alpha and beta by the frobenius-method
       g               =   0.0_dp
       g(1)            =   1.0_dp

       do n=1_i4,size(g)-1_i4
          do m=0_i4,n-1_i4
             g(n+1)  =   g(n+1)+(real(m,dp)*alpha(n+1-m)+s(laps_n)*beta(n+1-m))*g(m+1)
          end do
          g(n+1)      =   -1.0_dp/(real(n,dp)**2_i4)*g(n+1)
       end do

       !Define the series of the the first part of the second fundamental solution h with the aid of alpha, beta and g by the frob.-m.
       h               =   0.0_dp
       h(1)            =   1.0_dp

       do n=1_i4,size(h)-1_i4
          do m=0_i4,n-1_i4
             h(n+1)  =   h(n+1) + alpha(n+1-m)*(real(m,dp)*h(m+1)+g(m+1)) + s(laps_n)*beta(n+1-m)*h(m+1)
          end do
          h(n+1)      =   -1.0_dp/(real(n,dp)**2_i4)*(h(n+1)+2.0_dp*real(n,dp)*g(n+1))
       end do

       ! ------------------------------------------------------------------
       ! calculate the head
       ! ------------------------------------------------------------------

       coef                =   para(5)/(2.0_dp*PI_D*para(6)*para(7))

       sub3                =   (para(7)/para(9))*coef/s(laps_n)*bessk0(sqrt(para(8)**2_i4*para(4)*s(laps_n)/(para(6)*para(9))))
       !swap the nearfield conductivity with the reference-conductivity Kwell<->Kref in coef
       sub2                =   -coef/s(laps_n)

       sub1                =   (sub3-sub2*( poly(para(8),h)+log(para(8))*poly(para(8),g) ))/&
            poly(para(8),g)

       do rad_n=1_i4, size(r)
          lap_ext_theis3d_dp_d1(rad_n, laps_n)  =   (sub1+sub2*log(r(rad_n)))*poly(r(rad_n),g) + sub2*poly(r(rad_n),h)
       end do

    end do

  END FUNCTION lap_ext_theis3d_dp_d1

  FUNCTION ext_theis2d_sp(rad, time, params, inits, prop, grad, steh_n)

    use mo_laplace_inversion,    ONLY: NLInvSteh

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: rad
    REAL(sp),                           INTENT(IN)  :: time
    REAL(sp), DIMENSION(4),             INTENT(IN)  :: params               !params=(TG,var,corr,S)
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Ref,Qw)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality-factor [def=2]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=30]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]

    REAL(sp)                                        :: ext_theis2d_sp

    REAL(sp)                                        :: pr, Tref, C
    REAL(sp), DIMENSION(8)                          :: values
    INTEGER(i4)                                     :: stehlim, limit

    ! pr:               is the intern proportionality factor
    ! C:                abbreviation for the factor C=prop/corr
    ! Tref:             transmissivity at the reference point (Ref=inits(1))
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (le(rad      , 0.0_sp))      stop 'Radius must be positiv.'
    if (time .LT.     0.0_sp)       stop 'The time must be non-negativ.'
    if (le(inits(1) , 0.0_sp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_sp))      stop 'The transmissivity must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_sp))      stop 'The correlation length must be positiv.'

    if (present(prop)) then
       if (le(prop, 0.0_sp))       stop 'The proportionality factor must be positiv.'
       pr          =   prop
    else
       pr          =   2.0_sp
    endif

    if (present(steh_n)) then
       if (steh_n<=0_i4)            stop 'The Stehfest-boundary must be positiv.'
       stehlim      =   steh_n
    else
       stehlim      =   6_i4
    endif

    if (present(grad)) then
       if (grad<=3_i4)              stop 'The limit of the series expansion must be at least 4.'
       limit        =   grad
    else
       limit        =   30_i4
    endif

    C               =   pr/params(3)

    if (.not. (rad<C**(-1.0_sp)))   stop 'Radius out of convergence-area.'

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    Tref            =   params(1)*exp( -0.5_sp*(params(2)**2_i4)/(1.0_sp+(C*inits(1))**2_i4) )

    values          =   (/params(1), params(2)**2_i4, C, params(4), inits(1), inits(2), Tref, real(limit,sp)/)

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    if (time>0.0_sp) then
       ext_theis2d_sp  =   min(NLInvSteh(lap_ext_theis2d_sp, values, rad, time, stehlim),0.0_sp)
    else
       ext_theis2d_sp  =   0.0_sp
    end if

  END FUNCTION ext_theis2d_sp

  FUNCTION lap_ext_theis2d_sp(r, s, para)

    use mo_functions,       only: factorial
    use mo_nrutil,          only: poly
    use mo_combinatorics,   only: nextpart
    use mo_nr,              only: bessk0

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: r
    REAL(sp),                           INTENT(IN)  :: s
    REAL(sp), DIMENSION(:),             INTENT(IN)  :: para                     !para=(TG,var^2,C,S,Ref,Qw,Tref,limit)

    REAL(sp)                                        :: lap_ext_theis2d_sp

    !intern variables
    integer(i4), Parameter                          :: stdsiz=31_i4

    REAL(sp), &
         DIMENSION(max(nint(para(8),i4)+1,stdsiz))     :: g, h, alpha, beta, innerfunc   !expansions of the given functions

    integer(i4), &
         DIMENSION((size(beta)*(size(beta)-1_i4))/2_i4) :: partitions               !All partitionarrays in a row to calculate beta

    REAL(sp)                                        :: prod, sub1, sub2, sub3, coef, TH

    INTEGER(i4)                                     :: n,m

    !Calculate the series-expansion of alpha=1-3chi*(C*r)^2*(1+(C*r)^2)^(5/2)
    alpha           =   0.0_sp
    alpha(1)        =   1.0_sp

    do n=0_i4, min(Floor(para(8)/2.0_sp - 1.0_sp), ((stdsiz-1_i4)/2_i4 - 1_i4))
       alpha(2*(n+1)+1)=   ((-1.0_sp)**n)*para(2)*real(n+1_i4,sp)*para(3)**(2_i4*(n+1_i4))
    end do

    !Calculate the series-expansion of beta=-r^2*S/KCG(r)
    beta            =   0.0_sp
    partitions      =   0_i4

    !Set all partitions to the respectivly "biggest" one: (0,..,0,1)
    do n=1_i4,size(beta)-1_i4
       partitions((n*(n+1_i4))/2_i4)   =   1_i4
    end do

    !Define the series expansion of the inner function
    innerfunc       =   0.0_sp
    innerfunc(1)    =   0.5_sp*para(2)

    do n=1_i4, min(Floor(para(8)/2.0_sp), ((stdsiz-1_i4)/2_i4))
       innerfunc(2_i4*n+1_i4)  =   0.5_sp*((-1.0_sp)**n)*para(2)*para(3)**(2_i4*n)
    end do

    beta(3)         =   1.0_sp

    !Define the series expansion of beta by the rule of Faa di Bruno
    do n=2_i4,min(2_i4*Floor(para(8)/2.0_sp), stdsiz-1_i4),2_i4
       do
          prod=1.0_sp
          do m=1_i4, n

             if ((partitions(((n-1)*n)/2+m) .ne. 0_i4)) then
                if  (ne(innerfunc(m+1), 0.0_sp)) then
                   prod    =&
                        prod*(innerfunc(m+1)**real(partitions(((n-1)*n)/2+m),sp))/real(factorial(partitions(((n-1)*n)/2+m)),sp)
                else
                   prod    =   0.0_sp
                   exit
                end if
             end if

          end do

          if (.not. (isnan(prod) .or. (abs(prod) > huge(1.0_sp)) .or. (abs(prod) < tiny(1.0_sp))))&
               beta(n+3)=beta(n+3)+prod

          if (partitions(((n-1)*n)/2+1)==n) exit

          partitions(((n-1)*n)/2+1:((n+1)*n)/2)    =   nextpart(partitions(((n-1)*n)/2+1:((n+1)*n)/2))
       end do

    end do

    beta     =   -(para(4)/para(1))*exp(0.5_sp*para(2))*beta

    !Define the series expansion of the first fundamental solution g with the aid of alpha and beta by the frobenius-method
    g               =   0.0_sp
    g(1)            =   1.0_sp

    do n=1_i4,size(g)-1_i4
       do m=0_i4,n-1_i4
          g(n+1)  =   g(n+1)+(real(m,sp)*alpha(n+1-m)+s*beta(n+1-m))*g(m+1)
       end do
       g(n+1)      =   -1.0_sp/(real(n,sp)**2_i4)*g(n+1)
    end do

    !Define the series of the the first part of the second fundamental solution h with the aid of alpha, beta and g by the frobenius-m.
    h               =   0.0_sp
    h(1)            =   1.0_sp

    do n=1_i4,size(h)-1_i4
       do m=0_i4,n-1_i4
          h(n+1)  =   h(n+1) + alpha(n+1-m)*(real(m,sp)*h(m+1)+g(m+1)) + s*beta(n+1-m)*h(m+1)
       end do
       h(n+1)      =   -1.0_sp/(real(n,sp)**2_i4)*(h(n+1)+2.0_sp*real(n,sp)*g(n+1))
    end do

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    TH                  =   para(1)*exp(-para(2)*0.5_sp)                                !nearfield transmissivity

    coef                =   para(6)/(2.0_sp*PI*TH)

    sub3                =   (TH/para(7))*coef/s*bessk0(sqrt( para(5)**2_i4*para(4)*s/para(7) ))
    !swap the nearfield transmissivity with the reference-transmissivity TH<->Tref in coef
    sub2                =   -coef/s

    sub1                =   (sub3-sub2*( poly(para(5),h)+log(para(5))*poly(para(5),g) ))/&
         poly(para(5),g)

    lap_ext_theis2d_sp  =   (sub1+sub2*log(r))*poly(r,g) + sub2*poly(r,h)


  END FUNCTION lap_ext_theis2d_sp

  FUNCTION ext_theis2d_dp(rad, time, params, inits, prop, grad, steh_n)

    use mo_laplace_inversion,    ONLY: NLInvSteh

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: rad
    REAL(dp),                           INTENT(IN)  :: time
    REAL(dp), DIMENSION(4),             INTENT(IN)  :: params               !params=(TG,var,corr,S)
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Ref,Qw)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality-factor [def=2]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=30]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]

    REAL(dp)                                        :: ext_theis2d_dp

    REAL(dp)                                        :: pr, Tref, C !C=prop/corr
    REAL(dp), DIMENSION(8)                          :: values
    INTEGER(i4)                                     :: stehlim, limit

    ! pr:               is the intern proportionality factor
    ! Tref:             transmissivity at the reference point (Ref=inits(1))
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (le(rad      , 0.0_dp))      stop 'Radius must be positiv.'
    if (time .LT.     0.0_dp)       stop 'The time must be non-negativ.'
    if (le(inits(1) , 0.0_dp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_dp))      stop 'The transmissivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The correlation length must be positiv.'

    if (present(prop)) then
       if (le(prop, 0.0_dp))       stop 'The proportionality factor must be positiv.'
       pr          =   prop
    else
       pr          =   2.0_dp
    endif

    if (present(steh_n)) then
       if (steh_n<=0_i4)            stop 'The Stehfest-boundary must be positiv.'
       stehlim      =   steh_n
    else
       stehlim      =   6_i4
    endif

    if (present(grad)) then
       if (grad<=3_i4)              stop 'The limit of the series expansion must be at least 4.'
       limit        =   grad
    else
       limit        =   30_i4
    endif

    C               =   pr/params(3)

    if (.not. (rad<C**(-1.0_dp)))   stop 'Radius out of convergence-area.'

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    Tref            =   params(1)*exp( -0.5_dp*(params(2)**2_i4)/(1.0_dp+(C*inits(1))**2_i4) )

    values          =   (/params(1), params(2)**2_i4, C, params(4), inits(1), inits(2), Tref, real(limit,dp)/)

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    if (time>0.0_dp) then
       ext_theis2d_dp  =   min(NLInvSteh(lap_ext_theis2d_dp, values, rad, time, stehlim),0.0_dp)
    else
       ext_theis2d_dp  =   0.0_dp
    end if

  END FUNCTION ext_theis2d_dp

  FUNCTION lap_ext_theis2d_dp(r, s, para)

    use mo_functions,       only: factorial
    use mo_nrutil,          only: poly
    use mo_combinatorics,   only: nextpart
    use mo_nr,              only: bessk0

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: r
    REAL(dp),                           INTENT(IN)  :: s
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: para                     !para=(TG,var^2,C,S,Ref,Qw,Tref,limit)

    REAL(dp)                                        :: lap_ext_theis2d_dp

    !intern variables
    integer(i4), Parameter                          :: stdsiz=51_i4

    REAL(dp), &
         DIMENSION(max(nint(para(8),i4)+1,stdsiz))     :: g, h, alpha, beta, innerfunc   !expansions of the given functions

    integer(i4), &
         DIMENSION((size(beta)*(size(beta)-1_i4))/2_i4) :: partitions               !All partitionarrays in a row to calculate beta

    REAL(dp)                                        :: prod, sub1, sub2, sub3, coef, TH

    INTEGER(i4)                                     :: n,m

    !Calculate the series-expansion of alpha=1-3chi*(C*r)^2*(1+(C*r)^2)^(5/2)
    alpha           =   0.0_dp
    alpha(1)        =   1.0_dp

    do n=0_i4, min(Floor(para(8)/2.0_dp - 1.0_dp), ((stdsiz-1_i4)/2_i4 - 1_i4))
       alpha(2*(n+1)+1)=   ((-1.0_dp)**n)*para(2)*real(n+1_i4,dp)*para(3)**(2_i4*(n+1_i4))
    end do

    !Calculate the series-expansion of beta=-r^2*S/KCG(r)
    beta            =   0.0_dp
    partitions      =   0_i4

    !set all partitions to the "biggest" one
    do n=1_i4,size(beta)-1_i4
       partitions((n*(n+1_i4))/2_i4)   =   1_i4
    end do

    !Define the series expansion of the inner function
    innerfunc       =   0.0_dp
    innerfunc(1)    =   0.5_dp*para(2)

    do n=1_i4, min(Floor(para(8)/2.0_dp), ((stdsiz-1_i4)/2_i4))
       innerfunc(2_i4*n+1_i4)  =   0.5_dp*((-1.0_dp)**n)*para(2)*para(3)**(2_i4*n)
    end do

    beta(3)         =   1.0_dp

    !Define the series expansion of beta by the rule of Faa di Bruno
    do n=2_i4,min(2_i4*Floor(para(8)/2.0_dp), stdsiz-1_i4),2_i4
       do
          prod=1.0_dp
          do m=1_i4, n

             if ((partitions(((n-1)*n)/2+m) .ne. 0_i4)) then
                if  (ne(innerfunc(m+1), 0.0_dp)) then
                   prod    =&
                        prod*(innerfunc(m+1)**real(partitions(((n-1)*n)/2+m),dp))/real(factorial(partitions(((n-1)*n)/2+m)),dp)
                else
                   prod    =   0.0_dp
                   exit
                end if
             end if

          end do

          if (.not. (isnan(prod) .or. (abs(prod) > huge(1.0_dp)) .or. (abs(prod) < tiny(1.0_dp))))&
               beta(n+3)=beta(n+3)+prod

          if (partitions(((n-1)*n)/2+1)==n) exit

          partitions(((n-1)*n)/2+1:((n+1)*n)/2)    =   nextpart(partitions(((n-1)*n)/2+1:((n+1)*n)/2))
       end do

    end do

    beta     =   -(para(4)/para(1))*exp(0.5_dp*para(2))*beta

    !Define the series expansion of the first fundamental solution g with the aid of alpha and beta by the frobenius-method
    g               =   0.0_dp
    g(1)            =   1.0_dp

    do n=1_i4,size(g)-1_i4
       do m=0_i4,n-1_i4
          g(n+1)  =   g(n+1)+(real(m,dp)*alpha(n+1-m)+s*beta(n+1-m))*g(m+1)
       end do
       g(n+1)      =   -1.0_dp/(real(n,dp)**2_i4)*g(n+1)
    end do

    !Define the series of the the first part of the second fundamental solution h with the aid of alpha, beta and g by the frobenius-m.
    h               =   0.0_dp
    h(1)            =   1.0_dp

    do n=1_i4,size(h)-1_i4
       do m=0_i4,n-1_i4
          h(n+1)  =   h(n+1) + alpha(n+1-m)*(real(m,dp)*h(m+1)+g(m+1)) + s*beta(n+1-m)*h(m+1)
       end do
       h(n+1)      =   -1.0_dp/(real(n,dp)**2_i4)*(h(n+1)+2.0_dp*real(n,dp)*g(n+1))
    end do

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    TH                  =   para(1)*exp(-para(2)*0.5_dp)                                !nearfield transmissivity

    coef                =   para(6)/(2.0_dp*PI_D*TH)

    sub3                =   (TH/para(7))*coef/s*bessk0(sqrt( para(5)**2_i4*para(4)*s/para(7) ))
    !swap the nearfield transmissivity with the reference-transmissivity TH<->Tref in coef
    sub2                =   -coef/s

    sub1                =   (sub3-sub2*( poly(para(5),h)+log(para(5))*poly(para(5),g) ))/&
         poly(para(5),g)

    lap_ext_theis2d_dp  =   (sub1+sub2*log(r))*poly(r,g) + sub2*poly(r,h)


  END FUNCTION lap_ext_theis2d_dp


  FUNCTION ext_theis2d_sp_d1(rad, time, params, inits, prop, grad, steh_n)

    use mo_laplace_inversion,    ONLY: NLInvSteh

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(sp), DIMENSION(:),             INTENT(IN)  :: time
    REAL(sp), DIMENSION(4),             INTENT(IN)  :: params               !params=(TG,var,corr,S)
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Ref,Qw)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality-factor [def=2]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=30]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]

    REAL(sp), DIMENSION(size(rad),size(time))       :: ext_theis2d_sp_d1

    REAL(sp)                                        :: pr, Tref, C !C=prop/corr
    REAL(sp), DIMENSION(8)                          :: values
    INTEGER(i4)                                     :: stehlim, limit, n

    LOGICAL, DIMENSION(size(rad), size(time))       :: mask

    ! pr:               is the intern proportionality factor
    ! C:                abbreviation for the factor C=prop/corr
    ! Tref:             transmissivity at the reference point (Ref=inits(1))
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! mask:             is the mask for all non-zero time-points

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_sp))      stop 'Radius must be positiv.'
    if (any(time <    0.0_sp))      stop 'The time must be non-negativ.'
    if (le(inits(1) , 0.0_sp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_sp))      stop 'The transmissivity must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_sp))      stop 'The correlation length must be positiv.'

    if (present(prop)) then
       if (le(prop, 0.0_sp))       stop 'The proportionality factor must be positiv.'
       pr          =   prop
    else
       pr          =   2.0_sp
    endif

    if (present(steh_n)) then
       if (steh_n<=0_i4)            stop 'The Stehfest-boundary must be positiv.'
       stehlim      =   steh_n
    else
       stehlim      =   6_i4
    endif

    if (present(grad)) then
       if (grad<=3_i4)              stop 'The limit of the series expansion must be at least 4.'
       limit        =   grad
    else
       limit        =   30_i4
    endif

    C               =   pr/params(3)

    if (.not. all(rad<C**(-1.0_sp)))   stop 'Radius out of convergence-area.'

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    Tref            =   params(1)*exp( -0.5_sp*(params(2)**2_i4)/(1.0_sp+(C*inits(1))**2_i4) )

    values          =   (/params(1), params(2)**2_i4, C, params(4), inits(1), inits(2), Tref, real(limit,sp)/)

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do n=1_i4, size(rad)
       mask(n,:)       =   (time>0.0_sp)
    end do

    ext_theis2d_sp_d1   =   0.0_sp

    ext_theis2d_sp_d1   =&
         unpack( pack(NLInvSteh(lap_ext_theis2d_sp_d1, values, rad, pack(time, (time>0.0_sp)), stehlim), .true.),&
         mask,&
         ext_theis2d_sp_d1)

    where( ext_theis2d_sp_d1 .gt. 0.0_sp )  ext_theis2d_sp_d1   =   0.0_sp

  END FUNCTION ext_theis2d_sp_d1

  FUNCTION lap_ext_theis2d_sp_d1(r, s, para)

    use mo_functions,       only: factorial
    use mo_nrutil,          only: poly
    use mo_combinatorics,   only: nextpart
    use mo_nr,              only: bessk0

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: r
    REAL(sp), DIMENSION(:),             INTENT(IN)  :: s
    REAL(sp), DIMENSION(:),             INTENT(IN)  :: para                     !para=(TG,var^2,C,S,Ref,Qw,Tref,limit)

    REAL(sp), DIMENSION(size(r), size(s))           :: lap_ext_theis2d_sp_d1

    !intern variables
    integer(i4), Parameter                          :: stdsiz=31_i4

    REAL(sp), &
         DIMENSION(max(nint(para(8),i4)+1,stdsiz))     :: g, h, alpha, beta, innerfunc   !expansions of the given functions

    integer(i4), &
         DIMENSION((size(beta)*(size(beta)-1_i4))/2_i4) :: partitions               !All partitionarrays in a row to calculate beta

    REAL(sp)                                        :: prod, sub1, sub2, sub3, coef, TH

    INTEGER(i4)                                     :: n, m, rad_n, laps_n

    !Calculate the series-expansion of alpha=1-3chi*(C*r)^2*(1+(C*r)^2)^(5/2)
    alpha           =   0.0_sp
    alpha(1)        =   1.0_sp

    do n=0_i4, min(Floor(para(8)/2.0_sp - 1.0_sp), ((stdsiz-1_i4)/2_i4 - 1_i4))
       alpha(2*(n+1)+1)=   ((-1.0_sp)**n)*para(2)*real(n+1_i4,sp)*para(3)**(2_i4*(n+1_i4))
    end do

    !Calculate the series-expansion of beta=-r^2*S/KCG(r)
    beta            =   0.0_sp
    partitions      =   0_i4

    !set all partitions to the "biggest" one
    do n=1_i4,size(beta)-1_i4
       partitions((n*(n+1_i4))/2_i4)   =   1_i4
    end do

    !Define the series expansion of the inner function
    innerfunc       =   0.0_sp
    innerfunc(1)    =   0.5_sp*para(2)

    do n=1_i4, min(Floor(para(8)/2.0_sp), ((stdsiz-1_i4)/2_i4))
       innerfunc(2_i4*n+1_i4)  =   0.5_sp*((-1.0_sp)**n)*para(2)*para(3)**(2_i4*n)
    end do

    beta(3)         =   1.0_sp

    !Define the series expansion of beta by the rule of Faa di Bruno
    do n=2_i4,min(2_i4*Floor(para(8)/2.0_sp), stdsiz-1_i4),2_i4
       do
          prod=1.0_sp
          do m=1_i4, n

             if ((partitions(((n-1)*n)/2+m) .ne. 0_i4)) then
                if  (ne(innerfunc(m+1), 0.0_sp)) then
                   prod    =&
                        prod*(innerfunc(m+1)**real(partitions(((n-1)*n)/2+m),sp))/real(factorial(partitions(((n-1)*n)/2+m)),sp)
                else
                   prod    =   0.0_sp
                   exit
                end if
             end if

          end do

          if (.not. (isnan(prod) .or. (abs(prod) > huge(1.0_sp)) .or. (abs(prod) < tiny(1.0_sp))))&
               beta(n+3)=beta(n+3)+prod

          if (partitions(((n-1)*n)/2+1)==n) exit

          partitions(((n-1)*n)/2+1:((n+1)*n)/2)    =   nextpart(partitions(((n-1)*n)/2+1:((n+1)*n)/2))
       end do

    end do

    beta     =   -(para(4)/para(1))*exp(0.5_sp*para(2))*beta


    do laps_n=1_i4, size(s)

       !Define the series expansion of the first fundamental solution g with the aid of alpha and beta by the frobenius-method
       g               =   0.0_sp
       g(1)            =   1.0_sp

       do n=1_i4,size(g)-1_i4
          do m=0_i4,n-1_i4
             g(n+1)  =   g(n+1)+(real(m,sp)*alpha(n+1-m)+s(laps_n)*beta(n+1-m))*g(m+1)
          end do
          g(n+1)      =   -1.0_sp/(real(n,sp)**2_i4)*g(n+1)
       end do

       !Define the series of the the first part of the second fundamental solution h with the aid of alpha, beta and g by the frobenius-m.
       h               =   0.0_sp
       h(1)            =   1.0_sp

       do n=1_i4,size(h)-1_i4
          do m=0_i4,n-1_i4
             h(n+1)  =   h(n+1) + alpha(n+1-m)*(real(m,sp)*h(m+1)+g(m+1)) + s(laps_n)*beta(n+1-m)*h(m+1)
          end do
          h(n+1)      =   -1.0_sp/(real(n,sp)**2_i4)*(h(n+1)+2.0_sp*real(n,sp)*g(n+1))
       end do

       ! ------------------------------------------------------------------
       ! calculate the head
       ! ------------------------------------------------------------------

       TH                  =   para(1)*exp(-para(2)*0.5_sp)                                !nearfield transmissivity

       coef                =   para(6)/(2.0_sp*PI*TH)

       sub3                =   (TH/para(7))*coef/s(laps_n)*bessk0(sqrt( para(5)**2_i4*para(4)*s(laps_n)/para(7) ))
       !swap the nearfield transmissivity with the reference-transmissivity TH<->Tref in coef
       sub2                =   -coef/s(laps_n)

       sub1                =   (sub3-sub2*( poly(para(5),h)+log(para(5))*poly(para(5),g) ))/&
            poly(para(5),g)

       do rad_n=1_i4, size(r)
          lap_ext_theis2d_sp_d1(rad_n, laps_n)  =   (sub1+sub2*log(r(rad_n)))*poly(r(rad_n),g) + sub2*poly(r(rad_n),h)
       end do

    end do


  END FUNCTION lap_ext_theis2d_sp_d1

  FUNCTION ext_theis2d_dp_d1(rad, time, params, inits, prop, grad, steh_n)

    use mo_laplace_inversion,    ONLY: NLInvSteh

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: time
    REAL(dp), DIMENSION(4),             INTENT(IN)  :: params               !params=(TG,var,corr,S)
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Ref,Qw)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality-factor [def=2]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=30]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]

    REAL(dp), DIMENSION(size(rad),size(time))       :: ext_theis2d_dp_d1

    REAL(dp)                                        :: pr, Tref, C !C=prop/corr
    REAL(dp), DIMENSION(8)                          :: values
    INTEGER(i4)                                     :: stehlim, limit, n

    LOGICAL, DIMENSION(size(rad), size(time))       :: mask

    ! pr:               is the intern proportionality factor
    ! C:                abbreviation for the factor C=prop/corr
    ! Tref:             transmissivity at the reference point (Ref=inits(1))
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! mask:             is the mask for all non-zero time-points

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_dp))      stop 'Radius must be positiv.'
    if (any(time <    0.0_dp))      stop 'The time must be non-negativ.'
    if (le(inits(1) , 0.0_dp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_dp))      stop 'The transmissivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The correlation length must be positiv.'

    if (present(prop)) then
       if (le(prop, 0.0_dp))       stop 'The proportionality factor must be positiv.'
       pr          =   prop
    else
       pr          =   2.0_dp
    endif

    if (present(steh_n)) then
       if (steh_n<=0_i4)            stop 'The Stehfest-boundary must be positiv.'
       stehlim      =   steh_n
    else
       stehlim      =   6_i4
    endif

    if (present(grad)) then
       if (grad<=3_i4)              stop 'The limit of the series expansion must be at least 4.'
       limit        =   grad
    else
       limit        =   30_i4
    endif

    C               =   pr/params(3)

    if (.not. all(rad<C**(-1.0_dp)))   stop 'Radius out of convergence-area.'

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    Tref            =   params(1)*exp( -0.5_dp*(params(2)**2_i4)/(1.0_dp+(C*inits(1))**2_i4) )

    values          =   (/params(1), params(2)**2_i4, C, params(4), inits(1), inits(2), Tref, real(limit,dp)/)

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do n=1_i4, size(rad)
       mask(n,:)       =   (time>0.0_dp)
    end do

    ext_theis2d_dp_d1   =   0.0_dp

    ext_theis2d_dp_d1   =&
         unpack(&
         pack(NLInvSteh(lap_ext_theis2d_dp_d1, values, rad, pack(time, (time>0.0_dp)), stehlim), .true.),&
         mask,&
         ext_theis2d_dp_d1&
         )

    where( ext_theis2d_dp_d1 .gt. 0.0_dp )  ext_theis2d_dp_d1   =   0.0_dp

  END FUNCTION ext_theis2d_dp_d1

  FUNCTION lap_ext_theis2d_dp_d1(r, s, para)

    use mo_functions,       only: factorial
    use mo_nrutil,          only: poly
    use mo_combinatorics,   only: nextpart
    use mo_nr,              only: bessk0

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: r
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: s
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: para                     !para=(TG,var^2,C,S,Ref,Qw,Tref,limit)

    REAL(dp), DIMENSION(size(r), size(s))           :: lap_ext_theis2d_dp_d1

    !intern variables
    integer(i4), Parameter                          :: stdsiz=51_i4

    REAL(dp), &
         DIMENSION(max(nint(para(8),i4)+1,stdsiz))     :: g, h, alpha, beta, innerfunc   !expansions of the given functions

    integer(i4), &
         DIMENSION((size(beta)*(size(beta)-1_i4))/2_i4) :: partitions               !All partitionarrays in a row to calculate beta

    REAL(dp)                                        :: prod, sub1, sub2, sub3, coef, TH

    INTEGER(i4)                                     :: n, m, rad_n, laps_n

    !Calculate the series-expansion of alpha=1-3chi*(C*r)^2*(1+(C*r)^2)^(5/2)
    alpha           =   0.0_dp
    alpha(1)        =   1.0_dp

    do n=0_i4, min(Floor(para(8)/2.0_dp - 1.0_dp), ((stdsiz-1_i4)/2_i4 - 1_i4))
       alpha(2*(n+1)+1)=   ((-1.0_dp)**n)*para(2)*real(n+1_i4,dp)*para(3)**(2_i4*(n+1_i4))
    end do

    !Calculate the series-expansion of beta=-r^2*S/KCG(r)
    beta            =   0.0_dp
    partitions      =   0_i4

    !set all partitions to the "biggest" one
    do n=1_i4,size(beta)-1_i4
       partitions((n*(n+1_i4))/2_i4)   =   1_i4
    end do

    !Define the series expansion of the inner function
    innerfunc       =   0.0_dp
    innerfunc(1)    =   0.5_dp*para(2)

    do n=1_i4, min(Floor(para(8)/2.0_dp), ((stdsiz-1_i4)/2_i4))
       innerfunc(2_i4*n+1_i4)  =   0.5_dp*((-1.0_dp)**n)*para(2)*para(3)**(2_i4*n)
    end do

    beta(3)         =   1.0_dp

    !Define the series expansion of beta by the rule of Faa di Bruno
    do n=2_i4,min(2_i4*Floor(para(8)/2.0_dp), stdsiz-1_i4),2_i4
       do
          prod=1.0_dp
          do m=1_i4, n

             if ((partitions(((n-1)*n)/2+m) .ne. 0_i4)) then
                if  (ne(innerfunc(m+1), 0.0_dp)) then
                   prod    =&
                        prod*(innerfunc(m+1)**real(partitions(((n-1)*n)/2+m),dp))/real(factorial(partitions(((n-1)*n)/2+m)),dp)
                else
                   prod    =   0.0_dp
                   exit
                end if
             end if

          end do

          if (.not. (isnan(prod) .or. (abs(prod) > huge(1.0_dp)) .or. (abs(prod) < tiny(1.0_dp))))&
               beta(n+3)=beta(n+3)+prod

          if (partitions(((n-1)*n)/2+1)==n) exit

          partitions(((n-1)*n)/2+1:((n+1)*n)/2)    =   nextpart(partitions(((n-1)*n)/2+1:((n+1)*n)/2))
       end do

    end do

    beta     =   -(para(4)/para(1))*exp(0.5_dp*para(2))*beta


    do laps_n=1_i4, size(s)

       !Define the series expansion of the first fundamental solution g with the aid of alpha and beta by the frobenius-method
       g               =   0.0_dp
       g(1)            =   1.0_dp

       do n=1_i4,size(g)-1_i4
          do m=0_i4,n-1_i4
             g(n+1)  =   g(n+1)+(real(m,dp)*alpha(n+1-m)+s(laps_n)*beta(n+1-m))*g(m+1)
          end do
          g(n+1)      =   -1.0_dp/(real(n,dp)**2_i4)*g(n+1)
       end do

       !Define the series of the the first part of the second fundamental solution h with the aid of alpha, beta and g by the frob.-m.
       h               =   0.0_dp
       h(1)            =   1.0_dp

       do n=1_i4,size(h)-1_i4
          do m=0_i4,n-1_i4
             h(n+1)  =   h(n+1) + alpha(n+1-m)*(real(m,dp)*h(m+1)+g(m+1)) + s(laps_n)*beta(n+1-m)*h(m+1)
          end do
          h(n+1)      =   -1.0_dp/(real(n,dp)**2_i4)*(h(n+1)+2.0_dp*real(n,dp)*g(n+1))
       end do

       ! ------------------------------------------------------------------
       ! calculate the head
       ! ------------------------------------------------------------------

       TH                  =   para(1)*exp(-para(2)*0.5_dp)                                !nearfield transmissivity

       coef                =   para(6)/(2.0_dp*PI_D*TH)

       !swap the nearfield transmissivity with the reference-transmissivity TH<->Tref in coef
       sub3                =   (TH/para(7))*coef/s(laps_n)*bessk0(sqrt( para(5)**2_i4*para(4)*s(laps_n)/para(7) ))

       sub2                =   -coef/s(laps_n)

       sub1                =   (sub3-sub2*( poly(para(5),h)+log(para(5))*poly(para(5),g) ))/&
            poly(para(5),g)

       do rad_n=1_i4, size(r)
          lap_ext_theis2d_dp_d1(rad_n, laps_n)  =   (sub1+sub2*log(r(rad_n)))*poly(r(rad_n),g) + sub2*poly(r(rad_n),h)
       end do

    end do


  END FUNCTION lap_ext_theis2d_dp_d1

END MODULE mo_pumpingtests

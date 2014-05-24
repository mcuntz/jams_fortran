MODULE mo_groundwater

  ! This module serves different solutions for the groundwater-flow equation.

  ! Four kinds of solution are implemented:
  !     -Thiem's solution for steady-state and homogeneous aquifer                      :   thiem
  !     -Theis' solution for transient flow and homogenous aquifer                      :   theis
  !     -the extended thiem's solution for steady state and heterogeneous aquifer in 2D :   ext_thiem2d
  !     -the extended thiem's solution for steady state and heterogeneous aquifer and
  !      anisotropical impact in 3D                                                     :   ext_thiem3d
  !

  ! Written Sebastian Mueller, Apr 2014

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! If you use this routine in your work, you should cite the following reference
  ! Zech A, C L Schneider, and S Attinger (2012)
  !     The Extended Thiem’s solution: Including the impact of heterogeneity,
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
  USE mo_utils,     ONLY: eq, le

  IMPLICIT NONE

  PUBLIC :: thiem       ! thiem's solution for the steady-state and homogeneous aquifer
  PUBLIC :: theis       ! theis' solution for transient-flow and homogeneous aquifer
  PUBLIC :: ext_thiem2d ! the extended thiem's solution for the steady-state and heterogeneous aquifer in 2D
  PUBLIC :: ext_thiem3d ! the extended thiem's solution for the steady-state and heterogeneous aquifer in 3D


  ! ------------------------------------------------------------------

  !     NAME
  !         thiem

  !     PURPOSE
  !         Calculates the thiem solution h(r) for the steady-state and homogenous aquifers.

  !     CALLING SEQUENCE
  !         out = thiem(rad, params, inits)

  !     INTENT(IN)
  !         REAL(sp/dp)                 :: rad                          : input radius                  in m
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
  !         REAL(sp/dp)                 :: thiem        hydraulic head in m

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

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
  !         [1] Zech et al. 2012 The extended Thiem's solution -Including the Impact
  !             of heterogeneity, WRR

  !     HISTORY
  !         Written, Sebastian Mueller, Apr 2014

  INTERFACE thiem
     MODULE PROCEDURE thiem_sp, thiem_dp
  END INTERFACE thiem

  ! ------------------------------------------------------------------

  !     NAME
  !         theis

  !     PURPOSE
  !         Calculates the theis solution h(r) for transient flow and homogenous aquifers.

  !     CALLING SEQUENCE
  !         out = theis(rad, time, params, inits)

  !     INTENT(IN)
  !         REAL(sp/dp)                 :: rad          Radius          : input radius                  in m
  !         REAL(sp/dp)                 :: time         Time            : input time                    in s
  !         REAL(sp/dp), DIMENSION(3)   :: params       Parameters
  !                                                         params(1)   : specific storage coefficient  in m^3
  !                                                         params(2)   : hydraulic conductivity        in m/s
  !                                                         params(3)   : aquifer thickness             in m
  !         REAL(sp/dp), DIMENSION(1)   :: inits        Initial values
  !                                                         inits(1)    : pumping rate at the well      in m^3/s

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         REAL(sp/dp)                 :: theis        hydraulic head

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RESTRICTIONS
  !         1) Input values must be floating points.
  !         2) rad       >  0       [radius]
  !         2) time      >  0       [time]
  !         4) params(1) >  0       [storage coefficient]
  !         4) params(2) >  0       [conductivity]
  !         5) params(3) >  0       [aqifer thickness]

  !     EXAMPLE
  !         rad     = 0.5
  !         time    = 0.5
  !         inits   = (/ -0.001 /)
  !         params  = (/ 0.001 , 0.0001 , 1.0 /)

  !         h       = theis(rad, time, inits, params )

  !         -> see also example in test directory

  !     LITERATURE
  !         [1] Zech et al. 2012 The extended Thiem's solution -Including the Impact
  !             of heterogeneity, WRR

  !     HISTORY
  !         Written, Sebastian Mueller, Apr 2014

  INTERFACE theis
     MODULE PROCEDURE theis_sp, theis_dp
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
  !         REAL(sp/dp)                 :: rad          Radius          : input radius                  in m
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
  !         real(sp/dp)             :: ext_thiem2d      hydraulic head

  !     INTENT(IN), OPTIONAL
  !         REAL(sp)                :: prop             proportionality factor
  !                                                     default=2.0

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

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
     MODULE PROCEDURE ext_thiem2d_sp, ext_thiem2d_dp
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
  !         REAL(sp/dp)                 :: rad          Radius          : input radius                  in m
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
  !         REAL(sp/dp)                 :: ext_thiem3d  hydraulic head

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
     MODULE PROCEDURE ext_thiem3d_sp, ext_thiem3d_dp
  END INTERFACE ext_thiem3d

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

END MODULE mo_groundwater

MODULE mo_pumpingtests

  ! This module serves different solutions for the groundwater-flow equation,
  ! containing the classical Thiem and Theis solutions as well as
  ! the effektive well-flow solution for the effectiv Coarse-Graining Conductivity.

  ! These kinds of solution are implemented:
  !     -Thiem's solution for steady-state and homogeneous aquifer                      :   thiem
  !     -Theis' solution for transient flow and homogenous aquifer                      :   theis
  !     -the extended thiem's solution for steady state and heterogeneous aquifer in 2D :   ext_thiem2d
  !     -the extended thiem's solution for steady state and heterogeneous aquifer and
  !      anisotropical impact in 3D                                                     :   ext_thiem3d
  !     -the extended theis' solution for transient flow and heterogeneous aquifer and
  !      anisotropical impact in 2D                                                     :   ext_theis2d
  !     -the extended theis' solution for transient flow and heterogeneous aquifer and
  !      anisotropical impact in 3D                                                     :   ext_theis3d
  !

  ! Written Sebastian Mueller, Jun 2014 - Jun 2016

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2014-2016 Sebastian Mueller
  !
  ! Permission is hereby granted, free of charge, to any person obtaining a copy
  ! of this software and associated documentation files (the "Software"), to deal
  ! in the Software without restriction, including without limitation the rights
  ! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  ! copies of the Software, and to permit persons to whom the Software is
  ! furnished to do so, subject to the following conditions:
  !
  ! The above copyright notice and this permission notice shall be included in all
  ! copies or substantial portions of the Software.
  !
  ! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  ! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  ! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  ! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  ! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  ! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  ! SOFTWARE.

  ! If you use this routine in your work, you should cite the following reference
  ! Zech A, C L Schneider, and S Attinger (2012)
  !     The Extended Thiem's solution: Including the impact of heterogeneity,
  !     Water Resources Research 48, W10535, doi:10.1029/2012WR011852

  ! Note on Numerical Recipes License
  ! ---------------------------------
  ! Be aware that some code is under the Numerical Recipes License 3rd
  ! edition <http://www.nr.com/aboutNR3license.html>

  ! The Numerical Recipes Personal Single-User License lets you personally
  ! use Numerical Recipes code ("the code") on any number of computers,
  ! but only one computer at a time. You are not permitted to allow anyone
  ! else to access or use the code. You may, under this license, transfer
  ! precompiled, executable applications incorporating the code to other,
  ! unlicensed, persons, providing that (i) the application is
  ! noncommercial (i.e., does not involve the selling or licensing of the
  ! application for a fee), and (ii) the application was first developed,
  ! compiled, and successfully run by you, and (iii) the code is bound
  ! into the application in such a manner that it cannot be accessed as
  ! individual routines and cannot practicably be unbound and used in
  ! other programs. That is, under this license, your application user
  ! must not be able to use Numerical Recipes code as part of a program
  ! library or "mix and match" workbench.

  ! Businesses and organizations that purchase the disk or code download,
  ! and that thus acquire one or more Numerical Recipes Personal
  ! Single-User Licenses, may permanently assign those licenses, in the
  ! number acquired, to individual employees. Such an assignment must be
  ! made before the code is first used and, once made, it is irrevocable
  ! and can not be transferred.

  ! If you do not hold a Numerical Recipes License, this code is only for
  ! informational and educational purposes but cannot be used.

  USE mo_kind,      ONLY: i4, i8, sp, dp
  USE mo_constants, ONLY: PI_D, PI
  USE mo_utils,     ONLY: eq, le, ne, ge

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
  !         REAL(sp/dp), DIMENSION(1)   :: params       Parameters
  !                                                         params(1)   : hydraulic conductivity        in m/s
  !         REAL(sp/dp), DIMENSION(4)   :: inits        Initial values
  !                                                         inits(1)    : reference point radius R      in m
  !                                                         inits(2)    : reference head H              in m
  !                                                         inits(3)    : pumping rate at the well Q    in m^3/s
  !                                                         inits(4)    : aquifer thickness L           in m

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
  !         5) inits(4)  >  0       [aqifer thickness]

  !     EXAMPLE
  !         rad     = 0.5
  !         inits   = (/ 32.0 , 0.0 , -0.001 , 1.0 /)
  !         params  = (/ 0.0001 /)

  !         h       = thiem(rad, inits, params )

  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written, Sebastian Mueller, Apr 2014

  INTERFACE thiem
     MODULE PROCEDURE thiem_sp, thiem_dp, thiem_d1_sp, thiem_d1_dp
  END INTERFACE thiem
  PUBLIC :: thiem_sp, thiem_dp, thiem_d1_sp, thiem_d1_dp

  ! ------------------------------------------------------------------

  !     NAME
  !         theis

  !     PURPOSE
  !         Calculates the theis solution h(r) for transient flow and homogenous aquifers.
  !         Use either unstructured Radius-Time-points with "grid" or get a structured Radius-Time grid with "rad" and "time".

  !     CALLING SEQUENCE
  !         out = theis(rad, time, params, inits) or
  !         out = theis(grid, params, inits)

  !     INTENT(IN)
  !        [REAL(sp/dp),DIMENSION(:,2)  :: grid         Grid-Points     : input radius-time-points      in (m,s)]
  !         REAL(sp/dp),[DIMENSION(:)]  :: rad          Radius          : input radius                  in m
  !         REAL(sp/dp),[DIMENSION(:)]  :: time         Time            : input time                    in s
  !         REAL(sp/dp), DIMENSION(2)   :: params       Parameters
  !                                                         params(1)   : specific storage-coef. S      in m^3
  !                                                         params(2)   : hydraulic conductivity K      in m/s
  !         REAL(sp/dp), DIMENSION(2)   :: inits        Initial values
  !                                                         inits(1)    : pumping rate at the well Q    in m^3/s
  !                                                         inits(2)    : aquifer thickness L           in m

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
  !         6) inits(2)  >  0       [aqifer thickness]

  !     EXAMPLE
  !         rad     = 0.5
  !         time    = 0.5
  !         inits   = (/ -0.001 , 1.0 /)
  !         params  = (/ 0.001 , 0.0001 /)

  !         h       = theis(rad, time, inits, params )

  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written, Sebastian Mueller, Apr 2014

  INTERFACE theis
     MODULE PROCEDURE theis_sp, theis_dp, theis_d1_sp, theis_d1_dp, theis_grid_sp, theis_grid_dp
  END INTERFACE theis
  PUBLIC :: theis_sp, theis_dp, theis_d1_sp, theis_d1_dp, theis_grid_sp, theis_grid_dp

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
  !                                                         params(1)   : geom-mean transmissivity T_G  in m^2/s
  !                                                         params(2)   : variance of log-trans.
  !                                                         params(3)   : correlation length l          in m
  !         REAL(sp/dp), DIMENSION(3)   :: inits        Initial values
  !                                                         inits(1)    : reference point radius R      in m
  !                                                         inits(2)    : reference head H              in m
  !                                                         inits(3)    : pumping rate at the well Q    in m^3/s

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         REAL(sp/dp)                    :: prop         proportionality factor
  !                                                        default=2.0

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
  PUBLIC :: ext_thiem2d_sp, ext_thiem2d_dp, ext_thiem2d_d1_sp, ext_thiem2d_d1_dp

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
  !         REAL(sp/dp), DIMENSION(4)   :: params       Parameters
  !                                                         params(1)   : geom.-mean conductivity K_G  in m^2/s
  !                                                         params(2)   : variance of log-cond.
  !                                                         params(3)   : correlation length l          in m
  !                                                         params(4)   : anisotropy ratio e
  !         REAL(sp/dp), DIMENSION(4)   :: inits        Initial values
  !                                                         inits(1)    : reference point radius R      in m
  !                                                         inits(2)    : reference head H              in m
  !                                                         inits(3)    : pumping rate at the well Q    in m^3/s
  !                                                         inits(4)    : aquifer thickness L           in m

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
  !         8) inits(4)  >  0       [aquifer thickness]

  !     EXAMPLE
  !         rad     = 0.5
  !         inits   = (/ 32.0 , 0.0 , -0.001 , 1.0 /)
  !         params  = (/ 0.0001 , 1.0 , 10.0 , 1.0 /)

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
  PUBLIC :: ext_thiem3d_sp, ext_thiem3d_dp, ext_thiem3d_d1_sp, ext_thiem3d_d1_dp

  ! ------------------------------------------------------------------

  !     NAME
  !         ext_theis3d

  !     PURPOSE
  !         Calculates the extended Theis solution h(r) transient-state and heterogenous aquifers
  !         in 3 dimensions with anisotropical impact. Use either unstructured Radius-Time-points
  !         with "grid" or get a structured Radius-Time grid with "rad" and "time".

  !     CALLING SEQUENCE
  !         out = ext_theis3d(rad, time, params, inits, BC, prop, grad, steh_n, output)  or
  !         out = ext_theis3d(grid, params, inits, BC, prop, grad, steh_n, output)

  !     INTENT(IN)
  !        [REAL(sp/dp),DIMENSION(:,2)  :: grid         Grid-Points     : input radius-time-points      in (m,s)]
  !         REAL(sp/dp),[DIMENSION(:)]  :: rad          Radius          : input radius                  in m
  !         REAL(sp/dp),[DIMENSION(:)]  :: time         Time            : input radius                  in m
  !         REAL(sp/dp), DIMENSION(5)   :: params       Parameters
  !                                                         params(1)   : geometric mean conductivity   in m^2/s
  !                                                         params(2)   : variance of log-conductivity
  !                                                         params(3)   : correlation length            in m
  !                                                         params(4)   : specific storage coefficient  in m^3
  !                                                         params(5)   : anisotropy ratio
  !         REAL(sp/dp), DIMENSION(2)   :: inits        Initial values
  !                                                         inits(1)    : pumping rate at the well      in m^3/s
  !                                                         inits(2)    : aquifer thickness             in m

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
  !         LOGICAL                     :: output       gives information about the calculation during the computation
  !                                                     default=false
  !         INTEGER(i4)                 :: grad         limit of the series expansion for the analytical functions
  !                                                     default=40 (recommended)
  !         INTEGER(i4)                 :: steh_n       boundary-value for the stehfest-algorithm (numerical laplace-inversion)
  !                                                     default=6 (recommended, for fast computation)

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
  !         4) inits(1)  <  0       [pumping rate]
  !         5) inits(2)  >  0       [aquifer thickness]
  !         6) params(1) >  0       [conductivity]
  !         7) params(2) >  0       [variance]
  !         8) params(3) >  0       [correlation length]
  !         9) params(4) >  0       [storage coefficient]
  !        10) params(5) in [0,1]   [anisotropy ratio]

  !     EXAMPLE
  !         rad     = 0.5
  !         time    = 50.0
  !         inits   = (/-0.001, 1.0/)
  !         params  = (/0.0001, 1.0, 10.0, 0.001, 1.0/)

  !         h       = ext_theis3d(rad, time, inits, params )

  !         -> see also example in test directory

  !     LITERATURE
  !         [1] ?

  !     HISTORY
  !         Written, Sebastian Mueller, Oct 2015
  !         Modified, Matthias Cuntz, Sep 2016 - move g3D_dp outside of function

  INTERFACE ext_theis3d
     MODULE PROCEDURE ext_theis3d_sp, ext_theis3d_dp, ext_theis3d_sp_d1, ext_theis3d_dp_d1, ext_theis3d_grid_sp, ext_theis3d_grid_dp
  END INTERFACE ext_theis3d
  PUBLIC :: ext_theis3d_sp, ext_theis3d_dp, ext_theis3d_sp_d1, ext_theis3d_dp_d1, ext_theis3d_grid_sp, ext_theis3d_grid_dp

  ! ------------------------------------------------------------------

  !     NAME
  !         ext_theis2d

  !     PURPOSE
  !         Calculates the extended Theis solution h(r) transient-state and heterogenous aquifers in 2 dimensions.
  !         Use either unstructured Radius-Time-points with "grid" or get a structured Radius-Time grid with "rad" and "time".

  !     CALLING SEQUENCE
  !         out = ext_theis2d(rad, time, params, inits, prop, grad, steh_n, output) or
  !         out = ext_theis2d(grid, params, inits, prop, grad, steh_n, output)

  !     INTENT(IN)
  !        [REAL(sp/dp),DIMENSION(:,2)  :: grid         Grid-Points     : input radius-time-points      in (m,s)]
  !         REAL(sp/dp),[DIMENSION(:)]  :: rad          Radius          : input radius                  in m
  !         REAL(sp/dp),[DIMENSION(:)]  :: time         Time            : input radius                  in m
  !         REAL(sp/dp), DIMENSION(4)   :: params       Parameters
  !                                                         params(1)   : geometric mean transmissivity in m^2/s
  !                                                         params(2)   : variance of log-conductivity
  !                                                         params(3)   : correlation length            in m
  !                                                         params(4)   : specific storage coefficient  in m^3
  !         REAL(sp/dp), DIMENSION(1)   :: inits        Initial values
  !                                                         inits(1)    : pumping rate at the well      in m^3/s

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         REAL(sp/dp)                 :: prop         proportionality-factor
  !                                                     default=2.0
  !         LOGICAL                     :: output       gives information about the calculation during the computation
  !                                                     default=false
  !         INTEGER(i4)                 :: grad         limit of the series expansion for the analytical functions
  !                                                     default=40 (recommended)
  !         INTEGER(i4)                 :: steh_n       boundary-value for the stehfest-algorithm (numerical laplace-inversion)
  !                                                     default=6 (recommended, for fast computation)

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
  !         4) params(1) >  0       [transmissivity]
  !         5) params(2) >  0       [variance]
  !         6) params(3) >  0       [correlation length]
  !         7) params(4) >  0       [storage coefficient]

  !     EXAMPLE
  !         rad     = 0.5
  !         time    = 50.0
  !         inits   = (/ -0.001 /)
  !         params  = (/0.0001, 1.0, 10.0, 0.001/)

  !         h       = ext_theis2d(rad, time, inits, params )

  !         -> see also example in test directory

  !     LITERATURE
  !         [1] ?

  !     HISTORY
  !         Written, Sebastian Mueller, Oct 2015
  !         Modified, Matthias Cuntz, Sep 2016 - move g2D_dp outside of function

  INTERFACE ext_theis2d
     MODULE PROCEDURE ext_theis2d_sp, ext_theis2d_dp, ext_theis2d_sp_d1, ext_theis2d_dp_d1, ext_theis2d_grid_sp, ext_theis2d_grid_dp
  END INTERFACE ext_theis2d
  PUBLIC :: ext_theis2d_sp, ext_theis2d_dp, ext_theis2d_sp_d1, ext_theis2d_dp_d1, ext_theis2d_grid_sp, ext_theis2d_grid_dp

  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS


  !thiem solution

  FUNCTION thiem_sp(rad, params, inits)

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: rad
    REAL(sp), DIMENSION(1),             INTENT(IN)  :: params               !params=(K)
    REAL(sp), DIMENSION(4),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw,L)

    REAL(sp)                                        :: thiem_sp

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (le(rad, 0.0_sp))            stop 'Radius must be positiv!'
    if (le(inits(1), 0.0_sp))       stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_sp))      stop 'The conductivity must be positiv.'
    if (le(inits(4), 0.0_sp))       stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    thiem_sp        =   (-inits(3)/(2.0_sp*PI*inits(4)*params(1))) * log(rad/inits(1)) + inits(2)

  END FUNCTION thiem_sp

  FUNCTION thiem_dp(rad, params, inits)

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: rad
    REAL(dp), DIMENSION(1),             INTENT(IN)  :: params               !params=(K)
    REAL(dp), DIMENSION(4),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw,L)

    REAL(dp)                                        :: thiem_dp

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (le(rad, 0.0_dp))            stop 'Radius must be positiv!'
    if (le(inits(1), 0.0_dp))       stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_dp))      stop 'The conductivity must be positiv.'
    if (le(inits(4), 0.0_dp))       stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    thiem_dp        =   (-inits(3)/(2.0_dp*PI_D*inits(4)*params(1))) * log(rad/inits(1)) + inits(2)

  END FUNCTION thiem_dp

  FUNCTION thiem_d1_sp(rad, params, inits)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(sp), DIMENSION(1),             INTENT(IN)  :: params               !params=(K)
    REAL(sp), DIMENSION(4),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw,L)

    REAL(sp), DIMENSION(size(rad))                  :: thiem_d1_sp

    INTEGER(i4)                                     :: i

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_sp))      stop 'Radius must be positiv!'
    if (le(inits(1),  0.0_sp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_sp))      stop 'The conductivity must be positiv.'
    if (le(inits(4), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

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
    REAL(dp), DIMENSION(1),             INTENT(IN)  :: params               !params=(K)
    REAL(dp), DIMENSION(4),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw,L)

    REAL(dp), DIMENSION(size(rad))                  :: thiem_d1_dp

    INTEGER(i4)                                     :: i

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_dp))      stop 'Radius must be positiv!'
    if (le(inits(1),  0.0_dp))      stop 'Reference point must be at a positiv radius.'
    if (le(params(1), 0.0_dp))      stop 'The conductivity must be positiv.'
    if (le(inits(4), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do i=1_i4, size(rad)
       thiem_d1_dp(i)  =   thiem_dp(rad(i), params, inits)
    end do

  END FUNCTION thiem_d1_dp


  !theis solution

  FUNCTION theis_sp(rad, time, params, inits)

    USE mo_nr,                      ONLY: expint

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: rad, time
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: params               !params=(S,K)
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Qw,L)

    REAL(sp)                                        :: theis_sp

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (le(rad, 0.0_sp))            stop 'The radius must be positiv.'
    if (time .LT. 0.0_sp)           stop 'The time must be non-negativ.'
    if (le(params(1), 0.0_sp))      stop 'The specific storage coefficient must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The conductivity must be positiv.'
    if (le(inits(2), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    if (.not. eq(time, 0.0_sp)) then
       theis_sp    =   (inits(1)/(4.0_sp*PI*inits(2)*params(2))) &
            * expint(1_i4,((rad**2_i4)*params(1))/(4.0_sp*params(2)*time))
    else
       theis_sp    =   0.0_sp
    end if

  END FUNCTION theis_sp

  FUNCTION theis_dp(rad, time, params, inits)

    USE mo_nr,                      ONLY: expint

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: rad, time
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: params               !params=(S,K)
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Qw,L)

    REAL(dp)                                        :: theis_dp

    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (le(rad, 0.0_dp))            stop 'The radius must be positiv.'
    if (time .LT. 0.0_dp)           stop 'The time must be non-negativ.'
    if (le(params(1), 0.0_dp))      stop 'The specific storage coefficient must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The conductivity must be positiv.'
    if (le(inits(2), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    if (.not. eq(time, 0.0_dp)) then
       theis_dp    =   (inits(1)/(4.0_dp*PI_D*inits(2)*params(2))) &
            * expint(1_i4,((rad**2_i4)*params(1))/(4.0_dp*params(2)*time))
    else
       theis_dp    =   0.0_dp
    end if

  END FUNCTION theis_dp

  FUNCTION theis_d1_sp(rad, time, params, inits)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(sp), DIMENSION(:),             INTENT(IN)  :: time
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: params               !params=(S,K)
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Qw,L)

    REAL(sp), DIMENSION(size(rad),size(time))       :: theis_d1_sp

    INTEGER(i4)                                     :: i, j
    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_sp))      stop 'The radius must be positiv.'
    if (     any(time<0.0_sp))      stop 'The time must be non-negativ.'
    if (le(params(1), 0.0_sp))      stop 'The specific storage coefficient must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The conductivity must be positiv.'
    if (le(inits(2), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

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

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: time
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: params               !params=(S,K)
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Qw,L)

    REAL(dp), DIMENSION(size(rad),size(time))       :: theis_d1_dp

    INTEGER(i4)                                     :: i, j
    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_dp))      stop 'The radius must be positiv.'
    if (     any(time<0.0_dp))      stop 'The time must be non-negativ.'
    if (le(params(1), 0.0_dp))      stop 'The specific storage coefficient must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The conductivity must be positiv.'
    if (le(inits(2), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do i=1_i4, size(rad)
       do j=1_i4, size(time)
          theis_d1_dp(i,j)    =   theis_dp(rad(i), time(j), params, inits)
       end do
    end do

  END FUNCTION theis_d1_dp

  FUNCTION theis_grid_sp(grid, params, inits)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)  :: grid                 !Radius-Time points
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: params               !params=(S,K)
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Qw,L)

    REAL(sp), DIMENSION(size(grid(:,1)))            :: theis_grid_sp

    INTEGER(i4)                                     :: i
    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (size(grid(1,:)) .ne. 2_i4)      stop 'grid should contain 2 information: radius and time'
    if (.not. all(grid(:,1)>0.0_sp))    stop 'The radius must be positiv.'
    if (     any(grid(:,2)<0.0_sp))     stop 'The time must be non-negativ.'
    if (le(params(1), 0.0_sp))          stop 'The specific storage coefficient must be positiv.'
    if (le(params(2), 0.0_sp))          stop 'The conductivity must be positiv.'
    if (le(inits(2), 0.0_sp))           stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do i=1_i4, size(theis_grid_sp)
       theis_grid_sp(i)    =   theis_sp(grid(i,1), grid(i,2), params, inits)
    end do

  END FUNCTION theis_grid_sp

  FUNCTION theis_grid_dp(grid, params, inits)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)  :: grid                 !Radius-Time points
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: params               !params=(S,K)
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Qw,L)

    REAL(dp), DIMENSION(size(grid(:,1)))            :: theis_grid_dp

    INTEGER(i4)                                     :: i
    ! ------------------------------------------------------------------
    ! input tests
    ! ------------------------------------------------------------------

    if (size(grid(1,:)) .ne. 2_i4)      stop 'grid should contain 2 information: radius and time'
    if (.not. all(grid(:,1)>0.0_dp))    stop 'The radius must be positiv.'
    if (     any(grid(:,2)<0.0_dp))     stop 'The time must be non-negativ.'
    if (le(params(1), 0.0_dp))          stop 'The specific storage coefficient must be positiv.'
    if (le(params(2), 0.0_dp))          stop 'The conductivity must be positiv.'
    if (le(inits(2), 0.0_dp))           stop 'The aquifer thickness must be positiv.'

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do i=1_i4, size(theis_grid_dp)
       theis_grid_dp(i)    =   theis_dp(grid(i,1), grid(i,2), params, inits)
    end do

  END FUNCTION theis_grid_dp

  !extended thiems solution

  !2D

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

  !3D

  FUNCTION ext_thiem3d_sp(rad, params, inits, BC, prop)

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: rad
    REAL(sp), DIMENSION(4),             INTENT(IN)  :: params               !params=(KG,var,corr,e)
    REAL(sp), DIMENSION(4),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw,L)
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
    if (le(inits(4), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

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
            ( (1.0_sp-params(4)**2_i4)**(-1.5_sp) )*atan((1.0_sp/params(4)**2_i4 - 1.0_sp)**(0.5_sp)) &
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

    ext_thiem3d_sp  =   -(inits(3)/(2.0_sp*PI*inits(4)*Kefu))*&
         (&
         exp(-chi)*log(rad/inits(1)) + sinh(chi)*sub21 + (1.0_sp - cosh(chi))*sub22 &
         )&
         + inits(2)

  END FUNCTION ext_thiem3d_sp

  FUNCTION ext_thiem3d_dp(rad, params, inits, BC, prop)

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: rad
    REAL(dp), DIMENSION(4),             INTENT(IN)  :: params               !params=(KG,var,corr,e)
    REAL(dp), DIMENSION(4),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw,L)
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
    if (le(inits(4), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

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
            ( (1.0_dp-params(4)**2_i4)**(-1.5_dp) )*atan((1.0_dp/params(4)**2_i4 - 1.0_dp)**(0.5_dp)) &
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

    ext_thiem3d_dp  =   -(inits(3)/(2.0_dp*PI_D*inits(4)*Kefu))*&
         (&
         exp(-chi)*log(rad/inits(1)) + sinh(chi)*sub21 + (1.0_dp - cosh(chi))*sub22 &
         )&
         + inits(2)

  END FUNCTION ext_thiem3d_dp

  FUNCTION ext_thiem3d_d1_sp(rad, params, inits, BC, prop)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(sp), DIMENSION(4),             INTENT(IN)  :: params               !params=(KG,var,corr,e)
    REAL(sp), DIMENSION(4),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw,L)
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
    if (le(inits(4), 0.0_sp))      stop 'The aquifer thickness must be positiv.'

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
    REAL(dp), DIMENSION(4),             INTENT(IN)  :: params               !params=(KG,var,corr,e)
    REAL(dp), DIMENSION(4),             INTENT(IN)  :: inits                !inits=(Ref,href,Qw,L)
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
    if (le(inits(4), 0.0_dp))      stop 'The aquifer thickness must be positiv.'

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

  !3D

  FUNCTION ext_theis3d_sp(rad, time, params, inits, BC, prop, grad, steh_n, output)

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: rad
    REAL(sp),                           INTENT(IN)  :: time
    REAL(sp), DIMENSION(5),             INTENT(IN)  :: params               !params=(KG,sigma**2,corr,S,e)
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Qw,L)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=40]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: output               !enables output while computation

    REAL(sp)                                        :: ext_theis3d_sp

    REAL(sp)                                        :: pr
    REAL(sp), DIMENSION(1,1)                        :: temphead
    INTEGER(i4)                                     :: stehlim, limit
    LOGICAL                                         :: BCkind, output2


    ! aniso:            is the anisotropy-value depending on e
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


    if (.not. (rad>0.0_sp)   )      stop 'Radius must be positiv.'
    if (time <    0.0_sp     )      stop 'The time must be non-negativ.'
    if (le(inits(2) , 0.0_sp))      stop 'The aquifer thickness must be positiv.'
    if (ge(inits(1) , 0.0_sp))      stop 'The pumping-rate must be negativ.'
    if (le(params(1), 0.0_sp))      stop 'The geometric mean conductivity must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_sp))      stop 'The correlation length must be positiv.'
    if (le(params(4), 0.0_sp))      stop 'The storage must be positiv.'
    if (params(5) .LT. 0.0_sp)      stop 'The anisotropy ratio must be in [0,1].'
    if (params(5) .GT. 1.0_sp)      stop 'The anisotropy ratio must be in [0,1].'

    if (present(BC)) then
       BCkind       = BC
    else
       BCkind       = .TRUE.
    endif

    if (present(output)) then
       output2     = output
    else
       output2     = .false.
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
       limit        =   40_i4       !for limit>=44 you will get an error within factorial_i8 for floating-point exeption
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    temphead        =   ext_theis3d_sp_d1(&
         rad     =   (/rad/),&
         time    =   (/time/),&
         params  =   params,&
         inits   =   inits,&
         prop    =   pr,&
         grad    =   limit,&
         steh_n  =   stehlim,&
         BC      =   BCkind,&
         output  =   output2&
         )

    ext_theis3d_sp  =   temphead(1,1)

  END FUNCTION ext_theis3d_sp


  FUNCTION ext_theis3d_dp(rad, time, params, inits, BC, prop, grad, steh_n, output)

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: rad
    REAL(dp),                           INTENT(IN)  :: time
    REAL(dp), DIMENSION(5),             INTENT(IN)  :: params               !params=(KG,sigma**2,corr,S,e)
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Qw,L)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=40]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: output               !enables output while computation

    REAL(dp)                                        :: ext_theis3d_dp

    REAL(dp)                                        :: pr
    REAL(dp), DIMENSION(1,1)                        :: temphead
    INTEGER(i4)                                     :: stehlim, limit
    LOGICAL                                         :: BCkind, output2


    ! aniso:            is the anisotropy-value depending on e
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


    if (.not. (rad>0.0_dp)   )      stop 'Radius must be positiv.'
    if (    time <    0.0_dp )      stop 'The time must be non-negativ.'
    if (le(inits(2) , 0.0_dp))      stop 'The aquifer thickness must be positiv.'
    if (ge(inits(1) , 0.0_dp))      stop 'The pumping-rate must be negativ.'
    if (le(params(1), 0.0_dp))      stop 'The geometric mean conductivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The correlation length must be positiv.'
    if (le(params(4), 0.0_dp))      stop 'The storage must be positiv.'
    if (params(5) .LT. 0.0_dp)      stop 'The anisotropy ratio must be in [0,1].'
    if (params(5) .GT. 1.0_dp)      stop 'The anisotropy ratio must be in [0,1].'

    if (present(BC)) then
       BCkind       = BC
    else
       BCkind       = .TRUE.
    endif

    if (present(output)) then
       output2     = output
    else
       output2     = .false.
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
       limit        =   40_i4       !for limit>=44 you will get an error within factorial_i8 for floating-point exeption
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    temphead        =   ext_theis3d_dp_d1(&
         rad     =   (/rad/),&
         time    =   (/time/),&
         params  =   params,&
         inits   =   inits,&
         prop    =   pr,&
         grad    =   limit,&
         steh_n  =   stehlim,&
         BC      =   BCkind,&
         output  =   output2&
         )

    ext_theis3d_dp  =   temphead(1,1)

  END FUNCTION ext_theis3d_dp

  FUNCTION ext_theis3d_grid_sp(grid, params, inits, BC, prop, grad, steh_n, output)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)  :: grid
    REAL(sp), DIMENSION(5),             INTENT(IN)  :: params               !params=(KG,sigma**2,corr,S,e)
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Qw,L)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=40]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: output               !enables output while computation

    REAL(sp), DIMENSION(size(grid(:,1)),size(grid(:,1))) :: temphead
    REAL(sp), DIMENSION(size(grid(:,1)))            :: ext_theis3d_grid_sp

    REAL(sp)                                        :: pr
    INTEGER(i4)                                     :: stehlim, limit, i
    LOGICAL                                         :: BCkind, output2


    ! aniso:            is the anisotropy-value depending on e
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


    if (size(grid(1,:)) .ne. 2_i4)      stop 'grid should contain 2 information: radius and time'

    if (present(BC)) then
       BCkind       = BC
    else
       BCkind       = .TRUE.
    endif

    if (present(output)) then
       output2     = output
    else
       output2     = .false.
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
       limit        =   40_i4       !for limit>=44 you will get an error within factorial_i8 for floating-point exeption
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    temphead        =   real(ext_theis3d_dp_d1(&
         rad     =   real(grid(:,1),dp),&
         time    =   real(grid(:,2),dp),&
         params  =   real(params,dp),&
         inits   =   real(inits,dp),&
         prop    =   real(pr,dp),&
         grad    =   limit,&
         steh_n  =   stehlim,&
         BC      =   BCkind,&
         output  =   output2&
         ),sp)

    do i=1_i4, size(ext_theis3d_grid_sp)
       ext_theis3d_grid_sp(i)          =   temphead(i,i)
    end do

  END FUNCTION ext_theis3d_grid_sp

  FUNCTION ext_theis3d_grid_dp(grid, params, inits, BC, prop, grad, steh_n, output)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)  :: grid
    REAL(dp), DIMENSION(5),             INTENT(IN)  :: params               !params=(KG,sigma**2,corr,S,e)
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Qw,L)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=40]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: output               !enables output while computation

    REAL(dp), DIMENSION(size(grid(:,1)),size(grid(:,1))) :: temphead
    REAL(dp), DIMENSION(size(grid(:,1)))            :: ext_theis3d_grid_dp

    REAL(dp)                                        :: pr
    INTEGER(i4)                                     :: stehlim, limit, i
    LOGICAL                                         :: BCkind, output2


    ! aniso:            is the anisotropy-value depending on e
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


    if (size(grid(1,:)) .ne. 2_i4)      stop 'grid should contain 2 information: radius and time'

    if (present(BC)) then
       BCkind       = BC
    else
       BCkind       = .TRUE.
    endif

    if (present(output)) then
       output2     = output
    else
       output2     = .false.
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
       limit        =   40_i4       !for limit>=44 you will get an error within factorial_i8 for floating-point exeption
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    temphead        =   ext_theis3d_dp_d1(&
         rad     =   grid(:,1),&
         time    =   grid(:,2),&
         params  =   params,&
         inits   =   inits,&
         prop    =   pr,&
         grad    =   limit,&
         steh_n  =   stehlim,&
         BC      =   BCkind,&
         output  =   output2&
         )

    do i=1_i4, size(ext_theis3d_grid_dp)
       ext_theis3d_grid_dp(i)          =   temphead(i,i)
    end do

  END FUNCTION ext_theis3d_grid_dp


  FUNCTION ext_theis3d_sp_d1(rad, time, params, inits, BC, prop, grad, steh_n, output)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(sp), DIMENSION(:),             INTENT(IN)  :: time
    REAL(sp), DIMENSION(5),             INTENT(IN)  :: params               !params=(KG,sigma**2,corr,S,e)
    REAL(sp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Qw,L)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=40]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: output               !enables output while computation

    REAL(sp), DIMENSION(size(rad),size(time))       :: ext_theis3d_sp_d1

    REAL(sp)                                        :: pr
    INTEGER(i4)                                     :: stehlim, limit
    LOGICAL                                         :: BCkind, output2


    ! aniso:            is the anisotropy-value depending on e
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
    if (le(inits(2) , 0.0_sp))      stop 'The aquifer thickness must be positiv.'
    if (ge(inits(1) , 0.0_sp))      stop 'The pumping-rate must be negativ.'
    if (le(params(1), 0.0_sp))      stop 'The geometric mean conductivity must be positiv.'
    if (le(params(2), 0.0_sp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_sp))      stop 'The correlation length must be positiv.'
    if (le(params(4), 0.0_sp))      stop 'The storage must be positiv.'
    if (params(5) .LT. 0.0_sp)      stop 'The anisotropy ratio must be in [0,1].'
    if (params(5) .GT. 1.0_sp)      stop 'The anisotropy ratio must be in [0,1].'
    !    do n=2_i4, size(rad)
    !        if (rad(n) < rad(n-1))      stop 'Radial-points need to be sorted.'
    !    end do

    if (present(BC)) then
       BCkind       = BC
    else
       BCkind       = .TRUE.
    endif

    if (present(output)) then
       output2     = output
    else
       output2     = .false.
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
       limit        =   40_i4       !for limit>=44 you will get an error within factorial_i8 for floating-point exeption
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    ext_theis3d_sp_d1   =   real(ext_theis3d_dp_d1(&
         rad     =   real(rad,dp),&
         time    =   real(time,dp),&
         params  =   real(params,dp),&
         inits   =   real(inits,dp),&
         prop    =   real(pr,dp),&
         grad    =   limit,&
         steh_n  =   stehlim,&
         BC      =   BCkind,&
         output  =   output2&
         ),sp)

  END FUNCTION ext_theis3d_sp_d1


  FUNCTION ext_theis3d_dp_d1(rad, time, params, inits, BC, prop, grad, steh_n, output)

    use mo_laplace_inversion,    ONLY: NLInvSteh

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: time
    REAL(dp), DIMENSION(5),             INTENT(IN)  :: params               !params=(KG,sigma**2,corr,S,e)
    REAL(dp), DIMENSION(2),             INTENT(IN)  :: inits                !inits=(Qw,L)
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: BC                   !BC-kind (harmonic: true[default], arithmetic: false)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !prop-fac. for harm.case, arith.case autocorr. [def=1.6]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=40]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: output               !enables output while computation

    REAL(dp), DIMENSION(size(rad),size(time))       :: ext_theis3d_dp_d1

    REAL(dp)                                        :: pr, aniso, Kwell, Kefu, chi, output_real, C !C=prop/(e^(1/3)*corr)
    REAL(dp), DIMENSION(11)                         :: values
    INTEGER(i4)                                     :: stehlim, limit, n!, m
    LOGICAL                                         :: BCkind!, zero

    LOGICAL, DIMENSION(size(rad), size(time))       :: mask

    ! aniso:            is the anisotropy-value depending on e
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
    if (le(inits(2) , 0.0_dp))      stop 'The aquifer thickness must be positiv.'
    if (ge(inits(1) , 0.0_dp))      stop 'The pumping-rate must be negativ.'
    if (le(params(1), 0.0_dp))      stop 'The geometric mean conductivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The correlation length must be positiv.'
    if (le(params(4), 0.0_dp))      stop 'The storage must be positiv.'
    if (params(5) .LT. 0.0_dp)      stop 'The anisotropy ratio must be in [0,1].'
    if (params(5) .GT. 1.0_dp)      stop 'The anisotropy ratio must be in [0,1].'


    if (present(BC)) then
       BCkind       = BC
    else
       BCkind       = .TRUE.
    endif

    if (present(output)) then
       if (output) then
          output_real  = 1.0_dp
       else
          output_real  = 0.0_dp
       endif
    else
       output_real  = 0.0_dp
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
       limit        =   40_i4       !for limit=44 you will get an error within factorial_i8 for floating-point exeption
    endif

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    if (.NOT. BCkind)   pr = 0.5_dp*pr

    C               =   pr/( (params(5)**(1.0_dp/3.0_dp))*params(3) )

    if (eq(params(5), 1.0_dp)) then
       aniso       =   1.0_dp/3.0_dp
    else if (eq(params(5), 0.0_dp)) then
       aniso       =   0.0_dp
    else
       aniso       =   (params(5)/2.0_dp) * &
            (&
            ( (1.0_dp-params(5)**2_i4)**(-1.5_dp) )*atan((1.0_dp/params(5)**2_i4 - 1.0_dp)**(0.5_dp)) &
            - params(5)/(1.0_dp-params(5)**2_i4) &
            )
    endif

    Kefu            =   params(1)*exp(params(2)*(0.5_dp - aniso))

    if (BCkind) then
       Kwell       =   params(1)*exp(-params(2)/2.0_dp)
    else
       Kwell       =   params(1)*exp(params(2)/2.0_dp)
    endif

    chi             =   log(Kwell/Kefu)

    values          =   (/&
         params(1),          &! 1
         Kwell,              &! 2
         Kefu,               &! 3
         chi,                &! 4
         C,                  &! 5
         params(4),          &! 6
         inits(2),           &! 7
         inits(1),           &! 8
         real(limit,dp),     &! 9
         real(stehlim,dp),   &! 10
         output_real         &! 11
         /)

    !the optional output
    if (nint(output_real,i4)==1_i4) then
       write (*,*) ""
       write (*,*) "KG  : ",params(1)
       write (*,*) "sig2: ",params(2)
       write (*,*) "corr: ",params(3)
       write (*,*) "S   : ",params(4)
       write (*,*) "e   : ",params(5)
       write (*,*) "L   : ",inits(2)
       write (*,*) "Qw  : ",inits(1)
       write (*,*) "prop: ", pr
       if (BCkind) then
          write (*,*) "BC  : harmonic"
       else
          write (*,*) "BC  : arithmetic"
       endif
    end if

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    !set all heads to zero where time is zero
    do n=1_i4, size(rad)
       mask(n,:)       =   (time>0.0_dp)
    end do

    ext_theis3d_dp_d1   =   0.0_dp

    ext_theis3d_dp_d1   =&
         unpack(     pack(NLInvSteh(lap_ext_theis3d_dp_d1, values, C*rad, pack(time, (time>0.0_dp)), stehlim), .true.),&
         mask,&
         ext_theis3d_dp_d1)

    !correct numerical errors in the farfield through setting the head to 0.0
    where( ext_theis3d_dp_d1 .gt. 0.0_dp )          ext_theis3d_dp_d1   =   0.0_dp
    where( abs(ext_theis3d_dp_d1) .lt. 5e-4_dp )    ext_theis3d_dp_d1   =   0.0_dp

  END FUNCTION ext_theis3d_dp_d1


  FUNCTION lap_ext_theis3d_dp_d1(r, s, para)

    use mo_functions,       only: factorial
    use mo_nrutil,          only: poly
    use mo_combinatorics,   only: binomcoeffi, nextpart
    use mo_nr,              only: bessk0, bessk1
    use mo_interpol,        only: interpol
    use mo_ode_solver,      only: RBstiff
    use mo_linear_algebra,  only: solve_linear_equations

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: r
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: s
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: para              !para=(KG,Kwell,Kefu,chi,C,S,L,Qw,limit,stehlim,output)

    REAL(dp), DIMENSION(size(r), size(s))           :: lap_ext_theis3d_dp_d1

    !the standard limit for the series-expansion
    integer(i4), Parameter                          :: stdsiz=41_i4

    !expansions of the given function
    REAL(dp), &
         DIMENSION(max(nint(para(9),i4)+1,stdsiz))   :: dg, dh, g, h, alpha, beta, innerfunc

    !All partitionarrays in a row to calculate beta
    integer(i4), &
         DIMENSION((size(beta)*(size(beta)-1))/2)    :: partitions

    !functions for the rosenbrock-integration
    REAL(dp), DIMENSION(2)                          :: yi                   ! initial conditions
    REAL(dp), DIMENSION(:),     ALLOCATABLE         :: g1_xout, g2_xout     ! output time
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE         :: g1_out,  g2_out      ! output solutions
    REAL(dp), DIMENSION(1)                          :: g1temp,  g2temp      ! fundamental solutions evaluated at r(rad_n)
    REAL(dp), DIMENSION(5)                          :: g_para

    !variables for the rosenbrock-integration
    REAL(dp)                                        :: rinf,rmin,rmax,rlim  ! Boundary
    REAL(dp)                                        :: hstep                ! guessed stepsize
    REAL(dp)                                        :: hmin                 ! minimum allowed stepsize
    REAL(dp)                                        :: eps                  ! accurancy

    !intern variables
    REAL(dp)                                        :: prod, coef, T, t_temp, error, Kmid, ulim
    REAL(dp), DIMENSION(2,2)                        :: Mat
    REAL(dp), DIMENSION(2)                          :: LHS, RHS

    !loop-variables
    INTEGER(i4)                                     :: n, m, rad_n, laps_n, nn, mm

    !internt variable-names from para
    REAL(dp)                                        :: KG, Kwell, Kefu, chi, C, St, L, Qw, n_series
    INTEGER(i4)                                     :: n_steh
    LOGICAL                                         :: output

    !define the intern variables
    KG              =   para(1)
    Kwell           =   para(2)
    Kefu            =   para(3)
    chi             =   para(4)
    C               =   para(5)
    St              =   para(6)
    L               =   para(7)
    Qw              =   para(8)
    n_series        =   para(9)                 !is real
    n_steh          =   nint(para(10),i4)
    output          =   nint(para(11),i4)==1_i4


    !Calculate the series-expansion of alpha=1-3chi*(C*r)^2*(1+(C*r)^2)^(5/2)
    alpha           =   0.0_dp
    alpha(1)        =   1.0_dp

    do n=0_i4, min(Floor(n_series/2.0_dp - 1.0_dp), ((stdsiz-1_i4)/2_i4 - 1_i4))
       alpha(2*(n+1)+1)     =   -3.0_dp*chi*binomcoeffi(-2.5_dp,n)
    end do

    !Calculate the series-expansion of beta=-r^2*St/C/KCG(r)
    beta            =   0.0_dp
    partitions      =   0_i4

    !set all partitions to the "biggest" one
    do n=1_i4,size(beta)-1_i4
       partitions((n*(n+1_i4))/2_i4)   =   1_i4
    end do

    !Define the series expansion of the inner function in beta
    innerfunc       =   0.0_dp
    innerfunc(1)    =   -chi

    do n=1_i4, min(Floor(n_series/2.0_dp), ((stdsiz-1_i4)/2_i4))
       innerfunc(2_i4*n+1_i4)  =   -chi*binomcoeffi(-1.5_dp,n)
    end do

    beta(3)         =   1.0_dp

    !Define the series expansion of beta by the rule of Faa di Bruno
    do n=2_i4,min(2_i4*Floor(n_series/2.0_dp), stdsiz-1_i4)-2_i4,2_i4
       do
          prod=1.0_dp
          do m=1_i4, n
             if ((partitions(((n-1)*n)/2+m) .ne. 0_i4)) then
                if  (ne(innerfunc(m+1), 0.0_dp)) then
                   prod    =&
                        prod*(innerfunc(m+1)**partitions(((n-1)*n)/2+m))/real(factorial(int(partitions(((n-1)*n)/2+m),i8)),dp)
                else
                   prod    =   0.0_dp
                   exit
                end if
             end if
          end do

          beta(n+3)           =   beta(n+3)+prod

          if (partitions(((n-1)*n)/2+1)==n) exit

          partitions(((n-1)*n)/2+1:((n+1)*n)/2)    =   nextpart(partitions(((n-1)*n)/2+1:((n+1)*n)/2))

       end do
    end do

    beta            =   -(St/Kefu)*exp(-chi)/(C**2_i4)*beta

    !Define the boundaries for the different solutions (just reference-values)
    rinf            =   0.25_dp
    ulim            =   0.25_dp

    !define the upper limit for the solution such that the relativ error between KCG(r) and Kefu is less then 1% (error)
    error           =   0.01_dp

    if (chi>0.0_dp) then
       rlim        =   sqrt( (chi/log(1.0_dp+error))**(2.0_dp/3.0_dp) - 1.0_dp )
    elseif (chi<0.0_dp) then
       rlim        =   sqrt( (chi/log(1.0_dp-error))**(2.0_dp/3.0_dp) - 1.0_dp )
    else
       rlim        =   0.5_dp
    endif

    !define a weighted mean between Kwell and Kefu
    Kmid            =   (Kefu+2.0_dp*Kwell)/3.0_dp

    !parameters for the rosenbrock-integration with adaptive stepsize
    hstep           =   1e-12_dp                           ! incremental time step
    hmin            =   0.0_dp                             ! minimum allowed stepsize
    eps             =   1e-12_dp                           ! accuracy

    !optional output
    if (output) then
       write (*,*) ""
       write (*,*) "ulim: ", ulim
       write (*,*) "Rinf: ", rinf/C,   " Rinf*C: ", rinf
       write (*,*) "Rlim: ", rlim/C,   " Rlim*C: ", rlim
       write (*,*) "K_err:", error*100.0_dp,"percent"
       write (*,*) "C:    ", C
       write (*,*) "chi:  ", chi
       write (*,*) "KG:   ", KG
       write (*,*) "Kefu: ", Kefu     ," Kwell:  ", Kwell
       write (*,*) "Klim: ", Kefu*exp(chi/(sqrt((1.0_dp+(rlim)**2_i4))**3_i4))
       write (*,*) "Kmid: ", Kmid
       write (*,*) ""
    endif

    !loop over s(:)
    do nn=1_i4, size(s)/(n_steh)

       t_temp          =   1.0_dp/(s(nn)/log(2.0_dp))

       !Define the boundaries by comparisson to theis with the value U=r^2*St/(4*L*Kmid*t)<0.25
       rmax            =   min(sqrt(ulim*4.0_dp*L*t_temp*Kmid/St)*C, rlim)
       !        rmax            =   sqrt(ulim*4.0_dp*L*t_temp*Kmid/St)*C
       !        rmax            =   rlim
       rmin            =   min(rmax/2.0_dp, rinf)

       if (output) then
          write (*,*) "t_temp: ",t_temp," rmin: ",rmin/C," rmax: ",rmax/C, " eps: ",eps, " K_CG(rmax): ",&
               Kefu*exp(chi/(sqrt((1.0_dp+(rmax)**2_i4))**3_i4))
       endif

       do mm=0_i4, n_steh-1_i4
          laps_n = mm*size(s)/(n_steh) + nn

          !Define the series expansion of the first fundamental solution g with the aid of alpha and beta by the frobenius-method
          g               =   0.0_dp
          g(1)            =   1.0_dp
          dg              =   0.0_dp

          do n=1_i4,size(g)-1_i4
             do m=0_i4,n-1_i4
                g(n+1)  =   g(n+1)+(real(m,dp)*alpha(n+1-m) + s(laps_n)*beta(n+1-m))*g(m+1)
             end do
             g(n+1)      =   -1.0_dp/(real(n,dp)**2_i4)*g(n+1)
             dg(n)       =   real(n,dp)*g(n+1)
          end do


          !Define the series of the the first part of the second fundamental solution h with the aid of alph, beta & g by froben.
          h               =   0.0_dp
          h(1)            =   0.0_dp  !1.0_dp
          dh              =   0.0_dp

          do n=1_i4,size(h)-1_i4
             do m=0_i4,n-1_i4
                h(n+1)  =   h(n+1) + alpha(n+1-m)*(real(m,dp)*h(m+1)+g(m+1)) + s(laps_n)*beta(n+1-m)*h(m+1)
             end do
             h(n+1)      =   -1.0_dp/(real(n,dp)**2_i4)*(h(n+1)+2.0_dp*real(n,dp)*g(n+1))
             dh(n)       =   real(n,dp)*h(n+1)
          end do

          !Calculate 2 fundamental solutions in [rmin,rmax] with Rosenbrock-integration
          g_para          =   (/chi, s(laps_n), St, Kefu, C/)       !g_para=(chi, s, St, Kefu, C)

          coef            =   -Qw/(2.0_dp*PI_D*L*Kwell*s(laps_n))

          ! 1. solution
          yi              =   (/ poly(rmin,g),&
               poly(rmin,dg) /)

          call RBstiff( yi, rmin, rmax, hstep, g3D_dp, dg3D_dp, g_para, g1_xout, g1_out, hmin=hmin, eps=eps)

          ! 2. solution
          yi              =   (/ coef*(poly(rmin,h)  + log(rmin)*poly(rmin,g)),&
               coef*(poly(rmin,dh) + log(rmin)*poly(rmin,dg) + poly(rmin,g)/rmin) /)

          call RBstiff( yi, rmin, rmax, hstep, g3D_dp, dg3D_dp, g_para, g2_xout, g2_out, hmin=hmin, eps=eps)

          !clue the solution together with theis in the farfield (g\in C^1)
          Mat             =   0.0_dp
          RHS             =   0.0_dp

          !        T               =   sqrt(s(laps_n)*St/(Kefu*exp(chi/(sqrt((1.0_dp+(rmax)**2_i4))**3_i4))))/C
          T               =   sqrt(s(laps_n)*St/Kefu)/C

          Mat(1,1:1)      =   interpol(g1_out(:,1)  , g1_xout, (/rmax/))
          Mat(1,2)        =   -bessk0(T*rmax)
          Mat(2,1:1)      =   interpol(g1_out(:,2)  , g1_xout, (/rmax/))
          Mat(2,2)        =   T*bessk1(T*rmax)

          RHS(1:1)        =   -interpol(g2_out(:,1) , g2_xout, (/rmax/))
          RHS(2:2)        =   -interpol(g2_out(:,2) , g2_xout, (/rmax/))

          LHS             =   solve_linear_equations(Mat,RHS)

          !calculate the head in laplace-space
          do rad_n=1_i4, size(r)
             if  (r(rad_n)<rmin) then
                lap_ext_theis3d_dp_d1(rad_n, laps_n)  =   (coef*log(r(rad_n)) + LHS(1))*poly(r(rad_n),g) + coef*poly(r(rad_n),h)
             else if (r(rad_n)>rmax) then
                lap_ext_theis3d_dp_d1(rad_n, laps_n)  =   LHS(2)*bessk0(T*r(rad_n))
             else
                g1temp                                =   interpol(g1_out(:,1) , g1_xout, (/r(rad_n)/))
                g2temp                                =   interpol(g2_out(:,1) , g2_xout, (/r(rad_n)/))
                lap_ext_theis3d_dp_d1(rad_n, laps_n)  =   LHS(1)*g1temp(1) + g2temp(1)
             end if
          end do

       end do  !m
    end do  !n

  END FUNCTION lap_ext_theis3d_dp_d1

  subroutine g3D_dp( x, y, para1, dydx )   ! returns derivatives dydx at x for the first fundamental solution g

    implicit none

    real(dp), intent(in) :: x                       ! time
    real(dp), dimension(:), intent(in) :: y         ! unknowns of the equations
    real(dp), dimension(:), intent(in) :: para1      ! para1=(chi, s, S, Kefu, C)
    real(dp), dimension(:), intent(out) :: dydx     ! derivatives of y

    dydx(1)     =   y(2)
    dydx(2)     =   (3.0_dp*para1(1)*x/(sqrt(1+x**2_i4)**5_i4) - 1.0_dp/x)                               *   y(2) &
         + (exp(-para1(1)/(sqrt(1+x**2_i4)**3_i4))*para1(2)*para1(3)/para1(4)/((para1(5))**2_i4))   *   y(1)

  end subroutine g3D_dp

  subroutine dg3D_dp( x, y, para2, dfdx, dfdy )   ! returns derivatives dfdx and dfdy at x for the first fundamental solution g

    implicit none

    real(dp), intent(in) :: x                       ! time
    real(dp), dimension(:), intent(in) :: y         ! unknowns of the equations
    real(dp), dimension(:), intent(in) :: para2      ! para2=(chi, s, S, Kefu, C)
    real(dp), dimension(:), intent(out) :: dfdx     ! derivatives of f for x
    real(dp), dimension(:,:), intent(out) :: dfdy   ! derivatives of f for y

    dfdx(1)     =   0.0_dp
    dfdx(2)     =   ( 3.0_dp*para2(1)*(1.0_dp - 4.0_dp*x**2_i4)/(sqrt(1+x**2_i4)**7_i4) + 1.0_dp/(x**2_i4) )     *   y(2) &
         + ( exp(-para2(1)/(sqrt(1+x**2_i4)**3_i4))*para2(2)*para2(3)/para2(4)/((para2(5))**2_i4)*&
         3.0_dp*para2(1)*x/(sqrt(1+x**2_i4)**5_i4) )                                                *   y(1)

    dfdy(1,1)   =   0.0_dp
    dfdy(1,2)   =   1.0_dp
    dfdy(2,1)   =   (exp(-para2(1)/(sqrt(1+x**2_i4)**3_i4))*para2(2)*para2(3)/para2(4)/((para2(5))**2_i4))
    dfdy(2,2)   =   (3.0_dp*para2(1)*x/(sqrt(1+x**2_i4)**5_i4) - 1.0_dp/x)

  end subroutine dg3D_dp

  !2D

  FUNCTION ext_theis2d_sp(rad, time, params, inits, prop, grad, steh_n, output)

    IMPLICIT NONE

    REAL(sp),                           INTENT(IN)  :: rad
    REAL(sp),                           INTENT(IN)  :: time
    REAL(sp), DIMENSION(4),             INTENT(IN)  :: params               !params=(TG,sigma**2,corr,S)
    REAL(sp), DIMENSION(1),             INTENT(IN)  :: inits                !inits=(Qw)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality-factor [def=2]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=50]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: output               !enables output while computation

    REAL(sp)                                        :: ext_theis2d_sp

    REAL(sp)                                        :: pr
    REAL(sp), DIMENSION(1,1)                        :: temphead
    INTEGER(i4)                                     :: stehlim, limit

    LOGICAL                                         :: output2

    ! pr:               is the intern proportionality factor
    ! C:                abbreviation for the factor C=prop/corr
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! mask:             is the mask for all non-zero time-points

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (present(output)) then
       output2      = output
    else
       output2      = .false.
    endif

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
       limit        =   40_i4
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    temphead            =   real(ext_theis2d_dp_d1(&
         rad     =   (/real(rad,dp)/),&
         time    =   (/real(time,dp)/),&
         params  =   real(params,dp),&
         inits   =   real(inits,dp),&
         prop    =   real(pr,dp),&
         grad    =   limit,&
         steh_n  =   stehlim,&
         output  =   output2&
         ),sp)

    ext_theis2d_sp      =   temphead(1,1)

  END FUNCTION ext_theis2d_sp

  FUNCTION ext_theis2d_dp(rad, time, params, inits, prop, grad, steh_n, output)

    IMPLICIT NONE

    REAL(dp),                           INTENT(IN)  :: rad
    REAL(dp),                           INTENT(IN)  :: time
    REAL(dp), DIMENSION(4),             INTENT(IN)  :: params               !params=(TG,sigma**2,corr,S)
    REAL(dp), DIMENSION(1),             INTENT(IN)  :: inits                !inits=(Qw)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality-factor [def=2]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=50]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: output               !enables output while computation

    REAL(dp)                                        :: ext_theis2d_dp

    REAL(dp)                                        :: pr
    REAL(dp), DIMENSION(1,1)                        :: temphead
    INTEGER(i4)                                     :: stehlim, limit!, n, m

    LOGICAL                                         :: output2

    ! pr:               is the intern proportionality factor
    ! C:                abbreviation for the factor C=prop/corr
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! mask:             is the mask for all non-zero time-points

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (present(output)) then
       output2      = output
    else
       output2      = .false.
    endif

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
       limit        =   40_i4
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    temphead            =   ext_theis2d_dp_d1(&
         rad     =   (/rad/),&
         time    =   (/time/),&
         params  =   params,&
         inits   =   inits,&
         prop    =   pr,&
         grad    =   limit,&
         steh_n  =   stehlim,&
         output  =   output2&
         )

    ext_theis2d_dp      =   temphead(1,1)

  END FUNCTION ext_theis2d_dp

  FUNCTION ext_theis2d_grid_sp(grid, params, inits, prop, grad, steh_n, output)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:,:),           INTENT(IN)  :: grid
    REAL(sp), DIMENSION(4),             INTENT(IN)  :: params               !params=(TG,sigma**2,corr,S)
    REAL(sp), DIMENSION(1),             INTENT(IN)  :: inits                !inits=(Qw)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality-factor [def=2]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=50]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: output               !enables output while computation

    REAL(sp), DIMENSION(size(grid(:,1)),size(grid(:,1))) :: temphead
    REAL(sp), DIMENSION(size(grid(:,1)))            :: ext_theis2d_grid_sp

    REAL(sp)                                        :: pr
    INTEGER(i4)                                     :: stehlim, limit, i

    LOGICAL                                         :: output2

    ! pr:               is the intern proportionality factor
    ! C:                abbreviation for the factor C=prop/corr
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! mask:             is the mask for all non-zero time-points

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (size(grid(1,:)) .ne. 2_i4)      stop 'grid should contain 2 information: radius and time'

    if (present(output)) then
       output2      = output
    else
       output2      = .false.
    endif

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
       limit        =   40_i4
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    temphead            =   real(ext_theis2d_dp_d1(&
         rad     =   real(grid(:,1),dp),&
         time    =   real(grid(:,2),dp),&
         params  =   real(params,dp),&
         inits   =   real(inits,dp),&
         prop    =   real(pr,dp),&
         grad    =   limit,&
         steh_n  =   stehlim,&
         output  =   output2&
         ),sp)

    do i=1_i4, size(ext_theis2d_grid_sp)
       ext_theis2d_grid_sp(i)          =   temphead(i,i)
    end do

  END FUNCTION ext_theis2d_grid_sp

  FUNCTION ext_theis2d_grid_dp(grid, params, inits, prop, grad, steh_n, output)

    IMPLICIT NONE

    REAL(dp), DIMENSION(:,:),           INTENT(IN)  :: grid
    REAL(dp), DIMENSION(4),             INTENT(IN)  :: params               !params=(TG,sigma**2,corr,S)
    REAL(dp), DIMENSION(1),             INTENT(IN)  :: inits                !inits=(Qw)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality-factor [def=2]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=50]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: output               !enables output while computation

    REAL(dp), DIMENSION(size(grid(:,1)),size(grid(:,1))) :: temphead
    REAL(dp), DIMENSION(size(grid(:,1)))            :: ext_theis2d_grid_dp

    REAL(dp)                                        :: pr
    INTEGER(i4)                                     :: stehlim, limit, i

    LOGICAL                                         :: output2

    ! pr:               is the intern proportionality factor
    ! C:                abbreviation for the factor C=prop/corr
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! mask:             is the mask for all non-zero time-points

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (size(grid(1,:)) .ne. 2_i4)      stop 'grid should contain 2 information: radius and time'

    if (present(output)) then
       output2      = output
    else
       output2      = .false.
    endif

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
       limit        =   40_i4
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    temphead            =   ext_theis2d_dp_d1(&
         rad     =   grid(:,1),&
         time    =   grid(:,2),&
         params  =   params,&
         inits   =   inits,&
         prop    =   pr,&
         grad    =   limit,&
         steh_n  =   stehlim,&
         output  =   output2&
         )

    do i=1_i4, size(ext_theis2d_grid_dp)
       ext_theis2d_grid_dp(i)          =   temphead(i,i)
    end do

  END FUNCTION ext_theis2d_grid_dp

  FUNCTION ext_theis2d_sp_d1(rad, time, params, inits, prop, grad, steh_n, output)

    IMPLICIT NONE

    REAL(sp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(sp), DIMENSION(:),             INTENT(IN)  :: time
    REAL(sp), DIMENSION(4),             INTENT(IN)  :: params               !params=(TG,sigma**2,corr,S)
    REAL(sp), DIMENSION(1),             INTENT(IN)  :: inits                !inits=(Qw)
    REAL(sp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality-factor [def=2]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=50]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: output               !enables output while computation

    REAL(sp), DIMENSION(size(rad),size(time))       :: ext_theis2d_sp_d1

    REAL(sp)                                        :: pr
    INTEGER(i4)                                     :: stehlim, limit!, n, m

    LOGICAL                                         :: output2

    ! pr:               is the intern proportionality factor
    ! C:                abbreviation for the factor C=prop/corr
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! mask:             is the mask for all non-zero time-points

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (present(output)) then
       output2      = output
    else
       output2      = .false.
    endif

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
       limit        =   40_i4
    endif

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    ext_theis2d_sp_d1   =   real(ext_theis2d_dp_d1(&
         rad     =   real(rad,dp),&
         time    =   real(time,dp),&
         params  =   real(params,dp),&
         inits   =   real(inits,dp),&
         prop    =   real(pr,dp),&
         grad    =   limit,&
         steh_n  =   stehlim,&
         output  =   output2&
         ),sp)

  END FUNCTION ext_theis2d_sp_d1


  FUNCTION ext_theis2d_dp_d1(rad, time, params, inits, prop, grad, steh_n, output)

    use mo_laplace_inversion,    ONLY: NLInvSteh

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: rad
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: time
    REAL(dp), DIMENSION(4),             INTENT(IN)  :: params               !params=(TG,sigma**2,corr,S)
    REAL(dp), DIMENSION(1),             INTENT(IN)  :: inits                !inits=(Qw)
    REAL(dp),               OPTIONAL,   INTENT(IN)  :: prop                 !proportionality-factor [def=2]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: grad                 !limit of the series expansion [def=40]
    INTEGER(I4),            OPTIONAL,   INTENT(IN)  :: steh_n               !The Boundary for the stehfest-alg. [def=6]
    LOGICAL ,               OPTIONAL,   INTENT(IN)  :: output               !enables output while computation

    REAL(dp), DIMENSION(size(rad),size(time))       :: ext_theis2d_dp_d1

    REAL(dp)                                        :: pr, output_real, C !C=prop/corr
    REAL(dp), DIMENSION(8)                          :: values
    INTEGER(i4)                                     :: stehlim, limit, n!, m

    LOGICAL, DIMENSION(size(rad), size(time))       :: mask
    !    LOGICAL                                         :: zero

    ! pr:               is the intern proportionality factor
    ! C:                abbreviation for the factor C=prop/corr
    ! values:           Array of parameters to use the Stehfest-algorithm

    ! stehlim:          is the intern boundary for the stehfest-algorithm
    ! limit:            is the intern limit of the series expansion

    ! mask:             is the mask for all non-zero time-points

    ! ------------------------------------------------------------------
    ! Input tests
    ! ------------------------------------------------------------------

    if (.not. all(rad>0.0_dp))      stop 'Radius must be positiv.'
    if (any(time <    0.0_dp))      stop 'The time must be non-negativ.'
    if (ge(inits(1) , 0.0_dp))      stop 'The pumping-rate must negativ.'
    if (le(params(1), 0.0_dp))      stop 'The transmissivity must be positiv.'
    if (le(params(2), 0.0_dp))      stop 'The variance must be positiv.'
    if (le(params(3), 0.0_dp))      stop 'The correlation length must be positiv.'
    if (le(params(4), 0.0_dp))      stop 'The storage must be positiv.'
    !    do n=2_i4, size(rad)
    !        if (rad(n) < rad(n-1))      stop 'Radial-points need to be sorted.'
    !    end do

    if (present(output)) then
       if (output) then
          output_real  = 1.0_dp
       else
          output_real  = 0.0_dp
       endif
    else
       output_real  = 0.0_dp
    endif

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
       limit        =   40_i4
    endif

    ! ------------------------------------------------------------------
    ! setting of the intern variables
    ! ------------------------------------------------------------------

    C               =   pr/params(3)

    values          =   (/params(1), params(2), C, params(4), inits(1), real(limit,dp), real(stehlim,dp), output_real/)

    !the optional output
    if (nint(output_real,i4)==1_i4) then
       write (*,*) ""
       write (*,*) "TG  : ",params(1)
       write (*,*) "sig2: ",params(2)
       write (*,*) "corr: ",params(3)
       write (*,*) "S   : ",params(4)
       write (*,*) "Qw  : ",inits(1)
       write (*,*) "prop: ", pr
    end if

    ! ------------------------------------------------------------------
    ! calculate the head
    ! ------------------------------------------------------------------

    do n=1_i4, size(rad)
       mask(n,:)       =   (time>0.0_dp)
    end do

    ext_theis2d_dp_d1   =   0.0_dp

    ext_theis2d_dp_d1   =&
         unpack(&
         pack(NLInvSteh(lap_ext_theis2d_dp_d1, values, C*rad, pack(time, (time>0.0_dp)), stehlim), .true.),&
         mask,&
         ext_theis2d_dp_d1&
         )

    where( ext_theis2d_dp_d1 .gt. 0.0_dp )          ext_theis2d_dp_d1       =   0.0_dp
    where( abs(ext_theis2d_dp_d1) .lt. 1e-3_dp )    ext_theis2d_dp_d1       =   0.0_dp

    !    do n=1_i4, size(time)
    !        zero             =   .false.
    !        do m=1_i4, size(rad)
    !            if (eq(ext_theis2d_dp_d1(m,n),0.0_dp))  zero                    =   .true.
    !            if (zero)                               ext_theis2d_dp_d1(m,n)  =   0.0_dp
    !        end do
    !    end do

  END FUNCTION ext_theis2d_dp_d1


  FUNCTION lap_ext_theis2d_dp_d1(r, s, para)

    use mo_functions,       only: factorial
    use mo_nrutil,          only: poly
    use mo_combinatorics,   only: nextpart
    use mo_nr,              only: bessk0, bessk1
    use mo_interpol,        only: interpol
    use mo_ode_solver,      only: RBstiff
    use mo_linear_algebra,  only: solve_linear_equations

    IMPLICIT NONE

    REAL(dp), DIMENSION(:),             INTENT(IN)  :: r
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: s
    REAL(dp), DIMENSION(:),             INTENT(IN)  :: para                     !para=(TG,sig^2,C,S,Qw,limit, stehlim,output)

    REAL(dp), DIMENSION(size(r), size(s))           :: lap_ext_theis2d_dp_d1

    !the standard limit for the series-expansion
    integer(i4), Parameter                          :: stdsiz=41_i4

    !series-expansion of the given functions
    REAL(dp), &
         DIMENSION(max(nint(para(6),i4)+1,stdsiz))  :: dg, dh, g, h, alpha, beta, innerfunc

    !All partitionarrays in a row to calculate beta (notice that sum(i,1,n)=n*(n+1)/2)
    integer(i4), &
         DIMENSION((size(beta)*(size(beta)-1))/2)   :: partitions

    !functions for the rosenbrock-integration
    REAL(dp), DIMENSION(2)                          :: yi                   ! initial conditions
    REAL(dp), DIMENSION(:),     ALLOCATABLE         :: g1_xout, g2_xout     ! output time
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE         :: g1_out,  g2_out      ! output solutions
    REAL(dp), DIMENSION(1)                          :: g1temp,  g2temp      ! fundamental solutions evaluated at r(rad_n)
    REAL(dp), DIMENSION(5)                          :: g_para               ! parameters for the functions

    !variables for the rosenbrock-integration
    REAL(dp)                                        :: rinf,rmin,rmax,rlim  ! Boundary
    REAL(dp)                                        :: hstep                ! guessed stepsize
    REAL(dp)                                        :: hmin                 ! minimum allowed stepsize
    REAL(dp)                                        :: eps                  ! accurancy

    !intern variables
    REAL(dp)                                        :: prod, T, coef, TH, Tmid, t_temp, error, ulim
    REAL(dp), DIMENSION(2,2)                        :: Mat                  !Matrix for the LES to clue the solutions together
    REAL(dp), DIMENSION(2)                          :: LHS, RHS             !Left- and Righthandside of the LES

    !loop-variables
    INTEGER(i4)                                     :: n, m, rad_n, laps_n, nn, mm

    !intern variable-names from para
    REAL(dp)                                        :: TG, sig2, C, St, Qw, n_series
    INTEGER(i4)                                     :: n_steh
    LOGICAL                                         :: output

    !define the intern variables
    TG              =   para(1)
    sig2            =   para(2)
    C               =   para(3)
    St              =   para(4)
    Qw              =   para(5)
    n_series        =   para(6)                 !is real
    n_steh          =   nint(para(7),i4)
    output          =   nint(para(8),i4)==1_i4


    !Calculate the series-expansion of alpha
    alpha           =   0.0_dp
    alpha(1)        =   1.0_dp

    do n=0_i4, min(Floor(n_series/2.0_dp - 1.0_dp), ((stdsiz-1_i4)/2_i4 - 1_i4))
       alpha(2*(n+1)+1)        =   ((-1.0_dp)**n)*sig2*real(n+1_i4,dp)
    end do

    !Calculate the series-expansion of beta=-r^2*St/TCG(r)/C^2
    beta            =   0.0_dp
    partitions      =   0_i4

    !set all partitions to the "biggest" one
    do n=1_i4,size(beta)-1_i4
       partitions((n*(n+1_i4))/2_i4)   =   1_i4
    end do

    !Define the series expansion of the inner function in beta
    innerfunc       =   0.0_dp
    innerfunc(1)    =   0.5_dp*sig2

    do n=1_i4, min(Floor(n_series/2.0_dp), ((stdsiz-1_i4)/2_i4))
       innerfunc(2_i4*n+1_i4)  =   0.5_dp*((-1.0_dp)**n)*sig2
    end do

    beta(3)         =   1.0_dp

    !Define the series expansion of beta by the rule of Faa di Bruno
    do n=2_i4,min(2_i4*Floor(n_series/2.0_dp), stdsiz-1_i4)-2_i4,2_i4
       do
          prod=1.0_dp
          do m=1_i4, n
             if ((partitions(((n-1)*n)/2+m) .ne. 0_i4)) then
                if  (ne(innerfunc(m+1), 0.0_dp)) then
                   prod    =&
                        prod*(innerfunc(m+1)**partitions(((n-1)*n)/2+m))/real(factorial(int(partitions(((n-1)*n)/2+m),i8)),dp)
                else
                   prod    =   0.0_dp
                   exit
                end if
             end if
          end do

          beta(n+3)           =   beta(n+3)+prod

          if (partitions(((n-1)*n)/2+1)==n) exit

          partitions(((n-1)*n)/2+1:((n+1)*n)/2)    =   nextpart(partitions(((n-1)*n)/2+1:((n+1)*n)/2))

       end do
    end do

    beta                =   -(St/TG/C**2_i4)*exp(0.5_dp*sig2)*beta

    !Define the boundaries for the different solutions
    rinf                =   0.50_dp
    ulim                =   0.25_dp

    !define the upper limit for the solution such that the relativ error between TCG(r) and TG is less then 1% (error)
    error               =   0.01_dp
    rlim                =   sqrt( -0.5_dp*sig2/log(1.0_dp-error) - 1.0_dp )

    !define a weighted mean between TH and TG
    TH                  =   TG*exp(-sig2*0.5_dp)                                !nearfield transmissivity
    Tmid                =   (TG+2.0_dp*TH)/3.0_dp

    !parameters for the rosenbrock-integration with adaptive stepsize
    hstep               =   1e-8_dp                            ! incremental time step
    hmin                =   0.0_dp                             ! minimum allowed stepsize
    eps                 =   1e-8_dp                            ! accuracy

    !optional output
    if (output) then
       write (*,*) ""
       write (*,*) "ulim: "    , ulim
       write (*,*) "Rinf: "    , rinf/C, "Rinf*C: ", rinf
       write (*,*) "Rlim: "    , rlim/C, "Rlim*C: ", rlim
       write (*,*) "T_err:", error*100.0_dp,"percent"
       write (*,*) "C:    "    , C
       write (*,*) "sig^2:"    , sig2
       write (*,*) "TG:   "    , TG     , "TH:     ", TH
       write (*,*) "Tlim: " , TG*exp(-0.5_dp*(sig2)/(1.0_dp+(rlim)**2_i4))
       write (*,*) "Tmid: "    , Tmid
       write (*,*) ""
    endif

    !loop over s(:)
    do nn=1_i4, size(s)/n_steh
       !calculate the actual timepoint
       t_temp              =   1.0_dp/(s(nn)/log(2.0_dp))

       !Define the boundaries by comparisson to theis with the value U=r^2*St/(4*L*Tmid*t) < ulim
       rmax                =   min(sqrt(ulim*4.0_dp*t_temp*Tmid/St)*C, rlim)
       !        rmax                =   sqrt(ulim*4.0_dp*t_temp*Tmid/St)*C
       !        rmax                =   rlim

       rmin                =   min(rmax/2.0_dp, rinf)

       if (output) then
          write (*,*) "t_temp: ",     t_temp,&
               " rmin: ",      rmin/C,&
               " rmax: ",      rmax/C,&
               " eps: ",       eps,&
               " T(rmax): ",   TG*exp(-0.5_dp*(sig2)/(1.0_dp+(rmax)**2_i4))
       endif

       do mm=0_i4, n_steh-1_i4
          !go through s by every time-point
          laps_n              =   mm*size(s)/n_steh + nn

          !Define the series expansion of the first fundamental solution g with the aid of alpha and beta by the frobenius-method
          g                   =   0.0_dp
          g(1)                =   1.0_dp
          dg                  =   0.0_dp

          do n=1_i4,size(g)-1_i4
             do m=0_i4,n-1_i4
                g(n+1)      =   g(n+1)+(real(m,dp)*alpha(n+1-m)+s(laps_n)*beta(n+1-m))*g(m+1)
             end do
             g(n+1)          =   -1.0_dp/(real(n,dp)**2_i4)*g(n+1)
             dg(n)           =   real(n,dp)*g(n+1)
          end do

          !Define the series of the the first part of the second fundamental solution h with the aid of alpha, beta & g by froben.
          h                   =   0.0_dp
          h(1)                =   0.0_dp  !1.0_dp
          dh                  =   0.0_dp

          do n=1_i4,size(h)-1_i4
             do m=0_i4,n-1_i4
                h(n+1)      =   h(n+1) + alpha(n+1-m)*(real(m,dp)*h(m+1)+g(m+1)) + s(laps_n)*beta(n+1-m)*h(m+1)
             end do
             h(n+1)          =   -1.0_dp/(real(n,dp)**2_i4)*(h(n+1)+2.0_dp*real(n,dp)*g(n+1))
             dh(n)           =   real(n,dp)*h(n+1)
          end do

          !Calculate 2 fundamental solutions in [rmin,rmax] with Rosenbrock-integration
          g_para              =   (/sig2, s(laps_n), St, TG, C/)

          !boundary-condition at the well
          coef                =   -Qw/(2.0_dp*PI_D*TH*s(laps_n))

          ! 1. solution
          yi                  =   (/  poly(rmin,g),&
               poly(rmin,dg) /)

          call RBstiff( yi, rmin, rmax, hstep, g2D_dp, dg2D_dp, g_para, g1_xout, g1_out, hmin=hmin, eps=eps)

          ! 2. solution
          yi                  =   (/  coef*(poly(rmin,h)  + log(rmin)*poly(rmin,g)),&
               coef*(poly(rmin,dh) + log(rmin)*poly(rmin,dg) + poly(rmin,g)/rmin) /)

          call RBstiff( yi, rmin, rmax, hstep, g2D_dp, dg2D_dp, g_para, g2_xout, g2_out, hmin=hmin, eps=eps)

          !clue the solution together with theis in the farfield (modified bessel-functions)
          Mat                 =   0.0_dp
          RHS                 =   0.0_dp

          !        T                   =   sqrt(s(laps_n)*St/(TG*exp(-0.5_dp*(sig2)/(1.0_dp+(rmax)**2_i4))))/C
          T                   =   sqrt(s(laps_n)*St/(TG))/C

          !calculate the coefficients (As Bs Cs Ds)
          Mat(1,1:1)          =   interpol(g1_out(:,1) , g1_xout, (/rmax/))
          Mat(1,2)            =   -bessk0(T*rmax)
          Mat(2,1:1)          =   interpol(g1_out(:,2) , g1_xout, (/rmax/))
          Mat(2,2)            =   T*bessk1(T*rmax)

          RHS(1:1)            =   -interpol(g2_out(:,1) , g2_xout, (/rmax/))
          RHS(2:2)            =   -interpol(g2_out(:,2) , g2_xout, (/rmax/))

          LHS                 =   solve_linear_equations(Mat,RHS)

          !calculate the head in laplace-space
          do rad_n=1_i4, size(r)
             if  (r(rad_n)<rmin) then
                lap_ext_theis2d_dp_d1(rad_n, laps_n)  =   (coef*log(r(rad_n)) + LHS(1))*poly(r(rad_n),g) + coef*poly(r(rad_n),h)
             else if (r(rad_n)>rmax) then
                lap_ext_theis2d_dp_d1(rad_n, laps_n)  =   LHS(2)*bessk0(T*r(rad_n))
             else
                g1temp                                =   interpol(g1_out(:,1) , g1_xout, (/r(rad_n)/))
                g2temp                                =   interpol(g2_out(:,1) , g2_xout, (/r(rad_n)/))
                lap_ext_theis2d_dp_d1(rad_n, laps_n)  =   LHS(1)*g1temp(1) + g2temp(1)
             end if
          end do

          !laps_n
       end do  !m
    end do  !n

  END FUNCTION lap_ext_theis2d_dp_d1

  subroutine g2D_dp( x, y, para1, dydx )   ! returns derivatives dydx at x for the first fundamental solution g

    implicit none

    real(dp), intent(in) :: x                       ! time
    real(dp), dimension(:), intent(in) :: y         ! unknowns of the equations
    real(dp), dimension(:), intent(in) :: para1      ! para1=(sigma**2, s, S, TG, C)
    real(dp), dimension(:), intent(out) :: dydx     ! derivatives of y

    dydx(1)     =   y(2)
    dydx(2)     =   (-para1(1)*x/((1+x**2_i4)**2_i4) - 1.0_dp/x)                                 *   y(2) &
         + (exp(0.5_dp*para1(1)/(1+x**2_i4))*para1(2)*para1(3)/para1(4)/((para1(5))**2_i4)) *   y(1)

  end subroutine g2D_dp

  subroutine dg2D_dp( x, y, para2, dfdx, dfdy )   ! returns derivatives dydx at x for the first fundamental solution g

    implicit none

    real(dp), intent(in) :: x                       ! time
    real(dp), dimension(:), intent(in) :: y         ! unknowns of the equations
    real(dp), dimension(:), intent(in) :: para2      ! para2=(sigma**2, s, S, TG, C)
    real(dp), dimension(:), intent(out) :: dfdx     ! derivatives of f
    real(dp), dimension(:,:), intent(out) :: dfdy   ! derivatives of f

    dfdx(1)     =   0.0_dp
    dfdx(2)     =   ( -para2(1)*(1.0_dp - 3.0_dp*x**2_i4)/((1+x**2_i4)**3_i4) + 1.0_dp/(x**2_i4) )   *   y(2) &
         + ( exp(0.5_dp*para2(1)/(1+x**2_i4))*para2(2)*para2(3)/para2(4)/((para2(5))**2_i4)*&
         para2(1)*x/((1+x**2_i4)**2_i4) )                                               *   y(1)

    dfdy(1,1)   =   0.0_dp
    dfdy(1,2)   =   1.0_dp
    dfdy(2,1)   =   (exp(0.5_dp*para2(1)/(1+x**2_i4))*para2(2)*para2(3)/para2(4)/((para2(5))**2_i4))
    dfdy(2,2)   =   (-para2(1)*x/((1+x**2_i4)**2_i4) - 1.0_dp/x)

  end subroutine dg2D_dp

END MODULE mo_pumpingtests

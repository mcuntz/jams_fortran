!> \file mo_ansi_colors.f90

!> \brief Provide colouriser to output coloured text on terminal.

!> \details Function that surrounds string with ANSI color codes to output coloured text on terminal

!>  From Fortran Wiki: http://fortranwiki.org/fortran/show/ansi_colors

!> \authors Jason Blevins
!> \date April 2016
module mo_ansi_colors

  ! This module provides ANSI colors to output coloured terminal output.
  ! The module is originally from Jason Blevis, published at the Fortran Wiki:
  !     http://fortranwiki.org/fortran/show/ansi_colors
  ! and was adapted to the JAMS Fortran library.

  ! Written Matthias Cuntz, Feb 2020

  ! License
  ! -------
  ! This file is part of the JAMS Fortran package, distributed under the MIT License.
  !
  ! Copyright (c) 2016-2020 Jason Blevis, Matthias Cuntz - mc (at) macu (dot) de
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

  implicit none

  private

  public :: color ! Function returning str surrounded by color codes
  
  ! ansi colors black, red, green, yellow, blue, magenta cyan, and white
  character(len=*), parameter, public :: c_black   = '30'
  character(len=*), parameter, public :: c_red     = '31'
  character(len=*), parameter, public :: c_green   = '32'
  character(len=*), parameter, public :: c_yellow  = '33'
  character(len=*), parameter, public :: c_blue    = '34'
  character(len=*), parameter, public :: c_magenta = '35'
  character(len=*), parameter, public :: c_cyan    = '36'
  character(len=*), parameter, public :: c_white   = '37'

  character(len=1), parameter :: c_esc   = achar(27)
  character(len=2), parameter :: c_start = c_esc // '['
  character(len=1), parameter :: c_end   = 'm'
  character(len=*), parameter :: c_clear = c_start // '0' // c_end

contains

  ! ------------------------------------------------------------------
  !
  !     NAME
  !         color
  !
  !     PURPOSE
  !>        \brief Colorizer of strings
  !
  !>        \details Return string surrounded by color codes, which can be passsed to print or write statements
  !>        to produce coulered terminal output.
  !
  !     CALLING SEQUENCE
  !         out = color(str, code)
  !
  !     INTENT(IN)
  !>        \param[in] "character(len=*) :: str"   String
  !>        \param[in] "character(len=*) :: code"  Colour code string from mo_ansi_colors module.\n
  !>                                               Available codes are:\n
  !>                                               c_black, c_red, c_green, c_yellow, c_blue, c_magenta c_cyan, and c_white
  !
  !     RETURN
  !>       \return     character(len=len(str)+9) :: clor &mdash; str surrounded by color codes.
  !
  !     EXAMPLE
  !         write(*,*) color('o.k.', c_green)
  !
  !     HISTORY
  !>        \authors Jason Blevins, Matthias Cuntz
  !>        \date Apr 2016

  function color(str, code) result(out)

    implicit none

    character(len=*), intent(in)  :: str
    character(len=*), intent(in)  :: code
    character(len=:), allocatable :: out

    out = c_start // code // c_end // str // c_clear

  end function color

end module mo_ansi_colors

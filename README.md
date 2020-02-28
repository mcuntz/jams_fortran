# This is the JAMS Fortran library.

The JAMS Fortran library is a collection of general Fortran modules
offering miscellaneous functions in the categories
[Computational](#comp), [Date and Time](#date), [Input / Output](#io),
[Math / Numerics](#math), [Miscellaneous](#misc), and
[Screening, Sensitivity Analysis and Optimising / Fitting](#opti).

Created November 2011 by Matthias Cuntz and Juliane Mai  
while at the Department Computational Hydrosystems, Helmholtz Centre
for Environmental Research - UFZ, Permoserstr. 15, 04318 Leipzig,
Germany, and continued while M Cuntz is at the Institut National de
Recherche pour l’Agriculture, l’Alimentation et l’Environnement
(INRAE), Nancy, France, and J Mai is at the Department of Civil and
Environmental Engineering of the University of Waterloo, ON, Canada.

It is distributed under the MIT License (see LICENSE file and below).

Copyright (c) 2011-2020 Matthias Cuntz and Juliane Mai

Contact Matthias Cuntz - mc (at) macu (dot) de


---------------------------------------------------------------

### Installation

The library is maintained on a git repository at:

    https://github.com/mcuntz/jams_fortran/

To use it, checkout the git repository

    git clone https://github.com/mcuntz/jams_fortran.git

Copy the mo\_\*.f90 file of the chosen utility to your Fortran
project. _grep -i 'use' mo\_\*.f90_ to check for dependencies. All
Fortran modules use at least the kind definitions of _mo\_kind.f90_;
quite a few use also the general utilities in _mo\_utils.f90_.


---------------------------------------------------------------

### Documentation

Documentation of the functions is above each function code or above the
module procedure interface at the top of the file. There is a very
brief introduction to the modules at the beginning of each file.

Some but not all documentation is in doxygen format and can be
formatted with the appropriate tools, e.g. using the
[JAMS Makefile project](https://github.com/mcuntz/jams_makefile).

Each module has a test program in the directory
_test/test\_mo\_\*/_. These test programs try to test all
functionalities of the modules and provide hence examples of code
usage of (almost) all subroutines.


---------------------------------------------------------------

### Dependencies

The library is compatible with the Fortran 2003 standard. It was
tested with the GNU gfortran, Intel ifort, NAG nagfor, and
PGI pgfortran compilers in various revisions.

Interoperability with C is assured by the use of the intrinsic module
_iso\_c\_binding_ in _mo\_kind.f90_. The use of C functions were tested using
the C-header file _cfortran.h_, for example with the `qhull` library. See
the use of _cfortran.h_ in _qhull/qhull\_c.c_ and the corresponding
Fortran call in _mo\_qhull.f90_.

The library uses C preprocessor flags to deal with special compiler
behaviours, i.e. one has to compile with the -cpp or -fpp compiler
flag. For example, gfortran (revisions < 5) did not include an
intrinsic _ieee\_arithmetic_ module so that _mo\_utils_ has convenience
funtions such as _is\_nan(x)_, which uses either the gfortran extension
_isnan(x)_ or the Fortran standard _ieee\_is\_nan(x)_. The decision is
made during compile time, checking if the variable `__GFORTRAN__` is
defined. The following preprocessor variables are currently checked:  
`__ABSOFT__`, `__GFORTRAN__`, `__GFORTRAN41__`, `__NETCDF3__`,
`__pgiFortran__`, `__NAG__`, `__NAGf90Fortran__`.

Some modules use third-party packages such as LAPack or netCDF.  
The libraries [MinPack](https://www.netlib.org/minpack/),
[QHull](http://www.qhull.org), and
[netCDF3](https://www.unidata.ucar.edu/software/netcdf/release-notes-3.6.3.html)
are given for convenience with the JAMS library. It is always
recommended to use specific installations for the computer system, for
example for blas and lapack, which are highly optimised. The netcdf3
directory includes only the parts of the netCDF 3.6.3 version that are
necessary for Fortran projects. Read the README in the netcdf3
directory for two system-specific configuration files. The
[JAMS Makefile project](https://github.com/mcuntz/jams_makefile)
handles these kind of special configuration files.

The only other third-party softwares are
[LAPack](http://www.netlib.org/lapack/) and
[netCDF 4](https://www.unidata.ucar.edu/software/netcdf/), which have
to be installed independently. The latter can be installed by using,
for example the installation script _install\_netcdf_ of
[https://github.com/mcuntz/install_netcdf](https://github.com/mcuntz/install_netcdf).


---------------------------------------------------------------

###  Content

Modules are provided in the following categories:  
 * [Computational](#comp)  
 * [Date and Time](#date)  
 * [Input / Output](#io)  
 * [Math /  Numerics](#math)  
 * [Screening, Sensitivity Analysis and Optimising / Fitting](#opti)  
 * [Miscellaneous](#misc)  


---------------------------------------------------------------

### Modules ordered alphabetically

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | mo\_anneal | Simulated Annealing with estimation of initial temperature |
| | mo\_ansi\_colors | Provide colouriser to output coloured text on terminal |
| | mo\_append | Appends vectors and/or matrixes like \*nix *cat* and *paste* |
| | mo\_bootstrapping\_sensitivity\_analysis.f90 | Bootstrapping for sensitivity analysis |
| | mo\_boxcox | Box-Cox transformation, inverse transformation & estimating best exponent for transformation |
| | mo\_combinatorics | Combinatorial algorithms, e.g. binomial coefficient and k-subsets |
| | mo\_constants | Mathematical and physical constants |
| | mo\_corr | Correlation function with Fast Fourier Transform, or auto- & cross-correlations via direct calculation |
| | mo\_dds | Dynamically Dimensioned Search (DDS) and Modified DDS (MDDS) |
| | mo\_delsa | Distributed Evaluation of Local Sensitivity Analysis (DELSA) |
| | mo\_distributions.f90 | Continuous distributions density functions such as for the normal distribution |
| | mo\_elemeffects | Calculation of parameter's Elementary Effect on model/function output |
| | mo\_errormeasures | Distance or error measures between two datasets, e.g. bias, RMSE, ... |
| | mo\_file\_utils | Utilities for file handling, e.g. search free unit |
| | mo\_finish | Routine to end program gracefully |
| | mo\_fit | Linear & polynomial fitting, and general fit with singular value decomposition |
| | mo\_functional.f90 | Functional programming for modern Fortran |
| | mo\_functions | Special functions such as the Gamma function |
| | mo\_histo | Histogram of data (useable also for variogram) |
| | mo\_integrate | Integration routines |
| | mo\_interpol | Linear interpolation for irregular grids |
| | mo\_julian | Converts dates between calendars, e.g. Julian Day into Day, Month and Year, and vice versa; Standard and IMSL convention |
| | mo\_kernel | Kernel regression and kernel density estimation for PDF and CDF |
| | mo\_kind | Definition of numerical precisions |
| | mo\_laplace\_inversion.f90 | Numerical Laplace inversion |
| | mo\_linear\_algebra | Wrapper functions for LAPACK's F77 linear algebra routines + some convenience functions such as the diagonal of a matrix |
| | mo\_linfit | Fit a straight line with model I or model II (geometric mean) regression without error bars on input |
| | mo\_mad | Median absolute deviation test |
| | mo\_mcmc | Monte Carlo Markov Chain sampling of parameter distribution around optimum |
| | mo\_message | Write out message; works with num2str from mo\_string\_utils |
| | mo\_minpack | Optimization package minpack (F90 interfaces) |
| | mo\_moment | 1st to 4th moments, central and mixed central moments |
| | mo\_mpi\_stubs.f90 | Dummy functions for common MPI routines |
| | mo\_mtclim.f90 | Peter E Thornton's climate generator |
| | mo\_ncdump.f90 | Quick dumping of variables into netcdf files |
| | mo\_ncread | Reading nc files using the netcdf4 library |
| | mo\_ncwrite | Writing nc files using the netcdf4 library |
| | mo\_nelmin | Minimizes a function using the Nelder-Mead algorithm with the Applied Statistics algorithms No. 047 |
| | mo\_netcdf | David Schaefer's NetCDF Fortran 90 interface wrapper, with NETCDF3 support |
| | mo\_nml | Routines to handle namelist files |
| | mo\_nr | Main numerical recipes module containing the interfaces |
| | mo\_nrutil | Numerical recipes utilities module |
| | mo\_ode\_generator | Given N reactants, generates and solves all corresponding ODE systems |
| | mo\_ode\_solver | Iterative methods for the approximation of solutions of Ordinary Differential Equations (ODE) |
| | mo\_opt\_functions | Test functions for optimization routines |
| | mo\_orderpack | Orderpack 2.0 from Michel Olagnon provides order and unconditional, unique, and partial ranking, sorting, and permutation. Provides also convenience routines sort and sort\_index |
| | mo\_percentile | Median, Percentiles |
| | mo\_pi\_index | Parameter importance index PI or alternatively B index calculation |
| | mo\_poly | Tests if a given 2D point lies inside, outside, or on edge/vertex of a 2D polygon, compute area and center of mass |
| | mo\_pumpingtests.f90 | Thiem and Theis solutions for groundwater-flow equation, and effektive well-flow solution for the effective coarse-graining conductivity |
| | mo\_qhull.f90 | Wrapper of C program for convex hull calculation from number of points |
| | mo\_quicksort | Several different implementations of Quicksort including an OpenMP version |
| | mo\_random\_field.f90 | Generation of random fields with certain statistical properties, e.g. a given correlation length |
| | mo\_remap | Remaps a grid to another grid |
| | mo\_sampling | Random and Latin Hypercube Sampling for a set of parameters with Uniform(0,1) or Gaussian(0,1) Distribution |
| | mo\_sce | Shuffled Complex Evolution optimisation |
| | mo\_select\_distant\_entries.f90 | Select elements of matrix with largest distance to each other |
| | mo\_sobol | Sampling of parameters using Sobol sequences |
| | mo\_sobol\_index | Sobol indexes (main and total effects) |
| | mo\_sort | Quicksort arrays or indices |
| | mo\_spatialsimilarity.f90 | Routines for bias insensitive comparison of spatial patterns |
| | mo\_specan | Spectral analysis using FFT |
| | mo\_spline | Spline functions to approximate or interpolate data |
| | mo\_standard_score.f90 | Normalized (anomaly)/standard score/z score and deseasonalized (standard score on monthly basis) values of a time series |
| | mo\_string\_utils | Utilities for strings |
| | mo\_tee.f90 | Write out concatenated strings on standard out and to a given file or unit |
| | mo\_template | Module template demonstrating the coding standard of the library |
| | mo\_timer | CPU time routines, allows setting of multiple CPU timers |
| | mo\_utils | Provides general utilities such as comparisons of two reals, swapping of two elements in an array, is\_nan, etc. |
| | mo\_var2nc.f90 | Writing variables to netcdf files |
| | mo\_xor4096 | Generating Uniform or Gaussian Random Numbers using the xor4096 algorithm |
| | mo\_xor4096\_apps | Wrapper functions for random number generator xor4096 (Arrays of RNs, ranged RNs, Multivariate Normal Distribution) |


---------------------------------------------------------------

### Modules per category 

<a name="comp"></a>**Computational**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | mo\_append | Appends vectors and/or matrixes like \*nix *cat* and *paste* |
| | mo\_constants | Mathematical and physical constants |
| | mo\_functional.f90 | Functional programming for modern Fortran |
| | mo\_kind | Definition of numerical precisions |
| | mo\_message | Write out message; works with num2str from mo\_string\_utils |
| | mo\_mpi\_stubs.f90 | Dummy functions for common MPI routines |
| | mo\_nml | Routines to handle namelist files |
| | mo\_nr | Main numerical recipes module containing the interfaces |
| | mo\_nrutil | Numerical recipes utilities module |
| | mo\_orderpack | Orderpack 2.0 from Michel Olagnon provides order and unconditional, unique, and partial ranking, sorting, and permutation. Provides also convenience routines sort and sort\_index |
| | mo\_quicksort | Several different implementations of Quicksort including an OpenMP version |
| | mo\_remap | Remaps a grid to another grid |
| | mo\_sort | Quicksort arrays or indices |
| | mo\_string\_utils | Utilities for strings |
| | mo\_template | Module template demonstrating the coding standard of the library |
| | mo\_timer | CPU time routines, allows setting of multiple CPU timers |
| | mo\_utils | Provides general utilities such as comparisons of two reals, swapping of two elements in an array, is\_nan, etc. |
| | mo\_finish | Routine to end program gracefully |


<a name="date"></a>**Date and Time**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | mo\_julian | Converts dates between calendars, e.g. Julian Day into Day, Month and Year, and vice versa; Standard and IMSL convention |


<a name="io"></a>**Input / Output**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | mo\_ansi\_colors | Provide colouriser to output coloured text on terminal |
| | mo\_file\_utils | Utilities for file handling, e.g. search free unit |
| | mo\_ncdump.f90 | Quick dumping of variables into netcdf files |
| | mo\_ncread | Reading nc files using the netcdf4 library |
| | mo\_ncwrite | Writing nc files using the netcdf4 library |
| | mo\_netcdf | David Schaefer's NetCDF Fortran 90 interface wrapper, with NETCDF3 support |
| | mo\_tee.f90 | Write out concatenated strings on standard out and to a given file or unit |
| | mo\_var2nc.f90 | Writing variables to netcdf files |


<a name="math"></a>**Math / Numerics**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | mo\_boxcox | Box-Cox transformation, inverse transformation & estimating best exponent for transformation |
| | mo\_combinatorics | Combinatorial algorithms, e.g. binomial coefficient and k-subsets |
| | mo\_corr | Correlation function with Fast Fourier Transform, or auto- & cross-correlations via direct calculation |
| | mo\_distributions.f90 | Continuous distributions density functions such as for the normal distribution |
| | mo\_errormeasures | Distance or error measures between two datasets, e.g. bias, RMSE, ... |
| | mo\_functions | Special functions such as the Gamma function |
| | mo\_histo | Histogram of data (useable also for variogram) |
| | mo\_integrate | Integration routines |
| | mo\_interpol | Linear interpolation for irregular grids |
| | mo\_kernel | Kernel regression and kernel density estimation for PDF and CDF |
| | mo\_laplace\_inversion.f90 | Numerical Laplace inversion |
| | mo\_linear\_algebra | Wrapper functions for LAPACK's F77 linear algebra routines + some convenience functions such as the diagonal of a matrix |
| | mo\_linfit | Fit a straight line with model I or model II (geometric mean) regression without error bars on input |
| | mo\_mad | Median absolute deviation test |
| | mo\_moment | 1st to 4th moments, central and mixed central moments |
| | mo\_ode\_solver | Iterative methods for the approximation of solutions of Ordinary Differential Equations (ODE) |
| | mo\_percentile | Median, Percentiles |
| | mo\_poly | Tests if a given 2D point lies inside, outside, or on edge/vertex of a 2D polygon, compute area and center of mass |
| | mo\_qhull.f90 | Wrapper of C program for convex hull calculation from number of points |
| | mo\_sampling | Random and Latin Hypercube Sampling for a set of parameters with Uniform(0,1) or Gaussian(0,1) Distribution |
| | mo\_select\_distant\_entries.f90 | Select elements of matrix with largest distance to each other |
| | mo\_sobol | Sampling of parameters using Sobol sequences |
| | mo\_spatialsimilarity.f90 | Routines for bias insensitive comparison of spatial patterns |
| | mo\_spline | Spline functions to approximate or interpolate data |
| | mo\_standard_score.f90 | Normalized (anomaly)/standard score/z score and deseasonalized (standard score on monthly basis) values of a time series |
| | mo\_xor4096 | Generating Uniform or Gaussian Random Numbers using the xor4096 algorithm |
| | mo\_xor4096\_apps | Wrapper functions for random number generator xor4096 (Arrays of RNs, ranged RNs, Multivariate Normal Distribution) |


<a name="opti"></a>**Screening, Sensitivity Analysis and Optimising / Fitting**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | mo\_anneal | Simulated Annealing with estimation of initial temperature |
| | mo\_bootstrapping\_sensitivity\_analysis.f90 | Bootstrapping for sensitivity analysis |
| | mo\_dds | Dynamically Dimensioned Search (DDS) and Modified DDS (MDDS) |
| | mo\_delsa | Distributed Evaluation of Local Sensitivity Analysis (DELSA) |
| | mo\_elemeffects | Calculation of parameter's Elementary Effect on model/function output |
| | mo\_fit | Linear & polynomial fitting, and general fit with singular value decomposition |
| | mo\_mcmc | Monte Carlo Markov Chain sampling of parameter distribution around optimum |
| | mo\_minpack | Optimization package minpack (F90 interfaces) |
| | mo\_nelmin | Minimizes a function using the Nelder-Mead algorithm with the Applied Statistics algorithms No. 047 |
| | mo\_opt\_functions | Test functions for optimization routines |
| | mo\_pi\_index | Parameter importance index PI or alternatively B index calculation |
| | mo\_sce | Shuffled Complex Evolution optimisation |
| | mo\_sobol\_index | Sobol indexes (main and total effects) |


<a name="misc"></a>**Miscellaneous**

| | | |
| --- | --- | ------------------------------------------------------------------------------------------------------------- |
| | mo\_mtclim.f90 | Peter E Thornton's climate generator |
| | mo\_ode\_generator | Given N reactants, generates and solves all corresponding ODE systems |
| | mo\_pumpingtests.f90 | Thiem and Theis solutions for groundwater-flow equation, and effektive well-flow solution for the effective coarse-graining conductivity |
| | mo\_random\_field.f90 | Generation of random fields with certain statistical properties, e.g. a given correlation length |
| | mo\_specan | Spectral analysis using FFT |


---------------------------------------------------------------

###  License

This file is part of the JAMS Fortran library, distributed under the MIT License.

Copyright (c) 2011-2020 Matthias Cuntz, Juliane Mai, Stephan Thober - mc (at) macu (dot) de

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


---------------------------------------------------------------

## Note on Numerical Recipes License

Be aware that some code is under the Numerical Recipes License 3rd
edition
[http://www.nr.com/aboutNR3license.html](http://www.nr.com/aboutNR3license.html). 
You can identity it by the use of  mo\_nr or mo\_nrutil.

The Numerical Recipes Personal Single-User License lets you personally
use Numerical Recipes code ("the code") on any number of computers,
but only one computer at a time. You are not permitted to allow anyone
else to access or use the code. You may, under this license, transfer
precompiled, executable applications incorporating the code to other,
unlicensed, persons, providing that (i) the application is
noncommercial (i.e., does not involve the selling or licensing of the
application for a fee), and (ii) the application was first developed,
compiled, and successfully run by you, and (iii) the code is bound
into the application in such a manner that it cannot be accessed as
individual routines and cannot practicably be unbound and used in
other programs. That is, under this license, your application user
must not be able to use Numerical Recipes code as part of a program
library or "mix and match" workbench.

Businesses and organizations that purchase the disk or code download,
and that thus acquire one or more Numerical Recipes Personal
Single-User Licenses, may permanently assign those licenses, in the
number acquired, to individual employees. Such an assignment must be
made before the code is first used and, once made, it is irrevocable
and can not be transferred. 

If you do not hold a Numerical Recipes License, this code is only for
informational and educational purposes but cannot be used.

/*************************************************************************

 test program to demonstrate Fortran-C interoperability

*************************************************************************/

#include "mo_c.h"
#ifdef __CFORTRAN__
#include "cfortran.h"
#endif

/* Function that sums over some elements in the first dimension, and multiplying the result with a constant */
double testc(double *data, int d, int n, const double ref)
{
  int i;
  double out;

  out = 0.;
  for(i=d; i<=n; i++){
    out += data[i]*ref;
  }
  
  return out;
}
#ifdef __CFORTRAN__
/* Called as ctest in fortran; has 4 inputs.
   Put P in front of variable definition for pointers, i.e. INT -> int d, and PINT -> int *d
*/
FCALLSCFUN4(DOUBLE, testc, CTEST, ctest, PDOUBLE, INT, INT, DOUBLE)
#endif

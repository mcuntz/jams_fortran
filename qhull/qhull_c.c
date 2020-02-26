#include "qhull_a.h"
#ifdef __CFORTRAN__
#include "cfortran.h"
#endif

int qhull_c(char *flags, char *ofile, int dim, int numpoints, double *data) {
  
  coordT points[dim*numpoints]; /* array of coordinates for each point */
  boolT ismalloc = False;       /* True if qhull should free points in qh_freeqhull() or reallocation */
  FILE *outfile;                /* output from qh_produce_output(); use NULL to skip qh_produce_output() */
  FILE *errfile = stderr;       /* error messages from qhull code */
  int exitcode;                 /* 0 if no error from qhull */
  int curlong, totlong;         /* memory remaining after qh_memfreeshort */
  int i;

  /* Open output file */
  outfile = fopen(ofile, "w");
  
  /* Transform to coordT because cfortran knows only double */
  for (i=0; i<dim*numpoints; i++) points[i] = (coordT)data[i];

  /* Get the convex hull */
  exitcode = qh_new_qhull(dim, numpoints, points, ismalloc, flags, outfile, errfile);

  /* Clean up */
  qh_freeqhull(!qh_ALL);                   /* free long memory  */
  qh_memfreeshort (&curlong, &totlong);    /* free short memory and memory allocator */
  if (curlong || totlong)
    fprintf (errfile, "qhull_c: qhull internal warning (user_eg, #1): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
   
  /* Close output file */
  fclose(outfile);

  return exitcode;
}
#ifdef __CFORTRAN__
FCALLSCFUN5(INT, qhull_c, QHULL_F, qhull_f, STRING, STRING, INT, INT, PDOUBLE)
#endif

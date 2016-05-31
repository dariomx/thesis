#include <stdlib.h>
#include <stdbool.h>
#include "test_util.h"
#include "arpack.h"

int fiedler_arpack(int dim,
                   char * which,
                   int inev,
                   double isigma,
                   double itol,
                   void (matvec)(double *, double *),
                   void (solve)(double *, double *),
                   double * D,
                   double * Z,
                   int * iter);

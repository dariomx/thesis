#include "fiedler_arpack.h"

/**
 * Practical interface for ARPACK; simplified to run for standard
 * symmetric eigenvalue problem only.
 *
 * Returns the "exit code" of ARPACK; and if all went fine, it will
 * also write the requested eigenpairs (k) into lams and V (normally
 * we just want the first eigenpair, but this allows more flexibility).
 * allocated already). The output variable iter serves to save actual
 * number of iterations taken.
 *
 * Except for which, bmat, nev and tol; the rest of the parameters are
 * stolen from Scipy interface to ARPACK (arpack.py) 
 */
int fiedler_arpack(int dim,
                   char * which,
                   int inev,
                   double isigma,
                   double itol,
                   void (matvec)(double *, double *),
                   void (solve)(double *, double *),
                   double * D,
                   double * Z,
                   int * iter)
{
  int ido = 0;
  char * bmat = "I";
  int n = dim;
  int nev = inev;
  double tol = itol;
  double * resid = alloc_double(n);
  int ncv = 2*nev + 1;
  double * V = alloc_double(n * ncv);
  int ldv = n;
  int iparam[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  iparam[0] = 1;
  iparam[2] = n * 10;
  iparam[3] = 1;
  int mode = (isigma >= 0? 3 : 1); 
  iparam[6] = mode;
  int ipntr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double * workd = alloc_double(3*n);
  int lworkl = ncv * (ncv + 8);
  double * workl = alloc_double(lworkl);
  int info = 0;
  bool rvec = true;
  char * howmny = "A";
  bool * select = (bool *) malloc(nev * sizeof(bool));
  int ldz = n;
  double sigma = isigma;
  double *x, *y;
  
  while (TRUE)
  {
    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, V, &ldv,
           iparam, ipntr, workd, workl, &lworkl, &info);
    x = (double *) workd;
    y = ((double *) workd) + (ipntr[1] - 1);
    if ( ido == -1 || (ido == 1 && mode == 1) )
    {
      x += ipntr[0] - 1;
      matvec(x, y);
    }
    else if ( ido == 1 && mode == 3)
    {
      x += ipntr[2] - 1;
      solve(x, y);
    }
    else 
    {
      assert(ido == 99, "ARPACK iteration error in dsaupd: ido=%d", ido);
      break;
    }
  }

  dseupd_(&rvec, howmny, select, D, Z, &ldz, &sigma,
         bmat, &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam,
         ipntr, workd, workl, &lworkl, &info);
  assert(info == 0, "ARPACK deupd error: info=%d", info);
  *iter = iparam[3];
  return info;
}

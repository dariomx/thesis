#include "fiedler_arpack.h"

/*
 * tests fiedler_arpack + mult_matvec
 */

cholmod_sparse * L;
cholmod_dense * X;
cholmod_dense * Y;
cholmod_dense * FV;

void matvec(double * x, double * y)
{
  set_vec(x, X);
  mult_matvec(L, X, Y);
  copy_vec(Y, y);
}

int main(int argc, char * argv[])
{
  int n, iter, nev;
  double *lams, *V, ac, *fv, tol=1e-7, start, end;
  assert(argc == 3, "Usage: %s <lap_file> <fv_file>\n", argv[0]);
  test_start();
  L = load_mat(argv[1]);
  n = L->nrow;
  X = create_vec(n);
  Y = create_vec(n);
  nev = 2;
  lams = alloc_double(nev);
  V = alloc_double(nev * n);
  start = cputime();
  fiedler_arpack(n, "SM", nev, -1, tol, matvec, NULL, lams, V, &iter);
  end = cputime();
  ac = lams[1];
  fv = V + n;
  FV = create_vec(n);
  set_vec(fv, FV);
  save_vec(FV, argv[2]);
  printf("arpack took %10.8f (iter=%d), ac=%.3E\n", end-start, iter, ac);
  test_finish();
  return 0;
}

#include "fiedler_arpack.h"

/*
 * tests fiedler_arpack + cholesky fact/solve [cholmod library]
 */

cholmod_sparse * L;
cholmod_factor * factor;
cholmod_dense * FV;

void solve(double * x, double * y)
{
  lu_solve(L, factor, x, y);
}

int main(int argc, char * argv[])
{
  int n, iter, nev;
  double *lams, *V, ac, *fv, tol=1e-7, start, end, a;
  assert(argc == 3, "Usage: %s <lap_file> <fv_file>\n", argv[0]);
  test_start();
  L = load_mat(argv[1]);
  n = L->nrow;
  nev = 1;
  lams = alloc_double(nev);
  V = alloc_double(nev * n);
  start = cputime();
  factor = spec_upd_prec(L, &a);
  fiedler_arpack(n, "LM", nev, 0, tol, NULL, solve, lams, V, &iter);
  end = cputime();
  ac = lams[0] - a;
  fv = V;
  FV = create_vec(n);
  set_vec(fv, FV);
  save_vec(FV, argv[2]);
  printf("arpack took %10.8f (iter=%d), ac=%.3E\n", end-start, iter, ac);
  test_finish();
  return 0;
}

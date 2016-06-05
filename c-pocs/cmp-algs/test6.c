#include "fiedler_arpack.h"

/*
 * tests fiedler_arpack + cholesky fact/solve [cholmod library]
 */

cholmod_sparse * L;
cholmod_factor * factor;
cholmod_dense * b;
cholmod_dense * FV;
double solve_time, solvesl_time, solvecp_time;
int solve_calls;

void solve(double * x, double * y)
{
  double start, end;
  cholmod_dense * soln;
  solve_calls++;
  solvecp_time += take_time(set_vec(x, b));
  solvesl_time += take_time(soln = chol_solve(factor, b));
  solvecp_time += take_time(copy_vec(soln, y));
}

int main(int argc, char * argv[])
{
  int n, iter, nev;
  double *lams, *V, ac, *fv, tol=1e-7;
  double start, end, tot_start, tot_end, tot_time, a, sutime;
  assert(argc == 3, "Usage: %s <lap_file> <fv_file>\n", argv[0]);
  test_start();
  L = load_mat(argv[1]);
  n = L->nrow;
  nev = 1;
  lams = alloc_double(nev);
  V = alloc_double(nev * n);
  b = create_vec(n);
  solve_calls = solve_time = solvecp_time = solvesl_time = 0;
  tot_start = cputime();
  sutime = take_time(factor = spec_upd_prec(L, &a));
  fiedler_arpack(n, "LM", nev, 0, tol, NULL, solve, lams, V, &iter);
  tot_end = cputime();
  tot_time = tot_end - tot_start;
  solve_time = solvesl_time + solvecp_time;  
  ac = lams[0] - a;
  fv = V;
  FV = create_vec(n);
  set_vec(fv, FV);
  save_vec(FV, argv[2]);
  printf("arpack+su took %10.8f (iter=%d), ac=%.3E\n",
         tot_time, iter, ac);
  printf("spec upd stats: %.2f %10.8f\n",
         sutime/tot_time, sutime);
  printf("solve stats: %.2f, %d, %10.8f, %10.8f, %10.8f, %10.8f\n",
         solve_time / tot_time,
         solve_calls,
         solve_time,
         solve_time/(double) solve_calls,
         solvesl_time,
         solvecp_time);
  test_finish();
  return 0;
}

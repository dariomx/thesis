#include "test_util.h"

/*
 * tests chol fact/solver with a laplacian
 */

int main(int argc, char * argv[])
{
  double start, end;
  cholmod_sparse * L;
  cholmod_dense * B;
  cholmod_dense * X;
  cholmod_factor * factor;
  assert(argc == 4, "Usage: %s <L_file> <b_file> <x_file>\n", argv[0]);
  test_start();
  L = load_mat(argv[1]);
  start = cputime();
  diag_add(L, 1e-7);
  end = cputime();
  printf("diag_add took %10.8f\n", end-start);  
  B = load_vec(argv[2]);
  start = cputime();
  factor = chol_factor(L);
  X = chol_solve(factor, B);
  end = cputime();
  printf("chol fact/solve took %10.8f\n", end-start);
  save_vec(X, argv[3]);
  test_finish();
  return 0;
}

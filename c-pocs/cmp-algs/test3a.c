#include "test_util.h"
#include "test3a.h"

/*
 * tests LU fact/solver (constant arrays)
 */

int main(int argc, char * argv[])
{
  void *Symbolic, *Numeric ;
  double start, end;
  cholmod_dense * B;
  cholmod_dense * X;
  int n, ec;
  assert(argc == 3, "Usage: %s <b_file> <x_file>\n", argv[0]);
  test_start();
  n = 867;
  B = load_vec(argv[1]);
  X = create_vec(n);  
  start = cputime();
  ec = umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, null, null) ;
  assert (ec == 0, "failed symbolic factor %d", ec);
  ec = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
  assert (ec == 0, "failed numeric factor %d", ec);
  umfpack_di_free_symbolic (&Symbolic) ;
  ec = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, X->x, B->x, Numeric, null, null) ;
  assert (ec == 0, "failed solve %d", ec);  
  umfpack_di_free_numeric (&Numeric) ;
  end = cputime();
  printf("lu fact/solve took %10.8f\n", end-start);
  save_vec(X, argv[2]);
  test_finish();
  return 0;
}

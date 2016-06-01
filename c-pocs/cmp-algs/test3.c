#include "test_util.h"

/*
 * tests LU fact/solver
 */

int main(int argc, char * argv[])
{
  double start, end;
  cholmod_sparse * A;
  cholmod_dense * B;
  cholmod_dense * X;
  void * factor;
  int n;
  assert(argc == 4, "Usage: %s <A_file> <b_file> <x_file>\n", argv[0]);
  test_start();
  A = load_mat(argv[1]);
  {
    int i;
    Int * Ap , *Ai;
    double * Ax;
    Ap = A->p;
    Ai = A->i;
    Ax = A->x;
    eprint("\nA->p\n");
    for(i=0; i<A->ncol+1; i++) eprint("%d, ", Ap[i]);
    eprint("\nA->i\n");
    for(i=0; i<A->nzmax; i++) eprint("%d, ", Ai[i]);
    eprint("\nA->x\n");
    for(i=0; i<A->nzmax; i++) eprint("%.16f, ", Ax[i]);    
  }
  n = A->nrow;
  B = load_vec(argv[2]);
  X = create_vec(n);  
  start = cputime();
  factor = lu_factor(A);
  lu_solve(A, factor, (double *) B->x, (double *) X->x);
  end = cputime();
  printf("lu fact/solve took %10.8f\n", end-start);
  save_vec(X, argv[3]);
  test_finish();
  return 0;
}

#include "test_util.h"

cholmod_common CM;

void test_start()
{
  cholmod_start(&CM);
}

FILE * safe_open(char * fn, char * mode)
{
  FILE * f = fopen(fn, mode);
  assert(f != NULL, "Unable to file %s in mode %s\n", fn, mode);
  return f;
}

cholmod_sparse * load_mat(char * fn)
{
  FILE * f = safe_open(fn, "r");
  return cholmod_read_sparse(f, &CM);
  fclose(f);
}

void print_mat(cholmod_sparse * A, char * name)
{
  cholmod_print_sparse(A, name, &CM);
}

void save_mat(cholmod_sparse * A, char * fn)
{
  FILE * f = safe_open(fn, "w");
  cholmod_write_sparse(f, A, NULL, NULL, &CM);
}

/*
  Adds a positive sigma to L's diagonal, where L is assumed to be a
  weighted Laplacian (hence doing this will not alter the zero
  structure of the diagonal; which shall not have any zero actually);
  even isolated nodes will have at least the -self_similarity metric,
  which is assumed to have bigger magnitude than sigma.

  We also assume that assing sigma does not cause any overflow.
 */
void diag_add(cholmod_sparse * L, double sigma)
{
  int k,j;
  Int *Lp, *Li;
  double * Lx;
  Lp = L->p;
  Li = L->i;
  Lx = L->x;
  for(k=0, j=0; j < L->ncol; j++)
  {
    for(k=Lp[j]; Li[k] < j; k++);
    assert (Li[k] == j, "Could not find diagonal entry (%d,%d)", j, j);
    Lx[k] += sigma;
  }
}

cholmod_dense * load_vec(char * fn)
{
  FILE * f = safe_open(fn, "r");
  return cholmod_read_dense(f, &CM);
  fclose(f);
}

cholmod_dense * create_vec(int n)
{
  return cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, &CM);
}

void set_vec(double * x, cholmod_dense * v)
{
  memcpy(v->x, x, v->nrow * sizeof(double));
}

void copy_vec(cholmod_dense * v, double * x)
{
  memcpy(x, v->x, v->nrow * sizeof(double));
}

void print_vec(cholmod_dense * v, char * name)
{
  cholmod_print_dense(v, name, &CM);
}

void save_vec(cholmod_dense * v, char * fn)
{
  FILE * f = safe_open(fn, "w");
  cholmod_write_dense(f, v, NULL, &CM);
  fclose(f);
}

/* y = A * x */
void mult_matvec(cholmod_sparse * A, cholmod_dense * x, cholmod_dense * y)
{
  double one[2] = {1,0}, zero[2] = {0,0};
  cholmod_sdmult(A, 0, one, zero, x, y, &CM);
}

void * lu_factor(cholmod_sparse * A)
{
  void *symbolic, *numeric ;
  Int *Ap, *Ai;
  double * Ax;
  int n = A->nrow;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;  
  umfpack_di_symbolic(n, n, Ap, Ai, Ax, &symbolic, null, null) ;
  umfpack_di_numeric(Ap, Ai, Ax, symbolic, &numeric, null, null) ;
  umfpack_di_free_symbolic (&symbolic);
  return numeric;
}

void lu_solve(cholmod_sparse * A, void * factor, double * b, double * x)
{
  Int *Ap, *Ai;
  double * Ax;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;    
  umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, factor, null, null) ;  
}
                    
void test_finish()
{
  cholmod_finish(&CM);
}

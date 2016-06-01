#include "test_util.h"

cholmod_common CM;

cholmod_dense *X, *Y, *E;

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
  As a useful side calculation, it returns the minimum element of the
  diagonal (which for Laplacias is the minimum degree)
 */
double diag_add(cholmod_sparse * L, double sigma)
{
  int k,j;
  Int *Lp, *Li;
  double *Lx, min = L->ncol;
  Lp = L->p;
  Li = L->i;
  Lx = L->x;
  for(k=0, j=0; j < L->ncol; j++)
  {
    for(k=Lp[j]; Li[k] < j; k++);
    assert (Li[k] == j, "Could not find diagonal entry (%d,%d)", j, j);
    Lx[k] += sigma;
    if ( Lx[k] < min )
      min = Lx[k];
  }
  return min;
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

void scale_vec(cholmod_dense * v, double s)
{
  Int nz = MAX (1, v->nzmax);
  double *vx = v->x;
  int i;
  for(i=0; i<nz; i++)
    vx[i] *= s;
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

cholmod_factor * chol_factor(cholmod_sparse * A)
{
  cholmod_factor * factor = cholmod_analyze(A, &CM);
  int ec = cholmod_factorize(A, factor, &CM);
  assert(ec >= 0, "failed to analyze chol: %d\n", ec);
  if ( ec > 0 )
    eprint("warning in analyze chol: %d\n", ec);
  X = Y = E = NULL;
  return factor;
}

cholmod_dense * chol_solve(cholmod_factor * factor, cholmod_dense * b)
{
  int ec = cholmod_solve2 (CHOLMOD_A, factor, b, NULL, &X, NULL, &Y, &E, &CM) ;
  assert(ec >= 0, "failed to solve system using cholesky: %d\n", ec);
  if ( ec > 0 )
    eprint("warning in solve chol: %d\n", ec);
  return X;
}

/* spectral update preconditioning of the laplacian */
cholmod_factor * spec_upd_prec(cholmod_sparse * L, double * a_out)
{
  int ec;
  double n = (double) L->nrow;
  double a = 1e-2;
  double min_degree = diag_add(L, a);
  cholmod_factor * factor = chol_factor(L);  
  double ac_bound =  n/(n-1) * min_degree;
  double b = 1;
  double c = (ac_bound + b) / n;
  cholmod_sparse * v1 = cholmod_ones(L->nrow, 1, CHOLMODE_REAL, &CM);
  scale_vec(v1, c);
  ec = cholmod_updown(TRUE, v1, factor, &CM);
  assert(ec == 0, "error in spec_upd_prec: %d\n", ec);
  *a_out = a;
  return factor;
}

void test_finish()
{
  cholmod_finish(&CM);
}

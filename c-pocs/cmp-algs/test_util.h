#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cholmod_internal.h"
#include "cholmod.h"
#include "umfpack.h"

#define cputime SuiteSparse_time

#define TRUE 1

#define null ((double *) NULL)

#define eprint(fmt, args...) fprintf(stderr, fmt, ##args)

#define assert(cond, msg, args...) if (!(cond)) {eprint(msg , ##args); exit(1);}

#define alloc_double(n) (double *) malloc(n * sizeof(double))

#define take_time(s) (start=cputime(), s, end=cputime(), end-start)

void test_start();

cholmod_sparse * load_mat(char * fn);

void print_mat(cholmod_sparse * A, char * fn);

void save_mat(cholmod_sparse * A, char * fn);

double diag_add(cholmod_sparse * L, double sig);

cholmod_dense * load_vec(char * fn);

cholmod_dense * create_vec(int n);

void set_vec(double * x, cholmod_dense * v);

void copy_vec(cholmod_dense * v, double * x);

void print_vec(cholmod_dense * v, char * name);

void save_vec(cholmod_dense * v, char * fn);

void scale_vec(cholmod_dense * v, double s);

void mult_matvec(cholmod_sparse * A, cholmod_dense * x, cholmod_dense * y);

void * lu_factor(cholmod_sparse * A);

void lu_solve(cholmod_sparse * A, void * factor, double * b, double * x);

cholmod_factor * chol_factor(cholmod_sparse * A, double a);

cholmod_dense * chol_solve(cholmod_factor * factor, cholmod_dense * b);

cholmod_factor * spec_upd_prec(cholmod_sparse * L, double * a_out);

void test_finish();

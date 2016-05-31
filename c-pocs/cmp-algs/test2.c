#include "test_util.h"

/*
 * tests basic vector i/o functions and set_vec/copy_vec
 */
int main(int argc, char * argv[])
{
  int n, i;
  double *x, *y;
  assert(argc == 3, "Usage: %s <n> <out_file>\n", argv[0]);
  test_start();
  n = atoi(argv[1]);
  x = alloc_double(n);
  y = alloc_double(n);
  for(i=0; i<n; i++)
    x[i] = (double) i;
  cholmod_dense * v = create_vec(n);
  set_vec(x, v);
  print_vec(v, "v");
  save_vec(v, argv[2]);
  copy_vec(v, y);
  for(i=0; i<n; i++)
    printf("y[%d] = %f\n", i, y[i]);
  test_finish();
  return 0;
}

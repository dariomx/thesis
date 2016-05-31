#include "test_util.h"

/*
 * tests basic matrix i/o functions and diag_add
 */
int main(int argc, char * argv[])
{
  assert(argc == 4, "Usage: %s <lap_file> <sigma> <out_file>\n", argv[0]);
  test_start();
  cholmod_sparse * L = load_mat(argv[1]);
  double sig = atof(argv[2]);
  diag_add(L, sig);
  print_mat(L, "L");
  save_mat(L, argv[3]);
  test_finish();
  return 0;
}

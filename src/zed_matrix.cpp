#include "tetradiplomix.h"

//' Return a 6 x 4 matrix whose rows values of the zeds
//'
//' blah blah
//'
//' @keywords internal
//'
//' @return  an integer matrix of 0s and 1s
//'
//' @examples
//' example(tcf2param_list)
//' ale_glL <- geno_logL(ale_par_list)
//' @export
// [[Rcpp::export]]
IntegerMatrix zed_matrix() {
  IntegerMatrix ret(6, 4);

  ret(0, 2) = ret(0, 3) = 1L;
  ret(1, 1) = ret(1, 3) = 1L;
  ret(2, 1) = ret(2, 2) = 1L;
  ret(3, 0) = ret(3, 3) = 1L;
  ret(4, 0) = ret(4, 2) = 1L;
  ret(5, 0) = ret(5, 1) = 1L;

  return(ret);
}

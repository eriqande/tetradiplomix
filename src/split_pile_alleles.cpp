#include "tetradiplomix.h"

//' Count the number of alleles in G that are currently allocated to 0 or 1 according to Z
//'
//' blah blah
//'
//' @param G the N x 4 matrix of allele present in the N individuals.
//' @param Z a matrix of 0s and 1s parallel to G.
//' @param out a 2 x A NumericMatrix that counts the number of each allele allocated to 0 (first row)
//' and 1 (second row).  This is numeric so we can add a prior to it. A is the number of alleles.
//' You have to allocate memory to this outside the function.  The number of columns is important.
//'
//' @keywords internal
//'
//' @return  out is an output variable
//'
//' @examples
//' example(tcf2param_list)
//' ale_glL <- geno_logL(ale_par_list)
//' @export
// [[Rcpp::export]]
void split_pile_alleles(IntegerMatrix G, IntegerMatrix Z, NumericMatrix out) {
  int i,j;
  int A = out.ncol();
  int N = G.nrow();

  // initialize to prior to accumulate a sum on top of it
  for(i=0;i<A;i++) {
    out(0, i) = 1.0 / A;
    out(1, i) = 1.0 / A;
  }

  for(i=0;i<N;i++) {
    for(j=0;j<4;j++) {
      out(Z(i,j), G(i, j)) += 1.0;
    }
  }
}

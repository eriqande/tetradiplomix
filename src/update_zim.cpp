#include "tetradiplomix.h"


//' update the ZIM matrix
//'
//' This updates the allocations for each individuals of the gene copies
//' to locus 0 or locus 1.  It computes the full conditional for each of the
//' configurations in ZM, and then it samples from one of those and sets those
//' values in the row of ZIM.  This only works when called from with RCpp.
//' Otherwise ZIM does not get updated.
//' @export
// [[Rcpp::export]]
void update_zim(IntegerMatrix G,
                IntegerMatrix ZIM,
                IntegerMatrix ZM,
                NumericMatrix p) {

  int N = G.nrow();
  int C = G.ncol();  // will always be four for tetraploid...
  int nZ = ZM.nrow();  // the number of different zed configurations
  NumericVector probs(nZ);  // for storing the full conditionals
  int i,j,k;
  int zchosen;

  for(i=0;i<N;i++) {
    for(j=0;j<nZ;j++) {
      probs[j] = 1.0;  // to accumulate a product
      for(k=0;k<C;k++) {
        probs[j] *= p(ZM(j,k), G(i,k));
      }
    }
    probs = probs / sum(probs);

    // now, probs is the full conditional probilities and we want to
    // sample a zero-based index from that
    zchosen = samp_from_vec(probs);

    // now tranfer those values back to row i in ZIM
    for(k=0;k<C;k++) {
      ZIM(i,k) = ZM(zchosen, k);
    }

  }
}

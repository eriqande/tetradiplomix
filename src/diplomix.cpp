#include "tetradiplomix.h"


//' basic function to do mcmc on the model at a single locus
//'
//' more later
//'
//' @param G an N x 4 matrix (N is the number of individuals) giving
//' the inferred genotypes of each individual.  The alleles/haplotypes
//' are coded as integers from 0 up.  For example, an individual might
//' have the genotype 0031---two doses of the 0 allele, and one each of the
//' 1 and 3 alleles.  Note that the total number of alleles at the locus is
//' taken to be one greater than the max entry in the genos matrix.  There
//' can be no missing data in this matrix.
//'
//' @export
// [[Rcpp::export]]
List diplomix(IntegerMatrix G, int reps) {
  int A = max(G) + 1;  // the total number of distinct alleles
  int N = G.nrow();
  IntegerMatrix ZM;  // to hold the zed matrix
  IntegerMatrix ZIM(N, 4); // an N x 4 matrix that holds the actual Z_i's for each individual
                              // Note that we initialize them all to 0 for the first step
  NumericVector p0(A);
  NumericVector p1(A);
  NumericMatrix acnts(2, A);  // for counting up the number of alleles allocated to each locus
  NumericMatrix p(2, A); // for storing the current allele freqs
  List ret;
  int r;
  IntegerMatrix dump;
  List pout(reps), zimout(reps);  // preallocate space here

  // Get this for repeated use
  ZM = zed_matrix();


  // First get allele freqs with every allele allocated to 0
  split_pile_alleles(G, ZIM, acnts);
  p0 = acnts(0, _);

  // then simulate p1 and p0 to be a drifted version from that
  p0 = p0 / sum(p0);

  p1 = dirch_from_params(p0 * 20);  // 20 corresponds to drift of about F = 0.05
  p0 = dirch_from_params(p0 * 20);

  p(0, _) = p0;
  p(1, _) = p1;

  // Now that p0 and p1 are initialized, we do our reps
  for(r=0;r<reps;r++) {
    update_zim(G, ZIM, ZM, p);  // Gibbs sampling for the gene copy allocations

    // now, given those Zeds we update the allele frequencies via Gibbs sampling
    split_pile_alleles(G, ZIM, acnts);
    p(0, _) = dirch_from_params(acnts(0, _));
    p(1, _) = dirch_from_params(acnts(1, _));

    pout[r] = clone(p);
    zimout[r] = clone(ZIM);
  }


  ret = List::create(
    Named("N") = N,
    Named("A") = A,
    Named("ZM") = ZM,
    Named("ZIM") = ZIM,
    Named("p") = p,
    Named("acnts") = acnts,
    Named("ptrace") = pout,
    Named("zimtrace") = zimout
  );

  return(ret);
}

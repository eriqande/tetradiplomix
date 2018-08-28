#include <Rcpp.h>
using namespace Rcpp;




// function prototypes
IntegerMatrix zed_matrix();
void split_pile_alleles(IntegerMatrix G, IntegerMatrix Z, NumericMatrix out);
List diplomix(IntegerMatrix G, int reps);

void update_zim(IntegerMatrix G, IntegerMatrix ZIM, IntegerMatrix ZM, NumericMatrix p);

// in rcpp_utilities.cpp
double rgammadouble(int a, double b, double c);
NumericVector dirch_from_params(NumericVector C);
int samp_from_vec(NumericVector V);

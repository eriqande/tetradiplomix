# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Given a vector of parameters for different categories in 1...n,
#' simulate a Dirichlet random vector
#'
#' Takes a vector of counts for 1:n collections,
#' and returns a Dirichlet random variable generated by adding the prior to each
#' collection value, and simulating an alpha from a gamma distribution
#' with this shape parameter.
#'
#' The categories are labeled in C from 1 up to n.  n is the length of \code{lambda},
#' which is a vector of priors. Note that all elements of \code{lambda}
#' must be strictly greater than 0.
#' @keywords internal
#' @param C  a vector giving counts of categories
#' @param lambda priors for the categories
#' @export
dirch_from_params <- function(C) {
    .Call('_tetradiplomix_dirch_from_params', PACKAGE = 'tetradiplomix', C)
}

#' sample from 0:(n-1) given n probabilities
#'
#' boing
#' @export
#' @keywords internal
samp_from_vec <- function(V) {
    .Call('_tetradiplomix_samp_from_vec', PACKAGE = 'tetradiplomix', V)
}

#' basic function to do mcmc on the model at a single locus
#'
#' more later
#'
#' @param G an N x 4 matrix (N is the number of individuals) giving
#' the inferred genotypes of each individual.  The alleles/haplotypes
#' are coded as integers from 0 up.  For example, an individual might
#' have the genotype 0031---two doses of the 0 allele, and one each of the
#' 1 and 3 alleles.  Note that the total number of alleles at the locus is
#' taken to be one greater than the max entry in the genos matrix.  There
#' can be no missing data in this matrix.
#'
#' @export
diplomix <- function(G, reps) {
    .Call('_tetradiplomix_diplomix', PACKAGE = 'tetradiplomix', G, reps)
}

#' Count the number of alleles in G that are currently allocated to 0 or 1 according to Z
#'
#' blah blah
#'
#' @param G the N x 4 matrix of allele present in the N individuals.
#' @param Z a matrix of 0s and 1s parallel to G.
#' @param out a 2 x A NumericMatrix that counts the number of each allele allocated to 0 (first row)
#' and 1 (second row).  This is numeric so we can add a prior to it. A is the number of alleles.
#' You have to allocate memory to this outside the function.  The number of columns is important.
#'
#' @keywords internal
#'
#' @return  out is an output variable
#'
#' @examples
#' example(tcf2param_list)
#' ale_glL <- geno_logL(ale_par_list)
#' @export
split_pile_alleles <- function(G, Z, out) {
    invisible(.Call('_tetradiplomix_split_pile_alleles', PACKAGE = 'tetradiplomix', G, Z, out))
}

#' update the ZIM matrix
#'
#' This updates the allocations for each individuals of the gene copies
#' to locus 0 or locus 1.  It computes the full conditional for each of the
#' configurations in ZM, and then it samples from one of those and sets those
#' values in the row of ZIM.  This only works when called from with RCpp.
#' Otherwise ZIM does not get updated.
#' @export
update_zim <- function(G, ZIM, ZM, p) {
    invisible(.Call('_tetradiplomix_update_zim', PACKAGE = 'tetradiplomix', G, ZIM, ZM, p))
}

#' Return a 6 x 4 matrix whose rows values of the zeds
#'
#' blah blah
#'
#' @keywords internal
#'
#' @return  an integer matrix of 0s and 1s
#'
#' @examples
#' example(tcf2param_list)
#' ale_glL <- geno_logL(ale_par_list)
#' @export
zed_matrix <- function() {
    .Call('_tetradiplomix_zed_matrix', PACKAGE = 'tetradiplomix')
}


#' compute Fst from the allele frequencies at each rep
#'
#' To be applied to the output of diplomix.  This uses Nei 1973's expression
#' as grabbed from https://biology.stackexchange.com/questions/40756/calculating-pairwise-fst-from-allele-frequencies
#' @param p a list of matrices, each with 2 rows and A columns where A
#' is the number of alleles.
#' @return A vector of Fst values.
#' @export
fst_from_ptrace <- function(p) {
  tmp <- lapply(p, function(x) {

    HS <- mean(x^2)
    HT <- mean(colMeans(x) ^ 2)

   (HT - HS) / HT
  })

  unlist(tmp)
}

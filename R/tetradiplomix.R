#' tetradiplomix: inference of diploidizing regions of tetraploid genomes
#'
#'
#' @docType package
#' @name tetradiplomix
#' @importFrom Rcpp evalCpp
#' @useDynLib tetradiplomix
NULL


# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if(getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      "."
    )
  )
}

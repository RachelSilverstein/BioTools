#' Function Template
#'
#' Description of the function
#'
#' A more detailed description of what the function does.
#'
#' @param wt_sequnece Character string representing the RNA sequence to calculate RNA structure for (DNA also ok... will substitute U for T)
#' @param variant The variant
#'
#' @return The sum of x and y
#' @export
#'
#' @seealso \code{\link{another_function}} Info
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#' @examples
#'
#' example_function(2,3)
#' # [1] 5
#'
#'
#'
#'

varFold <- function(wt_sequence, variant = FALSE, start = NULL, stop = NULL, wt = NULL, mut = NULL) {
  if (!require(devtools)) {
    install.packages("devtools")
  }
  library(devtools)

  if (!existsFunction('mutateSeq')) {
    install_github("RachelSilverstein/BioTools")
  }
  library(BioTools)

  if (!variant) {
    sequence <- wt_sequence
  } else {
    if (is.null(start) | is.null(stop) | is.null(wt) | is.null(mut)) {
      stop("Variant details missing.")
    } else {
      sequence <- mutate(wt_sequence, start =  start, stop = stop, wt = wt, mut = mut)
    }
  }
}

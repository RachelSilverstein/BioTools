#' MutateSeq
#'
#' Mutate a DNA or RNA sequence
#'
#' Error if the wild type sequence at the postions provided does not match the wt seqment that you are trying to change.
#' Error if the mutant sequence provided does not match the length of the start and stop positions.
#'
#'
#' @param wt_sequence Character string representing the wild type DNA or RNA sequence to be mutated
#' @param wt Optional character string, the part of the sequence to be replaced by the mutant sequence. Included only for checking purposes to make sure the position to mutate is correct.
#' @param mut Character vector of mutant sequence that the indicated positions should be changed to.
#' @param start Integer nucleotide position that the mutation starts
#' @param stop Integer nucleotide position that the mutation ends (last mutant nucleotide)
#'
#' @return Character string of the mutant sequence.
#' @export
#'
#' @seealso \code{\link{another_function}} Info
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#' @examples
#'
#' mutateSeq("AAACATGGG", start = 4, stop = 6, wt = 'CAT', mut = 'GGC')
#' # [1] "AAAGGCGGG"
#'
#'

mutateSeq <- function(wt_sequence, start, stop, mut, wt = NULL) {
  sequence <- unlist(strsplit(wt_sequence, split = ""))
  if (!is.null(wt)) {
    wt <- unlist(strsplit(wt, split = ""))
    if (!all(sequence[c(start:stop)] == wt))
      stop("The sequence you are trying to mutate from does not match the wild type sequence at the given position.")
  }
  mut <- unlist(strsplit(mut, split = ""))
  if (!(length(mut) == stop - start + 1)) {
    stop("The length of the mutant sequence provided does not match start and stop positions.")
  }
  sequence[seq(start, stop)] <- mut
  sequence <- paste(sequence, sep = "", collapse = "")
  return(sequence)
}

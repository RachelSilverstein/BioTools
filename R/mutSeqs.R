#' Mutate Sequences
#'
#' Mutate a vector of DNA or RNA sequences; vectorized.
#'
#' Vectorized function: The ith sequences in the sequences vector will be mutated at the ith position in the start and stop vectors, and changed to the ith mutation provided in the mut vector.
#' Error if the wild type sequences at the postions provided does not match the wt seqment that you are trying to change.
#' Error if the mutant sequences provided does not match the length of the start and stop positions.
#'
#'
#' @param sequences Character vector representing the wild type DNA or RNA sequences to be mutated. Each sequences in a separate string.
#' @param wt Optional character vector, the part of the sequences to be replaced by the mutant sequences. In the same order as wt_sequences. Included only for checking purposes to make sure the position to mutate is correct.
#' @param mut Character vector of mutant sequences that the indicated positions should be changed to.
#' @param start Integer vector, ith elements represents the nucleotide positions that the ith mutation starts.
#' @param stop Integer vector, ith element represents the nucleotide position that the ith mutation ends (last mutant nucleotide)
#'
#' @return Character vector of the mutant sequences.
#' @export
#'
#'
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#' @examples
#'
#' mutSeqs(sequences = c('AAC', "GCC"), start = c(1,2), stop = c(2, 3), wt = c("AA", "CC"), mut = c('GG', 'AA'))
#'
#

mutSeqs <- function(sequences, start, stop, mut, wt = NULL) {
  # check that all the input vectors are the same length:
  l <- length(sequences)
  stopifnot((length(start) == l) & (length(stop) == l) & (length(mut) == l))

  sequences <- strsplit(sequences, split = "")
  # check that the wt provided matches the correct position in the sequences provided
  if (!is.null(wt)) { # if you provided a wt vector
    stopifnot(length(wt) == l)
    wt <- strsplit(wt, split = "")
    for (i in seq_along(sequences)) {
      cat("wt codon ", wt[[i]])
      cat("from sequence ", sequences[[i]][start[i]:stop[i]])
      stopifnot(sequences[[i]][start[i]:stop[i]] == wt[[i]])
    }
  }
  # check that the start and stop vectors match the length of the mutations provided
  mut <- strsplit(mut, split = "")
  for (i in seq_along(mut)) {
    stopifnot(length(mut[[i]]) == stop[i] - start[i] + 1)
  }

  # make the mutations
  for (i in seq_along(sequences)) {
    sequences[[i]][start[i]:stop[i]] <- mut[[i]]
  }
  f <- function (x) paste(x, collapse = "")
  sequences <- lapply(sequences, f)
  return(unlist(sequences))
}

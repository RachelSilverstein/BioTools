#' Variant RNAfold
#'
#' Get minimum free energy and ensemble free energy for specified variants of RNA sequences and their differences from wild type (calculated by \href{https://www.tbi.univie.ac.at/RNA/}{RNAfold}); vectorized.
#'
#' A more detailed description of what the function does.
#'
#' @param wt_sequences Character vector representing the wild type RNA sequences to (DNA also ok... will substitute U for T)
#' @param variants Optional logical. Indicates whether variants are provided that should be analyzed (TRUE) or whether the wild type sequence should be analyzed (FALSE). Default is FALSE.
#' @param start Optional integer vector same length as wt_sequences. Specifies the starting nucleotide positions of the variants. Must be provided if variants is TRUE.
#' @param stop Optional integer vector same length as wt_sequences. Specifies the ending nucleotide positions of the variants. Must be provided if variants is TRUE.
#' @param wt Optional character vector same length as wt_sequences. Not necessary to provide even if variant = TRUE. For error-checking purposes only.
#' @param mut Optional character vector same length as wt_sequences. Specifies what the sequence between start and stop nucleotides should be changed to. Must be provided if variants is TRUE.
#' @param MFEwt Optional numeric vector, the minimum free energies of the wild type sequences. This parameter is provided rather than calculated to avoid repetitive MFE calculations for vectors where all of the wild type sequences are identical.
#' @param ensFEwt Optional numeric vector, the ensemble free energies of the wild type sequences. This parameter is provided rather than calculated to avoid repetitive MFE calculations for vectors where all of the wild type sequences are identical.
#'
#' @return Data frame with 4 columns, MFE (minimum free energy), MFE_diff (minimum free energy difference, mutant minus wild type), ensFE (ensemble free energy), ensFE_diff (ensemble free energy difference, mutant minus wild type)
#' @export
#'
#' @seealso \code{\link{RNAfold}} Called by this function
#'
#' @seealso \code{\link{mutSeqs}} Called by this function
#'
#'
#' @examples
#'
#' # first get the wild type free energes using varfold without a variant
#' seq1 <- "AAAAAAGUUUGGGGGGCCCCCCTTTTTT"
#' seq2 <- "UUUUUUUUAAAAAAAUUUUUUUUAAAAAAAAUUUUUUUUAAAAAAAA"
#' wt_data <- varFold(wt_sequences = c(seq1, seq2), variants = FALSE)
#' MFEwt <- wt_data$MFE
#' ensFEwt <- wt_data$ensFE
#' # > wt_data
#' #     MFE                 MFE_diff  ensFE                 ensFE_diff
#' # 1  -7.4 No WT MFE data provided.  -8.35 No WT ensFE data provided.
#' # 2 -10.0 No WT MFE data provided. -11.54 No WT ensFE data provided.
#'
#' # now get the variant information with differences using the wild info just calculated
#' var_data <- varFold(wt_sequences = c(seq1, seq2), variants = TRUE, start = c(1, 2), stop = c(2, 3), wt = c("AA", "UU"), mut = c("GG", "CC"), MFEwt = MFEwt, ensFEwt = ensFEwt)
#' # > var_data
#' # MFE MFE_diff ensFE ensFE_diff
#' # MFE MFE_diff ensFE ensFE_diff
#' # 1 -7.6     -0.2 -8.49      -0.14
#' # 2 -8.2      1.8 -9.42       2.12
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#'

varFold <- function(wt_sequences,
                    variants = FALSE,
                    start = NULL,
                    stop = NULL,
                    wt = NULL,
                    mut = NULL,
                    MFEwt = NULL,
                    ensFEwt = NULL) {
  # chack for and install dependencies
  if (!require(devtools)) {
    install.packages("devtools")
  }
  library(devtools)

  if (!existsFunction('mutSeqs')) {
    install_github("RachelSilverstein/BioTools")
  }
  library(BioTools)

  if (!variants) { # keep the sequences as the wild type sequences
    sequences <- wt_sequences
  } else { # mutate them
    sequences <- BioTools::mutSeqs(wt_sequences, start =  start, stop = stop, wt = wt, mut = mut)
  }

  # Get the output of RNAfold for the sequences
  # the first item in () these brackets is the Minimum free energy
  # The second item in [] these brackets is the free energy of the ensemble
  raw_output <- BioTools::RNAfold(sequences)

  # parse the raw output
  library(stringr)
  # Get the parenthesis and what is inside
  MFE <- str_extract_all(raw_output, "\\([^()]+\\)$")
  MFE <- unlist(MFE)
  MFE <- gsub(pattern = "[()]", replacement = "", MFE) # remove parentheses
  MFE <- as.numeric(MFE) # convert to numeric

  # do the same thing for the free energy of ensemble
  ensFE <- str_extract_all(raw_output, "\\[[^()]+\\]$")
  ensFE <- unlist(ensFE)
  ensFE <- gsub(pattern = "\\[", replacement = "", ensFE)
  ensFE <- gsub(pattern = "\\]", replacement = "", ensFE)
  ensFE <- as.numeric(ensFE)

  # calculate the free energy differences from wild type if provided
  if (is.null(MFEwt)) {
    MFE_diff <- "No WT MFE data provided."
  } else {
    stopifnot(length(MFEwt) == length(MFE))
    MFE_diff <-  MFE - MFEwt
  }
  if (is.null(ensFEwt)) {
    ensFE_diff <- "No WT ensFE data provided."
  } else {
    stopifnot(length(ensFEwt) == length(ensFE))
    ensFE_diff <-  ensFE - ensFEwt
  }

  result <- data.frame(MFE, MFE_diff, ensFE, ensFE_diff)
  return(result)
}

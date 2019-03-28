#' Invoke RNAfold
#'
#' Input RNA sequences to be analyzed by RNAfold.
#'
#' Pass a character vector of RNA sequences to be analyzed (and optional additional parameters for RNA fold). Returns the output of RNA fold in string format.
#' Make sure to always keep the option --noPS. This funciton will also create dot.ps in the project directory (a dot plot). It will only create a single file if multiple inputs are supplied and continuously overwrite it. I don't want it to write this file but I can't seem to figure out how to turn it off, so the funcion jsut deletes this file before returning.
#' Note, RNAfold must be saved in /usr/local/bin/.
#'
#' @param RNAs Character vector of RNA sequences to be analyzed. A DNA sequence is also ok... T will be converted to U.
#' @param options Additional parameters to pass to RNAfold. Default is "--noPS -p --jobs='0'". (Calculates partition function and MFE structure). Warning: always use --noPS!
#'
#' @return The output of RNAfold as a character string.
#' @export
#'
#'
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#' @examples
#' RNAfold(c("AAAAUAUUUUUUUAUAUAUUUUUUUAAAAAAA", "CCCCCTTTTTCTATTTTAAAAAAAAUUUUUU"))
#'
#'
RNAfold <- function(RNAs, options = "--noPS -p --jobs='0'") {
  # Note that the file that is used for input for RNAfold must be in the
  # current working directory as RNAfold does not allow you to input an
  # entire file path, only a file name.

  # Create a temporary file in the system's working directory
  working_dir <- system("pwd", intern = TRUE)
  tf <- tempfile(pattern = "temp", tmpdir = working_dir)
  tf_name <- gsub(pattern = "^(.*/)", replacement = "", x = tf, perl = TRUE)

  # write the input for RNA fold to the temporary file
  # RNAfold takes a simple sequence file where each sequence is on a new line
  writeLines(RNAs, con = tf)

  # send command to RNA fold
  command <- paste(c("RNAfold ", options, " --infile='", tf_name, "'"), sep="", collapse = "")
  output <- system(command, intern = TRUE)

  # remove the temporary file
  file.remove(tf)
  # remove the dotplot file that RNAfold writes
  dot_file <- paste(c(working_dir, "/dot.ps"), collapse = "")
  file.remove(dot_file)
  return(output)
}

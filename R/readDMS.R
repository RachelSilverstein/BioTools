#' Read DMS Raw Data
#'
#' Read raw data from DMS experiment from a tsv file and format into a data frame.
#'
#' A more detailed description of what the function does.
#'
#' @param path Character, single path to the tsv file containing raw DMS data.
#' @param name Character, name of the gene.
#' @return Data frame containing the DMS data with the following headers: wt_aa pos mut_aa wt_codon mut_codon annotation   nonselect1   nonselect2
#' @export
#'
#' @seealso \code{\link{another_function}} Info
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#' @examples
#' # this example only works on my conputer
#' setwd("/Users/Rachel/desktop/Roth Lab/R_DMS")
#' path <- './data/rawData_UBE2I_solid.txt'
#' gene_name <- 'UBE2I'
#' readDMS(path, gene_name)
#'

readDMS <- function(path, name) {
  raw_data <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  raw_data$gene <- rep(name, length(raw_data$wt_aa))
  return(raw_data)
}

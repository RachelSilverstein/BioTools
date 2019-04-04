
#' Get SNPs
#'
#' @param data  data frame of mutations with a mut_codon and wt_codon colums
#' @raturn data frame plus 3 new columns, snp_pos, snp_from, snp_to
#'
#' @export

getSNPs <- function(data) {

  # take in a data frame of mutations with a mut_codon and wt_codon colums
  # return only rows that are SNPs and add 3 columns:
  # snp_pos
  # snp_from
  # snp_to
  SNPs <- data.frame()

  for (i in seq_along(data$fitness_score)) {
    row <- data[i, ]
    wt_codon <- unlist(strsplit(row$wt_codon, split = ""))
    mut_codon <- unlist(strsplit(row$mut_codon, split = ""))
    mismatches <- 0
    for (j in 1:3) {
      if (wt_codon[j] == mut_codon[j]) {
        NULL
      } else {
        pos <- j
        mismatches <- mismatches + 1
      }
    }
    if (mismatches == 1) {
      snp_pos <- pos
      snp_from <- wt_codon[pos]
      snp_to <- mut_codon[pos]
      row <- cbind(row, snp_from, snp_to, snp_pos)
      SNPs <- rbind(SNPs, row)
    }
  }
  return(SNPs)
}

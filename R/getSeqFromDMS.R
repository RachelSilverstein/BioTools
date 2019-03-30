


getSeqFromDMS <- function (DMS_data) {
  sequence <- character(length = max(DMS_data$pos))
  for (i in seq_along(DMS_data$wt_codon)) {
    pos <- DMS_data[i, "pos"]
    wt_codon <- DMS_data[i, "wt_codon"]
    sequence[pos] <- wt_codon
  }
  sequence[1] <- "ATG" # since this is not usually measured by DMS
  sequence <- paste(sequence, collapse = "")
  return(sequence)
}

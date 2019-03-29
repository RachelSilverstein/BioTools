#' Calculate Fitness Scores
#'
#' Calculate fitness scores of variants from raw DMS data.
#'
#' A more detailed description of what the function does.
#'
#' @param DMS_data Data frame FOR A SINGLE GENE of the form output by \code{readDMS}
#' @param nonselect_cutoff Filter out low confidence variants with less than this many reads per million in nonselective conditions. (Default is 200 reads per million)
#'
#' @return Input data frame plus the following columns: select (the average of select1 and select2), nonselect (the average of nonselect1 and nonselect2), fold_change, fitness_score
#' @export
#'
#' @seealso \code{\link{readDMS}} The output of this function should be used as input for calculateFitnessScores.
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#' @examples
#' # this example only works on my computer
#' setwd("/Users/Rachel/desktop/Roth Lab/R_DMS")
#' path <- './data/rawData_UBE2I_solid.txt'
#' gene_name <- 'UBE2I'
#' raw_data <- readDMS(path, gene_name)
#' head(calculateFitnessScores(raw_data))

calculateFitnessScores <- function(DMS_data, nonselect_cutoff = 200) {
  # Calclulates the FS of each mutation.
  # Fitness scores are calculated separately for each gene. MAKE SURE THAT
  # DMS_data ONLY CONTAINS DATA FOR A SINGLE GENE.
  # Filter out low nonselect mutations and add pseudocounts if the select is
  # less than 0.
  #
  # Parameters:
  #     DMS_data      data frame      contains the raw DMS data
  #                                   as output of read_in_raw_data()
  #                       !!!!!!!!!FOR ONLY ONE GENE AT A TIME!!!!!!!!!!
  #     nonselect_cutoff    numeric   filter out mutations with nonselect
  #                         lower than this number (in reads per million)
  # Value:
  #     DMS_data   data frame     the input DMS data with new columns
  # ---------------------------------------------------------------------
  # Check that data is for only one gene at a time. Throw error if not.
  if(!all(DMS_data$gene == DMS_data$gene[1])) {
    print("calculate_fitness_score() should be called for only one gene
          at a time!")
    stop()
  }
  # avarage the 2 trials together and subtract the control (reads without mutagenesis)
  # as the control represents the representation of that mutation due to sequening error
  DMS_data$select <-  ((DMS_data$select1-DMS_data$control1) + (DMS_data$select2-DMS_data$control2))/2
  DMS_data$nonselect <-  ((DMS_data$nonselect1-DMS_data$control1) + (DMS_data$nonselect2-DMS_data$control2))/2
  # add a pseudocount if select or nonselect is 0 to avoid dividing by 0 or taking log(0)
  DMS_data$select <- pmax(DMS_data$select, 1)
  DMS_data$nonselect <- pmax(DMS_data$nonselect, 1)
  # remove rows where nonselect is less than nonselect_cutoff
  keep <- DMS_data$nonselect > nonselect_cutoff
  DMS_data <- DMS_data[keep,]
  # calculate the fold change:
  DMS_data$fold_change = DMS_data$select/DMS_data$nonselect
  # calculate the mean fold change for stop and syn mutations:
  stops <- DMS_data[DMS_data$annotation == 'STOP',]
  # don't count the stops in the last 1/10 of the ORF in the median stop fold change
  # since these often have less effect
  stops <- stops[stops$annotation != 'STOP' | stops$pos < max(stops$pos)*(9/10),]
  syns <- DMS_data[DMS_data$annotation == 'SYN',]
  med_FC_stops <- median(stops$fold_change)
  med_FC_syns <- median(syns$fold_change)
  # calculate the fitness scores
  DMS_data$fitness_score <- (log(DMS_data$fold_change) - log(med_FC_stops))/
    (log(med_FC_syns) - log(med_FC_stops))
  return(DMS_data)
}

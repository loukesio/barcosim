#' Create barcode sequences from parent sequences
#' @description
#' The r_gpseq function allows you to replicate the parent DNA sequences, that have been
#' created from the gpseq function.
#'
#'
#' @param dna_seq A character, output from gpseq
#' @param num_replicates A numeric value
#' @param error_rate A numeric value from 0 to 1
#'
#' @return a data.frame
#'
#' @examples
#' dna_seq <- c("AAGA","AATC")
#' r_gpseq(dna_seq,3,0.1) #create 3 offspring sequences with random substitution
#'
#' @export
r_gpseq <- function(dna_seq, num_replicates, error_rate) {

  # Validate input: dna_seq
  if (!is.character(dna_seq) || any(!nzchar(dna_seq))) {
    stop("Error: dna_seq should be a non-empty string or vector of DNA sequences.")
  }

  # Validate input: num_replicates
  if (!is.integer(num_replicates) || num_replicates <= 0) {
    if (is.numeric(num_replicates) && num_replicates %% 1 == 0 && num_replicates > 0) {
      num_replicates <- as.integer(num_replicates)
    } else {
      stop("Error: num_replicates should be a positive integer.")
    }
  }

  # Validate input: error_rate
  if (!is.numeric(error_rate) || error_rate < 0 || error_rate > 1) {
    stop("Error: error_rate should be a number between 0 and 1.")
  }

  # Create a vector to store the replicated sequences
  replicated_dna <- vector(mode = "list", length = num_replicates * length(dna_seq))

  # Loop over the input DNA sequences
  for (i in seq_along(dna_seq)) {
    # Replicate the sequence num_replicates times
    for (j in seq_len(num_replicates)) {
      # Initialize a new sequence
      new_seq <- ""

      # Loop over each base in the original sequence
      for (k in seq_len(nchar(dna_seq[[i]]))) {
        # Generate a random number between 0 and 1
        rand_num <- stats::runif(1)

        # Check if the random number is less than the error rate
        if (rand_num < error_rate) {
          # If so, generate a random base to substitute
          new_base <- sample(c("A", "C", "G", "T",""), size = 1)
        } else {
          # Otherwise, keep the original base
          new_base <- substr(dna_seq[[i]], k, k)
        }

        # Append the new base to the new sequence
        new_seq <- paste0(new_seq, new_base)
      }

      # Store the replicated sequence in the output vector
      replicated_dna[[j + (i-1)*num_replicates]] <- list(parent = i, parent_seq = dna_seq[[i]], offspring = new_seq)
    }
  }

  # Convert the replicated sequences to a data frame
  replicated_df <- data.frame(do.call(rbind, replicated_dna))

  # Return the replicated sequences as a data frame
  return(replicated_df)
}

#' Create barcode sequences from parent sequences but with non random substitutions
#'
#' @description
#' #for each the parent sequence coming from the gpseq function, create n offspring sequences with a substitution rate of your interest.
#'
#' @param dna_seq A character vector of DNA sequences.
#' @param num_replicates An integer specifying the number of replicates for each DNA sequence.
#' @param error_rate A number between 0 and 1 representing the error rate for base substitution.
#' @param substitution_probs A list of length 5 containing the substitution probabilities for each base (A, C, G, T, and empty string).
#'
#' @return A data frame containing the replicated barcode sequences.
#'
#' @examples
#' # Example data
#' dna_seq <- c("AAGA", "AATC")
#' num_replicates <- 3
#' error_rate <- 0.1
#' substitution_probs <- list("A" = 0.1, "C" = 0.2, "G" = 0.3, "T" = 0.3, " " = 0.1)
#'
#' # Generate barcode sequences
#' barcode_sequences <- r_gpseq_csub(dna_seq, num_replicates, error_rate, substitution_probs)
#'
#' # Print the resulting barcode sequences
#' print(barcode_sequences)
#'
#' @export
#'

r_gpseq_csub  <- function(dna_seq, num_replicates, error_rate, substitution_probs) {

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
        # Choose the probability of substitution for the current base
        base <- substr(dna_seq[[i]], k, k)
        prob <- substitution_probs[[base]]

        # Generate a random number between 0 and 1
        rand_num <- runif(1)

        # Check if an error occurs
        if (rand_num < error_rate) {
          # If an error occurs, select a random base for substitution
          new_base <- sample(names(substitution_probs), size = 1, prob = unlist(substitution_probs))
        } else {
          # Otherwise, keep the original base
          new_base <- base
        }

        # Append the new base to the new sequence
        new_seq <- paste0(new_seq, new_base)
      }

      # Store the replicated sequence in the output vector
      replicated_dna[[j + (i - 1) * num_replicates]] <- list(parent = i, parent_seq = dna_seq[[i]], offspring = new_seq)
    }
  }

  # Convert the replicated sequences to a data frame
  replicated_df <- data.frame(do.call(rbind, replicated_dna))

  # Return the replicated sequences as a data frame
  return(replicated_df)
}



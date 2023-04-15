#' Create barcode sequences from parent sequences but with non random substitutions
#'
#' @description
#' #for each the parent sequence coming from the gpseq function, create n offspring sequences with a substitution rate of your interest.
#'
#' @param dna_seq A character, output from gpseq
#' @param num_replicates A numeric value
#' @param substitution_probs A list of probabilities for each based pair or an indel to appear
#'
#' @return A data.frame
#'
#' @examples
#' dna_seq <- c("AAGA","AATC")
#' substitution_probs <- list("A" = 0.1, "C" = 0.2, "G" = 0.3, "T" = 0.4, " " = 0.1)
#' r_gpseq_csub(dna_seq,3,substitution_probs)
#'
#' @export
r_gpseq_csub <- function(dna_seq, num_replicates, substitution_probs) {
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
        rand_num <- stats::runif(1)

        # Check if the random number is less than the substitution probability
        if (rand_num < prob) {
          # If so, substitute the base with a random base or an empty string
          new_base <- sample(c("A", "C", "G", "T", ""), size = 1, prob = c(rep(0.25, 4), 0.25))
        } else {
          # Otherwise, keep the original base
          new_base <- base
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

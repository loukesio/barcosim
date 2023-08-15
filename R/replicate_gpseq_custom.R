#' Create barcode sequences from parent sequences with non-random substitutions.
#'
#' @description
#' Given parent DNA sequences, this function generates offspring sequences with potential substitutions.
#' Substitutions are determined by an error rate and a specified set of substitution probabilities.
#'
#' @param dna_seq A character vector of parent DNA sequences.
#' @param num_replicates An integer specifying the number of offspring sequences per parent.
#' @param error_rate A numeric value (0-1) indicating the likelihood of a substitution at each base.
#' @param substitution_probs A named list indicating the substitution probabilities for each base (A, C, G, T, and '').
#' @param color A character flag. If set to "yes", the function will output color-coded sequences.
#'
#' @return A data frame containing the parent sequence, the offspring sequences, and parent index.
#'         If color is "yes", sequences are printed with color-coding and the function returns NULL.
#'
#' @examples
#' dna_seq <- c("AAGA", "AATC")
#' num_replicates <- 3
#' error_rate <- 0.1
#' substitution_probs <- list("A" = 0.1, "C" = 0.2, "G" = 0.3, "T" = 0.3, " " = 0.1)
#' barcode_sequences <- r_gpseq_csub(dna_seq, num_replicates, error_rate, substitution_probs, color = "yes")
#'
#' @export

r_gpseq_csub <- function(dna_seq, num_replicates, error_rate, substitution_probs, color = "no") {

  # Internal function to color-code DNA sequences for visualization.
  colorDNASequences <- function(dna_sequences) {
    # Colors for nucleotide bases.
    base_backgrounds <- c(A = "\033[41m", T = "\033[42m", G = "\033[44m", C = "\033[43m")

    for (dna_sequence in dna_sequences) {
      # Ensure DNA sequence is in uppercase.
      sequence <- toupper(dna_sequence)
      # Split sequence into individual bases.
      bases <- strsplit(sequence, "")[[1]]
      for (base in bases) {
        background_color <- base_backgrounds[base]
        # Print base with appropriate background color.
        cat(paste0("\033[30m", background_color, base, "\033[0m"))
      }
      # Move to the next line for the next sequence.
      cat("\n")
    }
  }

  # Input validation for dna_seq.
  if (!is.character(dna_seq) || any(!nzchar(dna_seq))) {
    stop("Error: dna_seq should be a non-empty string or vector of DNA sequences.")
  }

  # Input validation for num_replicates.
  if (!is.integer(num_replicates) || num_replicates <= 0) {
    if (is.numeric(num_replicates) && num_replicates %% 1 == 0 && num_replicates > 0) {
      num_replicates <- as.integer(num_replicates)
    } else {
      stop("Error: num_replicates should be a positive integer.")
    }
  }

  # Input validation for error_rate.
  if (!is.numeric(error_rate) || error_rate < 0 || error_rate > 1) {
    stop("Error: error_rate should be a number between 0 and 1.")
  }

  # Placeholder for generated sequences.
  replicated_dna <- vector(mode = "list", length = num_replicates * length(dna_seq))

  # Loop over parent sequences.
  for (i in seq_along(dna_seq)) {
    for (j in seq_len(num_replicates)) {
      new_seq <- ""  # Start with an empty sequence.

      # Examine each base for potential substitution.
      for (k in seq_len(nchar(dna_seq[[i]]))) {
        base <- substr(dna_seq[[i]], k, k)
        prob <- substitution_probs[[base]]
        rand_num <- runif(1)

        # Determine if a substitution should occur.
        if (rand_num < error_rate) {
          new_base <- sample(names(substitution_probs), size = 1, prob = unlist(substitution_probs))
        } else {
          new_base <- base
        }

        new_seq <- paste0(new_seq, new_base)
      }

      replicated_dna[[j + (i - 1) * num_replicates]] <- list(parent = i, parent_seq = dna_seq[[i]], offspring = new_seq)
    }
  }

  replicated_df <- data.frame(do.call(rbind, replicated_dna))

  # Output colored sequences if required.
  if (tolower(color) == "yes") {
    colorDNASequences(replicated_df$offspring)
    return(NULL)
  } else {
    return(replicated_df)
  }
}

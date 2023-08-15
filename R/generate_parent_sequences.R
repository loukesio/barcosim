#' Generate DNA sequences by merging initial and random sequences
#'
#' This function generates DNA sequences by merging an initial DNA sequence
#' with random sequences within a specified range. If `color` is set to "yes",
#' the sequences will be colored based on nucleotide type.
#'
#' @param num_sequences Number of sequences to generate (positive integer)
#' @param seq_length Length of each DNA sequence (positive integer)
#' @param range_start Start position of the range within the sequence (integer, 1 to seq_length)
#' @param range_end End position of the range within the sequence (integer, 1 to seq_length, greater than or equal to range_start)
#' @param color A character string, either "yes" or "no", indicating if the sequences should be colored.
#' @return A character vector of generated DNA sequences (when color="no") or printed colored sequences (when color="yes")
#'
#' @examples
#' # Create 4 DNA sequences, 150 bp each, variable in the area between 53-78, and colored.
#' gpseq(4,150,53,78, color="yes")
#'
#' @export

gpseq <- function(num_sequences, seq_length, range_start, range_end, color = "no") {
  # Convert inputs to integers if needed
  num_sequences <- as.integer(num_sequences)
  seq_length <- as.integer(seq_length)
  range_start <- as.integer(range_start)
  range_end <- as.integer(range_end)

  # Check if all inputs are integers
  if (!is.integer(num_sequences) || !is.integer(seq_length) || !is.integer(range_start) || !is.integer(range_end)) {
    stop("All inputs must be integers.")
  }

  # Check if the number of sequences is positive
  if (num_sequences <= 0) {
    stop("Number of sequences must be a positive integer.")
  }

  # Check if the sequence length is positive
  if (seq_length <= 0) {
    stop("Sequence length must be a positive integer.")
  }

  # Check if the range start and end positions are within the sequence length
  if (any(range_start < 1) || any(range_start > seq_length) || any(range_end < 1) || any(range_end > seq_length) || any(range_start > range_end)) {
    stop("Invalid range positions. Range start and end must be valid positions within the sequence length.")
  }

  # Generate initial DNA sequence of seq_length bp
  initial_seq <- paste(sample(c("A", "C", "G", "T"), seq_length, replace = TRUE), collapse = "")

  # Create multiple copies of the initial sequence
  initial_seqs <- rep(initial_seq, num_sequences)

  # Generate random sequences
  random_seqs <- replicate(num_sequences, paste(sample(c("A", "C", "G", "T"), seq_length, replace = TRUE), collapse = ""))

  # Merge initial and random sequences
  final_seqs <- character(num_sequences)
  for (i in 1:num_sequences) {
    seq <- initial_seqs[i]
    for (j in 1:length(range_start)) {
      start <- range_start[j]
      end <- range_end[j]
      seq <- paste(substr(seq, 1, start-1), substr(random_seqs[i], start, end), substr(seq, end+1, seq_length), sep = "")
    }
    final_seqs[i] <- seq
  }

  # If color is set to "yes", color the final sequences
  if (tolower(color) == "yes") {
    colorDNASequences(final_seqs)
  } else {
    return(final_seqs)
  }
}

colorDNASequences <- function(dna_sequences) {
  # Define background colors for bases
  base_backgrounds <- c(A = "\033[41m", T = "\033[42m", G = "\033[44m", C = "\033[43m") # Red, Green, Blue, Yellow background

  # Iterate over each DNA sequence
  for (dna_sequence in dna_sequences) {
    # Convert the sequence to uppercase to handle both upper and lower case letters
    sequence <- toupper(dna_sequence)

    # Iterate over each base in the sequence
    bases <- strsplit(sequence, "")[[1]]
    for (i in seq_along(bases)) {
      # Get the current base and background color
      base <- bases[i]
      background_color <- base_backgrounds[base]

      # Print the base surrounded by color with black foreground and background color
      cat(paste0("\033[30m", background_color, base, "\033[0m"))
    }

    cat("\n")  # Print a new line after each sequence
  }
}



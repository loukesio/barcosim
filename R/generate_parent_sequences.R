#' Generate DNA sequences by merging initial and random sequences
#'
#' This function generates DNA sequences by merging an initial DNA sequence
#' with random sequences within a specified range.
#'
#' @param num_sequences Number of sequences to generate (positive integer)
#' @param seq_length Length of each DNA sequence (positive integer)
#' @param range_start Start position of the range within the sequence (integer, 1 to seq_length)
#' @param range_end End position of the range within the sequence (integer, 1 to seq_length, greater than or equal to range_start)
#' @return A character vector of generated DNA sequences
#'
#' @examples
#' # Create 4 DNA sequneces, 150 bp each, variable in the are between 53-78.
#' gpseq(4,150,53,78)
#'
#' @export

gpseq <- function(num_sequences, seq_length, range_start, range_end) {

  # Convert inputs to integers
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
  if (range_start < 1 || range_start > seq_length || range_end < 1 || range_end > seq_length || range_start > range_end) {
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
    final_seqs[i] <- paste(substr(initial_seqs[i], 1, range_start-1), substr(random_seqs[i], range_start, range_end), substr(initial_seqs[i], range_end+1, seq_length), sep = "")
  }

  return(final_seqs)
}


gpseq(10, 10, 3, 6)


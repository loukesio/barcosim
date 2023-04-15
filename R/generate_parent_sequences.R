#' Generate parent sequences
#'
#'`gpseq` creates DNA sequences of certain length. These sequences are indetical apart from a specific range that refers to barcodes and its variable.
#'
#' @param num_sequences A numeric value
#' @param seq_length A numeric value
#' @param range_start A numeric value
#' @param range_end A numeric value
#'
#' @return A character vector, DNA sequences.
#' @examples
#' # Create 4 DNA sequneces, 150 bp each, variable in the are between 53-78.
#' gpseq(4,150,53,78)
#'
#' @export

gpseq <- function(num_sequences, seq_length, range_start, range_end) {

  # Generate initial DNA sequence of 150 bp
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

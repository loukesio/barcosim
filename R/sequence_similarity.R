#' Calculate sequence similarity
#'
#' This function calculates the similarity between sequences at each position.
#'
#' @param sequences A character vector of sequences to compare.
#' @return A numeric vector of similarity values for each position.
#'
#' @examples
#' # A dataset of DNA sequences
#'
#' sequences <- c("CGCCTCCACG", "CGAGCACACG", "CGGAGTCACG", "CGCCAACACG", "CGTAAGCACG",
#'                "CGCCAACACG", "CGTGCGCACG", "CGCTGTCACG", "CGGCCACACG", "CGGCGTCACG")
#'
#' # Checking at the similarity among them
#'
#' similarity <- calSeqSim(sequences)
#' similarity
#'
#' @export

calSeqSim <- function(sequences) {
  num_positions <- nchar(sequences[1])
  similarity <- numeric(num_positions)

  for (i in 1:num_positions) {
    position <- sapply(sequences, function(seq) substr(seq, i, i))
    unique_nucleotides <- unique(position)
    num_unique_nucleotides <- length(unique_nucleotides)

    if (num_unique_nucleotides == 1) {
      similarity[i] <- 100  # 100% similarity when all sequences have the same nucleotide at the position
    } else if (num_unique_nucleotides == length(sequences)) {
      similarity[i] <- 0  # 0% similarity when all sequences have different nucleotides at the position
    } else {
      count <- table(position)
      max_count <- max(count)
      similarity[i] <- max_count / sum(count) * 100
    }
  }

  return(similarity)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/replicate_gpseq.R
\name{r_gpseq}
\alias{r_gpseq}
\title{Create barcode sequences from parent sequences}
\usage{
r_gpseq(dna_seq, num_replicates, error_rate, color = "no")
}
\arguments{
\item{dna_seq}{A character, output from gpseq}

\item{num_replicates}{A numeric value}

\item{error_rate}{A numeric value from 0 to 1}

\item{color}{A character string, either "yes" or "no", indicating if the sequences should be colored.}
}
\value{
a data.frame
}
\description{
The r_gpseq function allows you to replicate the parent DNA sequences, that have been
created from the gpseq function.
}
\examples{
dna_seq <- c("AAGA","AATC")
r_gpseq(dna_seq,3,0.1, color="yes") #create 3 offspring sequences with random substitution

}

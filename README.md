[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/icon)](https://cran.r-project.org/package=icons)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Install the BarcoSim package
Install the package using the following commands  <img align="right" src="logo/BarcoSim_Logo.png" width=400>

```r
# for now you can install the developemental version of ltc
# first you need to install the devtools package 
# in case you have not already installed
install.packages("devtools") 
# and load it
library(devtools)

# then you can install the dev version of the ltc
devtools::install_github("loukesio/BarcoSim")
# and load it
library(BarcoSim)
```

### 1. Use the `gpseq` command to generate the parent sequences.

``` r
library(Biostrings) # Provides tools for working with biological sequences, such as DNA, RNA, and protein sequences
library(BarcoSim)   # BarcoSim: A package for simulating barcoded sequencing data
library(dplyr)      # A powerful package for data manipulation and transformation,


set.seed(123)       # sets the random seed to ensure the reproducibility of a random processes (generation of sequences)

#num_sequences = Number of sequences to generate 
#seq_length = Length of each DNA sequence
#range_start = Start position of the barcoded sequence
#range_end = End position of the barcoded sequence

# This function creates 10 parent sequences, each with 10 base pairs. The barcode ranges from base 3 to base 6.
df1 <- gpseq(num_sequences=10, seq_length=10, range_start=3, range_end=6)

df1 
  DNAStringSet()
#> DNAStringSet object of length 10:
#>      width seq
#>  [1]    10 TCGGGCCTGC
#>  [2]    10 TCCCCTCTGC
#>  [3]    10 TCCACCCTGC
#>  [4]    10 TCGTAGCTGC
#>  [5]    10 TCCCGTCTGC
#>  [6]    10 TCTTAGCTGC
#>  [7]    10 TCGAGCCTGC
#>  [8]    10 TCGTTGCTGC
#>  [9]    10 TCTATTCTGC
#> [10]    10 TCATTTCTGC
```

<sup>Created on 2023-04-15 with [reprex v2.0.2](https://reprex.tidyverse.org)</sup>

### 2. Use the `r_gpseq` command to replicate parent sequences and make a barcode data set.

``` r
library(Biostrings)
library(BarcoSim)
library(dplyr)

df1 <- gpseq(10, 10, 3, 6)
r_gpseq(df1, 4, 0.01) 
#>    parent parent_seq  offspring
#> 1       1 GTGGTGGTTT GTGGTGGTTT
#> 2       1 GTGGTGGTTT GTGGTGGTTT
#> 3       1 GTGGTGGTTT GTGGTGGTTT
#> 4       1 GTGGTGGTTT GTGGTGGTTT
#> 5       2 GTAGCCGTTT GTAGCCGTTT
#> 6       2 GTAGCCGTTT GTAGCCGTTT
#> 7       2 GTAGCCGTTT GTAGCCGTTT
#> 8       2 GTAGCCGTTT GTAGCCGTTT
#> 9       3 GTGGACGTTT GTGGACGTTT
#> 10      3 GTGGACGTTT GTGGACGTTT
#> 11      3 GTGGACGTTT GTGGACGTTT
#> 12      3 GTGGACGTTT GTGGACGTTT
#> 13      4 GTCACCGTTT GTCACCGTTT
#> 14      4 GTCACCGTTT GTCACCGTTT
#> 15      4 GTCACCGTTT GTCACCGTTT
#> 16      4 GTCACCGTTT GTCACCGTTT
#> 17      5 GTGGAGGTTT GTGGAGGTTT
#> 18      5 GTGGAGGTTT GTGGAGGTTT
#> 19      5 GTGGAGGTTT GTGGAGGTTT
#> 20      5 GTGGAGGTTT GTGGAGGTTT
#> 21      6 GTTTCAGTTT GTTTCAGTTT
#> 22      6 GTTTCAGTTT GTTTCAGTTT
#> 23      6 GTTTCAGTTT GTTTCAGTTT
#> 24      6 GTTTCAGTTT GTTTCAGTTT
#> 25      7 GTGGCTGTTT GTGGCTGTTT
#> 26      7 GTGGCTGTTT GTGGCTGTTT
#> 27      7 GTGGCTGTTT GTGGCTGTTT
#> 28      7 GTGGCTGTTT GTGGCTGTTT
#> 29      8 GTACATGTTT GTACATGTTT
#> 30      8 GTACATGTTT GTACATGTTT
#> 31      8 GTACATGTTT GTACATGTTT
#> 32      8 GTACATGTTT GTACATGTTT
#> 33      9 GTGCCGGTTT GTGCCGGTTT
#> 34      9 GTGCCGGTTT GTGCCGGTTT
#> 35      9 GTGCCGGTTT GTGCCGGTTT
#> 36      9 GTGCCGGTTT GTGCCGGTTT
#> 37     10 GTGAATGTTT GTGAATGTTT
#> 38     10 GTGAATGTTT GTGAATGTTT
#> 39     10 GTGAATGTTT GTGAATGTTT
#> 40     10 GTGAATGTTT GTGAATGTTT
```

<sup>Created on 2023-04-15 with [reprex v2.0.2](https://reprex.tidyverse.org)</sup>

### 3. Use the `r_gpseq_csub` command to replicate parent sequences with a certain error rate and a certain subsitution rate.

```
dna_seq <- c("AAGA","AATC")
substitution_probs <- list("A" = 0.1, "C" = 0.2, "G" = 0.3, "T" = 0.4, " " = 0.1)
r_gpseq_csub(dna_seq,3,substitution_probs)
```

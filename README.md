## Dashing Experiments

This repository contains code and data for reproducing experiments from the manuscript accompanying
the Dashing software application, available at https://github.com/dnbaker/dashing.

#### dsexp
dsexp contains experiments testing the performance of various data structures for Jaccard-coefficient calculation.

1. dsexp.cpp
  1. Performs numerical simulations comparing the error rates of bloom filters, minhash sketches, and hyperloglogs at specific sketch sizes for Jaccard-coefficient estimation of sets of varying sizes.
2. dsexp.Rmd
  1. Contains code for visualizing results from dsexp.cpp, which can be used to reproduce Fig. 1 and Supplementary Table 1 from the manuscript.

#### timing

1. `all_pairwise.py`
  1. Performs all pairwise comparisons between a set of genomes across varying sketch size and kmer length, comparing bindash, mash, and several estimation methods for HyperLogLogs.

####  accuracy

1. pairselector.py
  1. This script finds candidate genome pairs in specified ranges of Jaccard indices from a large, upper-triangular table of pairwise distances.
    1. We generated this table with `dashing dist`, with k=31 and p=16.
2. `pairwise_benchmark.cpp`
  1. For all pairs of genomes provided, calculate the exact Jaccard index with hash sets, and report the difference and errors between
     estimates and the true value for bindash, mash, and several estimation methods for HyperLogLogs.
  2. This can be rather memory-intensive due to the use of full hash sets; for this reason, we suggest omitting large genomes from the call generating the table used in pairselector.py.
3. `genomes_for_exp.txt`
  1. This is the set of genomes emitted by pairselector.py

# epiG
### *Statistical Inference of Epi-allelic Patterns from Whole-Genome Bisulphite Sequencing Data*

[![Travis-CI Build Status](https://travis-ci.org/vincent-dk/epiG.svg?branch=master)](https://travis-ci.org/vincent-dk/epiG)

[Epig R package manual](epiG-manual.pdf)

## Installation

Get the development version from github:

```R
# install.packages("devtools")
devtools::install_github("vincent-dk/epiG")
```

## Quick Start

### 1. Load the epiG R package
```R
# Load epiG package
library(epiG)
```

### 2. Create a epiG configuration
```R
# Specify site (GNAS)
chr <- "chr20"
start <- 57380000
end <- 57478000

# Create epiG configuration
config <- auto_config(
  bam_file = "SRX332736_sorted.bam",
  ref_file = "hg19_rCRSchrm.fa",
  alt_file = "dbsnp_135.hg19.fa",
  seq_type = "BSeq", # Std Bisulphite sequencing
  chr = chr,
  start = start,
  end = end)
```

### 3. Infer epigenomic haplotype patterns by running the epiG algorithm  

```R
# Run epiG
g <- epiG(max_threads = 3, config = config)

# Load additional data into model
g <- fetch_reads(g) # Load reads
g <- fetch_ref(g) # Load reference genome
g <- fetch_alt(g) # Load alternative nucleotides

# Print a summary of the inferred epigenomic haplotypes
g
```

### 4a. Get a overview of the inferred chains

```R
cinfo <- chai_info(g)
subset(cinfo, nreads > 10)
```

### 4b. Look at a specific position
```R
info <- position_info(g, 57415143 + c(0,1))
info[, c("position", "chain.id", "ref", "genotype", "fit.ratio",
 "methylated", "nreads.fwd", "nreads.rev")]
```

### 5. Make a plot
> More info coming soon...

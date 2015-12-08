############################################################
###
### LIBS ect
###

library(epiG)

raw.data.dir <- "/home/martin/Dropbox/epig/raw_data/"

############################################################
###
### SETUP
###

### bam file
# bam.file <- '/home/martin/seqdata/SRX729714/SRX729714_sorted.bam'
# name <- "SRX729714"

# bam.file <- '/home/martin/seqdata/SRX332737/SRX332737_sorted.bam'
# name <- "SRX332737"

bam.file <- '/home/martin/seqdata/SRX332736/SRX332736_sorted.bam'
name <- "SRX332736"

#! site

# GNAS
chr <- "chr20"
start <- 57380000
end <- 57478000
site.name <- "GNAS"

############################################################
###
### RUN epiG
###

### epiG configuration
config <- auto_config(
  bam_file = bam.file,
  ref_file = paste(raw.data.dir, "ref_data/hg19_rCRSchrm.fa", sep = ""),
  alt_file = paste(raw.data.dir, "ref_data/dbsnp_135.hg19.fa", sep = ""),
  chr = chr,
  start = start,
  end = end,
  use_paired_reads = TRUE,
  chunk_size = 150000,
  delta = 0.9)

fit <- epiG(bam.file, chr, start, end, max_threads = 3, config = config)

fit <- fetch.reads(fit)
fit <- fetch_ref(fit)
fit <- fetch_alt(fit)

### SAVE

save(fit, file = paste("/home/martin/Dropbox/epig/version12/code for plots/GNAS high cov plot/fit_", name,"_", site.name,".RData", sep=""))

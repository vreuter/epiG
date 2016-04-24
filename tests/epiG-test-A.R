library(epiG)

bam_file <- system.file("extdata", "GNAS_small.bam", package="epiG")
ref_file <- system.file("extdata", "hg19_rCRSchrm_chr20.fa", package="epiG")
alt_file <- system.file("extdata", "dbsnp_135.hg19_chr20.fa", package="epiG")

# Test file_info
file_info(bam_file)

# GNAS
chr <- "chr20"
start <- 57400000 
end <- 57400000 + 1000

### Build epiG configuration
config <- auto_config(
		bam_file = bam_file,
		ref_file = ref_file,
		alt_file = alt_file,
		chr = chr,
		start = start,
		end = end)

#### Run epiG
fit <- epiG(max_threads = 2, config = config)

#### Fetch
fit <- fetch_reads(fit)

fit <- fetch_ref(fit)

fit <- fetch_alt(fit)

# Test print
fit
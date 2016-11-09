library(epiG)

bam_file <- system.file("extdata", "GNAS_small.bam", package="epiG")
ref_file <- system.file("extdata", "hg19_GNAS.fa", package="epiG")
alt_file <- system.file("extdata", "dbsnp_135.hg19_GNAS.fa", package="epiG")

# Test file_info
file_info(bam_file)

# GNAS
chr <- "chr20"
start <- 57400000
end <- 57400000 + 1000

### Build epiG configuration
config_1 <- auto_config(
		bam_file = bam_file,
		ref_file = ref_file,
		ref_offset = 57380000,
		alt_file = alt_file,
		alt_offset = 57380000,
		chr = chr,
		start = start,
		end = end
  )

config_2 <- auto_config(
    bam_file = bam_file,
    ref_file = ref_file,
    ref_offset = 57380000,
    alt_file = alt_file,
    alt_offset = 57380000,
    chr = chr,
    start = end,
    end = end + 1000
  )

#### Run epiG.chunks
configs <- list(config_1, config_2)
fit <- epiG_chunks(configs, max_threads = 2)

#### Fetch
fit <- fetch_reads(fit)

fit <- fetch_ref(fit)

fit <- fetch_alt(fit)

# Test print
fit

# Test position info
dummy <- position_info(fit, start(fit[[1]]):end(fit[[1]]))
#TODO content tests

# Test read_info
dummy <- read_info(fit)
dummy <- read_info(fit, inc.symbols = TRUE)
#TODO content tests

# Test subregion
sub <- subregion(fit, start + 500, end + 500)
#TODO content tests

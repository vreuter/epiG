#
#     Description of this R script:
#     epiG-configuration scripts. 
#     This scripts contain functions need for creating epiG configurations
#     and bisulphite models
#
#     Intended for use with R.
#     Copyright (C) 2016 Martin Vincent
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http:#www.gnu.org/licenses/>
#


# Create bisulphite conversion model
.create_BS_model <- function(bisulphite_rate, bisulphite_inap_rate) {
	
	#TODO split up into 2 functions one for fwd model and one for rev model
	
	bisulphite_model <- list()
	
	bisulphite_model$fwd <- matrix(nrow = 4, ncol = 6)
	bisulphite_model$rev <- matrix(nrow = 4, ncol = 6)
	
	rownames(bisulphite_model$fwd) <- c('C', 'G', 'A', 'T')	
	rownames(bisulphite_model$rev) <- c('C', 'G', 'A', 'T')	
	colnames(bisulphite_model$rev) <- c('C', 'G', 'A', 'T', 'C^me', 'G_me')
	colnames(bisulphite_model$fwd) <- c('C', 'G', 'A', 'T', 'C^me', 'G_me')
	
	bisulphite_model$fwd[,] <- 0
	bisulphite_model$rev[,] <- 0
	
	bisulphite_model$fwd[1, 1] <- 1 - bisulphite_rate #C C
	bisulphite_model$fwd[1, 5] <- 1 - bisulphite_inap_rate #C c
	bisulphite_model$fwd[4, 5] <- bisulphite_inap_rate #T c
	bisulphite_model$fwd[2, 2] <- 1 #G G
	bisulphite_model$fwd[2, 6] <- 1 #G g
	bisulphite_model$fwd[3, 3] <- 1 #A A
	bisulphite_model$fwd[4, 1] <- bisulphite_rate #T C
	bisulphite_model$fwd[4, 4] <- 1 #T T
	
	bisulphite_model$rev[1, 1] <- 1 #C C
	bisulphite_model$rev[1, 5] <- 1 #C c
	bisulphite_model$rev[2, 2] <- 1 - bisulphite_rate #G G
	bisulphite_model$rev[2, 6] <- 1 - bisulphite_inap_rate #G g
	bisulphite_model$rev[3, 6] <- bisulphite_inap_rate #A g
	bisulphite_model$rev[3, 3] <- 1 #A A
	bisulphite_model$rev[3, 2] <- bisulphite_rate #A G
	bisulphite_model$rev[4, 4] <- 1 #T T
	
	return(bisulphite_model)
}

#' @title Bisulphite Conversion Model
#' @description
#' Create a bisulphite sequencing conversion model
#' 
#' @param bisulphite_rate bisulphite conversion rate (numeric in the range (0, 1])
#' @param bisulphite_inap_rate bisulphite inappropriate conversion rate (numeric in the range (0, 1])
#' @param Lmax maximal read length (integer)
#' @param ... ignored
#' @return an epiG conversion model
#' 
#' @author Martin Vincent
#' @export
#' @examples 
#' BSeq()
BSeq <- function(bisulphite_rate = .95, bisulphite_inap_rate = 0.05, Lmax = 110, ...) {
	
	bisulphite_rates <- rep(bisulphite_rate, Lmax)
	bisulphite_inap_rates <- rep(bisulphite_inap_rate, Lmax)
	
	model <- list()
	model$fwd <-lapply(1:Lmax, function(i) .create_BS_model(bisulphite_rates[i], bisulphite_inap_rates[i])$fwd) 
	model$rev <-lapply(1:Lmax, function(i) .create_BS_model(bisulphite_rates[i], bisulphite_inap_rates[i])$rev) 
	
	model$name <- "Bisulphite Sequencing Conversion Model"
	
	model$use_split <- FALSE
	
	class(model) <- "epiG.model"
	
	return(model)
}


#' @title NOMe-sequencing Conversion Model
#' @description
#' Create a NOMe-seq conversion model
#' 
#' @param bisulphite_rate bisulphite conversion rate (numeric in the range (0, 1])
#' @param bisulphite_inap_rate bisulphite inappropriate conversion rate (numeric in the range (0, 1])
#' @param Lmax maximal read length (integer)
#' @param ... ignored
#' @return an epiG conversion model
#' 
#' @author Martin Vincent
#' @export
#' @examples 
#' NOMeSeq()
NOMeSeq <- function(bisulphite_rate = .95, bisulphite_inap_rate = 0.05, Lmax = 110, ...) {
	
	model <- BSeq(
		bisulphite_rate = bisulphite_rate, 
		bisulphite_inap_rate = bisulphite_inap_rate, 
		Lmax = Lmax)

	
	ignor.me <- diag(6)
	rownames(ignor.me) <- c('C', 'G', 'A', 'T', 'C^me', 'G_me')
	colnames(ignor.me) <- c('C', 'G', 'A', 'T', 'C^me', 'G_me')
	ignor.me[c(1, 5), c(1, 5)] <- 1/2
	ignor.me[c(2, 6), c(2, 6)] <- 1/2
			
	model$fwd_HCGD <- model$fwd
	model$rev_HCGD <- model$rev
		
	model$fwd_DGCH <- model$fwd
	model$rev_DGCH <- model$rev
		
	model$fwd_CH <- lapply(model$fwd, function(x) x %*% ignor.me)
	model$rev_CH <- lapply(model$rev, function(x) x %*% ignor.me)
		
	model$fwd_C_G <- lapply(model$fwd, function(x) x %*% ignor.me)
	model$rev_C_G <- lapply(model$rev, function(x) x %*% ignor.me)
	
	
	model$name <- "NOMe-sequencing Conversion Model"
	
	model$use_split <- TRUE
	
	class(model) <- "epiG.model"
	
	return(model)
}


#' @title Create Standard Configuration
#' @description
#' Create a epiG configuration with standard parameters
#' 
#' @param ref_file genome reference file (path to .fa file)
#' @param alt_file alternative nucleotide file (path to .fa file)
#' @param bam_file bam file (path to .bam file)
#' @param chr reference name
#' @param start start position of region to processes 
#' @param end end position of region to processes 
#' @param seq_type sequencing type ("BSeq" for bisulphite sequencing, "NOMeSeq" for NOMe sequencing)
#' @param paired_reads should pair information be used (TRUE/FALSE or NULL, if NULL then paired information is used if pairs are present in bam file)
#' @param min_overlap minimum overlapping length (if NULL deafult value is used)
#' @param min_CG minimum overlapping CG positions (if NULL deafult value is used)
#' @param min_HCGD minimum overlapping HCGD positions (if NULL deafult value is used)
#' @param min_DGCH minimum overlapping DGCH positions (if NULL deafult value is used)
#' @param ... additional arguments (overrides default values)
#' 
#' @return An epiG configuration
#' 
#' @author Martin Vincent
#' @export
#' @examples 
#' 
#' # Retrieve paths to raw data files
#' bam_file <- system.file("extdata", "GNAS_small.bam", package="epiG")
#' ref_file <- system.file("extdata", "hg19_GNAS.fa", package="epiG")
#' alt_file <- system.file("extdata", "dbsnp_135.hg19_GNAS.fa", package="epiG")
#' 
#' # Specify region
#' chr <- "chr20"
#' start <- 57400000 
#' end <- 57400000 + 1000
#' 
#' # Build epiG configuration
#' config <- auto_config(
#' 		bam_file = bam_file,
#' 		ref_file = ref_file,
#' 		alt_file = alt_file,
#' 		chr = chr,
#' 		start = start,
#' 		end = end,
#' # If ref_file and alt_file contains the entire chromosome this is not needed
#' 		ref_offset = 57380000, 
#' 		alt_offset = 57380000)
#'
#' config #print a summary of the configuration
auto_config <- function(
		ref_file, 
		alt_file, 
		bam_file,
		chr, 
		start, 
		end, 
		seq_type = "BSeq",
		paired_reads = NULL,
		min_CG = NULL,
		min_HCGD = NULL,
		min_DGCH = NULL,
		min_overlap = NULL,
		...) {

	if(length(chr) == 0 || length(chr) != length(start) || length(chr) != length(end) || length(start) != length(end)) {
		stop("length of chr, start and end must be equal and nonzero")
	}
		
	if(length(chr) > 1) {
		
		confs <- lapply(1:length(chr), function(i) 
			auto_config(
				ref_file = ref_file, 
				alt_file = alt_file, 
				bam_file = bam_file, 
				chr = chr[i], 
				start = start[i], 
				end = end[i],
				seq_type = seq_type,
				paired_reads = paired_reads,
				min_CG = min_CG,
				min_HCGD = min_HCGD,
				min_DGCH = min_DGCH,
				min_overlap = min_overlap,
				...
			))
		
		return(confs)
	}
	
	# Fetch information about file
	
	reads <- fetch_read_info(bam_file, chr, start, end)
	n_reads <- nrow(reads)
	
	if(n_reads == 0) {
		warning("No reads found in region")
		return (NULL)
	}
	
	# Paired design
	if(is.null(paired_reads))  {
		if(max(table(reads$name[1:1000])) == 2) {
			paired_reads <- TRUE
		} else {
			paired_reads <- FALSE
		}
	}

	# Create conversion model

	if(seq_type == "NOMeSeq") {
		
		model <- NOMeSeq(Lmax = max(reads$length), ...)
		
		default_min_CG <- 0
		default_min_HCGD <- 0
		default_min_DGCH <- 2
		default_min_overlap <- 40
		
		
	} else if(seq_type == "BSeq") {
		
		model <- BSeq(Lmax = max(reads$length), ...)

		
		if(paired_reads) {
			default_min_CG <- 2
			default_min_overlap <- 50
			
		} else {
			default_min_CG <- 1
			default_min_overlap <- 40
		}
		
		default_min_HCGD <- 0
		default_min_DGCH <- 0
		
	} else {
		stop("Unknown seq_type")
	}
		
	
	

	# Create configuration
	
	config <- epiG_config(
			model = model,
			hard_limit = n_reads + 1000,
			chunk_size = n_reads,
			ref_file = ref_file,
			alt_file = alt_file,
			min_overlap = if(is.null(min_overlap)) default_min_overlap else min_overlap,
			min_CG = if(is.null(min_CG)) default_min_CG else min_CG,
			min_HCGD = if(is.null(min_HCGD)) default_min_HCGD else min_HCGD,
			min_DGCH = if(is.null(min_DGCH)) default_min_DGCH else min_DGCH,
			paired_reads = paired_reads,
			...
			)
	
	# Set run configuration
	config <- set_run_configuration(config, bam_file, chr, start, end)

	
	return(config)
}	

#TODO match C-side config names with R-side
#TODO examples epiG_config
#' @title Create an epiG Configuration
#' @description
#' Create a custom epiG configuration
#' 
#' @param model conversion model
#' @param ref_file genome reference file (path to .fa file)
#' @param alt_file alternative nucleotide file (path to .fa file)
#' @param min_overlap minimum overlapping length
#' @param min_CG minimum overlapping CG positions
#' @param min_HCGD minimum overlapping HCGD positions
#' @param min_DGCH minimum overlapping DGCH positions
#' @param paired_reads used pair information (reads with the same name in the bam file is paired and will be forced into the same haplotype chain)
#' @param ref_prior genotype prior parameter
#' @param structual_prior structural prior scaling
#' @param margin cut off margin
#' @param max_iterations maximal number of iterations
#' @param max_stages experimental stage optimization (if <= 1 then stage optimization is off)
#' @param chunk_size chunk size
#' @param chunk_method chunk method ('none' only one chunk, 'reads' chunks of approximately chunk_size reads, 'bases' chunks of chunk_size bases)
#' @param hard_limit maximal number of reads loaded per chunk (reads not loaded will be completely ignored)
#' @param ref_offset ref file offset (usually ref_offset = 1)
#' @param alt_offset alt file offset (usually alt_offset = 1)
#' @param quality_threshold discard reads with mean epsilon quality higher than quality_threshold
#' @param verbose show information while running
#' @param ... ignored
#' 
#' @return an epiG configuration
#' 
#' @author Martin Vincent
#' @export
epiG_config <- function(
		model, 
		ref_file, 
		alt_file, 
		min_overlap,
		min_CG,
		min_HCGD,
		min_DGCH,
		paired_reads,
		ref_prior = 1-1e-4, 
		structual_prior = 1,
		quality_threshold = 0.020,
		margin = 5,
		max_iterations = 1e5, 
		max_stages = 1,
		chunk_size = 5000, 
		chunk_method = "reads", 
		hard_limit = chunk_size + 1000,
		ref_offset = 1,
		alt_offset = 1,
		verbose = TRUE,
		...) {
	
	#TODO check config valid
	# 1) chunk_size < reads_hard_limit
	
	config <- list()
	
	config$ref_filename <- ref_file
	config$alt_filename <- alt_file
	config$ref_offset <- as.integer(ref_offset)-1L
	config$alt_offset <- as.integer(alt_offset)-1L 
	
	config$max_iterations <- as.integer(max_iterations)
	
	config$fwd_model <- model$fwd
	config$rev_model <- model$rev
	
	config$split_mode <- model$use_split

	config$model_name <- model$name
	
	if(config$split_mode) {
		
		config$fwd_DGCH_model <- model$fwd_DGCH
		config$rev_DGCH_model <- model$rev_DGCH
	
		config$fwd_HCGD_model <- model$fwd_HCGD 
		config$rev_HCGD_model <- model$rev_HCGD
		
		config$fwd_CH_model <- model$fwd_CH 
		config$rev_CH_model <- model$rev_CH
	
		config$fwd_C_G_model <- model$fwd_C_G 
		config$rev_C_G_model <- model$rev_C_G
		
	} else {
		
		config$fwd_DGCH_model <- list(matrix(0))
		config$rev_DGCH_model <- list(matrix(0))
		
		config$fwd_HCGD_model <- list(matrix(0))
		config$rev_HCGD_model <- list(matrix(0))
		
		config$fwd_CH_model <- list(matrix(0))
		config$rev_CH_model <- list(matrix(0))
		
		config$fwd_C_G_model <- list(matrix(0))
		config$rev_C_G_model <- list(matrix(0))
		
	}

	config$chunk_size <- chunk_size
	
	config$chunk_method <- chunk_method
	
	config$reads_hard_limit <- as.integer(hard_limit)
	
	config$quality_threshold <- quality_threshold
	
	config$ref_prior <- ref_prior
	
	config$min_overlap_length <- as.integer(min_overlap)
	
	config$min_CG_count <- as.integer(min_CG)
	config$min_HCGD_count <- as.integer(min_HCGD)
	config$min_DGCH_count <- as.integer(min_DGCH)
	config$margin <- as.integer(margin)
	
	config$max_stages <- as.integer(max_stages)
	
	config$structual_prior_scale <- structual_prior
		
	config$use_paired_reads <- paired_reads
				
	config$verbose <- verbose
	
	class(config) <- "epiG.config"
	
	return(config)
}

#TODO example set_run_configuration
#' @title Run Configuration
#' @description
#' set run configuration
#' 
#' @param config an epiG configuration 
#' @param filename bam file (path to .bam file)
#' @param refname reference name
#' @param start start position of region to processes 
#' @param end end position of region to processes 
#'
#' @return an epiG configuration
#'
#' @author Martin Vincent
#' @export
set_run_configuration <- function(config, filename, refname, start, end) {
	
	config$filename <- filename
	config$start <- start
	config$end <- end
	config$refname <- refname
	
	return(config)
}

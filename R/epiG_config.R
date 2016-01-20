# TODO: Add comment
# 
# Author: martin
###############################################################################
#' create_error_distributions
#' 
#' @param bisulfite_rate 
#' @param bisulfite_inap_rate 
#' @return bisulfite model
#' 
#' @author Martin Vincent
#' @export
create_error_distributions <- function(bisulfite_rate, bisulfite_inap_rate) {
	
	#TODO split up into 2 functions one for fwd model and one for rev model
	
	bisulfite_model <- list()
	
	bisulfite_model$fwd <- matrix(nrow = 4, ncol = 6)
	bisulfite_model$rev <- matrix(nrow = 4, ncol = 6)
	
	rownames(bisulfite_model$fwd) <- c('C', 'G', 'A', 'T')	
	rownames(bisulfite_model$rev) <- c('C', 'G', 'A', 'T')	
	colnames(bisulfite_model$rev) <- c('C', 'G', 'A', 'T', 'C^me', 'G_me')
	colnames(bisulfite_model$fwd) <- c('C', 'G', 'A', 'T', 'C^me', 'G_me')
	
	bisulfite_model$fwd[,] <- 0
	bisulfite_model$rev[,] <- 0
	
	bisulfite_model$fwd[1, 1] <- 1 - bisulfite_rate #C C
	bisulfite_model$fwd[1, 5] <- 1 - bisulfite_inap_rate #C c
	bisulfite_model$fwd[4, 5] <- bisulfite_inap_rate #T c
	bisulfite_model$fwd[2, 2] <- 1 #G G
	bisulfite_model$fwd[2, 6] <- 1 #G g
	bisulfite_model$fwd[3, 3] <- 1 #A A
	bisulfite_model$fwd[4, 1] <- bisulfite_rate #T C
	bisulfite_model$fwd[4, 4] <- 1 #T T
	
	bisulfite_model$rev[1, 1] <- 1 #C C
	bisulfite_model$rev[1, 5] <- 1 #C c
	bisulfite_model$rev[2, 2] <- 1 - bisulfite_rate #G G
	bisulfite_model$rev[2, 6] <- 1 - bisulfite_inap_rate #G g
	bisulfite_model$rev[3, 6] <- bisulfite_inap_rate #A g
	bisulfite_model$rev[3, 3] <- 1 #A A
	bisulfite_model$rev[3, 2] <- bisulfite_rate #A G
	bisulfite_model$rev[4, 4] <- 1 #T T
	
	return(bisulfite_model)
}

create_pcr_model <- function(rate) 1/3*rate + (1-4/3*rate)*diag(4)

#' create_bisulfite_model
#' 
#' @param bisulfite_rates 
#' @param bisulfite_inap_rate 
#' @param lambda 
#' @param Lmax 
#' @return ...
#' 
#' @author Martin Vincent
#' @export
create_bisulfite_model <- function(bisulfite_rates, bisulfite_inap_rate, Lmax) {
	
	#TODO postion decay
	#p <- exp_decay(lambda = lambda, Lmax = Lmax)
	#bisulfite_rates <- 1 - (1-bisulfite_rate)*Lmax*p
	
	bisulfite_rates <- rep(bisulfite_rates, Lmax)
	bisulfite_inap_rate <- rep(bisulfite_inap_rate, Lmax)
	
	model <- list()
	model$fwd <-lapply(1:Lmax, function(i) create_error_distributions(bisulfite_rates[i], bisulfite_inap_rate[i])$fwd) 
	model$rev <-lapply(1:Lmax, function(i) create_error_distributions(bisulfite_rates[i], bisulfite_inap_rate[i])$rev) 
	
	return(model)
}

create_haplo_prior <- function(delta, n_reads){
	# Compute log haplochain prior values
	return((0:n_reads)^2/n_reads)
} 

#TODO auto config for chunked mode
auto_config <- function(
		bam_file,
		ref_file, 
		alt_file, 
		chr, 
		start, 
		end, 
		min_overlap = NULL,
		...) {

	if(length(chr) != length(start) || length(chr) != length(end) || length(start) != length(end)) {
		stop("length of chr, start and end must be equal")
	}
	
	if(length(chr) > 1) {
		
		confs <- lapply(1:length(chr), function(i) auto_config(
							bam_file, 
							ref_file, 
							alt_file, 
							chr[i], 
							start[i], 
							end[i], 
							use_paired_reads,
							NOMEseq_mode,
							chunk_size,
							min_overlap))
		
		
		return(confs[ ! sapply(confs, is.null)])
	}
	
	reads <- fetch_reads_info(bam_file, chr, start, end)
	
	if(nrow(reads) == 0) {
		warning("No reads found") #TODO more info
		return (NULL)
	}
	
	n_reads <- nrow(reads)
	
	### Create bisulfite model
	#TODO auto detimen bisulfite rates
	model <- create_bisulfite_model(bisulfite_rate = .95, bisulfite_inap_rate = 0.05, Lmax = max(reads$length))
	
#	pcr.model <- create_pcr_model(rate = 0.01)
#	
#	model$fwd <- lapply(model$fwd, function(x) pcr.model %*% x)
#	model$rev <- lapply(model$rev, function(x) pcr.model %*% x)
	
	#TODO auto nomeseq mode
	model$fwd_HCGD <- model$fwd
	model$rev_HCGD <- model$rev
	
	model$fwd_DGCH <- model$fwd
	model$rev_DGCH <- model$rev 
	
	model$fwd_C_G <- model$fwd
	model$rev_C_G <- model$rev
	
	
	if(is.null(min_overlap)) {
	
		min_overlap <- max(mean(reads$length) - 
					quantile(diff(reads$start[seq(from = 1, to = n_reads, length.out = n_reads/2)]), p = 0.95), 25)
	}

	#TODO only if verbose=TRUE + nicer output use data.frame + rounding
	cat("Generating configuration with the following parameters:\n\n")
	cat(paste(" bisulphite conversion rate (fixed) =", 0.95,"\n"))
	cat(paste(" inappropriate bisulphite conversion rate (fixed) =", 0.05,"\n"))
	cat(paste(" Reference prior (fixed) =", 0.99,"\n"))
	cat(paste(" Min overlap length (computed) =", as.integer(min_overlap),"\n"))
	cat(paste(" Number of reads =", n_reads,"\n"))
	
	config <- epiG.algorithm.config(
			model = model,
			ref_prior = .99,
			min_overlap_length = min_overlap,
			reads_hard_limit = n_reads + 1000,
			chunk_size = n_reads,
			ref.file = ref_file,
			alt.file = alt_file,
			...
	)
	
	config <- add_run_configuration(config, bam_file, chr, start, end)
	
	return(config)
}	



#' Create a epiG configuration
#' 
#' @param ref.file 
#' @param alt.file 
#' @param max_iterations 
#' @param prior 
#' @param model 
#' @param sequence_quality_adjust 
#' @param haplo_prior modeled using a geometric distribution
#' @param ref_prior 
#' @param min_overlap_length 
#' @param chunk_size 
#' @param chunk_method 
#' @param reads_hard_limit 
#' @param verbose 
#' @return epiG configuration
#' 
#' @author martin
#' @export
epiG.algorithm.config <- function(
		model, 
		ref.file, 
		alt.file, 
		max_iterations = 3e5, 
		ref_prior = 0.9, 
		min_overlap_length = 50,
		min_CG_count = 1,
		min_HCGD_count = 0,
		min_DGCH_count = 0,
		max_stages = 1,
		structual_prior_scale = 1,
		chunk_size = 5000, 
		chunk_method = "reads", 
		reads_hard_limit = 7500,
		use_paired_reads = FALSE,
		NOMEseq_mode = FALSE,
		verbose = TRUE) {
	
	#TODO check config valid
	# 1) chunk_size < reads_hard_limit
	
	config <- list()
	
	config$ref.filename <- ref.file
	
	config$alt.filename <- alt.file
	
	config$max_iterations <- as.integer(max_iterations)
	
	config$fwd_model <- model$fwd
	config$rev_model <- model$rev
		
	config$fwd_DGCH_model <- model$fwd_DGCH
	config$rev_DGCH_model <- model$rev_DGCH
	
	config$fwd_HCGD_model <- model$fwd_HCGD 
	config$rev_HCGD_model <- model$rev_HCGD
	
	config$fwd_C_G_model <- model$fwd_C_G 
	config$rev_C_G_model <- model$rev_C_G
	
	config$chunk_size <- chunk_size
	
	config$chunk.method <- chunk_method
	
	config$reads_hard_limit <- as.integer(reads_hard_limit)
	
	config$ref_prior <- ref_prior
	
	config$min_overlap_length <- as.integer(min_overlap_length)
	
	config$min_CG_count <- as.integer(min_CG_count)
	config$min_HCGD_count <- as.integer(min_HCGD_count)
	config$min_DGCH_count <- as.integer(min_DGCH_count)
	
	config$max_stages <- as.integer(max_stages)
	
	config$structual_prior_scale <- structual_prior_scale
		
	config$use_paired_reads <- use_paired_reads
		
	config$NOMEseq_mode <- NOMEseq_mode
		
	config$verbose <- verbose
	
	config$fwd_GpC_model <- list(matrix(0))
	config$rev_GpC_model <- list(matrix(0))
	
	if(NOMEseq_mode) {
		config$fwd_GpC_model <- model$fwd_GpC
		config$rev_GpC_model <- model$rev_GpC
	} 
	
	class(config) <- "epiG.config"
	
	return(config)
}

add_run_configuration <- function(config, filename, refname, start, end) {
	
	config$filename <- filename
	config$start <- start
	config$end <- end
	config$refname <- refname
	
	return(config)
}

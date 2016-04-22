

.create_error_distributions <- function(bisulfite_rate, bisulfite_inap_rate) {
	
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

#' create_bisulfite_model
#' 
#' @param bisulfite_rates 
#' @param bisulfite_inap_rate 
#' @param Lmax 
#' @return ...
#' 
#' @author Martin Vincent
#' @export
create_bisulfite_model <- function(bisulfite_rate, bisulfite_inap_rate, Lmax) {
	
	bisulfite_rates <- rep(bisulfite_rate, Lmax)
	bisulfite_inap_rates <- rep(bisulfite_inap_rate, Lmax)
	
	model <- list()
	model$fwd <-lapply(1:Lmax, function(i) create_error_distributions(bisulfite_rates[i], bisulfite_inap_rates[i])$fwd) 
	model$rev <-lapply(1:Lmax, function(i) create_error_distributions(bisulfite_rates[i], bisulfite_inap_rates[i])$rev) 
	
	return(model)
}


#' auto_config
#' 
#' @param ref_file 
#' @param alt_file 
#' @param NOMEseq 
#' @param bam_file 
#' @param chr 
#' @param start 
#' @param end 
#' @param use_paired_reads 
#' @param min_overlap 
#' @param bisulfite_rate 
#' @param bisulfite_inap_rate 
#' @param quality_threshold 
#' @param ... 
#' 
#' @return ...
#' 
#' @author Martin Vincent
#' @export
auto_config <- function(
		ref_file, 
		alt_file, 
		NOMEseq = FALSE,
		bam_file,
		chr = NULL, 
		start = NULL, 
		end = NULL, 
		use_paired_reads = NULL,
		min_overlap = NULL,
		bisulfite_rate = .95, 
		bisulfite_inap_rate = 0.05, 
		quality_threshold = 0.020,
		...) {

	if( ! is.null(chr) && length(chr) > 1) {
		
		confs <- lapply(1:length(chr), function(i) auto_config(
							ref_file, 
							alt_file, 
							NOMEseq,
							bam_file, 
							chr[i], 
							start[i], 
							end[i], 
							use_paired_reads,
							min_overlap,
							bisulfite_rate = bisulfite_rate,
							bisulfite_inap_rate = bisulfite_inap_rate,
							quality_threshold = quality_threshold,
							...))
		
		return(confs)
	}
	
	# Fetch information about file
	
	reads <- fetch_reads_info(bam_file, chr, start, end)
	n_reads <- nrow(reads)
	
	if(n_reads == 0) {
		warning("No reads found") #TODO more info
		return (NULL)
	}
	
	# Paired design
	if(is.null(use_paired_reads))  {
		if(max(table(reads$name[1:1000])) == 2) {
			use_paired_reads <- TRUE
		} else {
			use_paired_reads <- FALSE
		}
	}

	### Create bisulfite model
	#TODO auto detimen bisulfite rates
	model <- create_bisulfite_model(
			bisulfite_rate = bisulfite_rate, 
			bisulfite_inap_rate = bisulfite_inap_rate, 
			Lmax = max(reads$length))

	if(NOMEseq) {
		
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
		
		min_CG_count <- 0
		min_HCGD_count <- 0
		min_DGCH_count <- 2
	
		structual_prior_scale <- 1
	
		if(is.null(min_overlap)) {
			min_overlap <- 40
		}
		
	} else {
		
		if(use_paired_reads) {
			min_CG_count <- 2
			
			if(is.null(min_overlap)) {
				min_overlap <- 50
			}
			
		} else {
			min_CG_count <- 1
			
			if(is.null(min_overlap)) {
				min_overlap <- 40
			}
			
		}
		
		min_HCGD_count <- 0
		min_DGCH_count <- 0
		
		structual_prior_scale <- 1
	}
	

	#TODO only if verbose=TRUE + nicer output use data.frame + rounding
	cat("Generating configuration with the following parameters:\n\n")
	cat(paste(" bisulphite conversion rate =", bisulfite_rate,"\n"))
	cat(paste(" inappropriate bisulphite conversion rate =", bisulfite_inap_rate,"\n"))
	cat(paste(" Min overlap length (computed) =", as.integer(min_overlap),"\n"))
	cat(paste(" Number of reads =", n_reads,"\n"))
	cat(paste(" Use paired reads =", use_paired_reads,"\n"))
	
	#TODO max chunk_size
	
	config <- epiG_config(
			model = model,
			reads_hard_limit = n_reads + 1000,
			chunk_size = n_reads,
			ref.file = ref_file,
			alt.file = alt_file,
			min_overlap_length = min_overlap,
			min_CG_count = min_CG_count,
			min_HCGD_count = min_HCGD_count,
			min_DGCH_count = min_DGCH_count,
			structual_prior_scale = structual_prior_scale,
			use_paired_reads = use_paired_reads,
			split_mode = NOMEseq,
			quality_threshold = quality_threshold,
			...
	)
	
	# Add run configuration
	if( ! is.null(chr)) {
		
		if(length(chr) != length(start) || length(chr) != length(end) || length(start) != length(end)) {
			stop("length of chr, start and end must be equal")
		}
	
		config <- add_run_configuration(config, bam_file, chr, start, end)
	}
	
	return(config)
}	



#' epiG_config
#' 
#' @param model 
#' @param ref_file 
#' @param alt_file 
#' @param max_iterations 
#' @param ref_prior 
#' @param min_overlap_length 
#' @param min_CG_count 
#' @param min_HCGD_count 
#' @param min_DGCH_count 
#' @param margin 
#' @param max_stages 
#' @param structual_prior_scale 
#' @param chunk_size 
#' @param chunk_method 
#' @param reads_hard_limit 
#' @param quality_threshold 
#' @param use_paired_reads 
#' @param split_mode 
#' @param verbose 
#' 
#' @return ...
#' 
#' @author Martin Vincent
#' @export
epiG_config <- function(
		model, 
		ref_file, 
		alt_file, 
		max_iterations = 1e5, 
		ref_prior = 1-1e-4, 
		min_overlap_length = 50,
		min_CG_count = 1,
		min_HCGD_count = 0,
		min_DGCH_count = 0,
		margin = 5,
		max_stages = 1,
		structual_prior_scale = 1,
		chunk_size = 5000, 
		chunk_method = "reads", 
		reads_hard_limit = 7500,
		quality_threshold = 0.020,
		use_paired_reads = FALSE,
		split_mode = FALSE,
		verbose = TRUE) {
	
	#TODO check config valid
	# 1) chunk_size < reads_hard_limit
	
	config <- list()
	
	config$ref_filename <- ref_file
	
	config$alt_filename <- alt_file
	
	config$max_iterations <- as.integer(max_iterations)
	
	config$fwd_model <- model$fwd
	config$rev_model <- model$rev
	
	if(split_mode) {
		
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
	
	config$reads_hard_limit <- as.integer(reads_hard_limit)
	
	config$quality_threshold <- quality_threshold
	
	config$ref_prior <- ref_prior
	
	config$min_overlap_length <- as.integer(min_overlap_length)
	
	config$min_CG_count <- as.integer(min_CG_count)
	config$min_HCGD_count <- as.integer(min_HCGD_count)
	config$min_DGCH_count <- as.integer(min_DGCH_count)
	config$margin <- as.integer(margin)
	
	config$max_stages <- as.integer(max_stages)
	
	config$structual_prior_scale <- structual_prior_scale
		
	config$use_paired_reads <- use_paired_reads
		
	config$split_mode <- split_mode
		
	config$verbose <- verbose
	
	class(config) <- "epiG.config"
	
	return(config)
}

#' add_run_configuration
#' 
#' @param config 
#' @param filename 
#' @param refname 
#' @param start 
#' @param end 
#' 
#' @author Martin Vincent
#' @export
add_run_configuration <- function(config, filename, refname, start, end) {
	
	config$filename <- filename
	config$start <- start
	config$end <- end
	config$refname <- refname
	
	return(config)
}

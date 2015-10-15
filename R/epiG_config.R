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
create_error_distributions <- function(bisulfite_rate = 0.94, bisulfite_inap_rate = 0.06) {
	
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

create_pcr_model <- function(rate = 0.01) 1/3*rate + (1-4/3*rate)*diag(4)

#' exp_decay
#' 
#' @param lambda 
#' @param Lmax 
#' @param x 
#' @return function values
#' 
#' @author martin
#' @export
exp_decay <- function(lambda = 0.1, Lmax = 100, x = 0:(Lmax-1)) {
	
	c <- (1-exp(-lambda))/(1-exp(-lambda*(Lmax)))
	
	return(c*exp(-lambda*x))
}

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
create_bisulfite_model <- function(bisulfite_rates = 0.94, bisulfite_inap_rate = 0.06, lambda = 0.1, Lmax = 100) {
	
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
		ref.file, 
		alt.file, 
		max_iterations = 1e5, 
		model = create_bisulfite_model(), 
		sequence_quality_adjust = 0.1, 
		haplo_prior = c(10, 2), 
		ref_prior = 0.9, 
		min_overlap_length = 1, 
		chunk_size = 5000, 
		chunk_method = "reads", 
		reads_hard_limit = 7500,
		use_paired_reads = FALSE,
		verbose = TRUE) {
	
	#TODO check config valid
	# 1) chunk_size < reads_hard_limit
	
	config <- list()
	
	config$ref.filename <- ref.file
	
	config$alt.filename <- alt.file
	
	config$max_iterations <- as.integer(max_iterations)
	
	config$fwd_model <- model$fwd
	
	config$rev_model <- model$rev
	
	config$sequence_quality_adjust <- sequence_quality_adjust
	
	config$chunk_size <- chunk_size
	
	config$chunk.method <- chunk_method
	
	config$reads_hard_limit <- as.integer(reads_hard_limit)
	
	config$ref_prior <- ref_prior
	
	config$min_overlap_length <- as.integer(min_overlap_length)
	
	config$use_paired_reads <- use_paired_reads
		
	config$verbose <- verbose
	
	# Compute log haplochain prior values
	# TODO move out in function
	c <- haplo_prior[1]
	N <- reads_hard_limit
	n <- 1:N
	r <- haplo_prior[2]
	hrel <- c*n*log(1+r/n)/(N*(log(1+r/N)))
	
	config$log_haplo_prior <- cumsum(hrel)
				
	class(config) <- "epiG.config"
	
	return(config)
}

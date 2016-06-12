#
#     Description of this R script:
#   TODO licens
#
#     Intended for use with R.
#     Copyright (C) 2014 Martin Vincent
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

#' @title Fit an epiG epigenotype model
#' @description
#' Fit an epiG epigenotype model
#' 
#' @param config epiG configuration
#' @param max_threads maximal number of threads to use
#' 
#' @return fitted model
#' 
#' @importFrom utils packageVersion
#' @export
#' @useDynLib epiG r_epiG_haplo_fit_filename
#' @useDynLib epiG r_epiG_haplo_fit_filename_chunks
#' @useDynLib epiG r_epiG_compute_chunk_positions
#' @author Martin Vincent
#'
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
#' # Run epiG
#' fit <- epiG(max_threads = 2, config = config)
#'
#' # Fetch additional information
#' fit <- fetch_reads(fit)
#' fit <- fetch_ref(fit)
#' fit <- fetch_alt(fit)
#'
#' # Information about fitted model
#' fit 
#' 
#' # Information about haplotype chains
#' chain_info(fit)
epiG <- function(config, max_threads = 2L) {
	
	start <- config$start
	end <- config$end
	refname <- config$refname
	filename <- config$filename
	
	if( ! is.character(refname)) {
		stop("refname must be a character string")
	}
	
	if(start >= end) {
		stop("end must be larger than start")
	}
	
	if(config$chunk_method == "none") {
		res <- .Call(r_epiG_haplo_fit_filename, filename, refname,  as.integer(start),  as.integer(end), as.integer(max_threads), as.integer(end-start+1), config)
		
	} else if(config$chunk_method == "bases") {
		res <- .Call(r_epiG_haplo_fit_filename, filename, refname,  as.integer(start),  as.integer(end), as.integer(max_threads), as.integer(config$chunk_size), config)
		
	} else if(config$chunk_method == "reads")  {
		s <- .Call(r_epiG_compute_chunk_positions, filename, refname, as.integer(start), as.integer(end), as.integer(config$chunk_size))
		
		if(length(s) >= 2) {
			chunks_start <- s[1:(length(s)-1)]
			chunks_end <- s[2:length(s)]-1L
			chunks_start[1] <- start
			refnames <- as.list(rep(refname, length(chunks_start)))
			configs <- replicate(length(chunks_start), config, simplify = FALSE) #TODO call epiG.chunks 
			
			res <- .Call(r_epiG_haplo_fit_filename_chunks, filename, refnames, as.integer(chunks_start), as.integer(chunks_end), as.integer(max_threads), configs)
		
		} else {
			res <- .Call(r_epiG_haplo_fit_filename, filename, refname,  as.integer(start),  as.integer(end), as.integer(max_threads), as.integer(end-start+1), config)
		}
		
	} else {
		stop("Unknown chunk method")
	}
	
	n_chunks <- res$number_of_chunks
	
	res.chunks <- list()
	
	for(i in 1:n_chunks) {
		
		res.chunks[[i]] <- list()
		
		res.chunks[[i]]$epiG_version <- packageVersion("epiG")
		res.chunks[[i]]$date <- date()
		
		res.chunks[[i]]$config <- config
		res.chunks[[i]]$filename <- filename
		res.chunks[[i]]$refname <- refname
		
		res.chunks[[i]]$offset <- res$chunks_start[i]
		
		res.chunks[[i]]$read_ids <- lapply(res$read_id[[i]], function(x) x + 1)
		res.chunks[[i]]$read_unique_chunk_ids <- sapply(res$read_unique_chunk_ids[[i]], function(x) x + 1)
		res.chunks[[i]]$read_names <- res$read_names
		res.chunks[[i]]$haplotype$chain <- as.integer(factor(res$haplotype[[i]]))
		res.chunks[[i]]$haplotype$start <- res$chain_start[[i]]
		res.chunks[[i]]$haplotype$end <- res$chain_end[[i]]
		
		res.chunks[[i]]$strands <- factor(res$strands[[i]])
		levels(res.chunks[[i]]$strands) <- c("fwd", "rev")
		
		res.chunks[[i]]$genotypes <-  lapply(res$genotypes[[i]], function(x) as.integer(x + 1))
		res.chunks[[i]]$loglikes <- res$loglikes[[i]]
				
		res.chunks[[i]]$length <- res$chunks_end[i] - res.chunks[[i]]$offset + 1
		
		class(res.chunks[[i]]) <- "epiG"
	}
	
	if(n_chunks == 1) {
		return(res.chunks[[1]])
	}	
	
	class(res.chunks) <- c("epiG", "chunks")
	return(res.chunks)
}

#' @title Fit epiG epigenotype models
#' @description
#' Fit epiG epigenotype models
#' 
#' Fit an epiG epigenotype model for each configuration in the list configs.
#' 
#' @param configs list of epiG configurations
#' @param max_threads maximal number of threads to use
#' 
#' @return list of fitted models
#' 
#' @importFrom utils packageVersion
#' @export
#' @useDynLib epiG r_epiG_haplo_fit_filename_chunks
#' @author Martin Vincent
epiG_chunks <- function(configs, max_threads = 8L) {
	
	if(class(configs) == "epiG.config") {
		stop("configs should be a list of configurations")
	}
	
	if(sum(sapply(configs, is.null)) > 0) {
		warning("removing null configurations")
		configs <- configs[ ! sapply(configs, is.null)]
	}
	
	filename = configs[[1]]$filename
	
	chunks_start <- sapply(configs, function(x) x$start)
	chunks_end <- sapply(configs, function(x) x$end)
	refnames <- lapply(configs, function(x) x$refname)

	
	res <- .Call(r_epiG_haplo_fit_filename_chunks, filename, refnames, as.integer(chunks_start),  as.integer(chunks_end), as.integer(max_threads), configs)
	
	n_chunks <- res$number_of_chunks
	
	res.chunks <- list()
			
	for(i in 1:n_chunks) {

		res.chunks[[i]] <- list()
		
		res.chunks[[i]]$epiG_version <- packageVersion("epiG")
		res.chunks[[i]]$date <- date()
		
		res.chunks[[i]]$config <- configs[[i]]
		res.chunks[[i]]$filename <- filename
		res.chunks[[i]]$refname <- refnames[[i]]
		
		res.chunks[[i]]$offset <- res$chunks_start[i]
		
		res.chunks[[i]]$read_ids <- lapply(res$read_id[[i]], function(x) x + 1)
		res.chunks[[i]]$read_unique_chunk_ids <- sapply(res$read_unique_chunk_ids[[i]], function(x) x + 1)
		res.chunks[[i]]$read_names <- res$read_names
		
		res.chunks[[i]]$haplotype$chain <- as.integer(factor(res$haplotype[[i]]))
		res.chunks[[i]]$haplotype$start <- res$chain_start[[i]]
		res.chunks[[i]]$haplotype$end <- res$chain_end[[i]]
		
		res.chunks[[i]]$strands <- factor(res$strands[[i]])
		levels(res.chunks[[i]]$strands) <- c("fwd", "rev")
		
		res.chunks[[i]]$genotypes <-  lapply(res$genotypes[[i]], function(x) as.integer(x + 1))
		res.chunks[[i]]$loglikes <- res$loglikes[[i]]
		
		res.chunks[[i]]$length <- as.integer(res$chunks_end[i] - res.chunks[[i]]$offset + 1)
		
		class(res.chunks[[i]]) <- "epiG"

	}
	
	if(n_chunks == 1) {
		return(res.chunks[[1]])
	}	
	
	class(res.chunks) <- c("epiG", "chunks")
	return(res.chunks)
}

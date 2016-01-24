#
#     Description of this R script:
#     TODO
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

#' Fit an epiG model
#' 
#' @param filename BAM file with aligned and ordered reads
#' @param refGenom_filename
#' @param altGenom_filename 
#' @param refname 
#' @param start 
#' @param end 
#' @param max_threads 
#' @param config 
#' @return fitted model
#' 
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_haplo_fit_filename
#' @useDynLib epiG r_epiG_haplo_fit_filename_chunks
epiG <- function(config, max_threads = 8L) {
	
	refGenom_filename = config$ref.filename
	altGenom_filename = config$alt.filename
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
	
	if(config$chunk.method == "none") {
		res <- .Call(r_epiG_haplo_fit_filename, filename, refGenom_filename, altGenom_filename, refname,  as.integer(start),  as.integer(end), as.integer(max_threads), as.integer(end-start+1), config)
		
	} else if(config$chunk.method == "bases") {
		res <- .Call(r_epiG_haplo_fit_filename, filename, refGenom_filename, altGenom_filename, refname,  as.integer(start),  as.integer(end), as.integer(max_threads), as.integer(config$chunk_size), config)
		
	} else if(config$chunk.method == "reads")  {
		s <- compute_chunk_positions(filename, refname, as.integer(start), as.integer(end), as.integer(config$chunk_size))
		
		if(length(s) >= 2) {
			chunks_start <- s[1:(length(s)-1)]
			chunks_end <- s[2:length(s)]-1L
			chunks_start[1] <- start
			refnames <- as.list(rep(refname, length(chunks_start)))
			configs <- replicate(length(chunks_start), config, simplify = FALSE) #TODO call epiG.chunks 
			
			res <- .Call(r_epiG_haplo_fit_filename_chunks, filename, refGenom_filename, altGenom_filename, refnames, as.integer(chunks_start), as.integer(chunks_end), as.integer(max_threads), configs)
		
		} else {
			res <- .Call(r_epiG_haplo_fit_filename, filename, refGenom_filename, altGenom_filename, refname,  as.integer(start),  as.integer(end), as.integer(max_threads), as.integer(end-start+1), config)
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

#' epiG.chunks
#' 
#' @param filename 
#' @param refGenom_filename 
#' @param altGenom_filename 
#' @param refnames 
#' @param chunks_start 
#' @param chunks_end 
#' @param max_threads 
#' @param config 
#' @return fitted models
#' 
#' @author martin
#' @export
#' @useDynLib epiG r_epiG_haplo_fit_filename_chunks
epiG.chunks <- function(configs, max_threads = 8L) {
	
	if(class(configs) == "epiG.config") {
		stop("configs should be a list of configurations")
	}
	
	refGenom_filename = configs[[1]]$ref.filename
	altGenom_filename = configs[[1]]$alt.filename
	filename = configs[[1]]$filename
	
	chunks_start <- sapply(configs, function(x) x$start)
	chunks_end <- sapply(configs, function(x) x$end)
	refnames <- lapply(configs, function(x) x$refname)

	
	res <- .Call(r_epiG_haplo_fit_filename_chunks, filename, refGenom_filename, altGenom_filename, refnames,  as.integer(chunks_start),  as.integer(chunks_end), as.integer(max_threads), configs)
	
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

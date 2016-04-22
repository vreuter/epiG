#
#     Description of this R script:
#     TODO
#
#     Intended for use with R.
#     Copyright (C) 2013 Martin Vincent
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#

#' fetch_read_info
#' 
#' @param filename 
#' @param refname 
#' @param start 
#' @param end 
#' @return info 
#' 
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_fetch_reads_info
fetch_read_info <- function(filename, refname, start, end) {

        res <- .Call(r_epiG_fetch_reads_info, filename, refname, as.integer(start), as.integer(end))

        info <- list()
        info$start <- as.integer(res$position)
        info$end <- as.integer(res$position + res$length - 1)
        info$length <- as.integer(res$length)
        info$nread <- length(res$position)

		# Return data.frame
		data.frame(
				name = res$names,
				start = as.integer(res$position),
				end = as.integer(res$position + res$length - 1),
				length = as.integer(res$length)
		)	
}

#' fetch_read_count
#' 
#' @param filename 
#' @param refname 
#' @param start 
#' @param end 
#' @return info 
#' 
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_fetch_reads_info
fetch_read_count <- function(filename, refname, start, end) {
	
	res <- .Call(r_epiG_fetch_read_count, filename, refname, as.integer(start), as.integer(end))
	
	data.frame(
		nreads = res$nreads,
		bp = res$bp.count
	)
}

#' fetch_reads_raw
#' 
#' @param filename 
#' @param refname 
#' @param start 
#' @param end 
#' @return reads
#' 
#' @author Martin Vincent
#' @export
fetch_reads_raw <- function(filename, refname, start, end) {
	
	reads <- .Call(r_epiG_fetch_reads_raw, filename, refname, as.integer(start), as.integer(end))
		
	class(reads) <- "epiG_reads"
	
	return(reads)
}

#' fetch_reads
#' 
#' @param object 
#' @return epiG model
#' 
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_fetch_reads
fetch_reads <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		reads <- .Call(r_epiG_fetch_reads, object$filename, object$refname, start(object), end(object), object$config$quality_threshold)

		if(length(reads$lengths) == 0) {
			
			object$reads <- reads
			
		} else {
		
			reads <- lapply(reads, function(x) x[object$read_unique_chunk_id])
			reads$positions <- reads$positions - object$offset
			object$reads <- reads
		}
	
	} else if(paste(class(object), collapse = ".") == "epiG.chunks") {
		
		object <- lapply(object, fetch_reads)
		class(object) <- c("epiG", "chunks")
		
	} else {
		stop("Unknown object -- object must be a epiG class")
	
	}
	
	return(object)
}

#' fetch_ref
#' 
#' @param object 
#' @return epiG model
#' 
#' @author Martin Vincent
#' @export
fetch_ref <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		object$ref <- read_fasta(object$config$ref.filename, object$refname, start(object), length(object))
		
	} else if(paste(class(object), collapse = ".") == "epiG.chunks") {
		object <- lapply(object, function(x) fetch_ref(x))
		class(object) <- c("epiG", "chunks")
		
	} else {
		stop("Unknown object -- object must be a epiG class")
	
	}
	
	return(object)
}

#' Read fasta
#' 
#' @param filename 
#' @param refname 
#' @param start 
#' @param len 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_read_fasta
read_fasta <- function(filename, refname, start, len) {
	return(.Call(r_epiG_read_fasta, 
					as.character(filename), 
					as.character(refname), 
					as.integer(start), 
					as.integer(len)))
}

#' fetch_alt
#' 
#' @param object 
#' @return epiG model
#' 
#' @author Martin Vincent
#' @export
fetch_alt <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		object$alt <- read_fasta(object$config$alt.filename, object$refname, start(object), length(object))
	
	} else if(paste(class(object), collapse = ".") == "epiG.chunks") {
		object <- lapply(object, function(x) fetch_alt(x))
		class(object) <- c("epiG", "chunks")
		
	} else {
		stop("Unknown object -- object must be a epiG class")
	
	}

	return(object)
}

#' header_info
#'  
#' @param filename 
#' 
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_fetch_header
header_info <- function(filename) {
 	
	tmp <- .Call(r_epiG_fetch_header, filename)
	
	data.frame(ref = tmp$refname, length = tmp$lengths, stringsAsFactors = FALSE)	
}

#' file_info
#' 
#' @param filename 
#' 
#' @author Martin Vincent
#' @export
file_info <- function(filename) {
	
	info <- fetch_header(filename)
	
	info$nreads <- NA
	info$mean_read_length <- NA
	
	for(i in 1:nrow(info)) {
		tmp <- fetch_reads_info(filename, info$ref[i], 0, info$length[i])
		info$nreads[i] <- tmp$nread
		info$mean_read_length[i] <- mean(tmp$length)
	}
	
	return(info)
}

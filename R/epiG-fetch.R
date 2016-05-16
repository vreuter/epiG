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

#' @title Fetch information about reads
#' @description 
#' Fetch information about reads overlapping the specified region
#'
#' @param file path to bam file
#' @param refname reference name
#' @param start start of region 
#' @param end end of region
#' @return data.frame with information abut reads. #TODO cols 
#'
#' @examples
#' # Retrieve paths to raw data files
#' bam_file <- system.file("extdata", "GNAS_small.bam", package="epiG")
#' 
#' fetch_read_info(bam_file, "chr20", 57400000, 57400000 + 100)
#'
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_fetch_reads_info
fetch_read_info <- function(file, refname, start, end) {

        res <- .Call(r_epiG_fetch_reads_info, file, refname, as.integer(start), as.integer(end))

		# Return data.frame
		data.frame(
				name = res$names,
				start = as.integer(res$position),
				end = as.integer(res$position + res$length - 1),
				length = as.integer(res$length)
		)	
}

#' @title Count reads overlapping specified region
#' @description 
#' 
#' @param file path to bam file
#' @param refname reference name
#' @param start start of region
#' @param end end of region
#' @return number of reads and total bps in reads overlapping region
#'
#' @examples
#' # Retrieve paths to raw data files
#' bam_file <- system.file("extdata", "GNAS_small.bam", package="epiG")
#' 
#' read_count(bam_file, "chr20", 57400000, 57400000 + 100)
#'
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_fetch_read_count
read_count <- function(file, refname, start, end) {
	
	res <- .Call(r_epiG_fetch_read_count, file, refname, as.integer(start), as.integer(end))
	
	data.frame(
		nreads = res$nreads,
		bp = res$bp.count
	)
}

#' @title Load reads 
#' @description 
#' 
#' @param file path to bam file
#' @param refname reference name
#' @param start start of region
#' @param end end of region
#' @param quality_threshold quality threshold
#' @param raw_quality_scores if TRUE raw quality score will be returned
#' @return TODO
#' 
#' @examples
#' # Retrieve paths to raw data files
#' bam_file <- system.file("extdata", "GNAS_small.bam", package="epiG")
#' 
#' load_reads(bam_file, "chr20", 57400000, 57400000 + 100)
#' 
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_fetch_reads
load_reads <- function(file, refname, start, end, quality_threshold = 1, raw_quality_scores = FALSE) {
	
	reads <- .Call(r_epiG_fetch_reads, file, refname, as.integer(start), as.integer(end), quality_threshold, ! raw_quality_scores)
		
	class(reads) <- "epiG.reads"
	
	return(reads)
}

#' @title Fetch reads
#' @description 
#'
#' Reads will be loaded and include in the epiG object
#' 
#' @param object epiG epigeotype model 
#'
#' @return model with reads included (this may increase the memory use) 
#' 
#' @author Martin Vincent
#' @export
fetch_reads <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		reads <- load_reads(object$filename, object$refname, start(object), end(object), object$config$quality_threshold, FALSE)

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

#' @title Fetch reference genom
#' @description 
#'
#' Load reference genom and include it in epiG object. 
#' 
#' @param object epiG epigenotype model
#' @return epiG epigenotype model
#' 
#' @author Martin Vincent
#' @export
fetch_ref <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		object$ref <- read_fasta(
						file = object$config$ref_filename, 
						refname = object$refname, 
						start = start(object),
						len = length(object), 
						offset = object$config$ref_offset)
		
	} else if(paste(class(object), collapse = ".") == "epiG.chunks") {
		object <- lapply(object, function(x) fetch_ref(x))
		class(object) <- c("epiG", "chunks")
		
	} else {
		stop("Unknown object -- object must be a epiG class")
	
	}
	
	return(object)
}

#' @title Read fasta file
#' @description 
#' 
#' @param file path to fasta file
#' @param refname referance name
#' @param start start of region
#' @param len length of region
#' @param offset file offset (position of first base in file)
#' @return TODO
#' 
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_read_fasta
read_fasta <- function(file, refname, start, len, offset = 0) {
	return(.Call(r_epiG_read_fasta, 
					as.character(file), 
					as.integer(offset),
					as.character(refname), 
					as.integer(start), 
					as.integer(len)))
}

#' @title fetch_alt
#' @description 
#' 
#' @param object 
#' @return epiG model
#' 
#' @author Martin Vincent
#' @export
fetch_alt <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		object$alt <- read_fasta(object$config$alt_filename, object$refname, start(object), length(object))
	
	} else if(paste(class(object), collapse = ".") == "epiG.chunks") {
		object <- lapply(object, function(x) fetch_alt(x))
		class(object) <- c("epiG", "chunks")
		
	} else {
		stop("Unknown object -- object must be a epiG class")
	
	}

	return(object)
}

#' @title Fetch header of bam file
#' @description 
#'  
#' @param file path to bam file
#' 
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_fetch_header
header_info <- function(file) {
 	
	tmp <- .Call(r_epiG_fetch_header, file)
	
	data.frame(ref = tmp$refname, length = tmp$lengths, stringsAsFactors = FALSE)	
}

#' @title Fetch information about bam file
#' @description 
#' 
#' @param file path to bam file
#' 
#' @author Martin Vincent
#' @export
file_info <- function(file) {
	
	info <- header_info(file)
	
	info$nreads <- NA
	info$mean_read_length <- NA
	
	for(i in 1:nrow(info)) {
		tmp <- fetch_read_info(file, info$ref[i], 0, info$length[i])
		info$nreads[i] <- nrow(tmp)
		
		if(nrow(tmp) > 0) {
			info$mean_read_length[i] <- mean(tmp$length)
		}
	}
	
	return(info)
}

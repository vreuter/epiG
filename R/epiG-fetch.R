#
#     Description of this R script:
#     R scripts for fetching information from bam files
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
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#

#' @title Information About Reads
#' @description 
#' Fetch information about reads overlapping the specified region
#'
#' @param file path to bam file
#' @param refname reference name
#' @param start start of region 
#' @param end end of region
#' @return data.frame with information abut reads. Columns:
#' \itemize{
#'  \item[\code{name}]{name of read}
#'  \item[\code{start}]{start position of read}
#'  \item[\code{end}]{end position of read}
#'  \item[\code{length}]{length of read}
#' }
#'
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_fetch_reads_info
#' @examples
#' # Retrieve paths to raw data files
#' bam_file <- system.file("extdata", "GNAS_small.bam", package="epiG")
#' 
#' fetch_read_info(bam_file, "chr20", 57400000, 57400000 + 100)
#'
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

#' @title Count Reads
#' @description 
#' Count reads overlapping the specified region
#' 
#' @param file path to bam file
#' @param refname reference name
#' @param start start of region
#' @param end end of region
#' @return number of reads and total bps in reads overlapping region
#'
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_fetch_read_count
#'
#' @examples
#' # Retrieve paths to raw data files
#' bam_file <- system.file("extdata", "GNAS_small.bam", package="epiG")
#' 
#' read_count(bam_file, "chr20", 57400000, 57400000 + 100)
read_count <- function(file, refname, start, end) {
	
	res <- .Call(r_epiG_fetch_read_count, file, refname, as.integer(start), as.integer(end))
	
	data.frame(
		nreads = res$nreads,
		bp = res$bp.count
	)
}

#' @title Load Reads 
#' @description 
#' Load the reads overlapping the specified region
#' 
#' @param file path to bam file
#' @param refname reference name
#' @param start start of region
#' @param end end of region
#' @param quality_threshold quality threshold
#' @param raw_quality_scores if TRUE raw quality score will be returned
#' @return a list with the following entries:
#' \itemize{
#'  \item[\code{reads}]{a list of reads (each read represented by a vector of bases)}
#'  \item[\code{quality}]{a list quality scores}
#'  \item[\code{positions}]{a vector of the start positions of the reads}
#'  \item[\code{lengths}]{a vector of the lengths of the reads}
#'  \item[\code{names}]{a vector of the names of the reads}
#' }
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_fetch_reads
#' @examples
#' # Retrieve paths to raw data files
#' bam_file <- system.file("extdata", "GNAS_small.bam", package="epiG")
#' 
#' info <- load_reads(bam_file, "chr20", 57400000, 57400000 + 100)
#' 
#' # Bases, qualities, start position, length and name of first read
#' info$reads[[1]]
#' info$quality[[1]]
#' info$position[1]
#' info$lengths[1]
#' info$names[1]
load_reads <- function(file, refname, start, end, quality_threshold = 1, raw_quality_scores = FALSE) {
	
	reads <- .Call(r_epiG_fetch_reads, file, refname, as.integer(start), as.integer(end), quality_threshold, ! raw_quality_scores)
		
	reads$reads <- lapply(reads$reads, .symbols)

	class(reads) <- "epiG.reads"
	
	return(reads)
}

#' @title Fetch Reads
#' @description
#' Reads will be loaded and include in the epiG object
#' 
#' @param object epiG epigenotype model 
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

#' @title Fetch Reference Genom
#' @description
#' Load reference genome and include it in epiG object. 
#' 
#' @param object epiG epigenotype model
#' @return epiG epigenotype model with reference genome included 
#' 
#' @author Martin Vincent
#' @export
fetch_ref <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		object$ref <- .read_fasta(
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

# Internal read_fast see description for read_fast
.read_fasta <- function(file, refname, start, len, offset) {
	return(.Call(r_epiG_read_fasta, 
					as.character(file), 
					as.integer(offset),
					as.character(refname), 
					as.integer(start), 
					as.integer(len)))
}

#' @title Read FASTA File
#' @description 
#' Read a raw FASTA file
#'
#' @param file path to fasta file
#' @param refname reference name
#' @param start start of region
#' @param len length of region
#' @param offset file offset (position of first base in file, usually offset = 1)
#' @return a vector of length \code{len} containing the bases 
#' 
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_read_fasta
#' @examples 
#' ref_file <- system.file("extdata", "hg19_GNAS.fa", package="epiG")
#'
#' # Specify region
#' chr <- "chr20"
#' start <- 57400000 
#' end <- 57400000 + 100
#' 
#' #Note that usually offset = 1
#' read_fasta(ref_file, chr, start, len = end-start+1, offset = 57380000)
read_fasta <- function(file, refname, start, len, offset = 1) {

	fas <- .read_fasta(file, refname, start, len, offset - 1L)
	fas <- sapply(fas, .symbols)

	return(fas)
}



#' @title Fetch Alternative Nucleotides
#' @description 
#' Load alternative nucleotides and include it in epiG object. 
#' 
#' @param object a epiG model
#' @return epiG model
#' 
#' @author Martin Vincent
#' @export
fetch_alt <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		object$alt <- .read_fasta(object$config$alt_filename, object$refname, start(object), length(object), object$config$alt_offset)
	
	} else if(paste(class(object), collapse = ".") == "epiG.chunks") {
		object <- lapply(object, function(x) fetch_alt(x))
		class(object) <- c("epiG", "chunks")
		
	} else {
		stop("Unknown object -- object must be a epiG class")
	
	}

	return(object)
}

#' @title Fetch Bam Header
#' @description 
#' Load bam file header
#'  
#' @param file path to bam file
#' @return a list of refnames and lengths associated with the file
#' 
#' @author Martin Vincent
#' @export
#' @useDynLib epiG r_epiG_fetch_header
#' @examples 
#' # Retrieve paths to raw data files
#' bam_file <- system.file("extdata", "GNAS_small.bam", package="epiG")
#'
#' header_info(bam_file)
header_info <- function(file) {
 	
	tmp <- .Call(r_epiG_fetch_header, file)
	
	data.frame(ref = tmp$refname, length = tmp$lengths, stringsAsFactors = FALSE)	
}

#' @title Fetch Bam File Information
#' @description 
#' Fetch information about bam file
#' 
#' @param file path to bam file
#' @return a data.frame with the following columns:
#' \itemize{
#'  \item[\code{ref}]{refname}
#'  \item[\code{length}]{length of ref}
#'  \item[\code{nreds}]{number of reads assigned to refname}
#'  \item[\code{mean_read_length}]{the mean read length of reads assigned to refname}
#' }
#' 
#' @author Martin Vincent
#' @export
#' @examples 
#' # Retrieve paths to raw data files
#' bam_file <- system.file("extdata", "GNAS_small.bam", package="epiG")
#' 
#' file_info(bam_file)
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

#
#     Description of this R script:
#     TODO
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

#' Search for pattern in integer vector
#' 
#' @param pattern integer vector 
#' @param x integer vector to search in
#' @return postion of pattern in x 
#'
#' @examples
#' TODO
#'
#' @author Martin Vincent
#' @useDynLib epiG r_epiG_locate
#' @export
vector_search <- function(pattern, x) {
	
	#TODO check input
	
	pattern <- as.integer(pattern)
	x <- as.integer(x)
		
	res <- .Call(r_epiG_locate, pattern, x)
	
	res <- res + 1L
	
	return(res)
}

#' locate C positions
#'
#' @param object epiG object
#' @return positions of C in object
#' 
#' @examples
#' todo
#'
#' @author Martin Vincent
#' @export
locate_C <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No ref genom found")
		}
		
		pos <- which(object$ref == 1)
		
		return(pos + object$offset - 1L)	
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(unlist(lapply(object, function(x) locate_C(x))))
	}
	
	stop("Unknown class")
}


#' locate GC positions
#' @param object 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
locate_GC <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No ref genom found")
		}
		
		pos <- vector_search(c(2,1), object$ref)
		
		return(pos + object$offset - 1L)	
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(unlist(lapply(object, function(x) locate_GC(x))))
	}
	
	stop("Unknown class")
}

#' locate CG positions
#' @param object 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
locate_CG <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No ref genom found")
		}
		
		pos <- vector_search(c(1,2), object$ref)
		
		return(pos + object$offset - 1L)	
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(unlist(lapply(object, function(x) locate_CG(x))))
	}
	
	stop("Unknown class")
}

#' locate_DGCH 
#' 
#' locate DGCH (isolated GpC) positions
#' 
#' @param object 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
locate_DGCH <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No ref genom found")
		}
		
		pos <- vector_search(c(2,1), object$ref)
		
		# Remove pos <= 1 
		if(any(pos == 1)) {
			pos <- pos[pos > 1]
		}
		
		# Remove pos > length(object$ref)-2
		if(any(pos > length(object$ref)-2)) {
			pos <- pos[pos < (length(object$ref)-2)]
		}
		
		# Remove GCG
		tmp <- vector_search(2, object$ref[pos+2])
		
		if(length(tmp) > 0) {
			pos <- pos[-tmp]
		}	
		
		# Remove CGC
		tmp <- vector_search(1, object$ref[pos-1])
		
		if(length(tmp) > 0) {
			pos <- pos[-tmp]
		}	
		
		return(pos + object$offset - 1L)	
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(unlist(lapply(object, function(x) locate_DGCH(x))))
	}
	
	stop("Unknown class")
}

#' locate_HCGD
#' 
#' locate HCGD (isolated CpG) positions
#' 
#' @param object 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
locate_HCGD <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No ref genom found")
		}
		
		pos <- vector_search(c(1, 2), object$ref)
		
		# Remove pos = 1 
		if(any(pos == 1)) {
			pos <- pos[pos > 1]
		}
		
		# Remove pos > length(object$ref)-2
		if(any(pos > length(object$ref)-2)) {
			pos <- pos[pos < (length(object$ref)-2)]
		}
		
		# Remove GCG
		tmp <- vector_search(2, object$ref[pos-1])
		
		if(length(tmp) > 0) {
			pos <- pos[-tmp]
		}	
		
		# Remove CGC
		tmp <- vector_search(1, object$ref[pos+2])
		
		if(length(tmp) > 0) {
			pos <- pos[-tmp]
		}	
		
		return(pos + object$offset - 1L)	
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(unlist(lapply(object, function(x) locate_HCGD(x))))
	}
	
	stop("Unknown class")
}

#' Locate positions where at least one chain has e genotype not matching with the reference.
#' 
#' @param object 
#' @return postions of mismatches
#' 
#' @author Martin Vincent
#' @export
locate_mismatch <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No ref genom found")
		}
		
		snp_pos <- as.integer()
		for(i in 1:length(object$genotype)) {
			pos <- (object$haplotype$start[i]):(object$haplotype$end[i]) - object$offset + 1
			sel <- pos >= 1 & pos <= length(object$ref)
			snp_pos <- c(snp_pos, pos[sel][object$ref[pos[sel]] != (object$genotypes[[i]][sel]-1)  %% 4 + 1] + object$offset - 1)
		}
		
		return(unique(snp_pos))
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(unlist(lapply(object, function(x) locate_nonref(x))))
	}
	
	stop("Unknown class")
}

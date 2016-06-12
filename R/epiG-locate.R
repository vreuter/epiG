#
#     Description of this R script:
#     TODO licens
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

#' @title Pattern Search
#' @description
#' Search for pattern in integer vector
#' 
#' @param pattern integer vector 
#' @param x integer vector to search in
#' @return position of pattern in x 
#'
#' @useDynLib epiG r_epiG_locate
#' @export
#' @author Martin Vincent
vector_search <- function(pattern, x) {
	
	#TODO check input
	
	pattern <- as.integer(pattern)
	x <- as.integer(x)
		
	res <- .Call(r_epiG_locate, pattern, x)
	
	res <- res + 1L
	
	return(res)
}

#' @title Locate C
#' @description
#' Locate C positions
#'
#' @param object epiG object
#' @return positions of C in object
#' 
#' @author Martin Vincent
#' @export
#'
#' @examples 
#' data(example)
#' 
#' locate_C(fit)
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

#' @title Locate GpC
#' @description
#' Locate GpC positions in the refrence genom
#' @param object epiG model
#' @return a vector of GpC positions 
#' 
#' @author Martin Vincent
#' @export
#'
#' @examples 
#' data(example)
#' 
#' locate_GC(fit)
locate_GC <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No refrence genom found in model")
		}
		
		pos <- vector_search(c(2,1), object$ref)
		
		return(pos + object$offset - 1L)	
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(unlist(lapply(object, function(x) locate_GC(x))))
	}
	
	stop("Unknown class")
}

#' @title Locate CpG
#' @description
#' Locate CpG positions in the refrence genom
#' @param object epiG model
#' @return a vector of CpG positions
#' 
#' @author Martin Vincent
#' @export
#'
#' @examples 
#' data(example)
#' 
#' locate_CG(fit)
locate_CG <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No refrence genom found in model")
		}
		
		pos <- vector_search(c(1,2), object$ref)
		
		return(pos + object$offset - 1L)	
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(unlist(lapply(object, function(x) locate_CG(x))))
	}
	
	stop("Unknown class")
}

#' @title Locate DGCH
#' @description
#' Locate DGCH (isolated GpC) positions in the refrence genom
#' 
#' @param object epiG model
#' @return a vector of isolated GpC positions
#' 
#' @author Martin Vincent
#' @export
#'
#' @examples 
#' data(example)
#' 
#' locate_DGCH(fit)
locate_DGCH <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No refrence genom found in model")
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

#' @title Locate HCGD 
#' @description
#' locate HCGD (isolated CpG) positions in the refrence genom
#' 
#' @param object epiG model
#' @return a vector of isolated CpG positions
#' 
#' @author Martin Vincent
#' @export
#'
#' @examples 
#' data(example)
#' 
#' locate_HCGD(fit)
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

#' @title Locate Mismatches
#' @description
#' Locate positions where at least one chain has a genotype not matching with the reference.
#' 
#' @param object an epiG model
#' @return a vector of postions of mismatches
#' 
#' @author Martin Vincent
#' @export
#'
#' @examples 
#' data(example)
#' 
#' locate_mismatch(fit)
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
		return(unlist(lapply(object, function(x) locate_mismatch(x))))
	}
	
	stop("Unknown class")
}

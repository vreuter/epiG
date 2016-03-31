# TODO: Add comment
# 
# Author: martin
###############################################################################


#' vector search
#' @param pattern 
#' @param text 
#' @return ??
#' @author Martin Vincent
#' @useDynLib epiG r_epiG_locate
#' @export
vector.search <- function(pattern, text) {
	
	#TODO check input
	
	pattern <- as.integer(pattern)
	text <- as.integer(text)
		
	res <- .Call(r_epiG_locate, pattern, text)
	
	res <- res + 1L
	
	return(res)
}

#' locate GC positions
#' @param object 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
locate.GC <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No ref genom found")
		}
		
		pos <- vector.search(c(2,1), object$ref)
		
		return(pos + object$offset - 1L)	
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(unlist(lapply(object, function(x) locate.GC(x))))
	}
	
	stop("Unknown class")
}

#' locate CpG positions
#' @param object 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
locate.CG <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No ref genom found")
		}
		
		pos <- vector.search(c(1,2), object$ref)
		
		return(pos + object$offset - 1L)	
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(unlist(lapply(object, function(x) locate.CG(x))))
	}
	
	stop("Unknown class")
}

#' locate DGCH positions
#' @param object 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
locate.DGCH <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No ref genom found")
		}
		
		pos <- vector.search(c(2,1), object$ref)
		
		# Remove pos <= 1 
		if(any(pos == 1)) {
			pos <- pos[pos > 1]
		}
		
		# Remove pos > length(object$ref)-2
		if(any(pos > length(object$ref)-2)) {
			pos <- pos[pos < (length(object$ref)-2)]
		}
		
		# Remove GCG
		tmp <- vector.search(2, object$ref[pos+2])
		
		if(length(tmp) > 0) {
			pos <- pos[-tmp]
		}	
		
		# Remove CGC
		tmp <- vector.search(1, object$ref[pos-1])
		
		if(length(tmp) > 0) {
			pos <- pos[-tmp]
		}	
		
		return(pos + object$offset - 1L)	
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(unlist(lapply(object, function(x) locate.DGCH(x))))
	}
	
	stop("Unknown class")
}

#' locate HCGD positions
#' @param object 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
locate.HCGD <- function(object) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(is.null(object[["ref"]])) {
			stop("No ref genom found")
		}
		
		pos <- vector.search(c(1, 2), object$ref)
		
		# Remove pos = 1 
		if(any(pos == 1)) {
			pos <- pos[pos > 1]
		}
		
		# Remove pos > length(object$ref)-2
		if(any(pos > length(object$ref)-2)) {
			pos <- pos[pos < (length(object$ref)-2)]
		}
		
		# Remove GCG
		tmp <- vector.search(2, object$ref[pos-1])
		
		if(length(tmp) > 0) {
			pos <- pos[-tmp]
		}	
		
		# Remove CGC
		tmp <- vector.search(1, object$ref[pos+2])
		
		if(length(tmp) > 0) {
			pos <- pos[-tmp]
		}	
		
		return(pos + object$offset - 1L)	
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(unlist(lapply(object, function(x) locate.HCGD(x))))
	}
	
	stop("Unknown class")
}

#' find positions of mutations 
#' @param object 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
locate.mutations <- function(object) {
	
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
		return(unlist(lapply(object, function(x) locate.mutations(x))))
	}
	
	stop("Unknown class")
}

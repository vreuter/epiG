# TODO: Add comment
# 
# Author: martin
###############################################################################

#' end position
#' 
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
end <- function(object, ... ) UseMethod("end")

#' start position
#' 
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
start <- function(object, ... ) UseMethod("start")

#' start
#' 
#' @param object 
#' @param ... 
#' @return numeric
#' 
#' @author Martin Vincent
#' @method start epiG
#' @export
start.epiG <- function(object, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(object$offset)
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sapply(object, start))
	}
	
	stop("Unknown class")
}



#' end
#' 
#' @param object 
#' @param ...
#' @return numeric
#' 
#' @author Martin Vincent
#' @method end epiG
#' @export
end.epiG <- function(object, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(start(object) + length(object)-1L)
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sapply(object, end))
	}
	
	stop("Unknown class")
}

#' Length of model in base pairs
#' 
#' @param x 
#' @return Length of model in base pairs
#' 
#' @author Martin Vincent
#' @method length epiG
#' @export
length.epiG <- function(x) {
	
	if(paste(class(x), collapse = ".") == "epiG") {
		return(x$length)
	}
	
	if(paste(class(x), collapse = ".") == "epiG.chunks") {
		return(sum(sapply(x, length)))
	}
	
	stop("Unknown class")
}

#' Number of reads in model
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
nread <- function(object, ... ) UseMethod("nread")

#' Number of reads in model
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @method nread epiG
#' @export
nread.epiG <- function(object, ...)  {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(length(unique(unlist(object$read_ids))))
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sapply(object, nread))
	}
	
	stop("Unknown class")
	
}

#' genotype
#' @param object 
#' @param pos 
#' @param remove.meth 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
genotype <- function(object, pos, remove.meth, ... ) UseMethod("genotype")

#' genotype
#' codeing C = 1, G = 2, A = 3, T = 4
#' @param object 
#' @param pos 
#' @param remove.meth 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @method genotype epiG
#' @export
genotype.epiG <- function(object, pos, remove.meth = FALSE, ...) {
	
	if(length(pos) > 1) {
		stop("pos must have length 1")
	}
	
	if(paste(class(object), collapse = ".") == "epiG") {
	
		if(start(object) > pos || end(object) < pos) {
			stop("Position not in range")
		}
	
		collected <- sapply(1:length(object$genotype), 
				function(i) if((pos - object$haplotype$start)[i] >= 0 && (object$haplotype$end[i] - pos) >= 0) 
						object$genotype[[i]][(pos - object$haplotype$start)[i]+1] else NA)
		
		names(collected) <- 1:length(object$genotype)
		
		if(remove.meth) {
			collected[collected %in% c(5,6)] <-  collected[collected %in% c(5,6)] %% 4
		}
		
		return(collected[!is.na(collected)])
	}
	
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		stop("Not yet implemented for chunks")
	}
	
	stop("Unknown class")
	
}

#' genotype
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
loglike <- function(object, pos, ... ) UseMethod("loglike")

#' genotype
#' codeing C = 1, G = 2, A = 3, T = 4
#' @param object 
#' @param pos 
#' @param remove.meth 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @method genotype epiG
#' @export
loglike.epiG <- function(object, pos, ...) {
	
	if(length(pos) > 1) {
		stop("pos must have length 1")
	}
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(start(object) > pos || end(object) < pos) {
			stop("Position not in range")
		}
		
		collected <- sapply(1:length(object$loglike), 
				function(i) if((pos - object$haplotype$start)[i] >= 0 && (object$haplotype$end[i] - pos) >= 0) 
						object$loglike[[i]][(pos - object$haplotype$start)[i]+1,] else rep(NA,6))
		
		colnames(collected) <- 1:length(object$loglike)
		
		rownames(collected) <- c("C", "G", "A", "T", "c", "g")
		
		return(collected[,apply(collected, 2, function(x) all(!is.na(x))), drop = FALSE])
	}
	
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		stop("Not yet implemented for chunks")
	}
	
	stop("Unknown class")
	
}

#' coverage
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
coverage <- function(object, pos = NULL, ... ) UseMethod("coverage")

#' coverage
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @method coverage epiG
#' @export
coverage.epiG <- function(object, pos = NULL, ...) {
	
	if(is.null(pos)) {
		return(sapply(start(object):end(object), function(pos) coverage(object, pos)))
	}
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(start(object) > pos || end(object) < pos) {
			stop("Position not in range")
		}
		
		if(length(object$read_ids) == 0) {
			return(0)
		}
		
		return(length(object$read_ids[[pos - start(object)+1]]))
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		
		stop("Not yet implemented")
		#TODO FIX
		tmp <- sapply(object, function(x) if(start(x) > pos || end(x) < pos) coverage(x) else NA)
		return(tmp[!is.na(tmp)])
	}
	
	stop("Unknown class")
	
}

#' position.info
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
position.info <- function(object, pos, ... ) UseMethod("position.info")

.methylation.status <- function(genotype.code, nfwd, nrev) {
	sapply(1:length(genotype.code), function(i) {
				if(genotype.code[i] %in% c(5,6)) return(TRUE)
				
				if(genotype.code[i] == 1 && nfwd[i] > 0) return(FALSE)
				
				if(genotype.code[i] == 2 && nrev[i] > 0) return(FALSE)
				
				return(NA)
			})
}

.ratio <- function(loglike, g) {
	
	if(g == 0) return(NA)
	
	if(g == 1) {
		return(2*(max(loglike[-c(1,5)]) - max(loglike[c(1,5)])))
	}
	
	if(g == 2) {
		return(2*(max(loglike[-c(2,6)]) - max(loglike[c(2,6)])))
	}
	
	
	if(g == 3) {
		return(2*(max(loglike[-3]) - max(loglike[3])))
	}

	if(g == 4) {
		return(2*(max(loglike[-4]) - max(loglike[4])))
	}
	
	# We should nerver get to this point
	stop("Internal errro")
}

#' position.info
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @method position.info epiG
#' @export
position.info.epiG <- function(object, pos, ...) {
	
	if(length(pos) == 0) {
		return(NULL)
	}
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(length(pos) > 1) {
			
			info.df <- NULL
			
			for(p in pos) {
				info.df <- rbind(info.df, position.info(object, p))
			}
			
			return(info.df)
		}
		
		if(coverage(object, pos) == 0) {
			#Return data.frame
			return(data.frame(
							position = pos, 
							chain.id = NA, 
							ref = NA, 
							alt = NA, 
							genotype = NA, 
							fit.ratio = NA,
							ref.ratio = NA,
							alt.ratio = NA,
							methylated = NA, 
							CpG = NA,
							nreads = NA, 
							nreads.fwd = NA, 
							nreads.rev = NA))
		}	
		
		chains <- object$haplotype$chain[object$read_ids[[pos - start(object)+1]]]
		strands <- object$strands[object$read_ids[[pos - start(object)+1]]]
		
		cid <- sort(unique(chains))
		nfwd <- sapply(cid, function(i) sum(strands[chains == i] == "fwd"))
		nrev <- sapply(cid, function(i) sum(strands[chains == i] == "rev"))
		g <- genotype(object, pos, remove.meth = TRUE)
		ll <- loglike(object, pos)
		
		info.df <- data.frame(
				position = pos, 
				chain.id = cid, 
				ref = NA, 
				alt = NA, 
				genotype = symbols(g)[as.character(cid)], 
				fit.ratio = sapply(as.character(cid), function(x) .ratio(ll[,x], g[x])),
				ref.ratio = NA,
				alt.ratio = NA,
				methylated =.methylation.status(genotype(object, pos, remove.meth = FALSE)[as.character(cid)], nfwd, nrev),  
				CpG = NA,
				nreads = sapply(cid, function(x) sum(chains == x)),
				nreads.fwd = nfwd,
				nreads.rev = nrev
		)
				
		if(!is.null(object[["ref"]])) {
			ref0 <- object$ref[pos - object$offset + 1]
			ref1 <- object$ref[pos - object$offset + 2]
			info.df$ref <- symbols(ref0)
			info.df$CpG <- ref0 == 1 & ref1 == 2
			info.df$ref.ratio = sapply(as.character(cid), function(x) .ratio(ll[,x], ref0))
			
		} 
		
		if(!is.null(object[["alt"]])) {
			alt <- object$alt[pos - object$offset + 1]
			info.df$alt <- symbols(alt)			
			info.df$alt.ratio = sapply(as.character(cid), function(x) .ratio(ll[,x], alt))
		} 
		
		return(info.df)
		
	}
	
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		
		tmp <- lapply(object, function(x) position.info(x, pos[pos %in% start(x):end(x)]), ...)		
		
		# Adjust chain.id
		chain.id.offset = 0
		for(i in 1:length(tmp)) {
			tmp[[i]]$chain.id <- tmp[[i]]$chain.id + chain.id.offset
			chain.id.offset <- chain.id.offset + max(object[[i]]$haplotype$chain, na.rm = TRUE)
		}
		
		return(do.call("rbind", tmp))
	}
	
	stop("Unknown class")
	
}

#' chain.info
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
chain.info <- function(object, ... ) UseMethod("chain.info")

#' chain.info
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @method chain.info epiG
#' @export
chain.info.epiG <- function(object, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		
		if(nchain(object) == 0) {
			stop("Contains no chains")
		}
		
		chains <- sort(unique(object$haplotype$chain))
		
		total.bp = sapply(chains, function(i) sum(object$reads$length[object$haplotype$chain == i]))
		depth <- sapply(1:length(chains), function(i) 
						sum(object$reads$length[(object$reads$positions + object$offset) %in% object$haplotype$start[i]:object$haplotype$end[i]])
						)
						
		return(data.frame(chain.id = chains, 
						start = object$haplotype$start, 
						end = object$haplotype$end, 
						length = object$haplotype$end - object$haplotype$start + 1,
						nreads = as.vector(table(object$haplotype$chain)), 
						nreads.fwd = sapply(chains, function(i) sum(object$strands[object$haplotype$chain == i] == "fwd")),
						nreads.rev = sapply(chains, function(i) sum(object$strands[object$haplotype$chain == i] == "rev")),
						depth.fraction = total.bp/depth
				))
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		tmp <- lapply(object, function(x) chain.info(x, ...))		
			
		# Adjust chain.id
		chain.id.offset = 0
		for(i in 1:length(tmp)) {
			tmp[[i]]$chain.id <- tmp[[i]]$chain.id + chain.id.offset
			chain.id.offset <- max(tmp[[i]]$chain.id, na.rm = TRUE)
		}
		
		return(do.call("rbind", tmp))
	}
	
	stop("Unknown class")
	
}

#' read.info
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
read.info <- function(object, ... ) UseMethod("read.info")

#' read.info
#' @param object 
#' @param  inc.symbols 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @method read.info epiG
#' @export
read.info.epiG <- function(object, inc.symbols = FALSE, ...) {
		
	if(paste(class(object), collapse = ".") == "epiG") {

		if(!("reads" %in% names(object))) {
			object <- fetch.reads(object)
		}
		
		info <- NULL
		
		if(inc.symbols) {
			for(idx in 1:nread(object)) {
				tmp <- data.frame(
						name = object$reads$name[idx], 
						position = object$reads$position[idx]:(object$reads$position[idx]+object$reads$length[idx]-1) + object$offset, 
						symbol = symbols(object$reads$reads[[idx]]), 
						read.id = idx, quality = object$reads$quality[[idx]],
						chain.id = object$haplotype$chain[idx], 
						strand = object$strands[idx])
			
				info <- rbind(info, tmp)	
			}
		} else {
			info <- data.frame(
						name = object$reads$name, 
						start = object$reads$position + object$offset,  
						end = object$reads$position+object$reads$length-1 + object$offset,
						length = object$reads$length,
						read.id = 1:nread(object),
						chain.id = object$haplotype$chain, 
						strand = object$strands
					)
		}
	
		return(info)
	}
	
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		tmp <- lapply(object, function(x) read.info(x, inc.symbols, ...))
		
		# Adjust chain.id read.id
		chain.id.offset = 0
		read.id.offset <- 0
		for(i in 1:length(tmp)) {
			
			tmp[[i]]$chain.id <- tmp[[i]]$chain.id + chain.id.offset
			chain.id.offset <- max(tmp[[i]]$chain.id, na.rm = TRUE)
			
			tmp[[i]]$read.id <- tmp[[i]]$read.id + read.id.offset
			read.id.offset <- max(tmp[[i]]$read.id)
			
			tmp[[i]]$chunk <- i
		}
		
		return(do.call("rbind", tmp))
	}
	
	stop("Unknown class")
	
}

read.info.epiG_reads <- function(object, inc.symbols = FALSE, ...) {
	
	info <- NULL
	
	if(inc.symbols) {
		for(idx in 1:length(object$reads)) {
			tmp <- data.frame(
					name = object$names[idx], 
					position = object$position[idx]:(object$position[idx]+object$length[idx]-1), 
					symbol = symbols(object$reads[[idx]]),
					quality = object$quality[[idx]],
					read.id = idx
					)
			
			info <- rbind(info, tmp)	
		}
	} else {
		info <- data.frame(
				name = object$names, 
				start = object$position,  
				end = object$position+object$length-1,
				length = object$length,
				read.id = 1:length(object$reads)
		)
	}
	
	return(info)
}


#' Number of chunks 
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
nchunks <- function(object, pos, ... ) UseMethod("nchunks")

#' Number of chunks 
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @method nchunks epiG
#' @export
nchunks.epiG <- function(object, ...) {
	if(paste(class(object), collapse = ".") == "epiG") {
		return(1)
	}
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sum(sapply(object, nchunks)))
	}
		
}

#' Number of chains
#' @param object 
#' @param pos 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
nchain <- function(object, pos, ... ) UseMethod("nchain")

#' Number of chains
#' @param object 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @method nchain epiG
#' @export
nchain.epiG <- function(object, ...) {
	
	if(paste(class(object), collapse = ".") == "epiG") {
		return(length(object$haplotype$start))
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		return(sapply(object, nchain))
	}
	
}

#' subregion
#' @param object 
#' @param start 
#' @param end 
#' @param chop.reads 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @export
subregion <- function(object, start, end, chop.reads = FALSE, ... ) UseMethod("subregion")

#' subregion
#' @param object 
#' @param start 
#' @param end 
#' @param chop.reads 
#' @param ... 
#' @return ??
#' 
#' @author Martin Vincent
#' @method subregion epiG
#' @export
subregion.epiG <- function(object, start, end, chop.reads = FALSE, ...) {
		
	if(paste(class(object), collapse = ".") == "epiG") {

		if(end < start(object) || start > end(object)) {
			stop("Out of range")
		}
		
		if(start < start(object) || end > end(object)) {
			stop("Out of range")
		}
				
		rel_start_pos <- max(start, start(object)) - object$offset
		rel_end_pos <- min(end, end(object)) - object$offset
		
		new_object <- NULL
		
		new_object$epiG_version <- object$epiG_version
		new_object$date <- object$date
		new_object$config <- object$config
		new_object$filename <- object$filename
		new_object$refname <- object$refname
		
		new_object$offset <- as.integer(start)
		new_object$length <- as.integer(rel_end_pos - rel_start_pos + 1)		
		
		new_object$read_ids <- object$read_ids[rel_start_pos:rel_end_pos+1]
		remaining_read_ids <- sort(unique(unlist(new_object$read_ids)))
	
		new_object$strands <- object$strands[remaining_read_ids]
		
		#FIXME add unique read ids + read names
		
		#Recalibrate read_ids
		new_id <- 1:length(remaining_read_ids)
		names(new_id) <- remaining_read_ids
		new_object$read_ids <- lapply(new_object$read_ids, function(x) as.integer(new_id[as.character(x)]))
		
		if(length(remaining_read_ids)) {
			warning("Region is empty")
		}
		
		# Ref and alt
		
		if("ref" %in% names(object)) {
			new_object$ref <-  object$ref[rel_start_pos:rel_end_pos+1]
		}
		
		if("alt" %in% names(object)) {
			new_object$alt <-  object$alt[rel_start_pos:rel_end_pos+1]
		}
		
		
		# Haplotype and genotype
		new_object$haplotype$chain <- object$haplotype$chain[remaining_read_ids]
		remaining_chains <- sort(unique(new_object$haplotype$chain))
		new_object$haplotype$chain <- as.integer(factor(new_object$haplotype$chain))
							
		new_object$haplotype$start <- sapply(object$haplotype$start[remaining_chains], function(x) max(start, x))
		new_object$haplotype$end <- sapply(object$haplotype$end[remaining_chains], function(x) min(end, x))
		
		new_object$genotypes <- list()
		new_object$loglikes <- list()
		for(i in remaining_chains) {
				chain_rel_start_pos <- max(start, object$haplotype$start[i]) - object$haplotype$start[i]
				chain_rel_end_pos <- min(end, object$haplotype$end[i]) - object$haplotype$start[i]
				
				new_object$genotypes[[length(new_object$genotypes)+1]] <- object$genotypes[[i]][chain_rel_start_pos:chain_rel_end_pos+1]		
				new_object$loglikes[[length(new_object$loglikes)+1]] <- object$loglikes[[i]][chain_rel_start_pos:chain_rel_end_pos+1,, drop = FALSE]				
		}
	
		# Reads
		if("reads" %in% names(object)) {
						
			new_object$reads$reads <- object$reads$reads[remaining_read_ids]
			new_object$reads$quality <- object$reads$quality[remaining_read_ids]
			new_object$reads$positions <- object$reads$positions[remaining_read_ids] - rel_start_pos
			new_object$reads$lengths <- object$reads$lengths[remaining_read_ids]			
			new_object$reads$names <- object$reads$names[remaining_read_ids]			
			
			if(chop.reads) {
			
				#FIXME bug in chop.reads
				
				for(i in 1:length(new_object$reads$lengths)) {
					
					read_rel_start_pos <- max(0,-new_object$reads$positions[i]) 
					read_length <- new_object$reads$lengths[i] - read_rel_start_pos - max(0, new_object$reads$positions[i] + new_object$reads$lengths[i] - new_object$length)
												
					new_object$reads$positions[i] <- new_object$reads$positions[i] + read_rel_start_pos
					new_object$reads$lengths[i] <- read_length
					new_object$reads$reads[[i]] <- new_object$reads$reads[[i]][(read_rel_start_pos+1):(read_rel_start_pos+read_length)]
					new_object$reads$quality[[i]] <- new_object$reads$quality[[i]][(read_rel_start_pos+1):(read_rel_start_pos+read_length)]
					
				}
			}
		}
		
		class(new_object) <- "epiG"
		return(new_object)
	}
	
	if(paste(class(object), collapse = ".") == "epiG.chunks") {
		
		new_object <- list()
		j <- 1
		
		for(i in 1:nchunks(object)) {
			if((start(object[[i]]) <= end && end <= end(object[[i]])) || (start(object[[i]]) <= start && start <= end(object[[i]]))) {
				new_object[[j]] <- subregion(object[[i]], max(start, start(object[[i]])), min(end, end(object[[i]])), chop.reads, ...)
				j <- j + 1
			} 		
		}
		
		if(j == 1) {
			stop("Out of range")
		}
		
		if(j == 2) {
			new_object <- new_object[[1]]
			class(new_object) <- c("epiG")
			return(new_object)
		}
		
		class(new_object) <- c("epiG", "chunks")
		return(new_object)
	}
}
	
#' symbols
#' 
#' @param g 
#' @return symbols
#' 
#' @author martin
#' @export
symbols <- function(g) {
		return(unlist(sapply(g, function(x) c("N", "C", "G", "A", "T", "c", "g")[x+1])))
}
	
	
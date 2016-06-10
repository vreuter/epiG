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


#' @title Print Information About an Fitted epiG Model
#' @description
#' Print Information About an Fitted epiG Model
#' 
#' @param x epiG model
#' @param ... ignored
#' 
#' @author Martin Vincent
#' @method print epiG
#' @export
#'
#' @examples 
#' #TODO examples print
print.epiG <- function(x, ...) {
	
        if("chunks" %in% class(x)) {
            message(paste("Epigenotype model (estimated by epiG version ", x[[1]]$epiG_version, ")", sep=""))
            message()
			#TODO list chunks
           # message("  ",paste(length(x), " chunks covering base position [", min(start(x)), ", ", max(end(x)), "]", " on ", x[[1]]$refname, ", ", sum(sapply(x, function(z) z$length)), " bases in total.", sep=""))
            message()
            message("  ",paste("Estimation was finalized ",x[[1]]$date, sep=""))
         } else {
            message(paste("Epigenotype model (estimated by epiG version ", x$epiG_version, ")", sep=""))
            message()
            message(paste("  Cover base position [", start(x), ", ", end(x), "]", " on ", x$refname, ", ", length(x), " bases in total.", sep=""))
            message()
            message("  ",paste("Estimation was finalized ",x$date, sep=""))
        }

#TODO print number of reads
}

#' @title Print Information About an epiG Configuration
#' @description
#' Print Information About an epiG Configuration
#' 
#' @param x epiG configuration
#' @param ... ignored
#' 
#' @author Martin Vincent
#' @method print epiG.config
#' @export
#'
#' @examples 
#' #TODO examples print.epiG.config
print.epiG.config <- function(x, ...) {
	
	message("\nepiG Config Object")
	
#	[1]         "ref_offset"           
#	[4] "alt_offset"            "max_iterations"        "fwd_model"            
#	[7] "rev_model"             "split_mode"            "model_name"           
#	[10] "fwd_DGCH_model"        "rev_DGCH_model"        "fwd_HCGD_model"       
#	[13] "rev_HCGD_model"        "fwd_CH_model"          "rev_CH_model"         
#	[16] "fwd_C_G_model"         "rev_C_G_model"         "chunk_size"           
#	[19] "chunk_method"          "reads_hard_limit"      "quality_threshold"    
#	[22] "ref_prior"                     "margin"               
#	[28] "max_stages"            "structual_prior_scale" "use_paired_reads"     
#	[31] "verbose"                             
	
	message("\n Files:")
	par <- c("filename", "ref_filename", "alt_filename") 
	tmp <- data.frame(unlist(x[par]))
	colnames(tmp) <- NULL
	print(tmp)
	
	message("\n Site:")
	par <- c("refname", "start" , "end") 
	tmp <- data.frame(unlist(x[par]))
	colnames(tmp) <- NULL
	print(tmp)
	
	message("\n Feasibility settings:")
	par <- c("min_overlap_length", "min_CG_count", "min_HCGD_count", "min_DGCH_count") 
	tmp <- data.frame(unlist(x[par]))
	colnames(tmp) <- NULL
	print(tmp)
}

#' @title Print Information About an epiG Conversion Model
#' @description
#' Print Information About an epiG Conversion Model
#' 
#' @param x epiG conversion model
#' @param ... ignored
#' 
#' @author Martin Vincent
#' @method print epiG.model
#' @export
#'
#' @examples 
#' #TODO examples print.epiG.model 
print.epiG.model <- function(x, ...) {
	
	cat(x$name)
	cat("\n\n")
	cat("Forward model:\n")
	print(x$fwd[[1]])
	cat("\n")
	cat("Reverse model:\n")
	print(x$rev[[1]])
}

#' @title Print Information About an epiG Reads Object
#' @description
#' Print information about an epiG reads object
#' 
#' @param x epiG reads object
#' @param ... ignored
#' 
#' @author Martin Vincent
#' @method print epiG.reads
#' @export
#'
#' @examples 
#' #TODO examples print.epiG.reads 
print.epiG.reads <- function(x, ...) {
	message("epiG reads object:")
}

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


#' @title print
#' @description
#' print
#' 
#' @param x 
#' @param ... 
#' 
#' @author Martin Vincent
#' @method print epiG
#' @export
#'
#' @examples 
#' #TODO
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

#' @title print config
#' @description
#' print config
#' 
#' @param x 
#' @param ... 
#' 
#' @author Martin Vincent
#' @method print epiG.config
#' @export
#'
#' @examples 
#' #TODO
print.epiG.config <- function(x, ...) {
	
	message("epiG config object:")
}

#' @title print model
#' @description
#' print model
#' 
#' @param x 
#' @param ... 
#' 
#' @author Martin Vincent
#' @method print epiG.model
#' @export
#'
#' @examples 
#' #TODO
print.epiG.model <- function(x, ...) {
	
	cat(x$name)
	cat("\n\n")
	cat("Forward model:\n")
	print(x$fwd[[1]])
	cat("\n")
	cat("Reverse model:\n")
	print(x$rev[[1]])
}

print.epiG.reads <- function(x, ...) {
	message("epiG reads object:")
}

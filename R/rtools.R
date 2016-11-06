#' @title Convert to sparse matrix
#'
#' @description
#' Convert sparse matrix retunred from .Call to sparseMatrix.
#'
#' @param x .Call retunred list
#' @author Martin Vincent
#' @export
sparseMatrix_from_C_format <- function(x) {
	sparseMatrix(p = x[[2]], i = x[[3]], x = x[[4]], dims = x[[1]], index1 = FALSE)
}

#' @title Prepare sparse matrix for .Call
#'
#' @description
#' Prepare sparse matrix for .Call
#'
#' @param x a spares matrix
#' @author Martin Vincent
#' @export
sparseMatrix_to_C_format <- function(x) {
	x <- as(x, "CsparseMatrix")
	return( list(dim(x), x@p, x@i, x@x) )
}

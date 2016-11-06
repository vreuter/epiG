#' @title Test internal rtools
#' @description
#' This function runs some internal tests and is not intended for users of the package.
#'
#' @author Martin Vincent
#' @export
#' @useDynLib epiG double_rtools_test
#' @useDynLib epiG field_double_rtools_test
#' @useDynLib epiG u32_rtools_test
#' @useDynLib epiG field_u32_rtools_test
#' @useDynLib epiG int_rtools_test
#' @useDynLib epiG field_int_rtools_test
#' @useDynLib epiG bool_rtools_test
#' @useDynLib epiG field_bool_rtools_test
#' @useDynLib epiG string_rtools_test
#' @useDynLib epiG field_string_rtools_test
#' @useDynLib epiG mat_double_rtools_test
#' @useDynLib epiG field_mat_double_rtools_test
#' @useDynLib epiG col_double_rtools_test
#' @useDynLib epiG field_col_double_rtools_test
#' @useDynLib epiG col_u32_rtools_test
#' @useDynLib epiG field_col_u32_rtools_test
#' @useDynLib epiG col_s32_rtools_test
#' @useDynLib epiG field_col_s32_rtools_test
#' @useDynLib epiG field_sp_mat_rtools_test
#' @useDynLib epiG sp_mat_rtools_test
#' @importFrom utils packageVersion
#' @importFrom stats runif

#' @import Matrix
test_rtools <- function() {

  errors <- NULL

  test_values <-  runif(100, -100, 100)
  if( ! all(sapply(test_values, function(x) x == .Call("double_rtools_test", x)))) {
    errors <- c(errors, "double test fail")
  }
  test_values <- as.list(test_values)
  if( ! all.equal(test_values, .Call("field_double_rtools_test", test_values))) {
    errors <- c(errors, "double test field fail")
  }

  test_values <-  runif(0, 1e5, 100)
  if( ! all(sapply(test_values, function(x) x == .Call("u32_rtools_test", x)))) {
    errors <- c(errors, "u32 test fail")
  }
  test_values <- as.list(test_values)
  if( ! all.equal(test_values, .Call("field_u32_rtools_test", test_values))) {
    errors <- c(errors, "u32 test field fail")
  }

  test_values <-  sample(-1e5:1e5, size = 100)
  if( ! all(sapply(test_values, function(x) x == .Call("int_rtools_test", x)))) {
    errors <- c(errors, "int test fail")
  }
  test_values <- as.list(test_values)
  if( ! all.equal(test_values, .Call("field_int_rtools_test", test_values))) {
    errors <- c(errors, "int test field fail")
  }

  test_values <-  c(TRUE, FALSE)
  if( ! all(sapply(test_values, function(x) x == .Call("bool_rtools_test", x)))) {
    errors <- c(errors, "bool test fail")
  }
  test_values <- as.list(test_values)
  if( ! all.equal(test_values, .Call("field_bool_rtools_test", test_values))) {
    errors <- c(errors, "bool test field fail")
  }

  test_values <- sapply(1:100, function(i)
    paste(sample(c(LETTERS, letters, ".", "\\", ","), size = i, replace = TRUE), collapse = ""))
  if( ! all(sapply(test_values, function(x) x == .Call("string_rtools_test", x)))) {
      errors <- c(errors, "string test fail")
  }
  test_values <- as.list(test_values)
  if( ! all.equal(test_values, .Call("field_string_rtools_test", test_values))) {
    errors <- c(errors, "string test field fail")
  }

  test_values <- replicate(25, matrix(runif(150, -1e5, 1e5), nrow = 10, ncol = 15), simplify = FALSE)
  if( ! all(sapply(test_values, function(x) all(x == .Call("mat_double_rtools_test", x))))) {
      errors <- c(errors, "mat_double test fail")
  }
  if( ! all.equal(test_values, .Call("field_mat_double_rtools_test", test_values))) {
    errors <- c(errors, "mat_double test field fail")
  }

  test_values <- replicate(25, runif(150, -1e5, 1e5), simplify = FALSE)
  if( ! all(sapply(test_values, function(x) all(x == .Call("col_double_rtools_test", x))))) {
      errors <- c(errors, "col_double test fail")
  }
  if( ! all.equal(test_values, .Call("field_col_double_rtools_test", test_values))) {
    errors <- c(errors, "col_double test field fail")
  }

  test_values <- replicate(25, sample(0:1e5, size = 500), simplify = FALSE)
  if( ! all(sapply(test_values, function(x) all(x == .Call("col_u32_rtools_test", x))))) {
      errors <- c(errors, "col_u32 test fail")
  }
  if( ! all.equal(test_values, .Call("field_col_u32_rtools_test", test_values))) {
    errors <- c(errors, "col_u32 test field fail")
  }

  test_values <- replicate(25, sample(-1e5:1e5, size = 500), simplify = FALSE)
  if( ! all(sapply(test_values, function(x) all(x == .Call("col_s32_rtools_test", x))))) {
      errors <- c(errors, "col_s32 test fail")
  }
  if( ! all.equal(test_values, .Call("field_col_s32_rtools_test", test_values))) {
    errors <- c(errors, "col_s32 test field fail")
  }

  test_values <- replicate(25,
    sparseMatrix(
      i = sample(1:100, size = 50),
      j = sample(1:100, size = 50),
      x = runif(50, -1e5, 1e5)
    ),
    simplify = FALSE
  )
  test_values <- lapply(test_values, sparseMatrix_to_C_format)
  if( ! all(sapply(test_values, function(x) all.equal(x, .Call("sp_mat_rtools_test", x))))) {
    errors <- c(errors, "sp_mat test fail")
  }
  if( ! all.equal(test_values, .Call("field_sp_mat_rtools_test", test_values))) {
    errors <- c(errors, "sp_mat test field fail")
  }

  if(length(errors) > 0) {
    stop(paste("\n", errors, collapse = ""))
  }
}

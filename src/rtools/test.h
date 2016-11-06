
extern "C" {
  SEXP double_rtools_test(SEXP  exp);
  SEXP int_rtools_test(SEXP  exp);
  SEXP char_rtools_test(SEXP  exp);
  SEXP u32_rtools_test(SEXP  exp);
  SEXP bool_rtools_test(SEXP  exp);
  SEXP string_rtools_test(SEXP  exp);
  SEXP mat_double_rtools_test(SEXP  exp);
  SEXP col_double_rtools_test(SEXP  exp);
  SEXP col_u32_rtools_test(SEXP  exp);
  SEXP col_s32_rtools_test(SEXP  exp);
  SEXP sp_mat_rtools_test(SEXP  exp);

  SEXP field_double_rtools_test(SEXP  exp);
  SEXP field_int_rtools_test(SEXP  exp);
  SEXP field_char_rtools_test(SEXP  exp);
  SEXP field_u32_rtools_test(SEXP  exp);
  SEXP field_bool_rtools_test(SEXP  exp);
  SEXP field_string_rtools_test(SEXP  exp);
  SEXP field_mat_double_rtools_test(SEXP  exp);
  SEXP field_col_double_rtools_test(SEXP  exp);
  SEXP field_col_u32_rtools_test(SEXP  exp);
  SEXP field_col_s32_rtools_test(SEXP  exp);
  SEXP field_sp_mat_rtools_test(SEXP  exp);
}

template<typename T>
SEXP rtools_test(SEXP  exp) {
  const T x = get_value<T>(exp);
  return rObject(x);
}

SEXP double_rtools_test(SEXP  exp) {
  return rtools_test<double>(exp);
}

SEXP int_rtools_test(SEXP  exp) {
  return rtools_test<int>(exp);
}

SEXP u32_rtools_test(SEXP  exp) {
  return rtools_test<arma::u32>(exp);
}

SEXP bool_rtools_test(SEXP  exp) {
  return rtools_test<bool>(exp);
}

SEXP string_rtools_test(SEXP  exp) {
  return rtools_test<std::string>(exp);
}

SEXP mat_double_rtools_test(SEXP  exp) {
  return rtools_test<arma::Mat<double> >(exp);
}

SEXP col_double_rtools_test(SEXP  exp) {
  return rtools_test<arma::Col<double> >(exp);
}

SEXP col_u32_rtools_test(SEXP  exp) {
  return rtools_test<arma::Col<arma::u32> >(exp);
}

SEXP col_s32_rtools_test(SEXP  exp) {
  return rtools_test<arma::Col<arma::s32> >(exp);
}

SEXP sp_mat_rtools_test(SEXP  exp) {
  return rtools_test<arma::sp_mat>(exp);
}

// TEST fields
template<typename T>
SEXP rtools_test_field(SEXP  exp) {
  const arma::field<T> x = get_field<T>(exp);
  return rObject(x);
}

SEXP field_double_rtools_test(SEXP  exp) {
  return rtools_test_field<double>(exp);
}

SEXP field_int_rtools_test(SEXP  exp) {
  return rtools_test_field<int>(exp);
}

SEXP field_u32_rtools_test(SEXP  exp) {
  return rtools_test_field<arma::u32>(exp);
}

SEXP field_bool_rtools_test(SEXP  exp) {
  return rtools_test_field<bool>(exp);
}

SEXP field_string_rtools_test(SEXP  exp) {
  return rtools_test_field<std::string>(exp);
}

SEXP field_mat_double_rtools_test(SEXP  exp) {
  return rtools_test_field<arma::Mat<double> >(exp);
}

SEXP field_col_double_rtools_test(SEXP  exp) {
  return rtools_test_field<arma::Col<double> >(exp);
}

SEXP field_col_u32_rtools_test(SEXP  exp) {
  return rtools_test_field<arma::Col<arma::u32> >(exp);
}

SEXP field_col_s32_rtools_test(SEXP  exp) {
  return rtools_test_field<arma::Col<arma::s32> >(exp);
}

SEXP field_sp_mat_rtools_test(SEXP  exp) {
  return rtools_test_field<arma::sp_mat>(exp);
}

//TODO rlist

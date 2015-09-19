/*
 Lightweight tools for R and c++ integration.
 Copyright (C) 2012 Martin Vincent

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

#ifndef ROBJECT_DEF_H_
#define ROBJECT_DEF_H_

//Copy constructor
rObject::rObject(rObject const& s) :
		exp(s.exp), number_of_protects(s.number_of_protects), unprotect_on_destruction(
				s.unprotect_on_destruction), exp_counter(s.exp_counter) {

	(*exp_counter)++;

}
rObject rObject::operator=(rObject const& s) {

	this->exp = s.getSEXP();
	const_cast<int&>(this->number_of_protects) = s.n_protects();
	const_cast<bool*&>(this->unprotect_on_destruction) =
			s.unprotect_on_destruction;
	const_cast<int*&>(this->exp_counter) = s.exp_counter;

	(*exp_counter)++;

	return *this;
}

rObject::rObject(SEXP exp) :
		exp(exp), number_of_protects(0), unprotect_on_destruction(new bool), exp_counter(
				new int) {

	*unprotect_on_destruction = false;
	*exp_counter = 1;
}

rObject::rObject(arma::u32 value) :
		number_of_protects(1), unprotect_on_destruction(new bool), exp_counter(
				new int) {

	*this->unprotect_on_destruction = true;
	*exp_counter = 1;

	PROTECT(exp = Rf_allocVector(INTSXP, 1));
	INTEGER(exp)[0] = value;
}

rObject::rObject(double value) :
		number_of_protects(1), unprotect_on_destruction(new bool), exp_counter(
				new int) {

	*this->unprotect_on_destruction = true;
	*exp_counter = 1;

	PROTECT(exp = Rf_allocVector(REALSXP, 1));
	REAL(exp)[0] = value;
}

rObject::rObject(int value) :
        number_of_protects(1), unprotect_on_destruction(new bool), exp_counter(
                new int) {

    *this->unprotect_on_destruction = true;
    *exp_counter = 1;

    PROTECT(exp = Rf_allocVector(INTSXP, 1));
    INTEGER(exp)[0] = value;
}

rObject::rObject(arma::Mat<double> const& m) :
		number_of_protects(2), unprotect_on_destruction(new bool), exp_counter(
				new int) {

	*this->unprotect_on_destruction = true;
	*exp_counter = 1;

	SEXP matrixDim;
	PROTECT(matrixDim = Rf_allocVector(INTSXP, 2));
	INTEGER(matrixDim)[0] = m.n_rows;
	INTEGER(matrixDim)[1] = m.n_cols;

	PROTECT(exp = Rf_allocVector(REALSXP, m.n_elem));

	//Copy data
    arma::arrayops::copy(REAL(exp), m.mem, m.n_elem);

	Rf_setAttrib(exp, R_DimSymbol, matrixDim);
}

rObject::rObject(arma::Mat<arma::u32> const& m) :
		number_of_protects(2), unprotect_on_destruction(new bool), exp_counter(
				new int) {

	*this->unprotect_on_destruction = true;
	*exp_counter = 1;

	SEXP matrixDim;
	PROTECT(matrixDim = Rf_allocVector(INTSXP, 2));
	INTEGER(matrixDim)[0] = m.n_rows;
	INTEGER(matrixDim)[1] = m.n_cols;

	PROTECT(exp = Rf_allocVector(INTSXP, m.n_rows * m.n_cols));

	//Copy data
	copy_cast(INTEGER(exp), m.mem, m.n_elem);

	Rf_setAttrib(exp, R_DimSymbol, matrixDim);
}

rObject::rObject(arma::Mat<int> const& m) :
        number_of_protects(2), unprotect_on_destruction(new bool), exp_counter(
                new int) {

    *this->unprotect_on_destruction = true;
    *exp_counter = 1;

    SEXP matrixDim;
    PROTECT(matrixDim = Rf_allocVector(INTSXP, 2));
    INTEGER(matrixDim)[0] = m.n_rows;
    INTEGER(matrixDim)[1] = m.n_cols;

    PROTECT(exp = Rf_allocVector(INTSXP, m.n_rows * m.n_cols));

    //Copy data
    copy_cast(INTEGER(exp), m.mem, m.n_elem);

    Rf_setAttrib(exp, R_DimSymbol, matrixDim);
}

rObject::rObject(arma::Mat<char> const& m) :
        number_of_protects(2), unprotect_on_destruction(new bool), exp_counter(
                new int) {

    *this->unprotect_on_destruction = true;
    *exp_counter = 1;

    SEXP matrixDim;
    PROTECT(matrixDim = Rf_allocVector(INTSXP, 2));
    INTEGER(matrixDim)[0] = m.n_rows;
    INTEGER(matrixDim)[1] = m.n_cols;

    PROTECT(exp = Rf_allocVector(INTSXP, m.n_rows * m.n_cols));

    //Copy data
    copy_cast(INTEGER(exp), m.mem, m.n_elem);

    Rf_setAttrib(exp, R_DimSymbol, matrixDim);
}

rObject::rObject(arma::Col<double> const& v) :
		number_of_protects(1), unprotect_on_destruction(new bool), exp_counter(
				new int) {

	*this->unprotect_on_destruction = true;
	*exp_counter = 1;

	PROTECT(exp = Rf_allocVector(REALSXP, v.n_elem));

	//Copy data
    arma::arrayops::copy(REAL(exp), v.mem, v.n_elem);

}

rObject::rObject(arma::Col<arma::u32> const& v) :
		number_of_protects(1), unprotect_on_destruction(new bool), exp_counter(
				new int) {

	*this->unprotect_on_destruction = true;
	*exp_counter = 1;

	PROTECT(exp = Rf_allocVector(INTSXP, v.n_elem));

	//Copy data
	copy_cast(INTEGER(exp), v.mem, v.n_elem);
}

rObject::rObject(arma::Col<int> const& v) :
        number_of_protects(1), unprotect_on_destruction(new bool), exp_counter(
                new int) {

    *this->unprotect_on_destruction = true;
    *exp_counter = 1;

    PROTECT(exp = Rf_allocVector(INTSXP, v.n_elem));

    //Copy data
    arma::arrayops::copy(INTEGER(exp), v.mem, v.n_elem);
}

rObject::rObject(arma::Col<char> const& v) :
        number_of_protects(1), unprotect_on_destruction(new bool), exp_counter(
                new int) {

    *this->unprotect_on_destruction = true;
    *exp_counter = 1;

    PROTECT(exp = Rf_allocVector(INTSXP, v.n_elem));

    //Copy data
    copy_cast(INTEGER(exp), v.mem, v.n_elem);
}

rObject::rObject(std::vector<int> const& v) :
        number_of_protects(1), unprotect_on_destruction(new bool), exp_counter(
                new int) {

    *this->unprotect_on_destruction = true;
    *exp_counter = 1;

    PROTECT(exp = Rf_allocVector(INTSXP, v.size()));

    //Copy data
    arma::arrayops::copy(INTEGER(exp), v.data(), v.size());
}

rObject::rObject(std::vector<std::string> const& v) :
		 number_of_protects(1),
		 unprotect_on_destruction(new bool),
		 exp_counter(new int) {

	 	*this->unprotect_on_destruction = true;
	    *exp_counter = 1;

	    PROTECT(exp = Rf_allocVector(STRSXP, v.size()));

	    //Copy data
	    for(int i = 0; i < v.size(); i++)  {
	    	SET_STRING_ELT(exp, i, mkChar(v[i].c_str()));
	    }
}


rObject::rObject(std::string const& v) :
		 number_of_protects(1),
		 unprotect_on_destruction(new bool),
		 exp_counter(new int) {

	 	*this->unprotect_on_destruction = true;
	    *exp_counter = 1;

	    PROTECT(exp = Rf_allocVector(STRSXP, 1));

	    SET_STRING_ELT(exp, 0, mkChar(v.c_str()));
}

rObject::rObject(arma::sp_mat const& m) :
		number_of_protects(5), unprotect_on_destruction(new bool), exp_counter(
				new int) {

	*this->unprotect_on_destruction = true;
	*exp_counter = 1;

	PROTECT(exp = Rf_allocVector(VECSXP, 4)); // Creating a list with 4 elements

	//TODO names on list

	SEXP dim;
	PROTECT(dim = Rf_allocVector(INTSXP, 2));
	SET_VECTOR_ELT(exp, 0, dim);
	INTEGER(dim)[0] = m.n_rows;
	INTEGER(dim)[1] = m.n_cols;

	SEXP col_ptrs;
	PROTECT(col_ptrs = Rf_allocVector(INTSXP, m.n_cols + 1));
	SET_VECTOR_ELT(exp, 1, col_ptrs);
	copy_cast(INTEGER(col_ptrs), m.col_ptrs, m.n_cols + 1);

	SEXP row_indices;
	PROTECT(row_indices = Rf_allocVector(INTSXP, m.n_nonzero));
	SET_VECTOR_ELT(exp, 2, row_indices);
	copy_cast(INTEGER(row_indices), m.row_indices, m.n_nonzero);

	SEXP values;
	PROTECT(values = Rf_allocVector(REALSXP, m.n_nonzero));
	SET_VECTOR_ELT(exp, 3, values);
    arma::arrayops::copy(REAL(values), m.values, m.n_nonzero);
}

template<typename T>
rObject::rObject(arma::field<T> const& field) :
        number_of_protects(1), unprotect_on_destruction(new bool), exp_counter(
                new int) {

    *this->unprotect_on_destruction = unprotect_on_destruction;
    *exp_counter = 1;

    PROTECT(exp = Rf_allocVector(VECSXP, field.n_elem)); // Creating a list with n_elem elements

    //Construct list
    unsigned int i;
    for (i = 0; i < field.n_elem; i++) {
        // attaching
        rObject tmp(field(i));
        SET_VECTOR_ELT(exp, i, tmp);
    }

}

rObject::rObject(rList const& list) :
		number_of_protects(2), unprotect_on_destruction(new bool), exp_counter(
				new int) {

	*this->unprotect_on_destruction = true;
	*exp_counter = 1;

	int number_of_elements = list.length();
	SEXP names;

	PROTECT(exp = Rf_allocVector(VECSXP, number_of_elements));
	PROTECT(names = Rf_allocVector(VECSXP, number_of_elements));

	//Create R list
	for (int i = 0; i < number_of_elements; ++i) {

		//Take over protection
		number_of_protects += list.get(i).takeover_protection();

		//Add to list
		SET_VECTOR_ELT(exp, i, list.get(i));
		SET_VECTOR_ELT(names, i, Rf_mkChar(list.getName(i).c_str()));
	}

	Rf_setAttrib(exp, R_NamesSymbol, names);

}

rObject::~rObject() {

	if (*exp_counter == 1) {

		//Last exp holder

		delete exp_counter;

		if (*unprotect_on_destruction) {
			UNPROTECT(number_of_protects);
		}
	}

	else {

		(*exp_counter)--;

	}

}

rObject::operator SEXP() const {
	return getSEXP();
}

SEXP rObject::getSEXP() const {
	return exp;
}

int rObject::n_protects() const {
	return number_of_protects;
}

//TODO remove -- we do not need this as no protection is need for elements of lists
int rObject::takeover_protection() {

	*unprotect_on_destruction = false;

	int r = number_of_protects;

	number_of_protects = 0;

	return r;
}

#endif /* ROBJECT_DEF_H_ */

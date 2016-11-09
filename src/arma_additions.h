/*
	Sgl template library for optimizing sparse group lasso penalized objectives.
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

#ifndef ARMA_ADDITIONS_H_
#define ARMA_ADDITIONS_H_

template<typename T>
bool is_zero(arma::Mat<T> obj) {

	typename arma::Mat<T>::const_iterator x = obj.begin();
	for(; x != obj.end(); ++x) {
		if(*x != 0) {
			return false;
		}
	}

	return true;
}

template<class E, class F>
  static void
  copy_cast(E* target, const F* source, arma::uword size)
  {

    for (arma::uword i = 0; i < size; ++i)
      {
        target[i] = static_cast<E>(source[i]);
      }

  }

#endif /* ARMA_ADDITIONS_H_ */

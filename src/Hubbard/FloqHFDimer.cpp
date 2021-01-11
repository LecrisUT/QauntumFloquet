//
// Created by Le Minh Cristian on 2020/12/07.
//

#include "FloqHFDimer.h"

using namespace QuanFloq;

template<typename T>
FloqHFDimer<T>::FloqHFDimer( int nFH, int nF_max, T v0, T v1, T U, T w, T t ) :
		vHamil<T>(2, std::is_same_v<T, float> || std::is_same_v<T, double> ? Sym : Her, false),
		vHFHamil<T>(2, 1, 2, false),
		vFloqHamil<T>(nFH, nF_max, true),
		vFloqHFHamil<T>(1, true) {
	this->indvH[0] = this->H;
	this->indvH[1] = this->H + 3;
	vFloqHamil<T>::setW(w, false);
	vFloqHFDimer<T>::setU(U, false);
	vFloqHFDimer<T>::setT(t, false);
	vFloqHFDimer<T>::setV0(v0, false);
	vFloqHFDimer<T>::setV1(v1, false);
	vFloqHFHamil<T>::Initialize();
	this->CalcSH();
}

// region Initialize templates
#ifdef BUILD_FLOAT
template
class QuanFloq::FloqHFDimer<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::FloqHFDimer<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::FloqHFDimer<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::FloqHFDimer<cdouble >;
#endif
// endregion
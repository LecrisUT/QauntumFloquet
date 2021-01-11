//
// Created by Le Minh Cristian on 2020/12/07.
//

#include "Dimer.h"

using namespace QuanFloq;

template<typename T>
Dimer<T>::Dimer( T v, T U, T t ) :
		vHamil<T>(3, std::is_same_v<T, float> || std::is_same_v<T, double> ? Sym : Her, true) {
//	T t2 = tt * sqrt((T)2.0f);
//	T M[9] = {tU+tv, t2, 0.0f,
//	          t2, 0.0f, t2,
//	          0.0f, t2, tU - tv};
//	this->setH(M);
	vDimer<T>::setT(t);
	vDimer<T>::setU(U);
	vDimer<T>::setV(v);
}

// region Initialize templates
#ifdef BUILD_FLOAT
template
class QuanFloq::Dimer<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::Dimer<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::Dimer<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::Dimer<cdouble >;
#endif
// endregion
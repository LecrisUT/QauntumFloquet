//
// Created by Le Minh Cristian on 2020/12/24.
//

#include "vDimer.h"

using namespace QuanFloq;

// region Constructor/Destructor
template<typename T>
vDimer<T>::vDimer() = default;
// endregion
template<typename T>
T vDimer<T>::getT() const {
	return t;
}
template<typename T>
void vDimer<T>::setT( T tt ) {
	t = tt;
	T t2 = tt * sqrt(T(2.0f));
	for (auto it:{1, 4})
		this->H[it] = t2;
}
template<typename T>
T vDimer<T>::getV() const {
	return v;
}
template<typename T>
void vDimer<T>::setV( T tv ) {
	v = tv;
	this->H[0] = tv + U;
	this->H[5] = -tv + U;
}
template<typename T>
T vDimer<T>::getU() const {
	return U;
}
template<typename T>
void vDimer<T>::setU( T tU ) {
	U = tU;
	this->H[0] = v + tU;
	this->H[5] = -v + tU;
}

// region Initialize templates
//#ifdef BUILD_VIRTUAL
#ifdef BUILD_FLOAT
template
class QuanFloq::vDimer<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::vDimer<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::vDimer<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::vDimer<cdouble >;
#endif
//#endif
// endregion
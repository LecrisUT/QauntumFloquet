//
// Created by Le Minh Cristian on 2020/12/07.
//

#include "vHFDimer.h"

using namespace QuanFloq;

// region Constructor/Destructor
template<typename T>
vHFDimer<T>::vHFDimer() = default;
// endregion

template<typename T>
T* vHFDimer<T>::MakeUEx( T U ) {
	auto UEx = new T[4][4]();
	UEx[0][0] = U;
	UEx[3][3] = U;
	return reinterpret_cast<T*>(UEx);
}

template<typename T>
void vHFDimer<T>::setT( T t ) {
	this->t = t;
	this->h[1] = t;
	this->H[1] = t;
}
template<typename T>
void vHFDimer<T>::setV( T v ) {
	this->v = v;
	auto v2 = v * T(.5f);
	this->h[0] = v2;
	this->H[0] = v2;
	this->vh[0] = v2;
	this->h[2] = -v2;
	this->H[2] = -v2;
	this->vh[1] = -v2;
}
template<typename T>
void vHFDimer<T>::setU( T U ) {
	this->U = U;
	this->tensUEx[0] = U;
	this->tensUEx[7] = U;
}

// region Initialize templates
//#ifdef BUILD_VIRTUAL
#ifdef BUILD_FLOAT
template
class QuanFloq::vHFDimer<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::vHFDimer<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::vHFDimer<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::vHFDimer<cdouble >;
#endif
//#endif
// endregion
#include "vFloqHFDimer.h"

using namespace QuanFloq;

// region Constructor/Destructor
template<typename T>
vFloqHFDimer<T>::vFloqHFDimer() = default;
// endregion
template<typename T>
void vFloqHFDimer<T>::setT( T t, bool CalcS ) {
	this->t = t;
	T ct;
	if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>)
		ct = t;
	else
		ct = std::conj(t);
	this->h[1] = t;
	this->H[1] = t;
	this->h[2] = ct;
	this->H[2] = ct;
	if (CalcS)
		this->CalcSH();
}
template<typename T>
void vFloqHFDimer<T>::setV0( T v0, bool CalcS ) {
	this->v = v0;
	auto v02 = v0 * T(.5f);
	this->h[0] = v02;
	this->H[0] = v02;
	this->vh[0] = v02;
	this->h[3] = -v02;
	this->H[3] = -v02;
	this->vh[1] = -v02;
	if (CalcS)
		this->CalcSH();
}
template<typename T>
void vFloqHFDimer<T>::setV1( T v1, bool CalcS ) {
	this->v1 = v1;
	auto v14 = v1 * T(.25f);
	this->h[4] = v14;
	this->H[4] = v14;
	this->vh[2] = v14;
	this->h[7] = -v14;
	this->H[7] = -v14;
	this->vh[3] = -v14;
	if (CalcS)
		this->CalcSH();
}
template<typename T>
void vFloqHFDimer<T>::setU( T U, bool CalcS ) {
	this->U = U;
	this->tensUEx[0] = U;
	this->tensUEx[7] = U;
	if (CalcS)
		this->CalcSH();
}
template<typename T>
void vFloqHFDimer<T>::setT( T t ) {
//	vFloqDimer<T>::setT(t);
	this->setT(t, true);
}
template<typename T>
void vFloqHFDimer<T>::setV( T v0 ) {
//	vFloqDimer<T>::setV(v0);
	this->setV0(v0, true);
}
template<typename T>
void vFloqHFDimer<T>::setU( T U ) {
//	vFloqDimer<T>::setU(U);
	this->setU(U, true);
}

// region Initialize templates
//#ifdef BUILD_VIRTUAL
#ifdef BUILD_FLOAT
template
class QuanFloq::vFloqHFDimer<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::vFloqHFDimer<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::vFloqHFDimer<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::vFloqHFDimer<cdouble >;
#endif
//#endif
// endregion
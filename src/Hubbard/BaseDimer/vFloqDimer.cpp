//
// Created by Le Minh Cristian on 2020/12/24.
//

#include "vFloqDimer.h"

using namespace QuanFloq;
// region Constructor/Destructor
template<typename T>
vFloqDimer<T>::vFloqDimer() = default;
template<typename T>
void vFloqDimer<T>::setT( T t, bool CalcS ) {
	this->t = t;
	auto t2 = t * sqrt(T(2.0f));
	for (auto it:{1, 5})
		this->H[it] = t2;
	for (auto it:{3, 7})
		if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>)
			this->H[it] = t2;
		else
			this->H[it] = std::conj(t2);
	if (CalcS)
		this->CalcSH();
}
template<typename T>
inline T vFloqDimer<T>::getV0() const {
	return this->v;
}
template<typename T>
void vFloqDimer<T>::setV0( T v0, bool CalcS ) {
	this->v = v0;
	this->H[0] = v0 + this->U;
	this->H[8] = -v0 + this->U;
	if (CalcS)
		this->CalcSH();
}
template<typename T>
inline T vFloqDimer<T>::getV1() const {
	return v1;
}
template<typename T>
void vFloqDimer<T>::setV1( T tv1, bool CalcS ) {
	v1 = tv1;
	auto v12 = tv1 / T(2.0f);
	this->H[9] = v12;
	this->H[17] = -v12;
	if (CalcS)
		this->CalcSH();
}
template<typename T>
void vFloqDimer<T>::setU( T U, bool CalcS ) {
	this->U = U;
	this->H[0] = this->v + U;
	this->H[8] = -this->v + U;
	if (CalcS)
		this->CalcSH();
}
template<typename T>
inline void vFloqDimer<T>::setT( T t ) {
	this->setT(t, true);
}
template<typename T>
inline void vFloqDimer<T>::setV( T v0 ) {
	this->setV0(v0, true);
}
template<typename T>
inline void vFloqDimer<T>::setU( T U ) {
	this->setU(U, true);
}

// region Initialize templates
//#ifdef BUILD_VIRTUAL
#ifdef BUILD_FLOAT
template
class QuanFloq::vFloqDimer<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::vFloqDimer<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::vFloqDimer<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::vFloqDimer<cdouble >;
#endif
//#endif
// endregion
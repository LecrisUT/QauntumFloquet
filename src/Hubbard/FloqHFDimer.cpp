//
// Created by Le Minh Cristian on 2020/12/07.
//

#include "FloqHFDimer.h"

using namespace QuanFloq;
template<typename T>
FloqHFDimer<T>::FloqHFDimer() = default;
template<typename T>
FloqHFDimer<T>::FloqHFDimer( int nFH, int nF_max, T v0, T v1, T U, T w, T t ) :
		vHamil<T>(2, std::is_same_v<T, float> || std::is_same_v<T, double> ? Sym : Her, false),
		vHFHamil<T>(2, 1, 2, false),
		vFloqHamil<T>(nFH, nF_max, true),
		vFloqHFHamil<T>(1, true) {
	this->indvH[0] = this->H;
	this->indvH[1] = this->H + 3;
	vFloqHamil<T>::setW(w, false);
	FloqHFDimer<T>::setU(U, false);
	FloqHFDimer<T>::setT(t, false);
	FloqHFDimer<T>::setV0(v0, false);
	FloqHFDimer<T>::setV1(v1, false);
	vFloqHFHamil<T>::Initialize();
	this->CalcSH();
}
template<typename T>
void FloqHFDimer<T>::setT( T t, bool CalcS ) {
	Dimer<T>::t = t;
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
void FloqHFDimer<T>::setV0( T v0, bool CalcS ) {
	Dimer<T>::v = v0;
	auto v02 = v0 / T(2.0f);
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
void FloqHFDimer<T>::setV1( T v1, bool CalcS ) {
	FloqDimer<T>::v1 = v1;
	auto v14 = v1 / T(4.0f);
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
void FloqHFDimer<T>::setU( T U, bool CalcS ) {
	Dimer<T>::U = U;
	this->tensUEx[0] = U;
	this->tensUEx[7] = U;
	if (CalcS)
		this->CalcSH();
}
template<typename T>
void FloqHFDimer<T>::setT( T t ) {
//	FloqDimer<T>::setT(t);
	this->setT(t, true);
}
template<typename T>
void FloqHFDimer<T>::setV( T v0 ) {
//	FloqDimer<T>::setV(v0);
	this->setV0(v0, true);
}
template<typename T>
void FloqHFDimer<T>::setU( T U ) {
//	FloqDimer<T>::setU(U);
	this->setU(U, true);
}

template
class QuanFloq::FloqHFDimer<float>;
template
class QuanFloq::FloqHFDimer<double>;
template
class QuanFloq::FloqHFDimer<cfloat >;
template
class QuanFloq::FloqHFDimer<cdouble >;
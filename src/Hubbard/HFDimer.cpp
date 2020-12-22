//
// Created by Le Minh Cristian on 2020/12/07.
//

#include "HFDimer.h"

using namespace QuanFloq;

template<typename T>
HFDimer<T>::HFDimer() = default;
template<typename T>
HFDimer<T>::HFDimer( T v, T U, T t ) :
		vHamil<T>(2, std::is_same_v<T, float> || std::is_same_v<T, double> ? Sym : Her, true),
		vHFHamil<T>(2, 1, 2, true) {
	this->indvH[0] = this->H;
	this->indvH[1] = this->H + 2;
	HFDimer<T>::setU(U);
	HFDimer<T>::setT(t);
	HFDimer<T>::setV(v);
}

template<typename T>
T* HFDimer<T>::MakeUEx( T U ) {
	auto UEx = new T[4][4]();
	UEx[0][0] = U;
	UEx[3][3] = U;
	return reinterpret_cast<T*>(UEx);
}

template<typename T>
void HFDimer<T>::setT( T t ) {
	Dimer<T>::t = t;
	this->h[1] = t;
	this->H[1] = t;
}
template<typename T>
void HFDimer<T>::setV( T v ) {
	Dimer<T>::v = v;
	auto v2 = v / T(2.0f);
	this->h[0] = v2;
	this->H[0] = v2;
	this->vh[0] = v2;
	this->h[2] = -v2;
	this->H[2] = -v2;
	this->vh[1] = -v2;
}
template<typename T>
void HFDimer<T>::setU( T U ) {
	Dimer<T>::U = U;
	this->tensUEx[0] = U;
	this->tensUEx[7] = U;
}


template
class QuanFloq::HFDimer<float>;
template
class QuanFloq::HFDimer<double>;
template
class QuanFloq::HFDimer<cfloat >;
template
class QuanFloq::HFDimer<cdouble >;
//
// Created by Le Minh Cristian on 2020/12/07.
//

#include "Dimer.h"

using namespace QuanFloq;

//static constexpr auto DimerSZ = HamilSize(3);
template<typename T>
Dimer<T>::Dimer()= default;
template<typename T>
Dimer<T>::Dimer( T v, T U, T t ) :
		vHamil<T>(3, std::is_same_v<T, float> || std::is_same_v<T, double> ? Sym : Her, true) {
//	T t2 = tt * sqrt((T)2.0f);
//	T M[9] = {tU+tv, t2, 0.0f,
//	          t2, 0.0f, t2,
//	          0.0f, t2, tU - tv};
//	this->setH(M);
	Dimer<T>::setT(t);
	Dimer<T>::setU(U);
	Dimer<T>::setV(v);
}
template<typename T>
T Dimer<T>::getT() const {
	return t;
}
template<typename T>
void Dimer<T>::setT( T tt ) {
	t = tt;
	T t2 = tt * sqrt(T(2.0f));
	for (auto it:{1, 4})
		this->H[it] = t2;
}
template<typename T>
T Dimer<T>::getV() const {
	return v;
}
template<typename T>
void Dimer<T>::setV( T tv ) {
	v = tv;
	this->H[0] = tv + U;
	this->H[5] = -tv + U;
}
template<typename T>
T Dimer<T>::getU() const {
	return U;
}
template<typename T>
void Dimer<T>::setU( T tU ) {
	U = tU;
	this->H[0] = v + tU;
	this->H[5] = -v + tU;
}

template
class QuanFloq::Dimer<float>;
template
class QuanFloq::Dimer<double>;
template
class QuanFloq::Dimer<cfloat >;
template
class QuanFloq::Dimer<cdouble >;


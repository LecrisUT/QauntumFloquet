//
// Created by Le Minh Cristian on 2020/12/07.
//

#include "FloqDimer.h"

using namespace QuanFloq;

//static constexpr auto DimerSZ = HamilSize(3);
template<typename T>
FloqDimer<T>::FloqDimer() = default;
template<typename T>
FloqDimer<T>::FloqDimer( int nF_max, T v0, T v1, T U, T w, T t ) :
		vHamil<T>(3, std::is_same_v<T, float> || std::is_same_v<T, double> ? Sym : Her, false),
		vFloqHamil<T>(1, nF_max, true) {
//	T t2 = tt * sqrt(2);
//	T M[9] = {tU+tv, t2, 0.0f,
//	          t2, 0.0f, t2,
//	          0.0f, t2, tU - tv};
//	this->setH(M);
	FloqDimer<T>::setT(t, false);
	FloqDimer<T>::setU(U, false);
	FloqDimer<T>::setV0(v0, false);
	FloqDimer<T>::setV1(v1, false);
	this->setW(w, false);
	this->CalcSH();
}
template<typename T>
void FloqDimer<T>::setT( T t, bool CalcS ) {
	Dimer<T>::t = t;
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
inline T FloqDimer<T>::getV0() const {
	return this->v;
}
template<typename T>
void FloqDimer<T>::setV0( T v0, bool CalcS ) {
	Dimer<T>::v = v0;
	this->H[0] = v0 + Dimer<T>::U;
	this->H[8] = -v0 + Dimer<T>::U;
	if (CalcS)
		this->CalcSH();
}
template<typename T>
inline T FloqDimer<T>::getV1() const {
	return v1;
}
template<typename T>
void FloqDimer<T>::setV1( T tv1, bool CalcS ) {
	v1 = tv1;
	auto v12 = tv1 / T(2.0f);
	this->H[9] = v12;
	this->H[17] = -v12;
	if (CalcS)
		this->CalcSH();
}
template<typename T>
void FloqDimer<T>::setU( T U, bool CalcS ) {
	Dimer<T>::U = U;
	this->H[0] = this->v + U;
	this->H[8] = -this->v + U;
	if (CalcS)
		this->CalcSH();
}
template<typename T>
inline void FloqDimer<T>::setT( T t ) {
	this->setT(t, true);
}
template<typename T>
inline void FloqDimer<T>::setV( T v0 ) {
	this->setV0(v0, true);
}
template<typename T>
inline void FloqDimer<T>::setU( T U ) {
	this->setU(U, true);
}

template
class QuanFloq::FloqDimer<float>;
template
class QuanFloq::FloqDimer<double>;
template
class QuanFloq::FloqDimer<cfloat >;
template
class QuanFloq::FloqDimer<cdouble >;
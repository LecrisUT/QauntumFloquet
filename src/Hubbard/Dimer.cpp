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

template
class QuanFloq::Dimer<float>;
template
class QuanFloq::Dimer<double>;
template
class QuanFloq::Dimer<cfloat >;
template
class QuanFloq::Dimer<cdouble >;


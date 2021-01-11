//
// Created by Le Minh Cristian on 2020/12/07.
//

#include "FloqDimer.h"

using namespace QuanFloq;

template<typename T>
FloqDimer<T>::FloqDimer( int nF_max, T v0, T v1, T U, T w, T t ):
		vHamil<T>(3, std::is_same_v<T, float> || std::is_same_v<T, double> ? Sym : Her, false),
		vFloqHamil<T>(1, nF_max, true) {
//	T t2 = tt * sqrt(2);
//	T M[9] = {tU+tv, t2, 0.0f,
//	          t2, 0.0f, t2,
//	          0.0f, t2, tU - tv};
//	this->setH(M);
	vFloqDimer<T>::setT(t, false);
	vFloqDimer<T>::setU(U, false);
	vFloqDimer<T>::setV0(v0, false);
	vFloqDimer<T>::setV1(v1, false);
	this->setW(w, false);
	this->CalcSH();
}
// endregion

template
class QuanFloq::FloqDimer<float>;
template
class QuanFloq::FloqDimer<double>;
template
class QuanFloq::FloqDimer<cfloat >;
template
class QuanFloq::FloqDimer<cdouble >;
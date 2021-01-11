//
// Created by Le Minh Cristian on 2020/12/24.
//

#include "HFDimer.h"

using namespace QuanFloq;

template<typename T>
HFDimer<T>::HFDimer( T v, T U, T t ):
		vHamil<T>(2, std::is_same_v<T, float> || std::is_same_v<T, double> ? Sym : Her, true),
		vHFHamil<T>(2, 1, 2, true) {
	this->indvH[0] = this->H;
	this->indvH[1] = this->H + 2;
	setU(U);
	setT(t);
	setV(v);
}
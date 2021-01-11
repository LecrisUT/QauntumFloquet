//
// Created by Le Minh Cristian on 2020/12/07.
//
#include "../../BaseHamil/vHamil.h"
#include "../vDimer.h"

#ifndef QF_DIMER_H
#define QF_DIMER_H

namespace QuanFloq {
	template<typename T>
	class Dimer;
	using sDimer = Dimer<float>;
	using dDimer = Dimer<double>;
	using cDimer = Dimer<cfloat >;
	using zDimer = Dimer<cdouble >;
	// endregion

	template <typename T>
	class Dimer final :
			public vDimer<T>{
	public:
		Dimer( T v, T U, T t = 1.0f );
	};
}

#endif //QF_DIMER_H

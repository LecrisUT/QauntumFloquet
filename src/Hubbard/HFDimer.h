//
// Created by Le Minh Cristian on 2020/12/24.
//

#ifndef QF_HFDIMER_H
#define QF_HFDIMER_H

#include "../vHFDimer.h"

namespace QuanFloq {
	template<typename T>
	class HFDimer;
	using sHFDimer = HFDimer<float>;
	using dHFDimer = HFDimer<double>;
	using cHFDimer = HFDimer<cfloat >;
	using zHFDimer = HFDimer<cdouble >;
	template<typename T>
	class HFDimer final :
			public vHFDimer<T> {
	public:
		HFDimer( T v, T U, T t = 1.0f );

	};
}
#endif //QF_HFDIMER_H
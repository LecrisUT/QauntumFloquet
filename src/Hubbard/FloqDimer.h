//
// Created by Le Minh Cristian on 2020/12/07.
//
#include "../vFloqDimer.h"

#ifndef QF_FLOQDIMER_H
#define QF_FLOQDIMER_H

namespace QuanFloq {
	template<typename T>
	class FloqDimer;
	using sFloqDimer = FloqDimer<float>;
	using dFloqDimer = FloqDimer<double>;
	using cFloqDimer = FloqDimer<cfloat >;
	using zFloqDimer = FloqDimer<cdouble >;
	// endregion

	template<typename T>
	class FloqDimer final :
			public vFloqDimer<T>{
	public:
		FloqDimer( int nF_max, T v0, T v1, T U, T w, T t = 1.0f );

	};
}

#endif //QF_FLOQDIMER_H

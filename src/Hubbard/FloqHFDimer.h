//
// Created by Le Minh Cristian on 2020/12/07.
//

#ifndef QF_FLOQHFDIMER_H
#define QF_FLOQHFDIMER_H

#include "../vFloqHFDimer.h"
namespace QuanFloq {
	template<typename T>
	class FloqHFDimer;
	using sFloqHFDimer = FloqHFDimer<float>;
	using dFloqHFDimer = FloqHFDimer<double>;
	using cFloqHFDimer = FloqHFDimer<cfloat >;
	using zFloqHFDimer = FloqHFDimer<cdouble >;
	// endregion

	template<typename T>
	class FloqHFDimer final :
			public vFloqHFDimer<T> {
	public:
		FloqHFDimer( int nFH, int nF_max, T v0, T v1, T u, T w, T t = 1.0f );

	};
}


#endif //QF_FLOQHFDIMER_H

//
// Created by Le Minh Cristian on 2020/12/07.
//
#include <vFloqDimer.h>

#ifndef QF_FLOQDIMER_H
#define QF_FLOQDIMER_H

namespace QuanFloq {
	// region Class Declaration
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

	// region Define the built Template instances
#ifdef BUILD_FLOAT
	extern template
	class FloqDimer<float>;
#endif
#ifdef BUILD_DOUBLE
	extern template
	class FloqDimer<double>;
#endif
#ifdef BUILD_CFLOAT
	extern template
	class FloqDimer<cfloat>;
#endif
#ifdef BUILD_CDOUBLE
	extern template
	class FloqDimer<cdouble>;
#endif
	// endregion
}

#endif //QF_FLOQDIMER_H

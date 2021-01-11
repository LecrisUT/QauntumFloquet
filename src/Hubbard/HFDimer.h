//
// Created by Le Minh Cristian on 2020/12/24.
//

#ifndef QF_HFDIMER_H
#define QF_HFDIMER_H

#include <vHFDimer.h>

namespace QuanFloq {
	// region Class Declarations
	template<typename T>
	class HFDimer;
	using sHFDimer = HFDimer<float>;
	using dHFDimer = HFDimer<double>;
	using cHFDimer = HFDimer<cfloat >;
	using zHFDimer = HFDimer<cdouble >;
	// endregion

	template<typename T>
	class HFDimer final :
			public vHFDimer<T> {
	public:
		HFDimer( T v, T U, T t = 1.0f );
	};

	// region Define the built Template instances
#ifdef BUILD_FLOAT
	extern template
	class HFDimer<float>;
#endif
#ifdef BUILD_DOUBLE
	extern template
	class HFDimer<double>;
#endif
#ifdef BUILD_CFLOAT
	extern template
	class HFDimer<cfloat>;
#endif
#ifdef BUILD_CDOUBLE
	extern template
	class HFDimer<cdouble>;
#endif
	// endregion
}
#endif //QF_HFDIMER_H
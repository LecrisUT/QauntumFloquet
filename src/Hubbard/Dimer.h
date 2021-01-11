//
// Created by Le Minh Cristian on 2020/12/07.
//
#include <vDimer.h>

#ifndef QF_DIMER_H
#define QF_DIMER_H

namespace QuanFloq {
	// region Class Declarations
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

	// region Define the built Template instances
#ifdef BUILD_FLOAT
	extern template
	class Dimer<float>;
#endif
#ifdef BUILD_DOUBLE
	extern template
	class Dimer<double>;
#endif
#ifdef BUILD_CFLOAT
	extern template
	class Dimer<cfloat>;
#endif
#ifdef BUILD_CDOUBLE
	extern template
	class Dimer<cdouble>;
#endif
	// endregion
}

#endif //QF_DIMER_H

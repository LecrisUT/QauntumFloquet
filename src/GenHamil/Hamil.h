//
// Created by Le Minh Cristian on 2020/11/24.
//
#include <vHamil.h>

#ifndef QF_HAMIL_H
#define QF_HAMIL_H

namespace QuanFloq {
	// region Class Declarations
	template<typename T>
	class Hamil;
	using sHamil = Hamil<float>;
	using dHamil = Hamil<double>;
	using cHamil = Hamil<cfloat >;
	using zHamil = Hamil<cdouble >;
	// endregion

	// region Class Definitions
	template<typename T>
	class Hamil final :
			public vHamil<T> {
	public:
		explicit Hamil( int nH, Hamil_Sym Sym );
		Hamil( int nH, Hamil_Sym Sym, T* H );
	};
	// endregion

	// region Define the built Template instances
#ifdef BUILD_FLOAT
	extern template
	class Hamil<float>;
#endif
#ifdef BUILD_DOUBLE
	extern template
	class Hamil<double>;
#endif
#ifdef BUILD_CFLOAT
	extern template
	class Hamil<cfloat>;
#endif
#ifdef BUILD_CDOUBLE
	extern template
	class Hamil<cdouble>;
#endif
	// endregion
}
#endif //QF_HAMIL_H

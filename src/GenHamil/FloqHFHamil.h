//
// Created by Le Minh Cristian on 2020/11/30.
//
#include <vFloqHFHamil.h>

#ifndef QF_FLOQHFHAMIL_H
#define QF_FLOQHFHAMIL_H

namespace QuanFloq {

	// region Class Declarations
	template<typename T>
	class FloqHFHamil;
	using sFloqHFHamil = FloqHFHamil<float>;
	using dFloqHFHamil = FloqHFHamil<double>;
	using cFloqHFHamil = FloqHFHamil<cfloat >;
	using zFloqHFHamil = FloqHFHamil<cdouble >;
	// endregion

	template<typename T>
	class FloqHFHamil final :
			vFloqHFHamil<T> {
	public:
		FloqHFHamil( int nH, Hamil_Sym SSym, int nFh, int nFH, int nF_max, int nElec, int nOrb, int nUEx );
		FloqHFHamil( int nH, Hamil_Sym SSym, int nFh, int nFH, int nF_max, int nElec, int nOrb, int nUEx, T* h, T w, T* UEx );
		FloqHFHamil( int nH, Hamil_Sym SSym, int nFh, int nFH, int nF_max, int nElec, int nOrb );
		FloqHFHamil( int nH, Hamil_Sym SSym, int nFh, int nFH, int nF_max, int nElec, int nOrb, T* h, T w, T* UEx );
	};

	// region Define the built Template instances
#ifdef BUILD_FLOAT
	extern template
	class FloqHFHamil<float>;
#endif
#ifdef BUILD_DOUBLE
	extern template
	class FloqHFHamil<double>;
#endif
#ifdef BUILD_CFLOAT
	extern template
	class FloqHFHamil<cfloat>;
#endif
#ifdef BUILD_CDOUBLE
	extern template
	class FloqHFHamil<cdouble>;
#endif
	// endregion
}

#endif //QF_FLOQHFHAMIL_H

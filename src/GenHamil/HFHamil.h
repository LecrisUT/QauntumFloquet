#include <vHFHamil.h>

#ifndef QF_HFHAMIL_H
#define QF_HFHAMIL_H

namespace QuanFloq {

	// region Class Declarations
	template<typename T>
	class HFHamil;
	using sHFHamil = HFHamil<float>;
	using dHFHamil = HFHamil<double>;
	using cHFHamil = HFHamil<cfloat >;
	using zHFHamil = HFHamil<cdouble >;
	// endregion

	template<typename T>
	class HFHamil final :
			virtual public vHFHamil<T> {
	public:
		HFHamil( int nH, Hamil_Sym Sym, int nElec, int nOrb, int nUEx );
		HFHamil( int nH, Hamil_Sym Sym, int nElec, int nOrb, int nUEx, T* h, T* UEx );
		HFHamil( int nH, Hamil_Sym Sym, int nElec, int nOrb );
		HFHamil( int nH, Hamil_Sym Sym, int nElec, int nOrb, T* h, T* UEx );
	};

	// region Define the built Template instances
#ifdef BUILD_FLOAT
	extern template
	class HFHamil<float>;
#endif
#ifdef BUILD_DOUBLE
	extern template
	class HFHamil<double>;
#endif
#ifdef BUILD_CFLOAT
	extern template
	class HFHamil<cfloat>;
#endif
#ifdef BUILD_CDOUBLE
	extern template
	class HFHamil<cdouble>;
#endif
	// endregion
}

#endif //QF_HFHAMIL_H

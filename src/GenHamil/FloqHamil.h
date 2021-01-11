#include <vFloqHamil.h>

#ifndef QF_FLOQHAMIL_H
#define QF_FLOQHAMIL_H
namespace QuanFloq {

	// region Class Declarations
	template<typename T>
	class FloqHamil;
	using sFloqHamil = FloqHamil<float>;
	using dFloqHamil = FloqHamil<double>;
	using cFloqHamil = FloqHamil<cfloat >;
	using zFloqHamil = FloqHamil<cdouble >;
	// endregion

	template<typename T>
	class FloqHamil final :
			virtual public vFloqHamil<T> {
	public:
		FloqHamil( int nH, Hamil_Sym SSym, int nFH, int nF_max );
		FloqHamil( int nH, Hamil_Sym SSym, int nFH, int nF_max, T* H, T w );
	};
	// region Define the built Template instances
#ifdef BUILD_FLOAT
	extern template
	class FloqHamil<float>;
#endif
#ifdef BUILD_DOUBLE
	extern template
	class FloqHamil<double>;
#endif
#ifdef BUILD_CFLOAT
	extern template
	class FloqHamil<cfloat>;
#endif
#ifdef BUILD_CDOUBLE
	extern template
	class FloqHamil<cdouble>;
#endif
	// endregion
}
#endif

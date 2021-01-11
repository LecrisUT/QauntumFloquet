#include "FloqHamil.h"

using namespace QuanFloq;

// region Constructor/Destructor
template<typename T>
FloqHamil<T>::FloqHamil( int nH, Hamil_Sym SSym, int nFH, int nF_max ) :
		vHamil<T>(nH, SSym, false),
		vFloqHamil<T>(nFH, nF_max, true) { }
template<typename T>
FloqHamil<T>::FloqHamil( int nH, Hamil_Sym SSym, int nFH, int nF_max, T* H, T w ) :
		FloqHamil<T>(nH, SSym, nFH, nF_max) {
	vFloqHamil<T>::Initialize(H, w);
}
// endregion

// region Initialize templates
#ifdef BUILD_FLOAT
template
class QuanFloq::FloqHamil<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::FloqHamil<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::FloqHamil<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::FloqHamil<cdouble >;
#endif
// endregion
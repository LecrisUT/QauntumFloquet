#include "Hamil.h"

using namespace QuanFloq;
// region Constructors
template<typename T>
Hamil<T>::Hamil( int nH, Hamil_Sym Sym ) :
		vHamil<T>(nH, Sym, true) { }
template<typename T>
Hamil<T>::Hamil( int nH, Hamil_Sym Sym, T* H ):
		Hamil<T>(nH, Sym) {
	vHamil<T>::Initialize(H);
}
// endregion

// region Initialize templates
#ifdef BUILD_FLOAT
template
class QuanFloq::Hamil<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::Hamil<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::Hamil<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::Hamil<cdouble >;
#endif
// endregion
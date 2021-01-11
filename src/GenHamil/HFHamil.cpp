#include "HFHamil.h"
#include <cstdlib>
#include <mkl.h>

using namespace QuanFloq;

// region Constructor/Destructor
template<typename T>
HFHamil<T>::HFHamil( int nH, Hamil_Sym Sym, int nElec, int nOrb, int nUEx ):
		vHamil<T>(nH, Sym),
		vHFHamil<T>(nElec, nOrb, nUEx, true) { }
template<typename T>
HFHamil<T>::HFHamil( int nH, Hamil_Sym Sym, int nElec, int nOrb, int nUEx, T* h, T* UEx ) :
		HFHamil<T>(nH, Sym, nElec, nOrb, nUEx) {
	vHFHamil<T>::Initialize(h, UEx);
}
template<typename T>
HFHamil<T>::HFHamil( int nH, Hamil_Sym Sym, int nElec, int nOrb ) :
		HFHamil<T>(nH, Sym, nElec, nOrb, nH * nH) { }
template<typename T>
HFHamil<T>::HFHamil( int nH, Hamil_Sym Sym, int nElec, int nOrb, T* h, T* UEx ) :
		HFHamil<T>(nH, Sym, nElec, nOrb, nH * nH, h, UEx) { }
// endregion

// region Initialize templates
#ifdef BUILD_FLOAT
template
class QuanFloq::HFHamil<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::HFHamil<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::HFHamil<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::HFHamil<cdouble >;
#endif
// endregion
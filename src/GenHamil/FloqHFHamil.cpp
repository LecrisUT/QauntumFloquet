//
// Created by Le Minh Cristian on 2020/11/30.
//

#include "FloqHFHamil.h"
using namespace QuanFloq;

// region Constructor/Destructor
template<typename T>
FloqHFHamil<T>::FloqHFHamil( int nH, Hamil_Sym SSym, int nFh, int nFH, int nF_max, int nElec, int nOrb, int nUEx ) :
		vHamil<T>(nH, SSym, false),
		vFloqHamil<T>(nFH, nF_max, true),
		vHFHamil<T>(nElec, nOrb, nUEx, false),
		vFloqHFHamil<T>(nFh, true) {

}
template<typename T>
FloqHFHamil<T>::FloqHFHamil( int nH, Hamil_Sym SSym, int nFh, int nFH, int nF_max, int nElec, int nOrb, int nUEx, T* h,
                             T w, T* UEx ) :
		FloqHFHamil(nH, SSym, nFh, nFH, nF_max, nElec, nOrb, nUEx) {
	vFloqHFHamil<T>::Initialize(h, w, UEx);
}
template<typename T>
FloqHFHamil<T>::FloqHFHamil( int nH, Hamil_Sym SSym, int nFh, int nFH, int nF_max, int nElec, int nOrb ) :
		FloqHFHamil(nH, SSym, nFh, nFH, nF_max, nElec, nOrb, nH * nH) {
}
template<typename T>
FloqHFHamil<T>::FloqHFHamil( int nH, Hamil_Sym SSym, int nFh, int nFH, int nF_max, int nElec, int nOrb, T* h, T w,
                             T* UEx ) :
		FloqHFHamil(nH, SSym, nFh, nFH, nF_max, nElec, nOrb, nH * nH, h, w, UEx) {
}
// endregion

// region Initialize templates
#ifdef BUILD_FLOAT
template
class QuanFloq::FloqHFHamil<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::FloqHFHamil<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::FloqHFHamil<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::FloqHFHamil<cdouble >;
#endif
// endregion

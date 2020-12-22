#include "Hamil.h"
#include <mkl.h>

using namespace QuanFloq;
// region Instantiations
// endregion

// TODO: Switch to RowMajor
// region Constructor/Destructor
template<typename T>
vHamil<T>::vHamil() :
		nH(), Sym() { }
template<typename T>
vHamil<T>::vHamil( int nH, Hamil_Sym Sym, bool initM ) :
		nH(nH), Sym(Sym),
		H(initM ? Sym == Full ? new T[nH * nH]() : new T[nH * (nH + 1) / 2]() : nullptr),
		Psi(initM ? new T[nH * nH]() : nullptr),
		E(initM ? new T[nH]() : nullptr),
		PreHPsi(), PostPsiHPsi(), PrePsiHPsi(), PostHPsi() { }
template<typename T>
vHamil<T>::vHamil( int nH, Hamil_Sym Sym, T* H, T* Psi, T* E ) :
		nH(nH), Sym(Sym), H(H), Psi(Psi), E(E) { }

template<typename T>
Hamil<T>::Hamil( int nH, Hamil_Sym Sym ) :
		vHamil<T>(nH, Sym, true) { }
template<typename T>
Hamil<T>::Hamil( int nH, Hamil_Sym Sym, T* H ):
		Hamil<T>(nH, Sym) {
	vHamil<T>::Initialize(H);
}

template<typename T, const HamilSize& Sz>
tHamil<T, Sz>::tHamil()
		: vHamil<T>(Sz.nH, false) {
	vHamil<T>::H = H;
	vHamil<T>::Psi = Psi;
	vHamil<T>::E = E;
}
template<typename T, const HamilSize& Sz>
tHamil<T, Sz>::tHamil( T* tH )
		: tHamil() {
	this->setH(tH);
}
// endregion

// region Get/Set
// TODO: Unclear if should be inlined or declared in header
template<typename T>
T* vHamil<T>::getH() const {
	return H;
}
template<typename T>
void vHamil<T>::setH( T* ttH ) {
	auto tH = reinterpret_cast<T(*)[nH]>(ttH);
	int ind = 0;
	for (int i = 0; i < nH; i++)
		for (int j = Sym == Full ? 0 : i; j < nH; j++)
			H[ind++] = tH[i][j];
	// Equivalent to:
//			H[ind++] = tH[i * nH + j];
}
template<typename T>
void vHamil<T>::setH( T** tH ) {
	int ind = 0;
	for (int i = 0; i < nH; i++)
		for (int j = Sym == Full ? 0 : i; j < nH; j++)
			H[ind++] = tH[i][j];
}

template<typename T>
inline const T* vHamil<T>::getPsi() const {
	return Psi;
}
template<typename T>
inline const T* vHamil<T>::getE() const {
	return E;
}
// endregion


// region Main methods
template<typename T>
void vHamil<T>::Initialize( T* tH ) {
	setH(tH);
}
// TODO: Unclear how inline affects library creation
template<typename T>
inline T* vHamil<T>::HPsi() {
	return HPsi(Psi);
}
template<typename T>
inline T* vHamil<T>::HPsi( T* tPsi ) {
	auto tHPsi = new T[nH]();
	HPsi(tPsi, tHPsi);
	return tHPsi;
}
// Main interface call
// TODO: Unclear if can inline
template<typename T>
void vHamil<T>::HPsi( T* tPsi, T* tHPsi ) {
	for (auto func : PreHPsi)
		func(this, tPsi);
	mHPsi(tPsi, tHPsi);
	for (auto func : PostHPsi)
		func(this, tPsi, tHPsi);
}
// TODO: Should be safe to inline because not public
// TODO: Make sure devirtualization occurs for inline to make sense
template<typename T>
inline void vHamil<T>::mHPsi( T* tPsi, T* tHPsi ) {
	if constexpr (std::is_same_v<T, float>)
		cblas_sspmv(CblasRowMajor, CblasUpper, nH,
		            1.0f, H, tPsi, 1, 0.0f, tHPsi, 1);
	else if constexpr (std::is_same_v<T, double>)
		cblas_dspmv(CblasRowMajor, CblasUpper, nH,
		            1.0f, H, tPsi, 1, 0.0f, tHPsi, 1);
	else if constexpr (std::is_same_v<T, cfloat >)
		cblas_chpmv(CblasRowMajor, CblasUpper, nH,
		            &cfloat1, H, tPsi, 1, &cfloat0, tHPsi, 1);
	else if constexpr (std::is_same_v<T, cdouble >)
		cblas_zhpmv(CblasRowMajor, CblasUpper, nH,
		            &cdouble1, H, tPsi, 1, &cdouble0, tHPsi, 1);
	else
			static_assert(sizeof(T) != sizeof(T), "Type not supported");
}
//// region Specific implementations: mHPsi
//template<>
//void vHamil<float>::mHPsi( float* tPsi, float* tHPsi ) {
//	cblas_sspmv(CblasRowMajor, CblasUpper, nH,
//	            1.0f, H, tPsi, 1, 0.0f, tHPsi, 1);
//}
//template<>
//void vHamil<double>::mHPsi( double* tPsi, double* tHPsi ) {
//	cblas_dspmv(CblasRowMajor, CblasUpper, nH,
//	            1.0f, H, tPsi, 1, 0.0f, tHPsi, 1);
//}
//template<>
//void vHamil<cfloat >::mHPsi( cfloat* tPsi, cfloat* tHPsi ) {
//	cblas_chpmv(CblasRowMajor, CblasUpper, nH,
//	            &cfloat1, H, tPsi, 1, &cfloat0, tHPsi, 1);
//}
//template<>
//void vHamil<cdouble >::mHPsi( cdouble* tPsi, cdouble* tHPsi ) {
//	cblas_zhpmv(CblasRowMajor, CblasUpper, nH,
//	            &cdouble1, H, tPsi, 1, &cdouble0, tHPsi, 1);
//}
//// endregion

template<typename T>
inline T vHamil<T>::PsiHPsi() {
	return PsiHPsi(Psi);
}
template<typename T>
inline T vHamil<T>::PsiHPsi( T* tPsi ) {
	T tE;
	auto tHPsi = new T[nH];
	PsiHPsi(tPsi, &tE, tHPsi);
	return tE;
}
template<typename T>
void vHamil<T>::PsiHPsi( T* tPsi, T* tE, T* tHPsi ) {
	for (auto func : PrePsiHPsi)
		func(this, tPsi);
	mPsiHPsi(tPsi, tE, tHPsi);
	for (auto func : PostPsiHPsi)
		func(this, tPsi, tE, tHPsi);
}
template<typename T>
inline void vHamil<T>::mPsiHPsi( T* tPsi, T* tE, T* tHPsi ) {
	mHPsi(tPsi, tHPsi);
	if constexpr (std::is_same_v<T, float>)
		*tE = cblas_sdot(nH, tPsi, 1, tHPsi, 1);
	else if constexpr (std::is_same_v<T, double>)
		*tE = cblas_ddot(nH, tPsi, 1, tHPsi, 1);
	else if constexpr (std::is_same_v<T, cfloat >)
		cblas_cdotc_sub(nH, tPsi, 1, tHPsi, 1, tE);
	else if constexpr (std::is_same_v<T, cdouble >)
		cblas_zdotc_sub(nH, tPsi, 1, tHPsi, 1, tE);
	else
			static_assert(sizeof(T) != sizeof(T), "Type not supported");
}
//// region Specific implementations: mPsiHPsi
//template<>
//void vHamil<float>::mPsiHPsi( float* tPsi, float* tE, float* tHPsi ) {
//	mHPsi(tPsi, tHPsi);
//	*tE = cblas_sdot(nH, tPsi, 1, tHPsi, 1);
//}
//template<>
//void vHamil<double>::mPsiHPsi( double* tPsi, double* tE, double* tHPsi ) {
//	mHPsi(tPsi, tHPsi);
//	*tE = cblas_ddot(nH, tPsi, 1, tHPsi, 1);
//}
//template<>
//void vHamil<cfloat >::mPsiHPsi( cfloat* tPsi, cfloat* tE, cfloat* tHPsi ) {
//	mHPsi(tPsi, tHPsi);
//	cblas_cdotc_sub(nH, tPsi, 1, tHPsi, 1, tE);
//}
//template<>
//void vHamil<cdouble >::mPsiHPsi( cdouble* tPsi, cdouble* tE, cdouble* tHPsi ) {
//	mHPsi(tPsi, tHPsi);
//	cblas_zdotc_sub(nH, tPsi, 1, tHPsi, 1, tE);
//}
//// endregion

template<typename T>
T vHamil<T>::Overlap( T* Bra, T* Ket ) {
	if constexpr (std::is_same_v<T, float>)
		return cblas_sdot(nH, Bra, 1, Ket, 1);
	else if constexpr (std::is_same_v<T, double>)
		return cblas_ddot(nH, Bra, 1, Ket, 1);
	else if constexpr (std::is_same_v<T, cfloat >) {
		cfloat val;
		cblas_cdotc_sub(nH, Bra, 1, Ket, 1, &val);
		return val;
	} else if constexpr (std::is_same_v<T, cdouble >) {
		cfloat val;
		cblas_cdotc_sub(nH, Bra, 1, Ket, 1, &val);
		return val;
	} else
			static_assert(sizeof(T) != sizeof(T), "Type not supported");
}
//// region Specific implementations: Overlap
//template<>
//float vHamil<float>::Overlap( float* Bra, float* Ket ) {
//	return cblas_sdot(nH, Bra, 1, Ket, 1);
//}
//template<>
//double vHamil<double>::Overlap( double* Bra, double* Ket ) {
//	return cblas_ddot(nH, Bra, 1, Ket, 1);
//}
//template<>
//cfloat vHamil<cfloat >::Overlap( cfloat* Bra, cfloat* Ket ) {
//	cfloat val;
//	cblas_cdotc_sub(nH, Bra, 1, Ket, 1, &val);
//	return val;
//}
//template<>
//cdouble vHamil<cdouble >::Overlap( cdouble* Bra, cdouble* Ket ) {
//	cdouble val;
//	cblas_zdotc_sub(nH, Bra, 1, Ket, 1, &val);
//	return val;
//}
//// endregion
template<typename T>
void vHamil<T>::NormalizePsi( T* tPsi, bool FlagNorm ) {
	if constexpr (std::is_same_v<T, float>) {
		// float type
		if (FlagNorm) {
			auto norm = Psi[0] < 0.0f ? -1.0f : 1.0f / cblas_snrm2(nH, tPsi, 1);
			cblas_sscal(nH, norm, tPsi, 1);
		} else if (Psi[0] < 0.0f)
			cblas_sscal(nH, -1.0f, tPsi, 1);
	} else if constexpr (std::is_same_v<T, double>) {
		// double type
		if (FlagNorm) {
			auto norm = Psi[0] < 0.0f ? -1.0f : 1.0f / cblas_dnrm2(nH, tPsi, 1);
			cblas_dscal(nH, norm, tPsi, 1);
		} else if (Psi[0] < 0.0f)
			cblas_dscal(nH, -1.0f, tPsi, 1);
	} else if constexpr (std::is_same_v<T, cfloat >) {
		// complex<float> type
		if (FlagNorm) {
			auto norm = Psi[0].real() < 0.0f ? -1.0f : 1.0f / cblas_scnrm2(nH, tPsi, 1);
			cblas_csscal(nH, norm, tPsi, 1);
		} else if (Psi[0].real() < 0.0f)
			cblas_csscal(nH, -1.0f, tPsi, 1);
	} else if constexpr (std::is_same_v<T, cdouble >) {
		// complex<double> type
		if (FlagNorm) {
			auto norm = Psi[0].real() < 0.0f ? -1.0f : 1.0f / cblas_dznrm2(nH, tPsi, 1);
			cblas_zdscal(nH, norm, tPsi, 1);
		} else if (Psi[0].real() < 0.0f)
			cblas_zdscal(nH, -1.0f, tPsi, 1);
	} else
			static_assert(sizeof(T) != sizeof(T), "Type not supported");
}
//// region Specific implementations: NormalizePsi
//template<>
//void vHamil<float>::NormalizePsi( float* tPsi, bool FlagNorm ) {
//	if (FlagNorm) {
//		auto norm = Psi[0] < 0.0f ? -1.0f : 1.0f / cblas_snrm2(nH, tPsi, 1);
//		cblas_sscal(nH, norm, tPsi, 1);
//	} else if (Psi[0] < 0.0f)
//		cblas_sscal(nH, -1.0f, tPsi, 1);
//}
//template<>
//void vHamil<double>::NormalizePsi( double* tPsi, bool FlagNorm ) {
//	if (FlagNorm) {
//		auto norm = Psi[0] < 0.0f ? -1.0f : 1.0f / cblas_dnrm2(nH, tPsi, 1);
//		cblas_dscal(nH, norm, tPsi, 1);
//	} else if (Psi[0] < 0.0f)
//		cblas_dscal(nH, -1.0f, tPsi, 1);
//}
//template<>
//void vHamil<cfloat >::NormalizePsi( cfloat* tPsi, bool FlagNorm ) {
//	if (FlagNorm) {
//		auto norm = Psi[0].real() < 0.0f ? -1.0f : 1.0f / cblas_scnrm2(nH, tPsi, 1);
//		cblas_csscal(nH, norm, tPsi, 1);
//	} else if (Psi[0].real() < 0.0f)
//		cblas_csscal(nH, -1.0f, tPsi, 1);
//}
//template<>
//void vHamil<cdouble >::NormalizePsi( cdouble* tPsi, bool FlagNorm ) {
//	if (FlagNorm) {
//		auto norm = Psi[0].real() < 0.0f ? -1.0f : 1.0f / cblas_dznrm2(nH, tPsi, 1);
//		cblas_zdscal(nH, norm, tPsi, 1);
//	} else if (Psi[0].real() < 0.0f)
//		cblas_zdscal(nH, -1.0f, tPsi, 1);
//}
//// endregion
// endregion
template
class QuanFloq::vHamil<float>;
template
class QuanFloq::vHamil<double>;
template
class QuanFloq::vHamil<cfloat >;
template
class QuanFloq::vHamil<cdouble >;
template
class QuanFloq::Hamil<float>;
template
class QuanFloq::Hamil<double>;
template
class QuanFloq::Hamil<cfloat >;
template
class QuanFloq::Hamil<cdouble >;
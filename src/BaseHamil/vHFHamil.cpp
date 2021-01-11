#include "HFHamil.h"
#include <cstdlib>
#include <mkl.h>

using namespace QuanFloq;

// region Instantiation
template<typename T>
void copyhtoH( T* h, T* H, int n );
//template<>
//void vHFHamil<float>::mCalcUEx( float* Bra, float* Ket );
//template<>
//void vHFHamil<double>::mCalcUEx( double* Bra, double* Ket );
//template<>
//void vHFHamil<cfloat >::mCalcUEx( cfloat* Bra, cfloat* Ket );
//template<>
//void vHFHamil<cdouble >::mCalcUEx( cdouble* Bra, cdouble* Ket );
//template<>
//void vHFHamil<float>::mCalcUEx_p1( float* BraKet, float* acc );
//template<>
//void vHFHamil<double>::mCalcUEx_p1( double* BraKet, double* acc );
//template<>
//void vHFHamil<cfloat >::mCalcUEx_p1( cfloat* BraKet, cfloat* acc );
//template<>
//void vHFHamil<cdouble >::mCalcUEx_p1( cdouble* BraKet, cdouble* acc );
// endregion

// region Constructor/Destructor
template<typename T>
vHFHamil<T>::vHFHamil() :
		nElec(), nOrb(), nUEx_max(), nUEx(),
		h(), tensUEx(), vh(), indvH() { }
template<typename T>
vHFHamil<T>::vHFHamil( int nElec, int nOrb, bool initM ) :
		nElec(nElec), nOrb(nOrb), nUEx_max(this->nH * this->nH), nUEx(nUEx_max),
		h(initM ? new T[this->nH * this->nH] : nullptr),
		tensUEx(initM ? new T[nUEx_max * this->nH * this->nH] : nullptr),
		vh(initM ? new T[nUEx_max] : nullptr),
		indvH(initM ? new T* [nUEx_max] : nullptr) {
	this->PreHPsi.push_front(&vHFHamil<T>::UpdateH);
	this->PrePsiHPsi.push_front(&vHFHamil<T>::UpdateH);
}
template<typename T>
vHFHamil<T>::vHFHamil( int nElec, int nOrb, int nUEx, bool initM ) :
		nElec(nElec), nOrb(nOrb), nUEx_max(nUEx), nUEx(nUEx_max),
		h(initM ? new T[this->nH * this->nH] : nullptr),
		tensUEx(initM ? new T[nUEx_max * this->nH * this->nH] : nullptr),
		vh(initM ? new T[nUEx_max] : nullptr),
		indvH(initM ? new T* [nUEx_max] : nullptr) {
	this->PreHPsi.push_front(&vHFHamil<T>::UpdateH);
	this->PrePsiHPsi.push_front(&vHFHamil<T>::UpdateH);
}
template<typename T>
vHFHamil<T>::vHFHamil( int nElec, int nOrb, int nUEx, T* h, T* tensUEx, T* vh, T** indvH ):
		nElec(nElec), nOrb(nOrb), nUEx_max(nUEx), nUEx(nUEx_max),
		h(h), tensUEx(tensUEx), vh(vh), indvH(indvH) {
	this->PreHPsi.push_front(&vHFHamil<T>::UpdateH);
	this->PrePsiHPsi.push_front(&vHFHamil<T>::UpdateH);
}

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

#if __cplusplus >= 202002L
template<typename T, HFHamilSize Sz>
tHFHamil<T, Sz>::tHFHamil() : tHamil<T,Sz>() {
	vHFHamil<T>::h = h;
	vHFHamil<T>::indvH = indvH;
	vHFHamil<T>::vh = vh;
	vHFHamil<T>::tensUEx = tensUEx;
}
template<typename T, HFHamilSize Sz>
tHFHamil<T, Sz>::tHFHamil( T* th ):tHFHamil() {
	std::copy(th, th + Sz.nH2t, h);
	copyhtoH<T>(h, this->H, this->nH);
}
#endif

// endregion

// region Get/Set
template<typename T>
T* vHFHamil<T>::geth() const {
	return h;
}
template<typename T>
void vHFHamil<T>::seth( T* th ) {
	int ind = 0;
	int iUEx = 0;
	for (int i = 0; i < this->nH; i++)
		for (int j = this->Sym == Full ? 0 : i; j < this->nH; j++) {
			auto val = th[i * this->nH + j];
			this->H[ind] = val;
			h[ind] = val;
			if (this->H + ind == indvH[iUEx]) {
				vh[iUEx] = val;
				iUEx++;
			}
			ind++;
		}
	assert(iUEx == nUEx);
}
template<typename T>
void vHFHamil<T>::seth( T** th ) {
	int ind = 0;
	int iUEx = 0;
	for (int i = 0; i < this->nH; i++)
		for (int j = this->Sym == Full ? 0 : i; j < this->nH; j++) {
			auto val = th[i][j];
			this->H[ind] = val;
			h[ind] = val;
			if (this->H + ind == indvH[iUEx]) {
				vh[iUEx] = val;
				iUEx++;
			}
			ind++;
		}
	assert(iUEx == nUEx);
}
template<typename T>
T* vHFHamil<T>::getUEx( T* tPsi ) {
	int nH2 = this->nH * (this->nH + 1) / 2;
	auto UEx = new T[nH2]();
	getUEx(UEx, tPsi);
	return UEx;
}
template<typename T>
void vHFHamil<T>::getUEx( T* tUEx, T* tPsi ) {
	if (tPsi != nullptr)
		CalcUEx(tPsi, tUEx);
	for (int i = 0; i < nUEx; i++) {
		ptrdiff_t pos = indvH[i] - this->H;
		tUEx[pos] = this->H[pos] - vh[i];
	}
}
template<typename T>
inline void copyhtoH( T* h, T* H, int n ) {
	std::copy(h, h + n * (n + 1) / 2, H);
}
// endregion

// region Main Methods
template<typename T>
void vHFHamil<T>::Initialize( T* th, T* UEx, double Tresh ) {
	auto nH = this->nH;
	auto nH2 = nH * nH;
	auto nH2t = this->Sym == Full ? nH2 : nH * (nH + 1) / 2;
	// TODO: Can reduce size to triangular from symmetry: e.g. [nH*(nH+1)/2][nH*(nH+1)/2] (Need class and transformers)
	auto tUEx = reinterpret_cast<T(*)[nH2]>(UEx);
	int tnUEx = 0;
	int iInd[nH2t];
	// region Temporary hack until symmetry is implemented
	if (this->Sym == Full) {
		for (int i = 0; i < nH2; i++)
			iInd[i] = i;
	} else {
		int ind = 0;
		for (int i = 0; i < nH; i++)
			for (int j = i; j < nH; j++)
				iInd[ind++] = i * nH + j;
		assert (ind == nH2t);
	}
	// endregion
	int ind = 0;
	for (int ii = 0; ii < nH2t; ii++) {
		int i = iInd[ii];
		for (int j = 0; j < nH2; j++)
//			if (CalcTensor_p1(UEx[i][j])) {
			if (std::abs(tUEx[i][j]) > Tresh) {
				std::copy(tUEx[i], tUEx[i] + nH2, tensUEx + ind);
				// If reduced by symmetry can be simplified
				indvH[tnUEx] = this->H + ii;
				tnUEx++;
				ind += nH2;
				break;
			}
	}
	// Maximum allowed nUEx is determined in the constructor
	assert(tnUEx <= nUEx_max);
	nUEx = tnUEx;
	seth(th);
}
template<typename T>
void vHFHamil<T>::HPsi( T* tPsi, T* tHPsi, bool tupdateH ) {
	if (tupdateH)
		UpdateH(tPsi);
	this->mHPsi(tPsi, tHPsi);
}
template<typename T>
void vHFHamil<T>::PsiHPsi( T* tPsi, T* tE, T* tHPsi, bool tupdateH ) {
	if (tupdateH)
		UpdateH(tPsi);
	this->mPsiHPsi(tPsi, tE, tHPsi);
}

template<typename T>
void vHFHamil<T>::UpdateH( T* tPsi ) {
//	T vhUEx[nUEx] = {vh};
	T vhUEx[nUEx];
	// TODO: Change to std::copy
	for (int i = 0; i < nUEx; i++)
		vhUEx[i] = vh[i];
//	std::memcpy(vhUEx,vh,nUEx);
	CalcUEx(tPsi, vhUEx);
	// NOTE: Might not be compiler optimized (no SIMD etc.)
	for (auto iH = 0; iH < nUEx; iH++)
		*indvH[iH] = vhUEx[iH];
//	auto ivhUEx = vhUEx;
//	auto ipMax = indvH + nUEx;
//	for (auto ipH = indvH; ipH < ipMax; ipH++) {
//		**ipH = *ivhUEx;
//		ivhUEx++;
//	}
}
template<typename T>
void vHFHamil<T>::UpdateH( vHamil<T>* Hamil, T* tPsi ) {
	// TODO: Test if virtualization works
	auto tHamil = dynamic_cast<vHFHamil<T>*>( Hamil);
//	auto tHamil = static_cast<vHFHamil<T>*>( BaseHamil);
//	auto tHamil = reinterpret_cast<vHFHamil<T>*>( BaseHamil);
//	auto tHamil = (vHFHamil<T>*)BaseHamil;
	tHamil->UpdateH(tPsi);
}

template<typename T>
void vHFHamil<T>::CalcUEx( T* tPsi, T* acc ) {
	for (int iBra = 0; iBra < nOrb; iBra++)
		for (int iKet = 0; iKet < nOrb; iKet++)
			// Currently only U interaction
			// TODO: Standardize Orbital spin -> include exchange interaction
			mCalcUEx(tPsi + this->nH * iBra, tPsi + this->nH * iKet, acc);
//	auto PsiEnd = tPsi + nOrb * this->nH2F_max;
//	for (auto pBra = tPsi; pBra < PsiEnd; pBra += this->nH2F_max)
//		for (auto pKet = tPsi; pKet < PsiEnd; pKet += this->nH2F_max)
//			mCalcUEx(pBra, pKet, acc);
}
template<typename T>
inline void vHFHamil<T>::mCalcUEx( T* Bra, T* Ket, T* acc ) {
	T BraKet[this->nH * this->nH];
	if constexpr (std::is_same_v<T, float>)
		cblas_sger(CblasRowMajor, this->nH, this->nH, 1.0f, Bra, 1, Ket, 1, BraKet, this->nH);
	else if constexpr (std::is_same_v<T, double>)
		cblas_dger(CblasRowMajor, this->nH, this->nH, 1.0f, Bra, 1, Ket, 1, BraKet, this->nH);
	else if constexpr (std::is_same_v<T, cfloat >)
		cblas_cgerc(CblasRowMajor, this->nH, this->nH, &cfloat1, Bra, 1, Ket, 1, BraKet, this->nH);
	else if constexpr (std::is_same_v<T, cdouble >)
		cblas_zgerc(CblasRowMajor, this->nH, this->nH, &cdouble1, Bra, 1, Ket, 1, BraKet, this->nH);
	else
			static_assert(sizeof(T) != sizeof(T), "Type not supported");
	mCalcUEx_p1(BraKet, acc);
}
//// region Specific Implementation mCalUEx
//template<>
//void vHFHamil<float>::mCalcUEx( float* Bra, float* Ket, float* acc ) {
//	float BraKet[nH * nH];
//	cblas_sger(CblasRowMajor, nH, nH, 1.0f, Bra, 1, Ket, 1, BraKet, nH);
//	mCalcUEx_p1(BraKet, acc);
//}
//template<>
//void vHFHamil<double>::mCalcUEx( double* Bra, double* Ket, double* acc ) {
//	double BraKet[nH * nH];
//	cblas_dger(CblasRowMajor, nH, nH, 1.0f, Bra, 1, Ket, 1, BraKet, nH);
//	mCalcUEx_p1(BraKet, acc);
//}
//template<>
//void vHFHamil<cfloat >::mCalcUEx( cfloat* Bra, cfloat* Ket, cfloat* acc ) {
//	cfloat BraKet[nH * nH];
//	cblas_cgerc(CblasRowMajor, nH, nH, &cfloat1, Bra, 1, Ket, 1, BraKet, nH);
//	mCalcUEx_p1(BraKet, acc);
//}
//template<>
//void vHFHamil<cdouble >::mCalcUEx( cdouble* Bra, cdouble* Ket, cdouble* acc ) {
//	cdouble BraKet[nH * nH];
//	cblas_zgerc(CblasRowMajor, nH, nH, &cdouble1, Bra, 1, Ket, 1, BraKet, nH);
//	mCalcUEx_p1(BraKet, acc);
//}
//// endregion
template<typename T>
inline void vHFHamil<T>::mCalcUEx_p1( T* BraKet, T* acc ) {
	auto nH2 = this->nH * this->nH;
	if constexpr (std::is_same_v<T, float>)
		cblas_sgemv(CblasRowMajor, CblasNoTrans, nUEx, nH2,
		            1.0f, tensUEx, nH2, BraKet, 1, 1.0f, acc, 1);
	else if constexpr (std::is_same_v<T, double>)
		cblas_dgemv(CblasRowMajor, CblasNoTrans, nUEx, nH2,
		            1.0f, tensUEx, nH2, BraKet, 1, 1.0f, acc, 1);
	else if constexpr (std::is_same_v<T, cfloat >)
		cblas_cgemv(CblasRowMajor, CblasNoTrans, nUEx, nH2,
		            &cfloat1, tensUEx, nH2, BraKet, 1, &cfloat1, acc, 1);
	else if constexpr (std::is_same_v<T, cdouble >)
		cblas_zgemv(CblasRowMajor, CblasNoTrans, nUEx, nH2,
		            &cdouble1, tensUEx, nH2, BraKet, 1, &cdouble1, acc, 1);
	else
			static_assert(sizeof(T) != sizeof(T), "Type not supported");
}
//// region Specific Implementation CalUEx
//template<>
//void vHFHamil<float>::mCalcUEx_p1( float* BraKet, float* acc ) {
//	cblas_sgemv(CblasRowMajor, CblasNoTrans, nUEx, nH * nH,
//	            1.0f, tensUEx, nH * nH, BraKet, 1, 1.0f, acc, 1);
//}
//template<>
//void vHFHamil<double>::mCalcUEx_p1( double* BraKet, double* acc ) {
//	cblas_dgemv(CblasRowMajor, CblasNoTrans, nUEx, nH * nH,
//	            1.0f, tensUEx, nH * nH, BraKet, 1, 1.0f, acc, 1);
//}
//template<>
//void vHFHamil<cfloat >::mCalcUEx_p1( cfloat* BraKet, cfloat* acc ) {
//	cblas_cgemv(CblasRowMajor, CblasNoTrans, nUEx, nH * nH,
//	            &cfloat1, tensUEx, nH * nH, BraKet, 1, &cfloat1, acc, 1);
//}
//template<>
//void vHFHamil<cdouble >::mCalcUEx_p1( cdouble* BraKet, cdouble* acc ) {
//	cblas_zgemv(CblasRowMajor, CblasNoTrans, nUEx, nH * nH,
//	            &cdouble1, tensUEx, nH * nH, BraKet, 1, &cdouble1, acc, 1);
//}
//// endregion
// endregion

template
class QuanFloq::vHFHamil<float>;
template
class QuanFloq::vHFHamil<double>;
template
class QuanFloq::vHFHamil<cfloat >;
template
class QuanFloq::vHFHamil<cdouble >;
template
class QuanFloq::HFHamil<float>;
template
class QuanFloq::HFHamil<double>;
template
class QuanFloq::HFHamil<cfloat >;
template
class QuanFloq::HFHamil<cdouble >;
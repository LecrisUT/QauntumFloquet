#include "vHFHamil.h"
#include <cstdlib>
#include <mkl.h>

using namespace QuanFloq;

// region Instantiation
template<typename T>
void copyhtoH( T* h, T* H, int n );
// endregion

// region Constructor/Destructor
template<typename T>
vHFHamil<T>::vHFHamil() :
		nElec(), nOrb(), nUEx_max(), nUEx(),
		h(), tensUEx(), vh(), indvH() { }
template<typename T>
vHFHamil<T>::vHFHamil( int nElec, int nOrb, bool initM ) :
		nElec(nElec), nOrb(nOrb), nUEx_max(this->nH * this->nH), nUEx(nUEx_max),
		h(initM ? this->Sym == Full ? new T[this->nH * this->nH]() : new T[this->nH * (this->nH + 1) / 2]() : nullptr),
		tensUEx(initM ? new T[nUEx_max * this->nH * this->nH]() : nullptr),
		vh(initM ? new T[nUEx_max] : nullptr),
		indvH(initM ? new T* [nUEx_max] : nullptr) {
	this->PreHPsi.push_front(&vHFHamil<T>::UpdateH);
	this->PrePsiHPsi.push_front(&vHFHamil<T>::UpdateH);
}
template<typename T>
vHFHamil<T>::vHFHamil( int nElec, int nOrb, int nUEx, bool initM ) :
		nElec(nElec), nOrb(nOrb), nUEx_max(nUEx), nUEx(nUEx_max),
		h(initM ? this->Sym == Full ? new T[this->nH * this->nH]() : new T[this->nH * (this->nH + 1) / 2]() : nullptr),
		tensUEx(initM ? new T[nUEx_max * this->nH * this->nH]() : nullptr),
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
const T* vHFHamil<T>::geth() const {
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
	// TODO: Fix dynamically allocated memory allocating random amounts
	// TODO: Struct template use pre-determined memory
	T vhUEx[nUEx];
//	auto vhUEx = new T[nUEx]();
//	auto vhUEx = (T*)calloc(nUEx, sizeof(T));
	std::copy(vh, vh + nUEx, vhUEx);
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
	const auto nH2 = this->nH * this->nH;
	T BraKet[nH2];
//	auto BraKet = new T[nH2]();
//	auto BraKet = (T*)calloc(nH2, sizeof(T));
//	static const auto T0 = T(0.0f);
	std::fill(BraKet, BraKet + nH2, vHamil<T>::T0);
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
// endregion

template<typename T>
const T* vHFHamil<T>::geth( [[maybe_unused]] int m, [[maybe_unused]] int n ) const {
	return h;
}
template<typename T>
const T* vHFHamil<T>::getH( [[maybe_unused]] int m, [[maybe_unused]] int n ) const {
	return this->H;
}
template<typename T>
const T* vHFHamil<T>::getPsi( [[maybe_unused]] int m, [[maybe_unused]] int n ) const {
	return this->Psi;
}
template<typename T>
const T* vHFHamil<T>::getE( [[maybe_unused]] int n ) const {
	return this->E;
}
template<typename T>
T vHFHamil<T>::Overlap( T* Bra, T* Ket, [[maybe_unused]] int n ) {
	return this->Overlap(Bra, Ket);
}
template<typename T>
void vHFHamil<T>::NormalizePsi( T* tPsi, [[maybe_unused]] int n, bool FlagNorm ) {
	this->NormalizePsi(tPsi, FlagNorm);
}
template<typename T>
void vHFHamil<T>::HPsi( T* tPsi, T* tHPsi, [[maybe_unused]] int m, [[maybe_unused]] int n ) {
	this->HPsi(tPsi, tHPsi);
}
template<typename T>
void vHFHamil<T>::PsiHPsi( T* tPsi, T* tE, T* tHPsi, [[maybe_unused]] int m, [[maybe_unused]] int n ) {
	this->PsiHPsi(tPsi, tE, tHPsi);
}
template<typename T>
void vHFHamil<T>::getUEx( T* tUEx, int m, int n ) {
	std::fill(tUEx, tUEx + m * n, vHamil<T>::T0);
	getUEx(tUEx, nullptr);
}
template<typename T>
void vHFHamil<T>::getUEx( T* tUEx, T* tPsi, int m, int n,
						  [[maybe_unused]] int mPsi, [[maybe_unused]] int nPsi ) {
	std::fill(tUEx, tUEx + m * n, vHamil<T>::T0);
	getUEx(tUEx, tPsi);
}
template<typename T>
void vHFHamil<T>::UpdateH( T* tPsi, [[maybe_unused]] int m, [[maybe_unused]] int n ) {
	UpdateH(tPsi);
}

// region Initialize templates
//#ifdef BUILD_VIRTUAL
#ifdef BUILD_FLOAT
template
class QuanFloq::vHFHamil<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::vHFHamil<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::vHFHamil<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::vHFHamil<cdouble >;
#endif
//#endif
// endregion
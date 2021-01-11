//
// Created by Le Minh Cristian on 2020/11/30.
//

#include "vFloqHFHamil.h"
#include <mkl.h>

using namespace QuanFloq;

// region Constructor/Destructor
template<typename T>
vFloqHFHamil<T>::vFloqHFHamil():
		nFh(), csr_nUEx(), csr_nUEx_max(), csr_ind_UEx() { }
template<typename T>
vFloqHFHamil<T>::vFloqHFHamil( int nFh, int csr_nUEx, bool initM ) :
		nFh(nFh), csr_nUEx_max(csr_nUEx), csr_nUEx(csr_nUEx_max),
		csr_ind_UEx(initM ? new int[csr_nUEx_max] : nullptr) {
	// TODO: Fix or throw error if H dimension is smaller than h
//	if (this->nFH < nFh)
//		this->nFH = nFh;
	assert(nFh <= this->nFH);
	if (initM) {
		int nH2 = this->nH * this->nH;
		this->h = new T[nH2 * (nFh + 1)]();
		this->vh = new T[this->nUEx_max * (this->nFH + 1)]();
		if (this->indvH == nullptr)
			this->indvH = new T* [this->nUEx_max];
		if (this->tensUEx == nullptr)
			this->tensUEx = new T[this->nUEx_max * this->nH * this->nH]();
	}
}
template<typename T>
vFloqHFHamil<T>::vFloqHFHamil( int nFh, bool initM ) :
		nFh(nFh), csr_nUEx_max(CalcCsr_nUEx_max()), csr_nUEx(csr_nUEx_max),
		csr_ind_UEx(initM ? new int[csr_nUEx_max] : nullptr) {
	// TODO: Fix or throw error if H dimension is smaller than h
//	if (this->nFH < nFh)
//		this->nFH = nFh;
	assert(nFh <= this->nFH);
	if (initM) {
		int nH2 = this->nH * this->nH;
		this->h = new T[nH2 * (nFh + 1)]();
		this->vh = new T[this->nUEx_max * (this->nFH + 1)]();
		if (this->indvH == nullptr)
			this->indvH = new T* [this->nUEx_max];
		if (this->tensUEx == nullptr)
			this->tensUEx = new T[this->nUEx_max * this->nH * this->nH]();
	}
}
template<typename T>
vFloqHFHamil<T>::vFloqHFHamil( int nFh, int csr_nUEx, int* csr_ind_UEx ) :
		nFh(nFh), csr_nUEx_max(csr_nUEx), csr_nUEx(csr_nUEx_max), csr_ind_UEx(csr_ind_UEx) {
}

// endregion

// region Get/Set
template<typename T>
void vFloqHFHamil<T>::seth( T* th ) {
	auto nH = this->nH;
	auto nH2 = nH * nH;
	auto nFH = this->nFH;
	std::copy(th, th + nH2 * (nFh + 1), this->h);
	std::copy(th, th + nH2 * (nFh + 1), this->H);
	if (nFH > nFh)
		std::fill(this->H + nH2 * (nFh + 1), this->H + nH2 * (nFH + 1), vHamil<T>::T0);
	auto vh = reinterpret_cast<T(*)[this->nUEx_max]>(this->vh);
	auto h = reinterpret_cast<T(*)[nH2]>(this->h);
	for (int ind = 0; ind < this->nUEx; ind++) {
		ptrdiff_t pos = this->indvH[ind] - this->H;
		for (int iF = 0; iF <= nFh; iF++)
			vh[iF][ind] = h[iF][pos];
		for (int iF = nFh + 1; iF <= nFH; iF++)
			vh[iF][ind] = vHamil<T>::T0;
	}
}
template<typename T>
void vFloqHFHamil<T>::seth( T** th ) {
	// TODO: Implement this variation
	assert(false);
}
template<typename T>
T* vFloqHFHamil<T>::getUEx( T* tPsi ) {
	auto nH2 = this->nH * this->nH;
	auto UEx = new T[nH2 * (this->nFH + 1)]();
	getUEx(UEx, tPsi);
	return UEx;
}
template<typename T>
void vFloqHFHamil<T>::getUEx( T* tUEx, T* tPsi ) {
	auto nUEx = this->nUEx;
	auto nFH = this->nFH;
	auto nH2 = this->nH * this->nH;
	if (tPsi != nullptr)
		CalcUEx(tPsi, tUEx);
	int indUEX = 0;
	int indH = 0;
	for (int iF = 0; iF <= nFh; iF++) {
		for (int i = 0; i < nUEx; i++) {
			ptrdiff_t pos = this->indvH[i] + indH - this->H;
			tUEx[pos] = this->H[pos] - this->vh[indUEX + i];
		}
		indH += nH2;
		indUEX += nUEx;
	}
	for (int iF = nFh + 1; iF <= nFH; iF++) {
		for (int i = 0; i < nUEx; i++) {
			ptrdiff_t pos = this->indvH[i] + indH - this->H;
			tUEx[pos] = this->H[pos];
		}
		indH += nH2;
	}
}
// endregion

// region Main Methods
template<typename T>
void vFloqHFHamil<T>::Initialize() {
	int nH2 = this->nH * this->nH;
	int indUEx = 0;
	auto csr_nvalues = this->csr_nvalues;
	auto nUEx = this->nUEx;
	auto csr_values_map = this->csr_values_map;
	auto H = this->H;
	auto indvH = this->indvH;
	for (int ind = 0; ind < csr_nvalues; ind++) {
		ptrdiff_t pos = (csr_values_map[ind] - H) % nH2;
		auto pposH = H + pos;
		for (int i = 0; i < nUEx; i++)
			if (indvH[i] == pposH) {
				csr_ind_UEx[indUEx] = ind;
				indUEx++;
				break;
			}
	}
	csr_nUEx = indUEx;
	assert(indUEx <= csr_nUEx_max);
}
template<typename T>
void vFloqHFHamil<T>::Initialize( T* h, T w, T* UEx, double Tresh ) {
	this->setW(w, false);
	vHFHamil<T>::Initialize(h, UEx, Tresh);
	vFloqHFHamil<T>::Initialize();
}
template<typename T>
int vFloqHFHamil<T>::CalcCsr_nUEx_max() {
	return CalcCsr_nUEx(this->nUEx_max, this->nFH, this->n2F_max);
}

int QuanFloq::CalcCsr_nUEx( int nUEx, int nFH, int n2F_max ) {
	return nUEx * ((nFH + 1) * n2F_max - nFH * (nFH + 1) / 2);
}
template<typename T>
void vFloqHFHamil<T>::CalcSH_UEx() {
	for (int ind = 0; ind < csr_nUEx; ind++) {
		int i = csr_ind_UEx[ind];
		this->csr_values_SH[i] = *this->csr_values_map[i];
		this->csr_values_SHf[i] = *this->csr_values_map[i];
	}
	this->CalcSH_part();
#ifdef UPDATESHF
	[[maybe_unused]] auto res = mkl_sparse_copy(this->SH, this->SDescr, &this->SHf);
	assert(res == SPARSE_STATUS_SUCCESS);
#endif
	this->CalcSHf();
}


template<typename T>
void vFloqHFHamil<T>::UpdateH( T* tPsi ) {
	// TODO: If Sparse matrix is optimized, need to recreate SH
	T vhUEx[this->nUEx * (this->nFH + 1)];
//	auto vhUEx = new T[this->nUEx * (this->nFH + 1)]();
//	auto vhUEx = (T*)calloc(this->nUEx * (this->nFH + 1), sizeof(T));
	auto nH2 = this->nH * this->nH;
	std::copy(this->vh, this->vh + this->nUEx * (nFh + 1), vhUEx);
	std::fill(vhUEx + this->nUEx * (nFh + 1), vhUEx + this->nUEx * (this->nFH + 1), vHamil<T>::T0);
	for (auto iF = 0; iF <= this->nFH; iF++)
		CalcUEx(tPsi, vhUEx + iF * this->nUEx, iF);
	// NOTE: Might not be compiler optimized (no SIMD etc.)
	for (auto iF = 0; iF <= this->nFH; iF++)
		for (auto iH = 0; iH < this->nUEx; iH++)
			*(this->indvH[iH] + iF * nH2) = vhUEx[iH + iF * this->nUEx];
	CalcSH_UEx();
}

template<typename T>
void vFloqHFHamil<T>::CalcUEx( T* tPsi, T* acc, int iF ) {
	auto iFH = this->nH * iF;
	auto nF = this->n2F_max - iF;
	for (int iBra = 0; iBra < this->nOrb; iBra++)
		for (int iKet = 0; iKet < this->nOrb; iKet++)
			// Currently only Coulomb U
			// TODO: Standardize spin orbitals to calculate Exchange
			mCalcUEx(tPsi + this->nH * iBra, tPsi + this->nH * iKet + iFH, acc, nF);
}
template<typename T>
void vFloqHFHamil<T>::mCalcUEx( T* Bra, T* lKet, T* acc, int nF ) {
	const auto nH2 = this->nH * this->nH;
	T BraKet[nH2];
//	auto BraKet = new T[nH2]();
//	auto BraKet = (T*)calloc(nH2, sizeof(T));
	std::fill(BraKet, BraKet + nH2, vHamil<T>::T0);
	int iFH = 0;
	for (int iF = 0; iF < nF; iF++) {
		if constexpr (std::is_same_v<T, float>)
			cblas_sger(CblasRowMajor, this->nH, this->nH, 1.0f, Bra + iFH, 1, lKet + iFH, 1, BraKet, this->nH);
		else if constexpr (std::is_same_v<T, double>)
			cblas_dger(CblasRowMajor, this->nH, this->nH, 1.0f, Bra + iFH, 1, lKet + iFH, 1, BraKet, this->nH);
		else if constexpr (std::is_same_v<T, cfloat >)
			cblas_cgerc(CblasRowMajor, this->nH, this->nH, &cfloat1, Bra + iFH, 1, lKet + iFH, 1, BraKet, this->nH);
		else if constexpr (std::is_same_v<T, cdouble >)
			cblas_zgerc(CblasRowMajor, this->nH, this->nH, &cdouble1, Bra + iFH, 1, lKet + iFH, 1, BraKet, this->nH);
		else
				static_assert(sizeof(T) != sizeof(T), "Type not supported");
		iFH += this->nH;
	}
	this->mCalcUEx_p1(BraKet, acc);
}
// endregion

template<typename T>
const T* vFloqHFHamil<T>::geth( [[maybe_unused]] int m, [[maybe_unused]] int n, [[maybe_unused]] int nf ) const {
	return this->h;
}
template<typename T>
const T* vFloqHFHamil<T>::geth( [[maybe_unused]] int m, [[maybe_unused]] int n ) const {
	return this->h;
}
template<typename T>
void vFloqHFHamil<T>::getUEx( T* tUEx, int m, int n ) {
	std::fill(tUEx, tUEx + m * n, vHamil<T>::T0);
	this->getUEx(tUEx, nullptr);
}
template<typename T>
void vFloqHFHamil<T>::getUEx( T* tUEx, T* tPsi, int m, int n,
                              [[maybe_unused]] int mPsi, [[maybe_unused]] int nPsi ) {
	std::fill(tUEx, tUEx + m * n, vHamil<T>::T0);
	this->getUEx(tUEx, tPsi);
}
template<typename T>
void vFloqHFHamil<T>::UpdateH( T* tPsi, [[maybe_unused]] int m, [[maybe_unused]] int n ) {
	this->UpdateH(tPsi);
}
template<typename T>
const T* vFloqHFHamil<T>::getH( [[maybe_unused]] int m, [[maybe_unused]] int n, [[maybe_unused]] int nf ) const {
	return this->H;
}
template<typename T>
const T* vFloqHFHamil<T>::getH( [[maybe_unused]] int m, [[maybe_unused]] int n ) const {
	return this->H;
}
template<typename T>
const T* vFloqHFHamil<T>::getPsi( [[maybe_unused]] int m, [[maybe_unused]] int n ) const {
	return this->Psi;
}
template<typename T>
const T* vFloqHFHamil<T>::getE( [[maybe_unused]] int n ) const {
	return this->E;
}
template<typename T>
T vFloqHFHamil<T>::Overlap( T* Bra, T* Ket, [[maybe_unused]] int n ) {
	return this->Overlap(Bra, Ket);
}
template<typename T>
void vFloqHFHamil<T>::NormalizePsi( T* tPsi, [[maybe_unused]] int n, bool FlagNorm ) {
	this->NormalizePsi(tPsi, FlagNorm);
}
template<typename T>
void vFloqHFHamil<T>::HPsi( T* tPsi, T* tHPsi, [[maybe_unused]] int m, [[maybe_unused]] int n ) {
	this->HPsi(tPsi, tHPsi);
}
template<typename T>
void vFloqHFHamil<T>::PsiHPsi( T* tPsi, T* tE, T* tHPsi, [[maybe_unused]] int m, [[maybe_unused]] int n ) {
	this->PsiHPsi(tPsi, tE, tHPsi);
}

// region Initialize templates
//#ifdef BUILD_VIRTUAL
#ifdef BUILD_FLOAT
template
class QuanFloq::vFloqHFHamil<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::vFloqHFHamil<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::vFloqHFHamil<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::vFloqHFHamil<cdouble >;
#endif
//#endif
// endregion
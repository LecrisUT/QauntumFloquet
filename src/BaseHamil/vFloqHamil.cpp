#include "vFloqHamil.h"

using namespace QuanFloq;

// region Constructor/Destructor
template<typename T>
vFloqHamil<T>::vFloqHamil():
		nFH(), nF_max(), n2F_max(), nH2F_max(),
		diag_ind(), partialt(),
		SH(), SHf(), SDescr(), CalcHf(), ActiveSH(), csr_values_Active(),
		csr_nvalues(), csr_values_SH(), csr_values_SHf(),
		csr_columns(), csr_rows(), csr_diagInd(), csr_values_map() { }
template<typename T>
vFloqHamil<T>::vFloqHamil( int tnFH, int tnF_max, bool initM ) :
		nFH(tnFH), nF_max(tnF_max), n2F_max(2 * nF_max + 1), nH2F_max(this->nH * n2F_max),
		diag_ind(initM ? new int[n2F_max * this->nH] : nullptr),
		partialt(initM ? new T[n2F_max * this->nH]() : nullptr),
		SH(), SHf(), SDescr(), CalcHf(true), ActiveSH(SHf), csr_values_Active(csr_values_SHf),
		csr_nvalues((n2F_max * (nFH + 1) - nFH * (nFH + 1) / 2) * this->nH * this->nH),
		csr_values_SH(initM ? new T[csr_nvalues] : nullptr), csr_values_SHf(initM ? new T[csr_nvalues] : nullptr),
		csr_columns(initM ? new int[csr_nvalues] : nullptr), csr_rows(initM ? new int[nH2F_max + 1] : nullptr),
		csr_diagInd(initM ? new int[n2F_max * this->nH] : nullptr),
		csr_values_map(initM ? new T* [csr_nvalues] : nullptr) {
	if (!initM)
		return;
	auto nH = this->nH;
	this->H = new T[nH * nH * (nFH + 1)]();
	this->Psi = new T[nH * n2F_max * nH]();
	this->E = new T[nH]();
	for (int i = 0; i < nH2F_max; i++)
		diag_ind[i] = i;
	Initialize();
}
template<typename T>
vFloqHamil<T>::vFloqHamil( int nFH, int nF_max, int csrNvalues, int* csrColumns, int* csrRows, int* csrDiagInd,
                           T* csrValuesSh, T* csrValuesSHf, T** csrValuesMap, int* diagInd, T* partialt ) :
		nFH(nFH), nF_max(nF_max), n2F_max(2 * nF_max + 1), nH2F_max(this->nH * n2F_max),
		diag_ind(diagInd), partialt(partialt),
		SH(), SHf(), SDescr(), CalcHf(true), ActiveSH(SHf),
		csr_nvalues(csrNvalues), csr_values_SH(csrValuesSh), csr_values_SHf(csrValuesSHf),
		csr_columns(csrColumns), csr_rows(csrRows), csr_diagInd(csrDiagInd), csr_values_map(csrValuesMap) {
	for (int i = 0; i < nH2F_max; i++)
		diag_ind[i] = i;
	Initialize();
}

#if __cplusplus >= 202002L
template<typename T, FloqHamilSize Sz>
tFloqHamil<T, Sz>::tFloqHamil()
		: vFloqHamil<T>(Sz.nH, Sz.nFH, Sz.nF_max, diag_ind, partialt) {
	vHamil<T>::H = H;
	vHamil<T>::Psi = Psi;
	vHamil<T>::E = E;
	vFloqHamil<T>::partialt = partialt;
	vFloqHamil<T>::diag_ind = diag_ind;
}
template<typename T, FloqHamilSize Sz>
tFloqHamil<T, Sz>::tFloqHamil( T* tH )
		: tFloqHamil() {
	this->setH(tH);
}
#endif
// endregion

// region Get/Set
template<typename T>
T vFloqHamil<T>::getW() const {
	return w;
}
template<typename T>
void vFloqHamil<T>::setW( T tw, bool CalcS ) {
	w = tw;
	if (CalcS)
		CalcSHf();
}
template<typename T>
void vFloqHamil<T>::setH( T* tH, bool CalcS ) {
	std::copy(tH, tH + this->nH * this->nH * this->nFH, this->H);
	if (CalcS)
		CalcSH();
}
template<typename T>
inline void vFloqHamil<T>::setH( T* tH ) {
	setH(tH, true);
}
template<typename T>
bool vFloqHamil<T>::isCalcHf() const {
	return CalcHf;
}
template<typename T>
void vFloqHamil<T>::setCalcHf( bool tCalcHf ) {
	CalcHf = tCalcHf;
	ActiveSH = tCalcHf ? (sparse_matrix*)SHf : (sparse_matrix*)SH;
	csr_values_Active = tCalcHf ? csr_values_SHf : csr_values_SH;
}
// endregion

// region Main Methods
template<typename T>
void vFloqHamil<T>::Initialize() {
	if (initialized)
		return;
	if (this->Sym == Her)
		SDescr.type = SPARSE_MATRIX_TYPE_HERMITIAN;
	else
		SDescr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
	SDescr.mode = SPARSE_FILL_MODE_UPPER;
	SDescr.diag = SPARSE_DIAG_NON_UNIT;
	const auto nF1 = nFH + 1;
	const auto nFH1 = nF1 * this->nH;
	const auto nF2 = 2 * nF_max - nFH;
	int trim = 0;
	csr_rows[0] = 0;
	int ind = 0;
	T* HM[this->nH][nFH1];
	const auto nH2 = this->nH * this->nH;
	// Reshape the matrix [N][N*(nFH+1)] = [N][N][nFH+1]
	for (int i = 0; i < this->nH; i++)
		for (int iF = 0; iF <= nFH; iF++)
			for (int j = 0; j < this->nH; j++)
				HM[i][iF * this->nH + j] = this->H + (iF * nH2 + i * this->nH + j);
	// Compute the CSR Data
	for (int iF = 0; iF < n2F_max; iF++) {
		if (iF > nF2)
			trim++;
		auto iFH = iF * this->nH;
		auto nCol = (nF1 - trim) * this->nH;
		for (int i = 0; i < this->nH; i++) {
			// Calculate the current rows_end
			csr_rows[iFH + i + 1] = csr_rows[iFH + i] + nCol;
			for (int j = 0; j < nCol; j++) {
				// Save the column index
				csr_columns[ind + j] = iFH + j;
				csr_values_map[ind + j] = HM[i][j];
			}
			csr_diagInd[iF * this->nH + i] = ind + i;
			ind += nCol;
		}
	}
	initialized = true;
}
template<typename T>
void vFloqHamil<T>::Initialize( T* tH, T tw ) {
	Initialize();
	setH(tH, false);
	setW(tw, false);
	CalcSH();
}
template<typename T>
void vFloqHamil<T>::CalcSH() {
	for (int i = 0; i < csr_nvalues; i++) {
		csr_values_SH[i] = *csr_values_map[i];
		csr_values_SHf[i] = *csr_values_map[i];
	}
	CalcSH_part();
	// TODO: mkl_sparse_order (and mkl_sparse_optimize) delink csr_* arrays from their sparse handle
//	mkl_sparse_order(SH);
//	mkl_sparse_set_mv_hint(SH,SPARSE_OPERATION_NON_TRANSPOSE,SDescr,400);
//	mkl_sparse_optimize(SH);
#ifdef UPDATESHF
	[[maybe_unused]] auto res = mkl_sparse_copy(SH, SDescr, &SHf);
	assert(res == SPARSE_STATUS_SUCCESS);
#endif
	CalcSHf();
}
template<typename T>
inline void vFloqHamil<T>::CalcSH_part() {
	[[maybe_unused]] sparse_status_t res;
	if constexpr (std::is_same_v<T, float>)
		res = mkl_sparse_s_create_csr(&SH, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
		                              csr_rows + 1, csr_columns, csr_values_SH);
	else if constexpr (std::is_same_v<T, double>)
		res = mkl_sparse_d_create_csr(&SH, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
		                              csr_rows + 1, csr_columns, csr_values_SH);
	else if constexpr (std::is_same_v<T, cfloat >)
		res = mkl_sparse_c_create_csr(&SH, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
		                              csr_rows + 1, csr_columns, csr_values_SH);
	else if constexpr (std::is_same_v<T, cdouble >)
		res = mkl_sparse_z_create_csr(&SH, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
		                              csr_rows + 1, csr_columns, csr_values_SH);
	else
			static_assert(sizeof(T) != sizeof(T), "Type not supported");
	assert(res == SPARSE_STATUS_SUCCESS);
}
template<typename T>
void vFloqHamil<T>::CalcSHf() {
	int ind = 0;
	for (int iF = -nF_max; iF <= nF_max; iF++) {
		auto val = (T)-iF * w;
		for (int i = 0; i < this->nH; i++) {
			csr_values_SHf[csr_diagInd[ind + i]] = csr_values_SH[csr_diagInd[ind + i]] + val;
			partialt[ind + i] = val;
		}
		ind += this->nH;
	}
#ifndef OLDMKL
	CalcSHf_part();
//	mkl_sparse_optimize(SHf);
	setCalcHf(CalcHf);
#endif
}

#ifndef OLDMKL
template<typename T>
inline void vFloqHamil<T>::CalcSHf_part() {
	[[maybe_unused]] sparse_status_t res;
	if constexpr (std::is_same_v<T, float>)
#ifdef UPDATESHF
#ifdef FUTUREMKL
		res = mkl_sparse_s_update_values(SHf, this->nH * n2F_max, diag_ind, diag_ind, partialt);
#else
		for (int i = 1; i < nH2F_max; i++) {
			res = mkl_sparse_s_set_value(SHf, i, i, csr_values_SH[csr_diagInd[i]] + partialt[i]);
			assert(res == SPARSE_STATUS_SUCCESS);
		}
#endif
#else
		res = mkl_sparse_s_create_csr(&SHf, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
		                              csr_rows + 1, csr_columns, csr_values_SHf);
#endif
	else if constexpr (std::is_same_v<T, double>)
#ifdef UPDATESHF
#ifdef FUTUREMKL
		res = mkl_sparse_d_update_values(SHf, this->nH * n2F_max, diag_ind, diag_ind, partialt);
#else
		for (int i = 1; i < nH2F_max; i++) {
			res = mkl_sparse_d_set_value(SHf, i, i, csr_values_SH[csr_diagInd[i]] + partialt[i]);
			assert(res == SPARSE_STATUS_SUCCESS);
		}
#endif
#else
		res = mkl_sparse_d_create_csr(&SHf, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
		                              csr_rows + 1, csr_columns, csr_values_SHf);
#endif
	else if constexpr (std::is_same_v<T, cfloat >)
#ifdef UPDATESHF
#ifdef FUTUREMKL
		res = mkl_sparse_c_update_values(SHf, this->nH * n2F_max, diag_ind, diag_ind, partialt);
#else
		for (int i = 1; i < nH2F_max; i++) {
			res = mkl_sparse_c_set_value(SHf, i, i, csr_values_SH[csr_diagInd[i]] + partialt[i]);
			assert(res == SPARSE_STATUS_SUCCESS);
		}
#endif
#else
		res = mkl_sparse_c_create_csr(&SHf, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
		                              csr_rows + 1, csr_columns, csr_values_SHf);
#endif
	else if constexpr (std::is_same_v<T, cdouble >)
#ifdef UPDATESHF
#ifdef FUTUREMKL
		res = mkl_sparse_z_update_values(SHf, this->nH * n2F_max, diag_ind, diag_ind, partialt);
#else
		for (int i = 1; i < nH2F_max; i++) {
			res = mkl_sparse_z_set_value(SHf, i, i, csr_values_SH[csr_diagInd[i]] + partialt[i]);
			assert(res == SPARSE_STATUS_SUCCESS);
		}
#endif
#else
		res = mkl_sparse_z_create_csr(&SHf, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
		                              csr_rows + 1, csr_columns, csr_values_SHf);
#endif
	else
			static_assert(sizeof(T) != sizeof(T), "Type not supported");
	assert(res == SPARSE_STATUS_SUCCESS);
}
#endif
template<typename T>
T* vFloqHamil<T>::HPsi( T* tPsi ) {
	auto tHPsi = new T[n2F_max * this->nH];
	HPsi(tPsi, tHPsi);
	return tHPsi;
}
template<typename T>
inline void vFloqHamil<T>::mHPsi( T* tPsi, T* tHPsi ) {
	[[maybe_unused]] sparse_status_t res;
#ifdef OLDMKL
	if constexpr (std::is_same_v<T, cfloat > || std::is_same_v<T, cdouble >)
		static const char matdescra[] = {'H', 'U', 'N', 'R'};
	else
		static const char matdescra[] = {'S', 'U', 'N', 'R'};
#endif
	if constexpr (std::is_same_v<T, float>)
#ifdef OLDMKL
		//	mkl_cspblas_scsrsymv(&charU, &nH2F_max, csr_values_Active, csr_rows, csr_columns, tPsi, tHPsi);
		mkl_scsrmv(&charN, &nH2F_max, &nH2F_max, &float1, matdescra,
				   csr_values_Active, csr_columns, csr_rows, csr_rows + 1, tPsi, &float0, tHPsi);
#else
		res = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE,
		                      float1, ActiveSH, SDescr, tPsi, float0, tHPsi);
#endif
	else if constexpr (std::is_same_v<T, double>)
#ifdef OLDMKL
		mkl_dcsrmv(&charN, &nH2F_max, &nH2F_max, &double1, matdescra,
				   csr_values_Active, csr_columns, csr_rows, csr_rows + 1, tPsi, &double0, tHPsi);
#else
		res = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,
		                      double1, ActiveSH, SDescr, tPsi, double0, tHPsi);
#endif
	else if constexpr (std::is_same_v<T, cfloat >)
#ifdef OLDMKL
		mkl_ccsrmv(&charN, &nH2F_max, &nH2F_max, &cfloat1, matdescra,
				   csr_values_Active, csr_columns, csr_rows, csr_rows + 1, tPsi, &cfloat0, tHPsi);
#else
		res = mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE,
		                      cfloat1, ActiveSH, SDescr, tPsi, cfloat0, tHPsi);
#endif
	else if constexpr (std::is_same_v<T, cdouble >)
#ifdef OLDMKL
		mkl_zcsrmv(&charN, &nH2F_max, &nH2F_max, &cdouble1, matdescra,
				   csr_values_Active, csr_columns, csr_rows, csr_rows + 1, tPsi, &cdouble0, tHPsi);
#else
		res = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE,
		                      cdouble1, ActiveSH, SDescr, tPsi, cdouble0, tHPsi);
#endif
	else
			static_assert(sizeof(T) != sizeof(T), "Type not supported");
#ifndef OLDMKL
	assert(res == SPARSE_STATUS_SUCCESS);
#endif
}
template<typename T>
inline T vFloqHamil<T>::PsiHPsi( T* tPsi ) {
	T tE;
	auto tHPsi = new T[n2F_max * this->nH];
	PsiHPsi(tPsi, &tE, tHPsi);
	return tE;
}
#ifndef OLDMKL
template<typename T>
inline void vFloqHamil<T>::mPsiHPsi( T* tPsi, T* tE, T* tHPsi ) {
	[[maybe_unused]] sparse_status_t res;
	if constexpr (std::is_same_v<T, float>)
		res = mkl_sparse_s_dotmv(SPARSE_OPERATION_NON_TRANSPOSE,
		                         float1, ActiveSH, SDescr, tPsi, float0, tHPsi, tE);
	else if constexpr (std::is_same_v<T, double>)
		res = mkl_sparse_d_dotmv(SPARSE_OPERATION_NON_TRANSPOSE,
		                         double1, ActiveSH, SDescr, tPsi, double0, tHPsi, tE);
	else if constexpr (std::is_same_v<T, cfloat >)
		res = mkl_sparse_c_dotmv(SPARSE_OPERATION_NON_TRANSPOSE,
		                         cfloat1, ActiveSH, SDescr, tPsi, cfloat0, tHPsi, tE);
	else if constexpr (std::is_same_v<T, cdouble >)
		res = mkl_sparse_z_dotmv(SPARSE_OPERATION_NON_TRANSPOSE,
		                         cdouble1, ActiveSH, SDescr, tPsi, cdouble0, tHPsi, tE);
	else
			static_assert(sizeof(T) != sizeof(T), "Type not supported");
	assert(res == SPARSE_STATUS_SUCCESS);
}
#endif
template<typename T>
T vFloqHamil<T>::Overlap( T* Bra, T* Ket ) {
	if constexpr (std::is_same_v<T, float>)
		return cblas_sdot(nH2F_max, Bra, 1, Ket, 1);
	else if constexpr (std::is_same_v<T, double>)
		return cblas_ddot(nH2F_max, Bra, 1, Ket, 1);
	else if constexpr (std::is_same_v<T, cfloat >) {
		cfloat val;
		cblas_cdotc_sub(nH2F_max, Bra, 1, Ket, 1, &val);
		return val;
	} else if constexpr (std::is_same_v<T, cdouble >) {
		cdouble val;
		cblas_zdotc_sub(nH2F_max, Bra, 1, Ket, 1, &val);
		return val;
	} else
			static_assert(sizeof(T) != sizeof(T), "Type not supported");
}
template<typename T>
void vFloqHamil<T>::NormalizePsi( T* tPsi, bool FlagNorm ) {
	if constexpr (std::is_same_v<T, float>) {
		// float type
		if (FlagNorm) {
			auto norm = tPsi[this->nH * nF_max] < 0.0f ? -1.0f : 1.0f / cblas_snrm2(nH2F_max, tPsi, 1);
			cblas_sscal(nH2F_max, norm, tPsi, 1);
		} else if (tPsi[this->nH * nF_max] < 0.0f)
			cblas_sscal(nH2F_max, -1.0f, tPsi, 1);
	} else if constexpr (std::is_same_v<T, double>) {
		// double type
		if (FlagNorm) {
			auto norm = tPsi[this->nH * nF_max] < 0.0f ? -1.0f : 1.0f / cblas_dnrm2(nH2F_max, tPsi, 1);
			cblas_dscal(nH2F_max, norm, tPsi, 1);
		} else if (tPsi[this->nH * nF_max] < 0.0f)
			cblas_dscal(nH2F_max, -1.0f, tPsi, 1);
	} else if constexpr (std::is_same_v<T, cfloat >) {
		// cfloat type
		if (FlagNorm) {
			auto norm = tPsi[this->nH * nF_max].real() < 0.0f ? -1.0f : 1.0f / cblas_scnrm2(nH2F_max, tPsi, 1);
			cblas_csscal(nH2F_max, norm, tPsi, 1);
		} else if (tPsi[this->nH * nF_max].real() < 0.0f)
			cblas_csscal(nH2F_max, -1.0f, tPsi, 1);
	} else if constexpr (std::is_same_v<T, cdouble >) {
		// cdouble type
		if (FlagNorm) {
			auto norm = tPsi[this->nH * nF_max].real() < 0.0f ? -1.0f : 1.0f / cblas_dznrm2(nH2F_max, tPsi, 1);
			cblas_zdscal(nH2F_max, norm, tPsi, 1);
		} else if (tPsi[this->nH * nF_max].real() < 0.0f)
			cblas_zdscal(nH2F_max, -1.0f, tPsi, 1);
	} else
			static_assert(sizeof(T) != sizeof(T), "Type not supported");
}

template<typename T>
const T* vFloqHamil<T>::getH( [[maybe_unused]] int m, [[maybe_unused]] int n, [[maybe_unused]] int nf ) const {
	return this->H;
}
template<typename T>
const T* vFloqHamil<T>::getH( [[maybe_unused]] int m, [[maybe_unused]] int n ) const {
	return this->H;
}
template<typename T>
const T* vFloqHamil<T>::getPsi( [[maybe_unused]] int m, [[maybe_unused]] int n ) const {
	return this->Psi;
}
template<typename T>
const T* vFloqHamil<T>::getE( [[maybe_unused]] int n ) const {
	return this->E;
}
template<typename T>
T vFloqHamil<T>::Overlap( T* Bra, T* Ket, [[maybe_unused]] int n ) {
	return this->Overlap(Bra, Ket);
}
template<typename T>
void vFloqHamil<T>::NormalizePsi( T* tPsi, [[maybe_unused]] int n, bool FlagNorm ) {
	this->NormalizePsi(tPsi, FlagNorm);
}
template<typename T>
void vFloqHamil<T>::HPsi( T* tPsi, T* tHPsi, [[maybe_unused]] int m, [[maybe_unused]] int n ) {
	this->HPsi(tPsi, tHPsi);
}
template<typename T>
void vFloqHamil<T>::PsiHPsi( T* tPsi, T* tE, T* tHPsi, [[maybe_unused]] int m, [[maybe_unused]] int n ) {
	this->PsiHPsi(tPsi, tE, tHPsi);
}

// region Initialize templates
//#ifdef BUILD_VIRTUAL
#ifdef BUILD_FLOAT
template
class QuanFloq::vFloqHamil<float>;
#endif
#ifdef BUILD_DOUBLE
template
class QuanFloq::vFloqHamil<double>;
#endif
#ifdef BUILD_CFLOAT
template
class QuanFloq::vFloqHamil<cfloat >;
#endif
#ifdef BUILD_CDOUBLE
template
class QuanFloq::vFloqHamil<cdouble >;
#endif
//#endif
// endregion
#include "FloqHamil.h"

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

template<typename T>
FloqHamil<T>::FloqHamil( int nH, Hamil_Sym SSym, int nFH, int nF_max ) :
		vHamil<T>(nH, SSym, false),
		vFloqHamil<T>(nFH, nF_max, true) { }
template<typename T>
FloqHamil<T>::FloqHamil( int nH, Hamil_Sym SSym, int nFH, int nF_max, T* H, T w ) :
		FloqHamil<T>(nH, SSym, nFH, nF_max) {
	vFloqHamil<T>::Initialize(H, w);
}

template<typename T, const FloqHamilSize& Sz>
tFloqHamil<T, Sz>::tFloqHamil()
		: vFloqHamil<T>(Sz.nH, Sz.nFH, Sz.nF_max, diag_ind, partialt) {
	vHamil<T>::H = H;
	vHamil<T>::Psi = Psi;
	vHamil<T>::E = E;
	vFloqHamil<T>::partialt = partialt;
	vFloqHamil<T>::diag_ind = diag_ind;
}
template<typename T, const FloqHamilSize& Sz>
tFloqHamil<T, Sz>::tFloqHamil( T* tH )
		: tFloqHamil() {
	this->setH(tH);
}

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
// region Specific implementations: CalcSH
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
//template<>
//void vFloqHamil<float>::CalcSH_part() {
//	[[maybe_unused]] auto res = mkl_sparse_s_create_csr(&SH, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
//	                                                    csr_rows + 1, csr_columns, csr_values_SH);
//	assert(res == SPARSE_STATUS_SUCCESS);
//}
//template<>
//void vFloqHamil<double>::CalcSH_part() {
//	[[maybe_unused]] auto res = mkl_sparse_d_create_csr(&SH, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
//	                                                    csr_rows + 1, csr_columns, csr_values_SH);
//	assert(res == SPARSE_STATUS_SUCCESS);
//}
//template<>
//void vFloqHamil<cfloat >::CalcSH_part() {
//	[[maybe_unused]] auto res = mkl_sparse_c_create_csr(&SH, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
//	                                                    csr_rows + 1, csr_columns, csr_values_SH);
//	assert(res == SPARSE_STATUS_SUCCESS);
//}
//template<>
//void vFloqHamil<cdouble >::CalcSH_part() {
//	[[maybe_unused]] auto res = mkl_sparse_z_create_csr(&SH, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
//	                                                    csr_rows + 1, csr_columns, csr_values_SH);
//	assert(res == SPARSE_STATUS_SUCCESS);
//}
// endregion
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
// region Specific implementations: CalcSHf
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
//template<>
//void vFloqHamil<float>::CalcSHf_part() {
//#ifdef UPDATESHF
//#ifdef FUTUREMKL
//	[[maybe_unused]] auto res = mkl_sparse_s_update_values(SHf, this->nH * n2F_max, diag_ind, diag_ind, partialt);
//	assert(res == SPARSE_STATUS_SUCCESS);
//#else
//	for (int i = 1; i < nH2F_max; i++) {
//		[[maybe_unused]] auto res = mkl_sparse_s_set_value(SHf, i, i, csr_values_SH[csr_diagInd[i]] + partialt[i]);
//		assert(res == SPARSE_STATUS_SUCCESS);
//	}
//#endif
//#else
//	[[maybe_unused]] auto res = mkl_sparse_s_create_csr(&SHf, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
//	                                                    csr_rows + 1, csr_columns, csr_values_SHf);
//	assert(res == SPARSE_STATUS_SUCCESS);
//#endif
//}
//template<>
//void vFloqHamil<double>::CalcSHf_part() {
//#ifdef UPDATESHF
//#ifdef FUTUREMKL
//	[[maybe_unused]] auto res = mkl_sparse_d_update_values(SHf, this->nH * n2F_max, diag_ind, diag_ind, partialt);
//	assert(res == SPARSE_STATUS_SUCCESS);
//#else
//	for (int i = 1; i < nH2F_max; i++) {
//		[[maybe_unused]] auto res = mkl_sparse_d_set_value(SHf, i, i, csr_values_SH[csr_diagInd[i]] + partialt[i]);
//		assert(res == SPARSE_STATUS_SUCCESS);
//	}
//#endif
//#else
//	[[maybe_unused]] auto res = mkl_sparse_d_create_csr(&SHf, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
//	                                                    csr_rows + 1, csr_columns, csr_values_SHf);
//	assert(res == SPARSE_STATUS_SUCCESS);
//#endif
//}
//template<>
//void vFloqHamil<cfloat >::CalcSHf_part() {
//#ifdef UPDATESHF
//#ifdef FUTUREMKL
//	[[maybe_unused]] auto res = mkl_sparse_c_update_values(SHf, this->nH * n2F_max, diag_ind, diag_ind, partialt);
//	assert(res == SPARSE_STATUS_SUCCESS);
//#else
//	for (int i = 1; i < nH2F_max; i++) {
//		[[maybe_unused]] auto res = mkl_sparse_c_set_value(SHf, i, i, csr_values_SH[csr_diagInd[i]] + partialt[i]);
//		assert(res == SPARSE_STATUS_SUCCESS);
//	}
//#endif
//#else
//	[[maybe_unused]] auto res = mkl_sparse_c_create_csr(&SHf, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
//	                                                    csr_rows + 1, csr_columns, csr_values_SHf);
//	assert(res == SPARSE_STATUS_SUCCESS);
//#endif
//}
//template<>
//void vFloqHamil<cdouble >::CalcSHf_part() {
//#ifdef UPDATESHF
//#ifdef FUTUREMKL
//	[[maybe_unused]] auto res = mkl_sparse_z_update_values(SHf, this->nH * n2F_max, diag_ind, diag_ind, partialt);
//	assert(res == SPARSE_STATUS_SUCCESS);
//#else
//	for (int i = 1; i < nH2F_max; i++) {
//		[[maybe_unused]] auto res = mkl_sparse_z_set_value(SHf, i, i, csr_values_SH[csr_diagInd[i]] + partialt[i]);
//		assert(res == SPARSE_STATUS_SUCCESS);
//	}
//#endif
//#else
//	[[maybe_unused]] auto res = mkl_sparse_z_create_csr(&SHf, SPARSE_INDEX_BASE_ZERO, nH2F_max, nH2F_max, csr_rows,
//	                                                    csr_rows + 1, csr_columns, csr_values_SHf);
//	assert(res == SPARSE_STATUS_SUCCESS);
//#endif
//}
// endregion
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
//template<>
//void vFloqHamil<float>::mHPsi( float* tPsi, float* tHPsi ) {
//#ifdef OLDMKL
//	//	mkl_cspblas_scsrsymv(&charU, &nH2F_max, csr_values_Active, csr_rows, csr_columns, tPsi, tHPsi);
//		static const char matdescra[] = {'S', 'U', 'N', 'C'};
//		mkl_scsrmv(&charN, &nH2F_max, &nH2F_max, &float1, matdescra,
//				   csr_values_Active, csr_columns, csr_rows, csr_rows + 1, tPsi, &float0, tHPsi);
//#else
//	[[maybe_unused]] auto res = mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE,
//	                                            float1, ActiveSH, SDescr, tPsi, float0, tHPsi);
//	assert(res == SPARSE_STATUS_SUCCESS);
//#endif
//}
//template<>
//void vFloqHamil<double>::mHPsi( double* tPsi, double* tHPsi ) {
//#ifdef OLDMKL
//	static const char matdescra[] = {'S', 'U', 'N', 'C'};
//	mkl_dcsrmv(&charN, &nH2F_max, &nH2F_max, &double1, matdescra,
//			   csr_values_Active, csr_columns, csr_rows, csr_rows + 1, tPsi, &double0, tHPsi);
//#else
//	[[maybe_unused]] auto res = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,
//	                                            double1, ActiveSH, SDescr, tPsi, double0, tHPsi);
//	assert(res == SPARSE_STATUS_SUCCESS);
//#endif
//}
//template<>
//void vFloqHamil<cfloat >::mHPsi( cfloat* tPsi, cfloat* tHPsi ) {
//#ifdef OLDMKL
//	static const char matdescra[] = {'H', 'U', 'N', 'C'};
//	mkl_ccsrmv(&charN, &nH2F_max, &nH2F_max, &cfloat1, matdescra,
//			   csr_values_Active, csr_columns, csr_rows, csr_rows + 1, tPsi, &cfloat0, tHPsi);
//#else
//	[[maybe_unused]] auto res = mkl_sparse_c_mv(SPARSE_OPERATION_NON_TRANSPOSE,
//	                                            cfloat1, ActiveSH, SDescr, tPsi, cfloat0, tHPsi);
//	assert(res == SPARSE_STATUS_SUCCESS);
//#endif
//}
//template<>
//void vFloqHamil<cdouble >::mHPsi( cdouble* tPsi, cdouble* tHPsi ) {
//#ifdef OLDMKL
//	static const char matdescra[] = {'H', 'U', 'N', 'C'};
//	mkl_zcsrmv(&charN, &nH2F_max, &nH2F_max, &cdouble1, matdescra,
//			   csr_values_Active, csr_columns, csr_rows, csr_rows + 1, tPsi, &cdouble0, tHPsi);
//#else
//	[[maybe_unused]] auto res = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE,
//	                                            cdouble1, ActiveSH, SDescr, tPsi, cdouble0, tHPsi);
//	assert(res == SPARSE_STATUS_SUCCESS);
//#endif
//}
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
//// region Specific implementations: PsiHPsi
//template<>
//void vFloqHamil<float>::mPsiHPsi( float* tPsi, float* tE, float* tHPsi ) {
//	mkl_sparse_s_dotmv(SPARSE_OPERATION_NON_TRANSPOSE,
//	                   float1, ActiveSH, SDescr, tPsi, float0, tHPsi, tE);
//}
//template<>
//void vFloqHamil<double>::mPsiHPsi( double* tPsi, double* tE, double* tHPsi ) {
//	mkl_sparse_d_dotmv(SPARSE_OPERATION_NON_TRANSPOSE,
//	                   double1, ActiveSH, SDescr, tPsi, double0, tHPsi, tE);
//}
//template<>
//void vFloqHamil<cfloat >::mPsiHPsi( cfloat* tPsi, cfloat* tE, cfloat* tHPsi ) {
//	mkl_sparse_c_dotmv(SPARSE_OPERATION_NON_TRANSPOSE,
//	                   cfloat1, ActiveSH, SDescr, tPsi, cfloat0, tHPsi, tE);
//}
//template<>
//void vFloqHamil<cdouble >::mPsiHPsi( cdouble* tPsi, cdouble* tE, cdouble* tHPsi ) {
//	mkl_sparse_z_dotmv(SPARSE_OPERATION_NON_TRANSPOSE,
//	                   cdouble1, ActiveSH, SDescr, tPsi, cdouble0, tHPsi, tE);
//}
//// endregion
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
//// region Specific implementations: Overlap
//template<>
//float vFloqHamil<float>::Overlap( float* Bra, float* Ket ) {
//	return cblas_sdot(nH2F_max, Bra, 1, Ket, 1);
//}
//template<>
//double vFloqHamil<double>::Overlap( double* Bra, double* Ket ) {
//	return cblas_ddot(nH2F_max, Bra, 1, Ket, 1);
//}
//template<>
//cfloat vFloqHamil<cfloat >::Overlap( cfloat* Bra, cfloat* Ket ) {
//	cfloat val;
//	cblas_cdotc_sub(nH2F_max, Bra, 1, Ket, 1, &val);
//	return val;
//}
//template<>
//cdouble vFloqHamil<cdouble >::Overlap( cdouble* Bra, cdouble* Ket ) {
//	cdouble val;
//	cblas_zdotc_sub(nH2F_max, Bra, 1, Ket, 1, &val);
//	return val;
//}
//// endregion
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
//// region Specific implementations: NormalizePsi
//template<>
//void vFloqHamil<float>::NormalizePsi( float* tPsi, bool FlagNorm ) {
//	if (FlagNorm) {
//		auto norm = tPsi[nH * nF_max] < 0.0f ? -1.0f : 1.0f / cblas_snrm2(nH2F_max, tPsi, 1);
//		cblas_sscal(nH2F_max, norm, tPsi, 1);
//	} else if (tPsi[nH * nF_max] < 0.0f)
//		cblas_sscal(nH2F_max, -1.0f, tPsi, 1);
//}
//template<>
//void vFloqHamil<double>::NormalizePsi( double* tPsi, bool FlagNorm ) {
//	if (FlagNorm) {
//		auto norm = Psi[nH * nF_max] < 0.0f ? -1.0f : 1.0f / cblas_dnrm2(nH2F_max, tPsi, 1);
//		cblas_dscal(nH2F_max, norm, tPsi, 1);
//	} else if (Psi[nH * nF_max] < 0.0f)
//		cblas_dscal(nH2F_max, -1.0f, tPsi, 1);
//}
//template<>
//void vFloqHamil<cfloat >::NormalizePsi( cfloat* tPsi, bool FlagNorm ) {
//	if (FlagNorm) {
//		auto norm = Psi[nH * nF_max].real() < 0.0f ? -1.0f : 1.0f / cblas_scnrm2(nH2F_max, tPsi, 1);
//		cblas_csscal(nH2F_max, norm, tPsi, 1);
//	} else if (Psi[nH * nF_max].real() < 0.0f)
//		cblas_csscal(nH2F_max, -1.0f, tPsi, 1);
//}
//template<>
//void vFloqHamil<cdouble >::NormalizePsi( cdouble* tPsi, bool FlagNorm ) {
//	if (FlagNorm) {
//		auto norm = Psi[nH * nF_max].real() < 0.0f ? -1.0f : 1.0f / cblas_dznrm2(nH2F_max, tPsi, 1);
//		cblas_zdscal(nH, norm, tPsi, 1);
//	} else if (Psi[nH * nF_max].real() < 0.0f)
//		cblas_zdscal(nH, -1.0f, tPsi, 1);
//}
//// endregion

template
class QuanFloq::vFloqHamil<float>;
template
class QuanFloq::vFloqHamil<double>;
template
class QuanFloq::vFloqHamil<cfloat >;
template
class QuanFloq::vFloqHamil<cdouble >;
template
class QuanFloq::FloqHamil<float>;
template
class QuanFloq::FloqHamil<double>;
template
class QuanFloq::FloqHamil<cfloat >;
template
class QuanFloq::FloqHamil<cdouble >;
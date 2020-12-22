//
// Created by Le Minh Cristian on 2020/12/07.
//

#include "Test.h"
//#include <cstdio>
#include <iostream>
#include <iomanip>
#include <mkl.h>

using namespace QuanFloq;
using namespace std;

static const int termWidth = 40;

template<typename T>
void Test<T>::PrintFloqMatrix( T* Mtr, int n, int nF, bool flag_Full ) {
//	T* tMtr = (T*)H->getH();
	auto tMtr = reinterpret_cast<T(*)[n][n]>(Mtr);
	if (flag_Full) {
		// M[2*nF+1][n][n]
//		auto M = new T[n * n * (2 * nF + 1)];
		auto M = reinterpret_cast<T(*)[n][n]>(new T[n * n * (2 * nF + 1)]());
		std::copy(tMtr[0][0], (T*)tMtr + nF * n * n, (T*)M[nF]);
//		int ind = 0;
		for (int k = 0; k <= nF; k++)
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++) {
//					M[i * n + (nF + k) * n + j] = tMtr[ind++];
					M[nF + k][i][j] = tMtr[k][i][j];
				}
		PrintMatrix(reinterpret_cast<T*>(M), n * (2 * nF + 1), n, n, n);
	} else
		PrintMatrix(reinterpret_cast<T*>(tMtr), n * (nF + 1), n, n, n);
}
template<typename T>
void Test<T>::PrintTriangMatrix( T* Mtr, int n ) {
	auto Mfull = new T[n * n]();
	int ind = 0;
	for (int i = 0; i < n; i++)
		for (int j = i; j < n; j++)
			Mfull[i * n + j] = Mtr[ind++];
	PrintMatrix(Mfull, n, n);
}
template<typename T>
inline void Test<T>::PrintMatrix( T* M, int m, int n, bool Transpose ) {
	PrintMatrix(M, m, n, m, n, Transpose);
}

template<typename T>
void Test<T>::PrintMatrix( T* M, int m, int n, int mSep, int nSep, bool Transpose ) {
	for (int i = 0; i < m; i++) {
		if (!(i % mSep))
			cout << string(termWidth, '-') << endl;
		for (int j = 0; j < n; j++) {
			if (!(j % nSep))
				cout << " | ";
			if (Transpose)
				cout << setw(8) << M[j * m + i] << "  ";
			else
				cout << setw(8) << M[i * n + j] << "  ";
		}
		cout << " | ";
		cout << endl;
	}
	cout << string(termWidth, '-') << endl;
}

template<typename T>
inline void Test<T>::PrintH( vHamil<T>* H ) {
	PrintTriangMatrix(H->getH(), H->nH);
}
template<typename T>
void Test<T>::PrintH( vFloqHamil<T>* H, bool flag_Full ) {
	PrintFloqMatrix(H->getH(), H->nH, H->nFH, flag_Full);
}
template<typename T>
inline void Test<T>::Printh( vHFHamil<T>* H ) {
	PrintTriangMatrix(H->geth(), H->nH);
}
template<typename T>
void Test<T>::Printh( vFloqHFHamil<T>* H, bool flag_Full ) {
	PrintFloqMatrix(H->geth(), H->nH, H->nFh, flag_Full);
}
template<typename T>
void Test<T>::PrintUEx( vHFHamil<T>* H ) {
	PrintTriangMatrix(H->getUEx(), H->nH);
}
template<typename T>
void Test<T>::PrintUEx( vFloqHFHamil<T>* H, bool flag_Full ) {
	PrintFloqMatrix(H->getUEx(), H->nH, H->nFH, flag_Full);
}
template<typename T>
void Test<T>::PrintSH( vFloqHamil<T>* H ) {
	auto nH2F_max = H->nH2F_max;
	auto nH = H->nH;
	auto* M = new T[nH2F_max * nH2F_max]();
#ifdef OLDMKL
	ConvertCsrToDense(H,H->csr_values_SH,M);
#else
	ConvertSparseToDense(H->SH, nH2F_max, H->SDescr, M);
#endif
	PrintMatrix(M, nH2F_max, nH2F_max, nH, nH);
}
template<typename T>
void Test<T>::PrintSHf( vFloqHamil<T>* H ) {
	auto nH2F_max = H->nH2F_max;
	auto nH = H->nH;
	auto M = new T[nH2F_max * nH2F_max]();
#ifdef OLDMKL
	ConvertCsrToDense(H,H->csr_values_SHf,M);
#else
	ConvertSparseToDense(H->SHf, nH2F_max, H->SDescr, M);
#endif
	PrintMatrix(M, nH2F_max, nH2F_max, nH, nH);
}
#ifdef OLDMKL
template<typename T>
void Test<T>::ConvertCsrToDense( vFloqHamil<T>* H, T* values, T* M ) {
	static_assert(sizeof(T) != sizeof(T), "Type not supported");
}
template<>
void Test<float>::ConvertCsrToDense( vFloqHamil<float>* H, float* values, float* M ) {
	static const int job[]={1,0,0,2,0,0};
	auto nH2F_max = H->nH2F_max;
	int info;
	mkl_sdnscsr(job,&nH2F_max,&nH2F_max,M,&nH2F_max,values,H->csr_columns,H->csr_rows,&info);
	assert(info==0);
}
template<>
void Test<double>::ConvertCsrToDense( vFloqHamil<double>* H, double* values, double* M ) {
	static const int job[]={1,0,0,2,0,0};
	auto nH2F_max = H->nH2F_max;
	int info;
	mkl_ddnscsr(job,&nH2F_max,&nH2F_max,M,&nH2F_max,values,H->csr_columns,H->csr_rows,&info);
	assert(info==0);
}
template<>
void Test<cfloat>::ConvertCsrToDense( vFloqHamil<cfloat>* H, cfloat* values, cfloat* M ) {
	static const int job[]={1,0,0,2,0,0};
	auto nH2F_max = H->nH2F_max;
	int info;
	mkl_cdnscsr(job,&nH2F_max,&nH2F_max,M,&nH2F_max,values,H->csr_columns,H->csr_rows,&info);
	assert(info==0);
}
template<>
void Test<cdouble>::ConvertCsrToDense( vFloqHamil<cdouble>* H, cdouble* values, cdouble* M ) {
	static const int job[]={1,0,0,2,0,0};
	auto nH2F_max = H->nH2F_max;
	int info;
	mkl_zdnscsr(job,&nH2F_max,&nH2F_max,M,&nH2F_max,values,H->csr_columns,H->csr_rows,&info);
	assert(info==0);
}
#else
template<typename T>
void Test<T>::ConvertSparseToDense( sparse_matrix* S, int N, matrix_descr Descr, T* M ) {
	static_assert(sizeof(T) != sizeof(T), "Type not supported");
}
template<>
void Test<float>::ConvertSparseToDense( sparse_matrix* S, int N, matrix_descr Descr, float* M ) {
	Descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	Descr.mode = SPARSE_FILL_MODE_FULL;
	auto* eye = new float[N * N]();
	for (int i = 0; i < N; i++)
		eye[i * N + i] = 1.0f;
	auto res = mkl_sparse_s_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, S, Descr, SPARSE_LAYOUT_ROW_MAJOR,
	                           eye, N, N, 0.0f, M, N);
	assert(res == SPARSE_STATUS_SUCCESS);
}
template<>
void Test<double>::ConvertSparseToDense( sparse_matrix* S, int N, matrix_descr Descr, double* M ) {
	auto* eye = new double[N * N]();
	for (int i = 0; i < N; i++)
		eye[i * N + i] = 1.0f;
	auto res = mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1.0f, S, Descr, SPARSE_LAYOUT_ROW_MAJOR,
	                           eye, N, N, 0.0f, M, N);
	assert(res == SPARSE_STATUS_SUCCESS);
}
template<>
void Test<cfloat >::ConvertSparseToDense( sparse_matrix* S, int N, matrix_descr Descr, cfloat* M ) {
	auto* eye = new cfloat[N * N]();
	for (int i = 0; i < N; i++)
		eye[i * N + i] = cfloat1;
	auto res = mkl_sparse_c_mm(SPARSE_OPERATION_NON_TRANSPOSE, cfloat1, S, Descr, SPARSE_LAYOUT_ROW_MAJOR,
	                           eye, N, N, cfloat0, M, N);
	assert(res == SPARSE_STATUS_SUCCESS);
}
template<>
void Test<cdouble >::ConvertSparseToDense( sparse_matrix* S, int N, matrix_descr Descr, cdouble* M ) {
	auto* eye = new cdouble[N * N]();
	for (int i = 0; i < N; i++)
		eye[i * N + i] = cdouble1;
	auto res = mkl_sparse_z_mm(SPARSE_OPERATION_NON_TRANSPOSE, cdouble1, S, Descr, SPARSE_LAYOUT_ROW_MAJOR,
	                           eye, N, N, cdouble0, M, N);
	assert(res == SPARSE_STATUS_SUCCESS);
}
#endif
template
class QuanFloq::Test<float>;
template
class QuanFloq::Test<double>;
template
class QuanFloq::Test<cfloat >;
template
class QuanFloq::Test<cdouble >;
/*template<typename T>
void QuanFloq::PrintH( vFloqHamil<T> H, bool flag_Full ) {
	auto nFH = flag_Full ? 2 * H.nFH + 1 : H.nFH + 1;
	PrintMatrix(H.H, H.nH * nFH, H.nH, H.nH, H.nH);
}
template<typename T>
void QuanFloq::Printh( vHFHamil<T> H ) {
	PrintMatrix(H.h, H.nH, H.nH);
}
template<typename T>
void QuanFloq::Printh( vFloqHFHamil<T> H, bool flag_Full ) {
	auto nFh = flag_Full ? 2 * H.nFh + 1 : H.nFh + 1;
	PrintMatrix(H.h, H.nH * nFh, H.nH, H.nH, H.nH);
}
template<typename T>
void QuanFloq::PrintUEx( vHFHamil<T> H ) {
	PrintMatrix(H.getUEx(), H.nH, H.nH);
}
template<typename T>
void QuanFloq::PrintUEx( vFloqHFHamil<T> H, bool flag_Full ) {
	auto nFH = flag_Full ? 2 * H.nFH + 1 : H.nFH + 1;
	PrintMatrix(H.h, H.nH * nFH, H.nH, H.nH, H.nH);
}
template<typename T>
void QuanFloq::PrintSH( vFloqHamil<T> H ) {

}
template<typename T>
void QuanFloq::PrintSHf( vFloqHamil<T> H ) {

}
*/
#include "StaticHamil.h"
#include <mkl_cblas.h>
#include <cstdlib>
#include <cstring>

using namespace QuanFloq;

// Main Hamiltonian object
template<int N>
dStaticHamil<N>::dStaticHamil( StaticHamil* H ) : staticH(H), dHamil<N>(H) { }
template<int N>
dStaticHamil<N>::dStaticHamil( double* H ) : dStaticHamil(new StaticHamil(H)) { }

template<int NValue>
void dStaticHamil<NValue>::PsiHPsi( double* Psi, double* E, double* HPsi ) {
	this->HPsi(Psi, HPsi);
//	cblas_dspmv(CblasColMajor, CblasUpper, N, 1.0f, (double*)H, psi, 1, 0.0f, hPsi, 1);
	*E = cblas_ddot(this->N, Psi, 1, HPsi, 1);
}

template<int NValue>
void dStaticHamil<NValue>::HPsi( double* Psi, double* HPsi ) {
	cblas_dspmv(CblasColMajor, CblasUpper, this->N,
	            1.0f, (double*)this->H, Psi, 1, 0.0f, HPsi, 1);
}
template<int NValue>
void dStaticHamil<NValue>::HPsi( double* Psi, double* HPsi, int N ) {
	cblas_dsymm(CblasColMajor, CblasLeft, CblasUpper, this->N, N,
	            1.0f, (double*)this->H, dStaticHamil::N, Psi, dStaticHamil::N, 0.0f, HPsi, dStaticHamil::N);
}

template<int NValue>
double dStaticHamil<NValue>::Overlap( double* Bra, double* Ket ) {
	return cblas_ddot(this->N, Bra, 1, Ket, 1);
}
template<int NValue>
void dStaticHamil<NValue>::NormalizePsi( double* Psi, bool FlagNorm ) {
	if (FlagNorm) {
		auto norm = Psi[0] < 0.0f ? -1.0f : 1.0f / cblas_dnrm2(this->N, Psi, 1);
		cblas_dscal(this->N, norm, Psi, 1);
	} else if (Psi[0] < 0.0f)
		cblas_dscal(this->N, -1.0f, Psi, 1);
}

// Hamiltonian object
template< int NValue>
dStaticHamil<NValue>::StaticHamil::StaticHamil( double* H, bool copy, bool clearOriginal ) : dHamil<NValue>::Hamil(H) {
	if (!copy)
		return;
	this->H = (double*)malloc(this->N * this->N * sizeof(double));
//	cblas_dcopy(this->N * this->N, H, 1, this->H, 1);
	memcpy(this->H, H, this->N * this->N * sizeof(double));
	if (!clearOriginal)
		return;
	free(H);
}

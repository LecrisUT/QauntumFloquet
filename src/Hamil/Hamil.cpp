#include "Hamil.h"
#include <cstdlib>

using namespace QuanFloq;

// Base Hamiltonian object
template<int N>
dHamil<N>::dHamil( Hamil* H, double* Psi, double* E ) : H(H), Psi(Psi), E(E) { }
template<int N>
dHamil<N>::dHamil( Hamil* H ) : dHamil<N>(H, (double*)malloc(N * N * sizeof(double)),
                                          (double*)malloc((N * sizeof(double)))) { }
template<int N>
dHamil<N>::~dHamil() {
	free(Psi);
	free(E);
}

template<int N>
double dHamil<N>::PsiHPsi() {
	return PsiHPsi(Psi);
}
template<int N>
double dHamil<N>::PsiHPsi( double* Psi ) {
	double E;
	auto HPsi = (double*)malloc(N * sizeof(double));
	PsiHPsi(Psi, &E, HPsi);
	return E;
}

template<int N>
double* dHamil<N>::HPsi() {
	return HPsi(Psi);
}
template<int N>
double* dHamil<N>::HPsi( double* Psi ) {
	auto HPsi = (double*)malloc(N * sizeof(double));
	this->HPsi(Psi, HPsi);
	return HPsi;
}

// Base Hamiltonian data object
template<int N>
dHamil<N>::Hamil::Hamil() {
	H = (double*)malloc(N * N * sizeof(double));
}
template<int N>
dHamil<N>::Hamil::Hamil( double* H ) : H(H) { }
template<int N>
dHamil<N>::Hamil::~Hamil() {
	free(H);
}
template<int N>
dHamil<N>::Hamil::operator double*() const {
	return H;
}

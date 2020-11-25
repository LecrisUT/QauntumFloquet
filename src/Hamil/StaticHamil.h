#include "Hamil.h"

#ifndef QF_STATICHAMIL_H
#define QF_STATICHAMIL_H
namespace QuanFloq {
	/**
	 * Main Hamiltonian object
	 */
	template<int N>
	class dStaticHamil : dHamil<N> {
	public:
		class StaticHamil;

		//Fields
		const StaticHamil* const staticH;

		// Methods
		explicit dStaticHamil( dStaticHamil<N>::StaticHamil* H );
		explicit dStaticHamil( double* H );

		void PsiHPsi( double* Psi, double* E, double* HPsi ) override;
		void HPsi( double* Psi, double* HPsi ) override;
		void HPsi( double* Psi, double* HPsi, int N ) override;
		double Overlap( double* Bra, double* Ket ) override;
		void NormalizePsi( double* Psi, bool FlagNorm = false ) override;
	};

	/**
	 * Stored Hamiltonian object
	 */
	template<int N>
	class dStaticHamil<N>::StaticHamil : dHamil<N>::Hamil {
	public:
		explicit StaticHamil( double* H, bool copy = true, bool clearOriginal = false );
	};
}
#endif
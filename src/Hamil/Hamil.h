//
// Created by Le Minh Cristian on 2020/11/24.
//

#ifndef QF_HAMIL_H
#define QF_HAMIL_H

namespace QuanFloq {
	template<int NValue>
	class dHamil {
	public:
		// Nested class
		class Hamil;

		// Fields
		static const int N = NValue;
		const Hamil* const H;
		double* Psi;
		double* E;

		// Methods
	protected:
		dHamil( dHamil<NValue>::Hamil* H, double* Psi, double* E );
		explicit dHamil( dHamil<NValue>::Hamil* H );
		~dHamil();
	public:
		/** @defgroup PsiHPsi Energy expectation
		 * Expectation value of E
		 * @{
		 */
		/**
		 * @return Energy expectation
		 */
		virtual double PsiHPsi();
		/**
		 * @copydoc double PsiHPsi()
		 * @param Psi Wavefunction
		 */
		virtual double PsiHPsi( double* Psi );
		virtual void PsiHPsi( double* Psi, double* E, double* HPsi ) = 0;
		/** @} */
		virtual double* HPsi();
		virtual double* HPsi( double* Psi );
		virtual void HPsi( double* Psi, double* HPsi, int N ) = 0;
		virtual void HPsi( double* Psi, double* HPsi ) = 0;

		virtual double Overlap( double* Bra, double* Ket ) = 0;
		virtual void NormalizePsi( double* Psi, bool FlagNorm = false ) = 0;
	};

	template<int N>
	class dHamil<N>::Hamil {
	protected:
		double* H;
		Hamil();
		explicit Hamil( double* H );
		~Hamil();
	public:
		virtual explicit operator double*() const;
	};
}
#endif //QF_HAMIL_H

//
// Created by Le Minh Cristian on 2020/11/30.
//
#include <vHFHamil.h>
#include <vFloqHamil.h>

#ifndef QF_BASEFLOQHFHAMIL_H
#define QF_BASEFLOQHFHAMIL_H

namespace QuanFloq {

	// region Class Declarations
	template<typename T>
	class vFloqHFHamil;
	using svFloqHFHamil = vFloqHFHamil<float>;
	using dvFloqHFHamil = vFloqHFHamil<double>;
	using cvFloqHFHamil = vFloqHFHamil<cfloat >;
	using zvFloqHFHamil = vFloqHFHamil<cdouble >;
	// endregion

	template<typename T>
	class vFloqHFHamil :
			virtual public vHFHamil<T>,
			virtual public vFloqHamil<T> {
		// region Fields
	private:
		bool initialized = false;
	public:
		const int nFh;
	protected:
		const int csr_nUEx_max;
		int csr_nUEx;
		int* csr_ind_UEx;
		// endregion

		// region Methods
		// region Constructor/Destructor
	protected:
		vFloqHFHamil();
		explicit vFloqHFHamil( int nFh, bool initM = false );
		vFloqHFHamil( int nFh, int csr_nUEx, bool initM = false );
		vFloqHFHamil( int nFh, int csr_nUEx, int* csr_ind_UEx );
		// endregion

		// region Get/Set
	public:
		void seth( T* th ) override;
		void seth( T** th ) override;
		T* getUEx( T* tPsi = nullptr ) override;
		void getUEx( T* tUEx, T* tPsi = nullptr ) override;
		// endregion

	public:
		using vHamil<T>::HPsi;
		using vHamil<T>::PsiHPsi;
		using vHFHamil<T>::HPsi;
		using vHFHamil<T>::PsiHPsi;
		using vHFHamil<T>::CalcUEx;
		// Temporary due to MATLAB interface
		using vHFHamil<T>::UpdateH;
		using vHamil<T>::Overlap;
		using vHamil<T>::NormalizePsi;
		using vHFHamil<T>::geth;
		using vHFHamil<T>::getUEx;
		using vHamil<T>::getH;
		using vHamil<T>::getPsi;
		using vHamil<T>::getE;

		// region Main Methods
	public:
		int CalcCsr_nUEx_max();
		void Initialize();
		void Initialize( T* h, T w, T* UEx, double Tresh = 1E-10 );
		void UpdateH( T* tPsi ) override;
		void CalcSH_UEx();
	protected:
		virtual void CalcUEx( T* tPsi, T* acc, int iF );
		virtual void mCalcUEx( T* Bra, T* Ket, T* acc, int nFn );
		// endregion

	public:
		virtual const T* geth( int m, int n, int nf ) const;
		const T* geth(int m, int n) const override;
		void getUEx( T* tUEx, int m, int n ) override;
		void getUEx( T* tUEx, T* tPsi, int m, int n, int mPsi, int nPsi ) override;
		void UpdateH( T* tPsi, int m, int n ) override;
		const T* getH( int m, int n, int nf ) const override;
		const T* getH( int m, int n ) const override;
		const T* getPsi( int m, int n ) const override;
		const T* getE( int n ) const override;
		T Overlap( T* Bra, T* Ket, int n ) override;
		void NormalizePsi( T* tPsi, int n, bool FlagNorm = false ) override;
		void HPsi( T* tPsi, T* tHPsi, int m, int n ) override;
		void PsiHPsi( T* tPsi, T* tE, T* tHPsi, int m, int n ) override;
		// endregion
	};

	int CalcCsr_nUEx( int nUEx, int nFH, int n2F_max );
}

#endif //QF_BASEFLOQHFHAMIL_H

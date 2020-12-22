//
// Created by Le Minh Cristian on 2020/11/30.
//
#include "HFHamil.h"
#include "FloqHamil.h"

#ifndef QF_FLOQHFHAMIL_H
#define QF_FLOQHFHAMIL_H

namespace QuanFloq {
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
	public:
		vFloqHFHamil();
		explicit vFloqHFHamil( int nFh, bool initM = false );
		vFloqHFHamil( int nFh, int csr_nUEx, bool initM = false );
		vFloqHFHamil( int nFh, int csr_nUEx, int* csr_ind_UEx );
		// endregion

		// region Get/Set
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
		using vFloqHamil<T>::HPsi;
		using vFloqHamil<T>::PsiHPsi;

		// region Main Methods
	public:
		int CalcCsr_nUEx_max();

		void Initialize();
		void Initialize( T* h, T w, T* UEx, double Tresh = 1E-10 );
		void CalcSH_UEx();
		void UpdateH( T* tPsi ) override;
		virtual void CalcUEx( T* tPsi, T* acc, int iF );
	protected:
		virtual void mCalcUEx( T* Bra, T* Ket, T* acc, int nFn );
		// endregion
		// endregion
	};

	template<typename T>
	class FloqHFHamil final :
			vFloqHFHamil<T> {
	public:
		FloqHFHamil( int nH, Hamil_Sym SSym, int nFh, int nFH, int nF_max, int nElec, int nOrb, int nUEx );
		FloqHFHamil( int nH, Hamil_Sym SSym, int nFh, int nFH, int nF_max, int nElec, int nOrb, int nUEx, T* h, T w, T* UEx );
		FloqHFHamil( int nH, Hamil_Sym SSym, int nFh, int nFH, int nF_max, int nElec, int nOrb );
		FloqHFHamil( int nH, Hamil_Sym SSym, int nFh, int nFH, int nF_max, int nElec, int nOrb, T* h, T w, T* UEx );
	};

	int CalcCsr_nUEx( int nUEx, int nFH, int n2F_max );
}


#endif //QF_FLOQHFHAMIL_H

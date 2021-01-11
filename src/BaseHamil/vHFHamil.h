#include <vHamil.h>
//#include <type_traits>

#ifndef QF_BASEHFHAMIL_H
#define QF_BASEHFHAMIL_H

namespace QuanFloq {
	// region Class Declarations
	template<typename T>
	class vHFHamil;
	using svHFHamil = vHFHamil<float>;
	using dvHFHamil = vHFHamil<double>;
	using cvHFHamil = vHFHamil<cfloat >;
	using zvHFHamil = vHFHamil<cdouble >;
	// endregion

	template<typename T>
	class vHFHamil :
			virtual public vHamil<T> {
		// region Fields
	private:
		bool initialized = false;
	protected:
		const int nUEx_max;
	public:
		const int nElec;
		const int nOrb;
		int nUEx;

	protected:
		T* h;

		T* tensUEx;
		T* vh;
		T** indvH;
		// endregion

		// region Methods
		// region Constructor/Destructor
	protected:
		vHFHamil();
		vHFHamil( int nElec, int nOrb, bool initM = false );
		vHFHamil( int nElec, int nOrb, int nUEx, bool initM = false );
		vHFHamil( int nElec, int nOrb, int nUEx, T* h, T* tensUEx, T* vh, T** indvH );
		// endregion

		// region Get/Set
	public:
		[[nodiscard]] const T* geth() const;
		virtual void seth( T* th );
		virtual void seth( T** th );
		virtual T* getUEx( T* tPsi = nullptr );
		virtual void getUEx( T* tUEx, T* tPsi = nullptr );
		// endregion

		// region Inherited overloads
	public:
		using vHamil<T>::HPsi;
		using vHamil<T>::PsiHPsi;
		// Temporary due to MATLAB interface
		using vHamil<T>::Overlap;
		using vHamil<T>::NormalizePsi;
		using vHamil<T>::getH;
		using vHamil<T>::getPsi;
		using vHamil<T>::getE;
		// endregion

		// region Main Methods
	public:
		void Initialize( T* h, T* UEx, double Tresh = 1E-10 );
//		bool CalcTensor_p1( T val );

		virtual void HPsi( T* tPsi, T* tHPsi, bool tupdateH );
		virtual void PsiHPsi( T* tPsi, T* tE, T* tHPsi, bool tupdateH );

		virtual void CalcUEx( T* tPsi, T* acc );
		virtual void UpdateH( T* tPsi );
		static void UpdateH( vHamil <T>* Hamil, T* tPsi );

	protected:
		virtual void mCalcUEx( T* Bra, T* Ket, T* acc );
		virtual void mCalcUEx_p1( T* BraKet, T* acc );
		// endregion

	public:
		virtual const T* geth(int m, int n) const;
		virtual void getUEx( T* tUEx, int m, int n );
		virtual void getUEx( T* tUEx, T* tPsi, int m, int n, int mPsi, int nPsi );
		virtual void UpdateH( T* tPsi, int m, int n );
		const T* getH( int m, int n ) const override;
		const T* getPsi( int m, int n ) const override;
		const T* getE( int n ) const override;
		T Overlap( T* Bra, T* Ket, int n ) override;
		void NormalizePsi( T* tPsi, int n, bool FlagNorm = false ) override;
		void HPsi( T* tPsi, T* tHPsi, int m, int n ) override;
		void PsiHPsi( T* tPsi, T* tE, T* tHPsi, int m, int n ) override;
		// endregion
	};

#if __cplusplus >= 202002L
	struct HFHamilSize :
			public HamilSize {
		const int nUEx;
		constexpr explicit HFHamilSize( int tnH, int tnUEx = -1 ) : HamilSize(tnH),
																	nUEx(tnUEx < 0 ? tnH * tnH : tnUEx) { }
	} __attribute__((aligned(4))) __attribute__((packed));

	template<typename T, HFHamilSize Sz>
	class [[maybe_unused]] tHFHamil :
			virtual public vHFHamil<T>,
			virtual public tHamil<T, Sz> {
	public:
		T h[Sz.nH * Sz.nH];
		T tensUEx[Sz.nUEx * Sz.nH * Sz.nH];
		T vh[Sz.nUEx];
		T* indvH[Sz.nUEx];

	public:
		tHFHamil();
		explicit tHFHamil( T* th );
	};
#endif
}

#endif //QF_HFHAMIL_H

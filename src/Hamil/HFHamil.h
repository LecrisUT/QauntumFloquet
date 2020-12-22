#include "Hamil.h"
#include <type_traits>

#ifndef QF_HFHAMIL_H
#define QF_HFHAMIL_H

namespace QuanFloq {
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
	public:
		vHFHamil();
		vHFHamil( int nElec, int nOrb, bool initM = false );
		vHFHamil( int nElec, int nOrb, int nUEx, bool initM = false );
		vHFHamil( int nElec, int nOrb, int nUEx, T* h, T* tensUEx, T* vh, T** indvH );
		// endregion

		// region Get/Set
		[[nodiscard]] T* geth() const;
		virtual void seth( T* th );
		virtual void seth( T** th );
		virtual T* getUEx( T* tPsi = nullptr );
		virtual void getUEx( T* tUEx, T* tPsi = nullptr );
		// endregion

		// region Inherited overloads
	public:
		using vHamil<T>::HPsi;
		using vHamil<T>::PsiHPsi;
		// endregion

		// region Main Methods
	public:
		void Initialize( T* h, T* UEx, double Tresh = 1E-10 );
//		bool CalcTensor_p1( T val );

		virtual void HPsi( T* tPsi, T* tHPsi, bool tupdateH );
		virtual void PsiHPsi( T* tPsi, T* tE, T* tHPsi, bool tupdateH );

		virtual void CalcUEx( T* tPsi, T* acc );
		virtual void UpdateH( T* tPsi );
		static void UpdateH( vHamil<T>* Hamil, T* tPsi );

	protected:
		virtual void mCalcUEx( T* Bra, T* Ket, T* acc );
		virtual void mCalcUEx_p1( T* BraKet, T* acc );
		// endregion
		// endregion
	};

	template<typename T>
	class HFHamil final:
			virtual public vHFHamil<T> {
	public:
		HFHamil( int nH, Hamil_Sym Sym, int nElec, int nOrb, int nUEx );
		HFHamil( int nH, Hamil_Sym Sym, int nElec, int nOrb, int nUEx, T* h, T* UEx );
		HFHamil( int nH, Hamil_Sym Sym, int nElec, int nOrb );
		HFHamil( int nH, Hamil_Sym Sym, int nElec, int nOrb, T* h, T* UEx );
	};

	struct HFHamilSize :
			public HamilSize {
		const int nUEx;
		constexpr explicit HFHamilSize( int tnH, int tnUEx = -1 ) : HamilSize(tnH),
		                                                            nUEx(tnUEx < 0 ? tnH * tnH : tnUEx) { }
	} __attribute__((aligned(4))) __attribute__((packed));
//	constexpr SizeStruct test = SizeStruct(0);

//	template<typename T, HFHamilSize Sz>
	template<typename T, HFHamilSize const& Sz>
	class [[maybe_unused]] tHFHamil :
			virtual public vHFHamil<T>,
			virtual public tHamil<T, (HamilSize&)Sz> {
	public:
		T h[Sz.nH * Sz.nH];
		T tensUEx[Sz.nUEx * Sz.nH * Sz.nH];
		T vh[Sz.nUEx];
		T* indvH[Sz.nUEx];

	public:
		tHFHamil();
		explicit tHFHamil( T* th );
	};
}

#endif //QF_HFHAMIL_H

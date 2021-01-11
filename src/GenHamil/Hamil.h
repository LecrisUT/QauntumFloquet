//
// Created by Le Minh Cristian on 2020/11/24.
//
#include <complex>
#include <list>
#include "../Conf.h"

#ifndef QF_HAMIL_H
#define QF_HAMIL_H

// TODO: Add final classes to devirtualize and inline function calls
// TODO: Make compatible with inline
namespace QuanFloq {
	// region Const values
	static constexpr cfloat cfloat1 = cfloat(1.0f);
	static constexpr cfloat cfloat0 = cfloat(0.0f);
	static constexpr cdouble cdouble1 = cdouble(1.0f);
	static constexpr cdouble cdouble0 = cdouble(0.0f);
	static constexpr float float1 = float(1.0f);
	static constexpr float float0 = float(0.0f);
	static constexpr double double1 = double(1.0f);
	static constexpr double double0 = double(0.0f);
	// endregion

	// region Class Declarations
//	template<typename T>
//	class IHamil;
	template<typename T>
	class vHamil;
	template<typename T>
	class Hamil;
#if __cplusplus >= 202002L
	struct HamilSize;
	template<typename T, HamilSize Sz>
	class [[maybe_unused]] tHamil;
#endif
	using sHamil = vHamil<float>;
	using dHamil = vHamil<double>;
	using cHamil = vHamil<cfloat >;
	using zHamil = vHamil<cdouble >;
	// endregion

	// region Class Definitions
	enum Hamil_Sym {
		Full,
		Sym,
		Her,
	};
	/**
	 * Basic Hamiltonian object
	 */
	template<typename T>
//	class vHamil : virtual public IHamil<T> {
	class vHamil {
		// region Fields
	private:
		bool initialized = false;
	public:
		const int nH;
		const Hamil_Sym Sym;
	protected:
		T* H;
		T* Psi;
		T* E;
	public:
		std::list<void (*)( vHamil<T>*, T* )> PreHPsi;
		std::list<void (*)( vHamil<T>*, T*, T* )> PostHPsi;
		std::list<void (*)( vHamil<T>*, T* )> PrePsiHPsi;
		std::list<void (*)( vHamil<T>*, T*, T*, T* )> PostPsiHPsi;
		// endregion

		// region Methods
		// region Constructor/Destructor
	public:
		vHamil();
		explicit vHamil( int nH, Hamil_Sym Sym, bool initM = false );
		vHamil( int nH, Hamil_Sym Sym, T* H, T* Psi, T* E );
		// endregion

		// region Get/Set
		[[nodiscard]] T* getH() const;
		virtual void setH( T* H );
		// TODO: add setH (T(*)[Sz.nH]) for templated Sz, include old interface as well
		virtual void setH( T** H );
		[[nodiscard]] const T* getPsi() const;
		[[nodiscard]] const T* getE() const;
		// endregion

		// region Main methods
	public:
		void Initialize( T* tH );
	protected:
		virtual void mHPsi( T* tPsi, T* tHPsi );
		virtual void mPsiHPsi( T* tPsi, T* tE, T* tHPsi );

	public:
		virtual T* HPsi();
		virtual T* HPsi( T* tPsi );
		virtual void HPsi( T* tPsi, T* tHPsi );
		virtual T PsiHPsi();
		virtual T PsiHPsi( T* tPsi );
		virtual void PsiHPsi( T* tPsi, T* tE, T* tHPsi );

		virtual T Overlap( T* Bra, T* Ket );
		virtual void NormalizePsi( T* tPsi, bool FlagNorm = false );
		// endregion
		// endregion
	};

	// TODO: Implement final class and rename virtual class
	template<typename T>
	class Hamil :
			public vHamil<T> {
	public:
		explicit Hamil( int nH, Hamil_Sym Sym );
		Hamil( int nH, Hamil_Sym Sym, T* H );
	};

#if __cplusplus >= 202002L
	struct HamilSize {
		const int nH;
		const int nH2 = nH * nH;
		const int nH2t = nH * (nH + 1) / 2;
		constexpr explicit HamilSize( int tnH ) : nH(tnH) { };
	} __attribute__((aligned(16)));

	template<typename T, HamilSize Sz>
	class tHamil :
			public vHamil<T> {
	public:
		T H[Sz.nH * (Sz.nH + 1) / 2];
		T Psi[Sz.nH * Sz.nH];
		T E[Sz.nH];

	public:
		tHamil();
		explicit tHamil( T* tH );
	};
#endif
	// endregion

	// region template initialization
	extern template
	class Hamil<float>;
	extern template
	class Hamil<double>;
	extern template
	class Hamil<cfloat>;
	extern template
	class Hamil<cdouble>;
	// endregion
}

//template
//class QuanFloq::IHamil<double>;
//template
//class QuanFloq::IHamil<float>;
//template
//class QuanFloq::IHamil<cfloat >;
//template
//class QuanFloq::IHamil<cdouble >;
#endif //QF_HAMIL_H

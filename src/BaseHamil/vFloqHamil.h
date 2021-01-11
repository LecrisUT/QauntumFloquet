#include <vHamil.h>
#include <mkl.h>

#ifndef QF_BASEFLOQHAMIL_H
#define QF_BASEFLOQHAMIL_H
namespace QuanFloq {
	// region Const values
#ifdef OLDMKL
	static constexpr char charU = 'U';
	static constexpr char charN = 'N';
#endif
	// endregion

	// region Class Declarations
	template<typename T>
	class vFloqHamil;
	using svFloqHamil = vFloqHamil<float>;
	using dvFloqHamil = vFloqHamil<double>;
	using cvFloqHamil = vFloqHamil<cfloat >;
	using zvFloqHamil = vFloqHamil<cdouble >;
	// endregion
/**
 * Main Floquet vHamil object
 */
	template<typename T>
	class vFloqHamil :
			virtual public vHamil<T> {
		// region Fields
	private:
		bool initialized = false;
	public:
		const int nFH;
		const int nF_max;
	protected:
		const int n2F_max;
		const int nH2F_max;
		T w;
		bool CalcHf;
	public:
		struct matrix_descr SDescr;
		sparse_matrix_t SH;
		sparse_matrix_t SHf;
		sparse_matrix* ActiveSH;
	protected:
		const int csr_nvalues;
		T* csr_values_SH;
		T* csr_values_SHf;
		T* csr_values_Active;
		int* csr_columns;
		int* csr_rows;
		int* csr_diagInd;
		T** csr_values_map;
		int* diag_ind;
		T* partialt;
		// endregion

		// region Methods
		// region Constructor/Destructor
	protected:
		vFloqHamil();
		vFloqHamil( int tnFH, int tnF_max, bool initM = false );
		vFloqHamil( int nFh, int nFMax, int csrNvalues, int* csrColumns, int* csrRows, int* csrDiagInd,
		            T* csrValuesSh, T* csrValuesSHf, T** csrValuesMap, int* diagInd, T* partialt );
		// endregion

		// region Get/Set
	public:
		T getW() const;
		virtual void setW( T tw, bool CalcS = true );
		void setH( T* tH, bool CalcS );
		void setH( T* tH ) override;
		[[nodiscard]] bool isCalcHf() const;
		void setCalcHf( bool calcHf );
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
		void Initialize();
		void Initialize( T* tH, T tw );
	protected:
		void mHPsi( T* tPsi, T* tHPsi ) override;
#ifndef OLDMKL
		void mPsiHPsi( T* tPsi, T* tE, T* tHPsi ) override;
#endif
	public:
		T* HPsi( T* tPsi ) override;
		T PsiHPsi( T* tPsi ) override;
		T Overlap( T* Bra, T* Ket ) override;
		void NormalizePsi( T* tPsi, bool FlagNorm ) override;

		virtual void CalcSH();
		virtual void CalcSHf();
	protected:
		virtual void CalcSH_part();
#ifndef OLDMKL
		virtual void CalcSHf_part();
#endif
		// endregion

	public:
		virtual const T* getH( int m, int n, int nf ) const;
		const T* getH( int m, int n ) const override;
		const T* getPsi( int m, int n ) const override;
		const T* getE( int n ) const override;
		T Overlap( T* Bra, T* Ket, [[maybe_unused]] [[maybe_unused]] int n ) override;
		void NormalizePsi( T* tPsi, int n, bool FlagNorm = false ) override;
		void HPsi( T* tPsi, T* tHPsi, int m, int n ) override;
		void PsiHPsi( T* tPsi, T* tE, T* tHPsi, int m, int n ) override;
		// endregion
	};

#if __cplusplus >= 202002L
	struct FloqHamilSize :
			public HamilSize {
		const int nF_max = 0, nFH = 0;
		const int n2F_max;
		const int nH2F_max;
		constexpr explicit FloqHamilSize( int tnH ) : HamilSize(tnH), n2F_max(2 * nF_max + 1),
		                                              nH2F_max(nH * n2F_max) { }
	} __attribute__((aligned(16))) __attribute__((packed));

	template<typename T, FloqHamilSize Sz>
	class [[maybe_unused]] tFloqHamil :
			public vFloqHamil<T> {
	public:
		T H[Sz.nH * Sz.nH * (Sz.nFH + 1)];
		T Psi[Sz.nH2F_max * Sz.nH];
		T E[Sz.nH];
		T partialt[Sz.nH2F_max];
		int diag_ind[Sz.nH2F_max];

	public:
		tFloqHamil();
		explicit tFloqHamil( T* tH );
	};
#endif
	// endregion
}
#endif

#include "Hamil.h"
#include <mkl.h>

#ifndef QF_FLOQHAMIL_H
#define QF_FLOQHAMIL_H
namespace QuanFloq {
	// region Const values
#ifdef OLDMKL
	static constexpr char charU = 'U';
	static constexpr char charN = 'N';
#endif
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
	protected:
		T w;
		bool CalcHf;
	public:
		struct matrix_descr SDescr;
		const int nFH;
		const int nF_max;
		const int n2F_max;
		const int nH2F_max;
		sparse_matrix_t SH;
		sparse_matrix_t SHf;
		sparse_matrix* ActiveSH;
		const int csr_nvalues;
	protected:
		T* csr_values_SH;
		T* csr_values_SHf;
		T* csr_values_Active;
		int* csr_columns;
		int* csr_rows;
		int* csr_diagInd;
		T** csr_values_map;
	protected:
		int* diag_ind;
		T* partialt;
		// endregion

		// region Methods
		// region Constructor/Destructor
	public:
		vFloqHamil();
		vFloqHamil( int tnFH, int tnF_max, bool initM = false );
		vFloqHamil( int nFh, int nFMax, int csrNvalues, int* csrColumns, int* csrRows, int* csrDiagInd,
		            T* csrValuesSh, T* csrValuesSHf, T** csrValuesMap, int* diagInd, T* partialt );
	public:
		// endregion
		// region Get/Set
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
		// endregion
	};

	template<typename T>
	class FloqHamil final :
			virtual public vFloqHamil<T> {
	public:
		FloqHamil( int nH, Hamil_Sym SSym, int nFH, int nF_max );
		FloqHamil( int nH, Hamil_Sym SSym, int nFH, int nF_max, T* H, T w );
	};

	struct FloqHamilSize :
			public HamilSize {
		const int nF_max = 0, nFH = 0;
		const int n2F_max;
		const int nH2F_max;
		constexpr explicit FloqHamilSize( int tnH ) : HamilSize(tnH), n2F_max(2 * nF_max + 1),
		                                              nH2F_max(nH * n2F_max) { }
	} __attribute__((aligned(16))) __attribute__((packed));
//	constexpr SizeStruct test = SizeStruct(0);

//	template<typename T, SizeStruct Sz>
	template<typename T, FloqHamilSize const& Sz>
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
	// endregion
}
#endif

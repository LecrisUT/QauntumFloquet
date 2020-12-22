//
// Created by Le Minh Cristian on 2020/12/07.
//
#include "../Hamil/HFHamil.h"
#include "Dimer.h"

#ifndef QF_HFDIMER_H
#define QF_HFDIMER_H

namespace QuanFloq {
	// region Class Declarations
	template<typename T>
	class HFDimer;
	using sHFDimer = HFDimer<float>;
	using dHFDimer = HFDimer<double>;
	using cHFDimer = HFDimer<cfloat >;
	using zHFDimer = HFDimer<cdouble >;
	// endregion

	template<typename T>
	class HFDimer :
			virtual public vHamil<T>,
			virtual public vHFHamil<T>,
			virtual public Dimer<T> {

	public:
		HFDimer();
		explicit HFDimer( T v, T U, T t = 1.0f );
	protected:
		T* MakeUEx( T U );

	public:
		void setT( T t ) override;
		void setV( T v ) override;
		void setU( T U ) override;
	};
}

#endif //QF_HFDIMER_H

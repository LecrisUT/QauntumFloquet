//
// Created by Le Minh Cristian on 2020/12/07.
//
#include "../Hamil/Hamil.h"

#ifndef QF_DIMER_H
#define QF_DIMER_H

namespace QuanFloq {
	// region Class Declarations
	template<typename T>
	class Dimer;
	using sDimer = Dimer<float>;
	using dDimer = Dimer<double>;
	using cDimer = Dimer<cfloat >;
	using zDimer = Dimer<cdouble >;
	// endregion

	template<typename T>
	class Dimer :
			virtual public vHamil<T> {
	protected:
		T t;
		T v;
		T U;

	public:
		Dimer();
		Dimer( T v, T U, T t = 1.0f );

		[[nodiscard]] T getT() const;
		virtual void setT( T t );
		[[nodiscard]] T getV() const;
		virtual void setV( T v );
		[[nodiscard]] T getU() const;
		virtual void setU( T U );
	};
}

#endif //QF_DIMER_H

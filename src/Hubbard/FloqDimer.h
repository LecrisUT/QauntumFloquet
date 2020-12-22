//
// Created by Le Minh Cristian on 2020/12/07.
//
#include "../Hamil/FloqHamil.h"
#include "Dimer.h"

#ifndef QF_FLOQDIMER_H
#define QF_FLOQDIMER_H

namespace QuanFloq {
	// region Class Declarations
	template<typename T>
	class FloqDimer;
	using sFloqDimer = FloqDimer<float>;
	using dFloqDimer = FloqDimer<double>;
	using cFloqDimer = FloqDimer<cfloat >;
	using zFloqDimer = FloqDimer<cdouble >;
	// endregion

	template<typename T>
	class FloqDimer :
			virtual public vHamil<T>,
			virtual public vFloqHamil<T>,
			virtual public Dimer<T> {
	protected:
		T v1;

	public:
		FloqDimer();
		FloqDimer( int nF_max, T v0, T v1, T U, T w, T t = 1.0f);


		virtual void setT( T t, bool CalcS = true );
		[[nodiscard]] T getV0() const;
		virtual void setV0( T v0, bool CalcS = true );
		[[nodiscard]] T getV1() const;
		virtual void setV1( T v1, bool CalcS = true );
		virtual void setU( T U, bool CalcS = true );
		void setT( T t ) override;
		void setV( T v0 ) override;
		void setU( T U ) override;
	};
}

#endif //QF_FLOQDIMER_H

//
// Created by Le Minh Cristian on 2020/12/07.
//
#include "../Hamil/FloqHFHamil.h"
#include "HFDimer.h"
#include "FloqDimer.h"

#ifndef QF_FLOQHFDIMER_H
#define QF_FLOQHFDIMER_H

namespace QuanFloq {
	// region Class Declarations
	template<typename T>
	class FloqHFDimer;
	using sFloqHFDimer = FloqHFDimer<float>;
	using dFloqHFDimer = FloqHFDimer<double>;
	using cFloqHFDimer = FloqHFDimer<cfloat >;
	using zFloqHFDimer = FloqHFDimer<cdouble >;
	// endregion

	template<typename T>
	class FloqHFDimer :
			virtual public vFloqHFHamil<T>,
			virtual public HFDimer<T>,
			virtual public FloqDimer<T>   {

	public:
		FloqHFDimer();
		FloqHFDimer( int nFH, int nF_max, T v0, T v1, T U, T w, T t = 1.0f );

	public:
		void setT( T t, bool CalcS ) override;
		void setV0( T v0, bool CalcS ) override;
		void setV1( T v1, bool CalcS ) override;
		void setU( T U, bool CalcS ) override;

		void setT( T t ) override;
		void setV( T v0 ) override;
		void setU( T U ) override;
	};
}


#endif //QF_FLOQHFDIMER_H

//
// Created by Le Minh Cristian on 2020/12/07.
//

#ifndef QF_BASEHFDIMER_H
#define QF_BASEHFDIMER_H

#include <vHFHamil.h>
#include <vDimer.h>

namespace QuanFloq {
	// region Class Declarations
	template<typename T>
	class vHFDimer;
	using svHFDimer = vHFDimer<float>;
	using dvHFDimer = vHFDimer<double>;
	using cvHFDimer = vHFDimer<cfloat >;
	using zvHFDimer = vHFDimer<cdouble >;
	// endregion

	template<typename T>
	class vHFDimer :
			virtual public vHFHamil<T>,
			virtual public vDimer<T> {

	protected:
		vHFDimer();

	protected:
		T* MakeUEx( T U );

	public:
		void setT( T t ) override;
		void setV( T v ) override;
		void setU( T U ) override;
	};
}

#endif //QF_BASEHFDIMER_H

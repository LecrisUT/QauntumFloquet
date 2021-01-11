//
// Created by Le Minh Cristian on 2020/12/07.
//

#ifndef QF_BASEHFDIMER_H
#define QF_BASEHFDIMER_H

#include "../BaseHamil/vHFHamil.h"
#include "vDimer.h"

namespace QuanFloq {
	// region Class Declarations
	template<typename T>
	class vHFDimer;
	// endregion

	template<typename T>
	class vHFDimer :
			virtual public vHamil<T>,
			virtual public vHFHamil<T>,
			virtual public vDimer<T> {

	public:
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

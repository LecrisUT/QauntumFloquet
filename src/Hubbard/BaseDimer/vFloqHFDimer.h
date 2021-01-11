//
// Created by Le Minh Cristian on 2020/12/24.
//

#ifndef QF_BASEFLOQHFDIMER_H
#define QF_BASEFLOQHFDIMER_H

#include "../BaseHamil/vFloqHFHamil.h"
#include "vHFDimer.h"
#include "vFloqDimer.h"

namespace QuanFloq {
	// region Class Declarations
	template<typename T>
	class vFloqHFDimer;
	template<typename T>
	class vFloqHFDimer :
			virtual public vFloqHFHamil<T>,
			virtual public vHFDimer<T>,
			virtual public vFloqDimer<T> {

	public:
		vFloqHFDimer();

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

#endif //QF_BASEFLOQHFDIMER_H
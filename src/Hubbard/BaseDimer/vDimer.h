//
// Created by Le Minh Cristian on 2020/12/24.
//

#ifndef QF_BASEDIMER_H
#define QF_BASEDIMER_H

#include "../BaseHamil/vHamil.h"
// region Class Declarations
namespace QuanFloq{
	template<typename T>
	class vDimer;
	template<typename T>
	class vDimer :
			virtual public vHamil<T> {
	protected:
		T t;
		T v;
		T U;

	public:
		vDimer();

		[[nodiscard]] T getT() const;
		virtual void setT( T t );
		[[nodiscard]] T getV() const;
		virtual void setV( T v );
		[[nodiscard]] T getU() const;
		virtual void setU( T U );
	};
}

#endif //QF_BASEDIMER_H

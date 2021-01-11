//
// Created by Le Minh Cristian on 2020/12/24.
//

#ifndef QF_BASEDIMER_H
#define QF_BASEDIMER_H

#include <vHamil.h>
namespace QuanFloq{
	// region Class Declarations
	template<typename T>
	class vDimer;
	using svDimer = vDimer<float>;
	using dvDimer = vDimer<double>;
	using cvDimer = vDimer<cfloat >;
	using zvDimer = vDimer<cdouble >;
	// endregion

	template<typename T>
	class vDimer :
			virtual public vHamil<T> {
	protected:
		T t;
		T v;
		T U;

	protected:
		vDimer();

	public:
		[[nodiscard]] T getT() const;
		virtual void setT( T t );
		[[nodiscard]] T getV() const;
		virtual void setV( T v );
		[[nodiscard]] T getU() const;
		virtual void setU( T U );
	};
}

#endif //QF_BASEDIMER_H

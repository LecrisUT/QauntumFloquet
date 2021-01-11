//
// Created by Le Minh Cristian on 2020/12/24.
//

#ifndef QF_BASEFLOQDIMER_H
#define QF_BASEFLOQDIMER_H

#include <vFloqHamil.h>
#include <vDimer.h>

namespace QuanFloq {
	// region Class Declaration
	template<typename T>
	class vFloqDimer;
	using svFloqDimer = vFloqDimer<float>;
	using dvFloqDimer = vFloqDimer<double>;
	using cvFloqDimer = vFloqDimer<cfloat >;
	using zvFloqDimer = vFloqDimer<cdouble >;
	// endregion

	template<typename T>
	class vFloqDimer :
			virtual public vFloqHamil<T>,
			virtual public vDimer<T> {
	protected:
		T v1;

	protected:
		vFloqDimer();

	public:
		virtual void setT( T t, bool CalcS );
		[[nodiscard]] T getV0() const;
		virtual void setV0( T v0, bool CalcS );
		[[nodiscard]] T getV1() const;
		virtual void setV1( T v1, bool CalcS );
		virtual void setU( T U, bool CalcS );
		void setT( T t ) override;
		void setV( T v0 ) override;
		void setU( T U ) override;
	};
}
#endif //QF_BASEFLOQDIMER_H
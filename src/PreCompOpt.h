#include "Hamil/PreCompOpt_Hamil.h"

#ifndef QF_PRECOMPOPT_H
#define QF_PRECOMPOPT_H

#include <cassert>

#ifndef cfloat
#if Use_mkl
#define cfloat MKL_Complex8
#else
#define cfloat std::complex<float>
#define MKL_Complex8 cfloat
#endif
#else
#define MKL_Complex8 cfloat
#endif

#ifndef cdouble
#if Use_mkl
#define cdouble MKL_Complex16
#else
#define cdouble std::complex<double>
#define MKL_Complex16 cdouble
#endif
#else
#define MKL_Complex16 cdouble
#endif

#endif //

//
// Created by Le Minh Cristian on 2020/12/07.
//
#include <vHamil.h>
#include <vFloqHamil.h>
#include <vHFHamil.h>
#include <vFloqHFHamil.h>

#ifndef QF_TEST_H
#define QF_TEST_H

namespace QuanFloq {
	template<typename T>
	class Test {
	public:
		static void PrintFloqMatrix (T* Mtr, int n, int nF, bool flag_Full = false);
		static void PrintTriangMatrix (T* Mtr, int n);
		static void PrintMatrix( T* M, int m, int n, bool Transpose = false );
		static void PrintMatrix( T* M, int m, int n, int mSep, int nSep, bool Transpose = false );
		static void PrintH( vHamil<T>* H );
		static void PrintH( vFloqHamil<T>* H, bool flag_Full = false );
		static void Printh( vHFHamil<T>* H );
		static void Printh( vFloqHFHamil<T>* H, bool flag_Full = false );
		static void PrintUEx( vHFHamil<T>* H );
		static void PrintUEx( vFloqHFHamil<T>* H, bool flag_Full = false );
		static void PrintSH( vFloqHamil<T>* H );
		static void PrintSHf( vFloqHamil<T>* H );
#ifdef OLDMKL
		static void ConvertCsrToDense (vFloqHamil<T>* H, T* values, T* M);
#else
		static void ConvertSparseToDense (sparse_matrix* S, int N, matrix_descr Descr, T* M);
#endif
	};
}

#endif //QF_TEST_H

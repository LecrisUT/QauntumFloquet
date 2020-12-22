//
// Created by Le Minh Cristian on 2020/12/07.
//
#include <cmath>
#include <iostream>
#include "../src/Hubbard/Dimer.h"
#include "../src/Hubbard/FloqDimer.h"
#include "../src/Hubbard/HFDimer.h"
#include "../src/Hubbard/FloqHFDimer.h"
#include "Test.h"

using namespace QuanFloq;

template<typename T>
void PerformTest();
template<typename T>
void PerformFloqTest();
template<typename T>
void PerformHFTest();
template<typename T>
void PerformFloqHFTest();

int main() {
	PerformTest<float>();
//	PerformTest<double>();
//	PerformTest<cfloat >();
//	PerformTest<cdouble >();
	PerformFloqTest<float>();
//	PerformFloqTest<double>();
//	PerformFloqTest<cfloat >();
//	PerformFloqTest<cdouble >();
	PerformHFTest<float>();
//	PerformHFTest<double>();
//	PerformHFTest<cfloat >();
//	PerformHFTest<cdouble >();
	PerformFloqHFTest<float>();
//	PerformFloqHFTest<double>();
//	PerformFloqHFTest<cfloat >();
//	PerformFloqHFTest<cdouble >();
	return 0;
}

template<typename T>
void PerformTest() {
	T Psi[] = {1.0f, 2.0f, 3.0f};
	T HPsi[3] = {};
	T E;
	auto H = Dimer<T>(1.0f, 2.0f);
	std::cout << "Testing for Dimer<" << typeid(T).name() << ">:" << std::endl;
	std::cout << "PreHPsi:" << H.PreHPsi.size();
	std::cout << " PostHPsi:" << H.PostHPsi.size();
	std::cout << " PrePsiHPsi:" << H.PrePsiHPsi.size();
	std::cout << " PostPsiHPsi:" << H.PostPsiHPsi.size() << std::endl;
	std::cout << "H" << std::endl;
	Test<T>::PrintH(&H);
	std::cout << "Psi" << std::endl;
	Test<T>::PrintMatrix(Psi, 3, 1);
	std::cout << "Init HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 3, 1);
	H.PsiHPsi(Psi, &E, HPsi);
	std::cout << "PsiHPsi: E=" << E << " HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 3, 1);
	H.HPsi(Psi, HPsi);
	std::cout << "HPsi: HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 3, 1);
}
template<typename T>
void PerformFloqTest() {
	T Psi[21] = {};
	Psi[9] = 1.0f;
	Psi[10] = 2.0f;
	Psi[11] = 3.0f;
	T HPsi[21] = {};
	T E;
	auto H = FloqDimer<T>(3, 2.0f, 3.0f, 5.0f, 1.5f);
	std::cout << "Testing for FloqDimer<" << typeid(T).name() << ">:" << std::endl;
	std::cout << "PreHPsi:" << H.PreHPsi.size();
	std::cout << " PostHPsi:" << H.PostHPsi.size();
	std::cout << " PrePsiHPsi:" << H.PrePsiHPsi.size();
	std::cout << " PostPsiHPsi:" << H.PostPsiHPsi.size() << std::endl;
	std::cout << "H" << std::endl;
	Test<T>::PrintH(&H);
	std::cout << "SH" << std::endl;
	Test<T>::PrintSH(&H);
	std::cout << "SHf" << std::endl;
	Test<T>::PrintSHf(&H);
	std::cout << "Psi" << std::endl;
	Test<T>::PrintMatrix(Psi, 3, 7, true);
	std::cout << "Init HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 3, 7, true);
	H.PsiHPsi(Psi, &E, HPsi);
	std::cout << "PsiHPsi: E=" << E << " HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 3, 7, true);
	H.HPsi(Psi, HPsi);
	std::cout << "HPsi: HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 3, 7, true);
}
template<typename T>
void PerformHFTest() {
	T Psi[] = {0.5f, std::sqrt(3.0f) * 0.5f};
	T HPsi[2] = {};
	T E;
	auto H = HFDimer<T>(1.0f, 2.0f, 1.0f);
	std::cout << "Testing for HFDimer<" << typeid(T).name() << ">:" << std::endl;
	std::cout << "PreHPsi:" << H.PreHPsi.size();
	std::cout << " PostHPsi:" << H.PostHPsi.size();
	std::cout << " PrePsiHPsi:" << H.PrePsiHPsi.size();
	std::cout << " PostPsiHPsi:" << H.PostPsiHPsi.size() << std::endl;
	std::cout << "H" << std::endl;
	Test<T>::PrintH(&H);
	std::cout << "h" << std::endl;
	Test<T>::Printh(&H);
	std::cout << "UEx" << std::endl;
	Test<T>::PrintUEx(&H);
	H.UpdateH(Psi);
	std::cout << "H" << std::endl;
	Test<T>::PrintH(&H);
	std::cout << "UEx" << std::endl;
	Test<T>::PrintUEx(&H);
	std::cout << "Psi" << std::endl;
	Test<T>::PrintMatrix(Psi, 2, 1);
	std::cout << "Init HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 2, 1);
	H.PsiHPsi(Psi, &E, HPsi, false);
	std::cout << "PsiHPsi: E=" << E << " HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 2, 1);
	H.PsiHPsi(Psi, &E, HPsi);
	std::cout << "PsiHPsi: E=" << E << " HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 2, 1);
	H.HPsi(Psi, HPsi, false);
	std::cout << "HPsi: HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 2, 1);
}
template<typename T>
void PerformFloqHFTest() {
	T Psi[14] = {};
	Psi[8] = 0.5f;
	Psi[9] = std::sqrt(3.0f) * 0.5f;
	T HPsi[14] = {};
	T E;
	auto H = FloqHFDimer<T>(2, 3, 2.0f, 3.0f, 5.0f, 1.5f);
	std::cout << "Testing for HFDimer<" << typeid(T).name() << ">:" << std::endl;
	std::cout << "PreHPsi:" << H.PreHPsi.size();
	std::cout << " PostHPsi:" << H.PostHPsi.size();
	std::cout << " PrePsiHPsi:" << H.PrePsiHPsi.size();
	std::cout << " PostPsiHPsi:" << H.PostPsiHPsi.size() << std::endl;
	std::cout << "H" << std::endl;
	Test<T>::PrintH(&H);
	std::cout << "h" << std::endl;
	Test<T>::Printh(&H);
	std::cout << "UEx" << std::endl;
	Test<T>::PrintUEx(&H);
	std::cout << "SH" << std::endl;
	Test<T>::PrintSH(&H);
	std::cout << "SHf" << std::endl;
	Test<T>::PrintSHf(&H);
	H.UpdateH(Psi);
	std::cout << "H" << std::endl;
	Test<T>::PrintH(&H);
	std::cout << "UEx" << std::endl;
	Test<T>::PrintUEx(&H);
	std::cout << "SH" << std::endl;
	Test<T>::PrintSH(&H);
	std::cout << "SHf" << std::endl;
	Test<T>::PrintSHf(&H);
	std::cout << "Psi" << std::endl;
	Test<T>::PrintMatrix(Psi, 2, 7, true);
	std::cout << "Init HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 2, 7, true);
	H.PsiHPsi(Psi, &E, HPsi, false);
	std::cout << "PsiHPsi: E=" << E << " HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 2, 7, true);
	H.PsiHPsi(Psi, &E, HPsi);
	std::cout << "PsiHPsi: E=" << E << " HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 2, 7, true);
	H.HPsi(Psi, HPsi, false);
	std::cout << "HPsi: HPsi" << std::endl;
	Test<T>::PrintMatrix(HPsi, 2, 7, true);
}
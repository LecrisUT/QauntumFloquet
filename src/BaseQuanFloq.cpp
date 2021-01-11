//
// Created by Le Minh Cristian on 2020/12/28.
//

#include "BaseQuanFloq.h"
#include <iostream>

using namespace QuanFloq;
void Version() {
	std::cout << "Version "
	          << QUANFLOQ_VERSION;
}
void Version( int* major ) {
	*major = QUANFLOQ_VERSION_MAJOR;
}
void Version( int* major, int* minor ) {
	*major = QUANFLOQ_VERSION_MAJOR;
	*minor = QUANFLOQ_VERSION_MINOR;
}
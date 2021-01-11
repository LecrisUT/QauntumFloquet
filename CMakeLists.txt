cmake_minimum_required(VERSION 3.16)
project(QuanFloq VERSION 0.1)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(MKL_USE_interface "lp64" CACHE STRING "for Intel(R)64 compatible arch: ilp64/lp64 or for ia32 arch: cdecl/stdcall")
find_package(MKL REQUIRED)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED True)

include(DefaultValues)
include(BuildDefs)

add_subdirectory(src)
add_subdirectory(test EXCLUDE_FROM_ALL)

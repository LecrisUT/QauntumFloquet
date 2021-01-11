add_library(Test INTERFACE)
target_sources(Test INTERFACE Test.cpp)
add_executable(testHamil TestHamil.cpp)

target_link_libraries(Test INTERFACE vHamil ${MKL_LIBRARIES})
target_link_libraries(testHamil PUBLIC Test;Dimer
		PUBLIC BaseQuanFloq
		)
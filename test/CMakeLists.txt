add_library(Test SHARED Test.cpp)
add_executable(testHamil TestHamil.cpp)

target_link_libraries(Test PUBLIC AllBaseHamils)
target_link_libraries(testHamil PUBLIC Test;Dimer;HFDimer;FloqDimer;FloqHFDimer)
# Generate base configuration
configure_file(Conf.h.in Conf.h)

add_library(Conf INTERFACE)
target_include_directories(Conf INTERFACE "${PROJECT_BINARY_DIR}/src")
install(FILES "${PROJECT_BINARY_DIR}/src/Conf.h" TYPE INCLUDE)

add_library(BaseQuanFloq INTERFACE)
target_sources(BaseQuanFloq INTERFACE BaseQuanFloq.cpp BaseQuanFloq.h)
target_link_libraries(BaseQuanFloq INTERFACE Conf)
target_include_directories(BaseQuanFloq INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
install(FILES BaseQuanFloq.h TYPE INCLUDE)

add_subdirectory(BaseHamil)
add_subdirectory(GenHamil)
add_subdirectory(Hubbard)

## Build the libraries
## TODO: Implement user option for STATIC/SHARED library
#add_library(vHamil Conf.cpp BaseHamil/vHamil.cpp BaseHamil/vFloqHamil.cpp BaseHamil/vHFHamil.cpp BaseHamil/vFloqHFHamil.cpp)
#add_library(vDimer Hubbard/BaseDimer/vDimer.cpp Hubbard/BaseDimer/vFloqDimer.cpp Hubbard/BaseDimer/vHFDimer.cpp Hubbard/BaseDimer/vFloqHFDimer.cpp)
## TODO: Resolve MATLAB not being able to load shared library
## Temporary hard-coded Static
#add_library(Hamil GenHamil/Hamil.cpp GenHamil/FloqHamil.cpp GenHamil/HFHamil.cpp GenHamil/FloqHFHamil.cpp)
#add_library(Dimer Hubbard/Dimer.cpp Hubbard/FloqDimer.cpp Hubbard/HFDimer.cpp Hubbard/FloqHFDimer.cpp)
#
## Using target, the definitions propagate to all including targets
#target_compile_definitions(vHamil PUBLIC ${DATA_DEFS})
#
#target_include_directories(vHamil PUBLIC "${PROJECT_BINARY_DIR}/src")
#target_include_directories(Hamil PUBLIC BaseHamil)
#target_include_directories(vDimer PUBLIC BaseHamil)
#target_include_directories(Dimer PUBLIC Hubbard/BaseDimer)
#
#target_link_libraries(vHamil PUBLIC ${MKL_LIBRARIES})
#target_link_libraries(Hamil PUBLIC vHamil)
#target_link_libraries(vDimer PUBLIC vHamil)
#target_link_libraries(Dimer PUBLIC vDimer)
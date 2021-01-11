add_subdirectory(BaseDimer)
if (MODULARITY EQUAL 0)
	# Each class is it's own library file
	# Define the library targets
	add_library(DumyBaseDimer INTERFACE)
	add_library(BaseDimer Dimer.cpp)
	add_library(FloqDimer FloqDimer.cpp)
	add_library(HFDimer HFDimer.cpp)
	add_library(FloqHFDimer FloqHFDimer.cpp)
	add_library(Dimer INTERFACE)

	# Various aliases for simple compatibility
	add_library(AllDimer ALIAS Dimer)
	add_library(BaseDimerSrc ALIAS BaseDimer)
	add_library(FloqDimerSrc ALIAS FloqDimer)
	add_library(HFDimerSrc ALIAS HFDimer)
	add_library(FloqHFDimerSrc ALIAS FloqHFDimer)
	add_library(DimerSrc ALIAS Dimer)

	# link the basic components
	target_link_libraries(Dimer INTERFACE BaseDimer FloqDimer HFDimer FloqHFDimer)
	target_link_libraries(BaseDimer PRIVATE vBaseDimerSrc PUBLIC DumyBaseDimer vBaseDimer)
	target_link_libraries(FloqDimer PRIVATE vFloqDimerSrc PUBLIC DumyBaseDimer vFloqDimer)
	target_link_libraries(HFDimer PRIVATE vHFDimerSrc PUBLIC DumyBaseDimer vHFDimer)
	target_link_libraries(FloqHFDimer PRIVATE vFloqHFDimerSrc PUBLIC DumyBaseDimer vFloqHFDimer)

	target_include_directories(DumyBaseDimer INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

	# Add to install targets
	install(TARGETS BaseDimer FloqDimer HFDimer FloqHFDimer)
else ()
	# All classes in one library file
	# Define the library target
	add_library(Dimer Dimer.cpp;
			FloqDimer.cpp;
			HFDimer.cpp;
			FloqHFDimer.cpp)

	# Various aliases for simple compatibility
	add_library(BaseDimer ALIAS Dimer)
	add_library(FloqDimer ALIAS Dimer)
	add_library(HFDimer ALIAS Dimer)
	add_library(FloqHFDimer ALIAS Dimer)
	add_library(AllDimer ALIAS Dimer)
	add_library(BaseDimerSrc ALIAS Dimer)
	add_library(FloqDimerSrc ALIAS Dimer)
	add_library(HFDimerSrc ALIAS Dimer)
	add_library(FloqHFDimerSrc ALIAS Dimer)
	add_library(DimerSrc ALIAS Dimer)

	# link to virtual components
	target_link_libraries(Dimer PUBLIC vDimer PRIVATE vDimerSrc)
	target_include_directories(Dimer PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

	# Add to install target
	install(TARGETS Dimer)
endif ()
install(FILES Dimer.h FloqDimer.h HFDimer.h FloqHFDimer.h
		TYPE INCLUDE
		)
if (BUILD_VIRTUAL)
	# Build the virtual classes to be linked to
	if (MODULARITY EQUAL 0)
		# Each class is it's own library file
		# Define the library targets
		add_library(vBaseDimer vDimer.cpp)
		add_library(vFloqDimer vFloqDimer.cpp)
		add_library(vHFDimer vHFDimer.cpp)
		add_library(vFloqHFDimer vFloqHFDimer.cpp)
		add_library(vDimer INTERFACE)

		# Various aliases for simple compatibility
		add_library(vAllDimer ALIAS vDimer)
		add_library(vBaseDimerSrc ALIAS vBaseDimer)
		add_library(vFloqDimerSrc ALIAS vFloqDimer)
		add_library(vHFDimerSrc ALIAS vHFDimer)
		add_library(vFloqHFDimerSrc ALIAS vFloqHFDimer)
		add_library(vDimerSrc ALIAS vDimer)

		# link the basic components
		target_link_libraries(vBaseDimer PUBLIC vBaseHamil)
		target_include_directories(vBaseDimer PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
		target_link_libraries(vFloqDimer PUBLIC vFloqHamil;vBaseDimer)
		target_link_libraries(vHFDimer PUBLIC vHFHamil;vBaseDimer)
		target_link_libraries(vFloqHFDimer PUBLIC vFloqHFHamil;vFloqDimer;vHFDimer)
		target_link_libraries(vDimer INTERFACE vBaseDimer vFloqDimer vHFDimer vFloqHFDimer)

		# Add to install targets
		install(TARGETS vBaseDimer vFloqDimer vHFDimer vFloqHFDimer)
	else ()
		# All classes in one library file
		# Define the library target
		add_library(vDimer vDimer.cpp;
				vFloqDimer.cpp;
				vHFDimer.cpp;
				vFloqHFDimer.cpp)

		# link the basic components
		target_link_libraries(vDimer PUBLIC vHamil)
		target_include_directories(vDimer PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

		# Various aliases for simple compatibility
		add_library(vBaseDimer ALIAS vDimer)
		add_library(vFloqDimer ALIAS vDimer)
		add_library(vHFDimer ALIAS vDimer)
		add_library(vFloqHFDimer ALIAS vDimer)
		add_library(vAllDimer ALIAS vDimer)
		add_library(vBaseDimerSrc ALIAS vDimer)
		add_library(vFloqDimerSrc ALIAS vDimer)
		add_library(vHFDimerSrc ALIAS vDimer)
		add_library(vFloqHFDimerSrc ALIAS vDimer)
		add_library(vDimerSrc ALIAS vDimer)

		# Add to install target
		install(TARGETS vDimer)
	endif ()
else ()
	# Link via interface. Only final classes are built.
	# Define all interfaces
	add_library(vBaseDimer INTERFACE)
	add_library(vFloqDimer INTERFACE)
	add_library(vHFDimer INTERFACE)
	add_library(vFloqHFDimer INTERFACE)
	add_library(vDimer INTERFACE)

	# Src-ed libraries used specifically when compiling final targets
	add_library(vBaseDimerSrc INTERFACE)
	add_library(vFloqDimerSrc INTERFACE)
	add_library(vHFDimerSrc INTERFACE)
	add_library(vFloqHFDimerSrc INTERFACE)
	add_library(vDimerSrc INTERFACE)
	add_library(vAllDimer ALIAS vDimer)

	target_sources(vBaseDimerSrc INTERFACE vDimer.cpp)
	target_sources(vFloqDimerSrc INTERFACE vFloqDimer.cpp)
	target_sources(vHFDimerSrc INTERFACE vHFDimer.cpp)
	target_sources(vFloqHFDimerSrc INTERFACE vFloqHFDimer.cpp)

	# link the basic components
	target_link_libraries(vBaseDimer INTERFACE vBaseHamil)
	target_include_directories(vBaseDimer INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

	# link derived libraries
	target_link_libraries(vFloqDimer INTERFACE vFloqHamil vBaseDimer)
	target_link_libraries(vHFDimer INTERFACE vHFHamil vBaseDimer)
	target_link_libraries(vFloqHFDimer INTERFACE vFloqHFHamil vFloqDimer vHFDimer)
	target_link_libraries(vBaseDimerSrc INTERFACE vBaseHamilSrc vBaseDimer )
	target_link_libraries(vFloqDimerSrc INTERFACE vFloqHamilSrc vFloqDimer)
	target_link_libraries(vHFDimerSrc INTERFACE vHFHamilSrc vHFDimer)
	target_link_libraries(vFloqHFDimerSrc INTERFACE vFloqHFHamilSrc vFloqHFDimer)

	target_link_libraries(vDimer INTERFACE vBaseDimer vFloqDimer vHFDimer vFloqHFDimer)
	target_link_libraries(vDimerSrc INTERFACE vDimer vBaseDimerSrc vFloqDimerSrc vHFDimerSrc vFloqHFDimerSrc)
endif ()
install(FILES vDimer.h vFloqDimer.h vHFDimer.h vFloqHFDimer.h
		TYPE INCLUDE
		)
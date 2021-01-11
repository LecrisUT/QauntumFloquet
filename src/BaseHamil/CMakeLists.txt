if (BUILD_VIRTUAL)
	# Build the virtual classes to be linked to
	if (MODULARITY EQUAL 0)
		# Each class is it's own library file
		# Define the library targets
		add_library(vBaseHamil vHamil.cpp)
		add_library(vFloqHamil vFloqHamil.cpp)
		add_library(vHFHamil vHFHamil.cpp)
		add_library(vFloqHFHamil vFloqHFHamil.cpp)
		add_library(vHamil INTERFACE)

		# Various aliases for simple compatibility
		add_library(vAllHamil ALIAS vHamil)
		add_library(vBaseHamilSrc ALIAS vBaseHamil)
		add_library(vFloqHamilSrc ALIAS vFloqHamil)
		add_library(vHFHamilSrc ALIAS vHFHamil)
		add_library(vFloqHFHamilSrc ALIAS vFloqHFHamil)
		add_library(vHamilSrc ALIAS vHamil)

		# link the basic components
		target_compile_definitions(vBaseHamil PUBLIC ${DATA_DEFS})
		target_link_libraries(vBaseHamil
				PUBLIC ${MKL_LIBRARIES}
				PUBLIC Conf
				)
		target_include_directories(vBaseHamil PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

		# link deived libraries
		target_link_libraries(vFloqHamil PUBLIC vBaseHamil)
		target_link_libraries(vHFHamil PUBLIC vBaseHamil)
		target_link_libraries(vFloqHFHamil PUBLIC vFloqHamil;vHFHamil)
		target_link_libraries(vHamil INTERFACE vBaseHamil vFloqHamil vHFHamil vFloqHFHamil)

		# Add to install targets
		install(TARGETS vBaseHamil vFloqHamil vHFHamil vFloqHFHamil)
	else ()
		# All classes in one library file
		# Define the library targets
		add_library(vHamil vHamil.cpp;
				vFloqHamil.cpp;
				vHFHamil.cpp;
				vFloqHFHamil.cpp)

		# Various aliases for simple compatibility
		add_library(vBaseHamil ALIAS vHamil)
		add_library(vFloqHamil ALIAS vHamil)
		add_library(vHFHamil ALIAS vHamil)
		add_library(vFloqHFHamil ALIAS vHamil)
		add_library(vAllHamil ALIAS vHamil)
		add_library(vBaseHamilSrc ALIAS vHamil)
		add_library(vFloqHamilSrc ALIAS vHamil)
		add_library(vHFHamilSrc ALIAS vHamil)
		add_library(vFloqHFHamilSrc ALIAS vHamil)
		add_library(vHamilSrc ALIAS vHamil)

		# link the basic components
		target_compile_definitions(vHamil PUBLIC ${DATA_DEFS})
		target_link_libraries(vHamil PUBLIC ${MKL_LIBRARIES} Conf)
		target_include_directories(vHamil PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

		install(TARGETS vHamil)
	endif ()
else ()
	# Link via interface. Only final classes are built.
	# Define all interfaces
	# TODO: Separate tHamil type and add them to the libraries
	# TODO: Check if infinite recursion occurs
	add_library(vBaseHamil INTERFACE)
	add_library(vFloqHamil INTERFACE)
	add_library(vHFHamil INTERFACE)
	add_library(vFloqHFHamil INTERFACE)
	add_library(vHamil INTERFACE)
	add_library(vAllHamil ALIAS vHamil)

	# Src-ed libraries used specifically when compiling final targets
	add_library(vBaseHamilSrc INTERFACE)
	add_library(vFloqHamilSrc INTERFACE)
	add_library(vHFHamilSrc INTERFACE)
	add_library(vFloqHFHamilSrc INTERFACE)
	add_library(vHamilSrc INTERFACE)

	# Add src names for specific compilation. Libraries will link Defs and Includes
	# TODO: MKL_Libraries are redundant for non-compilation part
	target_sources(vBaseHamilSrc INTERFACE vHamil.cpp)
	target_sources(vFloqHamilSrc INTERFACE vFloqHamil.cpp)
	target_sources(vHFHamilSrc INTERFACE vHFHamil.cpp)
	target_sources(vFloqHFHamilSrc INTERFACE vFloqHFHamil.cpp )

	# link the basic components
	target_compile_definitions(vBaseHamil INTERFACE ${DATA_DEFS})
	target_link_libraries(vBaseHamil INTERFACE Conf)
	target_include_directories(vBaseHamil INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

	target_link_libraries(vHamil INTERFACE vBaseHamil vFloqHamil vHFHamil vFloqHFHamil)
	target_link_libraries(vHamilSrc INTERFACE vHamil vBaseHamilSrc vFloqHamilSrc vHFHamilSrc vFloqHFHamilSrc)

	# link deived libraries
	target_link_libraries(vBaseHamilSrc INTERFACE vBaseHamil ${MKL_LIBRARIES} )
	target_link_libraries(vFloqHamilSrc INTERFACE vFloqHamil)
	target_link_libraries(vHFHamilSrc INTERFACE vHFHamil)
	target_link_libraries(vFloqHFHamilSrc INTERFACE vFloqHFHamil)
	# TODO: Sort-of redundant since source is in different library
	target_link_libraries(vFloqHamil INTERFACE vBaseHamil)
	target_link_libraries(vHFHamil INTERFACE vBaseHamil)
	target_link_libraries(vFloqHFHamil INTERFACE vFloqHamil;vHFHamil)
endif ()
list(APPEND RTLibrary ${MKL_LIBRARIES})
install(FILES vHamil.h vFloqHamil.h vHFHamil.h vFloqHFHamil.h
		TYPE INCLUDE
		)
# TODO: Separate tHamil types and define all the possible interface libraries
if (MODULARITY EQUAL 0)
	# Each class is it's own library file
	# Define the library targets
	add_library(DumyBaseHamil INTERFACE)
	add_library(BaseHamil Hamil.cpp)
	add_library(FloqHamil FloqHamil.cpp)
	add_library(HFHamil HFHamil.cpp)
	add_library(FloqHFHamil FloqHFHamil.cpp)

	# Various aliases for simple compatibility
	add_library(Hamil INTERFACE)
	add_library(AllHamil ALIAS Hamil)
	add_library(BaseHamilSrc ALIAS BaseHamil)
	add_library(FloqHamilSrc ALIAS FloqHamil)
	add_library(HFHamilSrc ALIAS HFHamil)
	add_library(FloqHFHamilSrc ALIAS FloqHFHamil)
	add_library(HamilSrc ALIAS Hamil)

	# link the basic components
	# For simplicity everything is PUBLIC: no need to write a lot of includes
	target_link_libraries(BaseHamil PUBLIC DumyBaseHamil vBaseHamil PRIVATE vBaseHamilSrc)
	target_link_libraries(FloqHamil PUBLIC DumyBaseHamil vFloqHamil PRIVATE vFloqHamilSrc)
	target_link_libraries(HFHamil PUBLIC DumyBaseHamil vHFHamil PRIVATE vHFHamilSrc)
	target_link_libraries(FloqHFHamil PUBLIC DumyBaseHamil vFloqHFHamil PRIVATE vFloqHFHamilSrc)
	target_link_libraries(Hamil INTERFACE BaseHamil FloqHamil HFHamil FloqHFHamil)

	target_include_directories(DumyBaseHamil INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

	# Add to install targets
	install(TARGETS BaseHamil FloqHamil HFHamil FloqHFHamil)
else ()
	# All classes in one library file
	# Define the library target
	add_library(Hamil Hamil.cpp;
			FloqHamil.cpp;
			HFHamil.cpp;
			FloqHFHamil.cpp)

	# Various aliases for simple compatibility
	add_library(BaseHamil ALIAS Hamil)
	add_library(FloqHamil ALIAS Hamil)
	add_library(HFHamil ALIAS Hamil)
	add_library(FloqHFHamil ALIAS Hamil)
	add_library(AllHamil ALIAS Hamil)
	add_library(BaseHamilSrc ALIAS Hamil)
	add_library(FloqHamilSrc ALIAS Hamil)
	add_library(HFHamilSrc ALIAS Hamil)
	add_library(FloqHFHamilSrc ALIAS Hamil)
	add_library(HamilSrc ALIAS Hamil)


	# link to virtual components
	target_link_libraries(Hamil PUBLIC vHamil PRIVATE vHamilSrc)
	target_include_directories(Hamil PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

	# Add to install target
	install(TARGETS Hamil)
endif ()
install(FILES Hamil.h FloqHamil.h HFHamil.h FloqHFHamil.h
		TYPE INCLUDE
		)
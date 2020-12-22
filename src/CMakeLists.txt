add_library(AllBaseHamils INTERFACE)
add_library(Hamil SHARED Hamil/Hamil.cpp)
add_library(FloqHamil SHARED Hamil/FloqHamil.cpp )
add_library(HFHamil SHARED Hamil/HFHamil.cpp)
add_library(FloqHFHamil SHARED Hamil/FloqHFHamil.cpp)

set(BaseHeaders "PreCompOpt.h;Hamil/PreCompOpt_Hamil.h;Hamil/Hamil.h;Hamil/FloqHamil.h;Hamil/HFHamil.h;Hamil/FloqHFHamil.h")
set(DimerHeaders "Hubbard/Dimer.h;Hubbard/FloqDimer.h;Hubbard/HFDimer.h;Hubbard/FloqHFDimer.h")
#target_link_libraries(AllBaseHamils INTERFACE Hamil FloqHamil HFHamil FloqHFHamil)
#target_link_libraries(DimerHamils INTERFACE Dimer FloqDimer HFDimer FloqHFDimer)

#add_library(AllBaseHamils INTERFACE)
#add_library(DimerHamils INTERFACE)

#install(TARGETS AllBaseHamils DESTINATION lib)
#install(TARGETS DimerHamils DESTINATION lib)

#install(FILES BaseHeaders PUBLIC_HEADER)
#install(FILES DimerHeaders PUBLIC_HEADER)

add_library(Dimer SHARED Hubbard/Dimer.cpp)
add_library(HFDimer SHARED Hubbard/HFDimer.cpp )
add_library(FloqDimer SHARED Hubbard/FloqDimer.cpp)
add_library(FloqHFDimer SHARED Hubbard/FloqHFDimer.cpp)

target_link_libraries(Hamil PUBLIC ${MKL_LIBRARIES})
target_link_libraries(FloqHamil PUBLIC Hamil)
target_link_libraries(HFHamil PUBLIC Hamil)
target_link_libraries(FloqHFHamil PUBLIC HFHamil)

target_link_libraries(Dimer PUBLIC Hamil)
target_link_libraries(HFDimer PUBLIC HFHamil;Dimer)
target_link_libraries(FloqDimer PUBLIC FloqHamil;Dimer)
target_link_libraries(FloqHFDimer PUBLIC FloqHFHamil;HFDimer;FloqDimer)

install(TARGETS Hamil FloqHamil HFHamil FloqHFHamil)
install(TARGETS Dimer FloqDimer HFDimer FloqHFDimer)
install(FILES ${BaseHeaders} DESTINATION include)
install(FILES ${DimerHeaders} DESTINATION include)
cmake_minimum_required(VERSION 2.6)
project(lect1)


add_library(lect1 INTERFACE)
add_executable(SYMM    utils/symm.cpp)
add_executable(CMPVEC  utils/cmpvec.cpp)
add_executable(CG      utils/cg.cpp)
add_executable(PCG     utils/pcg.cpp)
add_executable(Convert utils/conv.cpp)

include_directories(${PROJECT_SOURCE_DIR})


install(TARGETS lect1 EXPORT lect1-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect1-targets DESTINATION cmake)

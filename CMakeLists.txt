cmake_minimum_required(VERSION 2.6)
project(lect4)


add_library(lect4 INTERFACE)
add_executable(PCG_AMG_CHEB        utils/pcg_amg_cheb.cpp)

include_directories(${PROJECT_SOURCE_DIR})


install(TARGETS lect4 EXPORT lect4-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect4-targets DESTINATION cmake)

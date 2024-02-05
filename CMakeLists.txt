cmake_minimum_required(VERSION 2.6)
project(lect5)


add_library(lect5 INTERFACE)
add_executable(PCG_AMG_GS          utils/pcg_amg_gs.cpp)
add_executable(PCG_AMG_JAC         utils/pcg_amg_jac.cpp)

include_directories(${PROJECT_SOURCE_DIR})

install(TARGETS lect5 EXPORT lect5-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect5-targets DESTINATION cmake)

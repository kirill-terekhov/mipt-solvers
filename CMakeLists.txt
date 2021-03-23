cmake_minimum_required(VERSION 2.6)
project(lect5)


add_library(lect5 INTERFACE)
add_executable(PCG_AMG_CHEB        utils/pcg_amg_cheb.cpp)
add_executable(PCG_AMG_GS          utils/pcg_amg_gs.cpp)
add_executable(PCG_AMG_JAC         utils/pcg_amg_jac.cpp)
add_executable(BICGSTAB_AMG_CHEB   utils/bicgstab_amg_cheb.cpp)
add_executable(BICGSTAB_AMG_GS     utils/bicgstab_amg_gs.cpp)
add_executable(BICGSTAB_AMG_JAC    utils/bicgstab_amg_jac.cpp)

include_directories(${PROJECT_SOURCE_DIR})


install(TARGETS lect5 EXPORT lect5-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect5-targets DESTINATION cmake)

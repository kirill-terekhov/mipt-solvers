cmake_minimum_required(VERSION 2.6)
project(lect2)


add_library(lect2 INTERFACE)
add_executable(PCG_ILDUC        utils/pcg.cpp)
add_executable(BICGSTAB_ILDUC   utils/bicgstab.cpp)

include_directories(${PROJECT_SOURCE_DIR})


install(TARGETS lect2 EXPORT lect2-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect2-targets DESTINATION cmake)

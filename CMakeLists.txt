
cmake_minimum_required(VERSION 2.6)
project(lect3)


add_library(lect3 INTERFACE)
add_executable(MG        mg.cpp)

include_directories(${PROJECT_SOURCE_DIR})


install(TARGETS lect3 EXPORT lect3-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)

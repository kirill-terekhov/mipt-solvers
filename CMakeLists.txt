cmake_minimum_required(VERSION 2.6)
project(lect6)


add_library(lect6 INTERFACE)
add_executable(MLILDUCS        utils/mlilducS.cpp)

include_directories(${PROJECT_SOURCE_DIR})


install(TARGETS lect6 EXPORT lect6-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect6-targets DESTINATION cmake)

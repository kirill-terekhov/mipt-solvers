cmake_minimum_required(VERSION 2.6)
project(lect12)


add_library(lect12 INTERFACE)
add_executable(Deadoil          utils/deadoil_nln.cpp)

include_directories(${PROJECT_SOURCE_DIR})

install(TARGETS lect12 EXPORT lect12-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect12-targets DESTINATION cmake)

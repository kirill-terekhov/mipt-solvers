cmake_minimum_required(VERSION 2.6)
project(lect13)


add_library(lect13 INTERFACE)
add_executable(Carriers  carriers.cpp)

include_directories(${PROJECT_SOURCE_DIR})

install(TARGETS lect13 EXPORT lect13-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect13-targets DESTINATION cmake)

cmake_minimum_required(VERSION 2.6)
project(lect7)


add_library(lect7 INTERFACE)
add_executable(MLILDUC        utils/mlilduc.cpp)

include_directories(${PROJECT_SOURCE_DIR})


install(TARGETS lect7 EXPORT lect7-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect7-targets DESTINATION cmake)

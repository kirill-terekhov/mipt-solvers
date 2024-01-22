cmake_minimum_required(VERSION 2.6)
project(lect10)


add_library(lect10 INTERFACE)
add_executable(BPCG_AMG_GS          utils/bpcg_amg_gs.cpp)
add_executable(VANKA                utils/vanka.cpp)


include_directories(${PROJECT_SOURCE_DIR})

install(TARGETS lect10 EXPORT lect10-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect10-targets DESTINATION cmake)

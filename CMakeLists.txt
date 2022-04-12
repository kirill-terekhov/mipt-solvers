cmake_minimum_required(VERSION 2.6)
project(lect6)


add_library(lect6 INTERFACE)
add_executable(MLILDUC        utils/mlilduc.cpp)
add_executable(SYM_MLILDUC    utils/sym_mlilduc.cpp)
add_executable(MPT_MLILDUC    utils/mpt_mlilduc.cpp)
add_executable(SINK_MLILDUC   utils/sink_mlilduc.cpp)

include_directories(${PROJECT_SOURCE_DIR})


install(TARGETS lect6 EXPORT lect6-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect6-targets DESTINATION cmake)

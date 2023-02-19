cmake_minimum_required(VERSION 2.6)
project(lect9)

add_library(lect9 INTERFACE)

add_executable(MLILDUCS            utils/mlilducS.cpp)
add_executable(MLILDUCM            utils/mlilducM.cpp)
add_executable(MPT_WRCM_MLILDUCS   utils/mpt_wrcm_mlilducS.cpp)
add_executable(MPT_WRCM_MLILDUCM   utils/mpt_wrcm_mlilducM.cpp)


include_directories(${PROJECT_SOURCE_DIR})

install(TARGETS lect9 EXPORT lect9-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect9-targets DESTINATION cmake)

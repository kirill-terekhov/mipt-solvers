cmake_minimum_required(VERSION 2.6)
project(lect11)


add_library(lect11 INTERFACE)
add_executable(CPR_TS_ILU_AMG          utils/cpr_ts_ilu_amg.cpp)
add_executable(CPR_TS_GS_AMG           utils/cpr_ts_gs_amg.cpp)
add_executable(CPR_TSGS_ILU_AMG        utils/cpr_tsgs_ilu_amg.cpp)
add_executable(CPR_TSGS_GS_AMG         utils/cpr_tsgs_gs_amg.cpp)
add_executable(FIXED_STRESS_AMG        utils/fixed_stress_amg.cpp)


include_directories(${PROJECT_SOURCE_DIR})

install(TARGETS lect11 EXPORT lect11-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect11-targets DESTINATION cmake)

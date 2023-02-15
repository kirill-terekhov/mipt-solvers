cmake_minimum_required(VERSION 2.6)
project(lect8)


option(USE_METIS "Compile with METIS" OFF)


add_library(lect8 INTERFACE)

add_executable(ILDUC           utils/ilduc.cpp)
add_executable(SYM_ILDUC       utils/sym_ilduc.cpp)
add_executable(MPT_ILDUC       utils/mpt_ilduc.cpp)
add_executable(SINK_ILDUC      utils/sink_ilduc.cpp)
add_executable(MPT_WRCM_ILDUC  utils/mpt_wrcm_ilduc.cpp)


#after multilevel lect!
#add_executable(MLILDUC            utils/mlilduc.cpp)
#add_executable(SYM_MLILDUC        utils/sym_mlilduc.cpp)
#add_executable(MPT_MLILDUC        utils/mpt_mlilduc.cpp)
#add_executable(SINK_MLILDUC       utils/sink_mlilduc.cpp)
#add_executable(MPT_WRCM_MLILDUC       utils/mpt_wrcm_mlilduc.cpp)


if(USE_METIS)
	set(LIB_DOWNLOAD_PATH "${CMAKE_BINARY_DIR}/Libraries/download")
	include(ExternalProject)

	ExternalProject_Add(GKLIB
			#GIT_REPOSITORY "https://github.com/KarypisLab/GKlib.git"
			GIT_REPOSITORY "https://github.com/kirill-terekhov/GKlib.git"
			#GIT_TAG "5ed4280f0238e39e4ca141c9e8ad466caf8a537d"
			GIT_TAG "master"
			UPDATE_DISCONNECTED 1
			PREFIX "${LIB_DOWNLOAD_PATH}"
			CMAKE_ARGS  "-DCMAKE_INSTALL_PREFIX=${LIB_DOWNLOAD_PATH} "
				"-DCMAKE_INSTALL_INCLUDEDIR=${LIB_DOWNLOAD_PATH}/include/ "
				"-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} "
				"-DCMAKE_C_FLAGS='-D_POSIX_C_SOURCE=199309L' "
				)

	ExternalProject_Get_Property(GKLIB install_dir)
	ExternalProject_Get_Property(GKLIB source_dir)
	add_library(gklib STATIC IMPORTED )
	message("location: ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}GKlib${CMAKE_STATIC_LIBRARY_SUFFIX}")
	set_target_properties(gklib PROPERTIES IMPORTED_LOCATION
		${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}GKlib${CMAKE_STATIC_LIBRARY_SUFFIX})
	add_dependencies(gklib GKLIB)
	set(GKLIB_LIBRARIES gklib) 
	set(GKLIB_DIR ${install_dir})
	set(GKLIB_INCLUDE_DIRS ${install_dir}/include)
	find_package_handle_standard_args(GKLIB DEFAULT_MSG GKLIB_LIBRARIES GKLIB_INCLUDE_DIRS)
	include_directories(${GKLIB_INCLUDE_DIRS})

	ExternalProject_Add(METIS
			#GIT_REPOSITORY "https://github.com/KarypisLab/METIS.git"
			GIT_REPOSITORY "https://github.com/kirill-terekhov/METIS.git"
			#GIT_TAG "e2e633ac7624ea286abfe1cd302e6ff16d9073f9"
			GIT_TAG "master"
			#URL "http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz"
			UPDATE_DISCONNECTED 1
			PREFIX "${LIB_DOWNLOAD_PATH}"
			CMAKE_ARGS  "-DCMAKE_INSTALL_PREFIX=${LIB_DOWNLOAD_PATH} "
				"-DCMAKE_INSTALL_INCLUDEDIR=${LIB_DOWNLOAD_PATH}/include/ "
				"-DGKLIB_PATH=${GKLIB_DIR} "
				"-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} "
				"-DCMAKE_C_FLAGS='-DIDXTYPEWIDTH=32 -DREALTYPEWIDTH=64' "
				)

	ExternalProject_Get_Property(METIS install_dir)
	ExternalProject_Get_Property(METIS source_dir)
	add_library(metis STATIC IMPORTED )

	message("location: ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}metis${CMAKE_STATIC_LIBRARY_SUFFIX}")
	set_target_properties(metis PROPERTIES IMPORTED_LOCATION
		${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}metis${CMAKE_STATIC_LIBRARY_SUFFIX})

	add_dependencies(METIS GKLIB)
	add_dependencies(metis METIS)
	add_dependencies(metis GKLIB)

	set(METIS_LIBRARIES metis gklib) 
	set(METIS_DIR ${install_dir})
	set(METIS_INCLUDE_DIRS ${source_dir}/include)

	find_package_handle_standard_args(METIS DEFAULT_MSG METIS_LIBRARIES METIS_INCLUDE_DIRS)
	
	include_directories(${METIS_INCLUDE_DIRS})
	add_definitions(-DUSE_METIS)

	add_executable(ORDER  utils/order.cpp)
	add_executable(MPT_METIS_ILDUC  utils/mpt_metis_ilduc.cpp)
	#add_executable(MPT_METIS_MLILDUC  utils/mpt_metis_mlilduc.cpp)
	target_link_libraries(ORDER PUBLIC ${METIS_LIBRARIES})
	target_link_libraries(MPT_METIS_ILDUC PUBLIC ${METIS_LIBRARIES})
	#target_link_libraries(MPT_METIS_MLILDUC PUBLIC ${METIS_LIBRARIES})
endif(USE_METIS)

include_directories(${PROJECT_SOURCE_DIR})


install(TARGETS lect8 EXPORT lect8-targets
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)
install(EXPORT lect8-targets DESTINATION cmake)

#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 3.12)

#Name your project here
project(multilayer)

#set the module directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

#default to release mode
if(NOT CMAKE_BUILD_TYPE) 
	set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

#build the executable in the binary directory on MS Visual Studio
if ( MSVC )
	SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${OUTPUT_DIRECTORY}")
	SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${OUTPUT_DIRECTORY}")
	SET( LIBRARY_OUTPUT_DIRECTORY_DEBUG "${OUTPUT_DIRECTORY}")
	SET( LIBRARY_OUTPUT_DIRECTORY_RELEASE "${OUTPUT_DIRECTORY}")
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
endif ( MSVC )


#find packages-----------------------------------
#find the pthreads package
find_package(Threads)

#find the X11 package
find_package(X11)

#find CUDA, mostly for LA stuff using cuBLAS
find_package(CUDA REQUIRED)

#find Boost
#find_package(Boost)

#find the STIM library
find_package(STIM REQUIRED)

#find LAPACK and supporting link_libraries
find_package(clapack CONFIG REQUIRED)
find_package(OpenBLAS CONFIG REQUIRED)

#include include directories
include_directories(${CUDA_INCLUDE_DIRS}
			${STIM_INCLUDE_DIRS}
)

#Assign source files to the appropriate variables to easily associate them with executables
file(GLOB SRC "src/*.cpp")

#-----------------------------Create the executable--------------------------
#-----------------------------Show all four examples-------------------------
add_executable(multilayer
	       ${SRC}
)
link_directories(${CUDA_BIN_DIRS})
target_link_libraries(multilayer ${CUDA_LIBRARIES}
			${CUDA_CUBLAS_LIBRARIES}
			${CUDA_cusparse_LIBRARY}
			${CUDA_cusolver_LIBRARY}
			${CUDA_CUFFT_LIBRARIES}
			OpenBLAS::OpenBLAS
			f2c lapack
)

		
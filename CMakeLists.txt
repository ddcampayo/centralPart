# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( cPart )

cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package

set(EIGEN3_INCLUDE_DIR "/usr/include/eigen3")
find_package(Eigen3)

if(EIGEN3_FOUND)
message(STATUS "NOTICE: Eigen library found.")
include( ${EIGEN3_USE_FILE} )
else()
message(STATUS "NOTICE: Eigen library is not found.")
endif()

set(CMAKE_MODULE_PATH "/usr/local/include/eigen/cmake/;${CMAKE_MODULE_PATH}")
 #set(CHOLMOD_LIBRARIES "/usr/include/suitesparse/")
find_package( Cholmod )
find_library(LAPACK_LIB NAMES lapack)
find_library(BLAS_LIB NAMES blas)
find_library(SS_LIB NAMES suitesparseconfig)
include_directories( ${CHOLMOD_INCLUDES})



# Creating entries for target: pPart
# ############################

add_executable( cPart
		create.cpp
		linear_aux.cpp
		move.cpp
		volumes.cpp
		draw.cpp
		linear_DD_scalar_prod.cpp
		linear_fill_matrices.cpp
		linear_vect_to_field.cpp
		linear_ps_equation.cpp
		linear_p_equation.cpp
		linear_s_equation.cpp
		number.cpp
		fields.cpp
		lloyds.cpp
# choice of main file:
		cParticles_p_s.cpp
#		cParticles_s.cpp
)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS pPart )

# Link the executable to CGAL and third-party libraries
target_link_libraries(cPart   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )


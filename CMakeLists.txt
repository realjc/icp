cmake_minimum_required(VERSION 3.0)
project(my_icp)


SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

include_directories("include/")
find_package(Eigen3 REQUIRED)

find_package(PCL 1.7.2 REQUIRED COMPONENTS common io features visualization)
add_definitions(${PCL_DEFINITIONS})

include_directories (
  ${EIGEN3_INCLUDE_DIR}
  ${PCL_INCLUDE_DIRS}
  include
  )

# Set header files 
set(
	HEADER_FILES
	include/point_to_point.hpp
	)

# Set source files
set(
	SOURCE_FILES 
    src/point_to_point.cpp 
    test/test_icp.cpp 
    )

set(LIBRARIES
    ${EIGEN3_LIBRARY}
    ${PCL_LIBRARIES}
    ${PCL_COMMON_LIBRARIES}
    ${PCL_IO_LIBRARIES}
    ${PCL_FEATURES_LIBRARIES}
    ${PCL_VISUALIZATION_LIBRARIES}
    )
# executable files
add_executable(test_icp ${SOURCE_FILES} ${HEADER_FILES})
# target libraries 
target_link_libraries(test_icp ${PCL_LIBRARY_DIRS} ${LIBRARIES})


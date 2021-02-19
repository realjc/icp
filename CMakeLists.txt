cmake_minimum_required (VERSION 2.8)
project (P2P)
include_directories(include)

find_package(Sophus REQUIRED)
include_directories(${Sophus_INCLUDE_DIRS})

# find_package(Glog REQUIRED)
# include_directories(BEFORE ${GLOG_INCLUDE_DIRS}) 

set(SRC_FILES 
    test/test_icp_p2p.cpp
    src/icp_p2p.cpp
    src/load_data.cpp)
set(HEAD_FILES
    include/data_type.hpp
    include/icp_p2p.hpp
    include/load_data.hpp)

add_executable(point2point ${SRC_FILES} ${HEAD_FILES})
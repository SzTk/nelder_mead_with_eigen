cmake_minimum_required (VERSION 3.0)
project (amoeba)

# Add cmake too path
set (CMAKE_MODULE_PATH $(PROJECT_SOURCE_DIR)/cmaketools)

# Search Eigen3
find_package (Eigen3 REQUIRED NO_MODULE)

# Define target
add_executable (amoeba amoeba_try.cpp readMatrix.cpp)

# Add include
include_directories ( ${EIGEN3_INCLUDE_DIR} ß)
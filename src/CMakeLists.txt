cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

cmake_policy(SET CMP0074 NEW)

project(pyMBA
    LANGUAGES CXX
)

# set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")


set(CMAKE_BUILD_TYPE "Release")
# set(CMAKE_MACOSX_RPATH 1)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(pybind11 REQUIRED)
find_package(mba REQUIRED)

pybind11_add_module(pyMBA pyMBA.cpp)
target_link_libraries(pyMBA PRIVATE mba::mba Eigen3::Eigen)
target_link_options(pyMBA PUBLIC
    $<$<CXX_COMPILER_ID:GNU>:-Wno-undef>>
)
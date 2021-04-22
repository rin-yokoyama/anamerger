# CMakeLists.txt for PspmtAnalysis package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(merger)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install/ CACHE PATH "install dir" FORCE)
set(PAASS_DIR_PATH $ENV{PAASSSYS} CACHE STRING "path to paass directory" )
set( MERGER_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib/ CACHE STRING "path to install lib")
set(MERGER_WORK_DIR $ENV{HOME}/merger_work/ CACHE STRING "path to working directory")
set(BRIKEN_TOOLS_DIR ~/opt/ur17b/BrikenTools CACHE PATH "path to BrikenTools dir")
set(BRIKEN_TOOLS_MERGER_DIR ${BRIKEN_TOOLS_DIR}/Merger/MergerSoft)
set(BRIKEN_TOOLS_BUILD_DIR ${BRIKEN_TOOLS_DIR}/Merger/build CACHE PATH "path to BrikenTools build dir")
set(HOME $ENV{HOME} )
configure_file(include/LibraryConfig.h.in ${CMAKE_BINARY_DIR}/generated/LibraryConfig.h)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net Proof)
find_package(yaml-cpp)
set(YAML_CPP_LIB_DIR ${YAML_CPP_CMAKE_DIR}/../../../lib CACHE PATH "yaml-cpp lib dir")

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${YAML_CPP_INCLUDE_DIR})
include_directories(${PAASS_DIR_PATH}/include)
include_directories(${CMAKE_BINARY_DIR}/generated)
include_directories(${BRIKEN_TOOLS_MERGER_DIR})
link_directories(${PAASS_DIR_PATH}/lib)
link_directories(${YAML_CPP_LIB_DIR})
link_directories(${ROOTSYS}/lib)
link_directories(${BRIKEN_TOOLS_BUILD_DIR})

#---Create  a main program using the library
add_subdirectory(source)
add_subdirectory(share)
add_subdirectory(scripts)
add_subdirectory(config)
add_subdirectory(mergerLists)
add_executable(anamerger_main source/AnamergerMain.cpp)

target_link_libraries(anamerger_main ${ROOT_LIBRARIES} ${YAML_CPP_LIBRARIES} data2TreeLib AnamergerLib )

install(TARGETS anamerger_main
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        )

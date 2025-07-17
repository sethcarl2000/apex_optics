#---------------------------------------------------------
#   I copied most of this CMakeLists.txt from: https://cliutils.gitlab.io/modern-cmake/examples/root-dict/README.html 
cmake_minimum_required(VERSION 3.11)

#---------------------------------------------------------
#   Define the project 
project(ApexOpticsLibrary VERSION 0.1
        LANGUAGES CXX)  

#---------------------------------------------------------
#   Which include directories to use. 
set(CMAKE_CXX_STANDARD
    20
    CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)

#---------------------------------------------------------
#   Find ROOT (necessary)
find_package(ROOT 6.20 REQUIRED COMPONENTS Core RIO Physics)

#---------------------------------------------------------
#   Generate the ROOT dictionary
root_generate_dictionary(G__ApexOptics NPoly.h LINKDEF LinkDef.h)

#---------------------------------------------------------
#   Add Libraries
add_library(ApexOptics SHARED NPoly.h NPoly.cxx G__ApexOptics.cxx)
target_include_directories(ApexOptics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(ApexOptics PUBLIC ROOT::Core)

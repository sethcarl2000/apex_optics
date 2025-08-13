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
#   Define sources and headers
set(sources 
    NPoly.cxx 
    ApexOptics.cxx 
    RMatrix.cxx 
    NPolyArray.cxx
    MultiLayerPerceptron.cxx
)

set(headers 
    NPoly.h 
    ApexOptics.h 
    RMatrix.h 
    NPolyArray.h
    MultiLayerPerceptron.h
)

#---------------------------------------------------------
#   Generate the ROOT dictionary
root_generate_dictionary(G__ApexOptics ${headers} LINKDEF LinkDef.h)

#---------------------------------------------------------
#   Add Libraries
add_library(ApexOptics SHARED ${headers} ${sources} G__ApexOptics.cxx)
target_include_directories(ApexOptics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(ApexOptics PUBLIC ROOT::Core)

#---------------------------------------------------------
#   Add application executable
add_executable(pick_sieveholes pick_sieveholes.cxx G__ApexOptics.cxx ${headers} ${sources})

target_include_directories(pick_sieveholes PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${ROOT_INCLUDE_DIRS}"
)

target_link_libraries(pick_sieveholes PRIVATE
    ${ROOT_LIBRARIES}
    ApexOptics
)

cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)

project(FEGrowth)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(default_build_type "Release")

function(findLib libDir libOut libName)
    # Optional arguments ARGV3-ARVG6 can be used as alernative names for the library
    find_library(TEMP NAMES ${libName} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6}
        PATHS ${${libDir}} NO_DEFAULT_PATH)
    
    if(TEMP)
        set(${libOut} ${TEMP} PARENT_SCOPE)
        unset(TEMP CACHE)
    else()
        if(WIN32)
            message(SEND_ERROR "Could not find ${libName}.lib. Check ${libDir}.")
        elseif(APPLE)
            message(SEND_ERROR "Could not find lib${libName}.so, lib${libName}.a, or lib${libName}.dylib Check ${libDir}.")
        else()
            message(SEND_ERROR "Could not find lib${libName}.so, or lib${libName}.a. Check ${libDir}")
        endif()
        unset(TEMP CACHE)
    endif()
endfunction()

if(APPLE)
    set(CMAKE_OSX_ARCHITECTURES "x86_64")

    find_path(FE_SDK
        sdk/include/FECore/stdafx.h
        PATHS /Applications/FEBioStudio/
        DOC "FEBio SDK path")

    if(NOT FE_SDK)
        message(FATAL_ERROR "Unable to find FEBio SDK")
    endif()

    set(FE_SDK_INC "${FE_SDK}/sdk/include/")
    set(FE_SDK_LIB "${FE_SDK}/sdk/lib")
    message(STATUS "FE_SDK_LIB: ${FE_SDK_LIB}")

    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY omp)
    find_package(OpenMP REQUIRED)
    
    set(FE_LIBS "")
    set(FE_LIB_NAMES fecore febiomech)
    foreach(name IN LISTS FE_LIB_NAMES)
        findLib(FE_SDK_LIB FE_LIB ${name})
        list(APPEND FE_LIBS ${FE_LIB})
        unset(FE_LIB)
    endforeach(name)
    
    add_definitions(-D__APPLE__)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
    add_compile_options(-Xpreprocessor -fopenmp)
    include_directories(${OpenMP_C_INCLUDE_DIRS})

    include_directories(${FE_SDK_INC})
    
    file(GLOB INC_FILES "src/*.h")
    file(GLOB SRC_FILES "src/*.cpp")
    
    add_library(FEGrowth SHARED ${INC_FILES} ${SRC_FILES})
    target_link_libraries(FEGrowth ${FE_LIBS})
else() 
    message(FATAL_ERROR "Building supported on macOS only")
endif()

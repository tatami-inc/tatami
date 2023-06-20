cmake_minimum_required(VERSION 3.14)

project(tatami
    VERSION 1.0.0
    DESCRIPTION "A C++ interface for assorted matrix representations"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(tatami INTERFACE)

# Switch between include directories depending on whether the downstream is
# using the build directly or is using the installed package.
target_include_directories(tatami 
    INTERFACE 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# Defining how the package is to be installed.
#install(TARGETS tatami EXPORT tatami_targets)
#
#install(DIRECTORY include/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
#
#install(EXPORT tatami_targets
#    FILE tatami_targets.cmake
#    NAMESPACE tatami::
#    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami
#)

# Building the test-related machinery, if we are compiling this library directly.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(BUILD_TATAMI_GALLERY "compiling the gallery" OFF)
    if (BUILD_TATAMI_GALLERY)
        add_subdirectory(gallery)
    endif()

    option(INSTALL_GTEST OFF)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif() 
endif()

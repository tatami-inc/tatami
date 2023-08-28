cmake_minimum_required(VERSION 3.24)

project(tatami
    VERSION 2.0.2
    DESCRIPTION "A C++ interface for assorted matrix representations"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(tatami INTERFACE)
add_library(tatami::tatami ALIAS tatami)

include(GNUInstallDirs)
target_include_directories(tatami INTERFACE 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/tatami_tatami>"
)

# Building the test-related machinery, if we are compiling this library directly.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(TATAMI_TESTS "Build tatami's test suite." ON)
else()
    option(TATAMI_TESTS "Build tatami's test suite." OFF)
endif()

option(BUILD_TATAMI_GALLERY "compiling the gallery" OFF)
if(TATAMI_TESTS)
    if (BUILD_TATAMI_GALLERY)
        add_subdirectory(gallery)
    endif()

    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif() 
endif()

# Installing for find_package.
include(CMakePackageConfigHelpers)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tatami_tatami)

install(TARGETS tatami 
    EXPORT tatamiTargets)

install(EXPORT tatamiTargets
    FILE tatami_tatamiTargets.cmake
    NAMESPACE tatami::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatamiConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatamiConfig.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami)

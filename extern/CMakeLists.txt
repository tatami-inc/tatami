cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  byteme 
  GIT_REPOSITORY https://github.com/LTLA/byteme
  GIT_TAG 8a13a7ea0ea0712c2efb31c05f96d10f9efbfd41
)

FetchContent_MakeAvailable(byteme)

cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  byteme 
  GIT_REPOSITORY https://github.com/LTLA/byteme
  GIT_TAG 2e2185a52dd474d28cab89eb79329b467b8dac8e 
)

FetchContent_MakeAvailable(byteme)

include(FetchContent)

FetchContent_Declare(
  subpar
  GIT_REPOSITORY https://github.com/LTLA/subpar
  GIT_TAG master # ^3.0.0
)

FetchContent_MakeAvailable(subpar)

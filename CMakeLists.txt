cmake_minimum_required(VERSION 3.21)

project(ContingencyTable VERSION 1.0.0 LANGUAGES CXX)

include(FetchContent)

option(CONTINGENCYTABLE_BUILD_EXAMPLES "Build ContingencyTable examples" OFF)

# Keep DataTable dependency lightweight.
set(DATATABLE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(DATATABLE_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
  DataTable
  GIT_REPOSITORY https://github.com/nathanjrussell/DataTable.git
  GIT_TAG v1.1.0
)
FetchContent_MakeAvailable(DataTable)

add_library(ContingencyTable src/ContingencyTable.cpp)
add_library(ContingencyTable::ContingencyTable ALIAS ContingencyTable)

target_compile_features(ContingencyTable PUBLIC cxx_std_17)
target_include_directories(
  ContingencyTable
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(ContingencyTable PUBLIC DataTable::DataTable)

if(CONTINGENCYTABLE_BUILD_EXAMPLES)
  add_executable(parse_csv examples/parse_csv.cpp)
  target_link_libraries(parse_csv PRIVATE DataTable::DataTable)

  add_executable(contingency_example examples/contingency_example.cpp)
  target_link_libraries(contingency_example PRIVATE ContingencyTable::ContingencyTable)
endif()


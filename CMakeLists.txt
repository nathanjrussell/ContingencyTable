cmake_minimum_required(VERSION 3.21)

project(ContingencyTable VERSION 2.0.0 LANGUAGES CXX)

include(FetchContent)

option(CONTINGENCYTABLE_BUILD_EXAMPLES "Build ContingencyTable examples" OFF)
option(CONTINGENCYTABLE_BUILD_TESTS "Build ContingencyTable tests" ON)

# Keep DataTable dependency lightweight
set(DATATABLE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(DATATABLE_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
  DataTable
  GIT_REPOSITORY https://github.com/nathanjrussell/DataTable.git
  GIT_TAG v2.0.0
)
FetchContent_MakeAvailable(DataTable)

add_library(ContingencyTable
  src/ContingencyTable.cpp
  src/FeatureSelector.cpp
)
add_library(ContingencyTable::ContingencyTable ALIAS ContingencyTable)

target_compile_features(ContingencyTable PUBLIC cxx_std_17)
target_include_directories(
  ContingencyTable
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(ContingencyTable PUBLIC DataTable::DataTable)

if(CONTINGENCYTABLE_BUILD_EXAMPLES)
  add_executable(parse_csv examples/parse_csv.cpp)
  target_link_libraries(parse_csv PRIVATE DataTable::DataTable)

  add_executable(contingency_example examples/contingency_example.cpp)
  target_link_libraries(contingency_example PRIVATE ContingencyTable::ContingencyTable)

  add_executable(find_partition examples/find_partition.cpp)
  target_link_libraries(find_partition PRIVATE ContingencyTable::ContingencyTable)

  add_executable(verify_rowfilter examples/verify_rowfilter.cpp)
  target_link_libraries(verify_rowfilter PRIVATE ContingencyTable::ContingencyTable)

  add_executable(feature_selector examples/feature_selector.cpp)
  target_link_libraries(feature_selector PRIVATE ContingencyTable::ContingencyTable)
endif()

if(CONTINGENCYTABLE_BUILD_TESTS)
  enable_testing()
  include(${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt)
endif()



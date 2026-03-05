# Tests are opt-in from the parent CMake via CONTINGENCYTABLE_BUILD_TESTS.

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -O3")

include(FetchContent)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

add_executable(ContingencyTableTests
  ${CMAKE_CURRENT_LIST_DIR}/unit_tests.cpp
)

# Keep a single top-level build output layout.
set_target_properties(ContingencyTableTests PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

target_link_libraries(ContingencyTableTests
  PRIVATE
    ContingencyTable::ContingencyTable
    DataTable::DataTable
    GTest::gtest_main
)

target_compile_features(ContingencyTableTests PRIVATE cxx_std_17)

target_compile_definitions(ContingencyTableTests
  PRIVATE
    CONTINGENCYTABLE_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

include(GoogleTest)
gtest_discover_tests(ContingencyTableTests)

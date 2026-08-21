cmake_minimum_required(VERSION 3.15)
project(xpid_moorhen) 

include(CTest)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Enable Release mode and highest level optimization (Crucial for performance)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

# 1. Find all Gemmi source files (.cpp) automatically
file(GLOB GEMMI_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/gemmi/src/*.cpp")

# 2. Compile bundled Gemmi once and share it with the CLI and contract tests.
add_library(xpid_gemmi STATIC ${GEMMI_SOURCES})

# 3. Tell CMake where to find the headers
target_include_directories(xpid_gemmi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/gemmi/include)

# 4. Find and link the system ZLIB library (required by Gemmi for reading files)
find_package(ZLIB REQUIRED)
target_link_libraries(xpid_gemmi PUBLIC ZLIB::ZLIB)

add_executable(xhpi_test main.cpp)
target_link_libraries(xhpi_test PRIVATE xpid_gemmi)

if(BUILD_TESTING)
    add_executable(xpid_contract_tests test/contract_tests.cpp)
    target_link_libraries(xpid_contract_tests PRIVATE xpid_gemmi)
    target_include_directories(xpid_contract_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_definitions(
        xpid_contract_tests PRIVATE XPID_TEST_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME xpid_contract_tests COMMAND xpid_contract_tests)
endif()

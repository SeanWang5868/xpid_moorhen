cmake_minimum_required(VERSION 3.15)
project(xpid_moorhen) 

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Enable Release mode and highest level optimization (Crucial for performance)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

# 1. Find all Gemmi source files (.cpp) automatically
file(GLOB GEMMI_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/gemmi/src/*.cpp")

# 2. Compile our main.cpp ALONGSIDE all Gemmi source files
add_executable(xhpi_test main.cpp ${GEMMI_SOURCES})

# 3. Tell CMake where to find the headers
target_include_directories(xhpi_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/gemmi/include)

# 4. Find and link the system ZLIB library (required by Gemmi for reading files)
find_package(ZLIB REQUIRED)
target_link_libraries(xhpi_test PRIVATE ZLIB::ZLIB)
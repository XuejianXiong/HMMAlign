cmake_minimum_required(VERSION 3.15)
project(hmmalign LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find pybind11
find_package(pybind11 REQUIRED)

# Define the python module (_core)
# The output will be a .so or .pyd file that python can import
pybind11_add_module(_core src/viterbi.cpp)

# Optimization flags for high-performance bioinformatics
if(MSVC)
    target_compile_options(_core PRIVATE /O2 /arch:AVX2)
else()
    target_compile_options(_core PRIVATE -O3 -march=native -ffast-math)
endif()

# Ensure the module/compiled library is placed where the python package expects it
install(TARGETS _core DESTINATION hmmalign)

# Force the build output into the package folder during development
set_target_properties(_core PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/hmmalign")
# Make a basic CMakeLists.txt file



# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required(VERSION 3.5)

# Set the project name
project (CEM)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    src/main.cpp
    src/surface.cpp
)
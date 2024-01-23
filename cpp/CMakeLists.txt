
cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED False)
project (CEM)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    main.cpp
        parse_args.h
)

add_executable(CEM ${SOURCES})
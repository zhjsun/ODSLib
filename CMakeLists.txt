cmake_minimum_required (VERSION 3.11)
project (ODSLib_Demo)

# set the compiler and options
set (CMAKE_CXX_COMPILER clang++)
set (CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -g -ggdb")
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall")

include_directories (ODSLib)

aux_source_directory(src DIR_SRCS)

add_executable(runner ${DIR_SRCS})

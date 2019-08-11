cmake_minimum_required (VERSION 3.1)
project (ODSLib_Demo)

if (UNIX)
    # set the compiler and options
    set (CMAKE_CXX_COMPILER clang++)
    set (CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -g -ggdb")
    set (CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall")

    include_directories(BEFORE .)
elseif (WIN32)
    set (CMAKE_CXX_FLAGS "/utf-8")
    include_directories(BEFORE .)
    include_directories(AFTER D:\\local\\include)
    link_directories(D:\\local\\lib D:\\local\\bin)
endif ()

aux_source_directory(./src DIR_SRCS)

add_executable(runner ${DIR_SRCS})
target_link_libraries (runner dace)

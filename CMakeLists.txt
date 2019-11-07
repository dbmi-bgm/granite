cmake_minimum_required(VERSION 3.15)
project(novoCaller)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(novoCaller
        definitions.h
        novoCaller.cpp
        tests.h)

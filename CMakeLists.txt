cmake_minimum_required(VERSION 3.7)
project(bdd_steiner_tree_enumeration)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "")
set(SOURCE_FILES
        bdd.cpp
        bdd.h
        frontier.cpp
        frontier.h
        graph.cpp
        graph.h
        main.cpp
        state.cpp
        state.h
        utility.h)

add_executable(bdd_steiner_tree_enumeration ${SOURCE_FILES})
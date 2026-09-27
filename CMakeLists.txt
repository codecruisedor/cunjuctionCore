cmake_minimum_required(VERSION 3.20)
project(collisionengine)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)


FetchContent_Declare(
    sgp4
    GIT_REPOSITORY "https://github.com/dnwrnr/sgp4.git"
    GIT_TAG        4a100d5efc11ebf472472215658cf4f5eec0b9d7 # Release v3.0.0
)
FetchContent_MakeAvailable(sgp4)
add_executable(collisionengine collisionengine.cpp)

target_link_libraries(collisionengine PRIVATE sgp4)

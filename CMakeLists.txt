cmake_minimum_required(VERSION 3.20)
project(collisionengine)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)


FetchContent_Declare(
    sgp4
    GIT_REPOSITORY "https://github.com/dnwrnr/sgp4.git"
    GIT_TAG        4a100d5efc11ebf472472215658cf4f5eec0b9d7 # Release v3.0.0
)
FetchContent_MakeAvailable(sgp4)
get_target_property(sgp4_includes sgp4 INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "sgp4 include dirs: ${sgp4_includes}")
add_executable(collisionengine collisionengine.cpp)

target_link_libraries(collisionengine PRIVATE sgp4)

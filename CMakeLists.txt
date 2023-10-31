cmake_minimum_required(VERSION 3.26)
project(ProbabilisticStop)

set(CMAKE_CXX_STANDARD 17)

add_executable(main main.cpp src/ProbabilisticStop.cpp)
target_include_directories(main PRIVATE include/)
target_link_libraries(main PRIVATE armadillo)

add_library(ProbabilisticStop src/ProbabilisticStop.cpp)
target_include_directories(ProbabilisticStop PRIVATE include/)
target_link_libraries(ProbabilisticStop PRIVATE armadillo)

enable_testing ()
add_subdirectory (test)

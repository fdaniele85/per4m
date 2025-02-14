cmake_minimum_required(VERSION 3.26)
project(ProbabilisticStop)

set(CMAKE_CXX_STANDARD 17)

add_library(ProbabilisticStop ProbabilisticStop.cpp KDE.cpp ProbabilisticFilter.cpp CircularBuffer.cpp)

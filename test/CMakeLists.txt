Include(FetchContent)

FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

FetchContent_MakeAvailable(googletest)

add_executable(tests test.cpp ../src/ProbabilisticStop.cpp)
target_link_libraries(
        tests
        GTest::gtest_main
)

target_include_directories(tests PRIVATE ../include)

include(GoogleTest)
gtest_discover_tests(tests)

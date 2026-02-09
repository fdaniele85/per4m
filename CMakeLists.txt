cmake_minimum_required(VERSION 3.26)

execute_process(
        COMMAND git describe --tags --always
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE PER4M_GIT_DESCRIBE
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
)

# Use the git describe string for display, but keep a semver for CMake's version logic.
set(PER4M_VERSION_STRING "${PER4M_GIT_DESCRIBE}")
set(PER4M_VERSION_SEMVER "${PER4M_GIT_DESCRIBE}")
if(PER4M_GIT_DESCRIBE)
        set(PER4M_VERSION_STRING "${PER4M_GIT_DESCRIBE}")
        string(REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" PER4M_VERSION_SEMVER "${PER4M_GIT_DESCRIBE}")
        if(NOT PER4M_VERSION_SEMVER)
                set(PER4M_VERSION_SEMVER "0.1.0")
        endif()
endif()

project(per4m
        VERSION ${PER4M_VERSION_SEMVER}
        LANGUAGES CXX
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

add_library(per4m
        src/ProbabilisticStop.cpp
        src/KDE.cpp
        src/ProbabilisticFilter.cpp
        src/CircularBuffer.cpp
)

add_library(per4m::per4m ALIAS per4m)

target_compile_features(per4m PUBLIC cxx_std_17)

# Se hai una cartella include/per4m/... metti quella.
# Se invece gli header sono nella root, cambia BUILD_INTERFACE di conseguenza.
target_include_directories(per4m
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# ---- Install target + header ----
install(TARGETS per4m
        EXPORT per4mTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# ---- Export targets ----
install(EXPORT per4mTargets
        FILE per4mTargets.cmake
        NAMESPACE per4m::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/per4m
)

# ---- Config package ----
configure_package_config_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/per4mConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/per4mConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/per4m
)

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/per4mConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/per4mConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/per4mConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/per4m
)

# (Opzionale) per usare find_package dal build tree senza install
export(EXPORT per4mTargets
        FILE ${CMAKE_CURRENT_BINARY_DIR}/per4mTargets.cmake
        NAMESPACE per4m::
)

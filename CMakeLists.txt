# This file is part of the ARG-Needle genealogical inference and
# analysis software suite.
# Copyright (C) 2023-2025 ARG-Needle Developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.16)
message(STATUS "Using CMake version ${CMAKE_VERSION}")

project(arg_needle LANGUAGES CXX VERSION 1.0.3)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ARG_NEEDLE_PYTHON_BINDINGS "Whether to build the python bindings" OFF)
if (ARG_NEEDLE_PYTHON_BINDINGS)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif ()

# Project settings including default build type
include(cmake/ProjectSettings.cmake)

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# Sanitiser options if supported by compiler
include(cmake/Sanitisers.cmake)
enable_sanitisers(project_settings)

# allow for static analysis options
include(cmake/StaticAnalysers.cmake)

option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" OFF)

set(arg_needle_testdata_dir ${CMAKE_CURRENT_SOURCE_DIR}/testdata)
add_definitions(-DARG_NEEDLE_TESTDATA_DIR=\"${arg_needle_testdata_dir}\")


# Python bindings
if (ARG_NEEDLE_PYTHON_BINDINGS)
    option(ARG_NEEDLE_BUILDING_FROM_PYPROJECT "Are we building from pyproject.toml (pip install)?" OFF)

    if(ARG_NEEDLE_BUILDING_FROM_PYPROJECT)
        message(STATUS "Using pybind11 from pyproject.toml build environment")
        find_package(pybind11 REQUIRED)
    else()
        message(STATUS "Using FetchContent to get pybind11")
        include(FetchContent)
        FetchContent_Declare(
                pybind11
                GIT_REPOSITORY https://github.com/pybind/pybind11.git
                GIT_TAG f5fbe867d2d26e4a0a9177a51f6e568868ad3dc8 # Version 3.0.1
        )
        FetchContent_MakeAvailable(pybind11)
    endif()
endif ()


add_subdirectory(src)

option(ARG_NEEDLE_TESTING "Enable ARG Needle unit testing" ON)
if(ARG_NEEDLE_TESTING)
    Include(FetchContent)
    FetchContent_Declare(
            Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG b3fb4b9feafcd8d91c5cb510a4775143fdbef02f # Version 3.11.0
    )

    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    include(CTest)
    include(Catch)
    add_subdirectory(test/cpp)
endif()

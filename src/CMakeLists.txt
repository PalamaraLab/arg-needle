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

find_package(Boost COMPONENTS iostreams REQUIRED)

set(
        arg_needle_hashing_src
        hashing/FileUtils.cpp
        hashing/HapData.cpp
)

set(
        arg_needle_hashing_hdr
        hashing/FileUtils.hpp
        hashing/HapData.hpp
        hashing/utils.hpp
)

add_library(arg_needle_hashing STATIC ${arg_needle_hashing_src} ${arg_needle_hashing_hdr})

target_include_directories(
        arg_needle_hashing INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/hashing>
        $<INSTALL_INTERFACE:>
)

set_target_properties(arg_needle_hashing PROPERTIES PUBLIC_HEADER "${arg_needle_hashing_hdr}")

target_link_libraries(arg_needle_hashing PRIVATE Boost::headers Boost::iostreams)
target_link_libraries(arg_needle_hashing PRIVATE project_warnings)

if (PYTHON_BINDINGS)
    set_target_properties(arg_needle_hashing PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    pybind11_add_module(arg_needle_hashing_pybind hashing/pybind.cpp)
    target_link_libraries(arg_needle_hashing_pybind PRIVATE arg_needle_hashing)
endif ()

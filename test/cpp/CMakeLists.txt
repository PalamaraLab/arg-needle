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

set(
        test_files
        test_file_utils.cpp
)

add_executable(cpp_tests ${test_files})
target_link_libraries(cpp_tests PRIVATE arg_needle_hashing Catch2::Catch2WithMain)

catch_discover_tests(cpp_tests)

set(ARG_NEEDLE_TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(ARG_NEEDLE_RESOURCES_DIR "${CMAKE_SOURCE_DIR}/src/arg_needle/resources")

add_definitions(-DARG_NEEDLE_TEST_DIR=\"${ARG_NEEDLE_TEST_DIR}\")
add_definitions(-DARG_NEEDLE_RESOURCES_DIR=\"${ARG_NEEDLE_RESOURCES_DIR}\")

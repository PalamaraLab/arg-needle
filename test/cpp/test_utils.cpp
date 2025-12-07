/*
  This file is part of the ARG-Needle genealogical inference and
  analysis software suite.
  Copyright (C) 2023-2025 ARG-Needle Developers.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "HapData.hpp"
#include "utils.hpp"

using Catch::Matchers::ContainsSubstring;

void test_throw() {
  throw std::logic_error(MAKE_ERROR("Something went wrong"));
}


TEST_CASE("make_error", "[utils]") {

  REQUIRE_THROWS_WITH(test_throw(),
                      ContainsSubstring( "test_utils.cpp:" ) && ContainsSubstring( "Something went wrong" ));

  REQUIRE_THROWS_WITH(HapData("banana", ""),
                      ContainsSubstring( "HapData.cpp:" ) && ContainsSubstring( "Mode not recognized" ));
}


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

#include "FileUtils.hpp"


TEST_CASE( "FileUtils::fileExists", "[test_file_utils]" ) {
  REQUIRE(FileUtils::fileExists(ARG_NEEDLE_TEST_DIR "/CMakeLists.txt") == true);
  REQUIRE(FileUtils::fileExists(ARG_NEEDLE_TEST_DIR "/file_that_does_not_exist") == false);
}

TEST_CASE( "FileUtils::AutoGzIfstream", "[test_file_utils]")
{
  SECTION("open and close gz file")
  {
    REQUIRE(FileUtils::fileExists(ARG_NEEDLE_RESOURCES_DIR "/30-100-2000_CEU.decodingQuantities.gz") == true);
    FileUtils::AutoGzIfstream gz_file;
    gz_file.openOrExit(ARG_NEEDLE_RESOURCES_DIR "/30-100-2000_CEU.decodingQuantities.gz");
    gz_file.close();
  }

  SECTION("count lines in file")
  {
    REQUIRE(FileUtils::fileExists(ARG_NEEDLE_RESOURCES_DIR "/30-100-2000_CEU.decodingQuantities.gz") == true);
    REQUIRE(FileUtils::AutoGzIfstream::lineCount(ARG_NEEDLE_RESOURCES_DIR "/30-100-2000_CEU.decodingQuantities.gz") == 35245);
  }

  SECTION("extract line from a file")
  {
    REQUIRE(FileUtils::fileExists(ARG_NEEDLE_RESOURCES_DIR "/30-100-2000_CEU.decodingQuantities.gz") == true);

    FileUtils::AutoGzIfstream gz_file;
    gz_file.openOrExit(ARG_NEEDLE_RESOURCES_DIR "/30-100-2000_CEU.decodingQuantities.gz");

    std::string first_line;
    FileUtils::getline(gz_file, first_line);
    gz_file.close();

    REQUIRE(first_line == "TransitionType");
  }
}
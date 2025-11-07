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

// The below code was copied and modified from the Eagle software file
// https://github.com/poruloh/Eagle/blob/master/src/FileUtils.hpp
// developed by Po-Ru Loh and released under the GNU General Public
// License v3.0 (GPLv3).
//
// The license file can be found at 3rd_party/Eagle/COPYING from the
// root of this repository.

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "FileUtils.hpp"

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace FileUtils {

using std::cerr;
using std::endl;
using std::string;
using std::vector;

bool fileExists(const std::string& name) {
  std::ifstream f(name.c_str());
  return f.good();
}

int AutoGzIfstream::lineCount(const std::string& file) {
  AutoGzIfstream fin;
  fin.openOrExit(file);
  int ctr = 0;
  string line;
  while (getline(fin, line))
    ctr++;
  return ctr;
}

void AutoGzIfstream::openOrExit(const std::string& file, std::ios_base::openmode mode) {
  fin.open(file.c_str(), mode);
  if (!fin) {
    cerr << "ERROR: Unable to open file: " << file << endl;
    exit(1);
  }
  if ((int) file.length() > 3 && file.substr(file.length() - 3) == ".gz")
    boost_in.push(boost::iostreams::gzip_decompressor());
  boost_in.push(fin);
}

void AutoGzIfstream::close() {
  fin.close();
  boost_in.reset();
}

AutoGzIfstream::operator bool() const {
  return !boost_in.fail();
}

AutoGzIfstream& getline(AutoGzIfstream& in, std::string& s) {
  std::getline(in.boost_in, s);
  return in;
}

} // namespace FileUtils

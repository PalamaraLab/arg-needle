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

#include "FileUtils.hpp"

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

namespace FileUtils {

    struct AutoGzIfstream::Impl {
        boost::iostreams::filtering_istream boost_in;
        std::ifstream fin;
    };

    AutoGzIfstream::AutoGzIfstream() : pimpl(std::make_unique<Impl>()) {
    }

    AutoGzIfstream::~AutoGzIfstream() noexcept = default;

    bool fileExists(const std::filesystem::path &file) {
        std::ifstream f(file.c_str());
        return f.good();
    }

    int AutoGzIfstream::lineCount(const std::filesystem::path &file) {
        AutoGzIfstream fin;
        fin.openOrExit(file);
        int ctr = 0;
        std::string line;
        while (getline(fin, line)) {
            ctr++;
        }
        fin.close();
        return ctr;
    }

    void AutoGzIfstream::openOrExit(const std::filesystem::path &file, std::ios_base::openmode mode) {
        pimpl->fin.open(file.c_str(), mode);
        if (!pimpl->fin) {
            std::cerr << "ERROR: Unable to open file: " << file << std::endl;
            exit(1);
        }
        if (file.extension() == ".gz") {
            pimpl->boost_in.push(boost::iostreams::gzip_decompressor());
        }
        pimpl->boost_in.push(pimpl->fin);
    }

    void AutoGzIfstream::close() {
        pimpl->fin.close();
        pimpl->boost_in.reset();
    }

    AutoGzIfstream::operator bool() const noexcept {
        return !pimpl->boost_in.fail();
    }

    AutoGzIfstream &getline(AutoGzIfstream &in, std::string &s) {
        std::getline(in.pimpl->boost_in, s);
        return in;
    }

    AutoGzIfstream &AutoGzIfstream::operator>>(std::string &x) {
        pimpl->boost_in >> x;
        return *this;
    }
} // namespace FileUtils

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

#ifndef FILEUTILS_HPP
#define FILEUTILS_HPP

#include <filesystem>
#include <memory>
#include <string>

namespace FileUtils {
    /**
     * @brief Check whether a given file exists on disk.
     *
     * @param file Path to the file to check.
     * @return true if the file exists, false otherwise.
     */
    bool fileExists(const std::filesystem::path &file);

    /**
     * @class AutoGzIfstream
     * @brief Stream wrapper that transparently reads either plain-text or gzip-compressed files.
     *
     * AutoGzIfstream detects whether an input file is compressed (.gz) and automatically
     * opens it appropriately. It behaves similarly to std::ifstream but supports reading
     * gzip-compressed streams without requiring explicit decompression by the caller.
     *
     * Internally uses a pimpl to hide implementation details and avoid exposing boost
     * libraries at the interface level.
     */
    class AutoGzIfstream {
        struct Impl;
        std::unique_ptr<Impl> pimpl;

    public:
        /**
         * @brief Construct an unopened AutoGzIfstream.
         */
        AutoGzIfstream();

        /**
         * @brief Destructor closes the stream if open and releases internal resources.
         */
        ~AutoGzIfstream() noexcept;

        /**
         * @brief Count the number of lines in a file (supports gzipped and plain files).
         *
         * @param file Path to the file whose line count will be computed.
         * @return Number of lines in the file.
         */
        [[nodiscard]] static int lineCount(const std::filesystem::path &file);

        /**
         * @brief Open a file for reading or exit the program if opening fails.
         *
         * Automatically detects gzip compression based on file contents.
         *
         * @param file Path to the file to open.
         * @param mode Stream opening mode (defaults to std::ios::in).
         */
        void openOrExit(const std::filesystem::path &file,
                        std::ios_base::openmode mode = std::ios::in);

        /**
         * @brief Close the underlying stream.
         */
        void close();

        /**
         * @brief Read whitespace-delimited input into a string via the extraction operator.
         *
         * @param x Output string that will receive the parsed token.
         * @return Reference to this stream.
         */
        AutoGzIfstream &operator>>(std::string &x);

        /**
         * @brief Boolean conversion indicating whether the stream is currently valid.
         *
         * Allows usage in conditions such as:
         * @code
         * if (stream) { ... }
         * @endcode
         *
         * @return true if the stream is open and in a good state, false otherwise.
         */
        [[nodiscard]] explicit operator bool() const noexcept;

        /**
         * @brief Friend declaration enabling getline(AutoGzIfstream&, ...).
         */
        friend AutoGzIfstream &getline(AutoGzIfstream &in, std::string &s);
    };

    /**
     * @brief Read a full line from an AutoGzIfstream into a string.
     *
     * Supports both compressed and uncompressed input sources.
     *
     * @param in Stream to read from.
     * @param s Output string receiving the line (without delimiter).
     * @return Reference to the stream.
     */
    AutoGzIfstream &getline(AutoGzIfstream &in, std::string &s);
} // namespace FileUtils

#endif

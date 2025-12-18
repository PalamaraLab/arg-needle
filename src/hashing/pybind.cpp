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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>
#include <stdexcept>
#include <string>

#include "HapData.hpp"
#include "utils.hpp"

namespace py = pybind11;
using std::string;

PYBIND11_MODULE(arg_needle_hashing_pybind, m) {
  py::class_<HapData>(m, "HapData")
      .def(py::init<string, string, unsigned int, string, bool>(), "Initialize HapData",
           py::arg("mode"), py::arg("file_root_path"), py::arg("word_size") = 64,
           py::arg("map_file_path") = "", py::arg("fill_sites") = true)
      .def_readonly("num_haps", &HapData::num_haps)
      .def_readonly("num_sites", &HapData::num_sites)
      .def_readonly("word_size", &HapData::word_size)
      .def_readonly(
          "hashed_hap_ids", &HapData::hashed_hap_ids) // conversion from unordered_set to set
      .def_readonly(
          "physical_positions", &HapData::physical_positions) // conversion from vector to list
      .def_readonly(
          "genetic_positions", &HapData::genetic_positions) // conversion from vector to list
      .def_readonly("site_mafs", &HapData::site_mafs)       // conversion from vector to list
      .def("add_to_hash", &HapData::add_to_hash, py::arg("hap_id"))
      .def("get_closest_cousins", &HapData::get_closest_cousins, py::arg("hap_id"), py::arg("k"),
           py::arg("tolerance") = 0, py::arg("window_size_genetic") = 0,
           "Get K closest cousins to this one using hashing.")
      .def("print_hap", &HapData::print_hap, py::arg("hap_id"))
      .def("print_hashes", &HapData::print_hashes)
      .def("print_word_match_diagram", &HapData::print_word_match_diagram, py::arg("hap_id1"),
           py::arg("hap_id2"))
      .def(
          "physical_position_at",
          [](const HapData& data, size_t site) {
            if (site >= data.num_sites) {
              throw std::logic_error(MAKE_ERROR("Out of bounds site."));
            }
            return data.physical_positions[site];
          },
          py::arg("site"))
      .def(
          "genetic_position_at",
          [](const HapData& data, size_t site) {
            if (site >= data.num_sites) {
              throw std::logic_error(MAKE_ERROR("Out of bounds site."));
            }
            return data.genetic_positions[site];
          },
          py::arg("site"))
      .def(
          "site_maf_at",
          [](const HapData& data, size_t site) {
            if (site >= data.num_sites) {
              throw std::logic_error(MAKE_ERROR("Out of bounds site."));
            }
            return data.site_mafs[site];
          },
          py::arg("site"))
      .def(
          "sample_name",
          [](const HapData& data, size_t hap_id) {
            if (hap_id >= data.num_haps) {
              throw std::logic_error(MAKE_ERROR("Out of bounds hap_id."));
            }
            return data.sample_names[hap_id];
          },
          py::arg("hap_id"))
      .def("__repr__", [](const HapData& data) {
        std::ostringstream oss;
        oss << data;
        return oss.str();
      });
}

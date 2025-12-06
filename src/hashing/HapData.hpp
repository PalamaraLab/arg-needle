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

#ifndef ARG_NEELE_HAP_DATA_HPP
#define ARG_NEELE_HAP_DATA_HPP

#include <iostream>
#include <stdint.h>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using std::ostream;
using std::pair;
using std::string;
using std::tuple;
using std::unordered_map;
using std::unordered_set;
using std::vector;

struct Window {
  size_t start, end, index; // end is inclusive
  friend bool operator<(const Window& a, const Window& b) {
    if (a.start == b.start) {
      return a.end < b.end;
    }
    else {
      return a.start < b.start;
    }
  }
};

enum class HapDataMode { sequence, array };

class HapData {

public:
  typedef uint64_t word_type;
  unsigned int num_haps = 0;
  unsigned int num_sites = 0;
  unsigned int word_size;
  HapDataMode data_mode;
  vector<unsigned long> physical_positions;
  vector<double> genetic_positions;
  vector<float> site_mafs;
  vector<string> sample_names;
  vector<vector<bool>> sites;
  vector<vector<word_type>> words;

  vector<unordered_map<word_type, vector<size_t>>> hashes;
  unordered_set<size_t> hashed_hap_ids;

  HapData(string mode, string file_root_path, unsigned int _word_size = 64,
          string map_file_path = "", bool fill_sites = true);
  ~HapData();
  void add_to_hash(size_t hap_id);
  vector<tuple<size_t, size_t, vector<pair<size_t, double>>>>
  get_closest_cousins(size_t hap_id, unsigned int k, unsigned int tolerance = 0,
                      double window_size_genetic = 0);
  void print_hap(size_t hap_id);
  void print_hashes();
  void print_word_match_diagram(size_t hap_id1, size_t hap_id2);
  friend ostream& operator<<(ostream& os, const HapData& data);
};

#endif // ARG_NEELE_HAP_DATA_HPP

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

#include <algorithm>
#include <cassert>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "FileUtils.hpp"
#include "HapData.hpp"
#include "utils.hpp"


HapData::HapData(std::string mode, std::string file_root_path, unsigned int _word_size, std::string map_file_path,
                 bool fill_sites)
    : word_size(_word_size) {
  if (mode == "sequence") {
    data_mode = HapDataMode::sequence;
  }
  else if (mode == "array") {
    data_mode = HapDataMode::array;
  }
  else {
    throw std::logic_error(make_error("Mode not recognized."));
  }

  if (sizeof(word_type) != 8) {
    throw std::logic_error(make_error("Expected word_type to be 8 bytes (64 bits)."));
  }
  if (sizeof(1ull) < 8) {
    throw std::logic_error(
        make_error("Expected unsigned long long to be at least 8 bytes (64 bits)."));
  }
  if (word_size > 64 || word_size <= 0) {
    throw std::logic_error(make_error("Out of bounds word size."));
  }

  std::string line;
  std::stringstream ss;

  // read in .sample[s] file
  FileUtils::AutoGzIfstream file_samples;
  if (FileUtils::fileExists(file_root_path + ".samples")) {
    file_samples.openOrExit(file_root_path + ".samples");
  }
  else if (FileUtils::fileExists(file_root_path + ".sample")) {
    file_samples.openOrExit(file_root_path + ".sample");
  }
  else {
    std::cerr << "ERROR. Could not find sample file in " + file_root_path + ".sample[s]" << std::endl;
    exit(1);
  }

  while (getline(file_samples, line)) {
    std::vector<std::string> splitStr;
    std::istringstream iss(line);
    std::string buf;
    while (iss >> buf)
      splitStr.push_back(buf);

    // Skip first two lines (header) if present
    if ((splitStr[0] == "ID_1" && splitStr[1] == "ID_2" && splitStr[2] == "missing") ||
        (splitStr[0] == "0" && splitStr[1] == "0" && splitStr[2] == "0")) {
      continue;
    }

    sample_names.push_back(splitStr[0]);
    sample_names.push_back(splitStr[1]);
  }
  num_haps = sample_names.size();
  file_samples.close();

  // Parse .map[.gz] file
  FileUtils::AutoGzIfstream file_map;
  if (!map_file_path.empty()) {
    // Attempt to read in .map[.gz] file
    if (FileUtils::fileExists(map_file_path)) {
      file_map.openOrExit(map_file_path);
    }
    else {
      std::cerr << "ERROR. Could not open map file " + map_file_path + ", no such file" << std::endl;
      exit(1);
    }
  }
  else {
    // If no map file is specified, default to file_root_path.map[.gz]
    if (FileUtils::fileExists(file_root_path + ".map.gz")) {
      file_map.openOrExit(file_root_path + ".map.gz");
    }
    else if (FileUtils::fileExists(file_root_path + ".map")) {
      file_map.openOrExit(file_root_path + ".map");
    }
    else {
      std::cerr << "ERROR. Could not find map file in " + file_root_path + ".map.gz or " +
                  file_root_path + ".map"
           << std::endl;
      exit(1);
    }
  }
  std::string map_field[4];
  while (getline(file_map, line)) {
    ss.clear();
    ss.str(line);
    ss >> map_field[0] >> map_field[1] >> map_field[2] >> map_field[3];
    genetic_positions.push_back(stod(map_field[2]));
    physical_positions.push_back(stoul(map_field[3]));
  }
  num_sites = genetic_positions.size();
  file_map.close();

  // read in .hap[s][.gz] file
  FileUtils::AutoGzIfstream file_hap;
  if (FileUtils::fileExists(file_root_path + ".hap.gz")) {
    file_hap.openOrExit(file_root_path + ".hap.gz");
  }
  else if (FileUtils::fileExists(file_root_path + ".hap")) {
    file_hap.openOrExit(file_root_path + ".hap");
  }
  else if (FileUtils::fileExists(file_root_path + ".haps.gz")) {
    file_hap.openOrExit(file_root_path + ".haps.gz");
  }
  else if (FileUtils::fileExists(file_root_path + ".haps")) {
    file_hap.openOrExit(file_root_path + ".haps");
  }
  else {
    std::cerr << "ERROR. Could not find hap file in " + file_root_path + ".hap.gz, " + file_root_path +
                ".hap, " + ".haps.gz, or " + file_root_path + ".haps"
         << std::endl;
    exit(1);
  }

  if (fill_sites) {
    sites = std::vector<std::vector<bool>>(num_haps, std::vector<bool>());
  }
  words = std::vector<std::vector<word_type>>(num_haps, std::vector<word_type>());
  std::string marker_id;
  unsigned long int marker_pos;
  char al[2], inp;
  unsigned int site_id = 0u;
  while (getline(file_hap, line)) {
    // read the meta data
    ss.clear();
    ss.str(line);
    ss >> map_field[0] >> marker_id >> marker_pos >> al[0] >> al[1];
    if (map_field[0].empty()) {
      continue;
    }

    if (site_id % word_size == 0) {
      for (size_t hap_id = 0; hap_id < num_haps; ++hap_id) {
        words[hap_id].push_back(0);
      }
    }

    int maf_ctr = 0;
    if (fill_sites) {
      for (size_t hap_id = 0; hap_id < num_haps; ++hap_id) {
        ss >> inp;
        if (inp == '1') {
          ++maf_ctr;
          sites[hap_id].push_back(true);
          // important to use 1ull, not just 1!
          words[hap_id][site_id / word_size] ^= (1ull << (site_id % word_size));
        }
        else {
          sites[hap_id].push_back(false);
        }
      }
    }
    else {
      for (size_t hap_id = 0; hap_id < num_haps; ++hap_id) {
        ss >> inp;
        if (inp == '1') {
          ++maf_ctr;
          // important to use 1ull, not just 1!
          words[hap_id][site_id / word_size] ^= (1ull << (site_id % word_size));
        }
      }
    }
    float maf = static_cast<float>(maf_ctr) / static_cast<float>(num_haps);
    if (maf > 0.5f) {
      maf = 1.f - maf;
    }
    site_mafs.push_back(maf);
    ++site_id;
  }
  file_hap.close();
}

void HapData::add_to_hash(size_t hap_id) {
  if (hashed_hap_ids.find(hap_id) != hashed_hap_ids.end()) {
    throw std::logic_error(make_error("This haplotype has already been hashed."));
  }
  if (hap_id >= num_haps) {
    throw std::logic_error(make_error("Haplotype ID out of bounds."));
  }

  if (hashes.empty()) {
    hashes = std::vector<std::unordered_map<word_type, std::vector<size_t>>>(
        words[hap_id].size(), std::unordered_map<word_type, std::vector<size_t>>());
  }

  for (size_t i = 0; i < words[hap_id].size(); ++i) {
    std::vector<size_t>& hash_value =
        hashes[i][words[hap_id][i]]; // creates if not present, only hashes once
    hash_value.push_back(hap_id);
  }

  hashed_hap_ids.insert(hap_id);
}

void HapData::print_hap(size_t hap_id) {
  if (hap_id >= num_haps) {
    throw std::logic_error(make_error("Haplotype ID out of bounds."));
  }
  std::cout << "Bits for hap_id = " << hap_id << std::endl;
  for (size_t site_id = 0; site_id < num_sites; ++site_id) {
    std::cout << sites[hap_id][site_id];
    if (site_id % word_size == word_size - 1) {
      std::cout << " ";
    }
  }
  std::cout << std::endl;

  std::cout << "Words (hex) for hap_id = " << hap_id << std::endl;
  std::cout << std::hex << std::showbase;
  for (auto const& word : words[hap_id]) {
    std::cout << word << " ";
  }
  std::cout << std::endl;
  std::cout << std::dec << std::noshowbase;

  // cout << "Words (decimal)" << endl;
  // for (auto const &word : words[hap_id]) {
  //   cout << word << " ";
  // }
  // cout << endl;
}

void HapData::print_hashes() {
  for (size_t i = 0; i < hashes.size(); ++i) {
    std::cout << "Hash for word " << i << " of " << hashes.size() << std::endl;
    for (auto const& map_entry : hashes[i]) {
      unsigned long num_bits = word_size;
      if (i == hashes.size() - 1) {
        num_bits = ((num_sites - 1ul) % word_size) + 1ul;
      }
      for (size_t j = 0; j < num_bits; ++j) {
        std::cout << ((map_entry.first >> j) & 1);
      }
      std::cout << ":";
      for (const size_t id : map_entry.second) {
        std::cout << " " << id;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void HapData::print_word_match_diagram(size_t hap_id1, size_t hap_id2) {
  if (hap_id1 >= num_haps || hap_id2 >= num_haps) {
    throw std::logic_error(make_error("Haplotype ID out of bounds."));
  }
  for (size_t i = 0; i < words[hap_id1].size(); ++i) {
    if (i != 0) {
      if (i % 100 == 0) {
        std::cout << std::endl;
      }
      if (i % 25 == 0) {
        std::cout << std::endl;
      }
      else if (i % 5 == 0) {
        std::cout << " ";
      }
    }
    if (words[hap_id1][i] == words[hap_id2][i]) {
      std::cout << "x";
    }
    else {
      std::cout << "_";
    }
  }
  std::cout << std::endl;
}

std::vector<std::tuple<size_t, size_t, std::vector<std::pair<size_t, double>>>>
HapData::get_closest_cousins(size_t hap_id, unsigned int k, unsigned int tolerance,
                             double window_size_genetic) {

  // find the windows
  std::vector<Window> windows; // Window defined in HapData.hpp
  size_t num_words = words[hap_id].size();
  if (window_size_genetic <= 0) {
    // make a new window for each and every word
    for (size_t j = 0; j < num_words; ++j) {
      Window w{};
      w.start = j;
      w.end = j + 1;
      w.index = j;
      windows.push_back(w);
    }
  }
  else {
    size_t start_word = 0;
    double start_genetic = genetic_positions[0];
    size_t window_index = 0;
    for (size_t j = 0; j < num_words; ++j) {
      size_t last_word_site = std::min<size_t>((j + 1) * word_size - 1, num_sites - 1);
      // explanation: we need to leave enough room for the last window
      if (j == num_words - 1 ||
          (genetic_positions[last_word_site] - start_genetic >= window_size_genetic &&
           genetic_positions[num_sites - 1] - genetic_positions[last_word_site + 1] >=
               window_size_genetic)) {
        Window w{};
        w.start = start_word;
        w.end = j + 1;
        w.index = window_index;
        windows.push_back(w);
        window_index += 1;
        start_word = j + 1;
        if (last_word_site + 1 < num_sites) {
          start_genetic = genetic_positions[last_word_site + 1];
        }
      }
    }
  }

  std::vector<size_t> words_to_windows;
  for (size_t i = 0; i < windows.size(); ++i) {
    Window w = windows[i];
    for (size_t j = w.start; j < w.end; ++j) {
      words_to_windows.push_back(i);
    }
  }

  // how high each sample scores in each window
  // we only record samples that have matched
  std::vector<std::unordered_map<size_t, size_t>> window_scores(
      windows.size(), std::unordered_map<size_t, size_t>());
  // stretches of matching material separated by 2*k + 1 fillers, where k is the number
  // of mismatches, max size defined by 2*tolerance + 1
  std::vector<std::deque<std::pair<size_t, size_t>>> stretches(hap_id, std::deque<std::pair<size_t, size_t>>());

  // size_t num_overall_matches = 0;
  for (size_t i = 0; i < num_words; ++i) {
    // in some cases, the word does not yet exist in the hashmap
    if (hashes[i].find(words[hap_id][i]) != hashes[i].end()) {
      const std::vector<size_t>& matches = hashes[i].find(words[hap_id][i])->second;
      for (auto v : matches) {
        // check the end of stretches to figure out what to do
        if (stretches[v].empty()) {
          stretches[v].emplace_back(i, i + 1); // end is exclusive
        }
        else {
          std::pair<size_t, size_t>& back_pair = stretches[v].back();
          if (back_pair.second == i) {
            back_pair.second = i + 1; // end is exclusive
          }
          else {
            // we transform the number of mismatches, k, to 2*k - 1
            // then if the total number of mismatches is T, the total length L
            // of the stretches vector satisfies L = 2*T + 1 (assuming we start
            // and end with a match)
            // k = tolerance + 1 is the maximum we go up to since 2*k - 1 = 2*tolerance + 1
            // (actually all we need is 2*tolerance, but yeah better safe than sorry)
            size_t num_mismatches = std::min<size_t>(tolerance + 1, i - back_pair.second);
            size_t num_to_push = 2 * num_mismatches - 1;
            for (size_t push_reps = 0; push_reps < num_to_push; ++push_reps) {
              stretches[v].emplace_back(0, 0); // "NaN" value
            }
            stretches[v].emplace_back(i, i + 1); // end is exclusive
            // note: an improved algorithm is to pop first and process, then push
            // that way we only ever have to look at the front and the back
          }
        }

        // pop_front to get to size 2*tolerance + 1
        while (stretches[v].size() > 2 * tolerance + 1) {
          std::pair<size_t, size_t>& item = stretches[v].front();
          if (item.second != 0) {
            size_t range_start = item.first;
            // old version was buggy
            // size_t range_end = stretches[v].back().second;
            // new version is not ideal for complexity if tolerance is large
            size_t range_end = 0;
            for (size_t j = 0; j < 2 * tolerance + 1; ++j) {
              range_end = std::max(range_end, stretches[v][j].second);
            }
            size_t range_size = range_end - range_start;

            // we're given a half-open range [range_start, range_end)
            // we want to get the windows that overlap with this range
            // let's say windows are [0, 5), [5, 10), [10, 15), [15, 20)
            // if our range is [4, 14), we want [0, 5) to [10, 15) inclusive
            // if our range is [5, 15), we want [5, 10) to [10, 15) inclusive
            // if our range is [6, 16), we want [5, 10) to [15, 20) inclusive
            for (size_t window_index = words_to_windows[range_start];
                 window_index <= words_to_windows[range_end - 1]; ++window_index) {
              size_t& best_len =
                  window_scores[window_index][v]; // creates if not present, only hashes once
              if (range_size > best_len) {
                best_len = range_size;
              }
            }
          }
          stretches[v].pop_front();
        }
      }
    }
  }

  // go over all the stretches and pop_front
  for (size_t v = 0; v < hap_id; ++v) {
    while (!stretches[v].empty()) {
      std::pair<size_t, size_t> item = stretches[v].front();
      if (item.second != 0) {
        size_t range_start = item.first;
        // old version is buggy in general, but should work in this case
        size_t range_end = stretches[v].back().second;
        // new version is not ideal for complexity if tolerance is large
        // size_t range_end = 0;
        // for (size_t j = 0; j < stretches[v].size(); ++j) {
        //   range_end = std::max(range_end, stretches[v][j].second);
        // }
        size_t range_size = range_end - range_start;

        // we're given a half-open range [range_start, range_end)
        // we want to get the windows that overlap with this range
        // let's say windows are [0, 5), [5, 10), [10, 15), [15, 20)
        // if our range is [4, 14), we want [0, 5) to [10, 15) inclusive
        // if our range is [5, 15), we want [5, 10) to [10, 15) inclusive
        // if our range is [6, 16), we want [5, 10) to [15, 20) inclusive
        for (size_t window_index = words_to_windows[range_start];
             window_index <= words_to_windows[range_end - 1]; ++window_index) {
          size_t& hash_value =
              window_scores[window_index][v]; // creates if not present, only hashes once
          if (range_size > hash_value) {
            hash_value = range_size;
          }
        }
      }
      stretches[v].pop_front();
    }
  }

  // take the values in window_scores and sort to find top k
  std::vector<std::tuple<size_t, size_t, std::vector<std::pair<size_t, double>>>> results;
  for (const Window& w : windows) {
    size_t window_start_site = w.start * word_size;
    size_t window_end_site = std::min<size_t>(w.end * word_size - 1, num_sites - 1);

    std::vector<std::pair<double, size_t>> stats;
    for (const auto& map_entry : window_scores[w.index]) {
      size_t map_entry_hap_id = map_entry.first;
      auto score = static_cast<double>(map_entry.second);
      stats.emplace_back(score, map_entry_hap_id);
    }
    size_t actual_k = std::min<size_t>(k, stats.size());
    // use this if we want sorted
    std::partial_sort(
        stats.begin(), stats.begin() + static_cast<ptrdiff_t>(actual_k), stats.end(), std::greater<std::pair<double, size_t>>());
    // use this if we don't care about sorted
    // std::nth_element(stats.begin(), stats.begin() + actual_k, stats.end(),
    // std::greater<pair<double, size_t>>());

    // append to results
    results.emplace_back(window_start_site, window_end_site, std::vector<std::pair<size_t, double>>());
    for (size_t stats_idx = 0; stats_idx < actual_k; ++stats_idx) {
      std::get<2>(results[results.size() - 1])
          .emplace_back(stats[stats_idx].second, stats[stats_idx].first);
    }
  }

  return results;
}

std::ostream& operator<<(std::ostream& os, const HapData& data) {
  os << "HapData with " << data.num_haps << " haplotypes and " << data.num_sites;
  os << " sites, word size = " << data.word_size << " bits";
  return os;
}

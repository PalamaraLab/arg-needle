# This file is part of the ARG-Needle genealogical inference and
# analysis software suite.
# Copyright (C) 2023 ARG-Needle Developers.

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


"""ASMC decoding
"""

# General Python imports
from datetime import datetime
import logging
import math
import numpy as np
import os
import psutil; process = psutil.Process(os.getpid())
import shutil

# our packages
from asmc.asmc import DecodingParams, ASMC
import arg_needle_hashing_pybind as arg_needle_hashing
from .utils import btime

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


class ASMCDecoder(object):
    def __init__(self, site_positions, asmc_obj,
                 hasher=None, backup_hasher=None, asmc_pad_cm=100.0):
        self.site_positions = site_positions
        self.asmc_obj = asmc_obj
        if hasher is not None:
            self.hasher = hasher
        self.backup_hasher = backup_hasher
        self.cache_i = -1
        self.cache_map = None
        self.cache_mean = None
        self.expected_times = np.array(self.asmc_obj.get_expected_times())
        self.num_hash_queries = 0
        self.chunk_start = 0
        self.chunk_end = len(self.site_positions)
        self.asmc_pad_cm = asmc_pad_cm

    def __call__(self, i, j):
        """Decodes pair (i, j) between self.chunk_start and self.chunk_end.

        Performs batched queries and caches the results.
        """
        if j >= i:
            raise ValueError("Batching API expects first index > second index")

        if self.cache_i != i:
            self.asmc_obj.decode_pairs(
                [i for index in range(i)],
                [index for index in range(i)],
                self.chunk_start, self.chunk_end, self.asmc_pad_cm)
            results_struct = self.asmc_obj.get_ref_of_results()
            self.cache_mean = results_struct.per_pair_posterior_means
            self.cache_map = self.expected_times[results_struct.per_pair_MAPs]
            assert self.cache_mean.shape == (i, self.chunk_end - self.chunk_start)
            assert self.cache_map.shape == (i, self.chunk_end - self.chunk_start)
            self.cache_i = i

        return (self.site_positions,
                self.cache_map[j],
                self.cache_mean[j])

    def reset_chunked_upgma(self, start=None, end=None):
        self.cache_i = -1
        self.cache_map = None
        self.cache_mean = None
        if start is None:
            self.chunk_start = 0
        else:
            self.chunk_start = start
        if end is None:
            self.chunk_end = len(self.site_positions)
        else:
            self.chunk_end = end

    def compute_with_hashing(self, i, hash_topk, hash_cm, tmrca_mean=None,
                             tolerance=0, verbose=False, time_dict=None):
        """Given sample i, decodes for the hash_topk closest samples in each region

        See Supplementary Note 1 of our paper for more details.
        """
        if time_dict is None:
            time_dict = {"hash": 0, "asmc": 0, "smooth": 0, "thread": 0}
        indices = np.zeros(len(self.site_positions), dtype=np.int32)
        times = np.zeros(len(self.site_positions))

        def foo_func(x):
            time_dict["hash"] += x
        with btime(foo_func):
            closest_cousins = self.hasher.get_closest_cousins(i, hash_topk, tolerance, hash_cm)

            nothing_found = False
            low_score = False
            first_bad_id = 0
            for closest_cousin_id, entry in enumerate(closest_cousins):
                other_ids = [thing[0] for thing in entry[2]]
                # Heuristic 1: if a window finds nothing, then definitely refine
                if len(other_ids) == 0:
                    nothing_found = True
                    first_bad_id = closest_cousin_id
                else:
                    window_num_words = (entry[1] - entry[0]) // self.hasher.word_size + 1
                    sum_scores = sum([score_entry[1] for score_entry in entry[2]])
                    # Heuristic 2: if the sum of scores in a window is low, then also refine
                    if sum_scores < math.ceil(math.sqrt(hash_topk)) * window_num_words:
                        low_score = True
                        first_bad_id = closest_cousin_id

            if self.backup_hasher is not None and (nothing_found or low_score):
                if verbose:
                    logging.info("Refining hashing on ID {} using backup hasher".format(i))
                mid_bit = (closest_cousins[first_bad_id][0] + closest_cousins[first_bad_id][1]) // 2
                closest_cousins = self.backup_hasher.get_closest_cousins(i, hash_topk, tolerance, hash_cm)

                nothing_found = False
                for entry in closest_cousins:
                    other_ids = [thing[0] for thing in entry[2]]
                    if len(other_ids) == 0:
                        nothing_found = True

                if i > 2*hash_topk and nothing_found:
                    logging.warning("Warning: no cousins found for ID {}, even after refining. Consider decreasing the backup hash word size.".format(i))
            else:
                if i > 2*hash_topk:
                    if nothing_found:
                        logging.warning(
                            "Warning: no cousins found for ID {}. If you see this >10% of the time, decrease the hash word size.".format(i))
                    elif low_score:
                        logging.warning(
                            "Warning: few cousins found for ID {}. If you see this >10% of the time, decrease the hash word size.".format(i))

            if verbose:
                if self.num_hash_queries % 100 == 0:
                    logging.info("Selected hash results for {}:".format(i))
                    whole_length = len(closest_cousins)
                    logging.info("Length = {}, middle = {}".format(whole_length, whole_length // 2))
                    if whole_length > 2:
                        logging.info(closest_cousins[0])
                    if whole_length > 0:
                        logging.info(closest_cousins[whole_length // 2])
                    if whole_length > 2:
                        logging.info(closest_cousins[whole_length - 1])

            self.num_hash_queries += 1

        for entry in closest_cousins:
            def foo_func(x):
                time_dict["asmc"] += x
            with btime(foo_func):
                from_pos = entry[0]
                to_pos = entry[1] + 1
                other_ids = [thing[0] for thing in entry[2]]

                # TODO: be smarter by excluding already sampled IDs, and no replacement
                other_ids.extend(np.random.randint(0, i, size=(hash_topk - len(other_ids))))
                other_ids = np.array(other_ids, dtype=np.int32)
                assert len(other_ids) == hash_topk

                self.asmc_obj.decode_pairs(
                    [i for index in range(hash_topk)],
                    other_ids,
                    from_pos, to_pos, self.asmc_pad_cm)
                results_struct = self.asmc_obj.get_ref_of_results()
                batch_mean = results_struct.per_pair_posterior_means
                batch_map = self.expected_times[results_struct.per_pair_MAPs]
                assert batch_mean.shape == (hash_topk, to_pos - from_pos)
                assert batch_map.shape == (hash_topk, to_pos - from_pos)

            def foo_func(x):
                time_dict["smooth"] += x
            with btime(foo_func):
                tmrca_mean[other_ids, from_pos:to_pos] = batch_mean
                foo = np.argmin(batch_mean, axis=0)
                indices[from_pos:to_pos] = other_ids[foo]
                times[from_pos:to_pos] = batch_map[foo, range(to_pos - from_pos)]

        return indices, times


def make_asmc_decoder_simulation(
    simulation, base_tmp_dir, decoding_quant_file, recomb_rate=None, mapfile=None,
    asmc_tmp_string="asmc", hash_word_size=64, asmc_pad_cm=100.0,
    use_hashing=False, snp_ids=None, mode="array", verbose=False):

    assert mode in ["array", "sequence"]

    if simulation.num_samples % 2 != 0:
        raise ValueError("Even number of samples expected")

    # ASMC expects at least 300 samples by default to get an accurate CSFS
    # If we're less than 300, we augment by repeating data
    if simulation.num_samples < 300:
        augment_factor = math.ceil(300 / simulation.num_samples)
        if verbose:
            logging.info("Number of samples is less than 300, at {}, so augmenting by a factor of {}".format(
                simulation.num_samples, augment_factor))
    else:
        augment_factor = 1

    in_prefix = "in"
    time_string = datetime.now().strftime("day-%Y-%m-%d-time-%H:%M:%S.%f") # ensure unique temp directories!
    asmc_tmp_dir = os.path.join(base_tmp_dir, asmc_tmp_string + "_" + mode + "_" + time_string + "/")
    if verbose:
        logging.info("Using temporary directory " + asmc_tmp_dir)
    os.makedirs(asmc_tmp_dir, exist_ok=True)

    # Write .hap file
    out_path = os.path.join(asmc_tmp_dir, in_prefix + ".hap")
    site_positions = []
    with open(out_path, 'w') as out_file:
        last = -1
        snp_ids_index = 0
        for i, variant in enumerate(simulation.variants()):
            if snp_ids is None or snp_ids[snp_ids_index] == i:
                pos = int(variant.site.position)
                # in contrived cases, this might overflow the bounds of the genome
                if pos <= last:
                    pos = last + 1
                site_positions.append(pos)
                last = pos
                row_list = ['1', '.', str(pos), '0', '1']
                for i in range(augment_factor):
                    row_list += [str(entry) for entry in variant.genotypes]
                out_file.write(' '.join(row_list))
                out_file.write('\n')

                snp_ids_index += 1
                if snp_ids is not None and snp_ids_index == len(snp_ids):
                    break

    # Write .samples file
    site_positions = np.array(site_positions)  # this is good, don't change to np.float32
    out_path = os.path.join(asmc_tmp_dir, in_prefix + ".samples")
    with open(out_path, 'w') as out_file:
        out_file.write('\t'.join(["ID_1", "ID_2", "missing"]) + '\n')
        out_file.write('\t'.join(["0", "0", "0"]) + '\n')
        for i in range(augment_factor * simulation.num_samples // 2):
            out_file.write('\t'.join(["sample_" + str(i), "sample_" + str(i), "0"]) + '\n')

    # Write .map file
    if mapfile is not None:
        # In this case, read the cM information from the map file
        genetic_map_positions = []
        genetic_map_cms = []
        with open(mapfile, 'r') as fin:
            num_read = 0
            for line in fin:
                if num_read != 0:
                    chrom_string, position, rate, cm = line.strip('\n').split(' ')
                    genetic_map_positions.append(int(position))
                    genetic_map_cms.append(float(cm))
                num_read += 1
        genetic_map_positions = np.array(genetic_map_positions)
        genetic_map_cms = np.array(genetic_map_cms)

        # Linear interpolation
        indices = np.searchsorted(genetic_map_positions, site_positions, side="right")
        print(genetic_map_positions.shape)
        print(np.min(indices), np.max(indices))
        assert np.min(indices) > 0
        assert np.max(indices) < genetic_map_positions.shape[0]
        weights = (site_positions - genetic_map_positions[indices - 1]) / \
            (genetic_map_positions[indices] - genetic_map_positions[indices - 1])
        site_cms = weights * genetic_map_cms[indices] + (1-weights) * genetic_map_cms[indices - 1]

        # Write out result
        out_path = os.path.join(asmc_tmp_dir, in_prefix + ".map")
        with open(out_path, "w") as out_file:
            for i in range(site_positions.shape[0]):
                out_file.write('\t'.join([chrom_string, "SNP_" + str(site_positions[i]),
                    str(site_cms[i]), str(site_positions[i])]) + '\n')
    else:
        # In this case, just multiply by the rate
        assert recomb_rate is not None
        chrom_string = "chr"
        # Write out result
        out_path = os.path.join(asmc_tmp_dir, in_prefix + ".map")
        with open(out_path, "w") as out_file:
            for site_pos in site_positions:
                site_cm = (recomb_rate * 1e8) * site_pos / 1e6
                out_file.write('\t'.join([chrom_string, "SNP_" + str(site_pos),
                    str(site_cm), str(site_pos)]) + '\n')

    haps_file_root = os.path.join(asmc_tmp_dir, in_prefix)

    decoder_with_hasher = make_asmc_decoder(
        haps_file_root,
        decoding_quant_file,
        "" if mapfile is None else mapfile,
        mode,
        hash_word_size,
        backup_hash_word_size=0,
        asmc_pad_cm=asmc_pad_cm,
        use_hashing=use_hashing,
        verbose=verbose
        )
    shutil.rmtree(asmc_tmp_dir)
    return decoder_with_hasher


def make_asmc_decoder(
    haps_file_root, decoding_quant_file, mapfile="", mode="array",
    hash_word_size=64, backup_hash_word_size=0, asmc_pad_cm=100.0,
    use_hashing=False, verbose=False):

    # start to set up ASMC object
    noBatches = False
    hasher = None
    if use_hashing:
        if verbose:
            logging.info("Making HapData object")
        hasher = arg_needle_hashing.HapData(
            mode, haps_file_root, hash_word_size, mapfile, fill_sites=False)
        logging.info("Hashing data is {} by {}".format(hasher.num_haps, hasher.num_sites))

    backup_hasher = None
    if use_hashing and backup_hash_word_size > 0:
        if verbose:
            logging.info("Making backup HapData object")
        backup_hasher = arg_needle_hashing.HapData(
            mode, haps_file_root, backup_hash_word_size,
            map_file_path=mapfile, fill_sites=False)
        logging.info("Backup hashing data is {} by {}".format(hasher.num_haps, hasher.num_sites))

    if mode == "sequence":
        params = DecodingParams(
            in_file_root=haps_file_root,
            dq_file=decoding_quant_file,
            map_file=mapfile,
            compress=True, skip_CSFS_distance=float('nan'),
            decoding_mode_string="sequence", use_ancestral=False, no_batches=noBatches)
    elif mode == "array":
        params = DecodingParams(
            in_file_root=haps_file_root,
            dq_file=decoding_quant_file,
            map_file=mapfile,
            compress=False, skip_CSFS_distance=0,
            decoding_mode_string="array", use_ancestral=False, no_batches=noBatches)
    else:
        raise ValueError("Unrecognized mode, must be one of sequence or array")

    if verbose:
        logging.info("DecodingParams:")
        for field in ["foldData", "usingCSFS", "skipCSFSdistance", "compress",
                      "decodingMode", "useAncestral"]:
            logging.info("    {}: {}".format(field, getattr(params, field)))

    if verbose:
        logging.info("Memory: {}".format(process.memory_info().rss))
    logging.info("Making ASMCDecoder")
    asmc_obj = ASMC(params)
    asmc_obj.set_store_per_pair_posterior_mean()
    asmc_obj.set_store_per_pair_map()
    site_positions = np.array(asmc_obj.get_physical_positions())
    genetic_positions = np.array(asmc_obj.get_genetic_positions())
    if not np.all(np.diff(site_positions) > 0):
        raise ValueError("Physical positions from .map file must be strictly increasing")
    if not np.all(np.diff(genetic_positions) >= 0):
        raise ValueError("Genetic positions from .map file must not be decreasing")

    decoder = ASMCDecoder(site_positions, asmc_obj,
                          hasher=hasher, backup_hasher=backup_hasher,
                          asmc_pad_cm=asmc_pad_cm)
    return decoder

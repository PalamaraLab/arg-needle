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


"""ARG-Needle, ASMC-clust, and ARG normalization.
"""

# General Python imports
import logging
import math
import numpy as np
import os
import psutil; process = psutil.Process(os.getpid())

# Third-party packages
from fastcluster import linkage
import tskit

# our packages
import arg_needle_lib
from .decoders import make_asmc_decoder_simulation, make_asmc_decoder
from .simulator import Simulator # for ARG normalization
from .utils import btime, collect_garbage


logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

import arg_needle
DEFAULT_DECODING_QUANTITIES = os.path.join(
    os.path.dirname(arg_needle.__file__),
    "resources/30-100-2000_CEU.decodingQuantities.gz")


def add_default_arg_building_arguments(parser):
    # Number of samples
    parser.add_argument("--num_snp_samples", action="store", default=0, type=int,
        help="Number of haploid array samples at end of inference (default=0)")
    parser.add_argument("--num_sequence_samples", action="store", default=0, type=int,
        help="Number of haploid sequence samples at end of inference (default=0)")
    # ASMC and ASMC-clust parameters
    parser.add_argument("--asmc_clust", action="store", default=0, type=int,
        help="Whether to use ASMC-clust instead of ARG-Needle (default=0, nonzero means true)")
    parser.add_argument("--asmc_clust_chunk_sites", action="store", default=-1, type=int,
        help="Number of sites per chunk for memory-efficient ASMC-clust (default=-1 meaning not used)")
    parser.add_argument("--asmc_decoding_file", action="store", default=DEFAULT_DECODING_QUANTITIES,
        help=f"Where to find decoding quantities (default={DEFAULT_DECODING_QUANTITIES})")
    parser.add_argument("--asmc_pad_cm", action="store", default=2.0, type=float,
        help="ASMC padding in cM (default=2.0)")
    parser.add_argument("--asmc_tmp_string", action="store", default="asmc",
        help="String for prefixing ASMC temp files (default=asmc)")
    # Hashing parameters
    parser.add_argument("--sequence_hash_cm", action="store", default=0.3, type=float,
        help="Sequence mode hashing window size in cM (default=0.3)")
    parser.add_argument("--snp_hash_cm", action="store", default=0.5, type=float,
        help="Array mode hashing window size in cM (default=0.5)")
    parser.add_argument("--hash_tolerance", action="store", default=1, type=int,
        help="Hashing algorithm tolerance (default=1)")
    parser.add_argument("--hash_topk", action="store", default=64, type=int,
        help="How many closest samples to return using hashing (default=64)")
    parser.add_argument("--hash_word_size", action="store", default=16, type=int,
        help="Hashing word size, must be between 1 and 64 (default=16)")
    parser.add_argument("--backup_hash_word_size", action="store", default=8, type=int,
        help="Backup hashing word size (must be between 0 and 64, 0 means no backup for real data inference, default=8)")

def check_hash_word_sizes(args):
    if args.hash_word_size > 64 or args.hash_word_size <= 0:
        raise ValueError("hash_word_size must be between 1 and 64")
    if args.backup_hash_word_size > 64 or args.backup_hash_word_size < 0:
        raise ValueError("backup_hash_word_size must be between 0 and 64")


# To use this function, must also define args.rho or args.mapfile before passing in
# For constant recombination rate, set args.mapfile to None
def build_arg_simulation(args, simulation, base_tmp_dir, snp_indices=None,
                         mode="array", verbose=False, time_dict=None):
    check_hash_word_sizes(args)
    assert isinstance(base_tmp_dir, str)
    assert args.num_snp_samples > 0 or args.num_sequence_samples > 0
    use_asmc_clust = (args.asmc_clust != 0)
    if not mode in ["array", "sequence", "both"]:
        raise ValueError('mode must be one of "array", "sequence", or "both"')
    max_memory = 0

    arg = None
    if mode == "array":
        if args.num_snp_samples == 0:
            raise ValueError(f'Expected num_snp_samples > 0 for mode == "{mode}"')
        # Handle the SNP samples
        use_hashing = (not use_asmc_clust) and (args.hash_topk > 0)
        pairwise_decoder = make_asmc_decoder_simulation(
            simulation, base_tmp_dir, args.asmc_decoding_file, args.rho, args.mapfile,
            args.asmc_tmp_string, args.hash_word_size, args.asmc_pad_cm,
            use_hashing=use_hashing, snp_ids=snp_indices, mode="array", verbose=verbose)

        if use_asmc_clust:
            logging.info(f"Running ASMC-clust on {args.num_snp_samples} samples")
            arg = upgma_chunk(pairwise_decoder, args.num_snp_samples, simulation.sequence_length,
                              args.asmc_clust_chunk_sites, verbose)
        else:
            arg = thread_samples(arg, pairwise_decoder, args.num_snp_samples,
                                 args.hash_topk, args.snp_hash_cm,
                                 simulation.sequence_length,
                                 args.hash_tolerance, verbose, time_dict)
    else:
        if args.num_sequence_samples == 0:
            raise ValueError(f'Expected num_sequence_samples > 0 for mode == "{mode}"')

        # First handle the sequencing samples
        logging.info("Subsampling to only keep the sequence samples")
        sequence_sample_ids = [i for i in range(args.num_sequence_samples)]
        simulation_sequence = simulation.simplify(samples=sequence_sample_ids)
        logging.info("Done subsampling")

        use_hashing = (not use_asmc_clust) and (args.hash_topk > 0)
        pairwise_decoder = make_asmc_decoder_simulation(
            simulation_sequence, base_tmp_dir, args.asmc_decoding_file, args.rho, args.mapfile,
            args.asmc_tmp_string, args.hash_word_size, args.asmc_pad_cm,
            use_hashing=use_hashing, snp_ids=None, mode="sequence", verbose=verbose)

        if use_asmc_clust:
            logging.info(f"Running ASMC-clust on {args.num_sequence_samples} samples")
            arg = upgma_chunk(pairwise_decoder, args.num_sequence_samples,
                              simulation_sequence.sequence_length, args.asmc_clust_chunk_sites, verbose)
        else:
            arg = thread_samples(arg, pairwise_decoder, args.num_sequence_samples,
                                 args.hash_topk, args.sequence_hash_cm,
                                 simulation_sequence.sequence_length,
                                 args.hash_tolerance, verbose, time_dict)

        # Now thread the SNP samples (not possible to add using ASMC-clust)
        if mode == "both":
            if args.num_snp_samples == 0:
                raise ValueError(f'Expected num_snp_samples > 0 for mode == "{mode}"')

            arg_needle_lib.serialize_arg(
                arg, f'{base_tmp_dir}/{args.asmc_tmp_string}.step1.argn')
            arg = arg_needle_lib.deserialize_arg(
                f'{base_tmp_dir}/{args.asmc_tmp_string}.step1.argn',
                reserved_samples=args.num_sequence_samples + args.num_snp_samples)

            max_memory = max(max_memory, process.memory_info().rss)
            del pairwise_decoder
            if verbose:
                logging.info("Deleted decoder (and possible hasher)")
                collect_garbage()
 
            use_hashing = args.hash_topk > 0

            pairwise_decoder = make_asmc_decoder_simulation(
                simulation, base_tmp_dir, args.asmc_decoding_file, args.rho, args.mapfile,
                args.asmc_tmp_string, args.hash_word_size, args.asmc_pad_cm,
                use_hashing=use_hashing, snp_ids=snp_indices, mode="array", verbose=verbose)

            arg = thread_samples(arg, pairwise_decoder, args.num_snp_samples,
                                 args.hash_topk, args.snp_hash_cm,
                                 simulation.sequence_length,
                                 args.hash_tolerance, verbose, time_dict)

    max_memory = max(max_memory, process.memory_info().rss)
    del pairwise_decoder
    if verbose:
        logging.info("Deleted decoder (and possible hasher)")
        collect_garbage()

    assert isinstance(arg, arg_needle_lib.ARG)
    if arg.num_samples() != simulation.num_samples:
        logging.warning("Warning: number of samples in inferred ARG does not match simulation")
    logging.info("Done with inference")

    return arg, max_memory


def build_arg(args, haps_file_root, map_file, mode="array",
              verbose=False, time_dict=None):
    check_hash_word_sizes(args)
    if not mode in ["array", "sequence"]:
        raise ValueError('mode must be one of "array" or "sequence"')
    logging.info(f"Using {mode} mode")
    use_asmc_clust = (args.asmc_clust != 0)

    if mode == "array":
        if args.num_snp_samples == 0:
            raise ValueError(f'Expected num_snp_samples > 0 for mode == "{mode}"')
        num_samples = args.num_snp_samples
        hash_cm = args.snp_hash_cm
    else:
        if args.num_sequence_samples == 0:
            raise ValueError(f'Expected num_sequence_samples > 0 for mode == "{mode}"')
        num_samples = args.num_sequence_samples
        hash_cm = args.sequence_hash_cm

    arg = None
    max_memory = 0

    # Create ASMC decoder
    use_hashing = (not use_asmc_clust) and (args.hash_topk > 0)
    pairwise_decoder = make_asmc_decoder(
        haps_file_root, args.asmc_decoding_file, map_file,
        mode=mode, hash_word_size=args.hash_word_size,
        backup_hash_word_size=args.backup_hash_word_size,
        asmc_pad_cm=args.asmc_pad_cm, use_hashing=use_hashing,
        verbose=verbose)

    # using pairwise_decoder positions, figure out good arg_start and arg_end parameters
    first_pos = pairwise_decoder.site_positions[0]
    last_pos = pairwise_decoder.site_positions[-1]
    mean_spacing = (last_pos - first_pos) / (len(pairwise_decoder.site_positions) - 1)
    assert mean_spacing > 2
    arg_start = math.floor(first_pos - mean_spacing / 2)
    if arg_start < 0:
        logging.warning("Warning: predicted start position is negative, rewriting to 0")
        arg_start = 0
    arg_end = math.ceil(last_pos + mean_spacing / 2)

    if use_asmc_clust:
        logging.info("Running ASMC-clust on {} samples".format(num_samples))
        arg = upgma_chunk(pairwise_decoder, num_samples, (arg_start, arg_end),
                          args.asmc_clust_chunk_sites, verbose)
    else:
        arg = thread_samples(arg, pairwise_decoder, num_samples,
                             args.hash_topk, hash_cm,
                             (arg_start, arg_end),
                             args.hash_tolerance, verbose, time_dict)

    max_memory = max(max_memory, process.memory_info().rss)
    del pairwise_decoder
    if verbose:
        logging.info("Deleted decoder (and possible hasher)")
        collect_garbage()

    assert arg.num_samples() == num_samples
    logging.info("Done with ARG building")

    return arg, max_memory


def extend_arg(args, initial_arg, haps_file_root, map_file,
               mode="array", verbose=False, time_dict=None):
    check_hash_word_sizes(args)
    if not mode in ["array", "sequence"]:
        raise ValueError('mode must be one of "array" or "sequence"')
    logging.info(f"Using {mode} mode")
    if (args.asmc_clust != 0):
        raise ValueError("Cannot use ASMC-clust to extend an ARG, set asmc_clust = 0")

    if mode == "array":
        if args.num_snp_samples == 0:
            raise ValueError(f'Expected num_snp_samples > 0 for mode == "{mode}"')
        num_samples = args.num_snp_samples
        hash_cm = args.snp_hash_cm
    else:
        if args.num_sequence_samples == 0:
            raise ValueError(f'Expected num_sequence_samples > 0 for mode == "{mode}"')
        num_samples = args.num_sequence_samples
        hash_cm = args.sequence_hash_cm

    max_memory = 0

    # Do threading for the SNP samples
    if initial_arg is None:
        initial_arg_num_samples = 0
    else:
        initial_arg_num_samples = initial_arg.num_samples()
    if num_samples == initial_arg_num_samples:
        logging.info("Initial ARG already has {} samples, returning".format(initial_arg_num_samples))
        return initial_arg, max_memory
    if num_samples < initial_arg_num_samples:
        raise ValueError(
            "Trying to thread to {} haploid samples, but initial ARG already has {} haploid samples".format(
            num_samples, initial_arg_num_samples))
    num_thread_samples = num_samples - initial_arg_num_samples

    # Check that we have enough reserved samples for how much we want to thread
    if initial_arg is not None:
        if initial_arg.reserved_samples < num_samples:
            raise ValueError(
                "Please set reserved_samples to at least {}".format(num_samples))

    # Create ASMC decoder
    use_hashing = args.hash_topk > 0
    pairwise_decoder = make_asmc_decoder(
        haps_file_root, args.asmc_decoding_file, map_file,
        mode=mode, hash_word_size=args.hash_word_size,
        backup_hash_word_size=args.backup_hash_word_size,
        asmc_pad_cm=args.asmc_pad_cm, use_hashing=use_hashing,
        verbose=verbose)

    # using pairwise_decoder positions, figure out good arg_start and arg_end parameters
    first_pos = pairwise_decoder.site_positions[0]
    last_pos = pairwise_decoder.site_positions[-1]
    mean_spacing = (last_pos - first_pos) / (len(pairwise_decoder.site_positions) - 1)
    assert mean_spacing > 2
    arg_start = math.floor(first_pos - mean_spacing / 2)
    if arg_start < 0:
        logging.warning("Warning: predicted start position is negative, rewriting to 0")
        arg_start = 0
    arg_end = math.ceil(last_pos + mean_spacing / 2)

    arg = thread_samples(initial_arg, pairwise_decoder, num_thread_samples,
                         args.hash_topk, hash_cm,
                         (arg_start, arg_end),
                         args.hash_tolerance, verbose, time_dict)

    max_memory = max(max_memory, process.memory_info().rss)
    del pairwise_decoder
    if verbose:
        logging.info("Deleted decoder (and possible hasher)")
        collect_garbage()

    assert arg.num_samples() == num_samples
    logging.info("Done with ARG building")

    return arg, max_memory


def thread_samples(initial_arg, pairwise_decoder, num_next_samples,
                   hash_topk=0, hash_cm=0.1, arg_bounds=0,
                   tolerance=0, verbose=False, time_dict=None):
    """Takes an existing `arg_needle_lib.ARG` and threads additional samples

    Arguments:
        initial_arg: the initial `arg_needle_lib.ARG`, or `None`.
        pairwise_decoder: a function that takes two IDs and returns a tuple of
            sites, MAP, and posterior mean (all `numpy` arrays).
        num_next_samples: how many more samples to thread. If 0, no action is
            taken.
        hash_topk: how many closest samples to return using hashing
        hash_cm: hashing window size in cM
        arg_bounds: used to initialize with initial_arg = None, otherwise unused.
            If it's a tuple, it means both start and end. Otherwise, it just
            means end. Start must be a nonnegative int, otherwise an error is
            raised.
        tolerance: tolerance parameter
        verbose: whether to verbosely log information
        time_dict: None or a dict with keys "hash", "asmc", "smooth", "thread"
            for timing the respective steps
    """
    if time_dict is None:
        time_dict = {"hash": 0, "asmc": 0, "smooth": 0, "thread": 0}

    if num_next_samples == 0:
        return initial_arg

    if initial_arg is None:
        num_end_samples = num_next_samples
        start_thread_id = 0

        # Process arg_bounds
        if isinstance(arg_bounds, int) or isinstance(arg_bounds, float):
            arg_start = 0
            arg_length = arg_bounds
        else:
            assert isinstance(arg_bounds, tuple) or isinstance(arg_bounds, list)
            assert len(arg_bounds) == 2
            arg_start = arg_bounds[0]
            arg_length = arg_bounds[1] - arg_bounds[0]
        if not (isinstance(arg_start, int) and arg_start >= 0):
            raise ValueError("Start of ARG must be given as an int")
        if arg_length <= 0:
            raise ValueError("ARG length must be positive")
        if arg_start > pairwise_decoder.site_positions[0]:
            raise ValueError("Start of ARG must not exceed first site position")
        if arg_start + arg_length < pairwise_decoder.site_positions[-1]:
            raise ValueError("Start of ARG must not be before last site position")

        arg = arg_needle_lib.ARG(0, arg_length, reserved_samples=num_next_samples)
        arg.set_offset(arg_start)
    else:
        if not isinstance(initial_arg, arg_needle_lib.ARG):
            raise ValueError("Expecting an arg_needle_lib.ARG")
        num_end_samples = num_next_samples + initial_arg.num_samples()
        start_thread_id = initial_arg.num_samples()
        arg = initial_arg

    haploid_sample_size = pairwise_decoder.asmc_obj.get_haploid_sample_size()
    if num_end_samples > haploid_sample_size:
        raise ValueError(
            "Trying to thread to {} haploid samples, but data only has {} haploid samples".format(
                num_end_samples, haploid_sample_size))

    # Read ARG offset if it exists. Use this instead of arg_start.
    offset = 0
    if hasattr(arg, "offset"):
        offset = arg.offset

    logging.info("About to thread {} samples".format(num_next_samples))
    if hash_topk > 0:
        assert hasattr(pairwise_decoder, 'hasher')
        logging.info("Hashing parameters:")
        logging.info("  K for top K: {}".format(hash_topk))
        logging.info("  Word size: {}".format(pairwise_decoder.hasher.word_size))
        logging.info("  Window cM size: {}".format(hash_cm))
        logging.info("  Tolerance: {}".format(tolerance))

    if hash_topk > 0:
        hasher = pairwise_decoder.hasher
        id_set = hasher.hashed_hap_ids
        for i in range(start_thread_id):
            if i not in id_set:
                hasher.add_to_hash(i)
                if pairwise_decoder.backup_hasher is not None:
                    pairwise_decoder.backup_hasher.add_to_hash(i)
        assert len(hasher.hashed_hap_ids) == start_thread_id

    posterior_phys_pos = None
    for i in range(start_thread_id, start_thread_id + num_next_samples):
        arg.add_sample()
        if i == 0:
            if hash_topk != 0:
                hasher.add_to_hash(0)
                if pairwise_decoder.backup_hasher is not None:
                    pairwise_decoder.backup_hasher.add_to_hash(0)
            continue
        if i == 1 or i % 100 == 0 or i == start_thread_id + num_next_samples - 1:
            logging.info("Threading sample {}".format(i))
            if verbose:
                logging.info("Memory: {}".format(process.memory_info().rss))
        if verbose:
            if i % 10 == 0:
                logging.info(i)
        if hash_topk == 0:
            tmrca_map = None
            tmrca_mean = None
            for j in range(0, i):
                res = pairwise_decoder(i, j)  # a tuple containing the return values
                if posterior_phys_pos is None:
                    posterior_phys_pos = res[0]
                    posterior_phys_pos -= offset
                    # from posterior_phys_pos, compute midpoints which are fed into
                    # the threading procedure
                    if True:
                        threading_midpoints = posterior_phys_pos * 0.5
                        threading_midpoints[1:] += posterior_phys_pos[:-1] * 0.5
                        threading_midpoints[0] = 0
                    else:
                        threading_midpoints = posterior_phys_pos
                        threading_midpoints[0] = 0
                if tmrca_map is None:
                    tmrca_map = np.zeros((i, len(posterior_phys_pos)), dtype=np.float32)
                    tmrca_mean = np.zeros((i, len(posterior_phys_pos)), dtype=np.float32)
                tmrca_map[j] = res[1]
                tmrca_mean[j] = res[2]

            # Set up averaging of posterior mean using min MAP to determine boundaries
            indices = np.argmin(tmrca_mean, axis=0)
            times = tmrca_map[indices, range(len(posterior_phys_pos))]
        else:
            if posterior_phys_pos is None:
                posterior_phys_pos = pairwise_decoder.site_positions
                posterior_phys_pos -= offset
                # from posterior_phys_pos, compute midpoints which are fed into
                # the threading procedure
                if True:
                    threading_midpoints = posterior_phys_pos * 0.5
                    threading_midpoints[1:] += posterior_phys_pos[:-1] * 0.5
                    threading_midpoints[0] = 0
                else:
                    threading_midpoints = posterior_phys_pos
                    threading_midpoints[0] = 0

            if i == start_thread_id or i == 1:
                tmrca_mean = np.full((start_thread_id + num_next_samples - 1, len(posterior_phys_pos)), np.nan)

            # modifies tmrca_mean
            indices, times = pairwise_decoder.compute_with_hashing(
                i, hash_topk, hash_cm, tmrca_mean, tolerance, verbose, time_dict)
            def foo_func(x):
                time_dict["hash"] += x
            with btime(foo_func):
                pairwise_decoder.hasher.add_to_hash(i)
                if pairwise_decoder.backup_hasher is not None:
                    pairwise_decoder.backup_hasher.add_to_hash(i)

        def foo_func(x):
            time_dict["smooth"] += x
        with btime(foo_func):
            index_breaks = np.where(np.diff(indices) != 0)[0]
            time_breaks = np.where(np.diff(times) != 0)[0]

            # Two pointer method is O(N), but doesn't use numpy. This is O(N log N)
            c = np.unique(np.concatenate(([0], index_breaks+1, time_breaks+1)))
            if verbose:
                if i % 100 == 1:
                    if len(c) == 1:
                        logging.info("Threading {}, 1 smooth interval of size {} based on time {}".format(
                            i, len(posterior_phys_pos), times[0]))
                    else:
                        max_diff_idx = np.argmax(np.diff(c))
                        logging.info("Threading {}, {} smooth intervals, maximum size {} from {} to {} based on time {}".format(
                            i,
                            c.shape[0], np.max(np.diff(c)), c[max_diff_idx], c[max_diff_idx + 1],
                            times[c[max_diff_idx]]))

            # Averaging of posterior mean using min MAP to determine boundaries
            mean_times = np.zeros(c.shape[0], dtype=np.float32)
            num_nans = 0
            # could possibly replace for loop with numpy
            for k in range(c.shape[0]):
                begin = c[k]
                if k == c.shape[0] - 1:
                    end = None
                else:
                    end = c[k+1]
                mean_times[k] = np.mean(tmrca_mean[indices[begin], begin:end])
                if np.isnan(mean_times[k]):
                    num_nans += 1
            times_to_use = mean_times
            if num_nans > 0:
                logging.error("Got {} NaNs out of {} time values when threading sample {}".format(
                    num_nans, c.shape[0], i))
                raise ValueError("NaNs in threading times")

        def foo_func(x):
            time_dict["thread"] += x
        with btime(foo_func):
            # add some noise for good measure to ensure unique times when threading
            arg.thread_sample(
                threading_midpoints[c],
                indices[c],
                times_to_use * (1 + 1e-6*np.random.randn(c.shape[0])))

    return arg


def upgma_chunk(pairwise_decoder, num_samples, arg_bounds, num_chunk_sites, verbose=False):
    n = num_samples
    num_pairs = n * (n - 1) // 2

    # Process arg_bounds
    if isinstance(arg_bounds, int) or isinstance(arg_bounds, float):
        arg_start = 0
        arg_length = arg_bounds
    else:
        assert isinstance(arg_bounds, tuple) or isinstance(arg_bounds, list)
        assert len(arg_bounds) == 2
        arg_start = arg_bounds[0]
        arg_length = arg_bounds[1] - arg_bounds[0]
    if not (isinstance(arg_start, int) and arg_start >= 0):
        raise ValueError("Start of ARG must be given as an int")
    if arg_length <= 0:
        raise ValueError("ARG length must be positive")
    if arg_start > pairwise_decoder.site_positions[0]:
        raise ValueError("Start of ARG must not exceed first site position")
    if arg_start + arg_length < pairwise_decoder.site_positions[-1]:
        raise ValueError("Start of ARG must not be before last site position")

    tables = tskit.TableCollection(sequence_length=arg_length)
    for i in range(n):
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0.0)

    num_sites = len(pairwise_decoder.site_positions)
    if num_chunk_sites < 0:
        assert num_chunk_sites == -1
        num_chunk_sites = num_sites
    num_chunks = math.ceil(num_sites / num_chunk_sites)

    boundaries = np.zeros(num_sites + 1) # this was previously dtype=np.float32 but that led to errors
    midpoints = (pairwise_decoder.site_positions[:-1] + pairwise_decoder.site_positions[1:]) / 2
    boundaries[1:-1] = midpoints - arg_start
    boundaries[-1] = arg_length

    start_site = 0
    node_offset = 0
    distance_matrix_unrolled = None
    for k in range(num_chunks):
        end_site = min(num_sites, start_site + num_chunk_sites)
        pairwise_decoder.reset_chunked_upgma(start_site, end_site)
        logging.info("Chunk {} of {} (1-indexed), site {} to {}".format(k+1, num_chunks, start_site, end_site))
        del distance_matrix_unrolled
        if verbose:
            collect_garbage()
        distance_matrix_unrolled = np.ones((num_pairs, end_site - start_site), dtype=np.float32)
        for i in range(num_samples):
            if verbose:
                if i % 100 == 0 or i == num_samples - 1 or i == 10:
                    logging.info("Computing distance against sample {}".format(i))
                    logging.info("Memory: {}".format(process.memory_info().rss))
            index = i - 1
            for j in range(0, i):
                res = pairwise_decoder(i, j)  # a tuple containing the return values
                distance_matrix_unrolled[index, :] = res[2]
                # We need to fill this in along the columns
                index += (n - j - 2)
        if verbose:
            logging.info("Memory: {}".format(process.memory_info().rss))

        for j in range(end_site - start_site):
            if verbose:
                if j % 1000 == 0:
                    logging.info(j)
            clustering = linkage(distance_matrix_unrolled[:, j], method="average")
            last_distance = 0
            for i in range(clustering.shape[0]):
                row = clustering[i]
                height_to_use = row[2]
                if height_to_use == 0:
                    height_to_use = 0.001
                    raise RuntimeError("Bad: {} {} {}".format(row[0], row[1], row[2]))
                if height_to_use <= last_distance:
                    height_to_use = last_distance * 1.00001
                last_distance = height_to_use
                tables.nodes.add_row(flags=(not tskit.NODE_IS_SAMPLE), time=height_to_use)
                for child_id in [int(row[0]), int(row[1])]:
                    if child_id < n:
                        tables.edges.add_row(
                            left=boundaries[j+start_site], right=boundaries[j+1+start_site],
                            parent=node_offset + n + i, child=child_id)
                    else:
                        tables.edges.add_row(
                            left=boundaries[j+start_site], right=boundaries[j+1+start_site],
                            parent=node_offset + n + i, child=node_offset + child_id)
            node_offset += clustering.shape[0]
        start_site += num_chunk_sites

    # Need to do this in case we try to use the decoder afterwards for threading
    pairwise_decoder.reset_chunked_upgma()

    tables.sort()
    ts_to_return = tables.tree_sequence()
    arg_to_return = arg_needle_lib.tskit_to_arg(ts_to_return)
    # Now set the offset to arg_start
    arg_to_return.set_offset(arg_start)
    if verbose:
        logging.info("Memory: {}".format(process.memory_info().rss))
    return arg_to_return


def normalize_ts(ts, demofile, num_seeds=1000, start_seed=1, assume_haploid=False, verbose=False):
    num_samples = ts.num_samples
    logging.info("Computing ARG normalization")
    if verbose:
        logging.info("{} edges, {} nodes, {} trees".format(ts.num_edges, ts.num_nodes, ts.num_trees))
        logging.info("Reading edge data")
    heights_to_span = {}
    # could rewrite this using numpy
    for edge in ts.edges():
        height = ts.node(edge.parent).time
        if height not in heights_to_span:
            heights_to_span[height] = 0
        heights_to_span[height] += edge.right - edge.left # edge.span in tskit 0.3.2

    if verbose:
        logging.info("Done reading edge data")
    numpy_stuff = np.array(sorted(heights_to_span.items()))
    if verbose:
        logging.info(np.sum(numpy_stuff[:, 1]))

    cumsum = np.cumsum(numpy_stuff[:, 1])
    quantiles = (cumsum - 0.5 * numpy_stuff[:, 1]) / cumsum[-1]
    del heights_to_span
    del cumsum

    # set up simulator with sequence + SNP samples
    mapfile = None
    simulator = Simulator(mapfile, demofile, sample_size=num_samples,
        mu=0, rho=0, assume_haploid=assume_haploid)
    logging.info(f"Running {num_seeds} msprime simulations for ARG normalization")
    if not verbose:
        # set a stricter threshold for msprime logging before running simulations
        logger = logging.getLogger("msprime")
        logger.setLevel(logging.WARNING)

    node_times = []
    highest_times = []
    for seed_offset in range(num_seeds):
        seed = seed_offset + start_seed
        np.random.seed(seed)

        if verbose:
            if seed_offset == 0 or seed % 10 == 0:
                logging.info("Starting simulation " + str(seed))
        # 1e6 doesn't matter here as rho = 0
        simulation = simulator.simulation(1e6, random_seed=seed)
        highest = 0
        for node in simulation.nodes():
            if node.time > 0:
                node_times.append(node.time)
                if node.time > highest:
                    highest = node.time
        highest_times.append(highest)

    highest_times.sort()
    node_times.sort()
    if verbose:
        collect_garbage() # collects garbage and logs memory usage

    logging.info(f"Applying ARG normalization")
    if verbose:
        logging.info("Largest and smallest node times out of {}: {} {}".format(
            len(node_times), node_times[:5], node_times[-5:]))
        logging.info("Highest times mean: {}".format(np.mean(highest_times)))
    del highest_times

    node_times = np.array(node_times)
    node_times_pad = np.zeros(len(node_times) + 2)
    node_times_pad[1:-1] = np.array(node_times)
    node_times_pad[0] = node_times[0] * 0  # 0 is tunable, as long as it's < 1
    node_times_pad[-1] = node_times[-1] * 1.05  # 1.05 is tunable, as long as it's > 1
    sim_quantiles = np.linspace(0, 1, len(node_times_pad))
    del node_times

    corrected = np.interp(quantiles, sim_quantiles, node_times_pad)
    correction_dict = dict(zip(numpy_stuff[:, 0], corrected))
    del quantiles
    del sim_quantiles
    del node_times_pad
    del numpy_stuff
    del corrected
    if verbose:
        logging.info("Done with interpolation")
        collect_garbage() # collects garbage and logs memory usage

    # Note: because we reuse existing data, metadata is preserved
    tables = ts.dump_tables()
    ts_times = tables.nodes.time
    # Could be possible to use np.interp directly
    new_ts_times = np.zeros(ts_times.shape)
    for i, ts_time in enumerate(ts_times):
        if ts_time != 0:
            new_ts_times[i] = correction_dict[ts_time]

    del correction_dict
    if verbose:
        logging.info("Done with dict lookup")
        collect_garbage() # collects garbage and logs memory usage

    tables.nodes.set_columns(time=new_ts_times, flags=tables.nodes.flags)
    tables.sort()
    if verbose:
        logging.info("Done with rewrite and sort")
        collect_garbage() # collects garbage and logs memory usage

    new_ts = tables.tree_sequence()
    if verbose:
        logging.info("New ARG has {} edges, {} nodes, {} trees".format(
            new_ts.num_edges, new_ts.num_nodes, new_ts.num_trees))
        collect_garbage() # collects garbage and logs memory usage
    return new_ts


def normalize_arg(arg, demofile, num_seeds=1000, start_seed=1, assume_haploid=False, verbose=False):
    """Apply ARG normalization to an ARG and return the result

    Arguments:
        arg: `arg_needle_lib.ARG` object
        demofile: demography file, e.g. from
          https://github.com/PalamaraLab/ASMC_data/tree/main/demographies
        num_seeds: number of msprime simulations to perform using the demography
        start_seed: starting seed
        assume_haploid: whether the demography is in haploid sample sizes (default=False)
        verbose: whether to log verbose output (default=False)
    """
    ts = arg_needle_lib.arg_to_tskit(arg)
    new_ts = normalize_ts(ts, demofile, num_seeds, start_seed, assume_haploid, verbose)
    del ts
    new_arg = arg_needle_lib.tskit_to_arg(new_ts)
    logging.info("Done with ARG normalization")
    return new_arg


def trim_arg(arg, start_position=-1, end_position=-1, verbose=False):
    """Trim an ARG and return the result

    Arguments:
        arg: `arg_needle_lib.ARG` object
        start_position: start position for trimming in terms of genome coordinates, not
          relative to the ARG. -1 is a special value that means no start trimming (default=-1).
        end_position: end position for trimming in terms of genome coordinates, not
          relative to the ARG. -1 is a special value that means no start trimming (default=-1).
    """
    if arg.start != 0:
        raise ValueError("For trimming to work, ARG start should be zero (offset can be nonzero)")

    # We use the naming convention that start and end are in [0, arg.end], and
    # start_position, end_position are in [arg.offset, arg.offset + arg.end]
    if start_position == -1:
        start = 0
    else:
        start = start_position - arg.offset

    if end_position == -1:
        end = arg.end
    else:
        end = end_position - arg.offset

    assert 0 <= start < end <= arg.end

    # this automatically rewrites the offset as needed
    new_arg = arg_needle_lib.trim_arg(arg, start, end)
    logging.info("Done with ARG trimming")
    return new_arg

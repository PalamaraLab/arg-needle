# This file is part of the ARG-Needle genealogical inference and
# analysis software suite.
# Copyright (C) 2023-2025 ARG-Needle Developers.

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


# This test assumes you have installed the dev dependencies of arg-needle.
# In the root of this repository, run:
#
# pip install .[dev]

import subprocess
import tempfile
import pathlib
import h5py
import numpy as np
import shutil

DATA_DIR = pathlib.Path(__file__).resolve().parent / "data"

ARGN_KEYS = ['arg_file_version', 'chromosome', 'datetime_created', 'end', 'mutations', 'node_bounds', 'num_edges',
             'num_mutations', 'num_nodes', 'offset', 'start', 'threaded_samples']

def check_attr_keys(attrs):
    """
    Check if the keys are as expected in generated HDF5 file
    """
    keys_in_generated_file = sorted([str(x) for x in attrs.keys()])

    if keys_in_generated_file != ARGN_KEYS:
        print(f"Expected the following keys:\n{ARGN_KEYS}\n but got:\n{keys_in_generated_file}")
        return False

    return True

def check_attr_vals(attrs):
    """
    Check attr values are correct
    """

    # Deterministic values:
    assert attrs["arg_file_version"] == 2
    assert attrs["chromosome"] == 1
    assert np.isclose(attrs["start"], 0.0, rtol=1e-8)
    assert np.isclose(attrs["end"], 2000079.0, rtol=1e-8)
    assert attrs["mutations"] == False
    assert attrs["node_bounds"] == True
    assert attrs["offset"] == 10001457
    assert attrs["threaded_samples"] == 400

    # These values were calculated by running the example about 100 times
    nodes_mean = 17203.69792
    nodes_std = 127.8394651
    edges_mean = 93654.98958
    edges_std = 591.7562353

    # This should almost never fail
    assert attrs["num_nodes"] > nodes_mean - 3.0 * nodes_std
    assert attrs["num_nodes"] < nodes_mean + 3.0 * nodes_std

    assert attrs["num_edges"] > edges_mean - 3.0 * edges_std
    assert attrs["num_edges"] < edges_mean + 3.0 * edges_std

    return True


def test_script_output():

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = pathlib.Path(tmpdir)

        # Run scripts inside temporary directory
        subprocess.run(["prepare_example"], cwd=tmp, check=True)
        subprocess.run(["infer_args", "--normalize", "0"], cwd=tmp, check=True)

        # Check output file exists
        outfile = tmp / "example.argn"
        assert outfile.exists()

        # Compare to frozen output
        with h5py.File(outfile, "r") as arg_file:
            assert check_attr_keys(arg_file.attrs)
            assert check_attr_vals(arg_file.attrs)


if __name__ == '__main__':
    test_script_output()

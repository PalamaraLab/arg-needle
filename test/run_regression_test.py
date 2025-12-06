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

DATA_DIR = pathlib.Path(__file__).resolve().parent / "data"

def compare_attr_keys(a, b):
    if set(a.keys()) != set(b.keys()):
        return False
    return True

def compare_attr_vals(a, b):

    errors = False

    unchecked_keys = 'datetime_created'
    deterministic_keys = ['arg_file_version', 'chromosome', 'end', 'mutations', 'node_bounds', 'offset', 'start', 'threaded_samples']
    nondeterministic_keys = ['num_edges', 'num_nodes']

    for k in deterministic_keys:
        if a[k] != b[k]:
            print(f"Key {k} values ({a[k]} and {b[k]}) do not match")
            errors = True

    for k in nondeterministic_keys:
        if not np.isclose(a[k], b[k], rtol=0.01):
            print(f"Key {k} values ({a[k]} and {b[k]}) differ by more than 1%")
            errors = True

    if errors:
        return False
    return True


def check_attrs(h5file):

    return True

def hdf5_equal(file1, file2):
    with h5py.File(file1, "r") as f1, h5py.File(file2, "r") as f2:
        # compare file-level attributes and root group
        if not compare_attr_keys(f1.attrs, f2.attrs):
            return False
        if not compare_attr_vals(f1.attrs, f2.attrs):
            return False

        print("hooray!")
        return True

def test_script_output():

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = pathlib.Path(tmpdir)

        # Run scripts inside temporary directory
        # subprocess.run(["prepare_example"], cwd=tmp, check=True)
        # subprocess.run(["infer_args", "--normalize", "0"], cwd=tmp, check=True)
        #
        # # Check output file exists
        # outfile = tmp / "example.argn"
        # assert outfile.exists()

        # Compare to frozen output


        hdf5_equal(DATA_DIR / "example.argn", DATA_DIR / "arg_from_example.argn")



if __name__ == '__main__':
    test_script_output()

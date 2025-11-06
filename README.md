# arg-needle

A Python package implementing ARG-Needle and ASMC-clust genealogical inference.

## Installation

Use
```
pip install arg_needle
```
to install from PyPI, and
```
pip install .
```
to install from a local clone.

## Documentation

Please see the [ARG-Needle manual](https://palamaralab.github.io/software/argneedle/) for all usage instructions and documentation.

## For developers: making a release

- Bump the version number in [setup.py](setup.py) and [CMakeLists.txt](CMakeLists.txt)
- Update [RELEASE_NOTES.md](RELEASE_NOTES.md)
- Push changes and check that all [GitHub workflows](https://github.com/PalamaraLab/arg_needle/actions) pass
- Tag the commit in Git using syntax `vX.Y.Z`
- Make a release on GitHub, which should trigger a new build that will upload Python wheels to PyPI

## Acknowledgements

`arg-needle` has been developed by Brian C. Zhang with support from Fergus Cooper, Árni Freyr Gunnarsson, and Pier Francesco Palamara.

The files `src/hashing/FileUtils.cpp` and `src/hashing/FileUtils.hpp` are adapted from analogous files in the Eagle software, which can be found at https://github.com/poruloh/Eagle/tree/master/src. Eagle was developed by Po-Ru Loh and released under the GNU General Public License v3.0 (GPLv3). Its license can be found under the `3rd_party/` directory.

## Reference

If you use this software, please cite:

B. C. Zhang, A. Biddanda, Á. F. Gunnarsson, F. Cooper, P. F. Palamara, Biobank-scale inference of ancestral recombination graphs enables genealogical analysis of complex traits. [Nature Genetics, 2023](https://www.nature.com/articles/s41588-023-01379-x).

# arg-needle

This repository contains arg-needle, which implements the ARG inference algorithms ARG-Needle and ASMC-clust.
Prebuilt CPython wheels are available for Linux (compatible with glibc ≥ 2.28) and macOS (built on macOS 15 for x86_64 and macOS 14 for arm64).

| Platform \ CPython          | ≤3.8 | 3.9 | 3.10 | 3.11 | 3.12 | 3.13 | 3.14 |
|-----------------------------|------|-----|------|------|------|------|------|
| Linux x86_64                | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| Linux aarch64               | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| macOS Intel (x86_64)        | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| macOS Apple Silicon (arm64) | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |

## Quickstart

### Install the Python module from PyPI

The Python module can be installed with:

```bash
pip install arg-needle
```

### Documentation

Please see the [ARG-Needle manual](https://palamaralab.github.io/software/argneedle/) for all usage instructions and documentation.

## License

arg-needle is distributed under the GNU General Public License v3.0 (GPLv3). For any questions or comments on arg-needle, please contact Pier Palamara using `<lastname>@stats.ox.ac.uk`.

## Acknowledgements

arg-needle has been developed by Brian C. Zhang with support from Fergus Cooper, Árni Freyr Gunnarsson, and Pier Francesco Palamara.

## Reference

If you use this software, please cite:

B. C. Zhang, A. Biddanda, Á. F. Gunnarsson, F. Cooper, P. F. Palamara, Biobank-scale inference of ancestral recombination graphs enables genealogical analysis of complex traits. [Nature Genetics, 2023](https://www.nature.com/articles/s41588-023-01379-x).

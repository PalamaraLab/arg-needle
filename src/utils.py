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

from contextlib import contextmanager
import gc
import logging
import os
import psutil; process = psutil.Process(os.getpid())
import sys
import time

# https://www.pythoncentral.io/measure-time-in-python-time-time-vs-time-clock/
if sys.platform == 'win32':
    # On Windows, the best timer is time.clock
    default_timer = time.clock
else:
    # On most other platforms the best timer is time.time
    default_timer = time.time
# from timeit import default_timer  # we don't want this as it disables gc


logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def btime_default_f(time_in_seconds):
    print("Time elapsed (seconds):", time_in_seconds)


@contextmanager
def btime(f=btime_default_f):
    """Timer.

    Times a block of code and applies function f to the resulting time in seconds.

    Inspired by https://stackoverflow.com/a/30024601.
    """
    start = default_timer()
    yield
    end = default_timer()
    time_in_seconds = end - start
    f(time_in_seconds)


def collect_garbage(sleep_seconds=1.):
    gc.collect()
    time.sleep(sleep_seconds)
    logging.info("Ran garbage collection")
    memory_in_bytes = process.memory_info().rss
    logging.info("Memory: {}".format(process.memory_info().rss))
    return memory_in_bytes


if __name__ == "__main__":
    with btime():
        time.sleep(1)
    # Nesting example inspired by https://github.com/hector-sab/ttictoc
    with btime():
        time.sleep(1)
        with btime():
            time.sleep(1)
        time.sleep(1)

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

import logging
import math
import numpy as np
import pandas

import msprime

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


class Simulator:
    def __init__(self, mapfile, demofile, sample_size, mu, rho=1.2e-8,
        assume_haploid=False, check=False):
        self.mu = mu
        self.sample_size = sample_size
        if mapfile is not None:
            logging.info("Reading " + mapfile, "...")
            self.recomb_map = msprime.RecombinationMap.read_hapmap(mapfile)
        else:
            self.recomb_map = None
        self.demography = self.read_demo(demofile, assume_haploid)
        self.rho = rho
        self.configure_demography()
        if check:
            self.check_demographics()

    def read_demo(self, demofile, assume_haploid=False):
        df = pandas.read_csv(demofile, sep="\t", header=None)
        df.columns = ['generation', 'size']
        if assume_haploid:
            df['size'] = df['size'] / 2
            logging.info("Read in demography assuming haploid")
        else:
            logging.info("Read in demography assuming diploid")
        return df

    def configure_demography(self):
        self.demographic_events = []
        self.pc = [msprime.PopulationConfiguration(self.sample_size)]
        for index in self.demography.index:
            if index == self.demography.shape[0] - 1: break
            now_time = self.demography['generation'][index]
            now_size = self.demography['size'][index]
            self.demographic_events.append(
                msprime.PopulationParametersChange(now_time, now_size, growth_rate=0))

    def check_demographics(self):
        dp = msprime.DemographyDebugger(
            population_configurations=self.pc, 
            demographic_events=self.demographic_events)
        dp.print_history()

    def simulation(self, length, random_seed=10, output=None):
        if self.recomb_map is None:
            tree_seq = msprime.simulate(
                population_configurations = self.pc, 
                demographic_events = self.demographic_events,
                mutation_rate = self.mu, 
                length=length,
                recombination_rate=self.rho,
                random_seed=random_seed
            )
        else:
            tree_seq = msprime.simulate(
                population_configurations = self.pc, 
                demographic_events = self.demographic_events,
                mutation_rate = self.mu, 
                recombination_map=self.recomb_map,
                random_seed=random_seed
            )
        if output != None:
            with open(output + ".vcf", "w") as vcf_file:
                tree_seq.write_vcf(vcf_file, 2)
            tree_seq.dump(output + ".tree")
        return tree_seq

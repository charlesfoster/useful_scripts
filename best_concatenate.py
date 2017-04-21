#!/usr/bin/env python2.7
# Usage: python best_concatenate.py concatenated_filename.nex
# This script will concatenate all nexus files within your directory
# All nexus files must end with '.nex'
# If taxa are missing any genes, the gene sequence will be filled with '?'

from os import getcwd, listdir
from os.path import isfile, join
from Bio.Nexus import Nexus
from sys import argv

mypath = getcwd()

nex_files = [f for f in listdir(mypath) if isfile(join(mypath, f)) if f.endswith('.nex')]

nexi =  [(fname, Nexus.Nexus(fname)) for fname in nex_files]

combined = Nexus.combine(nexi);
combined.write_nexus_data(filename=argv[1]);

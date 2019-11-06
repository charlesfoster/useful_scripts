#!/usr/bin/env python2.7
# Usage: ./deconcatenate.py concatenated_filename.nex
# This script will deconcatenate a nexus file
# Requires biopython (and its dependencies)
# Modified from a script found at https://gist.github.com/brantfaircloth

from Bio.Nexus import Nexus
from sys import argv

aln = Nexus.Nexus()
aln.read(argv[1])

# get count of charsets:
print "You have %d charsets" %len(aln.charsets.keys())

# take a gander at the charsets:
print "Your charsets are named %s" %aln.charsets.keys()

# split the concatenated file into charsets files (prepending whatever text you place after filename='')
aln.write_nexus_data_partitions(filename='my', charpartition=aln.charsets)


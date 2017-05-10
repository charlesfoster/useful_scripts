#!/usr/bin/env python2.7
# Usage: ./best_concatenate.py fasta_suffix
# e.g. ./best_concatenate.py fasta
# This script will convert all fasta files within your directory into nexus format
# This can be preferable to the converter.py script if you want non-interactive conversions, such as for usage in scripts
# Requires biopython (and its dependencies)

from os import getcwd, listdir, rename
from os.path import isfile, join
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from sys import argv

file_suffix = argv[1]

print "\n*****************\nHello, and welcome to my terminal. I will be your converter for today.\n"

print "Your files end with: %s\n" % file_suffix

mypath = getcwd()

print "Your path is: %s\n"  % mypath
 
fas_files = [f for f in listdir(mypath) if isfile(join(mypath, f)) if f.endswith(argv[1])] 

print "Your fasta files to be converted are: %s\n" % fas_files

print "I will convert them now. Thank me later.\n"

for file in fas_files:
     count = SeqIO.convert(file, "fasta", file+".nex", "nexus", generic_dna)
     print "Converting: %s" % file
     print"Converted %i taxa\n" % count

new_files = [f for f in listdir(mypath) if isfile(join(mypath, f)) if f.endswith(argv[1]+".nex")]

for file in new_files:
    rename(file, file.replace(argv[1]+".nex", "nex"))
     
print "\nI am finished. Enjoy your new files, my friend!"

#!/usr/bin/env python2.7

### If you use this script, please place your hands together and bow in thanks to Charles Foster from MEEP, the University of Sydney ###
### Note: the script relies on your file extensions accurately mirroring the file contents
### Acceptable extensions: .nex (nexus), .fasta .fas .fa (fasta), .phy (phylip)
### Usage: ./converter_CLI.py <original_file_type> <desired_file_type>
### Requires biopython (and its dependencies)

from os import getcwd, listdir, rename
from os.path import isfile, join
from sys import argv
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import generic_dna

print "\n*****************\n\nHello, and welcome to my terminal. I will be your converter for today.\n"


mypath = getcwd()

print "Your path is: %s\n"  % mypath

fas_files = [f for f in listdir(mypath) if isfile(join(mypath, f)) if f.endswith("fasta") or f.endswith("fas") or f.endswith("fa")] 
nex_files = [f for f in listdir(mypath) if isfile(join(mypath, f)) if f.endswith("nex") or f.endswith("nexus")] 
phy_files = [f for f in listdir(mypath) if isfile(join(mypath, f)) if f.endswith("phy")] 

all_files = fas_files + nex_files + phy_files

print "You have the following files:\n\nFasta: %s\nNexus: %s\nPhylip: %s\n" % (fas_files,nex_files,phy_files)

option = argv[1] + "_" + argv[2]

if option == "fasta_nexus":
    print "\nGood choice! I'll do it now...\n"
    for file in fas_files:
         count = SeqIO.convert(file, "fasta", file+".nex", "nexus", generic_dna)
         print "\nConverting: %s" % file
         print "Converted %i taxa\n" % count

elif option == "fasta_phy":
    print "\nGood choice! I'll do it now...\n"
    for file in fas_files:
         count = SeqIO.convert(file, "fasta", file+".phy", "phylip-relaxed")
         print "\nConverting: %s" % file
         print "Converted %i taxa\n" % count

elif option == "nexus_fasta":
    print "\nGood choice! I'll do it now...\n"
    for file in nex_files:
         count = SeqIO.convert(file, "nexus", file+".fas", "fasta")
         print "\nConverting: %s" % file
         print "Converted %i taxa\n" % count

elif option == "nexus_phy":
    print "\nGood choice! I'll do it now...\n"
    for file in nex_files:
         count = SeqIO.convert(file, "nexus", file+".phy", "phylip-relaxed")
         print "\nConverting: %s" % file
         print "Converted %i taxa\n" % count

elif option == "phylip_fasta":
    print "\nGood choice! I'll do it now...\n"
    for file in phy_files:
         count = SeqIO.convert(file, "phylip-relaxed", file+".fasta", "fasta")
         print "\nConverting: %s" % file
         print "Converted %i taxa\n" % count

elif option == "phylip_nexus":
    print "\nGood choice! However, I should warn you that any gene boundaries in a concatenated file will be lost. I'll do the conversion now...\n"
    for file in phy_files:
         count = SeqIO.convert(file, "phylip-relaxed", file+".nex", "nexus", generic_dna)
         print "\nConverting: %s" % file
         print "Converted %i taxa\n" % count

elif not option == "fasta_nexus" or option == "fasta_phy" or option == "nexus_fasta" or option == "nexus_phy" or option == "phylip_fasta" or option == "phylip_nexus":
    print "That was not an option. This is why we can't have nice things. I quit!\a\a\a"
    quit()

for file in os.listdir('.'):
    if file.endswith('.fas.nex'):
        rename(file, file.replace(".fas.nex", ".nex"))
    elif file.endswith('.fasta.nex'):
        rename(file, file.replace(".fasta.nex", ".nex"))
    elif file.endswith('.fa.nex'):
        rename(file, file.replace(".fa.nex", ".nex"))
    elif file.endswith('.phy.nex'):
        rename(file, file.replace(".phy.nex", ".nex"))
    elif file.endswith('.fas.phy'):
        rename(file, file.replace(".fas.phy", ".phy"))
    elif file.endswith('.fasta.phy'):
        rename(file, file.replace(".fasta.phy", ".phy"))
    elif file.endswith('.nex.phy'):
        rename(file, file.replace(".nex.phy", ".phy"))
    elif file.endswith('.nex.fasta'):
        rename(file, file.replace(".nex.fasta", ".fasta"))
    elif file.endswith('.nex.fas'):
        rename(file, file.replace(".nex.fas", ".fasta"))
    elif file.endswith('.phy.fasta'):
        rename(file, file.replace(".phy.fasta", ".fasta"))
 
print "\nI am now finished converting your files. Have a great day!\n"

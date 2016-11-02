from os import getcwd, walk, listdir, path, rename, system, chdir
from os.path import isfile, join
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from glob import iglob
from fileinput import input, filelineno
from shutil import move, copyfile, copy2


mypath = getcwd()
fas_files = [f for f in listdir(mypath) if isfile(join(mypath, f))]

for file in fas_files:
    count = SeqIO.convert(file, "fasta", file+".nex", "nexus", generic_dna)
    print"Converted %i records" % count

from os import getcwd, walk, listdir, path, rename, system, chdir
from os.path import isfile, join
from Bio.Nexus import Nexus
from sys import argv

mypath = getcwd()

nex_files = [f for f in listdir(mypath) if isfile(join(mypath, f)) if f.endswith('.nex')]

nexi =  [(fname, Nexus.Nexus(fname)) for fname in nex_files]

combined = Nexus.combine(nexi);
combined.write_nexus_data(filename=argv[1]);
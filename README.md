# useful_scripts
Throughout my PhD I have used/written many different scripts, most of which are very simple. I thought this might be a good place to store some of them. Most of these are amalgamations or modifications of existing scripts I've found, but some are novel. All of these have been useful at some stage or another, but I'll give a short outline of the most useful ones below.

converter.py: This interactive script will convert all files format within a directory of a specified format to another specified format. For example, it will convert all fasta files to nexus files. This is useful when one has many separate genes within a directory.

best_concatenate.py: This script will concatenate all nexus file within a directory. The resulting concatenated file adds question marks if a taxon is missing within a particular gene. Additionally, gene boundaries are maintained within a nexus block at the end of the file.

CLC_assembler.sh: This is a pipeline for genome assembly using the CLC cell assembly command line tools. It is specific to a relatively novel method we trialled in our lab group in which we pooled DNA samples from multiple untagged organisms. Hence, the script carries out a de novo assembly, followed by modifications, then extracts contigs specific to reference organisms that we specify using BLAST.

prune_nex.sh: This script prunes taxa from all nexus files in a directory. It keeps all taxa specified in a separate text file.

I will update this with descriptions of other scripts when I get time.

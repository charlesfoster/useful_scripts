# useful_scripts
Throughout my PhD and postdoctoral career I have used/written many different scripts, most of which are very simple. I thought this might be a good place to store some of them. Most of these are amalgamations or modifications of existing scripts I've found, but some are novel. All of these have been useful at some stage or another, but I'll give a short outline of some of the more useful ones below.

converter.py: This interactive script will convert all files format within a directory of a specified format to another specified format. For example, it will convert all fasta files to nexus files. This is useful when one has many separate genes within a directory.

best_concatenate.py: This script will concatenate all nexus file within a directory. The resulting concatenated file adds question marks if a taxon is missing within a particular gene. Additionally, gene boundaries are maintained within a nexus block at the end of the file.

DE_GO_pipeline.r: This pipeline will read in raw counts for each transcript within a transcriptome, as estimated using Salmon, summarise the counts at the gene level, carry out differential expression using DESeq2, then carry out GO and KEGG enrichment using GOseq. The required input files are detailed within the script. Note that many steps are specific to the workflow I employed, but reading through the file might give some ideas.

subset_fasta.sh: a simple bash script to subset a fasta file based on sequence identifiers in another text file (one per line). At its heart the script employs GNU grep, its main benefit is avoiding having to remember the correct flags every time. The fasta seqs need to be linearised, but for that you can use...

linearise_fasta.sh: a script employing awk to linearise fasta files.


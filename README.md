# useful_scripts
Throughout my PhD and postdoctoral career I have used/written many different scripts, most of which are very simple. I thought this might be a good place to store some of them. Most of these are amalgamations or modifications of existing scripts I've found, but some are novel. All of these have been useful at some stage or another, but I'll give a short outline of some of the more useful ones below.

converter.py: This interactive script will convert all files format within a directory of a specified format to another specified format. For example, it will convert all fasta files to nexus files. This is useful when one has many separate genes within a directory.

best_concatenate.py: This script will concatenate all nexus file within a directory. The resulting concatenated file adds question marks if a taxon is missing within a particular gene. Additionally, gene boundaries are maintained within a nexus block at the end of the file.

CLC_assembler.sh: This is a pipeline for genome assembly using the CLC cell assembly command line tools. It is specific to a relatively novel method we trialled in our lab group in which we pooled DNA samples from multiple untagged organisms. Hence, the script carries out a de novo assembly, followed by modifications, then extracts contigs specific to reference organisms that we specify using BLAST.

prune_nex.sh: This script prunes taxa from all nexus files in a directory. It keeps all taxa specified in a separate text file.

DE_GO_pipeline.r: This pipeline will read in raw counts for each transcript within a transcriptome, as estimated using Salmon, summarise the counts at the gene level, carry out differential expression using DESeq2, then carry out GO and KEGG enrichment using GOseq. The required input files are detailed within the script.

I will update this with descriptions of other scripts when I get time.

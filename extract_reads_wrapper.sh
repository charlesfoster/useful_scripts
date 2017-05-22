#!/bin/bash
# Script credit: Charles Foster, University of Sydney
# This shell script generates a .pbs script to run analyses on the Artemis HPC at the University of Sydney
# The aim of the analyses is to map Illumina reads against a reference genome, and then extract the mapped 
# reads into new fastq files for assembly. The unmapped reads are also stored in a separate .bam file, just
# in case you want to examine them further
#
# Usage: ./extract_reads_wrapper.sh ReferenceName SampleNumber
# e.g.: ./extract_reads_wrapper.sh Plant 4
#
# After running this wrapper, you will need to submit your job to the queue with qsub

echo "#!/bin/bash" >> extract_reads_$1.pbs
echo "#PBS -q defaultQ -P RDS-FSC-Cockroaches_PBH-RW -l select=1:ncpus=8:mem=12GB -M charles.foster@sydney.edu.au -m ae" >> extract_reads_$1.pbs
echo "#PBS -l walltime=23:00:00" >> extract_reads_$1.pbs
echo "" >> extract_reads_$1.pbs
echo "# Load the modules" >> extract_reads_$1.pbs
echo "module load bwa/0.7.12" >> extract_reads_$1.pbs
echo "module load samtools/1.3.1" >> extract_reads_$1.pbs
echo "" >> extract_reads_$1.pbs
echo "# Change to the directory where you will be running the analyses" >> extract_reads_$1.pbs
echo "cd $PWD" >> extract_reads_$1.pbs
echo "" >> extract_reads_$1.pbs
echo "bwa index -p $1Index -a is $1_reference.fasta" >> extract_reads_$1.pbs
echo "bwa mem -t 8 -k 19 -w 100 -d 100 -r 1.5 -c 10000 -A 1 -B 4 -O 6 -E 1 -L 5 -U 9 -T 30 -v 3 $1Index /project/RDS-FSC-Cockroaches_PBH-RW/Macrogen2017/DNA$2_$1_R1.fastq.gz /project/RDS-FSC-Cockroaches_PBH-RW/Macrogen2017/DNA$2_$1_R2.fastq.gz | \\" >> extract_reads_$1.pbs
echo "samtools view -bS - > DNA$2_$1_bwa.bam" >> extract_reads_$1.pbs
echo "samtools view -b -F 4 -@ 8 DNA$2_$1_bwa.bam > DNA$2_$1_mapped.bam" >> extract_reads_$1.pbs
echo "samtools view -b -f 4 -@ 8 DNA$2_$1_bwa.bam > DNA$2_$1_unmapped.bam" >> extract_reads_$1.pbs
echo "" >> extract_reads_$1.pbs
echo "# The *_mapped.bam file created in the previous step contains only the mapped reads. In this next" >> extract_reads_$1.pbs
echo "# step you will extract them into forward and reverse fastq files. These can then be assembled" >> extract_reads_$1.pbs
echo "# using CLC on our local computers. Note: we also created a *_unmapped.bam file  in the " >> extract_reads_$1.pbs
echo "# previous step, which, you guessed it, contains only the unmapped reads. You might want to " >> extract_reads_$1.pbs
echo "# check these further in case there are novel genes in there from your organism. Up to you." >> extract_reads_$1.pbs
echo "# Also, and singletons (unpaired reads) will be saved in the next step. This can easily be" >> extract_reads_$1.pbs
echo "# modified so we keep those separately." >> extract_reads_$1.pbs
echo "samtools fastq -1 DNA$2_$1_R1.fastq -2 DNA$2_$1_R2.fastq DNA$2_$1_mapped.bam" >> extract_reads_$1.pbs

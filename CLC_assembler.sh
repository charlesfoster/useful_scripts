#!/bin/bash

# This is a script to carry out a de novo assembly of pooled DNA samples from different 
# organisms with no barcodes, and then assemble the genome of each organism separately
# using reference genomes. In this case, it assembles separate genomes based on plant, 
# mantis shrimp, and cockroach reference genomes

# Required installations: clc_assembly_cell (command line version of CLC) with valid 
# licence, BLAST suite. These commands will need to be in your PATH, or the absolute paths 
# for these programs can be added below. Please also note the required names for the 
# reference genomes - these can be updated easily.

# To run: ./CLC_assembler.sh

# Assemble the data
echo "##### STEP 1: De novo assembly of the reads #####"
echo ""
echo "I will now do a de novo assembly for each sample."
for i in *_1.fq;
do
db=`echo $i | awk -F '_' '{print $1}'`
echo "I am assembling: $db"
REVERSE=`echo ${i%_1.fq}_2.fq`
echo "Your forward reads are: $i"
echo "Your reverse reads are: $REVERSE"
clc_assembler -o $db.assembled.fas -p fb ss 150 600 -q -i $i $REVERSE

# Map the reads to the original assembly
ASSEMBLED=`echo ${i%_1.fq}.assembled.fas`
echo "I will now map the reads from $db to this assembly: $ASSEMBLED"
clc_mapper -o $db.map1 -p fb ss 150 600 -q -i $i $REVERSE -d $ASSEMBLED -g 1 -e 1
echo "I am finished working with $db for now."
done

# Use the mapping file to correct any mistakes that occurred during the assembly
echo "Your reads have been mapped to the de novo assembly."
echo "I will now use the mapping file to correct any mistakes that occurred during the assembly."
for i in *.map1;
do
clc_extract_consensus -o ${i%.map1}.assembled2.fas -a $i -w
done

echo ""
echo "STEP 1 is now finished."
echo ""

echo "##### STEP 2: Running BLAST #####"
echo ""

# Move the first assemblies into a separate directory
mkdir PreliminaryAssemblies
mv *.assembled.fas PreliminaryAssemblies/

# Copy this assembly into relevant files to get just the contigs for each organism
echo "I will now make an assembly file for each organism."
for i in *.assembled2.fas;
do
cp $i Plant_$i
cp $i Shrimp_nuc_$i
cp $i Shrimp_mito_$i
cp $i Cockroach_nuc_$i
cp $i Cockroach_mito_$i
rm $i
done

# Make a BLAST database for each organism reference file
echo "I will create a BLAST reference database for each assembly."

for i in *.assembled2.fas;
do
makeblastdb -in $i -dbtype nucl -parse_seqids -out $i.blastdb.aaa
touch $i.blastdb.aaa
done

# To blast against the local reference database you created
echo "I will now blast each of the organisms against their correct references."
echo "### Plant Assemblies ###"
for i in Plant*.blastdb.aaa;
do
echo "I am blasting: $i"
blastn -query Plant_reference.fasta -db $i -out $i.txt -evalue 1e-170 -num_alignments 0
echo "Done! Unless you received an error message, your results are in: $i.txt"
done 

echo "### Shrimp Assemblies ###"
for i in Shrimp_nuc*.blastdb.aaa;
do
echo "I am blasting: $i"
blastn -query Shrimp_nuc_reference.fasta -db $i -out $i.txt -evalue 1e-170 -num_alignments 0
echo "Done! Unless you received an error message, your results are in: $i.txt"
done 

echo "### Shrimp Assemblies ###"
for i in Shrimp_mito*.blastdb.aaa;
do
echo "I am blasting: $i"
blastn -query Shrimp_mito_reference.fasta -db $i -out $i.txt -evalue 1e-170 -num_alignments 0
echo "Done! Unless you received an error message, your results are in: $i.txt"
done 

echo "### Cockroach Assemblies ###"
for i in Cockroach_nuc*.blastdb.aaa;
do
echo "I am blasting: $i"
blastn -query Cockroach_nuc_reference.fasta -db $i -out $i.txt -evalue 1e-170 -num_alignments 0
echo "Done! Unless you received an error message, your results are in: $i.txt"
done 

echo "### Cockroach Assemblies ###"
for i in Cockroach_mito*.blastdb.aaa;
do
echo "I am blasting: $i"
blastn -query Cockroach_mito_reference.fasta -db $i -out $i.txt -evalue 1e-170 -num_alignments 0
echo "Done! Unless you received an error message, your results are in: $i.txt"
done 

echo ""
echo "STEP 2 is now finished."
echo ""

echo "##### STEP 3: Retrieving sequences #####"
echo ""

# Retrieve sequences from the Blast results
echo "I will now retrieve your sequences for each organism from the BLAST results"
for i in *.fas.blastdb.aaa;
do
db=`echo $i | awk '{print $1}'`
RESULTS=`echo $i | sed 's/.fas.blastdb.aaa/.fas.blastdb.aaa.txt/'`
echo "Your database is: $db"
echo "Your BLAST results being used here are in: $RESULTS"
blastdbcmd -db $i -entry_batch $RESULTS -outfmt "%f %s" -out $RESULTS.fasta
done

for i in *.assembled2.fas.blastdb.aaa.txt.fasta;
do
mv $i ${i%.assembled2.fas.blastdb.aaa.txt.fasta}_final.fasta
echo "Your results have been saved to: ${i%.assembled2.fas.blastdb.aaa.txt.fasta}_final.fasta"
done

mkdir BlastFiles
mv *.assembled2.fas.blastdb.aaa* BlastFiles

mkdir SecondaryAssemblies
for i in *.assembled2.fas;
do
mv $i SecondaryAssemblies
done

echo ""
echo "STEP 3 is now finished."
echo ""

mkdir Plant_results
mv Plant*final.fasta Plant_results/
mkdir ShrimpMito_results
mv Shrimp_mito*final.fasta ShrimpMito_results/
mkdir ShrimpNuc_results
mv Shrimp_nuc*final.fasta ShrimpNuc_results/
mkdir CockroachMito_results
mv Cockroach_mito*final.fasta CockroachMito_results/
mkdir CockroachNuc_results
mv Cockroach_nuc*final.fasta CockroachNuc_results/

echo "Enjoy your final assemblies. You can find them in their respective results folders."
echo ""

#!/bin/bash
# This script batch renames taxa in all files with a specified extension, based on a text
# file detailing all names to be changed. The original names should be in the first column
# of the text file, and the new names should be in the second column.
#
# Usage: ./batch_rename_taxa.sh fasta names.txt

while read i; 
do
	WRONG=`echo ${i} | cut -f1 -d " "` 
	RIGHT=`echo ${i} | cut -f2 -d " "`
	
	for gene in *.$1
	do
		sed -i '' "s|${WRONG}|${RIGHT}|g" $gene
	done
done < $2

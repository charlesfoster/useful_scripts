#!/bin/bash
# Usage: sh rename_for_assembly.sh <outfile_name>

ls *_R1_norm* | sort > fastq_list_R1.txt
ls *_R2_norm* | sort > fastq_list_R2.txt

while read L
do 
	printf "\nPreparing new names for ${L}\n"
	V1=`echo $L | cut -f1 -d " "`
	V2=`echo $L | cut -f2 -d " "`
	V3=`grep $V1 fastq_list_R1.txt`
	V4=`grep $V1 fastq_list_R2.txt`
	NAME1=$(echo ${V3}_${V2})
	NAME2=$(echo ${V4}_${V2})
	echo $NAME1 >> new_R1_names.txt
	echo $NAME2 >> new_R2_names.txt
done < renaming.txt

printf "\nNow I'll rename...\n"
awk -F'_' 'BEGIN{OFS="_"}{print($1,$(NF),$(NF-2),$(NF-1))}' new_R1_names.txt | sort > tmp1 && mv tmp1 new_R1_names.txt
awk -F'_' 'BEGIN{OFS="_"}{print($1,$(NF),$(NF-2),$(NF-1))}' new_R2_names.txt | sort > tmp2 && mv tmp2 new_R2_names.txt

paste fastq_list_R1.txt new_R1_names.txt -d " " | sed "s|^|mv |g" | parallel -j1 -v
paste fastq_list_R2.txt new_R2_names.txt -d " " | sed "s|^|mv |g" | parallel -j1 -v

printf "\nNow I'll prepare your sample metadata file for Trinity...\n"

ls *_R1_norm.fq.gz | cut -f2 -d "_" > tmp_group.txt
ls *_R1_norm.fq.gz | cut -f1,2 -d "_" > tmp_repl.txt
ls *_R1_norm.fq.gz | sed "s|^|${PWD}/|g" > tmp_R1.txt
ls *_R2_norm.fq.gz | sed "s|^|${PWD}/|g" > tmp_R2.txt
paste tmp_group.txt tmp_repl.txt tmp_R1.txt tmp_R2.txt > $1
rm tmp_group.txt tmp_repl.txt tmp_R1.txt tmp_R2.txt

printf "\nAll done!\n"


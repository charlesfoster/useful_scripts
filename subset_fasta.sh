#!/bin/bash
# Author: Charles Foster
# Purpose: Prune a large fasta file to only contain the sequences named in a text file

l_flag=''
f_flag=''
h_flag=''
o_flag=''

print_usage() {
  printf "Usage: ./subset_fasta.sh -l <list_of_identifiers> -f <fasta_file> -o <outname.fa>\n"
  printf "Notes: 
  - all fasta sequences must be linearised (on one line each).
  - the script uses GNU grep, installed as ggrep on my Mac. Update the program name as required.\n"
}

if [[ $# -eq 0 ]] ; then
	print_usage 
    exit 1
fi

while getopts 'l:f:ho:' flag; do
  case "${flag}" in
    l) NAMES="${OPTARG}" ;;
    f) FASTA="${OPTARG}" ;;
    h) print_usage 
           exit 1 ;;
    o) OUTNAME="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

ggrep -w -A 1 -f ${NAMES} ${FASTA} --no-group-separator >> ${OUTNAME}

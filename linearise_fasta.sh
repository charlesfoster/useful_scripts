#!/bin/bash
# Author: Charles Foster, modified from an example found on Stack Overflow :)
# Purpose: linearise a fasta file

f_flag=''
h_flag=''
o_flag=''

print_usage() {
  printf "Usage: ./linearise_fasta.sh -f <fasta_file> -o <outname.fa>\n"
}

if [[ $# -eq 0 ]] ; then
	print_usage 
    exit 1
fi

while getopts 'f:ho:' flag; do
  case "${flag}" in
    f) FASTA="${OPTARG}" ;;
    h) print_usage 
           exit 1 ;;
    o) OUTNAME="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' ${FASTA} > ${OUTNAME}

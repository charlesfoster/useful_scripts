#!/bin/bash
### Prunes taxa from a nexus file by specifying a list of those to be kept.
### Usage: ./prune_nex.sh taxa_to_keep.txt
for i in *.nex;
do
sed '/matrix/q' $i >> ${i%.nex}_final.nex
while read j; do grep "^$j" $i >> ${i%.nex}_final.nex;done < keep.txt;
echo ";" >> ${i%.nex}_final.nex
echo "end;" >> ${i%.nex}_final.nex
KEEP_LEN=`wc -l keep.txt | sed "s|^ *||g"| sed "s| keep.txt||g"`
NTAX_ORIG=`grep -o "ntax=[0-9]*" $i`
sed -i '' "s|${NTAX_ORIG}|ntax=${KEEP_LEN}|g" ${i%.nex}_final.nex
done

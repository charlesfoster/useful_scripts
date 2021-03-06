#!/bin/bash
### Prunes taxa from a nexus file by specifying a list of those to be kept.
### Usage: ./prune_nex.sh taxa_to_keep.txt
### Caution: experimental. Check the output to make sure it looks OK.
for i in *.nex;
do
sed '/matrix/q' $i >> ${i%.nex}_pruned.nex
while read j; do grep "^$j" $i >> ${i%.nex}_pruned.nex;done < $1;
echo ";" >> ${i%.nex}_pruned.nex
echo "end;" >> ${i%.nex}_pruned.nex
KEEP_LEN=`wc -l $1 | sed "s|^ *||g"| sed "s| $1||g"`
NTAX_ORIG=`grep -o "ntax=[0-9]*" $i`
sed -i '' "s|${NTAX_ORIG}|ntax=${KEEP_LEN}|g" ${i%.nex}_pruned.nex
done

#!/bin/bash
### Prunes taxa from a phylip file by specifying a list of those to be kept.
### Usage: ./prune_phy.sh taxa_to_keep.txt
for i in *.phy;
do
head -n1 $i >> ${i%.phy}_pruned.phy
while read j; do grep "^$j" $i >> ${i%.phy}_pruned.phy;done < $1;
KEEP_LEN=`wc -l $1 | sed "s|^ *||g"| sed "s| $1||g"`
NTAX_ORIG=`egrep -m1 -o "^ *[0-9]+" $i`
sed -i '' "s|${NTAX_ORIG}|${KEEP_LEN}|g" ${i%.phy}_pruned.phy
done

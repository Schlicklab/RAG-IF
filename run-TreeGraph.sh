#!/bin/bash


prefix=$1

#echo ${prefix}


if [ -f ${prefix}.bpseq ]; then
   python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/modified-treeGraph/treeGraphs.py  ${prefix}.bpseq > ${prefix}.tg_log 
fi


if [ -f ${prefix}.nupack.bpseq ]; then
   python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/modified-treeGraph/treeGraphs.py  ${prefix}.nupack.bpseq > ${prefix}.nupack.tg_log 
fi


if [ -f ${prefix}.mfe.bpseq ]; then
   python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/modified-treeGraph/treeGraphs.py  ${prefix}.mfe.bpseq > ${prefix}.mfe.tg_log 
fi



if [ -f ${prefix}.centroid.bpseq ]; then
   python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/modified-treeGraph/treeGraphs.py  ${prefix}.centroid.bpseq > ${prefix}.centroid.tg_log 
fi





#if [ ! -f /tmp/foo.txt ]; then
#    echo "File not found!"
#fi




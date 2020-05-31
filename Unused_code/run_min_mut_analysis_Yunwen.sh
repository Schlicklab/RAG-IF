#!/bin/bash
# to run analysis to find minimal mutations 
#
#

for i in Top*;do

cd  $i
if [  -f *.survivors ]; then
    echo "in folder:" $i
    python /Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/end_seq_order.py $1  > min_mut.analysis4   
fi
cd ..	

done 



#!/bin/bash
# to run analysis to find minimal mutations 
#
#

for i in Top*;do

cd  $i
if [  -f *.survivors ]; then
    echo "in folder:" $i
    python /Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/end_seq_order_Swati.py $1  > min_mut.analysis4_Swati   
fi
cd ..	

done 



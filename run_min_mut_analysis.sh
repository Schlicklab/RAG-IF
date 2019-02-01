#!/bin/bash
# to run analysis to find minimal mutations 
#
#

for i in Top*;do

cd  $i
if [  -f *.survivors ]; then
    echo "in folder:" $i
    python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/products-7_4/InitialDesign/IIIb/7_4-gaif-run1/6/Top29/result_analysis/end_seq_order.py 8_6  > min_mut.analysis4 &   
fi
cd ..	

done 



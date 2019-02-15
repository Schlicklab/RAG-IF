#!/bin/bash

cd $1 

for file1 in *.bpseq; do
  /Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/forForna_NuPACK.sh  ${file1}
done 



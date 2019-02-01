#!/bin/bash

cd BPSEQ 

for file1 in *.bpseq; do
  ../forForna_NuPACK.sh  ${file1}
done 



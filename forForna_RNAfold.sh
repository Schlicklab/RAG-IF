#!/bin/bash

# This script takes SEQ files first, then calculate secondary structures using RNAfold 


filename=`basename ${1} .seq`


cp ${filename}.seq ${filename}.fa

#sed -i '1s/^/>This is comment line\n/' ${filename}.fa
ex -sc '1i|>This is comment line' -cx ${filename}.fa



/opt/viennaRNA2.4.8/bin/RNAfold -p -d2 --noLP < ${filename}.fa > ${filename}.rnafold
rm *.ps
rm *.fa



# perl -nle 'print && exit if $. == 1' Top6.rnafold 

# perl -nle 'print && exit if $. == 3' Top6.rnafold |awk '{print $1}'

echo ">"${filename}-RNAfold-mfe > ${filename}.forna
perl -nle 'print && exit if $. == 2' ${filename}.rnafold >>  ${filename}.forna
perl -nle 'print && exit if $. == 3' ${filename}.rnafold |awk '{print $1}' >> ${filename}.forna

echo ">"${filename}-RNAfold-centroid >> ${filename}.forna
perl -nle 'print && exit if $. == 2' ${filename}.rnafold >>  ${filename}.forna
perl -nle 'print && exit if $. == 5' ${filename}.rnafold |awk '{print $1}' >> ${filename}.forna







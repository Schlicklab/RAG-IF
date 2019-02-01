#!/bin/bash
# This script is made to create inputs for "forna"


##
#cp Top1.bpseq.bak Top1.bpseq
##

filename=`basename ${1} .bpseq`
#echo $filename


cp ${filename}.bpseq  ${filename}.bpseq.bak
 


#sed -i '1s/^/#This is comment line\n/' $1

ex -sc '1i|#This is comment line' -cx $1
# insert a first line

#module load python/2

rm -rf ${filename}.dotbracket
rm -rf ${filename}.fa
rm -rf ${filename}.in


python /Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/bp2dotfa/bpseq2dotfa.py $1
 
#~/nupack/bin/mfe -material rna

tail -n +2 ${filename}.fa  > ${filename}.in

/opt/nupack3.2.2/bin/mfe  -material rna ${filename}
# output file is ${filename}.mfe
cp ${filename}.fa ${filename}.nupack.fa

#sed -i "1 s|$|-nupack|" ${filename}.nupack.fa 
#
sed ' 1 s/.*/&-nupack/' ${filename}.nupack.fa > tmp1
mv tmp1 ${filename}.nupack.fa




answer1=`grep -q "\.\." ${filename}.mfe ; echo $?`
answer2=`grep -q "((" ${filename}.mfe ; echo $?`
answer3=`grep -q "))" ${filename}.mfe ; echo $?`


if [ $answer1 = 0 ]; then
   grep  "\.\." ${filename}.mfe > ${filename}.nupack.dotbracket  
fi

if [ $answer2 = 0 ]; then
   grep  "((" ${filename}.mfe > ${filename}.nupack.dotbracket  
fi

if [ $answer3 = 0 ]; then
   grep  "))" ${filename}.mfe > ${filename}.nupack.dotbracket  
fi



# convert Fasta + Dot-Bracket to bpseq for nupack

python  /Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/dotfa2bp/dotfa2bpseq.py  ${filename}.nupack.fa ${filename}.nupack.dotbracket
 





# generate input for Forna

cat ${filename}.nupack.fa  ${filename}.nupack.dotbracket ${filename}.fa ${filename}.dotbracket > ${filename}.forna

# generate 


# Clean-up
rm -rf  ${filename}.in  ${filename}.fa ${filename}.dotbracket ${filename}.nupack.fa ${filename}.nupack.dotbracket ${filename}.mfe

#  
mv ${filename}.bpseq.bak ${filename}.bpseq


 

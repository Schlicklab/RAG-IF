# coding: utf-8

# modified on 08/21/2018

# This script is designed to run GAIF program to carry out the mutation generation
# and it is also designed for verify the mutated sequence candidates.


import os 
import tempfile
import sys

from checkBPSEQ import * # modified on 2018-07-16



#RNAinverseDir="/Users/yt34/NYU_Drive_Google/Work/RNA-projects/GAIF/"#"/opt/viennaRNA2.4.8/bin/"
#NuPACKDir="/opt/nupack3.2.2/bin/"
#TGpath="/Users/yt34/NYU_Drive_Google/Work/RNA-projects/modified-treeGraph/"
#fd2seqDir="/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/dotfa2bp/"
RNAinverseDir="/Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/"#"/opt/viennaRNA2.4.8/bin/"
NuPACKDir="/opt/nupack3.2.2/bin/"
TGpath="/Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/modified-treeGraph/"
fd2seqDir="/Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/"


#inpf = "Top49.mfe.rnaInverseInp"
#tmpf = "Top49.seq"
inpf = sys.argv[1] 
tmpf = sys.argv[2] # newly added for template file  

prefix = inpf.split(".")[0]


# Determine the number of sequences to be generated:
# if < 200, 4^N
# else 200



Rvalue = 500
with open(inpf,"r") as myfile:
     data = myfile.read().replace('\n','')
     freq = data.count('N')
    

#print freq 

if 4**freq < 500:
   Rvalue = 4**freq

# Run Commands !!!!!!!!!!


#os.system(RNAinverseDir+"RNAinverse -Fm −−noGU −−noClosingGU -f 0.5 -R"+str(Rvalue)+"  -d2 <"+inpf+" >"+prefix+".mut_seq1" )
os.system("python "+RNAinverseDir+"gaif.py "+ inpf +" " + tmpf)
   # output file is "heaven.txt"



# Remove duplicate generated sequences 
seq1 = []
with open("heaven.txt") as pointer:
     for line in pointer:
         #
         seq1.append( line.split()[0].upper() )
         
#print len(seq1)

seq2 = list(set(seq1))

print "Number of candidate sequences for this case:",len(seq2) # list(set(seq1)))

#seq2 = list(set(seq1))

#print seq2

# Verify these generated sequences with NuPACK and RNAfold
# one by one stored in seq2


## NuPACK
### NuPACK input - [prefix].in





final_seq = []
final_label = []
# # # # # # # # # # # # # # # # # # 
# Loop over all sequnces in seq2 # 
# # # # # # # # # # # # # # # # # 

#ex_seq = seq2[0]

process_bar = ShowProcess(len(seq2), 'OK')

for i in range(len(seq2)):
  ex_seq = seq2[i]

  good_seq = 'y'

  # Get the targeted topology
  with open(prefix+".tg_log") as pointer:
       for line in pointer:
           #
           if "Graph ID" in line:
               #print line
               target_top = line.split()[-1] 




  # Verify with RNAfold

  f1 = open("tmpRNAfold.in","w")
  f1.write(ex_seq)
  f1.close()

  os.system("/opt/viennaRNA2.3.5/bin/RNAfold -p -d2 --noLP < tmpRNAfold.in > tmpRNAfold.out" )
  os.system("rm -rf tmpRNAfold.in rna.ps dot.ps")

  f1 = open("tmpRNAfold-mfe.fa","w")
  f1.write(">test\n")
  f1.write(ex_seq)
  f1.close()

  os.system("cp tmpRNAfold-mfe.fa tmpRNAfold-centroid.fa")

  # get dot-bracket and then convert to bpseq
  # 1. mfe - 2nd line
  # 2. centroid - 4th line

  counter = 0
  with open("tmpRNAfold.out") as pointer:
       for line in pointer:
           counter = counter + 1 
           if counter == 2:
              tmpRNA_mfe_ss = line.split()[0]
           if counter == 4: 
              tmpRNA_centroid_ss = line.split()[0] 
  os.system("rm -rf tmpRNAfold.out")



  f1 = open("tmpRNAfold-mfe.dotbracket","w")
  f1.write(tmpRNA_mfe_ss)
  f1.close()
  os.system("rm -rf tmpRNAfold-mfe.bpseq")
  os.system("python "+fd2seqDir+"dotfa2bpseq.py tmpRNAfold-mfe.fa tmpRNAfold-mfe.dotbracket"+">/dev/null ")
  os.system("rm -rf tmpRNAfold-mfe.fa tmpRNAfold-mfe.dotbracket")

  f1 = open("tmpRNAfold-centroid.dotbracket","w")
  f1.write(tmpRNA_centroid_ss)
  f1.close()
  os.system("rm -rf tmpRNAfold-centroid.bpseq")
  os.system("python "+fd2seqDir+"dotfa2bpseq.py tmpRNAfold-centroid.fa tmpRNAfold-centroid.dotbracket"+">/dev/null ")
  os.system("rm -rf tmpRNAfold-centroid.fa tmpRNAfold-centroid.dotbracket ")


  # mfe - treegraph
  os.system("tail -n +2 tmpRNAfold-mfe.bpseq > tmp1")
  os.system("mv tmp1 tmpRNAfold-mfe.bpseq")
  os.system("python "+TGpath+"treeGraphs.py tmpRNAfold-mfe.bpseq > tmpRNAfold-mfe.tg_log")
  with open("tmpRNAfold-mfe.tg_log") as pointer:
       for line in pointer:
           if "Graph ID" in line:
               rnafold_mfe_top = line.split()[-1]

  os.system("rm -rf tmpRNAfold-mfe.bpseq tmpRNAfold-mfe.tg_log")
  ###print rnafold_mfe_top

  # centroid - treegraph
  os.system("tail -n +2 tmpRNAfold-centroid.bpseq > tmp1")
  os.system("mv tmp1 tmpRNAfold-centroid.bpseq")
  #os.system("python "+TGpath+"treeGraphs.py tmpRNAfold-centroid.bpseq > tmpRNAfold-centroid.tg_log")
  # modified on 2018-07-18

  ifSkip = chkBPSEQ("tmpRNAfold-centroid.bpseq") # modified on 2018-07-16
  if ifSkip == 'n':
     os.system("python "+TGpath+"treeGraphs.py tmpRNAfold-centroid.bpseq > tmpRNAfold-centroid.tg_log") 
     with open("tmpRNAfold-centroid.tg_log") as pointer:
          for line in pointer:
              if "Graph ID" in line:
                  rnafold_centroid_top = line.split()[-1]
     os.system("rm -rf tmpRNAfold-centroid.bpseq tmpRNAfold-centroid.tg_log")
  ###print rnafold_centroid_top
  else:
     rnafold_centroid_top = 'invalid'


  if (rnafold_centroid_top != target_top ) and (rnafold_mfe_top != target_top):
     good_seq = 'n'

  # Verify with NuPack

  f1 = open("tmpNupack.in","w")
  f1.write(ex_seq)
  f1.close()

  os.system("/opt/nupack3.2.2/bin/mfe -material rna tmpNupack") # 2
  os.system("rm tmpNupack.in")

  f1 = open("tmpNupack.fa","w")
  f1.write(">test\n")
  f1.write(ex_seq)
  f1.close()

  # Get dot-bracket from "tmpNupack.mfe" 
  with open("tmpNupack.mfe") as pointer:
       for line in pointer:
           #
           if "((" in line:
              target = line
           if ".." in line:
              target = line
           if "))" in line:
              target = line

  f1 = open("tmpNupack.dotbracket","w")
  f1.write(target)
  f1.close()

  os.system("rm -rf tmpNupack.bpseq")
  os.system("python "+fd2seqDir+"dotfa2bpseq.py tmpNupack.fa tmpNupack.dotbracket"+">/dev/null ")
  os.system("tail -n +2 tmpNupack.bpseq > tmp1")
  os.system("mv tmp1 tmpNupack.bpseq")


  # double check the NuPACK verification result
  ifSkipB = 'n'
  ifSkipB = chkTermi("tmpNupack.bpseq") # 07-18


  os.system("python "+TGpath+"treeGraphs.py tmpNupack.bpseq > tmpNupack.tg_log"+" 2>/dev/null") #modified 2018-07-18
  # clean up
  os.system("rm -rf tmpNupack.mfe tmpNupack.fa tmpNupack.dotbracket tmpNupack.bpseq")

  # get its top. graph id
  unpack_top = "not_assigned"
  with open("tmpNupack.tg_log") as pointer:
       for line in pointer:
           if "Graph ID" in line:
               nupack_top = line.split()[-1]

  os.system("rm -rf tmpNupack.tg_log ")

  if ifSkipB == "y": # 07-18
     nupack_top = 'verified_failed'


  if target_top != nupack_top:
     good_seq = "n"


  #print good_seq

  if good_seq == 'y':
     final_seq.append( ex_seq )

  # print out status
  stat_1 = (target_top == nupack_top)
  stat_2 = (target_top == rnafold_mfe_top )
  stat_3 = (target_top == rnafold_centroid_top)
  #print stat_1, stat_2, stat_3  #,"->",(stat_1 and stat_2 and stat_3)
  if stat_2 and stat_3:
     final_label.append("B")
  if stat_2 and (not stat_3):
     final_label.append("M")
  if (not stat_2) and stat_3:
     final_label.append("C")
  
  #
  process_bar.show_process()


print "Number of successful sequences for this case:",len(final_seq)


# MFE only -> M
# CENTROID ONLY -> C
# Both -> B

# export successful sequences 
f1 = open(prefix+".survivors","w")
for i in range(len(final_seq)):
    f1.write( final_seq[i].strip() )
    f1.write(" "+final_label[i] )
    f1.write("\n")
f1.close()





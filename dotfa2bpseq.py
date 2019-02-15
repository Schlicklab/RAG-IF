# Convert Dot-Bracket + FASTA to BPSEQ 
# Y. Tao - 2018/06/12
# modified on 2018/06/17 for the pairing problem in Bot-Bracket file 

import sys
import os
from dotPairing import  inputDotBracket 


#inpf1 = "1F1T.fa"
#inpf2 = "1F1T.dotbracket"

inpargs = sys.argv
#print inpargs
if "-h" in inpargs or "--help" in inpargs:
   print "usage: python dotfa2bpseq.py {INPUT}.fa {INPUT}.dotbracket"
   sys.exit()

if len(inpargs) != 3:
   print "Please check input.\nAbort."
   print "usage: python dotfa2bpseq.py {INPUT}.fa {INPUT}.dotbracket"
   sys.exit()

inpf1 = inpargs[1]
inpf2 = inpargs[2]


inpfilename1, inpfilesuffix1 = os.path.splitext(inpf1)
inpfilename2, inpfilesuffix2 = os.path.splitext(inpf2)

if inpfilename1 != inpfilename2: 
   print "Input file names not match.\nAbort."
   sys.exit()

if (inpfilesuffix1.lower()) != ".fa":
   print "Please check input FASTA file suffix."
   print "Abort."
   sys.exit()

if (inpfilesuffix2.lower()) != ".dotbracket":
   print "Please check input dot-bracket file suffix."
   print "Abort."
   sys.exit()

#lineN = -1
first = -1
NTlabel = []
NTln = []
with open(inpf1) as pointer:
   for line in pointer:
       strpline = line.strip()
       if len(strpline) > 0:
          first = first + 1
          if first == 0:
             if strpline[0] != ">":
                print "Please check FASTA input file format.\nAbort."
                sys.exit()
             continue
          strpline = strpline.replace(" ","")
          NTln = list(strpline)
          #print NTln
          for i in range(len(NTln)):
              if NTln[i] not in ["A","U","C","G"]:
                 print "Unknown NT for RNA.\n\Abort."
                 print NTln[i]
                 sys.exit()
              NTlabel.append( NTln[i] )  

# Flexible format for FASTA supported 
#  1. space in between
#  2. multiple lines
                
#print NTlabel

NumNT = len(NTlabel) 

# dot-bracket file
first = -1
DB = []
with open(inpf2) as pointer:
  for line in pointer:
      strpline = line.strip()
      if len(strpline) > 0:
         first = first + 1  
         DB = list(strpline)
         for i in range(len(DB)):
             if DB[i] not in ["(",")","."]:
                print "Unknown symbol for dot-bracket notation.\n\Abort."
                print BD[i]
                sys.exit()
      if first == 0:
         break 

NumNT2 = len(DB)
if NumNT != NumNT2:
   print "Inconsistent length.\nAbort."
   sys.exit()

Lp = DB.count("(")
Rp = DB.count(")")
if Lp != Rp :
   print "Incorrect pair numbers.\nAbort."
   sys.exit()


# We need to modify starting here

pairs = inputDotBracket(DB)
 
#for i in range(len(DB)):
#    for j in range(len(DB)):
#        if DB[i] == "(": 
#           if DB[j] == ")":
#              if i > j: 
#                 print "Error in pairing.\nAbort."
#                 sys.exit()
#
#
#DB2 = DB[:]
#DB2.reverse()
#pairs = [] 
#for i in range(len(DB)):
#    for j in range(len(DB2)):
#        #if i < j :
#        if DB[i] == "(" and DB2[j] == ")":
#           DB[i] = "P"
#           DB2[j]= "P"
#           #print i+1, NumNT-j
#           pairs.append( [i+1,NumNT-j] )
#           break 





       
# generate BPSEQ file
if os.path.exists(inpfilename1+".bpseq"):
   print inpfilename1+".bpseq"+" already exists. Please check whether you still need it."
   print "Aborting."
   sys.exit() 


outf = open (inpfilename1+".bpseq","w")
outf.write("#Comment line\n")
for i in range(NumNT):
    outf.write(str(i+1)+" ")
    outf.write(NTlabel[i]+" ")
    thirdC = 0
    for j in range(len(pairs)):
        if pairs[j][0] == (i+1):
           thirdC = pairs[j][1]
        if pairs[j][1] == (i+1):
           thirdC = pairs[j][0]
    outf.write(str(thirdC)+"\n")
outf.close()
print "BPSEQ file written to "+inpfilename1+".bpseq"



         



  













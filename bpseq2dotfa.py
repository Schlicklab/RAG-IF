# Convert BQSEQ format into DotBracket + FASTA
# Y. Tao - 2018/06/12
# modified for correcting the pairing in Dot-Bracket file - 2018/06/17
    

import os
import sys

#inpf = "1F1T.bpseq"
inpargs = sys.argv

if "-h" in inpargs or "--help" in inpargs:
   print "usage: python bpseq2dotfa.py {INPUT}.bpseq"
   sys.exit()

if len(inpargs) != 2:
   print "Please check input.\nAbort."
   print "usage: python bpseq2dotfa.py {INPUT}.bpseq"
   sys.exit()

inpf = inpargs[1] 

inpfilename, inpfilesuffix = os.path.splitext(inpf)

#print inpfilename
#print inpfilesuffix
#print inpfilesuffix.lower()

if (inpfilesuffix.lower()) != ".bpseq":
  print "Please check input file suffix.\nAborting."
  sys.exit()  

NumNT = 0
first = -1
NTlabel = []
NumUP = 0
with open(inpf) as pointer:
   for line in pointer:
       strpline = line.strip()
       if len(strpline) > 0:
          first = first + 1
          if first == 0: 
             continue
          #print first 
          if str(first) != strpline.split()[0]:
             print "Error found in first column of input file.\nAborting." 
             print strpline
             sys.exit()
          NTlabel.append(strpline.split()[1])
          if strpline.split()[2] == "0":
             NumUP = NumUP + 1
         

#print NTlabel


# generate FASTA file
if os.path.exists(inpfilename+".fa"):
   print inpfilename+".fa"+" already exists. Please check whether you still need it."
   print "Aborting."
   sys.exit() 

faf = open( inpfilename+".fa" ,"w")
faf.write(">FASTA file for "+inpfilename+"\n")
for i in range(len(NTlabel)):
   faf.write(NTlabel[i])
faf.write("\n")
faf.close()
print "FASTA file written to "+inpfilename+".fa"


# generate dot-bracket file
#print first
if (first-NumUP)%2 != 0:
   print "Odd number of paired NTs.\nAbort."
   sys.exit()
NPair = (first-NumUP)/2 # Number of base pairs 


# We need to modify starting here 
DotBracket=[]
for i in range(first): # Initialize the Dot-Bracket notation 
   DotBracket.append(".") 


first=-1
pairs = []
#PairC = 0
with open(inpf) as pointer:
   for line in pointer:
       strpline = line.strip()
       if len(strpline) > 0:
#          #print strpline 
          first = first + 1
          if first == 0: 
             continue
          if strpline.split()[2] != "0":
             pairs.append([ int(strpline.split()[0]), int(strpline.split()[2]) ])  
#             PairC = PairC + 1
#             if PairC <= NPair:
#                DotBracket.append("(")
#             else:
#                DotBracket.append(")")
#          else:
#             DotBracket.append(".")

if (2*NPair) != len(pairs): 
   print "Please double check number of NT pairs.\nAbort."
   sys.exit()

while NPair >= 1: 
   #for i in range(len(pairs)):
   a = pairs[0][0]
   b = pairs[0][1]  
   DotBracket[a-1]="("
   DotBracket[b-1]=")"
   pairs.remove([a,b])
   pairs.remove([b,a])
   NPair = NPair - 1
     

 


#print DotBracket
if os.path.exists(inpfilename+".dotbracket"):
   print inpfilename+".dotbracket"+" already exists. Please check whether you still need it."
   print "Aborting."
   sys.exit() 

dbf = open( inpfilename+".dotbracket" ,"w")
for i in range(len(DotBracket)):
   dbf.write(DotBracket[i])
dbf.write("\n")
dbf.close()
print "Dot-Bracket file written to "+inpfilename+".dotbracket"







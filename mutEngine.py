# This mutation engine script is designed to generate the input file for 
# RNAinverse program from ViennaRNA package 
# Y. Tao - 2018/07/02


import os 
import sys 


# We need to make this script command-line compatiable 
inpargs = sys.argv

if "-h" in inpargs or "--help" in inpargs:
   print "usage: python mutEngine.py  [desired_ss].bpseq  *-ifmutate.log "
   sys.exit()

if len(inpargs) != 3:
   print "Please check input.\nAbort."
   print "usage: python mutEngine.py  [desired_ss].bpseq  *-ifmutate.log "
   sys.exit()



# Set up the location where the conversion script is 
#bpseq2dotfaDir = "/Users/yt34/NYU_Drive_Google/Work/RNA-projects/myOwnScripts/bp2dotfa/"
bpseq2dotfaDir = "/Users/sj78/Documents/labwork/MutationsForDesign/RAG-IF_Code/"


# BPSEQ file containing desired secondary structure
#inpf1 = "Top151.bpseq"
inpf1 = inpargs[1]
#### BPSEQ file containing current secondary structure 
#### inpf2 =
#### This file might not be necessary after testing with RNAinverse 


# Output file from s2_script.py 
#inpf2 = "Top151-ifmutate.log"
inpf2 = inpargs[2]



# Task I - convert inpf1 into dot-bracket file 
if bpseq2dotfaDir[-1] != '/':
   bpseq2dotfaDir = bpseq2dotfaDir + '/' 

#os.system ( "sed -i '1s/^/#This is comment line\n/'" + inpf1 ) # Add the first line 
os.system( "ex -sc '1i|#This is comment line' -cx " + inpf1 ) # adjusted for Mac OS 

os.system("rm -rf *.fa")
os.system("rm -rf *.dotbracket")
os.system("python "+ bpseq2dotfaDir + "bpseq2dotfa.py "+ inpf1 )

os.system ("tail -n +2 "+ inpf1 + "> tmpf1")  # delete the first line
os.system ("mv tmpf1 " + inpf1  )             # delete the first line 
os.system("rm -rf *.fa")

inpf1filename,inpf1suffix = os.path.splitext(inpf1)
# the generate dot-bracket file is {inpf1filename}.dotbracket


# Task II - generate the sequence with lower-cases as constant residues 

flag = -1
seqStr = ''
with open(inpf2) as pointer:
     for line in pointer:
         if "Print IFMUTATE info" in line:
            flag = 0
            continue
         if flag >= 0:
            flag = flag + 1 
            if len(line.strip()) != 0:
               if len( line.split() ) == 5 :
                  #
                  if line.split()[-1] == "no": # Easy to make mistakes here
                     seqStr = seqStr + line.split()[1].lower()
                  else:
                     seqStr = seqStr + "N"   #line.split()[1] 

         
#print seqStr

# Task III - generate the input file for RNAinverse

os.system( "mv "+inpf1filename+".dotbracket "+inpf1filename+".rnaInverseInp"  )
with open(inpf1filename+".rnaInverseInp","a") as myfile:
    myfile.write(seqStr+"\n")

num_lines = sum(1 for line in open( inpf1filename+".rnaInverseInp" ))
if num_lines == 2:
   print "RNAinverse input file written to",inpf1filename+".rnaInverseInp"








# This script is designed to analyze the difference between two RNA 2D structures of the same RNA molecule
# as the basis for automatic mutation.
# Y. Tao - 2018/06/18
# S. Jain - 2019/03/22

import sys
from analyzer_util import *

def read_input(file,string):

    topo=""
    info=[]
    start = -1
    with open(file) as pointer:
        for line in pointer:
            strpline = line.strip()
            if len(strpline) > 0:
                if "Graph ID" in strpline:
                    topo = strpline.split()[-1]
                    print "RAG ID for " + string + " RNA is",topo.strip()
                    start = 0
                    continue
                if start > -1 :
                    if int(strpline.split()[0]) != start+1:
                        print "Warning: Inconsistency found in first column.\nAbort."
                        #print strpline
                        #sys.exit()
                        pass
                    nt = residue(int(strpline.split()[1]),strpline.split()[2],int(strpline.split()[3]),strpline.split()[4],int(strpline.split()[5]))
                    info.append(nt)
                    start = start + 1
    return info
#-----------------

#------------------#
#------MAIN--------#
#------------------#

if "-h" in sys.argv or "--help" in sys.argv:
   print "usage: python analyzer.py {current_topology}.log {target_topology}.log"
   sys.exit()

if len(sys.argv) != 3:
   print "Please check input.\nAbort."
   print "usage: python analyzer.py {current_topology}.log {target_topology}.log"
   sys.exit()

# I. Read input file
#
#   Input: # The input files are generate by TreeGraph code (modified)
#

# I.1 Read in file A
print "*Basic Info."
info_A = read_input(sys.argv[1],"first") # reading the input files and extracting residue informations
start_A = len(info_A)+1

# I.2 Read in file B
info_B = read_input(sys.argv[2],"second")
start_B = len(info_B)+1

#print start_A, start_B
if start_A == start_B:
   print "Length of RNA:",start_A
else:
   print "Length of RNA does not match in two files.\nAbort."
   sys.exit()

# II. Find common vertices
#
#   Input: info_A, info_B
#

# II.1 Find all vertices and related residues 

#NumVerA = 0
#NumVerB = 0

VerA = extractVer( info_A )
VerB = extractVer( info_B )
#print VerA
#print VerB
NumVerA = len(VerA) - 1 # as VerA and VerB contain residues with loopID = 0, which is the loop label for all helic residues
NumVerB = len(VerB) - 1

print "Number of vertices in first RNA:", NumVerA
print "Number of vertices in second RNA:",NumVerB

# II.2 Correlate vertices between two RNAs

a2b = verCorr1(VerA,VerB) # identifies how many loops of 2 are residues in loops of 1 are divided into
b2a = verCorr1(VerB,VerA)

#print "a2b, b2a 1"
#print a2b
#print b2a

# II.2.a Correlate identical vertices 
one2one = easyCorr(a2b,b2a) # This function may change a2b and b2a.

if len(one2one) != 0:
    print "*Comparative Info."
    print "Number of identical vertices:",len(one2one)
    print "Correlation result 1:"
    print one2one
#print "a2b, b2a 2"
#print a2b
#print b2a

# II.2.b Correlate enlarged/shrunk vertices
one2oneI,one2oneII = minorVariation(a2b,b2a)

if len(one2oneI) != 0 or len(one2oneII) != 0:
   print "Number of enlarged/shrunk vertices:",len(one2oneI)+len(one2oneII)
   if len(one2oneI) != 0:
      pass
      print "Correlation result 2a:"
      print one2oneI
   if len(one2oneII) != 0:
      pass
      print "Correlation result 2b:"
      print one2oneII

#print one2oneI
#print one2oneII
#print "a2b, b2a 3"
#print a2b
#print b2a

# II.2.c Find Two-to-One [first-2,second-1] vertex sets
two2one = easySplit(a2b,b2a)

if len(two2one) != 0:
   print "Number of Two-to-One [first-2,second-1] vertex sets:",len(two2one)
   for i in range(len(two2one)):
       m = two2one[i][0][0]
       n = two2one[i][0][1]
       s = two2one[i][1]
       print "Correlation result 3:"+"("+m+","+n+") -> "+s      

#print "a2b, b2a 4"
#print a2b
#print b2a
#print two2one

# II.2.d  Unsuccessful pairing for 5'- and 3'-termini
      
termis = checkTermini(start_A,a2b,b2a,VerA,VerB)
if len(termis) != 0:
   pass
   print "Unsuccessful pairing at termini detected."
   print "Correlation result 4:"
   print termis

#print "a2b, b2a, 5"
#print a2b
#print b2a

# II.2.e Find inter-changable vertex dimer (pair)

two2two = findInterchangablePair(a2b,b2a)
if len(two2two) != 0:
   print "Number of Two-to-Two [first-2,second-2] vertex pairs:", len(two2two)/2
   print "Correlation result 5:"
   print two2two
#print "a2b, b2a, 6"
#print a2b
#print b2a

if len(a2b) == 0 and len(b2a) == 0:
   print "All vertices have been analyzed."
else:
   print "The siuation is more complicated."

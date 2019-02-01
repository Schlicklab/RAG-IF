# This script is designed to analyze the difference between two RNA molecules 
# as the basis for automatic mutation.
#
# Y. Tao - 2018/06/18
#

import sys
from analyzer_util import *

# The input files are generate by TreeGraph code (modified)

inpargs = sys.argv

if "-h" in inpargs or "--help" in inpargs:
   print "usage: python analyzer.py {current_topology}.log {target_topology}.log"
   sys.exit()

if len(inpargs) != 3:
   print "Please check input.\nAbort."
   print "usage: python analyzer.py {current_topology}.log {target_topology}.log"
   sys.exit()




#inpf1 = "home-ii-b-current.log"
#inpf2 = "home-ii-b-target.log"

inpf1 = inpargs[1]
inpf2 = inpargs[2]

# I.1 Read in file A 
topo_A=""
info_A=[]
start_A = -1
with open(inpf1) as pointer: 
   for line in pointer: 
       strpline = line.strip()
       if len(strpline) > 0:
           if "Graph ID" in strpline:
               topo_A = strpline.split()[-1]
               print "*Basic Info."
               print "RAG ID for first RNA is",topo_A.strip()
               start_A = 0
               continue
           if start_A > -1 :
               if int(strpline.split()[0]) != start_A+1:
                    
                   print "Warning: Inconsistency found in first column.\nAbort."
                   #print strpline
                   #sys.exit()
                   pass
               nt = residue( int(strpline.split()[1]), \
                             strpline.split()[2], \
                             int(strpline.split()[3]), \
                             strpline.split()[4], \
                             int(strpline.split()[5]) )
               info_A.append(nt) 
               start_A = start_A + 1 

#print start_A
#


# I.2 Read in file B 
topo_B=""
info_B=[]
start_B = -1
with open(inpf2) as pointer: 
   for line in pointer: 
       strpline = line.strip()
       if len(strpline) > 0:
           if "Graph ID" in strpline:
               topo_B = strpline.split()[-1]
               print "RAG ID for second RNA is",topo_B.strip()
               start_B = 0
               continue
           if start_B > -1 :
               if int(strpline.split()[0]) != start_B+1:
                   print "Warning: Inconsistency found in first column.\nAbort."
                   #sys.exit()
                   #firstC = start_B + 1
                   pass
               else:
                   #firstC = int(strpline.split()[1])
                   pass

               nt = residue( int(strpline.split()[1]), \
                             strpline.split()[2], \
                             int(strpline.split()[3]), \
                             strpline.split()[4], \
                             int(strpline.split()[5]) )
               info_B.append(nt) 
               start_B = start_B + 1 

#print start_A, start_B
if start_A == start_B:
   print "Length of RNA:",start_A
else:
   print "Length of RNA not match in two files.\nAbort."
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
NumVerA = len(VerA) - 1
NumVerB = len(VerB) - 1

print "Number of vertices in first RNA:", NumVerA
print "Number of vertices in second RNA:",NumVerB

# II.2 Correlate vertices between two RNAs  

#print "VerA, VerB"
#print VerA
#print VerB

a2b = verCorr1(VerA,VerB)
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








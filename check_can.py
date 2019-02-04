# This script is designed to find potential candidate for successful design.# Y. Tao - 2018-07-05
# Changes by Swati Jain - 02/01/2019 to optimize the code and remove hard coded arguments
import sys

directory=sys.argv[1]
desiredTop = sys.argv[2]

inpf1 = directory + "/Nupack_Rank_Topo.txt"
inpf2 = directory + "/RNAfold_Rank_Topo.txt"

data1 = []
# read in file 1 - Nupack topologies
with open (inpf1) as pointer:
     for line in pointer:
         if len(line) > 1 and (len(line.split()) == 2)  :
            #print line

            dat = [ line.split()[0], line.split()[1] ]
            data1.append( dat )
#print data1

data2 = []
# read in file 2 - RNAfold topologies
with open (inpf2) as pointer:
     for line in pointer:
         if len(line) > 1 and (len(line.split())==3):
            dat = [ line.split()[0], line.split()[1], line.split()[2] ]
            data2.append( dat )
#print data2

# check I - successful candidates [loose]
# S.J. - fold onto same topology by both RNAfold (either mfe or centroid) and NUPACK
# check Ib - successful candidates [loose]
# S.J. - fold onto the target topology (desiredTop variable) by both RNAfold (either mfe or centroid) and NUPACK
# check II - successful candidates [tight]
# S.J. - fold onto the same topology by both RNAfold (both mfe and centroid) and NUPACK
# check IIb - successful candidates [tight]
# S.J. - fold onto the target topology (desiredTop variable) by both RNAfold (both mfe and centroid) and NUPACK
# Check IIIa - only one side is desired topology
# S.J. - NUPACK folds onto the target topology but RNAfold (both mfe and centroid) do not
# Check IIIb - only one side is desired topology
# S.J. - RNAfold (both mfe and centroid) folds onto the target topology but NUPACK does not
# S.J. - RNAfold (only mfe) folds onto the target topology but NUPACK does not
# S.J. - RNAfold (only centroid) folds onto the target topology but NUPACK does not

# S.J. combine all for loops into 1, don't need to read it again and again - 02/04/2019
for i in range(len(data1)):
    flag=0
    for j in range(len(data2)):
        if data1[i][0] == data2[j][0]: # same rank
            flag=1 # sequence found
                
            if (data1[i][1] == data2[j][1]) or (data1[i][1] == data2[j][2]): # NUPACK equal to either one of the RNAfold
                #print "Candidate Type I:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
                
            if (data1[i][1] == data2[j][1] == desiredTop) or (data1[i][1] == data2[j][2] == desiredTop):# NUPACK equal to either one of the RNAfold, but also with desired topology
                #print "Candidate Type Ib:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
                    
            if (data1[i][1] == data2[j][1]) and (data1[i][1] == data2[j][2]): # NUPACK equal to both of the RNAfold
                #print "Candidate Type II:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
                    
            if (data1[i][1] == data2[j][1] == desiredTop) and (desiredTop == data1[i][1] == data2[j][2]): # NUPACK equal to both of the RNAfold, but also with desired topology
                #print "Candidate Type IIb:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
                    
            if (data1[i][1] == desiredTop ) and (data2[j][1] != desiredTop) and (data2[j][2] != desiredTop): # NUPACK folds onto the target topology but RNAfold (both mfe and centroid) do not
                #print "Candidate Type IIIa:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
                    
            if (data1[i][1] != desiredTop ) and (data2[j][1] == desiredTop) and (data2[j][2] == desiredTop): # RNAfold (both mfe and centroid) folds onto the target topology but NUPACK does not
                print "Candidate Type IIIb:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
                    
            if (data1[i][1] != desiredTop ) and (data2[j][1] == desiredTop) and (data2[j][2] != desiredTop): # RNAfold (only mfe) folds onto the target topology but NUPACK does not
                print "Candidate Type IIIb:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
                    
            if (data1[i][1] != desiredTop ) and (data2[j][1] != desiredTop) and (data2[j][2] == desiredTop): # RNAfold (only centroid) folds onto the target topology but NUPACK does not
                print "Candidate Type IIIb:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
    
        if flag == 1: # sequence found and type assigned
            break

# check I - successful candidates [loose]
#for i in range(len(data1)):
#    for j in range(len(data2)):
#        if data1[i][0] == data2[j][0]:
#            if (data1[i][1] == data2[j][1]) or (data1[i][1] == data2[j][2]): # NUPACK equal to either one of the RNAfold
#              print "Candidate Type I:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
#print ""

# check Ib - successful candidates [loose]
#for i in range(len(data1)):
#    for j in range(len(data2)):
#        if data1[i][0] == data2[j][0]:
#           if (data1[i][1] == data2[j][1] == desiredTop) or (data1[i][1] == data2[j][2] == desiredTop):
#              print "Candidate Type Ib:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
#print ""

# check II - successful candidates [tight]
#for i in range(len(data1)):
#    for j in range(len(data2)):
#        if data1[i][0] == data2[j][0]:
#           if (data1[i][1] == data2[j][1]) and (data1[i][1] == data2[j][2]): # NUPACK equal to both of the RNAfold
#              print "Candidate Type II:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
#print ""

# check IIb - successful candidates [tight]
#for i in range(len(data1)):
#    for j in range(len(data2)):
#        if data1[i][0] == data2[j][0]:
#           if (data1[i][1] == data2[j][1] == desiredTop) and (desiredTop == data1[i][1] == data2[j][2]):
#              print "Candidate Type IIb:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
#print ""

# Check IIIa - only one side is desired topology
#for i in range(len(data1)):
#    for j in range(len(data2)):
#        if data1[i][0] == data2[j][0]:
#           if (data1[i][1] == desiredTop ) and (data2[j][1] != desiredTop) and (data2[j][2] != desiredTop) :
#              print "Candidate Type IIIa:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
#print ""

# Check IIIb - only one side is desired topology

#for i in range(len(data1)):
#    for j in range(len(data2)):
#        if data1[i][0] == data2[j][0]:
#           if (data1[i][1] != desiredTop ) and (data2[j][1] == desiredTop) and (data2[j][2] == desiredTop) :
#              print "Candidate Type IIIb:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]


#for i in range(len(data1)):
#    for j in range(len(data2)):
#        if data1[i][0] == data2[j][0]:
#           if (data1[i][1] != desiredTop ) and (data2[j][1] == desiredTop) and (data2[j][2] != desiredTop) :
#              print "Candidate Type IIIb:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]


#for i in range(len(data1)):
#    for j in range(len(data2)):
#        if data1[i][0] == data2[j][0]:
#           if (data1[i][1] != desiredTop ) and (data2[j][1] != desiredTop) and (data2[j][2] == desiredTop) :
#              print "Candidate Type IIIb:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
#print ""













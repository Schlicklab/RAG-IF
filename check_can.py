# This script is designed to find potential candidate for successful design.# Y. Tao - 2018-07-05



inpf1 = "Nupack_Rank_Topo.txt"
inpf2 = "RNAfold_Rank_Topo.txt"
desiredTop = "7_4" # change this according to the target topology


data1 = []
# read in file 1 
with open (inpf1) as pointer:
     for line in pointer:
         if len(line) > 1 and (len(line.split()) == 2)  :
            #print line

            dat = [ line.split()[0], line.split()[1] ]
            data1.append( dat )

#print data1 


data2 = []
# read in file 2 
with open (inpf2) as pointer:
     for line in pointer:
         if len(line) > 1 and (len(line.split())==3):
            dat = [ line.split()[0], line.split()[1], line.split()[2] ]
            data2.append( dat )

#print data2 


# check I - successful candidates [loose]   
for i in range(len(data1)):
    for j in range(len(data2)):
        if data1[i][0] == data2[j][0]:
           if (data1[i][1] == data2[j][1]) or (data1[i][1] == data2[j][2]):
              #
              print "Candidate Type I:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
print ""



# check Ib - successful candidates [loose]   
for i in range(len(data1)):
    for j in range(len(data2)):
        if data1[i][0] == data2[j][0]:
           if (data1[i][1] == data2[j][1] == desiredTop) or (data1[i][1] == data2[j][2] == desiredTop):
              #
              print "Candidate Type Ib:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
print ""





# check II - successful candidates [tight]   
for i in range(len(data1)):
    for j in range(len(data2)):
        if data1[i][0] == data2[j][0]:
           if (data1[i][1] == data2[j][1]) and (data1[i][1] == data2[j][2]):
              #
              print "Candidate Type II:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
print ""


# check IIb - successful candidates [tight]   
for i in range(len(data1)):
    for j in range(len(data2)):
        if data1[i][0] == data2[j][0]:
           if (data1[i][1] == data2[j][1] == desiredTop) and (desiredTop == data1[i][1] == data2[j][2]):
              #
              print "Candidate Type IIb:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
print ""


# Check IIIa - only one side is desired topology
for i in range(len(data1)):
    for j in range(len(data2)):
        if data1[i][0] == data2[j][0]:
           if (data1[i][1] == desiredTop ) and (data2[j][1] != desiredTop) and (data2[j][2] != desiredTop) :  
              #
              print "Candidate Type IIIa:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
print ""


# Check IIIb - only one side is desired topology
for i in range(len(data1)):
    for j in range(len(data2)):
        if data1[i][0] == data2[j][0]:
           if (data1[i][1] != desiredTop ) and (data2[j][1] == desiredTop) and (data2[j][2] == desiredTop) :  
              #
              print "Candidate Type IIIb:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]


for i in range(len(data1)):
    for j in range(len(data2)):
        if data1[i][0] == data2[j][0]:
           if (data1[i][1] != desiredTop ) and (data2[j][1] == desiredTop) and (data2[j][2] != desiredTop) :
              print "Candidate Type IIIb:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
           
           
for i in range(len(data1)):
    for j in range(len(data2)):
        if data1[i][0] == data2[j][0]:
           if (data1[i][1] != desiredTop ) and (data2[j][1] != desiredTop) and (data2[j][2] == desiredTop) :
              print "Candidate Type IIIb:",data1[i][0],data1[i][1],data2[j][1],data2[j][2]
print ""













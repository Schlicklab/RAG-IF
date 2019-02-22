# This script is designed to find stems connecting two
# two known vertices. These known vertices have been correlated
# between two topologies using the previous 'analyzer' script.
#
#
#
# Y. Tao - 2018-06-28
#        - modified: 2018-08-20 [current script takes three input files]
#        - modified: 2018-09-11 (two loops in 2->1 case should also be labeled as 'yes')
#       -> modified: 2018-09-18 [if there should be a loop in a stem connecting two known vertices]  
#


from s2_script_util import *
import sys


if "-h" in sys.argv or "--help" in sys.argv:
   print "usage: python s2_script.py {analyzer-output-file} {TreeGraph(modified)-output-file(left)} {TreeGraph(modified)-output-file(right)} "
   sys.exit()



if len(sys.argv) != 4:
   print "Please check input."
   print "usage: python s2_script.py {analyzer-output-file} {TreeGraph(modified)-output-file(left)} {TreeGraph(modified)-output-file(right)} "
   sys.exit()




inpf_analyzerOut = sys.argv[1] #"*.summary1" # analyzer output file as one input 
inpf_leftInfo = sys.argv[2]    #"*-nupack-tree.log" 
inpf_rightInfo = sys.argv[3]   # *-RNAfold-tree.log


print "S1 script output as input:",inpf_analyzerOut
print "Left topology information:", inpf_leftInfo
print "Right topology information:", inpf_rightInfo

# here 'cr' stands for "Correlation Result"
cr1 = 0
cr1str = ''
cr2a = 0
cr2astr = ''
cr2b = 0
cr2bstr = ''

cr3 = 0
cr3N = 0
cr3count = 999
cr3str = []


cr4 = 0
cr4str = ''

cr5 = 0
cr5str = ''


# Read in information
with open(inpf_analyzerOut) as pointer:
    for line in pointer:

        if "Correlation result 1" in line:
           cr1 = 1
           continue
        if cr1 == 1:
           cr1str = line.strip() 
           cr1 = 2

        if "Correlation result 2a" in line:
           cr2a = 1
           continue
        if cr2a == 1:
           cr2astr = line.strip()
           cr2a = 2

        if "Correlation result 2b" in line:
           cr2b = 1
           continue
        if cr2b == 1:
           cr2bstr = line.strip()
           cr2b = 2

        if "Number of Two-to-One" in line:
           cr3 = 1
           cr3N = int(line.split()[-1])
           #print "cr3N",cr3N
           cr3count = 0
           continue
        if cr3count < cr3N : # read in "Correlation result 3" line(s)
           txt = line.split(':')[1] 
           cr3str.append( txt.strip() )
           cr3count = cr3count + 1

        if "Correlation result 4" in line:
           cr4 = 1
           continue
        if cr4 == 1:
           cr4str = line.strip()
           cr4 = 2
         
        if "Correlation result 5" in line:
           cr5 = 1
           continue
        if cr5 ==1:
           cr5str = line.strip()
           cr5 = 2 





# Extract related vertices 
LeftV = []
RightV = []



##  Correlation result 1
cr1dic =  procDicStr(cr1str)
LeftV,RightV = extDicList( cr1dic,LeftV,RightV )


## Correlation result 2a
cr2adic = procDicStr( cr2astr )
LeftV,RightV = extDicList( cr2adic, LeftV, RightV)


## Correlation result 2a
cr2bdic = procDicStr( cr2bstr )
LeftV,RightV = extDicList( cr2bdic, LeftV, RightV)

## Correlation result 3 

toBreak = []
#print cr3str
LeftV, RightV,toBreak = procCr3StrList( cr3str, LeftV, RightV, toBreak )




## Correlation result 4 
cr4dic = procDicStr(cr4str) 
LeftV,RightV = extDicList( cr4dic,LeftV,RightV )


## Correlation 5 
cr5dic = procDicStr(cr5str)
LeftV,RightV = extDicList( cr5dic,LeftV,RightV )

print "Left known vertices:",LeftV
print "Right known vertices:",RightV
if len(toBreak) != 0:
   print "Two-to-One Left side:",toBreak # this should be labeled as 'yes' - 2018-09-11 


# Find all correlated vertices 
# Then find all edges and those edges that need to be mutated 

## Read in the left topology information

info_left = []
start_A = -1 
with open(inpf_leftInfo) as pointer:
   for line in pointer:
       strpline = line.strip()
       if len(strpline) > 0:
          if "Graph ID" in strpline:
              start_A = 0
              continue
          if start_A > -1:
             start_A = start_A + 1 # counter 
             nt = residue( int(start_A), \
                           int(strpline.split()[1]), \
                           strpline.split()[2], \
                           strpline.split()[3], \
                           strpline.split()[4], \
                           strpline.split()[5]      )  
             #print nt.ifmutate
             info_left.append( nt )


## Read in the right topology information
info_right = []
start_A = -1 
with open(inpf_rightInfo) as pointer:
   for line in pointer:
       strpline = line.strip()
       if len(strpline) > 0:
          if "Graph ID" in strpline:
              start_A = 0
              continue
          if start_A > -1:
             start_A = start_A + 1 # counter 
             nt = residue( int(start_A), \
                           int(strpline.split()[1]), \
                           strpline.split()[2], \
                           strpline.split()[3], \
                           strpline.split()[4], \
                           strpline.split()[5]      )  
             #print nt.ifmutate
             info_right.append( nt )




## Label those residues on the loops which are not correlated  as 'yes'

#for i in range(len(info_left)):
#    info_left[i].prtMut()

#printMutInf(info_left)


for i in range(len(info_left)):
    if info_left[i].loopID != '0':
       if info_left[i].loopID not in LeftV:
          info_left[i].ifmutate = 'yes'
          #print i+1 

#printMutInf(info_left)

## scan all stems

stemPool = []
for i in range(len(info_left)):
    if info_left[i].stemID != '0':
       stemPool.append( info_left[i].stemID )

stemPool = list(set(stemPool))

print "Stem list:",stemPool


for i in range(len(stemPool)): # Loop over all stems one by one 
    aStem = stemPool[i] # Current stem  
    thisStem = []
    thatStem = []
    for j in range(len(info_left)): # Loop over all NTs 
        if info_left[j].stemID == aStem:
           thisStem.append([j, info_left[j].loopID] ) # collect these NTs into a list            
           pass


    #print thisStem
    toDel = []
    for j in range(len(thisStem)):
        if thisStem[j][1] == '0': # leave out the NTs in the middle of stems 
           toDel.append( thisStem[j] )
    for j in toDel:
        thisStem.remove( j )


    #print thisStem
    if len(thisStem) != 4:
       print "Unexpected stem!\nAbort."
       sys.exit()

    secondPart = []
    for j in range(len(thisStem)):
        secondPart.append( thisStem[j][1] ) # secondPart collects the loop IDs

    #print secondPart
    secondPart = list(set(secondPart))  # 4 -> 2
    #print secondPart

    if len(secondPart) != 2:
       print "Unexpected stem!\nAbort."
       sys.exit()

    # new issue - 09/18/2018
    # we need to consider this situation: there should be a loop in the stem which connects two known vertices
    # in this case, this loop was not considered and should be added
    if (secondPart[0] in LeftV) and (secondPart[1] in LeftV):
       for j in range(len(info_left)):
           if info_left[j].stemID != '0':
              if info_left[j].stemID == aStem:
                  if info_left[j].loopID == '0':
                     #if info_left[j].loopID not in LeftV: 
                        if info_right[j].loopID != '0':
                           if info_right[j].loopID not in RightV: # a loop which is not correlated in target 
                              info_left[j].ifmutate = 'yes'        

                 #  

       pass 


    # the following part should be corrected ... - 08/20/2018 
    # secondPart collects the loop IDs
    if (secondPart[0] not in LeftV) or (secondPart[1] not in LeftV): # at least one end is unknown loop
                                    # changed 'or' into 'and' - 08/20/2018 -> both ends are unknown
                                    # changed back to 'or'
        for j in range(len(info_left)):
            if info_left[j].stemID != '0':
                if info_left[j].stemID == aStem:
                   info_left[j].ifmutate = 'yes'
        pass

    # sometimes, if one end loop is correlated while the other end loop is not correlated,
    # the stem in between may stay the same. In this situation, the residues on such a stem
    # should stay the same and not to be mutated.
    if ((secondPart[0]     in LeftV) and (secondPart[1] not in LeftV)) or \
       ((secondPart[0] not in LeftV) and (secondPart[1]     in LeftV)):
       # compare this stem in left and right RNAs 
       # input: thatStem (not thisStem), info_left, info_right
       # output: info_left[x].ifmutate 
       for j in range(len(info_left)): # Loop over all NTs 
           if info_left[j].stemID == aStem:
              thatStem.append([j, info_left[j].loopID] ) # collect these NTs into a list            
              # thatStem is identical to the original thisStem 

       for p in range(len(thatStem)):
           j = thatStem[p][0]
           if info_left[j].pairNum != 0: 
              if info_left[j].pairNum == info_right[j].pairNum:
                 info_left[j].ifmutate = 'no' 
                 pass

       # 




       pass         
    #





    # Check into Two-to-One cases
    if len(toBreak) != 0:
       for j in range(len(toBreak)):
           a = toBreak[j][0]
           b = toBreak[j][1]
           if (secondPart[0] == a and secondPart[1] == b) or \
              (secondPart[0] == b and secondPart[1] == a):
              #
              for k in range(len(info_left)):
                  if info_left[k].stemID != '0':
                     if info_left[k].stemID == aStem:
                        info_left[k].ifmutate = 'yes'


    # END OF LOOPING OVER STEMS 



# over-write the two loops of '2->1' as "yes"
# added on 2018-09-11
for i in range(len(toBreak)):
    a = toBreak[i][0]
    b = toBreak[i][1]
    for k in range(len(info_left)):
        if info_left[k].loopID == a or info_left[k].loopID == b:
           if info_left[k].pairNum == 0:
              info_left[k].ifmutate = 'yes'





#
print "Print IFMUTATE info.:"
printMutInf(info_left)



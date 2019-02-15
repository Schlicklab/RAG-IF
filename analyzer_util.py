class residue:
    def __init__(self,pairNum,NTlabel,stemID,ss,loopID):
       self.pairNum = pairNum
       self.NTlabel = NTlabel
       self.stemID = stemID
       self.ss = ss
       self.loopID = loopID



def extractVer( objList ): 
    # Take the list of residue objects, return dictionary
    num = len(objList)
    dic0 = {}# "1":[2,3,7]
    for i in range(num):
        if str(objList[i].loopID) not in dic0: 
           dic0[str( objList[i].loopID )] = [i+1] #[objList[i].loopID]
        else:
           dic0[str( objList[i].loopID )].append(i+1) #( objList[i].loopID ) 
     
    return dic0


def findwhichVer(ntLabel,dic0):
     result = ""
     for key in dic0:
        if ntLabel in dic0[key]:
            result = key

     return result



def verCorr1(dicA,dicB):

     result = {} # correlation dictionary

     #keysA = []
     for key in dicA:
        if key != "0": #specify one vertex     
           ntLabels = dicA[key][:]
           #print ntLabels
           result[key] = []
           for i in range(len(ntLabels)): # loop over all NT in one V.
               inKeyB = findwhichVer(ntLabels[i],dicB)
               if inKeyB not in result[key]:
                   result[key].append (inKeyB)
               else:
                   continue

     #print result
     return result




def verCorr2(dicA,dicB): # Reverse correlation 
     return verCorr1(dicB,dicA)


def easyCorr(a2b,b2a):
    # Find one-to-one correlation pair
    pass
    result = {}
    toDelA = []
    toDelB = []
    for key in a2b:
        if len(a2b[key]) == 1:
            b0 = a2b[key][0] 
            if b0 in b2a: 
              if len(b2a[b0]) == 1:
                if b2a[b0][0] == key:
                   #print key, b0
                   result[key] = b0
                   #del a2b[key]
                   #del b2a[b0]
                   toDelA.append(key)
                   toDelB.append(b0)

    for key in toDelA:
        del a2b[key]
    for b0 in toDelB:
        del b2a[b0]
    return result 

def minorVariation(a2b,b2a):
    # Find the situations of [a] -> [0,a] and [0,b] -> [b]
    #                    and [0,c] -> [c,0] 

    pass
    result1 = {}
    result2 = {}
    toDelA = []
    toDelB = []

    # [a] -> [0,a] 
    for key in a2b:
        if len(a2b[key]) == 1:
            b0 = a2b[key][0]
            if b0 in b2a: 
              if len(b2a[b0]) == 2:
                if (key in b2a[b0]) and ("0" in b2a[b0]):
                   result1[key] = b0
                   toDelA.append(key)
                   toDelB.append(b0)
    #print result1


    # [0,b] -> [b]
    for b0 in b2a:
        if len(b2a[b0]) == 1:
            key = b2a[b0][0]
            if key in a2b: 
              if len(a2b[key]) == 2:
                if (b0 in a2b[key]) and ("0" in a2b[key]):
                   result2[key] = b0 # This part has not been tested 
                   toDelA.append(key)
                   toDelB.append(b0)
    #print result2


    #print "minor"
    # [0,c] -> [c,0]
    for key in a2b:
        if len(a2b[key]) == 2:
           go1 = 0 
           if a2b[key][0] == '0':
                 b0 = a2b[key][1]
                 go1 = 1
           if a2b[key][1] == '0':
                 b0 = a2b[key][0]
                 go1 = 1
           if go1 == 1:
              #print key,b0
              if b0 in b2a:
                 if len(b2a[b0]) == 2:
                    if '0' in b2a[b0]:
                       if key in b2a[b0]:
                          result2[key] = b0 # no bother result3 
                          toDelA.append(key)
                          toDelB.append(b0)

    for key in toDelA:
        del a2b[key]
    for b0 in toDelB:
        del b2a[b0]
    
    return result1, result2  
 
def easySplit(a2b,b2a):
    #print a2b
    #print b2a
    # Find the situation of [a] -> [c], [b] -> [c]
    # Improvement: the initial implementation only considers the situation of "c":[a,b]
    #              the improvement also considers the situation of "c":[a,b,0] -> enlarged loop in target RNA  
    # Search in one-way, instead of two-way
    pass 
    result = [] # Note here is a list instead of dic. 
    toDelA = []
    toDelB = []
    for b0 in b2a:
        if len(b2a[b0]) == 3: # Added on 2018-06-21 
           pass
           if '0' in b2a[b0]:
              three2two = b2a[b0][:]
              three2two.remove('0')
              key1 = three2two[0]
              key2 = three2two[1]
              if len(a2b[key1]) == 1:
                if a2b[key1][0] == b0:
                  if len(a2b[key2]) == 1:
                    if a2b[key2][0] == b0:
                       result.append([[key1,key2],b0]) 
                       toDelA.append(key1)
                       toDelA.append(key2)
                       toDelB.append(b0)
 
       
        if len(b2a[b0]) == 2:
           key1 = b2a[b0][0]
           key2 = b2a[b0][1]
           if (key1 in a2b) and (key2 in a2b):
             if len(a2b[key1]) == 1:
               if a2b[key1][0] == b0:
                 if len(a2b[key2]) == 1:
                   if a2b[key2][0] == b0:
                      result.append([[key1,key2],b0]) 
                      toDelA.append(key1)
                      toDelA.append(key2)
                      toDelB.append(b0)

    for key in toDelA:
        del a2b[key]
    for b0 in toDelB:
        del b2a[b0]

    return result 
     


def checkTermini(rnaLength,a2b,b2a,VerA,VerB):
    # to find out unsuccessful pairing at the first and second position of 5',3'-terminus. 
    # Assumption: in our current design procedure, the 5' terminus and 3' terminus should be paired 

    # situations to be covered -> [1-o,2-c], [1-o,2-o],  

    result = {}
    leftKey=""
    rightKey=""
    flag = 0
 
    for k in VerA:
        #1-o, 2-c
        if len(VerA[k]) == 2:
           if (rnaLength-1) in VerA[k]:
              if 2 in VerA[k]:
                 leftKey = k
                 flag = 1

        #1-o,2-o
        # to be added here
        if len(VerA[k]) == 2:
           if (rnaLength-2) in VerA[k]:
              if 3 in VerA[k]:
                leftKey = k
                flag = 1


    #print leftKey

 
    for k in VerB:
        if len(VerB[k]) == 2:
           if rnaLength in VerB[k]:
              if 1 in VerB[k]:
                 rightKey = k

    # We can skip to check the results in a2b and b2a
    # We should expect {leftKey:['0']} for a2b, {rightKey:['0']} for b2a


    #
    #print leftKey,rightKey 
    if leftKey != "" and rightKey != "":
       #print "Hi" 
       result[leftKey] = rightKey
       del a2b[leftKey]
       del b2a[rightKey]
 
    if flag != 0:
       pass
       #print "Unsuccessful pairing at termini detected."    

    return result


def findInterchangablePair(a2b,b2a):
    # Consider the situation of 
    # a2b -> 'a':['c'], 'b':['c','d']
    # b2a -> 'd':['b'], 'c':['a','b']

    # '0' may occur anywhere in the 4 loops

    result = {}
    toDelA = []
    toDelB = []

    go5 = 0

    # 
    for a in a2b:
        if (len(a2b[a]) == 1 or len(a2b[a]) == 2 ):
           go1 = 0
           if (len(a2b[a]) == 1):
               if a2b[a][0] != '0':
                  go1 = 1
                  c = a2b[a][0]
           else:
               if a2b[a][0] == '0':
                  go1 = 1
                  c = a2b[a][1]
               if a2b[a][1] == '0':
                  go1 = 1
                  c = a2b[a][0]
           #
           if go1 == 1:
              for b in a2b:
                if (len(a2b[b]) == 2 or len(a2b[b]) == 3 ):
                   go2 = 0
                   if len(a2b[b]) == 2:
                       if '0' not in a2b[b]:
                          if a2b[b][0] == c:
                             go2 = 1
                             d = a2b[b][1]
                          if a2b[b][1] == c:
                             go2 = 1
                             d = a2b[b][0]  
                   if len(a2b[b]) == 3:
                       if '0' in a2b[b]:
                          if a2b[b][0] == c and a2b[b][1] == '0':
                             go2 = 1
                             d = a2b[b][2]
                          if a2b[b][0] == '0' and a2b[b][1] == c:
                             go2 = 1
                             d = a2b[b][2]
   
                          if a2b[b][0] == c and a2b[b][2] == '0':
                             go2 = 1
                             d = a2b[b][1]
                          if a2b[b][0] == '0' and a2b[b][2] == c:
                             go2 = 1
                             d = a2b[b][1]
 
                          if a2b[b][1] == c and a2b[b][2] == '0':
                             go2 = 1 
                             d = a2b[b][0]
                          if a2b[b][1] == '0' and a2b[b][2] == c:
                             go2 = 1 
                             d = a2b[b][0]
                         
                   if go2 == 1:
                      go3 = 0
                      if d in b2a:
                         if len(b2a[d]) == 1 or len(b2a[d]) == 2:
                            if len(b2a[d]) == 1:
                               if b2a[d][0] == b:
                                  go3 = 1
                            if len(b2a[d]) == 2:
                               if b in b2a[d]:
                                  if '0' in b2a[d]:
                                     go3 = 1
 
                      if go3 == 1:
                         go4 = 0
                         if c in b2a:
                            if len(b2a[c]) == 2 or len(b2a[c]) == 3:
                               if len(b2a[c]) == 2:
                                  if a in b2a[c]:
                                     if b in b2a[c]:
                                        go4 = 1
                               if len(b2a[c]) == 3: 
                                  if a in b2a[c]:
                                     if b in b2a[c]:
                                        if '0' in b2a[c]:
                                           go4 = 1
                         if go4 == 1:
                            go5 = 1
                            toDelA.append(a)
                            toDelA.append(b)
                            toDelB.append(c)
                            toDelB.append(d)
                            result[a] = c
                            result[b] = d                   
   
    if go5 == 1: 
       for key in toDelA:
         del a2b[key]
       for b0 in toDelB:
         del b2a[b0]

    return result 
  
                    
 
              
                



    

#    for l1 in a2b:
#        if ( len(a2b[l1]) == 1 and a2b[l1][0] == l1 ) or \
#           ( len(a2b[l1]) == 2 and (l1 in a2b[l1]) and ('0' in a2b[l1]) ):
#           for l2 in a2b:
#               if ( len(a2b[l2]) == 2 and (l1 in a2b[l2]) and (l2 in a2b[l2]) ) or \
#                  ( len(a2b[l2]) == 3 and (l1 in a2b[l2]) and (l2 in a2b[l2]) and ('0' in a2b[l2]) ):
#                  # Switch to b2a
#                  if ( len(b2a[l2]) == 1 and b2a[l2][0] == l2 ) or \
#                     ( len(b2a[l2]) == 2 and (l2 in b2a[l2]) and ('0' in b2a[l2]) ):
#                     if ( len(b2a[l1]) == 2 and (l1 in b2a[l1]) and (l2 in b2a[l1]) ) or \
#                        ( len(b2a[l1]) == 3 and (l1 in b2a[l1]) and (l2 in b2a[l1]) and ('0' in b2a[l1]) ):
                           
 
    



  
    pass
  






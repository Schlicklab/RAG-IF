# utility for s2_script

import ast
import sys


def procDicStr(dicstring):
    # process dictionary string

    dic0 = {}

    if dicstring == "":
        pass
    else:
        dic0 = ast.literal_eval( dicstring )
    return dic0    




def extDicList(dic0,leftL,rightL):
    # convert dictionary to list information 
    checker = []
    if type(leftL) != type(checker):
       print "Error in extDicList.\nAbort."
       sys.exit()
    if type(rightL) != type(checker):
       print "Error in extDicList.\nAbort."
       sys.exit()
    
    #
    for left in dic0:
        right = dic0[left]
        if right not in rightL:
           rightL.append( right )
        if left not in leftL:
           leftL.append( left )


    return leftL,rightL
   

def procCr3StrList( cr3str, LeftV, RightV, toBreak ):
    if cr3str == []:
       pass
    else:
        for txt in cr3str:
            l1 = txt.split(',')[0].split('(')[1] 
            l2 = txt.split(',')[1].split(')')[0]
            r  = txt.split(',')[1].split('->')[1].strip()

            toBreak.append( [l1,l2] )

            if l1 not in LeftV:
               LeftV.append( l1 )
            if l2 not in LeftV:
               LeftV.append( l2 )
            if r  not in RightV:
               RightV.append( r )
    #
    return LeftV,RightV,toBreak



class residue:
    def __init__(self,no,pairNum,NTlabel,stemID,ss,loopID):
       self.no = no # Number 
       self.pairNum = pairNum # With which NT it is paired 
       self.NTlabel = NTlabel # The NT label - A U C G
       self.stemID = stemID   # Stem ID
       self.ss = ss           #
       self.loopID = loopID   # Loop ID
       self.ifmutate = 'no'

    def prtMut(self):
        print self.no, self.NTlabel  , self.stemID, self.loopID , self.ifmutate 



def printMutInf(resList):
    for i in range(len(resList)):
        resList[i].prtMut()

    return 0




#a="{'0':'1'}"



#b=procDicStr(a)
#print b['0']




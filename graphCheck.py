#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 22:31:46 2020

@author: qiyaozhu
"""
import random
import os
import os.path
import sys
import time
from functools import partial
import multiprocessing
from ClassesFunctions import *
from dualGraphs import *


def nupackCheck(survivor):
    
    seq = survivor[0]
    jobID = survivor[1]
    
    with open("tmpRNAfold"+jobID+".in","w") as f:
        f.write(">84 MT246482.1 0 13405 84\n")
        f.write(seq)
    
    os.system("mfe -pseudo -material rna tmpRNAfold"+jobID+" 2>/dev/null ")
    with open("tmpRNAfold"+jobID+".mfe", 'r') as f:
        fold = f.readlines()[16]
    with open("tmpRNAfold"+jobID+".mfe", 'w') as f:
        f.write(">84 MT246482.1 0 13405 84\n")
        f.write(seq + "\n")
        f.write(fold)
    os.system("dot2ct tmpRNAfold"+jobID+".mfe tmpRNAfold"+jobID+".ct")
    
    RNA = getCTInfo("tmpRNAfold"+jobID+".ct")
    os.system("rm -rf tmpRNAfold"+jobID+".in tmpRNAfold"+jobID+".ct tmpRNAfold"+jobID+".mfe")
    
    countHelices(RNA) 
    changeHelices(RNA)
    RNA.makeMatrices()
    connectHelices(RNA)
    for i in range(0,len(RNA.adjMatrix)): # S.J. 07/11/2018 - to keep track of vertexOrder
        vertexOrder.append(0)
        
    success, graph = calcEigen(RNA)
    correctHNumbers(RNA)
    if len(RNA.adjMatrix)==1 or len(RNA.adjMatrix)>9:
        print ("No matching graph exists because vertex number is either 1 or greater than 10.")
        return None
    elif success == 0: # no graph ID was assigned as eigen values not in the library S.J. 11/09/2017
        print ("No matching graph exists (even if the vertex number is between 2 and 9).")
        return None
    else:
        print(graph)
        return graph


def getSeqJob(heaven):
    
    Survivors = []
    
    with open(heaven, 'r') as f:
        lines = f.readlines()
        
    for i in range(len(lines)):
        if lines[i][0] != '>':
            continue
        else:
            seq = lines[i+1].split(" ")[0]            
            jobID = str(random.randint(10000,99999))
            jobID = jobID+str(time.time()).split('.')[1]
            list1= ['a','b','c','d','e','f','g','h','i','j']
            list2= [1,2,3,4,5,6,7,8,9,0]
            jobID = jobID + random.choice(list1)+ str(random.choice(list2))
            jobID = jobID + random.choice(list1)+ str(random.choice(list2))
            jobID = jobID + random.choice(list1)+ str(random.choice(list2))
            Survivors.append([seq, jobID])
            
    return Survivors
 
           
def main(heaven):
    
    Survivors = getSeqJob(heaven)
    pool = multiprocessing.Pool(4)
    result = pool.map( nupackCheck, Survivors )
    pool.close()
    pool.join()
    
    
if __name__== "__main__":

   if len(sys.argv) < 2:
      print("incorrect excution...")
      sys.exit()
   
   heaven = sys.argv[1]
   if not os.path.isfile(heaven):
      print("Heaven file not exist...")
      sys.exit()

   main(heaven)
    
            
            
    
    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 12:22:05 2020

@author: qiyaozhu
"""


import os

def withSlippery(inpf):
    
    design = inpf.split("i")[0]
    name = design + "Forna"
    
    with open("COV_PK.out", 'r') as f:
        seq = f.readlines()[1]
        
    with open(inpf, 'r') as f:
        lines = f.readlines()
        fold = lines[0]
        mutseq = lines[1]
    start = mutseq.find('N')
    num = mutseq.count('N')
        
    with open(name, 'w') as f:
        f.write(">Mutation design "+design+"\n")
        f.write(seq)
        f.write(fold)
        for i in range(num):
            start = start + 1
            f.write(str(start)+':red ')

def withoutSlippery(inpf):

    design = inpf.split("i")[0]
    name = design + "Forna"
    
    with open("COV_PK_noSlippery.out", 'r') as f:
        seq = f.readlines()[1]
        
    with open(inpf, 'r') as f:
        lines = f.readlines()
        fold = lines[0]
        mutseq = lines[1]
    start = mutseq.find('N')
    num = mutseq.count('N')
        
    with open(name, 'w') as f:
        f.write(">Mutation design "+design+"\n")
        f.write(seq)
        f.write(fold)
        for i in range(num):
            start = start + 1
            f.write(str(start)+':red ')
            
withoutSlippery('3_5_3inpf')

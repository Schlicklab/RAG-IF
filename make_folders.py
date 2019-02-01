# This script is designed to make folders after the second filter 
# Y. Tao -2018/07/12
#
import os
import os.path
import sys


FatherFolder = "7_4_test"
# this needs to be changed according to the target and the folder - this will be created

inpf = "result-filter.txt"

if not os.path.isdir( FatherFolder ):
   os.system("mkdir "+FatherFolder)

with open(inpf) as pointer:
     for line in pointer:
         #
         if len(line) > 1: 
            #
            subf = line.split()[-1]
            subsubf = line.split()[0]
            
            #
            if not os.path.isdir( FatherFolder+"/"+subf ):
               os.system("mkdir "+FatherFolder+"/"+subf)
            os.system("mkdir "+FatherFolder+"/"+subf+"/"+subsubf)

             




















